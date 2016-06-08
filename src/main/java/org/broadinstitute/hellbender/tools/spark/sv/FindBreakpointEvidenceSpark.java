package org.broadinstitute.hellbender.tools.spark.sv;

import com.esotericsoftware.kryo.DefaultSerializer;
import com.esotericsoftware.kryo.Kryo;
import com.esotericsoftware.kryo.io.Input;
import com.esotericsoftware.kryo.io.Output;
import com.google.cloud.dataflow.sdk.options.PipelineOptions;
import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.*;
import org.apache.commons.collections4.iterators.SingletonIterator;
import org.apache.spark.HashPartitioner;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.VoidFunction;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchHashSet;
import org.broadinstitute.hellbender.tools.spark.utils.MapPartitioner;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import scala.Tuple2;

import java.io.*;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Tool to describe reads that support a hypothesis of a genomic breakpoint.
 */
@CommandLineProgramProperties(summary="Find reads that evidence breakpoints.",
        oneLineSummary="Dump FASTQs for local assembly of putative genomic breakpoints.",
        programGroup = SparkProgramGroup.class)
public final class FindBreakpointEvidenceSpark extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    // The minimum mapping quality for reads used to gather evidence of breakpoints.
    private static final int MIN_MAPQ = 20;

    // The minimum length of the matched portion of an interesting alignment.  Reads that don't match at least
    // this many reference bases won't be used in gathering evidence.
    private static final int MIN_MATCH_LEN = 45;

    // Intervals with more than this much coverage are filtered out, because the reads mapped to that interval are
    // clearly not exclusively local to the interval.
    private static final int MAX_INTERVAL_COVERAGE = 1000;

    // A guess about how many reads will be in the final assembly, given the number reads mapped to the interval.
    // Just used to size a hashmap to a capacity that's not wildly inappropriate.  This was determined empirically.
    private static final float ASSEMBLY_TO_MAPPED_SIZE_RATIO = 7.f;

    @Argument(doc = "directory for fastq output", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private String outputDir;

    @Argument(doc = "directory for evidence output", fullName = "breakpointEvidenceDir", optional = true)
    private String evidenceDir;

    @Argument(doc = "file for breakpoint intervals output", fullName = "breakpointIntervals", optional = true)
    private String intervalFile;

    @Argument(doc = "file for mapped qname intervals output", fullName = "qnameIntervalsMapped", optional = true)
    private String qNamesMappedFile;

    @Argument(doc = "file for high frequency kmers output", fullName = "highFrequencyKmers", optional = true)
    private String highFrequencyKmersFile;

    @Argument(doc = "file for kmer intervals output", fullName = "kmerIntervals", optional = true)
    private String kmerFile;

    @Argument(doc = "file for mapped qname intervals output", fullName = "qnameIntervalsForAssembly", optional = true)
    private String qNamesAssemblyFile;

    /**
     * This is a path to a file of kmers that appear too frequently in the reference to be usable as probes to localize
     * reads.  We don't calculate it here, because it depends only on the reference.
     * The program FindBadGenomicKmersSpark can produce such a list for you.
     */
    @Argument(doc = "file containing ubiquitous kmer list", fullName = "kmersToIgnore")
    private String kmersToIgnoreFilename;

    @Override
    public boolean requiresReads()
    {
        return true;
    }

    @Override
    protected void runTool( final JavaSparkContext ctx ) {
        if ( getHeaderForReads().getSortOrder() != SAMFileHeader.SortOrder.coordinate ) {
            throw new GATKException("The reads must be coordinate sorted.");
        }

        final Locations locations =
                new Locations(evidenceDir, intervalFile, qNamesMappedFile,
                                highFrequencyKmersFile, kmerFile, qNamesAssemblyFile);
        final PipelineOptions pipelineOptions = getAuthenticatedGCSOptions();
        final JavaRDD<GATKRead> unfilteredReads = getUnfilteredReads();
        final JavaRDD<GATKRead> allPrimaryLines =
                unfilteredReads.filter(read -> !read.isSecondaryAlignment() && !read.isSupplementaryAlignment());
        final int[] nIntervals = new int[1];

        // develop evidence, intervals, and, finally, a set of template names for each interval
        final HopscotchHashSet<QNameAndInterval> qNamesSet =
                getMappedQNamesSet(ctx, getHeaderForReads(), unfilteredReads, locations, pipelineOptions, nIntervals);

        if ( nIntervals[0] == 0 ) return;

        // supplement the template names with other reads that share kmers
        addAssemblyQNames(ctx, kmersToIgnoreFilename, qNamesSet, allPrimaryLines, locations, pipelineOptions);

        // write a FASTQ file for each interval
        final String outDir = outputDir;
        generateFastqs(ctx, qNamesSet, nIntervals[0], allPrimaryLines,
                intervalAndFastqBytes -> writeFastq(intervalAndFastqBytes, outDir));
        log("Wrote FASTQs for assembly.");
    }

    /**
     * Find the breakpoint evidence,
     * cluster the evidence into intervals,
     * find the template names mapped into each interval,
     * kmerize each of these templates,
     * clean up by removing some intervals that are bogus as evidenced by ubiquitous kmers,
     * and return a set of template names and the intervals to which they belong.
     */
    private HopscotchHashSet<QNameAndInterval> getMappedQNamesSet(
            final JavaSparkContext ctx,
            final SAMFileHeader header,
            final JavaRDD<GATKRead> unfilteredReads,
            final Locations locations,
            final PipelineOptions pipelineOptions,
            final int[] nIntervalsArr )
    {
        final JavaRDD<GATKRead> mappedReads =
                unfilteredReads.filter(read ->
                        !read.isDuplicate() && !read.failsVendorQualityCheck() && !read.isUnmapped());
        final ReadMetadata readMetadata = getMetadata(header);
        log("Metadata retrieved.");
        final Broadcast<ReadMetadata> broadcastMetadata = ctx.broadcast(readMetadata);
        final List<Interval> intervals = getIntervals(broadcastMetadata, header, mappedReads, locations);

        final int nIntervals = intervals.size();
        nIntervalsArr[0] = intervals.size();
        log("Discovered " + nIntervals + " intervals.");

        if ( nIntervals == 0 ) return null;

        final List<Interval> goodIntervals =
                removeHighCoverageIntervals(ctx, broadcastMetadata, intervals, mappedReads);

        final int nKilledIntervals = intervals.size() - goodIntervals.size();
        log("Killed " + nKilledIntervals + " intervals that had >" + MAX_INTERVAL_COVERAGE + "x coverage.");

        // record the intervals
        if ( locations.intervalFile != null ) {
            try (final OutputStreamWriter writer = new OutputStreamWriter(new BufferedOutputStream(
                    BucketUtils.createFile(locations.intervalFile, pipelineOptions)))) {
                final List<SAMSequenceRecord> contigs = header.getSequenceDictionary().getSequences();
                for (final Interval interval : goodIntervals) {
                    final String seqName = contigs.get(interval.getContig()).getSequenceName();
                    writer.write(seqName + " " + interval.getStart() + " " + interval.getEnd() + "\n");
                }
            } catch (final IOException ioe) {
                throw new GATKException("Can't write intervals file " + locations.intervalFile, ioe);
            }
        }

        final HopscotchHashSet<QNameAndInterval> qNamesSet =
                getQNames(ctx, broadcastMetadata, goodIntervals, mappedReads);
        broadcastMetadata.destroy();

        if ( locations.qNamesMappedFile != null ) {
            dumpQNames(locations.qNamesMappedFile, pipelineOptions, qNamesSet);
        }
        log("Discovered " + qNamesSet.size() + " mapped template names.");

        return qNamesSet;
    }

    /**
     * Kmerize each read mapped into a breakpoint interval,
     * get the template names of all reads sharing these kmers (regardless of where or if they're mapped),
     * and add these template names to the set of names for each interval.
     */
    private void addAssemblyQNames(
            final JavaSparkContext ctx,
            final String kmersToIgnoreFilename,
            final HopscotchHashSet<QNameAndInterval> qNamesSet,
            final JavaRDD<GATKRead> allPrimaryLines,
            final Locations locations,
            final PipelineOptions pipelineOptions )
    {
        final JavaRDD<GATKRead> goodPrimaryLines =
                allPrimaryLines.filter(read -> !read.isDuplicate() && !read.failsVendorQualityCheck());

        qNamesSet.addAll(
                getAssemblyQNames(
                        ctx,
                        getKmerAndIntervalsSet(ctx, kmersToIgnoreFilename, qNamesSet, goodPrimaryLines, locations, pipelineOptions),
                        goodPrimaryLines));

        if ( locations.qNamesAssemblyFile != null ) {
            dumpQNames(locations.qNamesAssemblyFile, pipelineOptions, qNamesSet);
        }

        log("Discovered "+qNamesSet.size()+" unique template names for assembly.");
    }

    /**
     * Kmerize reads having template names in a given set,
     * filter out kmers that appear too often in this read set or in the genome to be helpful in localizing reads,
     * and return the set of kmers that appear in each interval.
     */
    private HopscotchHashSet<KmerAndInterval> getKmerAndIntervalsSet(
            final JavaSparkContext ctx,
            final String kmersToIgnoreFilename,
            final HopscotchHashSet<QNameAndInterval> qNamesSet,
            final JavaRDD<GATKRead> goodPrimaryLines,
            final Locations locations,
            final PipelineOptions pipelineOptions )
    {
        final Set<SVKmer> kmerKillSet = SVKmer.readKmersFile(kmersToIgnoreFilename, pipelineOptions);
        log("Ignoring " + kmerKillSet.size() + " genomically common kmers.");
        final List<SVKmer> kmerKillList = getHighCountKmers(goodPrimaryLines, locations, pipelineOptions);
        log("Ignoring " + kmerKillList.size() + " common kmers in the reads.");
        kmerKillSet.addAll(kmerKillList);
        log("Ignoring a total of " + kmerKillSet.size() + " unique common kmers.");

        final HopscotchHashSet<KmerAndInterval> kmerAndIntervalSet = new HopscotchHashSet<>(
                getKmerIntervals(ctx, qNamesSet, kmerKillSet, goodPrimaryLines, locations, pipelineOptions));
        log("Discovered " + kmerAndIntervalSet.size() + " kmers.");

        return kmerAndIntervalSet;
    }

    /**
     * Transform all the reads for a supplied set of template names in each inverval into FASTQ records
     * for each interval, and do something with the list of FASTQ records for each interval (like write it to a file).
     */
    @VisibleForTesting static void generateFastqs(final JavaSparkContext ctx,
                                       final HopscotchHashSet<QNameAndInterval> qNamesSet,
                                       final int nIntervals,
                                       final JavaRDD<GATKRead> reads,
                                       final VoidFunction<Tuple2<Integer, List<byte[]>>> fastqWriter) {
        final Broadcast<HopscotchHashSet<QNameAndInterval>> broadcastQNamesSet = ctx.broadcast(qNamesSet);
        final int nPartitions = reads.partitions().size();

        reads
            .mapPartitionsToPair(readItr ->
                    new ReadsForQNamesFinder(broadcastQNamesSet.value(), nIntervals).call(readItr), false)
            .combineByKey(x -> x, FindBreakpointEvidenceSpark::combineLists, FindBreakpointEvidenceSpark::combineLists,
                        new HashPartitioner(nPartitions), false, null)
            .foreach(fastqWriter);

        broadcastQNamesSet.destroy();
    }

    /** Concatenate two lists. */
    private static List<byte[]> combineLists( final List<byte[]> list1, final List<byte[]> list2 ) {
        final List<byte[]> result = new ArrayList<>(list1.size() + list2.size());
        result.addAll(list1);
        result.addAll(list2);
        return result;
    }

    /**
     * Sort a list of FASTQ records.
     * This puts them into proper order for an interleaved FASTQ file (since the FASTQ record begins with the
     * template name).
     */
    @VisibleForTesting static byte[][] sortFastqRecs( final List<byte[]> fastqRecs ) {
        final byte[][] fastqsArray = new byte[fastqRecs.size()][];
        fastqRecs.toArray(fastqsArray);

        // since the fastq bytes start with the read name, this will interleave pairs
        Arrays.sort(fastqsArray, FindBreakpointEvidenceSpark::compareByteArrays);

        return fastqsArray;
    }

    /** write a FASTQ file for an assembly */
    private static void writeFastq( final Tuple2<Integer, List<byte[]>> intervalAndFastqs,
                                    final String outputDir ) {
        final byte[][] fastqsArray = sortFastqRecs(intervalAndFastqs._2);

        // we're not going to try to marshal PipelineOptions for now -- not sure it's a good idea, anyway
        final PipelineOptions pipelineOptions = null;
        final String fileName = outputDir + "/assembly" + intervalAndFastqs._1 + ".fastq";
        try ( final OutputStream writer = new BufferedOutputStream(BucketUtils.createFile(fileName, pipelineOptions)) ) {
            for ( final byte[] fastqBytes : fastqsArray ) {
                writer.write(fastqBytes);
            }
        } catch ( final IOException ioe ) {
            throw new GATKException("Can't write "+fileName, ioe);
        }
    }

    /** Sure seems like this should be available in the standard Java library somewhere, but I can't find it. */
    private static int compareByteArrays( final byte[] arr1, final byte[] arr2 ) {
        final int len = Math.min(arr1.length, arr2.length);
        for ( int idx = 0; idx != len; ++idx ) {
            final int result = Integer.compare(arr1[idx] & 0xff, arr2[idx] & 0xff);
            if ( result != 0 ) return result;
        }
        return Integer.compare(arr1.length, arr2.length);
    }

    /**
     * Grab template names for all reads that contain kmers associated with a given breakpoint.
     */
    @VisibleForTesting static List<QNameAndInterval> getAssemblyQNames(
            final JavaSparkContext ctx,
            final HopscotchHashSet<KmerAndInterval> kmerAndIntervalSet,
            final JavaRDD<GATKRead> reads ) {
        final Broadcast<HopscotchHashSet<KmerAndInterval>> broadcastKmerAndIntervalSet =
                ctx.broadcast(kmerAndIntervalSet);

        final List<QNameAndInterval> qNames =
            reads
                .mapPartitions(readItr ->
                        new MapPartitioner<>(readItr,
                                new QNamesForKmersFinder(broadcastKmerAndIntervalSet.value())), false)
                .collect();

        broadcastKmerAndIntervalSet.destroy();

        return qNames;
    }

    /** find kmers for each interval */
    @VisibleForTesting static List<KmerAndInterval> getKmerIntervals(
            final JavaSparkContext ctx,
            final HopscotchHashSet<QNameAndInterval> qNamesSet,
            final Set<SVKmer> kmerKillSet,
            final JavaRDD<GATKRead> reads,
            final Locations locations,
            final PipelineOptions pipelineOptions ) {

        final Broadcast<Set<SVKmer>> broadcastKmerKillSet = ctx.broadcast(kmerKillSet);
        final Broadcast<HopscotchHashSet<QNameAndInterval>> broadcastQNameAndIntervalsSet =
                ctx.broadcast(qNamesSet);

        // given a set of template names with interval IDs and a kill set of ubiquitous kmers,
        // produce a set of interesting kmers for each interval ID
        final List<KmerAndInterval> kmerIntervals =
            reads
                .mapPartitionsToPair(readItr ->
                        new MapPartitioner<>(readItr,
                            new QNameKmerizer(broadcastQNameAndIntervalsSet.value(),
                                            broadcastKmerKillSet.value())), false)
                .reduceByKey(Integer::sum)
                .mapPartitions(itr -> new KmerCleaner().call(itr))
                .collect();

        broadcastQNameAndIntervalsSet.destroy();
        broadcastKmerKillSet.destroy();

        // record the kmers with their interval IDs
        if ( locations.kmerFile != null ) {
            try (final OutputStreamWriter writer = new OutputStreamWriter(new BufferedOutputStream(
                    BucketUtils.createFile(locations.kmerFile, pipelineOptions)))) {
                for (final KmerAndInterval kmerAndInterval : kmerIntervals) {
                    writer.write(kmerAndInterval.toString() + "\n");
                }
            } catch (final IOException ioe) {
                throw new GATKException("Can't write kmer intervals file " + locations.kmerFile, ioe);
            }
        }

        return kmerIntervals;
    }

    /** Find all the kmers in this read set that occur with high (KmerReducer.MAX_KMER_COUNT) frequency. */
    @VisibleForTesting static List<SVKmer> getHighCountKmers(
            final JavaRDD<GATKRead> reads,
            final Locations locations,
            final PipelineOptions pipelineOptions ) {
        final int nPartitions = reads.partitions().size();
        final List<SVKmer> kmers =
                reads
                    .mapPartitions(KmerCounter::new, false)
                    .mapToPair(kmerAndCount -> new Tuple2<>(kmerAndCount, null))
                    .partitionBy(new HashPartitioner(nPartitions))
                    .map(tuple -> tuple._1)
                    .mapPartitions(KmerReducer::new)
                    .map(SVKmer::new)
                    .collect();

        if ( locations.highFrequencyKmersFile != null ) {
            final String fileName = locations.highFrequencyKmersFile;
            try ( final Writer writer =
                          new BufferedWriter(new OutputStreamWriter(BucketUtils.createFile(fileName, pipelineOptions))) ) {
                for ( final SVKmer kmer : kmers ) {
                    writer.write(kmer.toString(SVConstants.KMER_SIZE));
                    writer.write('\n');
                }
            }
            catch ( final IOException ioe ) {
                throw new GATKException("Unable to write high frequency kmers file.", ioe);
            }
        }
        return kmers;
    }

    /** figure out the coverage for each interval and filter out those with ridiculously high coverage */
    @VisibleForTesting static List<Interval> removeHighCoverageIntervals(
            final JavaSparkContext ctx,
            final Broadcast<ReadMetadata> broadcastMetadata,
            final List<Interval> intervals,
            final JavaRDD<GATKRead> reads) {
        final Broadcast<List<Interval>> broadcastIntervals = ctx.broadcast(intervals);
        final Map<Integer, Integer> intervalCoverage =
                reads
                    .mapPartitionsToPair(readItr ->
                            new IntervalCoverageFinder(broadcastMetadata.value(),broadcastIntervals.value(),readItr))
                    .reduceByKey(Integer::sum)
                    .collectAsMap();
        final List<Interval> result =
                IntStream.range(0, intervals.size())
                .filter(idx -> intervalCoverage.getOrDefault(idx, 0)/intervals.get(idx).getLength() <= MAX_INTERVAL_COVERAGE)
                .mapToObj(intervals::get)
                .collect(Collectors.toList());
        broadcastIntervals.destroy();
        return result;
    }

    /** find template names for reads mapping to each interval */
    @VisibleForTesting static HopscotchHashSet<QNameAndInterval> getQNames(
                                    final JavaSparkContext ctx,
                                    final Broadcast<ReadMetadata> broadcastMetadata,
                                    final List<Interval> intervals,
                                    final JavaRDD<GATKRead> reads ) {
        final Broadcast<List<Interval>> broadcastIntervals = ctx.broadcast(intervals);
        final List<QNameAndInterval> qNameAndIntervalList =
                reads
                        .mapPartitions(readItr ->
                                new MapPartitioner<>(readItr,
                                        new QNameFinder(broadcastMetadata.value(),
                                                broadcastIntervals.value())), false)
                        .collect();
        broadcastIntervals.destroy();

        final HopscotchHashSet<QNameAndInterval> result =
                new HopscotchHashSet<>((int)(ASSEMBLY_TO_MAPPED_SIZE_RATIO*qNameAndIntervalList.size()));
        result.addAll(qNameAndIntervalList);
        return result;
    }

    /** write template names and interval IDs to a file. */
    private static void dumpQNames( final String qNameFile,
                                    final PipelineOptions pipelineOptions,
                                    final Collection<QNameAndInterval> qNames ) {
        try (final OutputStreamWriter writer = new OutputStreamWriter(new BufferedOutputStream(
                BucketUtils.createFile(qNameFile, pipelineOptions)))) {
            for (final QNameAndInterval qnameAndInterval : qNames) {
                writer.write(qnameAndInterval.toString() + "\n");
            }
        } catch (final IOException ioe) {
            throw new GATKException("Can't write qname intervals file " + qNameFile, ioe);
        }
    }

    /**
     * Identify funky reads that support a hypothesis of a breakpoint in the vicinity, group the reads,
     * and declare a breakpoint interval where there is sufficient density of evidence.
     */
    @VisibleForTesting static List<Interval> getIntervals( final Broadcast<ReadMetadata> broadcastMetadata,
                                                final SAMFileHeader header,
                                                final JavaRDD<GATKRead> reads,
                                                final Locations locations ) {
        // find all breakpoint evidence, then filter for pile-ups
        final int maxFragmentSize = broadcastMetadata.value().getMaxMedianFragmentSize();
        final List<SAMSequenceRecord> contigs = header.getSequenceDictionary().getSequences();
        final int nContigs = contigs.size();
        final JavaRDD<BreakpointEvidence> evidenceRDD =
                reads
                    .filter(read ->
                            read.getMappingQuality() >= MIN_MAPQ &&
                            read.getCigar().getCigarElements()
                                    .stream()
                                    .filter(ele -> ele.getOperator()==CigarOperator.MATCH_OR_MISMATCH ||
                                                    ele.getOperator()==CigarOperator.EQ ||
                                                    ele.getOperator()==CigarOperator.X)
                                    .mapToInt(CigarElement::getLength)
                                    .sum() >= MIN_MATCH_LEN)
                    .mapPartitions(readItr ->
                            new MapPartitioner<>(readItr, new ReadClassifier(broadcastMetadata.value())), true)
                    .mapPartitions(evidenceItr ->
                            new MapPartitioner<>(evidenceItr, new BreakpointClusterer(2*maxFragmentSize)), true)
                    .mapPartitions(evidenceItr ->
                            new MapPartitioner<>(evidenceItr,
                                    new WindowSorter(3*maxFragmentSize), new BreakpointEvidence(nContigs)), true);

        // record the evidence
        if ( locations.evidenceDir != null ) {
            evidenceRDD.cache();
            evidenceRDD.saveAsTextFile(locations.evidenceDir);
        }

        // find discrete intervals that contain the breakpoint evidence
        final Iterator<Interval> intervalItr =
                evidenceRDD
                        .mapPartitions(evidenceItr ->
                                new MapPartitioner<>(evidenceItr,
                                        new EvidenceToIntervalMapper(maxFragmentSize),
                                        new BreakpointEvidence(nContigs)), true)
                        .collect()
                        .iterator();

        // coalesce overlapping intervals (can happen at partition boundaries)
        final List<Interval> intervals = new ArrayList<>();
        if ( intervalItr.hasNext() ) {
            Interval prev = intervalItr.next();
            while ( intervalItr.hasNext() ) {
                final Interval next = intervalItr.next();
                if ( prev.isDisjointFrom(next) ) {
                    intervals.add(prev);
                    prev = next;
                } else {
                    prev = prev.join(next);
                }
            }
            intervals.add(prev);
        }

        if ( locations.evidenceDir != null ) evidenceRDD.unpersist();

        return intervals;
    }

    /** gather some interesting factoids about the reads in aggregate */
    @VisibleForTesting static ReadMetadata getMetadata( final SAMFileHeader header ) {
        final List<SAMReadGroupRecord> groups = header.getReadGroups();
        final int nGroups = groups.size();
        final ReadMetadata.ReadGroupFragmentStatistics groupStats =
                //TODO: get real data
                new ReadMetadata.ReadGroupFragmentStatistics(400.f, 75.f);
        final List<ReadMetadata.ReadGroupFragmentStatistics> stats = new ArrayList<>(nGroups);
        for ( int idx = 0; idx != nGroups; ++idx ) {
            stats.add(groupStats);
        }
        return new ReadMetadata(header, stats, groupStats);
    }

    private void log( final String message ) {
        logger.info(message);
    }

    @VisibleForTesting static class Locations {
        public final String evidenceDir;
        public final String intervalFile;
        public final String qNamesMappedFile;
        public final String highFrequencyKmersFile;
        public final String kmerFile;
        public final String qNamesAssemblyFile;

        public Locations( final String evidenceDir, final String intervalFile,
                          final String qNamesMappedFile, final String highFrequencyKmersFile,
                          final String kmerFile, final String qNamesAssemblyFile ) {
            this.evidenceDir = evidenceDir;
            this.intervalFile = intervalFile;
            this.qNamesMappedFile = qNamesMappedFile;
            this.highFrequencyKmersFile = highFrequencyKmersFile;
            this.kmerFile = kmerFile;
            this.qNamesAssemblyFile = qNamesAssemblyFile;
        }
    }

    /**
     * A class that acts as a filter for breakpoint evidence.
     * It passes only that evidence that is part of a putative cluster.
     */
    private static final class BreakpointClusterer
            implements Function<BreakpointEvidence, Iterator<BreakpointEvidence>> {
        private final int staleEventDistance;
        private final SortedMap<BreakpointEvidence, Boolean> locMap = new TreeMap<>();
        private final List<Map.Entry<BreakpointEvidence, Boolean>> reportableEntries = new ArrayList<>(2*MIN_EVIDENCE);
        private final Iterator<BreakpointEvidence> noEvidence = Collections.emptyIterator();
        private int currentContig = -1;

        private static final int MIN_EVIDENCE = 15; // minimum evidence count in a cluster

        BreakpointClusterer( final int staleEventDistance ) {
            this.staleEventDistance = staleEventDistance;
        }

        public Iterator<BreakpointEvidence> apply( final BreakpointEvidence evidence ) {
            if ( evidence.getContigIndex() != currentContig ) {
                currentContig = evidence.getContigIndex();
                locMap.clear();
            }

            locMap.put(evidence, true);

            final int locusStart = evidence.getEventStartPosition();
            final int locusEnd = evidence.getContigEnd();
            final int staleEnd = locusStart - staleEventDistance;
            int evidenceCount = 0;
            reportableEntries.clear();
            final Iterator<Map.Entry<BreakpointEvidence, Boolean>> itr = locMap.entrySet().iterator();
            while ( itr.hasNext() ) {
                final Map.Entry<BreakpointEvidence, Boolean> entry = itr.next();
                final BreakpointEvidence evidence2 = entry.getKey();
                final int contigEnd = evidence2.getContigEnd();
                if ( contigEnd <= staleEnd ) itr.remove();
                else if ( evidence2.getEventStartPosition() >= locusEnd ) break;
                else if ( contigEnd > locusStart ) {
                    evidenceCount += 1;
                    if ( entry.getValue() ) reportableEntries.add(entry);
                }
            }

            if ( evidenceCount >= MIN_EVIDENCE ) {
                return reportableEntries.stream()
                        .map(entry -> { entry.setValue(false); return entry.getKey(); })
                        .iterator();
            }
            return noEvidence;
        }
    }

    /**
     * Class to fully sort a stream of nearly sorted BreakpointEvidences.
     */
    private static final class WindowSorter
            implements Function<BreakpointEvidence, Iterator<BreakpointEvidence>> {
        private final SortedSet<BreakpointEvidence> recordSet = new TreeSet<>();
        private final List<BreakpointEvidence> reportableEvidence = new ArrayList<>();
        private final int windowSize;
        private int currentContig = -1;

        WindowSorter( final int windowSize ) {
            this.windowSize = windowSize;
        }

        public Iterator<BreakpointEvidence> apply( final BreakpointEvidence evidence ) {
            reportableEvidence.clear();
            if ( evidence.getContigIndex() != currentContig ) {
                reportableEvidence.addAll(recordSet);
                recordSet.clear();
                currentContig = evidence.getContigIndex();
            } else {
                final int reportableEnd = evidence.getEventStartPosition() - windowSize;
                final Iterator<BreakpointEvidence> itr = recordSet.iterator();
                while ( itr.hasNext() ) {
                    final BreakpointEvidence evidence2 = itr.next();
                    if ( evidence2.getEventStartPosition() >= reportableEnd ) break;
                    reportableEvidence.add(evidence2);
                    itr.remove();
                }
            }
            recordSet.add(evidence);
            return reportableEvidence.iterator();
        }
    }

    /**
     * Minimalistic simple interval.
     */
    @DefaultSerializer(Interval.Serializer.class)
    @VisibleForTesting static final class Interval {
        private final int contig;
        private final int start;
        private final int end;

        Interval( final int contig, final int start, final int end ) {
            this.contig = contig;
            this.start = start;
            this.end = end;
        }

        private Interval( final Kryo kryo, final Input input ) {
            contig = input.readInt();
            start = input.readInt();
            end = input.readInt();
        }

        private void serialize( final Kryo kryo, final Output output ) {
            output.writeInt(contig);
            output.writeInt(start);
            output.writeInt(end);
        }

        public int getContig() { return contig; }
        public int getStart() { return start; }
        public int getEnd() { return end; }
        public int getLength() { return end-start; }

        public boolean isDisjointFrom( final Interval that ) {
            return this.contig != that.contig || this.end < that.start || that.end < this.start;
        }

        public int overlapLen( final Interval that ) {
            if ( this.contig != that.contig ) return 0;
            return Math.max(0, Math.min(this.end, that.end) - Math.max(this.start, that.start));
        }

        public Interval join( final Interval that ) {
            if ( this.contig != that.contig ) throw new GATKException("Joining across contigs.");
            return new Interval(contig, Math.min(this.start, that.start), Math.max(this.end, that.end));
        }

        @Override
        public boolean equals( final Object obj ) {
            if ( !(obj instanceof Interval) ) return false;
            final Interval that = (Interval)obj;
            return this.contig == that.contig && this.start == that.start && this.end == that.end;
        }

        @Override
        public int hashCode() {
            return 47*(47*(47*(47*contig)+start)+end);
        }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<Interval> {
            @Override
            public void write( final Kryo kryo, final Output output, final Interval interval ) {
                interval.serialize(kryo, output);
            }

            @Override
            public Interval read( final Kryo kryo, final Input input, final Class<Interval> klass ) {
                return new Interval(kryo, input);
            }
        }
    }

    /**
     * A class to examine a stream of BreakpointEvidence, and group it into Intervals.
     */
    private static final class EvidenceToIntervalMapper implements Function<BreakpointEvidence, Iterator<Interval>> {
        private final Iterator<Interval> noInterval = Collections.emptyIterator();
        private final int gapSize;
        private int contig = -1;
        private int start;
        private int end;

        EvidenceToIntervalMapper(final int gapSize ) {
            this.gapSize = gapSize;
        }

        public Iterator<Interval> apply( final BreakpointEvidence evidence ) {
            Iterator<Interval> result = noInterval;
            if ( evidence.getContigIndex() != contig ) {
                if ( contig != -1 ) {
                    result = new SingletonIterator<>(new Interval(contig, start, end));
                }
                contig = evidence.getContigIndex();
                start = evidence.getEventStartPosition();
                end = evidence.getContigEnd();
            } else if ( evidence.getEventStartPosition() >= end+gapSize ) {
                result = new SingletonIterator<>(new Interval(contig, start, end));
                start = evidence.getEventStartPosition();
                end = evidence.getContigEnd();
            } else {
                end = Math.max(end, evidence.getContigEnd());
            }
            return result;
        }
    }

    /**
     * Class to find the coverage of the intervals.
     */
    private static final class IntervalCoverageFinder implements Iterable<Tuple2<Integer, Integer>> {
        private final int[] basesInInterval;

        IntervalCoverageFinder( final ReadMetadata metadata,
                                final List<Interval> intervals,
                                final Iterator<GATKRead> readItr ) {
            basesInInterval = new int[intervals.size()];
            int intervalsIndex = 0;
            while ( readItr.hasNext() ) {
                final GATKRead read = readItr.next();
                final int readContigId = metadata.getContigID(read.getContig());
                final int readStart = read.getUnclippedStart();
                final int intervalsSize = intervals.size();
                while (intervalsIndex < intervalsSize) {
                    final Interval interval = intervals.get(intervalsIndex);
                    if (interval.getContig() > readContigId) break;
                    if (interval.getContig() == readContigId && interval.getEnd() > read.getStart()) break;
                    intervalsIndex += 1;
                }
                if (intervalsIndex >= intervalsSize) break;
                final Interval indexedInterval = intervals.get(intervalsIndex);
                final Interval readInterval = new Interval(readContigId, readStart, read.getUnclippedEnd());
                basesInInterval[intervalsIndex] += indexedInterval.overlapLen(readInterval);
            }
        }

        public Iterator<Tuple2<Integer,Integer>> iterator() {
            return IntStream
                    .range(0, basesInInterval.length)
                    .filter(idx -> basesInInterval[idx] > 0)
                    .mapToObj(idx -> new Tuple2<>(idx, basesInInterval[idx]))
                    .iterator();
        }
    }

    /**
     * A template name and an intervalId.
     * Note:  hashCode does not depend on intervalId, and that's on purpose.
     * This is actually a compacted (K,V) pair, and the hashCode is on K only.
     * Using a MultiMap of these is more memory-efficient than the more natural Tuple2(qName,List(intervalId)) because
     * most qNames are associated with a single interval.
     */
    @DefaultSerializer(QNameAndInterval.Serializer.class)
    @VisibleForTesting static final class QNameAndInterval {
        private final byte[] qName;
        private final int hashVal;
        private final int intervalId;

        QNameAndInterval( final String qName, final int intervalId ) {
            this.qName = qName.getBytes();
            this.hashVal = qName.hashCode();
            this.intervalId = intervalId;
        }

        private QNameAndInterval( final Kryo kryo, final Input input ) {
            final int nameLen = input.readInt();
            qName = input.readBytes(nameLen);
            hashVal = input.readInt();
            intervalId = input.readInt();
        }

        private void serialize( final Kryo kryo, final Output output ) {
            output.writeInt(qName.length);
            output.writeBytes(qName);
            output.writeInt(hashVal);
            output.writeInt(intervalId);
        }

        public String getQName() { return new String(qName); }
        public int getIntervalId() { return intervalId; }

        @Override
        public int hashCode() { return hashVal; }

        @Override
        public boolean equals( final Object obj ) {
            return obj instanceof QNameAndInterval && equals((QNameAndInterval) obj);
        }

        public boolean equals( final QNameAndInterval that ) {
            return this.intervalId == that.intervalId && Arrays.equals(this.qName, that.qName);
        }

        public boolean sameName( final byte[] name ) { return Arrays.equals(qName, name); }

        public String toString() { return new String(qName)+" "+intervalId; }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<QNameAndInterval> {
            @Override
            public void write( final Kryo kryo, final Output output, final QNameAndInterval qNameAndInterval ) {
                qNameAndInterval.serialize(kryo, output);
            }

            @Override
            public QNameAndInterval read( final Kryo kryo, final Input input, final Class<QNameAndInterval> klass ) {
                return new QNameAndInterval(kryo, input);
            }
        }
    }

    /**
     * Class to find the template names associated with reads in specified intervals.
     */
    private static final class QNameFinder implements Function<GATKRead, Iterator<QNameAndInterval>> {
        private final ReadMetadata metadata;
        private final List<Interval> intervals;
        private final Iterator<QNameAndInterval> noName = Collections.emptyIterator();
        private int intervalsIndex = 0;

        QNameFinder( final ReadMetadata metadata,
                     final List<Interval> intervals ) {
            this.metadata = metadata;
            this.intervals = intervals;
        }

        @Override
        public Iterator<QNameAndInterval> apply( final GATKRead read ) {
            final int readContigId = metadata.getContigID(read.getContig());
            final int readStart = read.getUnclippedStart();
            final int intervalsSize = intervals.size();
            while ( intervalsIndex < intervalsSize ) {
                final Interval interval = intervals.get(intervalsIndex);
                if ( interval.getContig() > readContigId ) break;
                if ( interval.getContig() == readContigId && interval.getEnd() > read.getStart() ) break;
                intervalsIndex += 1;
            }
            if ( intervalsIndex >= intervalsSize ) return noName;
            final Interval indexedInterval = intervals.get(intervalsIndex);
            final Interval readInterval = new Interval(readContigId, readStart, read.getUnclippedEnd());
            if ( indexedInterval.isDisjointFrom(readInterval) ) return noName;
            return new SingletonIterator<>(new QNameAndInterval(read.getName(), intervalsIndex));
        }
    }

    /**
     * A <Kmer,IntervalId> pair.
     * Note:  hashCode does not depend on intervalId, and that's on purpose.
     * This is actually a compacted (K,V) pair, and the hashCode is on K only.
     * Using a MultiMap of these is more memory-efficient than the more natural Tuple2(SVKmer,List(intervalId)) because
     * most kmers are associated with a single interval.
     */
    @DefaultSerializer(KmerAndInterval.Serializer.class)
    @VisibleForTesting final static class KmerAndInterval extends SVKmer {
        private final int intervalId;

        KmerAndInterval(final SVKmer kmer, final int intervalId ) {
            super(kmer);
            this.intervalId = intervalId;
        }

        private KmerAndInterval(final Kryo kryo, final Input input ) {
            super(kryo, input);
            intervalId = input.readInt();
        }

        protected void serialize( final Kryo kryo, final Output output ) {
            super.serialize(kryo, output);
            output.writeInt(intervalId);
        }

        @Override
        public boolean equals( final Object obj ) {
            return obj instanceof KmerAndInterval && equals((KmerAndInterval)obj);
        }

        public boolean equals( final KmerAndInterval that ) {
            return super.equals(that) && this.intervalId == that.intervalId;
        }

        public int getIntervalId() { return intervalId; }

        @Override
        public String toString() { return super.toString(SVConstants.KMER_SIZE)+" "+intervalId; }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<KmerAndInterval> {
            @Override
            public void write( final Kryo kryo, final Output output, final KmerAndInterval kmerAndInterval) {
                kmerAndInterval.serialize(kryo, output);
            }

            @Override
            public KmerAndInterval read(final Kryo kryo, final Input input,
                                        final Class<KmerAndInterval> klass ) {
                return new KmerAndInterval(kryo, input);
            }
        }
    }

    /**
     * Class that acts as a mapper from a stream of reads to a stream of KmerAndIntervals.
     * The template names of reads to kmerize, along with a set of kmers to ignore are passed in (by broadcast).
     */
    private static final class QNameKmerizer implements Function<GATKRead, Iterator<Tuple2<KmerAndInterval, Integer>>> {
        private final HopscotchHashSet<QNameAndInterval> qNameAndIntervalSet;
        private final Set<SVKmer> kmersToIgnore;
        private final ArrayList<Tuple2<KmerAndInterval, Integer>> tupleList = new ArrayList<>();

        QNameKmerizer( final HopscotchHashSet<QNameAndInterval> qNameAndIntervalSet,
                       final Set<SVKmer> kmersToIgnore ) {
            this.qNameAndIntervalSet = qNameAndIntervalSet;
            this.kmersToIgnore = kmersToIgnore;
        }

        public Iterator<Tuple2<KmerAndInterval, Integer>> apply( final GATKRead read ) {
            final String qName = read.getName();
            final int qNameHash = qName.hashCode();
            final byte[] qNameBytes = qName.getBytes();
            final Iterator<QNameAndInterval> names = qNameAndIntervalSet.bucketIterator(qName);
            tupleList.clear();
            while ( names.hasNext() ) {
                final QNameAndInterval qNameAndInterval = names.next();
                if ( qNameAndInterval.hashCode() == qNameHash && qNameAndInterval.sameName(qNameBytes) ) {
                    SVKmerizer.stream(read.getBases(), SVConstants.KMER_SIZE)
                            .map(kmer -> kmer.canonical(SVConstants.KMER_SIZE))
                            .filter(kmer -> !kmersToIgnore.contains(kmer))
                            .map(kmer -> new KmerAndInterval(kmer, qNameAndInterval.getIntervalId()))
                            .forEach(kmerCountAndInterval -> tupleList.add(new Tuple2<>(kmerCountAndInterval,1)));
                    break;
                }
            }
            return tupleList.iterator();
        }
    }


    /**
     * A <Kmer,count> pair.
     * Note:  hashCode and equality does not depend on count, and that's on purpose.
     * This is actually a compacted (K,V) pair, and the hashCode and equality is on K only.
     */
    @DefaultSerializer(KmerAndCount.Serializer.class)
    private final static class KmerAndCount extends SVKmer {
        private int count;

        KmerAndCount(final SVKmer kmer ) {
            super(kmer);
            this.count = 0;
        }

        private KmerAndCount(final Kryo kryo, final Input input ) {
            super(kryo, input);
            count = input.readInt();
        }

        protected void serialize( final Kryo kryo, final Output output ) {
            super.serialize(kryo, output);
            output.writeInt(count);
        }

        @Override
        public boolean equals( final Object obj ) {
            return obj instanceof KmerAndCount && equals((KmerAndCount)obj);
        }

        public boolean equals( final KmerAndCount that ) {
            return super.equals(that);
        }

        public int getCount() { return count; }
        public void incrementCount() { count += 1; }
        public void incrementCount( final int val ) { count += val; }

        @Override
        public String toString() { return super.toString(SVConstants.KMER_SIZE)+" "+count; }

        public static final class Serializer extends com.esotericsoftware.kryo.Serializer<KmerAndCount> {
            @Override
            public void write( final Kryo kryo, final Output output, final KmerAndCount kmerAndCount) {
                kmerAndCount.serialize(kryo, output);
            }

            @Override
            public KmerAndCount read(final Kryo kryo, final Input input, final Class<KmerAndCount> klass ) {
                return new KmerAndCount(kryo, input);
            }
        }
    }

    /**
     * Kmerizes reads into a big Set, ignores those that occur only once (i.e., < MIN_KMER_COUNT),
     * and returns a Kmer and occurrence count for the rest.
     */
    private final static class KmerCounter implements Iterable<KmerAndCount> {
        private static final int KMERS_PER_PARTITION_GUESS = 20000000;
        private static final int MIN_KMER_COUNT = 2;
        private final HopscotchHashSet<KmerAndCount> kmerSet;

        KmerCounter( final Iterator<GATKRead> readItr ) {
            kmerSet = new HopscotchHashSet<>(KMERS_PER_PARTITION_GUESS);
            while ( readItr.hasNext() ) {
                final Iterator<SVKmer> kmers = new SVKmerizer(readItr.next().getBases(),SVConstants.KMER_SIZE);
                while ( kmers.hasNext() ) {
                    kmerSet.put(new KmerAndCount(kmers.next().canonical(SVConstants.KMER_SIZE))).incrementCount();
                }
            }
            final Iterator<KmerAndCount> kmerItr = kmerSet.iterator();
            while ( kmerItr.hasNext() ) {
                if ( kmerItr.next().getCount() < MIN_KMER_COUNT ) kmerItr.remove();
            }
        }

        @Override
        public Iterator<KmerAndCount> iterator() { return kmerSet.iterator(); }
    }

    /**
     * Returns kmers that occur very frequently (>= MAX_KMER_COUNT).
     */
    private final static class KmerReducer implements Iterable<KmerAndCount> {
        private static final int KMERS_PER_PARTITION_GUESS = 250000;
        private static final int MAX_KMER_COUNT = 200;
        private final HopscotchHashSet<KmerAndCount> kmerSet;

        KmerReducer( final Iterator<KmerAndCount> kmerItr ) {
            kmerSet = new HopscotchHashSet<>(KMERS_PER_PARTITION_GUESS);
            while ( kmerItr.hasNext() ) {
                final KmerAndCount kmerAndCount = kmerItr.next();
                final KmerAndCount tableKmerAndCount = kmerSet.put(kmerAndCount);
                if ( kmerAndCount != tableKmerAndCount ) tableKmerAndCount.incrementCount(kmerAndCount.getCount());
            }
            final Iterator<KmerAndCount> kmerItr2 = kmerSet.iterator();
            while ( kmerItr2.hasNext() ) {
                if ( kmerItr2.next().getCount() < MAX_KMER_COUNT ) kmerItr2.remove();
            }
        }

        @Override
        public Iterator<KmerAndCount> iterator() { return kmerSet.iterator(); }
    }

    /**
     * Eliminates dups, and removes over-represented kmers.
     */
    private static final class KmerCleaner {
        private static final int MAX_INTERVALS = 3;
        private static final int MIN_KMER_COUNT = 3;
        private static final int MAX_KMER_COUNT = 125;
        private static final int MAX_KMERS_PER_REHASHED_PARTITION_GUESS = 600000;

        private final HopscotchHashSet<KmerAndInterval> kmerSet =
                new HopscotchHashSet<>(MAX_KMERS_PER_REHASHED_PARTITION_GUESS);

        public Iterable<KmerAndInterval> call(final Iterator<Tuple2<KmerAndInterval, Integer>> kmerCountItr ) {

            // remove kmers with extreme counts that won't help in building a local assembly
            while ( kmerCountItr.hasNext() ) {
                final Tuple2<KmerAndInterval, Integer> kmerCount = kmerCountItr.next();
                final int count = kmerCount._2;
                if ( count >= MIN_KMER_COUNT && count <= MAX_KMER_COUNT ) kmerSet.add(kmerCount._1);
            }

            // remove ubiquitous kmers that appear in too many intervals
            final int kmerSetCapacity = kmerSet.capacity();
            for ( int idx = 0; idx != kmerSetCapacity; ++idx ) {
                cleanBucket(idx);
            }

            return kmerSet;
        }

        /** clean up ubiquitous kmers that hashed to some given hash bucket */
        private void cleanBucket( final int bucketIndex ) {
            final int nKmers = SVUtils.iteratorSize(kmerSet.bucketIterator(bucketIndex));
            if ( nKmers <= MAX_INTERVALS ) return;

            // create an array of all the entries in the given bucket
            final KmerAndInterval[] kmers = new KmerAndInterval[nKmers];
            final Iterator<KmerAndInterval> itr = kmerSet.bucketIterator(bucketIndex);
            int idx = 0;
            while ( itr.hasNext() ) {
                kmers[idx++] = itr.next();
            }

            // sort the entries by kmer
            Arrays.sort(kmers);

            // remove entries sharing a kmer value that appear in more than MAX_INTERVALS intervals
            // readIdx starts a group, and testIdx is bumped until the kmer changes (ignoring interval id)
            int readIdx = 0;
            int testIdx = 1;
            while ( testIdx < nKmers ) {
                // first kmer in the group
                final SVKmer test = kmers[readIdx];

                // bump testIdx until the kmer changes (or we run out of entries)
                while ( testIdx < nKmers && test.equals(kmers[testIdx]) ) {
                    ++testIdx;
                }

                // if the number of entries with the same kmer is too big
                if ( testIdx-readIdx > MAX_INTERVALS ) {
                    // kill 'em
                    while ( readIdx != testIdx ) {
                        kmerSet.remove(kmers[readIdx++]);
                    }
                }
                readIdx = testIdx;
                testIdx += 1;
            }
        }
    }

    /**
     * Class that acts as a mapper from a stream of reads to a stream of <intervalId,read> pairs.
     * It knows which breakpoint(s) a read belongs to (if any) by kmerizing the read, and looking up each SVKmer in
     * a multi-map of SVKmers onto intervalId.
     */
    private static final class QNamesForKmersFinder implements Function<GATKRead, Iterator<QNameAndInterval>> {
        private final HopscotchHashSet<KmerAndInterval> kmerAndIntervalSet;
        private final Set<Integer> intervalIdSet = new HashSet<>();
        private final List<QNameAndInterval> qNameAndIntervalList = new ArrayList<>();
        private final Iterator<QNameAndInterval> emptyIterator = Collections.emptyIterator();

        QNamesForKmersFinder(final HopscotchHashSet<KmerAndInterval> kmerAndIntervalSet) {
            this.kmerAndIntervalSet = kmerAndIntervalSet;
        }

        public Iterator<QNameAndInterval> apply(final GATKRead read) {
            intervalIdSet.clear();
            SVKmerizer.stream(read.getBases(), SVConstants.KMER_SIZE)
                    .map( kmer -> kmer.canonical(SVConstants.KMER_SIZE) )
                    .forEach( kmer -> {
                        final Iterator<KmerAndInterval> itr = kmerAndIntervalSet.bucketIterator(kmer);
                        while ( itr.hasNext() ) {
                            final KmerAndInterval kmerAndInterval = itr.next();
                            if (kmer.equals(kmerAndInterval)) {
                                intervalIdSet.add(kmerAndInterval.getIntervalId());
                            }
                        }
                    });
            if (intervalIdSet.isEmpty()) return emptyIterator;
            qNameAndIntervalList.clear();
            final String qName = read.getName();
            intervalIdSet.stream().forEach(intervalId ->
                    qNameAndIntervalList.add(new QNameAndInterval(qName, intervalId)));
            return qNameAndIntervalList.iterator();
        }
    }

    /**
     * Find <intervalId,fastqBytes> pairs for interesting template names.
     */
    private static final class ReadsForQNamesFinder {
        private final HopscotchHashSet<QNameAndInterval> qNamesSet;
        private final int nIntervals;
        private final int nReadsPerInterval;

        @SuppressWarnings("unchecked")
        ReadsForQNamesFinder( final HopscotchHashSet<QNameAndInterval> qNamesSet, final int nIntervals ) {
            this.qNamesSet = qNamesSet;
            this.nIntervals = nIntervals;
            this.nReadsPerInterval = 2*qNamesSet.size()/nIntervals;
        }

        public Iterable<Tuple2<Integer, List<byte[]>>> call( final Iterator<GATKRead> readsItr ) {
            @SuppressWarnings({"unchecked", "rawtypes"})
            final List<byte[]>[] intervalReads = new List[nIntervals];
            int nPopulatedIntervals = 0;
            while ( readsItr.hasNext() ) {
                final GATKRead read = readsItr.next();
                final String readName = read.getName();
                final int readNameHash = readName.hashCode();
                final byte[] readNameBytes = readName.getBytes();
                final Iterator<QNameAndInterval> namesItr = qNamesSet.bucketIterator(readName);
                byte[] fastqBytes = null;
                while (namesItr.hasNext()) {
                    final QNameAndInterval nameAndInterval = namesItr.next();
                    if (nameAndInterval.hashCode() == readNameHash && nameAndInterval.sameName(readNameBytes)) {
                        if (fastqBytes == null) fastqBytes = fastqForRead(read).getBytes();
                        final int intervalId = nameAndInterval.getIntervalId();
                        if ( intervalReads[intervalId] == null ) {
                            intervalReads[intervalId] = new ArrayList<>(nReadsPerInterval);
                            nPopulatedIntervals += 1;
                        }
                        intervalReads[intervalId].add(fastqBytes);
                    }
                }
            }
            final List<Tuple2<Integer, List<byte[]>>> fastQRecords = new ArrayList<>(nPopulatedIntervals);
            if ( nPopulatedIntervals > 0 ) {
                for ( int idx = 0; idx != nIntervals; ++idx ) {
                    final List<byte[]> readList = intervalReads[idx];
                    if ( readList != null ) fastQRecords.add(new Tuple2<>(idx, readList));
                }
            }
            return fastQRecords;
        }

        private String fastqForRead( final GATKRead read ) {
            final String nameSuffix = read.isPaired() ? (read.isFirstOfPair() ? "/1" : "/2") : "";
            //final String mappedLocation = read.isUnmapped() ? "*" : read.getContig()+":"+read.getStart();
            return "@" + read.getName() + nameSuffix +
                    // "|" + mappedLocation +
                    "\n" +
                    read.getBasesString() + "\n" +
                    "+\n" +
                    ReadUtils.getBaseQualityString(read)+"\n";
        }
    }
}

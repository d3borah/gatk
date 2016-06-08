package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.utils.HopscotchHashSet;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import scala.Tuple2;

import java.io.*;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Set;

public final class FindBreakpointEvidenceSparkUnitTest extends BaseTest {
    private static final ReadMetadata.ReadGroupFragmentStatistics testStats =
            new ReadMetadata.ReadGroupFragmentStatistics(400.f, 75.f);
    private static final FindBreakpointEvidenceSpark.Interval testInterval =
            new FindBreakpointEvidenceSpark.Interval(1, 33143134, 33143478);

    private final String toolDir = getToolTestDataDir();
    private final String readsFile = toolDir+"SVBreakpointsTest.bam";
    private final String qNamesFile = toolDir+"SVBreakpointsTest.qnames";
    private final String kmersFile = toolDir+"SVBreakpointsTest.kmers";
    private final String asmQNamesFile = toolDir+"SVBreakpointsTest.assembly.qnames";
    private final String fastqFile = toolDir+"SVBreakpointsTest.assembly";

    private final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();
    private final ReadsSparkSource readsSource = new ReadsSparkSource(ctx);
    private final SAMFileHeader header = readsSource.getHeader(readsFile, null, null);
    private final JavaRDD<GATKRead> reads = readsSource.getParallelReads(readsFile, null, null, 0L);
    private final JavaRDD<GATKRead> mappedReads = reads.filter(read -> !read.isUnmapped());
    private final ReadMetadata readMetadataExpected =
            new ReadMetadata(header, Arrays.asList(testStats, testStats, testStats), testStats);
    private final Broadcast<ReadMetadata> broadcastMetadata = ctx.broadcast(readMetadataExpected);
    private final FindBreakpointEvidenceSpark.Locations locations =
        new FindBreakpointEvidenceSpark.Locations(null, null, null, null, null, null);
    private final HopscotchHashSet<FindBreakpointEvidenceSpark.QNameAndInterval> expectedQNames = loadExpectedQNames(qNamesFile);
    private final HopscotchHashSet<FindBreakpointEvidenceSpark.QNameAndInterval> expectedAssemblyQNames = loadExpectedQNames(asmQNamesFile);
    private final List<FindBreakpointEvidenceSpark.Interval> expectedIntervalList = Collections.singletonList(testInterval);

    @Test(groups = "spark")
    public void getMetadataTest() {
        final ReadMetadata readMetadataActual = FindBreakpointEvidenceSpark.getMetadata(header);
        Assert.assertEquals(readMetadataActual, readMetadataExpected);
    }

    @Test(groups = "spark")
    public void getIntervalsTest() {
        final List<FindBreakpointEvidenceSpark.Interval> actualIntervals =
                FindBreakpointEvidenceSpark.getIntervals(broadcastMetadata, header, mappedReads, locations);
        Assert.assertEquals(actualIntervals, expectedIntervalList);
    }

    @Test(groups = "spark")
    public void getQNamesTest() {
        final HopscotchHashSet<FindBreakpointEvidenceSpark.QNameAndInterval> qNamesSetActual =
                FindBreakpointEvidenceSpark.getQNames(ctx, broadcastMetadata, expectedIntervalList, mappedReads);
        Assert.assertEquals(qNamesSetActual, expectedQNames);
    }

    @Test(groups = "spark")
    public void getKmerIntervalsTest() {
        final HopscotchHashSet<SVKmer> killSet = new HopscotchHashSet<>(2);
        killSet.add(SVKmerizer.toKmer("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
        killSet.add(SVKmerizer.toKmer("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"));
        final List<SVKmer> highCountKmers = FindBreakpointEvidenceSpark.getHighCountKmers(reads, locations, null);
        Assert.assertEquals(highCountKmers.size(), 2);
        killSet.addAll(highCountKmers);
        Assert.assertEquals(killSet.size(), 2);

        final Set<FindBreakpointEvidenceSpark.KmerAndInterval> actualKmerAndIntervalSet =
                new HopscotchHashSet<>(
                        FindBreakpointEvidenceSpark.getKmerIntervals(ctx, expectedQNames, killSet, reads, locations, null));
        final Set<SVKmer> expectedKmers = SVKmer.readKmersFile(kmersFile, null);
        Assert.assertEquals(actualKmerAndIntervalSet.size(), expectedKmers.size());
        for ( final FindBreakpointEvidenceSpark.KmerAndInterval kmerAndInterval : actualKmerAndIntervalSet ) {
            Assert.assertTrue(expectedKmers.contains(new SVKmer(kmerAndInterval)));
        }
    }

    @Test(groups = "spark")
    public void getAssemblyQNamesTest() throws FileNotFoundException {
        final Set<SVKmer> expectedKmers = SVKmer.readKmersFile(kmersFile, null);
        final HopscotchHashSet<FindBreakpointEvidenceSpark.KmerAndInterval> kmerAndIntervalSet =
                new HopscotchHashSet<>(expectedKmers.size());
        for ( final SVKmer kmer : expectedKmers ) {
            kmerAndIntervalSet.add(new FindBreakpointEvidenceSpark.KmerAndInterval(kmer, 0));
        }
        final HopscotchHashSet<FindBreakpointEvidenceSpark.QNameAndInterval> actualAssemblyQNames =
                new HopscotchHashSet<>(FindBreakpointEvidenceSpark.getAssemblyQNames(ctx, kmerAndIntervalSet, reads));
        Assert.assertEquals(expectedAssemblyQNames.size(), actualAssemblyQNames.size());
    }

    @Test(groups = "spark")
    public void generateFastqsTest() {
        final String expectedFile = fastqFile;
        FindBreakpointEvidenceSpark.generateFastqs(ctx, expectedAssemblyQNames, 1, reads,
                intervalAndFastqBytes -> compareFastqs(intervalAndFastqBytes, expectedFile));
    }

    private static void compareFastqs(
            final Tuple2<Integer, List<byte[]>> intervalAndFastqBytes,
            final String fastqFile ) throws IOException {
        final byte[][] fastqs = FindBreakpointEvidenceSpark.sortFastqRecs(intervalAndFastqBytes._2);
        final byte[] concatenatedFastqs = new byte[Arrays.stream(fastqs).mapToInt(fastq -> fastq.length).sum()];
        int idx = 0;
        for ( final byte[] fastq : fastqs ) {
            System.arraycopy(fastq, 0, concatenatedFastqs, idx, fastq.length);
            idx += fastq.length;
        }
        final ByteArrayInputStream actualStream = new ByteArrayInputStream(concatenatedFastqs);
        try( InputStream expectedStream = new BufferedInputStream(new FileInputStream(fastqFile)) ) {
            int val;
            while ( (val = expectedStream.read()) != -1 )
                Assert.assertEquals(actualStream.read(), val);
            Assert.assertEquals(actualStream.read(), -1);
        }
    }

    private HopscotchHashSet<FindBreakpointEvidenceSpark.QNameAndInterval> loadExpectedQNames( final String fileName ) {
        final HopscotchHashSet<FindBreakpointEvidenceSpark.QNameAndInterval> expectedSet = new HopscotchHashSet<>(300);
        try( BufferedReader rdr = new BufferedReader(new FileReader(fileName)) ) {
            for ( String line = rdr.readLine(); line != null; line = rdr.readLine() ) {
                expectedSet.add(new FindBreakpointEvidenceSpark.QNameAndInterval(line, 0));
            }
        } catch ( final IOException ioe ) {
            throw new GATKException("Unable to read " + fileName);
        }
        return expectedSet;
    }
}

package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.metrics.Header;
import org.broadinstitute.hellbender.exceptions.UserException;

import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import htsjdk.samtools.util.IterableAdapter;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.AuthHolder;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.metrics.*;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.R.RScriptExecutorException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Serializable;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * Worker class to collect insert size metrics, add metrics to file, and provides
 * accessors to stats of groups of different level.
 */
public final class InsertSizeMetricsCollectorSparkTemp
        implements MetricsCollectorSpark<InsertSizeMetricsArgs>,
                Serializable
{
    private static final long serialVersionUID = 1L;

    // TODO: This should really go into the constructor of the XMetricsCollector along with the SAMGroupRecords
    MetricsFile<InsertSizeMetrics, Integer> metricsFile = null;

    private InsertSizeMetricsArgs inputArgs = null;
    private InsertSizeMetricsCollectorTemp collector = null;

    // path to Picard R script for producing histograms in PDF files.
    private static final String R_SCRIPT = "insertSizeHistogram.R";

    List<SAMReadGroupRecord> rgRecords = null;

    //TODO: rationalize what is passed into the constructors
    // initializeCollector; initializeCollector should only take input args defined
    // for the collector (ie.e. command line args); things like SamFileHeader, default headers
    // metrics file, etc should go in the cosntructor ??
    public InsertSizeMetricsCollectorSparkTemp(final List<SAMReadGroupRecord> samRgRecords) {
        this.rgRecords = samRgRecords;
    }

    /**
     * Initialize the collector with input arguments;
     */
    @Override
    public void initializeCollector(
            final InsertSizeMetricsArgs inputArgs,
            final List<Header> defaultHeaders)
    {
        collector = new InsertSizeMetricsCollectorTemp(inputArgs, rgRecords);
        metricsFile = new MetricsFile<InsertSizeMetrics, Integer>();
        defaultHeaders.stream().forEach(h -> metricsFile.addHeader(h));
        this.inputArgs = inputArgs;
    }

    /**
     * Return a read filter to be used for this collector.
     *
     * @param samFileHeader
     * @return ReadFilter
     */
    @Override
    public ReadFilter getCollectorReadFilter(final SAMFileHeader samFileHeader) {

        final InsertSizeMetricsReadFilter svFilter = new InsertSizeMetricsReadFilter(
                inputArgs.useEnd,
                inputArgs.filterNonProperlyPairedReads,
                !inputArgs.useDuplicateReads,
                !inputArgs.useSecondaryAlignments,
                !inputArgs.useSupplementaryAlignments,
                inputArgs.MQPassingThreshold);

        return new WellformedReadFilter(samFileHeader).and(svFilter);
    }

    @Override
    public void collectMetrics(
            final JavaRDD<GATKRead> filteredReads,
            final SAMFileHeader samHeader,
            final String inputBaseName,
            final AuthHolder authHolder)
    {
        if (filteredReads.isEmpty()) {
            throw new GATKException("No valid reads found in input file.");
        }

        InsertSizeMetricsCollectorTemp summary =
                filteredReads//.repartition(2)
                        .mapPartitions((Iterator<GATKRead> it) -> {
                            // we could allocate a collector here instead of having this class "is-a" collector ??
                            it.forEachRemaining(r -> collector.acceptRecord(r.convertToSAMRecord(samHeader), null));
                            return Arrays.asList(collector);
                        }).reduce(collector::combine);

//        InsertSizeMetricsCollectorTemp summary =
//                filteredReads.repartition(2)
//                        .aggregate(
//                                collector,
//                                (acc, read) -> { acc.acceptRecord(read.convertToSAMRecord(samHeader), null); return acc;},
//                                collector::combine
//                        );

        summary.finish();
        summary.addAllLevelsToFile(metricsFile);
        MetricsUtils.saveMetrics(metricsFile, inputArgs.output.getAbsolutePath(), authHolder);
        writeHistogramPDF(inputBaseName);
    }

    /**
     * Customized serializable reads filter, based on cmd line arguments provided
     */
    private static final class InsertSizeMetricsReadFilter implements ReadFilter{
        private static final long serialVersionUID = 1L;

        private final ReadFilter combinedReadFilter;

        public InsertSizeMetricsReadFilter(final InsertSizeMetricsArgs.EndToUse whichEnd,
                                           final boolean filterNonProperlyPairedReads,
                                           final boolean filterDuplicatedReads,
                                           final boolean filterSecondaryAlignments,
                                           final boolean filterSupplementaryAlignments,
                                           final int     MQThreshold){

            final InsertSizeMetricsArgs.EndToUse endVal = whichEnd;

            ReadFilter tempFilter = ReadFilterLibrary.MAPPED;
            tempFilter = tempFilter.and(GATKRead::isPaired);
            tempFilter = tempFilter.and(read -> 0!=read.getFragmentLength());
            tempFilter = tempFilter.and(read -> endVal == (read.isFirstOfPair() ?
                    InsertSizeMetricsArgs.EndToUse.FIRST :
                    InsertSizeMetricsArgs.EndToUse.SECOND));

            if(filterNonProperlyPairedReads)  { tempFilter = tempFilter.and(GATKRead::isProperlyPaired); }
            if(filterDuplicatedReads)         { tempFilter = tempFilter.and(read -> !read.isDuplicate()); }
            if(filterSecondaryAlignments)     { tempFilter = tempFilter.and(read -> !read.isSecondaryAlignment()); }
            if(filterSupplementaryAlignments) { tempFilter = tempFilter.and(read -> !read.isSupplementaryAlignment()); }

            if(0!=MQThreshold)  { tempFilter = tempFilter.and(read -> read.getMappingQuality() >= MQThreshold);}

            combinedReadFilter = tempFilter;
        }

        @Override
        public boolean test(final GATKRead read){
            return combinedReadFilter.test(read);
        }
    }

    /**
     * Calls R script to plot histogram(s) in PDF.
     * @throws RScriptExecutorException
     */
    @VisibleForTesting
    void writeHistogramPDF(String inputName) throws RScriptExecutorException {

        final File histogramPlotPDF = inputArgs.histogramPlotFile;
        IOUtil.assertFileIsWritable(histogramPlotPDF);

        final RScriptExecutor executor = new RScriptExecutor();
        executor.addScript(new Resource(R_SCRIPT, org.broadinstitute.hellbender.tools.spark.pipelines.metrics.InsertSizeMetricsCollectorSpark.class));
        executor.addArgs(inputArgs.output,                                // text-based metrics file
                histogramPlotPDF.getAbsolutePath(),    // PDF graphics file
                inputName);                  // input bam file
        executor.exec();
    }

}
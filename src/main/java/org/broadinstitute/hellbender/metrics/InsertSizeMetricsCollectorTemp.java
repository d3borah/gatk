package org.broadinstitute.hellbender.metrics;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMReadGroupRecord;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamPairUtil;
import htsjdk.samtools.metrics.Header;
import htsjdk.samtools.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.util.IOUtil;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.AuthHolder;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.filters.WellformedReadFilter;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.pipelines.metrics.InsertSizeMetricsCollectorSparkArgs;
import org.broadinstitute.hellbender.tools.spark.pipelines.metrics.MetricsCollectorSpark;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.R.RScriptExecutorException;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import scala.Tuple2;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

public class InsertSizeMetricsCollectorTemp
        extends MultiLevelCollectorTemp<InsertSizeMetrics, Integer, InsertSizeMetricsCollectorArgs, PerUnitInsertSizeMetricsCollector>
        implements Serializable
{
    private static final long serialVersionUID = 1L;

    // path to Picard R script for producing histograms in PDF files.
    private static final String R_SCRIPT = "insertSizeHistogram.R";

    InsertSizeMetricsArgs inputArgs = null;
    MetricsFile<InsertSizeMetrics, Integer> metricsFile = null;


    public InsertSizeMetricsCollectorTemp(
            final InsertSizeMetricsArgs inputArgs,
            final List<SAMReadGroupRecord> samRgRecords)
    {
        super();
        this.inputArgs = inputArgs;
        setup(inputArgs.metricAccumulationLevel, samRgRecords);
    }

    // We will pass insertSize and PairOrientation with the DefaultPerRecordCollectorArgs passed to the record collectors
    // This method is called once Per samRecord
    @Override
    protected InsertSizeMetricsCollectorArgs makeArg(SAMRecord samRecord, ReferenceSequence refSeq) {
        final int insertSize = Math.abs(samRecord.getInferredInsertSize());
        final SamPairUtil.PairOrientation orientation = SamPairUtil.getPairOrientation(samRecord);

        return new InsertSizeMetricsCollectorArgs(insertSize, orientation);
    }

    /** Make an InsertSizeCollector with the given arguments */
    @Override
    protected PerUnitInsertSizeMetricsCollector makeChildCollector(
            final String sample,
            final String library,
            final String readGroup) {
        return new PerUnitInsertSizeMetricsCollector(
                sample,
                library,
                readGroup,
                inputArgs.minimumPct,
                inputArgs.maxMADTolerance,
                inputArgs.histogramWidth);
    }

    @Override
    public void acceptRecord(final SAMRecord record, final ReferenceSequence refSeq) {
        //TODO: we shouldn't mix filtering with accept/(process) since for Spark
        //they're separate

        if (!record.getReadPairedFlag() ||
                record.getReadUnmappedFlag() ||
                record.getMateUnmappedFlag() ||
                record.getFirstOfPairFlag() ||
                record.isSecondaryOrSupplementary() ||
                record.getDuplicateReadFlag() ||
                record.getInferredInsertSize() == 0) {
            return;
        }

        super.acceptRecord(record, refSeq);
    }

    public InsertSizeMetricsCollectorTemp combine(InsertSizeMetricsCollectorTemp source, InsertSizeMetricsCollectorTemp target) {

        // TODO: this mutates source!
        source.allReadCollector = combineUnit(source.allReadCollector, target.allReadCollector);
        super.combine(source, target);
        return source;
    }

    public PerUnitInsertSizeMetricsCollector combineUnit(
            PerUnitInsertSizeMetricsCollector sourceCollector,
            PerUnitInsertSizeMetricsCollector targetCollector) {
        return sourceCollector.combine(targetCollector);
    }

    /**
     * Calls R script to plot histogram(s) in PDF.
     *
     * @throws RScriptExecutorException
     */
    @VisibleForTesting
    void writeHistogramPDF(String inputName) throws RScriptExecutorException {

        //final File histogramPlotPDF = new File(inputArgs.histogramPlotFile);
        IOUtil.assertFileIsWritable(inputArgs.histogramPlotFile);

        final RScriptExecutor executor = new RScriptExecutor();
        executor.addScript(new Resource(R_SCRIPT, org.broadinstitute.hellbender.tools.spark.pipelines.metrics.InsertSizeMetricsCollectorSpark.class));
        executor.addArgs(inputArgs.output,                                // text-based metrics file
                inputArgs.histogramPlotFile.getAbsolutePath(),    // PDF graphics file
                inputName);                  // input bam file
        executor.exec();
    }

}

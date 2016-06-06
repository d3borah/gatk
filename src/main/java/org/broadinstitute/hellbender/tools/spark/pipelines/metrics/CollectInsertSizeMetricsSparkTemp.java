package org.broadinstitute.hellbender.tools.spark.pipelines.metrics;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.metrics.Header;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.AuthHolder;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.metrics.InsertSizeMetricsArgs;
import org.broadinstitute.hellbender.metrics.InsertSizeMetricsCollectorTemp;
import org.broadinstitute.hellbender.utils.read.GATKRead;

import java.io.ByteArrayOutputStream;
import java.io.ObjectOutputStream;
import java.util.List;

/**
 * Spark tool for collecting insert size metrics.
 */
@CommandLineProgramProperties(
        summary        = "Program to collect insert size distribution information in SAM/BAM/CRAM file(s)",
        oneLineSummary = "Collect Insert Size Distribution on Spark",
        programGroup   = SparkProgramGroup.class)
public final class CollectInsertSizeMetricsSparkTemp
        extends MetricsCollectorToolSpark<InsertSizeMetricsArgs> {

    private static final long serialVersionUID = 1L;

    @ArgumentCollection
    InsertSizeMetricsArgs insertSizeArgs = new InsertSizeMetricsArgs();

    InsertSizeMetricsCollectorSparkTemp insertSizeCollector = null;

    @Override
    public InsertSizeMetricsArgs gatherInputArguments() {
        return insertSizeArgs;
    }

    public void initializeCollector(
            final InsertSizeMetricsArgs inputArgs,
            final List<Header> defaultHeaders)
    {
        insertSizeCollector = new InsertSizeMetricsCollectorSparkTemp(
                getHeaderForReads().getReadGroups());
        insertSizeCollector.initializeCollector(inputArgs, defaultHeaders);
    }

    /**
     * Expose the read filter required for this collector
     */
    @Override
    public ReadFilter getCollectorReadFilter(final SAMFileHeader samHeader) {
        //Delegate actual collection to InsertSizeMetricCollector
        return insertSizeCollector.getCollectorReadFilter(samHeader);
    }

    @Override
    public void collectMetrics(
            final JavaRDD<GATKRead> filteredReads,
            final  SAMFileHeader samHeader,
            final String inputBaseName,
            final AuthHolder authHolder)
    {
        insertSizeCollector.collectMetrics(
                filteredReads,
                samHeader,
                inputBaseName,
                authHolder
        );
    }
}

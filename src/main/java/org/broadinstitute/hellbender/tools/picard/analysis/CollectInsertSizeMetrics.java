package org.broadinstitute.hellbender.tools.picard.analysis;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.metrics.MetricsFile;
//import org.broadinstitute.hellbender.metrics.MetricsFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.IOUtil;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.QCProgramGroup;
import org.broadinstitute.hellbender.metrics.*;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.io.Resource;

import java.io.File;
import java.util.EnumSet;
import java.util.Set;

/**
 * Command line program to read non-duplicate insert sizes, create a Histogram
 * and report distribution statistics.
 *
 * @author Doug Voet (dvoet at broadinstitute dot org)
 */
@CommandLineProgramProperties(
        summary = "Reads a SAM/BAM/CRAM file and writes a file containing metrics about " +
                "the statistical distribution of insert size (excluding duplicates) " +
                "and generates a Histogram plot.",
        oneLineSummary = "Produces metrics for insert size distribution for a SAM/BAM/CRAM file",
        programGroup = QCProgramGroup.class
)
public final class CollectInsertSizeMetrics extends SinglePassSamProgram {
    private static final Logger log = LogManager.getLogger(CollectInsertSizeMetrics.class);

    private static final String R_SCRIPT = "insertSizeHistogram.R";

    @ArgumentCollection
    public InsertSizeMetricsArgs inputArgs = new InsertSizeMetricsArgs();

    // Calculates InsertSizeMetrics for all METRIC_ACCUMULATION_LEVELs provided
    private InsertSizeMetricsCollectorTemp multiCollector;

    /**
     * Put any custom command-line validation in an override of this method.
     * clp is initialized at this point and can be used to print usage and access argv.
     * Any options set by command-line parser can be validated.
     *
     * @return null if command line is valid.  If command line is invalid, returns an array of error message
     *         to be written to the appropriate place.
     */
    @Override
    protected String[] customCommandLineValidation() {
         if (inputArgs.minimumPct < 0 || inputArgs.minimumPct > 0.5) {
             return new String[]{"MINIMUM_PCT was set to " + inputArgs.minimumPct +
                     ". It must be between 0 and 0.5 so all data categories don't get discarded."};
         }

         return super.customCommandLineValidation();
    }

    @Override
    protected boolean usesNoRefReads() { return false; }

    @Override
    protected void setup(final SAMFileHeader header, final File samFile) {
        IOUtil.assertFileIsWritable(inputArgs.output);
        IOUtil.assertFileIsWritable(inputArgs.histogramPlotFile);

        //Delegate actual collection to InsertSizeMetricCollector
        multiCollector = new InsertSizeMetricsCollectorTemp(
                inputArgs,
                header.getReadGroups());
    }

    @Override
    protected void acceptRead(final SAMRecord record, final ReferenceSequence ref) {
        multiCollector.acceptRecord(record, ref);
    }

    @Override
    protected void finish() {
        multiCollector.finish();

        //TODO hack
//        final MetricsFile<InsertSizeMetrics, Integer> file = getMetricsFile();
        //TODO hack
        final MetricsFile<InsertSizeMetrics, Integer> file = new MetricsFile<InsertSizeMetrics, Integer>();
        multiCollector.addAllLevelsToFile(file);

        if(file.getNumHistograms() == 0) {
            //can happen if user sets MINIMUM_PCT = 0.5, etc.
            log.warn("All data categories were discarded because they contained < " + inputArgs.minimumPct +
                     " of the total aligned paired data.");
            final PerUnitInsertSizeMetricsCollector allReadsCollector = multiCollector.getAllReadsCollector();
            log.warn("Total mapped pairs in all categories: " + (allReadsCollector == null ? allReadsCollector : allReadsCollector.getTotalInserts()));
        }
        else  {
            file.write(inputArgs.output);
            if (inputArgs.producePlot){
                final RScriptExecutor executor = new RScriptExecutor();
                executor.addScript(new Resource(R_SCRIPT, CollectInsertSizeMetrics.class));
                executor.addArgs(inputArgs.output.getAbsolutePath(), inputArgs.histogramPlotFile.getAbsolutePath(), INPUT.getName());
                if (inputArgs.histogramWidth != null) executor.addArgs(String.valueOf(inputArgs.histogramWidth));
                executor.exec();
            }
        }
    }
}

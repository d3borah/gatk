package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.base.Stopwatch;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import scala.Tuple2;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.concurrent.TimeUnit;

@CommandLineProgramProperties(summary="Align assembled contigs to the reference and call breakpoints from them.",
        oneLineSummary="Align contigs and call breakpoints",
        programGroup = SparkProgramGroup.class)
public class AlignContigsAndCallBreakpointsSpark extends GATKSparkTool {

    @Argument(doc = "file for assembled breakpoint output", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String output;

    @Argument(doc = "Input file of assembled contigs", shortName = StandardArgumentDefinitions.INPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME, optional = false)
    private String input;

    private static final Logger log = LogManager.getLogger(AlignContigsAndCallBreakpointsSpark.class);

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {

        final JavaPairRDD<String, String> contigCollectionByBreakpointId = ctx.textFile(input).mapToPair(this::splitAssemblyLine);

        final JavaPairRDD<String, RunSGAViaProcessBuilderOnSpark.ContigsCollection> breakpointIdsToContigsCollection =
                contigCollectionByBreakpointId.mapValues(v -> new RunSGAViaProcessBuilderOnSpark.ContigsCollection(Arrays.asList(v.split("|"))));

        final JavaPairRDD<String, ContigAligner.AssembledBreakpoint> assembledBreakpoints = breakpointIdsToContigsCollection.mapPartitionsToPair(iter -> () -> {
            try {
                final String ref = referenceArguments.getReferenceFileName();
                final File localRef = BucketUtils.isHadoopUrl(ref) ? localizeReferenceAndBwaIndexFiles(ref) : new File(ref);
                final ContigAligner contigAligner = new ContigAligner();
                final List<Tuple2<String, ContigAligner.AssembledBreakpoint>> results = new ArrayList<>();
                iter.forEachRemaining(cc -> {
                    final List<ContigAligner.AssembledBreakpoint> breakpointAssemblies = contigAligner.alignContigs(cc._2, localRef);
                    breakpointAssemblies.forEach(b -> {
                        results.add(new Tuple2<>(cc._1, b));
                    });
                });
                return results.iterator();
            } catch (IOException e) {
                throw new GATKException("Cannot run BWA-MEM", e);
            }

        });

        assembledBreakpoints.saveAsTextFile(output);

    }

    private Tuple2<String, String> splitAssemblyLine(final String assemblyLine) {
        final String[] split = assemblyLine.split("\t");
        return new Tuple2<>(split[0], split[1]);
    }

    private File localizeReferenceAndBwaIndexFiles(String referenceFilename) throws IOException {
        Stopwatch downloadRefStopwatch = Stopwatch.createStarted();
        File localRef = File.createTempFile("ref", ".fa");
        if (!localRef.delete()) {
            throw new IOException("Cannot delete temporary file for reference: " + localRef);
        }
        Path localPath = localRef.toPath();
        Path remotePath = IOUtils.getPath(referenceFilename);
        Files.copy(remotePath, localPath);
        for (String extension : new String[] { ".amb", ".ann", ".bwt", ".pac", ".sa" }) {
            Files.copy(remotePath.resolveSibling(remotePath.getFileName() + extension),
                    localPath.resolveSibling(localPath.getFileName() + extension));
        }
        downloadRefStopwatch.stop();
        log.info("Time to download reference: " + downloadRefStopwatch.elapsed(TimeUnit.SECONDS) + "s");
        return localRef;
    }

}
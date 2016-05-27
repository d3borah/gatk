package org.broadinstitute.hellbender.tools.spark.sv;

import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariationSparkProgramGroup;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import scala.Tuple2;

import java.io.BufferedOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.HashSet;
import java.util.List;

/**
 * Created by shuang on 5/27/16.
 */
@CommandLineProgramProperties(
        summary        = "Temporary program to remove duplicated fastq records in a fastq file. ",
        oneLineSummary = "Perform SGA-based local assembly on fasta files on Spark.",
        programGroup   = StructuralVariationSparkProgramGroup.class)
public final class DeduplicateFASTQFiles extends GATKSparkTool {
    private static final long serialVersionUID = 1L;

    @Argument(doc       = "An URI to the directory where all interleaved FASTQ files for putative breakpoints are located.",
            shortName = "in",
            fullName  = "inDir",
            optional  = false)
    public String pathToAllInterleavedFASTQFiles = null;

    @Override
    public void runTool(final JavaSparkContext ctx){

        // first load RDD of pair that has path to FASTQ file path as its first and FASTQ file contents as its second
        final JavaPairRDD<String, String> fastqContentsForEachBreakpoint = ctx.wholeTextFiles(pathToAllInterleavedFASTQFiles);

        final JavaPairRDD<String, String> result = fastqContentsForEachBreakpoint.mapToPair(DeduplicateFASTQFiles::dedup);

        result.foreach(DeduplicateFASTQFiles::writeFastq);
    }

    private static final Tuple2<String, String> dedup(final Tuple2<String, String> duplicatedFASTQContents){

        final String contentsOneLiner = duplicatedFASTQContents._2();
        final String[] contents = contentsOneLiner.split("\n");

        final StringBuilder dedupedContentsBuilder = new StringBuilder();

        final HashSet<String> uniqueTemplateNames = new HashSet<>();
        final int sz = contents.length;
        for(int i=0; i<sz; i+=4){
            if(!uniqueTemplateNames.contains(contents[i])){
                dedupedContentsBuilder.append(contents[i]);
                dedupedContentsBuilder.append(contents[i+1]);
                dedupedContentsBuilder.append(contents[i+2]);
                dedupedContentsBuilder.append(contents[i+3]);
                uniqueTemplateNames.add(contents[i]);
            }
        }

        final String originalSpace = duplicatedFASTQContents._1();
        final String targetSpace = originalSpace.replace(".sav", ".sav.dedup");

        return new Tuple2<>(targetSpace, dedupedContentsBuilder.toString());
    }

    private static final void writeFastq( final Tuple2<String, String> pair) {
        final String fileName = pair._1();
        try ( final OutputStream writer = new BufferedOutputStream(BucketUtils.createNonGCSFile(fileName)) ) {
            writer.write(pair._2().getBytes());
        } catch ( final IOException ioe ) {
            throw new GATKException("Can't write "+fileName, ioe);
        }
    }
}

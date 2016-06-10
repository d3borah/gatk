package org.broadinstitute.hellbender.tools.spark.sv;

import com.github.lindenb.jbwa.jni.AlnRgn;
import com.github.lindenb.jbwa.jni.BwaIndex;
import com.github.lindenb.jbwa.jni.BwaMem;
import com.github.lindenb.jbwa.jni.ShortRead;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bwa.BWANativeLibrary;
import scala.Tuple2;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class ContigAligner {

    public ContigAligner() {
        BWANativeLibrary.load();
    }

    public List<AssembledBreakpoint> alignContigs (RunSGAViaProcessBuilderOnSpark.ContigsCollection contigsCollection, File referenceFile) {
        List<AssembledBreakpoint> assembledBreakpoints = new ArrayList<>();

        try {
            final BwaIndex index;

            index = new BwaIndex(referenceFile);
            final BwaMem bwaMem = new BwaMem(index);
            try {
                for(final Tuple2<RunSGAViaProcessBuilderOnSpark.ContigsCollection.ContigID, RunSGAViaProcessBuilderOnSpark.ContigsCollection.ContigSequence> contigInfo : contigsCollection.getContents()) {
                    final String contigId = contigInfo._1.toString();
                    final byte[] sequence = contigInfo._2.toString().getBytes();
                    final AlnRgn[] alnRgns = bwaAlignSequence(bwaMem, contigId, sequence);

                    // todo: parse cigar for internal indels
                    if (alnRgns.length > 1) {
                        // todo: order these by location in the query sequence? Does BWA do this?
                        // todo: how does BWA order secondary alignments? Should I assume they are all at the end?
                        for (int i = 0; i < alnRgns.length; i++) {
                            if (alnRgns[i].getSecondary() != -1) {
                                break;
                            }
                            if (i > 0) {
                                final AssembledBreakpoint assembledBreakpoint = new AssembledBreakpoint();
                                assembledBreakpoint.contigId = contigId;

                                final AlnRgn alnRgn1 = alnRgns[i-1];
                                final Cigar cigar1 = TextCigarCodec.decode(alnRgn1.getCigar());
                                assembledBreakpoint.region1 = new AlignmentRegion();
                                final int pos1 = (int) alnRgn1.getPos() + 1;
                                assembledBreakpoint.region1.interval = new SimpleInterval(alnRgn1.getChrom(), pos1, pos1 + cigar1.getReferenceLength());
                                assembledBreakpoint.region1.forwardStrand = ('+' == alnRgn1.getStrand());
                                assembledBreakpoint.region1.mqual = alnRgn1.getMQual();

                                final AlnRgn alnRgn2 = alnRgns[i];
                                final Cigar cigar2 = TextCigarCodec.decode(alnRgn2.getCigar());
                                assembledBreakpoint.region2 = new AlignmentRegion();
                                final int pos2 = (int) alnRgn2.getPos() + 1;
                                assembledBreakpoint.region2.interval = new SimpleInterval(alnRgn2.getChrom(), pos2, pos2 + cigar2.getReferenceLength());
                                assembledBreakpoint.region2.forwardStrand = ('+' == alnRgn2.getStrand());
                                assembledBreakpoint.region2.mqual = alnRgn2.getMQual();

                                //todo: add inserted sequence and microhomology detection
                                assembledBreakpoint.homology = "";
                                assembledBreakpoint.insertedSequence = "";

                                assembledBreakpoints.add(assembledBreakpoint);
                            }
                        }
                    }
                }
            } finally {
                bwaMem.dispose();
            }
        } catch (IOException e) {
            throw new GATKException("could not execute BWA");
        }

        return assembledBreakpoints;
    }

    private AlnRgn[] bwaAlignSequence(final BwaMem bwaMem, final String contigId, final byte[] sequence) throws IOException {
        final ShortRead contigShortRead = new ShortRead(contigId, sequence, qualSequence(sequence.length));
        AlnRgn[] alignmentRegions = bwaMem.align(contigShortRead);
        return alignmentRegions;
    }

    private byte[] qualSequence(final int length) {
        final byte[] quals = new byte[length];
        for (int i = 0; i < quals.length; i++) {
            quals[i] = 'A';
        }
        return quals;
    }

    static class AssembledBreakpoint {
        String contigId;
        AlignmentRegion region1;
        AlignmentRegion region2;
        String insertedSequence;
        String homology;
    }

    static class AlignmentRegion {
        boolean forwardStrand;
        SimpleInterval interval;
        int mqual;
    }
}

package org.broadinstitute.hellbender.tools.spark.sv;

import com.github.lindenb.jbwa.jni.AlnRgn;
import com.github.lindenb.jbwa.jni.BwaIndex;
import com.github.lindenb.jbwa.jni.BwaMem;
import com.github.lindenb.jbwa.jni.ShortRead;
import com.sun.xml.internal.ws.util.Pool;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.bwa.BWANativeLibrary;
import scala.Tuple2;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;

public class ContigAligner {

    public ContigAligner() {
        BWANativeLibrary.load();
    }

    public List<AssembledBreakpoint> alignContigs (final RunSGAViaProcessBuilderOnSpark.ContigsCollection contigsCollection, final File referenceFile) {
        final List<AssembledBreakpoint> assembledBreakpoints = new ArrayList<>();

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

                        final List<AlignmentRegion> alignmentRegionList = Arrays.stream(alnRgns)
                                .map(AlignmentRegion::new)
                                .sorted(Comparator.comparing(a -> a.startInContig))
                                .collect(Collectors.toList());

                        for (int i = 0; i < alignmentRegionList.size(); i++) {
                            if (alnRgns[i].getSecondary() != -1) {
                                break;
                            }
                            if (i > 0) {
                                final AssembledBreakpoint assembledBreakpoint = new AssembledBreakpoint();
                                assembledBreakpoint.contigId = contigId;

                                assembledBreakpoint.region1 = alignmentRegionList.get(i-1);
                                assembledBreakpoint.region2 = alignmentRegionList.get(i);

                                assembledBreakpoint.homology = "";
                                if (assembledBreakpoint.region1.endInContig >= assembledBreakpoint.region2.startInContig) {
                                    assembledBreakpoint.homology = getSequence(contigsCollection, contigId, assembledBreakpoint.region2.startInContig, assembledBreakpoint.region1.endInContig);
                                }
                                //todo: add inserted sequence
                                assembledBreakpoint.insertedSequence = "";

                                assembledBreakpoints.add(assembledBreakpoint);
                            }
                        }
                    }
                }
            } finally {
                bwaMem.dispose();
            }
        } catch (final IOException e) {
            throw new GATKException("could not execute BWA");
        }

        return assembledBreakpoints;
    }

    private String getSequence(final RunSGAViaProcessBuilderOnSpark.ContigsCollection contigsCollection, final String contigId, final int startInContig, final int endInContig) {
        for (final Tuple2<RunSGAViaProcessBuilderOnSpark.ContigsCollection.ContigID, RunSGAViaProcessBuilderOnSpark.ContigsCollection.ContigSequence> contig : contigsCollection.getContents()) {
            if (contig._1.equals(new RunSGAViaProcessBuilderOnSpark.ContigsCollection.ContigID(contigId))) {
                final String sequence = contig._2.toString();
                return sequence.substring(startInContig, endInContig);
            }
        }
        throw new GATKException("Could not find contig " + contigId);
    }

    private AlnRgn[] bwaAlignSequence(final BwaMem bwaMem, final String contigId, final byte[] sequence) throws IOException {
        final ShortRead contigShortRead = new ShortRead(contigId, sequence, qualSequence(sequence.length));
        return bwaMem.align(contigShortRead);
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

        final Cigar cigar;
        final boolean forwardStrand;
        final SimpleInterval interval;
        final int mqual;
        final int startInContig;
        final int endInContig;
        final int contigLength;


        public AlignmentRegion(final AlnRgn alnRgn) {
            this.forwardStrand = alnRgn.getStrand() == '+';
            this.cigar = TextCigarCodec.decode(alnRgn.getCigar());
            this.interval = new SimpleInterval(alnRgn.getChrom(), (int) alnRgn.getPos() + 1, (int) (alnRgn.getPos() + 1 + cigar.getReferenceLength()));
            this.mqual = alnRgn.getMQual();
            this.contigLength = cigar.getReadLength();
            this.startInContig = startOfAlignmentInContig();
            this.endInContig = endOfAlignmentInContig();
        }

        private int startOfAlignmentInContig() {
            return getClippedBases(forwardStrand);
        }

        private int endOfAlignmentInContig() {
            return contigLength - getClippedBases(!forwardStrand);
        }

        private int getClippedBases(final boolean forwardStrand) {
            int posInContig = 0;
            if (forwardStrand) {
                int i = 0;
                while (cigar.getCigarElement(i).getOperator() == CigarOperator.H || cigar.getCigarElement(i).getOperator() == CigarOperator.S) {
                    posInContig += cigar.getCigarElement(i).getLength();
                    i = i + 1;
                }
            } else {
                int i = cigar.getCigarElements().size() - 1;
                while (cigar.getCigarElement(i).getOperator() == CigarOperator.H || cigar.getCigarElement(i).getOperator() == CigarOperator.S) {
                    posInContig += cigar.getCigarElement(i).getLength();
                    i = i - 1;
                }
            }
            return posInContig;
        }

    }
}

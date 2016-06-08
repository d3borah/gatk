package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.broadinstitute.hellbender.tools.walkers.genotyper.GLVectorSizeCache;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * Handles addressing general ploidy genotype likelihoods. We store likelihoods as a map
 * of form int[] -> double (to be more precise, IntArrayWrapper -> Double).
 * For a given ploidy (chromosome count) and number of alleles, we need a form to iterate deterministically
 * across all possible allele conformations, that is, all possible ways in which N whole numbers will sum to P,
 * where N is number of alleles and P is number of chromosomes.
 * For example, with P=2,N=2, ploidy = 2 iterator will produce
 * alleleCounts = [2 0 ] [1 1] [ 0 2]
 *
 *
 */
public final class GenotpyeIterator {
    private int[] alleleCounts;
    private final int ploidy;
    private final int numAlleles;
    private boolean hasNext = true;
    private int linearIndex = 0;

    /**
     * Shortcut constructor for common use case: iterator will produce
     * all vectors of length numAlleles whose sum = ploidy
     * @param numAlleles              Number of alleles
     * @param ploidy          Ploidy
     */
    public GenotpyeIterator(final int numAlleles, final int ploidy) {
        Utils.validateArg(numAlleles > 0, "must have at least one allele");
        Utils.validateArg(ploidy >= 0, "negative ploidy forbidden");
        this.numAlleles = numAlleles;
        this.ploidy = ploidy;
        alleleCounts = new int[numAlleles];
        alleleCounts[0] = ploidy;
    }

    /**
     * To get the next genotype:
     *
     * if current state has no ref alleles, (1) find the first allele with non-zero count
     * and give all of it counts to ref (if this allele is ref, do nothing); (2) increment the allele after this one,
     * decrementing the ref,For example [0 0 2 1] -> [2 0 0 1] (step 1) -> [1 0 0 2] (step 2)
     *
     * If the last alt allele has all the counts eg [0 0 0 3] we stop and set {@code hasNext} to {@code false}.
     *
     * Thus we have, eg for 3 alleles and ploidy = 2, the ordering
     * [2 0 0] -> [1 1 0] -> [0 2 0] -> [1 0 1] -> [0 1 1] -> [0 0 2]
     */
    public void next() {
        if (ploidy == 0) {
            hasNext = false;
            return;
        }
        final int firstNonEmptyAllele = IntStream.range(0, numAlleles).filter(a -> alleleCounts[a] > 0).findFirst().getAsInt();
        if (firstNonEmptyAllele == numAlleles - 1) {
            hasNext = false;
            return;
        }
        final int firstNonEmptyAlleleCount = alleleCounts[firstNonEmptyAllele];
        alleleCounts[0] += firstNonEmptyAlleleCount;
        alleleCounts[firstNonEmptyAllele]-= firstNonEmptyAlleleCount;
        alleleCounts[firstNonEmptyAllele + 1]++;
        alleleCounts[0]--;
        linearIndex++;
    }

    public void reset() {
        Arrays.fill(alleleCounts, 0);
        alleleCounts[0] = ploidy;
        linearIndex = 0;
        hasNext = true;
    }

    public int[] getAlleleCounts() { return alleleCounts; }
    public int getLinearIndex() { return linearIndex; }
    public boolean hasNext() { return hasNext; }

    //TODO: maybe this class is the wrong place.  It kind of makes sense because this is the inverse of getAlleleCounts
    //TODO: in a sense.
    public static int getLinearIndex(final int[] alleleCounts, final int numAlleles, final int ploidy) {
        if (ploidy <= 0) {
            return 0;
        }

        int numGenotypesWithLowerIndex = 0;
        int remainingPloidy = ploidy;
        for (int allele = numAlleles-1; allele >=1; allele--) {
            final int countAtThisAllele = alleleCounts[allele];
            // how many blocks are before current position
            if (countAtThisAllele == 0) {
                continue;
            }
            for (int smallerCountAtThisAllele=0; smallerCountAtThisAllele < countAtThisAllele; smallerCountAtThisAllele++) {
                numGenotypesWithLowerIndex += GLVectorSizeCache.getNumLikelihoodElements(allele, remainingPloidy - smallerCountAtThisAllele);
            }

            remainingPloidy -= countAtThisAllele;
        }

        return numGenotypesWithLowerIndex;

    }
}

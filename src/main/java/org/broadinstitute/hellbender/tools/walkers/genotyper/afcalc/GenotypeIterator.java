package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import java.util.Arrays;

/**
 * For a given ploidy (chromosome count) and number of alleles, we need a form to iterate deterministically
 * across all possible allele conformations.
 * Problem equivalent to listing in deterministic order all possible ways in which N integers will sum to P,
 * where N is number of alleles and P is number of chromosomes.
 * For example, with P=2,N=2, iterator will produce
 * [2 0 ] [1 1] [ 0 2]
 */
public final class GenotypeIterator {
    private final int ploidy;   // this is always > 0 in calling code
    private final int numAlleles;
    private int[] alleleCounts;
    private boolean hasNext;
    private int linearIndex;

    public GenotypeIterator(final int numAlleles, final int ploidy) {
        this.ploidy = ploidy;
        this.numAlleles = numAlleles;
        alleleCounts = new int[numAlleles];
        reset();
    }

    public void next() {
        hasNext = false;
        for (int allele=1; allele < numAlleles; allele++) {
            if (alleleCounts[0] == 0) {
                alleleCounts[0] += alleleCounts[allele];
                alleleCounts[allele] = 0;
            }
            else {
                alleleCounts[allele]++;
                alleleCounts[0]--;
                hasNext = true;
                break;
            }
        }

        if (hasNext) {
            linearIndex++;
        }
    }

    public void reset() {
        Arrays.fill(alleleCounts, 0);
        alleleCounts[0] = ploidy;
        hasNext = true;
        linearIndex = 0;
    }
    public int[] getCurrentVector() {
        return alleleCounts;
    }
    public int getLinearIndex() {
        return linearIndex;
    }
    public boolean hasNext() {return hasNext;}
}

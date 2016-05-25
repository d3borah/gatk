package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.broadinstitute.hellbender.utils.MathUtils;

import java.util.Arrays;

/**
 * Class that produces allele-frequency priors.
 */
public class AFPriorProvider {

    private double[][] priorByTotalPloidy;

    private final double heterozygosity;
    private final double log10Heterozygosity;

    /**
     * Construct a new provider given the heterozygosity value.
     * @param heterozygosity must be a valid heterozygosity between larger than 0 and smaller than 1.
     * @throws IllegalArgumentException if {@code heterozygosity} is not valid one in the interval (0,1).
     */
    public AFPriorProvider(final double heterozygosity) {
        if (heterozygosity <= 0) {
            throw new IllegalArgumentException("the heterozygosity must be greater than 0");
        }
        if (heterozygosity >= 1) {
            throw new IllegalArgumentException("the heterozygosity must be less than 1");
        }
        if (Double.isNaN(heterozygosity)) {
            throw new IllegalArgumentException("the heterozygosity cannot be a NaN");
        }
        this.heterozygosity = heterozygosity;
        this.log10Heterozygosity = Math.log10(heterozygosity);
    }

    /**
     * Returns the priors given a total-ploidy (the total number of genome copies across all samples).
     *
     * <p>
     *     For performance sake the returned value is a direct reference ot the cached prior, so the client code
     *     must not modify its content.
     * </p>
     *
     * @param totalPloidy the requested total ploidy.
     *
     * @return never {@code null}, please do not modify its content. An array of {@code totalPloidy + 1} positions where
     *  the ith position is the log10(prior) probability of the an alternative allele AC to be exactly <i>i</i> copies in
     *  a draw of {@code totalPloidy} elements.
     */
    public double[] forTotalPloidy(final int totalPloidy) {
        if (totalPloidy < 0) {
            throw new IllegalArgumentException("the total-ploidy cannot be negative");
        }
        ensureCapacity(totalPloidy);
        final double[] cachedResult = priorByTotalPloidy[totalPloidy];
        if (cachedResult == null) {
            return priorByTotalPloidy[totalPloidy] = buildPriors(totalPloidy);
        } else {
            return cachedResult;
        }
    }

    /**
     * Make sure that structures have enough capacity to hold the information up to the given total-ploidy.
     * @param totalPloidy
     */
    protected void ensureCapacity(final int totalPloidy) {
        if (totalPloidy < 0) {
            throw new IllegalArgumentException("the total-ploidy cannot be negative");
        }
        if (priorByTotalPloidy == null) {
            priorByTotalPloidy = new double[totalPloidy + 1][];  // just enough for those cases in where we have a fix total-ploidy.
        } else if (priorByTotalPloidy.length - 1 < totalPloidy) {
            priorByTotalPloidy = Arrays.copyOf(priorByTotalPloidy, Math.max(priorByTotalPloidy.length << 1, totalPloidy + 1));
        }
    }

    /**
     * Given a total ploidy construct the allele prior probabilities array.
     * @param totalPloidy the target total-ploidy. Code can assume that is a non-negative number.
     *
     * @return never {@code null}, an array of exactly {@code totalPloidy + 1} positions that satisifed the
     *  contract {@link #forTotalPloidy(int) forTotalPloidy(totalPloidy)}.
     */
    protected double[] buildPriors(final int totalPloidy) {
        final double[] result = new double [totalPloidy + 1];
        Arrays.fill(result, log10Heterozygosity);
        result[0] = Double.NEGATIVE_INFINITY;
        for (int i = 1; i <= totalPloidy; i++) {
            result[i] -= MathUtils.log10(i);
        }
        final double log10Sum = MathUtils.approximateLog10SumLog10(result);
        if (log10Sum >= 0) {
            throw new IllegalArgumentException("heterozygosity " + heterozygosity + " is too large of total ploidy " + totalPloidy);
        }
        result[0] = MathUtils.log10OneMinusPow10(log10Sum);
        return result;
    }

}

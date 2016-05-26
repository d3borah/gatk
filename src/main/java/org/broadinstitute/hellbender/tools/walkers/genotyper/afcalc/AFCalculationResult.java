package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;

/**
 * Describes the results of the AFCalc
 *
 * Only the bare essentials are represented here, as all AFCalc models must return meaningful results for
 * all of these fields.
 *
 * Note that all of the values -- i.e. priors -- are checked now that they are meaningful, which means
 * that users of this code can rely on the values coming out of these functions.
 */
public final class AFCalculationResult {
    private static final int AF0 = 0;
    private static final int AF1p = 1;

    private final double[] log10PosteriorsOfAC;

    private final Map<Allele, Double> log10pRefByAllele;

    /**
     * The AC values for all ALT alleles at the MLE
     */
    private final int[] alleleCountsOfMLE;

    /**
     * The list of alleles actually used in computing the AF
     */
    private final List<Allele> allelesUsedInGenotyping;

    /**
     * Create a results object capability of storing results for calls with up to maxAltAlleles
     */
    public AFCalculationResult(final int[] alleleCountsOfMLE,
                               final List<Allele> allelesUsedInGenotyping,
                               final double[] log10LikelihoodsOfAC,
                               final double[] log10PriorsOfAC,
                               final Map<Allele, Double> log10pRefByAllele) {
        this.alleleCountsOfMLE = Utils.nonNull(alleleCountsOfMLE, "alleleCountsOfMLE cannot be null").clone();
        this.allelesUsedInGenotyping = Collections.unmodifiableList(new ArrayList<>(allelesUsedInGenotyping));
        this.log10PosteriorsOfAC = computePosteriors(log10LikelihoodsOfAC, log10PriorsOfAC);
        this.log10pRefByAllele = Collections.unmodifiableMap(new LinkedHashMap<>(log10pRefByAllele));
    }

    /**
     * Returns the AC of allele a la #getAlleleCountsOfMLE
     *
     * @param allele the allele whose AC we want to know.  Error if its not in allelesUsedInGenotyping
     * @throws IllegalStateException if allele isn't in allelesUsedInGenotyping
     * @return the AC of allele
     */
    public int getAlleleCountAtMLE(final Allele allele) {
        Utils.nonNull(allele);
        return alleleCountsOfMLE[altAlleleIndex(allele)];
    }

    /**
     * Get the list of alleles actually used in genotyping.
     *
     * Due to computational / implementation constraints this may be smaller than
     * the actual list of alleles requested
     *
     * @return a non-empty list of alleles used during genotyping, the first of which is the reference allele
     */
    public List<Allele> getAllelesUsedInGenotyping() {
        return allelesUsedInGenotyping;
    }

    /**
     * Get the log10 normalized -- across all ACs -- posterior probability of AC == 0 for all alleles
     */
    public double getLog10PosteriorOfAFEq0() {
        return log10PosteriorsOfAC[AF0];
    }

    /**
     * Get the log10 normalized -- across all ACs -- posterior probability of AC > 0 for any alleles
     */
    public double getLog10PosteriorOfAFGT0() {
        return log10PosteriorsOfAC[AF1p];
    }

    /**
     * Are we sufficiently confident in being non-ref that the site is considered polymorphic?
     *
     * We are non-ref if the probability of being non-ref > the emit confidence (often an argument).
     * Suppose posterior AF > 0 is log10: -5 => 10^-5
     * And that log10minPNonRef is -3.
     * We are considered polymorphic since 10^-5 < 10^-3 => -5 < -3
     *
     * Note that log10minPNonRef is really the minimum confidence, scaled as an error rate, so
     * if you want to be 99% confidence, then log10PNonRef should be log10(0.01) = -2.
     *
     * @param log10minPNonRef the log10 scaled min pr of being non-ref to be considered polymorphic
     *
     * @return true if there's enough confidence (relative to log10minPNonRef) to reject AF == 0
     */
    public boolean isPolymorphic(final Allele allele, final double log10minPNonRef) {
        Utils.nonNull(allele);
        return getLog10PosteriorOfAFEq0ForAllele(allele) < log10minPNonRef;
    }

    /**
     * Same as #isPolymorphic but takes a phred-scaled quality score as input
     */
    public boolean isPolymorphicPhredScaledQual(final Allele allele, final double minPNonRefPhredScaledQual) {
        Utils.nonNull(allele);
        if ( minPNonRefPhredScaledQual < 0 ) {
            throw new IllegalArgumentException("phredScaledQual " + minPNonRefPhredScaledQual + " < 0 ");
        }
        final double log10Threshold = minPNonRefPhredScaledQual / -10;
        return isPolymorphic(allele, log10Threshold);
    }

    /**
     * Returns the log10 probability that allele is not segregating
     *
     * Note that this function is p not segregating so that we can store
     * internally the log10 value of AF == 0, which grows very quickly
     * negative and yet has sufficient resolution for high confidence tests.
     * For example, if log10pRef == -100, not an unreasonably high number,
     * if we tried to store log10pNonRef we'd be looking at 1 - 10^-100, which
     * quickly underflows to 1.  So the logic here is backward from what
     * you really want (the p of segregating) but we do that for numerical
     * reasons
     *
     * Unlike the sites-level annotation, this calculation is specific to allele, and can be
     * used to separately determine how much evidence there is that allele is independently
     * segregating as opposed to the site being polymorphic with any allele.  In the bi-allelic
     * case these are obviously the same but for multiple alt alleles there can be lots of
     * evidence for one allele but not so much for any other allele
     *
     * @param allele the allele we're interested in, must be in getAllelesUsedInGenotyping
     * @return the log10 probability that allele is not segregating at this site
     */
    public double getLog10PosteriorOfAFEq0ForAllele(final Allele allele) {
        Utils.nonNull(allele);
        final Double log10pNonRef = log10pRefByAllele.get(allele);
        Utils.nonNull(log10pNonRef, "Unknown allele " + allele);
        return log10pNonRef;
    }

    /**
     * Returns the log10 normalized posteriors given the log10 likelihoods and priors
     *
     * @param log10LikelihoodsOfAC
     * @param log10PriorsOfAC
     *
     * @return freshly allocated log10 normalized posteriors vector
     */
    private static double[] computePosteriors(final double[] log10LikelihoodsOfAC, final double[] log10PriorsOfAC) {
        final double[] log10UnnormalizedPosteriors = new double[log10LikelihoodsOfAC.length];
        for ( int i = 0; i < log10LikelihoodsOfAC.length; i++ ) {
            log10UnnormalizedPosteriors[i] = log10LikelihoodsOfAC[i] + log10PriorsOfAC[i];
        }
        return MathUtils.normalizeFromLog10(log10UnnormalizedPosteriors, true, false);
    }

    /**
     * Computes the offset into linear vectors indexed by alt allele for allele
     *
     * Things like our MLE allele count vector are indexed by alt allele index, with
     * the first alt allele being 0, the second 1, etc.  This function computes the index
     * associated with allele.
     *
     * @param allele the allele whose alt index we'd like to know
     * @throws IllegalArgumentException if allele isn't in allelesUsedInGenotyping
     * @return an index value greater than 0 suitable for indexing into the MLE and other alt allele indexed arrays
     */
    private int altAlleleIndex(final Allele allele) {
        if ( allele.isReference() ) {
            throw new IllegalArgumentException("Cannot get the alt allele index for reference allele " + allele);
        }
        final int index = allelesUsedInGenotyping.indexOf(allele);
        if ( index == -1 ) {
            throw new IllegalArgumentException("could not find allele " + allele + " in " + allelesUsedInGenotyping);
        } else {
            return index - 1;
        }
    }
}
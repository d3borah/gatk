package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.GenotypesContext;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.List;


/**
 * Generic interface for calculating the probability of alleles segregating given priors and genotype likelihoods
 */
public abstract class AFCalculator {

    protected static final Logger logger = LogManager.getLogger(AFCalculator.class);

    /**
     * Compute the probability of the alleles segregating given the genotype likelihoods of the samples in vc
     *
     * @param vc the VariantContext holding the alleles and sample information.  The VariantContext
     *           must have at least 1 alternative allele
     * @param log10AlleleFrequencyPriors a prior vector nSamples x 2 in length indicating the Pr(AF = i)
     * @return result (for programming convenience)
     */
    public AFCalculationResult getLog10PNonRef(final VariantContext vc, final int defaultPloidy, final int maximumAlternativeAlleles, final double[] log10AlleleFrequencyPriors) {
        Utils.nonNull(vc, "VariantContext cannot be null");
        Utils.nonNull(log10AlleleFrequencyPriors, "priors vector cannot be null");
        if ( vc.getNAlleles() == 1 ) {
            throw new IllegalArgumentException("VariantContext has only a single reference allele, but getLog10PNonRef requires at least one at all " + vc);
        }

        final VariantContext vcWorking = reduceScope(vc,defaultPloidy, maximumAlternativeAlleles);

        return computeLog10PNonRef(vcWorking, defaultPloidy, log10AlleleFrequencyPriors);
    }

    // ---------------------------------------------------------------------------
    //
    // Abstract methods that should be implemented by concrete implementations
    // to actually calculate the AF
    //
    // ---------------------------------------------------------------------------

    /**
     * Look at VC and perhaps return a new one of reduced complexity, if that's necessary
     *
     * Used before the call to computeLog10PNonRef to simply the calculation job at hand,
     * if vc exceeds bounds.  For example, if VC has 100 alt alleles this function
     * may decide to only genotype the best 2 of them.
     *
     * @param vc the initial VC provided by the caller to this AFcalculation
     * @return a potentially simpler VC that's more tractable to genotype
     */
    protected abstract VariantContext reduceScope(final VariantContext vc, final int defaultPloidy, final int maximumAlternativeAlleles);

    /**
     * Actually carry out the log10PNonRef calculation on vc, storing results in results
     *
     * @param vc                                variant context with alleles and genotype likelihoods,
     *                                          must have at least one alt allele
     * @param log10AlleleFrequencyPriors        priors
     * @return a AFCalcResult object describing the results of this calculation
     */
    protected abstract AFCalculationResult computeLog10PNonRef(final VariantContext vc, final int defaultPloidy,
                                                        final double[] log10AlleleFrequencyPriors);

    /**
     * Subset VC to the just allelesToUse, updating genotype likelihoods
     *
     * Must be overridden by concrete subclasses
     *
     * @param vc                                variant context with alleles and genotype likelihoods
     * @param defaultPloidy                     default ploidy to assume in case {@code vc} does not indicate it for a sample.
     * @param allelesToUse                      alleles to subset
     * @param assignGenotypes
     * @return GenotypesContext object
     */
    public abstract GenotypesContext subsetAlleles(final VariantContext vc,
                                                   final int defaultPloidy,
                                                   final List<Allele> allelesToUse,
                                                   final boolean assignGenotypes);
}
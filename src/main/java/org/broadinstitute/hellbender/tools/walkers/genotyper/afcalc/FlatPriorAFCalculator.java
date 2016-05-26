package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeAlleleCounts;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculator;
import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotypeLikelihoodCalculators;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * AFCalculator in which genotype posterior probabilities are equal to their likelihoods, normalized.
 *
 * For now I'm just using this as a simple implementation to hook up to the GenotypingEngine in order to delete all
 * the other implementations and their associated classes like StateTracker and ExactACSet.
 */
public class FlatPriorAFCalculator extends AFCalculator {
    @Override
    protected AFCalculationResult computeLog10PNonRef(final VariantContext vc, final int defaultPloidy) {
        Utils.nonNull(vc, "vc is null");

        final int[] alleleCountsOfMLE = null;
        final List<Allele> allelesUsedInGenotyping = vc.getAlleles();
        final double log10PosteriorOfAFEq0 = -1.0;    //TODO: placeholder
        final Map<Allele, Double> log10pRefByAllele = null;

        final GenotypesContext GLs = vc.getGenotypes();
        final int numAlleles = vc.getNAlleles();

        for (final Genotype genotype : GLs.iterateInSampleNameOrder()) {
            if (!genotype.hasPL()) {
                continue;
            }
            final double[] gls = genotype.getLikelihoods().getAsVector();
            final double[] genotypePosteriors = MathUtils.normalizeFromLog10(gls, false);   // real, not log

        }

        return new AFCalculationResult(alleleCountsOfMLE, allelesUsedInGenotyping, log10PosteriorOfAFEq0, log10pRefByAllele);
    }


    /**
     * A WHOLE BUNCH OF CODE COPIED FROM ExactAFCalculator and GeneralPloidyExactAFCalculator to implement reduceScope()
     *
     * I copied for the sake of expedience since this logic is seprate from the exact model.  I'm sure that I will want to destroy this.
     */

    protected static final int PL_INDEX_OF_HOM_REF = 0;

    /**
     * Sorts {@link LikelihoodSum} instances where those with higher likelihood are first.
     */
    protected static final Comparator<LikelihoodSum> LIKELIHOOD_SUM_COMPARATOR = Comparator.<LikelihoodSum>comparingDouble(o->o.sum).reversed();

    /**
     * Sorts {@link LikelihoodSum} instances where those with higher likelihood are first but make sure that
     * NON_REF alleles are last.
     */
    protected static final Comparator<LikelihoodSum> LIKELIHOOD_NON_REF_THEN_SUM_COMPARATOR = (o1, o2) -> {
        if (o1.allele == GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE) {
            return 1;
        } else if (o2.allele == GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE) {
            return -1;
        } else {
            return o1.compareTo(o2);
        }
    };
    /**
     * Sorts {@link LikelihoodSum} instances where those with lower alternative allele index are first regardless of
     * the likelihood sum.
     */
    protected static final Comparator<LikelihoodSum> LIKELIHOOD_INDEX_COMPARATOR = Comparator.comparingInt(o->o.index);



    /**
     * Look at VC and perhaps return a new one of reduced complexity, if that's necessary
     *
     * Used before the call to computeLog10PNonRef to simply the calculation
     *
     * @param vc the initial VC provided by the caller to this AFcalculation
     * @return a potentially simpler VC that's more tractable to genotype
     */
    @Override
    protected VariantContext reduceScope(final VariantContext vc, final int defaultPloidy, final int maximumAlternativeAlleles) {
        Utils.nonNull(vc, "vc is null");
        final List<Allele> inputAltAlleles = vc.getAlternateAlleles();
        final List<Allele> outputAltAlleles = reduceScopeAlleles(vc, defaultPloidy, maximumAlternativeAlleles);

        if (inputAltAlleles.size() == outputAltAlleles.size()) {
            return vc;
        } else if (inputAltAlleles.size() < outputAltAlleles.size()) {
            throw new IllegalStateException("unexpected: reduction increased the number of alt. alleles!");
        }
        logger.warn("this tool is currently set to genotype at most " + maximumAlternativeAlleles
                + " alternate alleles in a given context, but the context at " + vc.getContig() + ":" + vc.getStart()
                + " has " + (vc.getAlternateAlleles().size())
                + " alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument");

        final List<Allele> alleles = Stream.concat(Stream.of(vc.getReference()), outputAltAlleles.stream()).collect(Collectors.toList());
        final GenotypesContext genotypes = subsetAlleles(vc,defaultPloidy,alleles,false);
        return new VariantContextBuilder(vc).alleles(alleles).genotypes(genotypes).make();
    }



    /**
     * Returns a the new set of alleles to use.
     * @param vc target variant context.
     * @param numAllelesToChoose number of alleles to keep.
     * @return the list of alternative allele to keep.
     */
    protected List<Allele> reduceScopeAlleles(final VariantContext vc, final int defaultPloidy, final int numAllelesToChoose) {
        Utils.nonNull(vc, "vc is null");

        final int numOriginalAltAlleles = vc.getAlternateAlleles().size();

        final int nonRefAltAlleleIndex = GATKVariantContextUtils.indexOfAltAllele(vc, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE, false);
        final boolean nonRefAltAllelePresent = nonRefAltAlleleIndex >= 0;

        // <NON_REF> should not be considered in the downsizing
        final int numProperOriginalAltAlleles = numOriginalAltAlleles - (nonRefAltAllelePresent ? 1 : 0);

        // Avoid pointless allele reduction:
        if (numAllelesToChoose >= numProperOriginalAltAlleles) {
            return vc.getAlternateAlleles();
        }

        final double[] likelihoodSums = new double[numOriginalAltAlleles + 1]; //likelihood sums for all alleles including ref
        for ( final Genotype genotype : vc.getGenotypes().iterateInSampleNameOrder() ) {
            final double[] gls = genotype.getLikelihoods().getAsVector();
            if (gls == null || GATKVariantContextUtils.likelihoodsAreUninformative(gls)) {
                continue;
            }

            final int indexOfBestGL = MathUtils.maxElementIndex(gls);

            final double GLDiffBetweenRefAndBest = gls[indexOfBestGL] - gls[PL_INDEX_OF_HOM_REF];
            final int ploidy = genotype.getPloidy() > 0 ? genotype.getPloidy() : defaultPloidy;

            final int[] acCount = getAlleleCountFromPLIndex(1 + numOriginalAltAlleles, ploidy, indexOfBestGL);
            // by convention, first count coming from getAlleleCountFromPLIndex comes from reference allele
            for (int k=1; k < acCount.length; k++) {
                if (acCount[k] > 0) {
                    likelihoodSums[k] += acCount[k] * GLDiffBetweenRefAndBest;
                }
            }
        }
        Collections.sort(likelihoodSums, nonRefAltAllelePresent ? LIKELIHOOD_NON_REF_THEN_SUM_COMPARATOR : LIKELIHOOD_SUM_COMPARATOR);

        // Keep track of original index order of likelihoods via a heap
        final PriorityQueue<LikelihoodSum> mostLikelyAllelesHeapByIndex = new PriorityQueue<>(numOriginalAltAlleles, LIKELIHOOD_INDEX_COMPARATOR);
        mostLikelyAllelesHeapByIndex.addAll(likelihoodSums.subList(0, numAllelesToChoose));


        // guaranteed no to have been added at this point thanks for checking on whether reduction was
        // needed in the first place.
        if (nonRefAltAllelePresent) {
            mostLikelyAllelesHeapByIndex.add(likelihoodSums.get(nonRefAltAlleleIndex));
        }

        final List<Allele> orderedBestAlleles = new ArrayList<>(numAllelesToChoose);

        while (!mostLikelyAllelesHeapByIndex.isEmpty()) {
            orderedBestAlleles.add(mostLikelyAllelesHeapByIndex.remove().allele);
        }

        return orderedBestAlleles;
    }

    /**
     * Given a scalar index, what's the alelle count conformation corresponding to it?
     * @param nAlleles                    Number of alleles
     * @param numChromosomes              Ploidy
     * @param PLindex                     Index to query
     * @return                            Allele count conformation, according to iteration order from SumIterator
     */
    private static int[] getAlleleCountFromPLIndex(final int nAlleles, final int numChromosomes, final int PLindex) {
        final GenotypeLikelihoodCalculator calculator = new GenotypeLikelihoodCalculators().getInstance(numChromosomes, nAlleles);
        final GenotypeAlleleCounts alleleCounts = calculator.genotypeAlleleCountsAt(PLindex);
        return alleleCounts.alleleCountsByIndex(nAlleles - 1);
    }

    /**
     * Wrapper class that compares two likelihoods associated with two alleles
     */
    protected static final class LikelihoodSum implements Comparable<LikelihoodSum> {
        public double sum = 0.0;
        public final Allele allele;
        public final int index;

        public LikelihoodSum(final Allele allele, final int index) { this.allele = allele; this.index = index; }

        public int compareTo(final LikelihoodSum other) {
            final double diff = Double.compare(sum, other.sum);
            return ( diff < 0.0 ) ? 1 : (diff > 0.0 ) ? -1 : 0;
        }
    }

    /**
     * END OF A WHOLE BUNCH OF CODE COPIED FROM ExactAFCalculator and GeneralPloidyExactAFCalculator to implement reduceScope()
     *
     */

    /**
     * A WHOLE BUNCH OF CODE COPIED FROM ExactAFCalculator and GeneralPloidyExactAFCalculator to implement subsetAlleles()
     *
     */


    static final int MAX_LENGTH_FOR_POOL_PL_LOGGING = 100; // if PL vectors longer than this # of elements, don't log them
    /**
     * From a given variant context, extract a given subset of alleles, and update genotype context accordingly,
     * including updating the PL's, and assign genotypes accordingly
     * @param vc                                variant context with alleles and genotype likelihoods
     * @param defaultPloidy                     ploidy to assume in case that {@code vc} does not contain that information
     *                                          for a sample.
     * @param allelesToUse                      alleles to subset
     * @param assignGenotypes                   true: assign hard genotypes, false: leave as no-call
     * @return                                  GenotypesContext with new PLs
     */
    public GenotypesContext subsetAlleles(final VariantContext vc, final int defaultPloidy, final List<Allele> allelesToUse,
                                          final boolean assignGenotypes) {
        Utils.nonNull(vc, "vc is null");
        Utils.nonNull(allelesToUse, "allelesToUse is null");

        // we need to determine which of the alternate alleles (and hence the likelihoods) to use and carry forward
        final int numOriginalAltAlleles = vc.getAlternateAlleles().size();
        final int numNewAltAlleles = allelesToUse.size() - 1;


        // create the new genotypes
        final GenotypesContext newGTs = GenotypesContext.create();
        for ( final Genotype g : vc.getGenotypes().iterateInSampleNameOrder() ) {
            final int ploidy = g.getPloidy() > 0 ? g.getPloidy() : defaultPloidy;
            if ( !g.hasLikelihoods() ) {
                newGTs.add(GenotypeBuilder.create(g.getSampleName(), GATKVariantContextUtils.noCallAlleles(ploidy)));
                continue;
            }

            // create the new likelihoods array from the alleles we are allowed to use
            // Optimization: if no new alt alleles (pure ref call), keep original likelihoods and skip normalization and subsetting
            final double[] originalLikelihoods = g.getLikelihoods().getAsVector();
            double[] newLikelihoods = numOriginalAltAlleles == numNewAltAlleles || numNewAltAlleles == 0 ? originalLikelihoods :
                    MathUtils.normalizeFromLog10(subsetToAlleles(originalLikelihoods, ploidy, vc.getAlleles(), allelesToUse), false, true);

            // if there is no mass on the (new) likelihoods, then just no-call the sample
            if ( GATKVariantContextUtils.likelihoodsAreUninformative(newLikelihoods) ) {
                newGTs.add(GenotypeBuilder.create(g.getSampleName(), GATKVariantContextUtils.noCallAlleles(ploidy)));
            }
            else {
                final GenotypeBuilder gb = new GenotypeBuilder(g);

                //TODO: is it really necessary to handle numNewAltAlleles == 0 separately?
                gb.PL(numNewAltAlleles == 0 ? null : newLikelihoods);

                // if we weren't asked to assign a genotype, then just no-call the sample
                if ( !assignGenotypes || MathUtils.sum(newLikelihoods) > GATKVariantContextUtils.SUM_GL_THRESH_NOCALL ) {
                    gb.alleles(GATKVariantContextUtils.noCallAlleles(ploidy));
                } else {
                    assignGenotype(gb, newLikelihoods, allelesToUse, ploidy);
                }
                newGTs.add(gb.make());
            }
        }

        return newGTs;

    }

    /**
     * Given set of alleles with corresponding vector of likelihoods, subset to a new set of alleles
     *
     * @param oldLikelihoods        Vector of PL's corresponding to original alleles
     * @param numChromosomes        Ploidy (number of chromosomes describing PL's)
     * @param originalAlleles       List of original alleles
     * @param allelesToSubset       Alleles to subset
     * @return                      Vector of new PL's, ordered accorrding to SumIterator's ordering
     */
    private static double[] subsetToAlleles(final double[] oldLikelihoods, final int numChromosomes,
                                            final List<Allele> originalAlleles, final List<Allele> allelesToSubset) {

        final int newPLSize = getNumLikelihoodElements(allelesToSubset.size(), numChromosomes);
        final double[] newPLs = new double[newPLSize];

        // First fill boolean array stating whether each original allele is present in new mapping
        final List<Boolean> allelePresent = originalAlleles.stream().map(a -> allelesToSubset.contains(a)).collect(Collectors.toList());

        // compute mapping from old idx to new idx
        // Example. Original alleles: {T*,C,G,A}. New alleles: {G,C}. Permutation key = [2,1]
        final int[] oldIndicesOfNewAlleles = allelesToSubset.stream().mapToInt(a -> originalAlleles.indexOf(a)).toArray();

        final SumIterator iterator = new SumIterator(originalAlleles.size(),numChromosomes);
        while (iterator.hasNext()) {
            // for each entry in logPL table, associated originally with allele count stored in vec[],
            // see if this allele count conformation will be present in new logPL table.
            // For entry to be present, elements in dimensions not present in requested allele list have to have count = 0
            final int[] pVec = iterator.getCurrentVector();
            final double pl = oldLikelihoods[iterator.getLinearIndex()];

            boolean keyPresent = true;
            for (int k=0; k < allelePresent.size(); k++) {
                if (pVec[k] > 0 && !allelePresent.get(k)) {
                    keyPresent = false;
                }
            }

            if (keyPresent) {// skip to next entry in logPLs if this conformation is not present in subset

                final int[] newCount = new int[allelesToSubset.size()];

                // map from old allele mapping count to new allele mapping
                // In pseudo-Matlab notation: newCount = vec[oldIndicesOfNewAlleles] for oldIndicesOfNewAlleles vector
                for (int idx = 0; idx < newCount.length; idx++) {
                    newCount[idx] = pVec[oldIndicesOfNewAlleles[idx]];
                }

                // get corresponding index from new count
                final int outputIdx = getLinearIndex(newCount, allelesToSubset.size(), numChromosomes);
                newPLs[outputIdx] = pl;
            }
            iterator.next();
        }

        return  newPLs;
    }
    /*
* a cache of the PL vector sizes as a function of # of alleles and ploidy
*/
    @VisibleForTesting
    static int getNumLikelihoodElements(final int numAlleles, final int ploidy) {
        return GLVECTORSIZES[numAlleles][ploidy];
    }

    private static final int MAX_NUM_ALLELES_TO_CACHE = 20;
    private static final int MAX_NUM_SAMPLES_PER_POOL = 1000;

    //Note: this is shared state but it's not modified at runtime
    private static final int[][] GLVECTORSIZES = fillGLVectorSizeCache(MAX_NUM_ALLELES_TO_CACHE, 2*MAX_NUM_SAMPLES_PER_POOL);

    private static int[][] fillGLVectorSizeCache(final int maxAlleles, final int maxPloidy) {
        final int[][] cache = new int[maxAlleles][maxPloidy];
        for (int numAlleles=1; numAlleles < maxAlleles; numAlleles++) {
            for (int ploidy=0; ploidy < maxPloidy; ploidy++) {
                cache[numAlleles][ploidy] = numAlleles == 1 ? 1 : Arrays.stream(cache[numAlleles - 1], 0, ploidy + 1).sum();
            }
        }
        return cache;
    }

    private static int getLinearIndex(final int[] vectorIdx, final int numAlleles, final int ploidy) {
        if (ploidy <= 0) {
            return 0;
        }

        int linearIdx = 0;
        int cumSum = ploidy;
        for (int k = numAlleles - 1; k >= 1; k--) {
            final int idx = vectorIdx[k];
            // how many blocks are before current position
            if (idx == 0) {
                continue;
            }
            for (int p=0; p < idx; p++) {
                linearIdx += getNumLikelihoodElements(k, cumSum - p);
            }

            cumSum -= idx;
        }

        return linearIdx;

    }


    /**
     * Assign genotypes (GTs) to the samples in the Variant Context greedily based on the PLs
     *
     * @param newLikelihoods       the PL array
     * @param allelesToUse         the list of alleles to choose from (corresponding to the PLs)
     * @param numChromosomes        Number of chromosomes per pool
     */
    private static void assignGenotype(final GenotypeBuilder gb, final double[] newLikelihoods, final List<Allele> allelesToUse,
                                       final int numChromosomes) {
        final int numNewAltAlleles = allelesToUse.size() - 1;

        // find the genotype with maximum likelihoods
        final int PLindex = numNewAltAlleles == 0 ? 0 : MathUtils.maxElementIndex(newLikelihoods);
        final GenotypeLikelihoodCalculator calculator = new GenotypeLikelihoodCalculators().getInstance(numChromosomes, allelesToUse.size());
        final GenotypeAlleleCounts alleleCounts = calculator.genotypeAlleleCountsAt(PLindex);
        gb.alleles(alleleCounts.asAlleleList(allelesToUse));

        // remove PLs if necessary
        if (newLikelihoods.length > MAX_LENGTH_FOR_POOL_PL_LOGGING) {
            gb.noPL();
        }

        if ( numNewAltAlleles > 0 ) {
            gb.log10PError(GenotypeLikelihoods.getGQLog10FromLikelihoods(PLindex, newLikelihoods));
        }
    }

    /**
     * END OF A WHOLE BUNCH OF CODE COPIED FROM ExactAFCalculator and GeneralPloidyExactAFCalculator to implement subsetAlleles()
     *
     */
}

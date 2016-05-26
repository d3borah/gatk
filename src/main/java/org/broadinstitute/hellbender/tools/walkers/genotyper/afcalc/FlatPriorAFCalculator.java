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

/**
 * AFCalculator in which genotype posterior probabilities are equal to their likelihoods, normalized.
 *
 * For now I'm just using this as a simple implementation to hook up to the GenotypingEngine in order to delete all
 * the other implementations and their associated classes like StateTracker and ExactACSet.
 */
public class FlatPriorAFCalculator extends AFCalculator {


    @Override
    protected AFCalculationResult computeLog10PNonRef(final VariantContext vc, final int defaultPloidy,
                                                      final double[] log10AlleleFrequencyPriors) {
        Utils.nonNull(vc, "vc is null");
        Utils.nonNull(log10AlleleFrequencyPriors, "log10AlleleFrequencyPriors is null");

        final int[] alleleCountsOfMLE = null;
        final List<Allele> allelesUsedInGenotyping = vc.getAlleles();
        final double[] log10LikelihoodsOfAC = null; //refactor AFCalculationResult to not use this!!
        final double[] log10PriorsOfAC = null;
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

        return new AFCalculationResult(alleleCountsOfMLE, allelesUsedInGenotyping, log10LikelihoodsOfAC, log10PriorsOfAC, log10pRefByAllele);
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
     * Used before the call to computeLog10PNonRef to simply the calculation job at hand,
     * if vc exceeds bounds.  For example, if VC has 100 alt alleles this function
     * may decide to only genotype the best 2 of them.
     *
     * @param vc the initial VC provided by the caller to this AFcalculation
     * @return a potentially simpler VC that's more tractable to genotype
     */
    @Override
    protected VariantContext reduceScope(final VariantContext vc, final int defaultPloidy, final int maximumAlternativeAlleles) {
        Utils.nonNull(vc, "vc is null");
        // don't try to genotype too many alternate alleles
        final List<Allele> inputAltAlleles = vc.getAlternateAlleles();
        final List<Allele> outputAltAlleles = reduceScopeAlleles(vc, defaultPloidy, maximumAlternativeAlleles);

        // only if output allele has reduced from the input alt allele set size we should care.
        final int altAlleleReduction = inputAltAlleles.size() - outputAltAlleles.size();

        if (altAlleleReduction == 0) {
            return vc;
        }
        logger.warn("this tool is currently set to genotype at most " + maximumAlternativeAlleles
                + " alternate alleles in a given context, but the context at " + vc.getContig() + ":" + vc.getStart()
                + " has " + (vc.getAlternateAlleles().size())
                + " alternate alleles so only the top alleles will be used; see the --max_alternate_alleles argument");

        final List<Allele> alleles = new ArrayList<>(maximumAlternativeAlleles + 1);
        alleles.add(vc.getReference());
        alleles.addAll(reduceScopeAlleles(vc, defaultPloidy, maximumAlternativeAlleles));
        final VariantContextBuilder builder = new VariantContextBuilder(vc);
        builder.alleles(alleles);
        builder.genotypes(reduceScopeGenotypes(vc, defaultPloidy, alleles));
        if (altAlleleReduction < 0) {
            throw new IllegalStateException("unexpected: reduction increased the number of alt. alleles!: " + -altAlleleReduction + " " + vc + " " + builder.make());
        }
        return builder.make();
    }



    /**
     * Returns a the new set of alleles to use.
     * @param vc target variant context.
     * @param numAllelesToChoose number of alleles to keep.
     * @return the list of alternative allele to keep.
     */
    protected List<Allele> reduceScopeAlleles(final VariantContext vc, final int defaultPloidy, final int numAllelesToChoose) {
        Utils.nonNull(vc, "vc is null");

        // Look  for the <NON_REF> allele to exclude it from the pruning if present.
        final int numOriginalAltAlleles = vc.getAlternateAlleles().size();

        final int nonRefAltAlleleIndex = GATKVariantContextUtils.indexOfAltAllele(vc, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE, false);
        final boolean nonRefAltAllelePresent = nonRefAltAlleleIndex >= 0;

        // <NON_REF> should not be considered in the downsizing, so we need to count it out when
        // considering if alt. allele downsizing is required.
        final int numProperOriginalAltAlleles = numOriginalAltAlleles - (nonRefAltAllelePresent ? 1 : 0);

        // Avoid pointless allele reduction:
        if (numAllelesToChoose >= numProperOriginalAltAlleles) {
            return vc.getAlternateAlleles();
        }

        final LikelihoodSum[] likelihoodSums = new LikelihoodSum[numOriginalAltAlleles];
        for ( int i = 0; i < numOriginalAltAlleles; i++ ) {
            final Allele allele = vc.getAlternateAllele(i);
            likelihoodSums[i] = new LikelihoodSum(allele,i);
        }

        // Calculate the allele likelihood sums.
        reduceScopeCalculateLikelihoodSums(vc, defaultPloidy, likelihoodSums);

        // sort them by probability mass and choose the best ones
        // Make sure that the <NON_REF> allele is last if present.
        Collections.sort(Arrays.asList(likelihoodSums), nonRefAltAllelePresent ? LIKELIHOOD_NON_REF_THEN_SUM_COMPARATOR : LIKELIHOOD_SUM_COMPARATOR);

        // We need to return the best likelihood alleles in the original alternative allele index order.
        // This heap will keep track of that index order.
        final PriorityQueue<LikelihoodSum> mostLikelyAllelesHeapByIndex = new PriorityQueue<>(numOriginalAltAlleles, LIKELIHOOD_INDEX_COMPARATOR);

        for ( int i = 0; i < numAllelesToChoose; i++ ) {
            mostLikelyAllelesHeapByIndex.add(likelihoodSums[i]);
        }

        // guaranteed no to have been added at this point thanks for checking on whether reduction was
        // needed in the first place.
        if (nonRefAltAllelePresent) {
            mostLikelyAllelesHeapByIndex.add(likelihoodSums[nonRefAltAlleleIndex]);
        }

        final List<Allele> orderedBestAlleles = new ArrayList<>(numAllelesToChoose);

        while (!mostLikelyAllelesHeapByIndex.isEmpty()) {
            orderedBestAlleles.add(mostLikelyAllelesHeapByIndex.remove().allele);
        }

        return orderedBestAlleles;
    }

    protected GenotypesContext reduceScopeGenotypes(final VariantContext vc, final int defaultPloidy, final List<Allele> allelesToUse) {
        Utils.nonNull(vc, "vc is null");
        Utils.nonNull(allelesToUse, "allelesToUse is null");
        return subsetAlleles(vc,defaultPloidy,allelesToUse,false);
    }

    protected void reduceScopeCalculateLikelihoodSums(final VariantContext vc, final int defaultPloidy, final LikelihoodSum[] likelihoodSums) {
        Utils.nonNull(vc, "vc is null");
        Utils.nonNull(likelihoodSums, "likelihoodSums is null");

        final int numOriginalAltAlleles = likelihoodSums.length;
        final GenotypesContext genotypes = vc.getGenotypes();
        for ( final Genotype genotype : genotypes.iterateInSampleNameOrder() ) {
            if (!genotype.hasPL()) {
                continue;
            }
            final double[] gls = genotype.getLikelihoods().getAsVector();
            if (MathUtils.sum(gls) >= GATKVariantContextUtils.SUM_GL_THRESH_NOCALL) {
                continue;
            }

            final int PLindexOfBestGL = MathUtils.maxElementIndex(gls);

            final double bestToHomRefDiffGL = PLindexOfBestGL == PL_INDEX_OF_HOM_REF ? 0.0 : gls[PLindexOfBestGL] - gls[PL_INDEX_OF_HOM_REF];
            final int declaredPloidy = genotype.getPloidy();
            final int ploidy = declaredPloidy <= 0 ? defaultPloidy : declaredPloidy;

            final int[] acCount = getAlleleCountFromPLIndex(1 + numOriginalAltAlleles, ploidy, PLindexOfBestGL);
            // by convention, first count coming from getAlleleCountFromPLIndex comes from reference allele
            for (int k=1; k < acCount.length;k++) {
                if (acCount[k] > 0) {
                    likelihoodSums[k - 1].sum += acCount[k] * bestToHomRefDiffGL;
                }
            }
        }
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
    public GenotypesContext subsetAlleles(final VariantContext vc,
                                          final int defaultPloidy,
                                          final List<Allele> allelesToUse,
                                          final boolean assignGenotypes) {
        Utils.nonNull(vc, "vc is null");
        Utils.nonNull(allelesToUse, "allelesToUse is null");

        // the genotypes with PLs
        final GenotypesContext oldGTs = vc.getGenotypes();

        // samples
        final List<String> sampleIndices = oldGTs.getSampleNamesOrderedByName();

        // the new genotypes to create
        final GenotypesContext newGTs = GenotypesContext.create();

        // we need to determine which of the alternate alleles (and hence the likelihoods) to use and carry forward
        final int numOriginalAltAlleles = vc.getAlternateAlleles().size();
        final int numNewAltAlleles = allelesToUse.size() - 1;


        // create the new genotypes
        for ( int k = 0; k < oldGTs.size(); k++ ) {
            final Genotype g = oldGTs.get(sampleIndices.get(k));
            final int declaredPloidy = g.getPloidy();
            final int ploidy = declaredPloidy <= 0 ? defaultPloidy : declaredPloidy;
            if ( !g.hasLikelihoods() ) {
                newGTs.add(GenotypeBuilder.create(g.getSampleName(), GATKVariantContextUtils.noCallAlleles(ploidy)));
                continue;
            }

            // create the new likelihoods array from the alleles we are allowed to use
            final double[] originalLikelihoods = g.getLikelihoods().getAsVector();
            double[] newLikelihoods;

            // Optimization: if # of new alt alleles = 0 (pure ref call), keep original likelihoods so we skip normalization
            // and subsetting
            if ( numOriginalAltAlleles == numNewAltAlleles || numNewAltAlleles == 0) {
                newLikelihoods = originalLikelihoods;
            } else {
                newLikelihoods = subsetToAlleles(originalLikelihoods, ploidy, vc.getAlleles(), allelesToUse);

                // might need to re-normalize
                newLikelihoods = MathUtils.normalizeFromLog10(newLikelihoods, false, true);
            }

            // if there is no mass on the (new) likelihoods, then just no-call the sample
            if ( MathUtils.sum(newLikelihoods) > GATKVariantContextUtils.SUM_GL_THRESH_NOCALL ) {
                newGTs.add(GenotypeBuilder.create(g.getSampleName(), GATKVariantContextUtils.noCallAlleles(ploidy)));
            }
            else {
                final GenotypeBuilder gb = new GenotypeBuilder(g);

                if ( numNewAltAlleles == 0 ) {
                    gb.noPL();
                } else {
                    gb.PL(newLikelihoods);
                }

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


        int idx = 0;
        // First fill boolean array stating whether each original allele is present in new mapping
        final boolean [] allelePresent = new boolean[originalAlleles.size()];
        for ( final Allele allele : originalAlleles ) {
            allelePresent[idx++] = allelesToSubset.contains(allele);
        }


        // compute mapping from old idx to new idx
        // This might be needed in case new allele set is not ordered in the same way as old set
        // Example. Original alleles: {T*,C,G,A}. New alleles: {G,C}. Permutation key = [2,1]

        final int[] permutationKey = new int[allelesToSubset.size()];
        for (int k=0; k < allelesToSubset.size(); k++)
        // for each allele to subset, find corresponding index in original allele list
        {
            permutationKey[k] = originalAlleles.indexOf(allelesToSubset.get(k));
        }


        final SumIterator iterator = new SumIterator(originalAlleles.size(),numChromosomes);

        while (iterator.hasNext()) {
            // for each entry in logPL table, associated originally with allele count stored in vec[],
            // see if this allele count conformation will be present in new logPL table.
            // For entry to be present, elements in dimensions not present in requested allele list have to have count = 0
            final int[] pVec = iterator.getCurrentVector();
            final double pl = oldLikelihoods[iterator.getLinearIndex()];

            boolean keyPresent = true;
            for (int k=0; k < allelePresent.length; k++) {
                if (pVec[k] > 0 && !allelePresent[k]) {
                    keyPresent = false;
                }
            }

            if (keyPresent) {// skip to next entry in logPLs if this conformation is not present in subset

                final int[] newCount = new int[allelesToSubset.size()];

                // map from old allele mapping count to new allele mapping
                // In pseudo-Matlab notation: newCount = vec[permutationKey] for permutationKey vector
                for (idx = 0; idx < newCount.length; idx++) {
                    newCount[idx] = pVec[permutationKey[idx]];
                }

                // get corresponding index from new count
                final int outputIdx = getLinearIndex(newCount, allelesToSubset.size(), numChromosomes);
                newPLs[outputIdx] = pl;
            }
            iterator.next();
        }

        return  newPLs;
    }


    private static int getLinearIndex(final int[] vectorIdx, final int numAlleles, final int ploidy) {

        if (ploidy <= 0) {
            return 0;
        }

        int linearIdx = 0;
        int cumSum = ploidy;
        for (int k=numAlleles-1;k>=1; k--) {
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

    /*
* a cache of the PL ivector sizes as a function of # of alleles and pool sizes
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
                if (numAlleles == 1) {
                    cache[numAlleles][ploidy] = 1;
                } else if (ploidy == 1) {
                    cache[numAlleles][ploidy] = numAlleles;
                } else {
                    int acc =0;
                    for (int k=0; k <= ploidy; k++ ) {
                        acc += cache[numAlleles - 1][ploidy - k];
                    }

                    cache[numAlleles][ploidy] = acc;
                }
            }
        }
        return cache;
    }

    /**
     * Crucial inner class that handles addressing elements of pool likelihoods. We store likelihoods as a map
     * of form int[] -> double (to be more precise, IntArrayWrapper -> Double).
     * For a given ploidy (chromosome count) and number of alleles, we need a form to iterate deterministically
     * across all possible allele conformations.
     * Problem equivalent to listing in determistic order all possible ways in which N integers will sum to P,
     * where N is number of alleles and P is number of chromosomes.
     * There's an option to list all integers so that sum will be UP to P.
     * For example, with P=2,N=2, restrictSumTo = 2 iterator will produce
     * [2 0 ] [1 1] [ 0 2]
     *
     *
     */
    private static final class SumIterator {
        private int[] currentState;
        private final int[] finalState;
        private final int restrictSumTo;
        private final int dim;
        private boolean hasNext;
        private int linearIndex;
        private int currentSum;

        /**
         * Default constructor. Typical use case: restrictSumTo = -1 if there's no sum restriction, or will generate int[]
         * vectors so that all add to this value.
         *
         * @param finalState                    End state - typically we should set value to (P,P,P,...)
         * @param restrictSumTo                 See above
         */
        public SumIterator(final int[] finalState,final int restrictSumTo) {
            this.finalState = finalState;
            this.dim = finalState.length;
            this.restrictSumTo = restrictSumTo;
            currentState = new int[dim];
            reset();

        }

        /**
         * Shortcut constructor for common use case: iterator will produce
         * all vectors of length numAlleles whose sum = numChromosomes
         * @param numAlleles              Number of alleles
         * @param numChromosomes          Ploidy
         */
        public SumIterator(final int numAlleles, final int numChromosomes) {
            this(getInitialStateVector(numAlleles, numChromosomes), numChromosomes);
        }


        private static int[] getInitialStateVector(final int nAlleles, final int numChromosomes) {
            final int[] initialState = new int[nAlleles];
            Arrays.fill(initialState, numChromosomes);
            return initialState;
        }

        public void next() {
            final int initialDim = (restrictSumTo > 0)?1:0;
            hasNext = next(finalState, initialDim);
            if (hasNext) {
                linearIndex++;
            }
        }

        private boolean next(final int[] finalState, final int initialDim) {
            boolean hasNextState = false;
            for (int currentDim=initialDim; currentDim < finalState.length; currentDim++) {
                final int x = currentState[currentDim]+1;

                if (x > finalState[currentDim] || (currentSum >= restrictSumTo && initialDim > 0)) {
                    // update vector sum, and reset position
                    currentSum -= currentState[currentDim];
                    currentState[currentDim] = 0;
                    if (currentDim >= dim-1) {
                        hasNextState = false;
                        break;
                    }
                }
                else {
                    currentState[currentDim] = x;
                    hasNextState = true;
                    currentSum++;
                    break;
                }
            }
            if (initialDim > 0) {
                currentState[0] = restrictSumTo - currentSum;
            }
            return hasNextState;
        }

        public void reset() {
            Arrays.fill(currentState, 0);
            if (restrictSumTo > 0) {
                currentState[0] = restrictSumTo;
            }
            hasNext = true;
            linearIndex = 0;
            currentSum = 0;
        }
        public int[] getCurrentVector() {
            return currentState;
        }

        public int getLinearIndex() {
            return linearIndex;
        }

        public boolean hasNext() {
            return hasNext;
        }
    }


    /**
     * Assign genotypes (GTs) to the samples in the Variant Context greedily based on the PLs
     *
     * @param newLikelihoods       the PL array
     * @param allelesToUse         the list of alleles to choose from (corresponding to the PLs)
     * @param numChromosomes        Number of chromosomes per pool
     */
    private static void assignGenotype(final GenotypeBuilder gb,
                                       final double[] newLikelihoods,
                                       final List<Allele> allelesToUse,
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

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
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * AFCalculator in which genotype posterior probabilities are equal to their likelihoods, normalized.
 *
 * For now I'm just using this as a simple implementation to hook up to the GenotypingEngine in order to delete all
 * the other implementations and their associated classes like StateTracker and ExactACSet.
 */
public class FlatPriorAFCalculator extends AFCalculator {

    private static final int INDEX_OF_HOM_REF = 0;
    @Override
    protected AFCalculationResult computeLog10PNonRef(final VariantContext vc, final int defaultPloidy) {
        Utils.nonNull(vc, "vc is null");

        final int[] alleleCountsOfMLE = new int[vc.getNAlleles()];
        final List<Allele> allelesUsedInGenotyping = vc.getAlleles();
        final double log10PosteriorOfAFEq0 = -1.0;    //TODO: placeholder
        final Map<Allele, Double> log10pRefByAllele = null;

        final GenotypesContext GLs = vc.getGenotypes();
        final int numAlleles = vc.getNAlleles();

        double log10ProbAllHomRef = 0.0;
        double[] log10ProbAlleleAbsent = new double[vc.getNAlleles()];
        for (final Genotype genotype : GLs.iterateInSampleNameOrder()) {
            if (!genotype.hasPL()) {
                continue;
            }
            final double[] gls = genotype.getLikelihoods().getAsVector();
            final double[] genotypePosteriors = MathUtils.normalizeFromLog10(gls, false);   // real, not log
            log10ProbAllHomRef += genotypePosteriors[0];   // log10PAllHomRef = product of individual samples being hom ref


        }


        return new AFCalculationResult(alleleCountsOfMLE, allelesUsedInGenotyping, log10PosteriorOfAFEq0, log10pRefByAllele);
    }


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

        final List<Allele> outputAlleles = Stream.concat(Stream.of(vc.getReference()), outputAltAlleles.stream()).collect(Collectors.toList());
        final GenotypesContext outputGenotypes = subsetAlleles(vc,defaultPloidy,outputAlleles,false);
        return new VariantContextBuilder(vc).alleles(outputAlleles).genotypes(outputGenotypes).make();
    }

    /**
     * Returns a the new set of alleles to use.
     * @param vc target variant context.
     * @param numAllelesToChoose number of alleles to keep.
     * @return the list of alternative alleles to keep.
     */
    protected List<Allele> reduceScopeAlleles(final VariantContext vc, final int defaultPloidy, final int numAllelesToChoose) {
        Utils.nonNull(vc, "vc is null");
        final int nonRefAltAllele = GATKVariantContextUtils.indexOfAltAllele(vc, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE, false);
        final int[] properAltAlleles = IntStream.range(1, vc.getNAlleles()).filter(n -> n != nonRefAltAllele).toArray();

        // Avoid pointless allele reduction:
        if (numAllelesToChoose >= properAltAlleles.length) {
            return vc.getAlternateAlleles();
        }

        final double[] likelihoodSums = new double[vc.getNAlleles()];
        for ( final Genotype genotype : vc.getGenotypes().iterateInSampleNameOrder() ) {
            final double[] gls = genotype.getLikelihoods().getAsVector();
            if (gls == null || GATKVariantContextUtils.likelihoodsAreUninformative(gls)) {
                continue;
            }

            final int indexOfMostLikelyGenotype = MathUtils.maxElementIndex(gls);
            final double GLDiffBetweenRefAndBest = gls[indexOfMostLikelyGenotype] - gls[INDEX_OF_HOM_REF];
            final int ploidy = genotype.getPloidy() > 0 ? genotype.getPloidy() : defaultPloidy;

            final int[] alleleCounts = getAlleleCountsFromIndex(vc.getNAlleles(), ploidy, indexOfMostLikelyGenotype);
            for (int allele = 1; allele < alleleCounts.length; allele++) {
                likelihoodSums[allele] += alleleCounts[allele] * GLDiffBetweenRefAndBest;
            }
        }

        final List<Double> properAltAlleleLikelihoodSums = Arrays.stream(properAltAlleles)
                .mapToObj(n -> likelihoodSums[n]).collect(Collectors.toList());
        Collections.sort(properAltAlleleLikelihoodSums, Collections.reverseOrder());
        final double likelihoodSumThreshold = properAltAlleleLikelihoodSums.get(numAllelesToChoose);
        return IntStream.range(1, vc.getNAlleles())
                .filter(n -> n == nonRefAltAllele || likelihoodSums[n] > likelihoodSumThreshold)
                .map(n -> n-1)  //go from allele index to alt allele index
                .mapToObj(vc::getAlternateAllele).collect(Collectors.toList());
    }

    /**
     * Given a scalar index, what's the alelle count conformation corresponding to it?
     * @param nAlleles                    Number of alleles
     * @param numChromosomes              Ploidy
     * @param index                     Index to query
     * @return                            Allele count conformation, according to iteration order from GenotypeIterator
     */
    private static int[] getAlleleCountsFromIndex(final int nAlleles, final int numChromosomes, final int index) {
        final GenotypeLikelihoodCalculator calculator = new GenotypeLikelihoodCalculators().getInstance(numChromosomes, nAlleles);
        final GenotypeAlleleCounts alleleCounts = calculator.genotypeAlleleCountsAt(index);
        return alleleCounts.alleleCountsByIndex(nAlleles - 1);
    }

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
                gb.PL(newLikelihoods);

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
     * Given set of alleles with corresponding vector of genotype likelihoods (in the canonical order)
     * output likelihoods in the canonical order for a subset of the alleles.
     *
     * @param oldLikelihoods        Vector of PL's corresponding to original alleles
     * @param numChromosomes        Ploidy (number of chromosomes describing PL's)
     * @param originalAlleles       List of original alleles
     * @param newAlleles            Alleles to subset
     * @return                      Vector of new PL's, ordered accorrding to GenotypeIterator's ordering
     */
    private static double[] subsetToAlleles(final double[] oldLikelihoods, final int numChromosomes,
                                            final List<Allele> originalAlleles, final List<Allele> newAlleles) {
        final int newPLSize = numGenotypes(newAlleles.size(), numChromosomes);
        final double[] newPLs = new double[newPLSize];

        // First fill boolean array stating whether each original allele is present in new mapping
        final List<Boolean> allelePresent = originalAlleles.stream().map(a -> newAlleles.contains(a)).collect(Collectors.toList());

        // compute mapping from old idx to new idx eg Original: {T*,C,G,A}, New: {G,C}. Permutation key = [2,1]
        final int[] oldIndicesOfNewAlleles = newAlleles.stream().mapToInt(a -> originalAlleles.indexOf(a)).toArray();

        // iterate over all genotypes containing the old alleles, and for each genotype that contains only new alleles
        // copy its GL to the appropriate index of the new GL array
        for (GenotypeIterator iterator = new GenotypeIterator(originalAlleles.size(),numChromosomes); iterator.hasNext(); iterator.next()) {
            // skip this genotype if it contains any alleles not in the new allele subset
            final int[] alleleCounts = iterator.getCurrentVector();
            if (IntStream.range(0, originalAlleles.size()).anyMatch(k -> alleleCounts[k] > 0 && !allelePresent.get(k))) {
                continue;
            }

            final int[] newAlleleCounts = IntStream.range(0, newAlleles.size())
                    .map(newAllele -> alleleCounts[oldIndicesOfNewAlleles[newAllele]]).toArray();

            // get corresponding index from new count
            final int outputIdx = getGenotypeIndex(newAlleleCounts, newAlleles.size(), numChromosomes);
            newPLs[outputIdx] = oldLikelihoods[iterator.getLinearIndex()];
        }

        return  newPLs;
    }

    @VisibleForTesting
    static int numGenotypes(final int numAlleles, final int ploidy) { return GL_ARRAY_SIZES[numAlleles][ploidy]; }

    private static final int MAX_NUM_ALLELES_TO_CACHE = 20;
    private static final int MAX_NUM_SAMPLES_PER_POOL = 1000;

    //Note: this is shared state but it's not modified at runtime
    private static final int[][] GL_ARRAY_SIZES = fillGLVectorSizeCache(MAX_NUM_ALLELES_TO_CACHE, 2*MAX_NUM_SAMPLES_PER_POOL);

    private static int[][] fillGLVectorSizeCache(final int maxAlleles, final int maxPloidy) {
        final int[][] cache = new int[maxAlleles][maxPloidy];
        for (int numAlleles=1; numAlleles < maxAlleles; numAlleles++) {
            for (int ploidy=0; ploidy < maxPloidy; ploidy++) {
                cache[numAlleles][ploidy] = numAlleles == 1 ? 1 : Arrays.stream(cache[numAlleles - 1], 0, ploidy + 1).sum();
            }
        }
        return cache;
    }

    // given vector of allele counts (i.e. number of copies of each allele a genotype contains)
    // find the genotype's index in the canonical ordering.
    // To do this total the number of genotypes that equal this one past a certain allele and have lower index below.
    private static int getGenotypeIndex(final int[] alleleCounts, final int numAlleles, final int ploidy) {
        if (ploidy == 0) {
            return 0;
        }
        int numGenotypesWithLowerIndex = 0;
        int remainingPloidy = ploidy;
        for (int allele = numAlleles - 1; allele > 0; allele--) {
            final int countAtThisAllele = alleleCounts[allele];
            for (int smallerCountAtThisAllele=0; smallerCountAtThisAllele < countAtThisAllele; smallerCountAtThisAllele++) {
                numGenotypesWithLowerIndex += numGenotypes(allele, remainingPloidy - smallerCountAtThisAllele);
            }
            remainingPloidy -= countAtThisAllele;
        }
        return numGenotypesWithLowerIndex;
    }


    //TODO: first, it's kind of ugly how this is a side effect of subsetting and not its own method
    //TODO: second, there should be a more principled way of assigning GTs that is not simply based on likelihoods
    /**
     * Assign genotypes (GTs) to the samples in the Variant Context greedily based on the PLs
     *
     * @param newLikelihoods       the PL array
     * @param allelesToUse         the list of alleles to choose from (corresponding to the PLs)
     * @param numChromosomes        Number of chromosomes per pool
     */
    private static void assignGenotype(final GenotypeBuilder gb, final double[] newLikelihoods, final List<Allele> allelesToUse,
                                       final int numChromosomes) {
        final int mostLikelyGenotypeIndex = MathUtils.maxElementIndex(newLikelihoods);
        final GenotypeLikelihoodCalculator calculator = new GenotypeLikelihoodCalculators().getInstance(numChromosomes, allelesToUse.size());
        final GenotypeAlleleCounts alleleCounts = calculator.genotypeAlleleCountsAt(mostLikelyGenotypeIndex);
        gb.alleles(alleleCounts.asAlleleList(allelesToUse));

        // remove PLs if necessary
        if (newLikelihoods.length > MAX_LENGTH_FOR_POOL_PL_LOGGING) {
            gb.noPL();
        }
        if ( allelesToUse.size() > 1 ) {
            gb.log10PError(GenotypeLikelihoods.getGQLog10FromLikelihoods(mostLikelyGenotypeIndex, newLikelihoods));
        }
    }
}

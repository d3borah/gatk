package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.*;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVCFConstants;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Created by davidben on 6/8/16.
 */
public class AlleleSubsettingUtils {
    private AlleleSubsettingUtils() {}  //prevent instantiation
    private static final int PL_INDEX_OF_HOM_REF = 0;
    private static final int MAX_LENGTH_FOR_PL_LOGGING = 100; // if PL vectors longer than this # of elements, don't log them

    /**
     * From a given variant context, extract a given subset of alleles, and update genotype context accordingly,
     * including updating the PL's, and assign genotypes accordingly
     * @param vc                                variant context with alleles and genotype likelihoods
     * @param defaultPloidy                     ploidy to use when {@code vc} is missing information for a sample
     * @param allelesToKeep                     subset of alleles in {@code vc} to keep
     * @param method                            whether to assign hard genotypes or leave as no-call
     * @return                                  GenotypesContext with new PLs
     */
    public static GenotypesContext subsetAlleles(final VariantContext vc, final int defaultPloidy, final boolean isDiploid,
                                                              final List<Allele> allelesToKeep, final GenotypeAssignmentMethod method) {
        return isDiploid ? subsetDiploidAlleles(vc, allelesToKeep, method) :
                subsetGeneralPloidyAlleles(vc, defaultPloidy, allelesToKeep, method);
    }

    /**
     * Restrict a variant context to its "best" alleles.
     * @param vc                    variant context with alleles and genotype likelihoods
     * @param defaultPloidy         ploidy to assume in case that {@code vc} does not contain that information
     * @param numAltAllelesToKeep
     * @param isDiploid
     * @return
     */
    public static VariantContext reduceNumberOfAlleles(final VariantContext vc, final int defaultPloidy,
                                                       final int numAltAllelesToKeep, final boolean isDiploid,
                                                       final GenotypeAssignmentMethod method) {
        Utils.nonNull(vc, "vc is null");
        final List<Allele> originalAltAlleles = vc.getAlternateAlleles();
        if (numAltAllelesToKeep >= originalAltAlleles.size()) {
            return vc;
        }
        final List<Allele> outputAlleles = calculateAllelesToKeep(vc, defaultPloidy, numAltAllelesToKeep);
        final GenotypesContext reducedGenotypes = subsetAlleles(vc, defaultPloidy, isDiploid, outputAlleles, method);
        return new VariantContextBuilder(vc).alleles(outputAlleles).genotypes(reducedGenotypes).make();
    }

    /**
     * Updates the PLs and AD of the Genotypes in the newly selected VariantContext to reflect the fact that some alleles
     * from the original VariantContext are no longer present.
     *
     * @param selectedVC  the selected (new) VariantContext
     * @param originalVC  the original VariantContext
     * @return a new non-null GenotypesContext
     */
    public static GenotypesContext updatePLsAndAD(final VariantContext selectedVC, final VariantContext originalVC) {
        final int numNewAlleles = selectedVC.getNAlleles();
        final int numOriginalAlleles = originalVC.getNAlleles();
        Utils.validateArg( numNewAlleles <= numOriginalAlleles, "Attempting to fix PLs and AD from what appears to be a *combined* VCF and not a selected one");

        final GenotypesContext oldGs = selectedVC.getGenotypes();
        return numNewAlleles == numOriginalAlleles ? oldGs :
                createDiploidGenotypesWithSubsettedLikelihoods(oldGs, originalVC, selectedVC.getAlleles(), GenotypeAssignmentMethod.DO_NOT_ASSIGN_GENOTYPES);
    }

    /**
     * Add the genotype call (GT) field to GenotypeBuilder using the requested {@link GenotypeAssignmentMethod}
     *
     * @param originalGT the original genotype calls, cannot be null
     * @param gb the builder where we should put our newly called alleles, cannot be null
     * @param assignmentMethod the method to use to do the assignment, cannot be null
     * @param newLikelihoods a vector of likelihoods to use if the method requires PLs, should be log10 likelihoods, cannot be null
     * @param allelesToUse the alleles we are using for our subsetting
     */
    public static void updateGenotypeCallAfterSubsetting(final List<Allele> originalGT,
                                                         final GenotypeBuilder gb,
                                                         final GenotypeAssignmentMethod assignmentMethod,
                                                         final double[] newLikelihoods,
                                                         final List<Allele> allelesToUse) {
        switch ( assignmentMethod ) {
            case DO_NOT_ASSIGN_GENOTYPES:
                break;
            case SET_TO_NO_CALL:
                gb.alleles(GATKVariantContextUtils.noCallAlleles(2)); // @TODO assumes diploid
                gb.noGQ();
                break;
            case USE_PLS_TO_ASSIGN:
                if ( newLikelihoods == null || !GATKVariantContextUtils.isInformative(newLikelihoods) ) {
                    // if there is no mass on the (new) likelihoods, then just no-call the sample
                    gb.alleles(GATKVariantContextUtils.noCallAlleles(2)); // @TODO assumes diploid
                    gb.noGQ();
                } else {
                    // find the genotype with maximum likelihoods
                    final int PLindex = MathUtils.maxElementIndex(newLikelihoods);
                    GenotypeLikelihoods.GenotypeLikelihoodsAllelePair alleles = GenotypeLikelihoods.getAllelePair(PLindex);
                    gb.alleles(Arrays.asList(allelesToUse.get(alleles.alleleIndex1), allelesToUse.get(alleles.alleleIndex2)));
                    gb.log10PError(GenotypeLikelihoods.getGQLog10FromLikelihoods(PLindex, newLikelihoods));
                }
                break;
            case BEST_MATCH_TO_ORIGINAL:
                final List<Allele> best = new LinkedList<>();
                final Allele ref = allelesToUse.get(0); // WARNING -- should be checked in input argument
                for ( final Allele originalAllele : originalGT ) {
                    best.add(allelesToUse.contains(originalAllele) ? originalAllele : ref);
                }
                gb.noGQ();
                gb.noPL();
                gb.alleles(best);
                break;
        }
    }

    /**
     * From a given variant context, extract a given subset of alleles, and update genotype context accordingly,
     * including updating the PL's, and assign genotypes accordingly
     * @param vc                                variant context with alleles and genotype likelihoods
     * @param defaultPloidy                     ploidy to assume in case that {@code vc} does not contain that information
     *                                          for a sample.
     * @param allelesToKeep                      alleles to subset
     * @param method                            whether to assign hard genotypes or leave as no-call
     * @return                                  GenotypesContext with new PLs
     */
    private static GenotypesContext subsetGeneralPloidyAlleles(final VariantContext vc, final int defaultPloidy,
                                                               final List<Allele> allelesToKeep, final GenotypeAssignmentMethod method) {
        Utils.nonNull(vc, "vc is null");
        Utils.nonNull(allelesToKeep, "allelesToKeep is null");
        if (allelesToKeep.size() == 1) {
            return GATKVariantContextUtils.subsetToRefOnly(vc, defaultPloidy);
        }
        final GenotypesContext newGTs = GenotypesContext.create();

        // create the new genotypes
        for ( final Genotype g : vc.getGenotypes() ) {
            final int ploidy = g.getPloidy() > 0 ? g.getPloidy() : defaultPloidy;
            if ( !g.hasLikelihoods() ) {
                newGTs.add(GenotypeBuilder.create(g.getSampleName(), GATKVariantContextUtils.noCallAlleles(ploidy)));
                continue;
            }

            final double[] originalLikelihoods = g.getLikelihoods().getAsVector();
            final double[] newLikelihoods = subsettedGenotypeLikelihoods(originalLikelihoods, ploidy, vc.getAlleles(), allelesToKeep);

            if (!GATKVariantContextUtils.isInformative(newLikelihoods) ) {
                newGTs.add(GenotypeBuilder.create(g.getSampleName(), GATKVariantContextUtils.noCallAlleles(ploidy)));
                continue;
            }

            final GenotypeBuilder gb = new GenotypeBuilder(g);
            gb.PL(newLikelihoods);


            // if we weren't asked to assign a genotype, then just no-call the sample
            if (method == GenotypeAssignmentMethod.DO_NOT_ASSIGN_GENOTYPES || method == GenotypeAssignmentMethod.SET_TO_NO_CALL) {
                gb.alleles(GATKVariantContextUtils.noCallAlleles(ploidy));
            } else {
                assignGenotypes(gb, newLikelihoods, allelesToKeep, ploidy);
            }
            newGTs.add(gb.make());
        }
        return newGTs;
    }

    /**
     * subset the Variant Context to the specific set of alleles passed in (pruning the PLs appropriately)
     *
     * @param vc                 variant context with genotype likelihoods
     * @param allelesToKeep       which alleles from the vc are okay to use; *** must be in the same relative order as those in the original VC ***
     * @param assignGenotypes    assignment strategy for the (subsetted) PLs
     * @return a new non-null GenotypesContext
     */
    private static GenotypesContext subsetDiploidAlleles(final VariantContext vc,
                                                        final List<Allele> allelesToKeep,
                                                        final GenotypeAssignmentMethod assignGenotypes) {
        Utils.nonNull(allelesToKeep, "allelesToKeep is null");
        Utils.validateArg(allelesToKeep.get(0).isReference(), "First allele must be the reference allele");
        if(allelesToKeep.size() == 1) {
            return GATKVariantContextUtils.subsetToRefOnly(vc, 2);
        }

        // optimization: if no input genotypes, just exit
        return vc.getGenotypes().isEmpty() ? GenotypesContext.create() :
                createDiploidGenotypesWithSubsettedLikelihoods(vc.getGenotypes(), vc, allelesToKeep, assignGenotypes);
    }

    /**
     * Assign genotypes (GTs) to the samples in the Variant Context greedily based on the PLs
     *
     * @param newLikelihoods       the PL array
     * @param allelesToKeep         the list of alleles to choose from (corresponding to the PLs)
     * @param ploidy               ploidy
     */
    private static void assignGenotypes(final GenotypeBuilder gb, final double[] newLikelihoods,
                                        final List<Allele> allelesToKeep, final int ploidy) {
        final int numNewAltAlleles = allelesToKeep.size() - 1;

        // find the genotype with maximum likelihoods
        final int indexOfLikeliestGenotype = MathUtils.maxElementIndex(newLikelihoods);
        final GenotypeLikelihoodCalculator calculator = new GenotypeLikelihoodCalculators().getInstance(ploidy, allelesToKeep.size());
        final GenotypeAlleleCounts alleleCounts = calculator.genotypeAlleleCountsAt(indexOfLikeliestGenotype);

        gb.alleles(alleleCounts.asAlleleList(allelesToKeep));

        // remove PLs if necessary in case of giant ploidy due to pooling
        if (newLikelihoods.length > MAX_LENGTH_FOR_PL_LOGGING) {
            gb.noPL();
        }

        if ( numNewAltAlleles > 0 ) {
            gb.log10PError(GenotypeLikelihoods.getGQLog10FromLikelihoods(indexOfLikeliestGenotype, newLikelihoods));
        }
    }

    /**
     * Given set of alleles with corresponding vector of likelihoods, subset likelihoods to new set of alleles
     *
     * @param oldLikelihoods        Vector of PL's corresponding to original alleles
     * @param ploidy                Ploidy (number of chromosomes describing PL's)
     * @param originalAlleles       List of original alleles
     * @param newAlleles            New alleles -- must be a subset of {@code originalAlleles}
     * @return                      Vector of new PL's, ordered accorrding to GenotpyeIterator's ordering
     */
    private static double[] subsettedGenotypeLikelihoods(final double[] oldLikelihoods, final int ploidy,
                                                         final List<Allele> originalAlleles, final List<Allele> newAlleles) {
        //Since this is a private method, we will *assume* that newAlleles is a subset of originalAlleles
        // Then we have shortcuts based on two trivial cases:
        if (newAlleles.size() == originalAlleles.size()) {
            return oldLikelihoods;  //no subsetting
        } else if (newAlleles.size() == 1) {
            return new double[] { oldLikelihoods[0]};   //hom ref is the only genotype
        }

        final int newPLSize = GLVectorSizeCache.getNumLikelihoodElements(newAlleles.size(), ploidy);
        final double[] newPLs = new double[newPLSize];

        // fill boolean array stating whether each original allele is present in new mapping
        final List<Boolean> allelePresent = originalAlleles.stream()
                .map(a -> newAlleles.contains(a)).collect(Collectors.toList());

        // compute old indices of new alleles in case new allele set is not ordered in the same way as old set
        // Example. Original alleles: {T*,C,G,A}. New alleles: {G,C}. oldIndexTable = [2,1]
        final int[] oldIndexTable = newAlleles.stream().mapToInt(a -> originalAlleles.indexOf(a)).toArray();
        for(final GenotpyeIterator iterator = new GenotpyeIterator(originalAlleles.size(),ploidy); iterator.hasNext(); iterator.next()) {
            final int[] oldAlleleCounts = iterator.getAlleleCounts();
            final boolean containsOnlyNewAlleles = IntStream.range(0, originalAlleles.size())
                    .filter(k -> oldAlleleCounts[k] > 0)
                    .allMatch(allelePresent::get);

            if (containsOnlyNewAlleles) {
                final int[] newAlleleCounts = IntStream.range(0, newAlleles.size())
                        .map(a -> oldAlleleCounts[oldIndexTable[a]]).toArray();
                final int newIndex = GenotpyeIterator.getLinearIndex(newAlleleCounts, newAlleles.size(), ploidy);
                newPLs[newIndex] = oldLikelihoods[iterator.getLinearIndex()];
            }
        }

        return  newPLs;
    }

    /**
     * Create the new GenotypesContext with the subsetted PLs and ADs
     *
     * @param originalGs               the original GenotypesContext
     * @param vc                       the original VariantContext
     * @param allelesToUse             the actual alleles to use with the new Genotypes
     * @param assignGenotypes          assignment strategy for the (subsetted) PLs
     * @return a new non-null GenotypesContext
     */
    private static GenotypesContext createDiploidGenotypesWithSubsettedLikelihoods(final GenotypesContext originalGs,
                                                                                  final VariantContext vc,
                                                                                  final List<Allele> allelesToUse,
                                                                                  final GenotypeAssignmentMethod assignGenotypes) {
        final List<Integer> likelihoodIndexesToUse = determineDiploidLikelihoodIndexesToUse(vc, allelesToUse);
        final GenotypesContext newGTs = GenotypesContext.create(originalGs.size());

        // make sure we are seeing the expected number of likelihoods per sample
        final int expectedNumLikelihoods = GenotypeLikelihoods.numLikelihoods(vc.getNAlleles(), 2);

        // create the new genotypes
        for (final Genotype g : originalGs.iterateInSampleNameOrder()) {
            final GenotypeBuilder gb = new GenotypeBuilder(g);

            // create the new likelihoods array from the alleles we are allowed to use
            double[] newLikelihoods = null;
            if (g.hasLikelihoods()) {
                final double[] originalLikelihoods = g.getLikelihoods().getAsVector();
                newLikelihoods = originalLikelihoods.length == expectedNumLikelihoods ?
                        MathUtils.normalizeFromLog10(likelihoodIndexesToUse.stream()
                                .mapToDouble(idx -> originalLikelihoods[idx]).toArray(), false, true) : null;
            }
            if ( newLikelihoods == null || !GATKVariantContextUtils.isInformative(newLikelihoods) )
                gb.noPL();
            else
                gb.PL(newLikelihoods);

            updateGenotypeCallAfterSubsetting(g.getAlleles(), gb, assignGenotypes, newLikelihoods, allelesToUse);
            newGTs.add(gb.make());
        }

        return subsetADToKeptAlleles(newGTs, vc, allelesToUse);
    }

    /**
     * Figure out which likelihood indexes to use for a selected down set of alleles
     *
     * For non-diploid samples a fancier algorithm is required.
     *
     * @param vc        the original VariantContext
     * @param allelesToKeep      the subset of alleles to use
     * @return a list of PL indexes to use or null if none
     */
    private static List<Integer> determineDiploidLikelihoodIndexesToUse(final VariantContext vc, final List<Allele> allelesToKeep) {
        final List<Boolean> alleleIndexesToUse = indicesOfAllelesToKeep(vc, allelesToKeep);

        final int numLikelihoods = GenotypeLikelihoods.numLikelihoods(vc.getNAlleles(), GATKVariantContextUtils.DEFAULT_PLOIDY);
        return IntStream.range(0, numLikelihoods).filter(idx -> {
                    final GenotypeLikelihoods.GenotypeLikelihoodsAllelePair alleles = GenotypeLikelihoods.getAllelePair(idx);
                    return alleleIndexesToUse.get(alleles.alleleIndex1) && alleleIndexesToUse.get(alleles.alleleIndex2);})
                .boxed().collect(Collectors.toList());
    }

    /**
     * Fix the AD for the GenotypesContext of a VariantContext that has been subset
     *
     * @param originalGenotypes       the original GenotypesContext
     * @param originalVC       the original VariantContext
     * @param allelesToKeep     the new (sub)set of alleles to use
     * @return a new non-null GenotypesContext
     */
    private static GenotypesContext subsetADToKeptAlleles(final GenotypesContext originalGenotypes, final VariantContext originalVC,
                                                          final List<Allele> allelesToKeep) {
        final List<Boolean> alleleIndexesToUse = indicesOfAllelesToKeep(originalVC, allelesToKeep);

        final GenotypesContext newGTs = GenotypesContext.create(originalGenotypes.size());
        for (final Genotype g : originalGenotypes.iterateInSampleNameOrder()) {
            if(!g.hasAD()) {
                newGTs.add(g);
            } else {
                final GenotypeBuilder builder = new GenotypeBuilder(g);

                final int[] oldAD = g.getAD();
                if ( oldAD.length != alleleIndexesToUse.size() ) {
                    builder.noAD();
                } else {
                    final int[] newAD = IntStream.range(0, oldAD.length)
                            .filter(alleleIndexesToUse::get)
                            .map(i -> oldAD[i]).toArray();
                    builder.AD(newAD);
                }
                newGTs.add(builder.make());
            }
        }

        return newGTs;
    }

    /**
     * Given an original VariantContext and a list of alleles from that VC to keep,
     * returns a bitset representing which allele indexes should be kept
     *
     * @param originalVC      the original VC
     * @param allelesToKeep   the list of alleles to keep
     * @return                a non-null List
     */
    private static List<Boolean> indicesOfAllelesToKeep(final VariantContext originalVC, final List<Allele> allelesToKeep) {
        return originalVC.getAlleles().stream()
                .map(a -> a.isReference() || allelesToKeep.contains(a)).collect(Collectors.toList());
    }

    /**
     * Returns the new set of alleles to use.
     * @param vc target variant context.
     * @param numAltAllelesToKeep number of alleles to keep.
     * @return the list of alleles to keep, including the reference
     */
    private static List<Allele> calculateAllelesToKeep(final VariantContext vc, final int defaultPloidy,
                                                       final int numAltAllelesToKeep) {
        Utils.nonNull(vc, "vc is null");
        final int nonRefAltAlleleIndex = GATKVariantContextUtils.indexOfAltAllele(vc, GATKVCFConstants.NON_REF_SYMBOLIC_ALLELE, false);
        final int[] properAltAlleleIndices = IntStream.range(1, vc.getNAlleles()).filter(n -> n != nonRefAltAlleleIndex).toArray();

        if (numAltAllelesToKeep >= properAltAlleleIndices.length) {
            return vc.getAlleles();
        }

        final double[] likelihoodSums = new double[vc.getNAlleles()];
        for ( final Genotype genotype : vc.getGenotypes().iterateInSampleNameOrder() ) {
            final double[] gls = genotype.getLikelihoods().getAsVector();
            if (gls == null || !GATKVariantContextUtils.isInformative(gls)) {
                continue;
            }
            final int indexOfMostLikelyGenotype = MathUtils.maxElementIndex(gls);
            final double GLDiffBetweenRefAndBest = gls[indexOfMostLikelyGenotype] - gls[PL_INDEX_OF_HOM_REF];
            final int ploidy = genotype.getPloidy() > 0 ? genotype.getPloidy() : defaultPloidy;

            final int[] alleleCounts = new GenotypeLikelihoodCalculators()
                    .getInstance(ploidy, vc.getNAlleles()).genotypeAlleleCountsAt(indexOfMostLikelyGenotype)
                    .alleleCountsByIndex(vc.getNAlleles() - 1);

            for (int allele = 1; allele < alleleCounts.length; allele++) {
                likelihoodSums[allele] += alleleCounts[allele] * GLDiffBetweenRefAndBest;
            }
        }

        final List<Double> properAltAlleleLikelihoodSums = Arrays.stream(properAltAlleleIndices)
                .mapToObj(n -> likelihoodSums[n]).collect(Collectors.toList());
        Collections.sort(properAltAlleleLikelihoodSums, Collections.reverseOrder());
        final double likelihoodSumThreshold = properAltAlleleLikelihoodSums.get(numAltAllelesToKeep);
        return IntStream.range(0, vc.getNAlleles()) //keep ref, non-ref, and alts above the threshold
                .filter(n -> n == 0 || n == nonRefAltAlleleIndex || likelihoodSums[n] > likelihoodSumThreshold)
                .mapToObj(n -> vc.getAlternateAllele(n-1)).collect(Collectors.toList());
    }
}

package org.broadinstitute.hellbender.tools.walkers.genotyper;

import htsjdk.variant.variantcontext.Allele;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;

/**
 * Encapsulates the data use to make the genotype calls.
 */
public final class GenotypingData<A extends Allele> implements SampleList, AlleleList<A> {

    private final PloidyList ploidyList;

    private final ReadLikelihoods<A> likelihoods;

    /**
     * Constructs a new genotyping-data collection providing the ploidy model to apply to the input model
     * and the read-likelihoods collection.
     *
     *
     * @param ploidyList the ploidy model.
     * @param likelihoods the read-likelihoods collection.
     *
     * @throws IllegalArgumentException if either {@code ploidyList} or {@code likelihoods} is {@code null},
     *   or they are not compatible in terms of the samples they contain; their lists must match.
     */
    public GenotypingData(final PloidyList ploidyList, final ReadLikelihoods<A> likelihoods) {
        Utils.nonNull(ploidyList, "the ploidy model cannot be null");
        Utils.nonNull(likelihoods, "the likelihood object cannot be null");
        this.ploidyList = ploidyList;
        this.likelihoods = likelihoods;
        if (!ploidyList.asListOfSamples().equals(likelihoods.asListOfSamples())) {
            throw new IllegalArgumentException("sample list are different between ploidy-model and read-likelihood collection, perhaps just the order");
        }
    }

    /**
     * Returns the ploidy model that corresponds to the data provided.
     * @return never {@code null}.
     */
    public PloidyList ploidyModel() {
        return ploidyList;
    }

    @Override
    public int numberOfSamples() {
        return ploidyList.numberOfSamples();
    }

    @Override
    public int indexOfSample(final String sample) {
        Utils.nonNull(sample);
        return ploidyList.indexOfSample(sample);
    }

    @Override
    public String getSample(final int sampleIndex) {
        Utils.validateArg(sampleIndex >= 0, "sampleIndex");
        return ploidyList.getSample(sampleIndex);
    }

    /**
     * Returns read-likelihoods to use for genotyping.
     * @return never {@code null}.
     */
    public ReadLikelihoods<A> readLikelihoods() {
        return likelihoods;
    }

    @Override
    public int numberOfAlleles() {
        return likelihoods.numberOfAlleles();
    }

    @Override
    public int indexOfAllele(final A allele) {
        Utils.nonNull(allele);
        return likelihoods.indexOfAllele(allele);
    }

    @Override
    public A getAllele(final int index) {
        Utils.validateArg(index >= 0, "index");
        return likelihoods.getAllele(index);
    }
}

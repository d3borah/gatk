package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.SampleList;

import java.util.Arrays;
import java.util.Collections;

public final class PloidyList implements SampleList{

    private final SampleList sampleList;

    private final int[] ploidies;

    private final int ploidySum;

    private final boolean isHomogeneous;

    public PloidyList(final SampleList sampleList, final int[] ploidies) {
        this.sampleList = Utils.nonNull(sampleList, "the sample list cannot be null");
        this.ploidies = Utils.nonNull(ploidies, "the ploidies cannot be null").clone();
        if (sampleList.numberOfSamples() != ploidies.length) {
            throw new IllegalArgumentException("sample-list and ploidy array length must match");
        }
        Arrays.stream(ploidies).forEach(p -> {if (p < 0) throw new IllegalArgumentException("no ploidy can be less than 0");});
        ploidySum = (int) MathUtils.sum(ploidies);
        isHomogeneous = ploidies.length == 0 || Arrays.stream(ploidies).allMatch(p -> p == ploidies[0]);
    }

    public PloidyList(final SampleList sampleList, final int ploidy) {
        this(sampleList, Collections.nCopies(sampleList.numberOfSamples(), ploidy).stream().mapToInt(n->n).toArray());
    }

    public int samplePloidy(final int sampleIndex) {
        Utils.validIndex(sampleIndex, ploidies.length);
        return ploidies[sampleIndex];
    }

    public boolean isHomogeneous() {
        return isHomogeneous;
    }

    public int totalPloidy() {
        return ploidySum;
    }

    public int numberOfSamples() {
        return ploidies.length;
    }

    public int indexOfSample(final String sample) {
        Utils.nonNull(sample);
        return sampleList.indexOfSample(sample);
    }

    public String getSample(final int sampleIndex) {
        Utils.validIndex(sampleIndex, numberOfSamples());
        return sampleList.getSample(sampleIndex);
    }
}
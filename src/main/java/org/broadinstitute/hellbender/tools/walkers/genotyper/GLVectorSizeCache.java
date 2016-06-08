package org.broadinstitute.hellbender.tools.walkers.genotyper;

import java.util.Arrays;

/**
 * a cache of the PL vector sizes (i.e. the number of possible unphased genotypes) as a function of # of alleles and ploidy
 *
 * Created by davidben on 6/9/16.
 */
public class GLVectorSizeCache {
    //TODO: THIS MAY BE REDUNDANT with a similar cache in GenotypeLikelihoods
    private GLVectorSizeCache() { } //prevent instantiation
    private static final int MAX_NUM_ALLELES_TO_CACHE = 20;
    private static final int MAX_PLOIDY = 2000; //TODO: this was originally for pooled sequencing -- is this still relevant?

    public static int getNumLikelihoodElements(final int numAlleles, final int ploidy) {
        return GL_VECTOR_SIZE_CACHE[numAlleles][ploidy];
    }

    private static final int[][] GL_VECTOR_SIZE_CACHE = fillGLVectorSizeCache(MAX_NUM_ALLELES_TO_CACHE, MAX_PLOIDY);

    private static int[][] fillGLVectorSizeCache(final int maxAlleles, final int maxPloidy) {
        final int[][] cache = new int[maxAlleles][maxPloidy];
        for (int numAlleles=1; numAlleles < maxAlleles; numAlleles++) {
            for (int ploidy=0; ploidy < maxPloidy; ploidy++) {
                cache[numAlleles][ploidy] = numAlleles == 1 ? 1 : Arrays.stream(cache[numAlleles - 1], 0, ploidy + 1).sum();
            }
        }
        return cache;
    }
}

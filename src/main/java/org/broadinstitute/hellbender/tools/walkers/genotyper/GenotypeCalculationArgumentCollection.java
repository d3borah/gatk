package org.broadinstitute.hellbender.tools.walkers.genotyper;

import org.broadinstitute.hellbender.cmdline.Advanced;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.HomoSapiensConstants;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public final class GenotypeCalculationArgumentCollection {

    /**
     * Creates a GenotypeCalculationArgumentCollection with default values.
     */
    public GenotypeCalculationArgumentCollection() {}

    /**
     * Creates a GenotypeCalculationArgumentCollection with the values from other
     *
     * @param other GenotypeCalculationArgumentCollection from which to copy values
     */
    public GenotypeCalculationArgumentCollection( final GenotypeCalculationArgumentCollection other ) {
        Utils.nonNull(other);
        this.ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED = other.ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED;
        this.STANDARD_CONFIDENCE_FOR_CALLING = other.STANDARD_CONFIDENCE_FOR_CALLING;
        this.STANDARD_CONFIDENCE_FOR_EMITTING = other.STANDARD_CONFIDENCE_FOR_EMITTING;
        this.MAX_ALTERNATE_ALLELES = other.MAX_ALTERNATE_ALLELES;
        this.inputPrior = new ArrayList<>(other.inputPrior);
        this.samplePloidy = other.samplePloidy;
    }

    /**
     * Depending on the value of the --max_alternate_alleles argument, we may genotype only a fraction of the alleles being sent on for genotyping.
     * Using this argument instructs the genotyper to annotate (in the INFO field) the number of alternate alleles that were originally discovered at the site.
     */
    @Argument(fullName = "annotateNDA", shortName = "nda", doc = "If provided, we will annotate records with the number of alternate alleles that were discovered (but not necessarily genotyped) at a given site", optional = true)
    public boolean ANNOTATE_NUMBER_OF_ALLELES_DISCOVERED = false;

    /**
     * The minimum phred-scaled Qscore threshold to separate high confidence from low confidence calls. Only genotypes with
     * confidence >= this threshold are emitted as called sites. A reasonable threshold is 30 for high-pass calling (this
     * is the default).
     */
    @Argument(fullName = "standard_min_confidence_threshold_for_calling", shortName = "stand_call_conf", doc = "The minimum phred-scaled confidence threshold at which variants should be called", optional = true)
    public double STANDARD_CONFIDENCE_FOR_CALLING = 30.0;

    /**
     * This argument allows you to emit low quality calls as filtered records.
     */
    @Argument(fullName = "standard_min_confidence_threshold_for_emitting", shortName = "stand_emit_conf", doc = "The minimum phred-scaled confidence threshold at which variants should be emitted (and filtered with LowQual if less than the calling threshold)", optional = true)
    public double STANDARD_CONFIDENCE_FOR_EMITTING = 30.0;

    /**
     * If there are more than this number of alternate alleles presented to the genotyper (either through discovery or GENOTYPE_GIVEN ALLELES),
     * then only this many alleles will be used.  Note that genotyping sites with many alternate alleles is both CPU and memory intensive and it
     * scales exponentially based on the number of alternate alleles.  Unless there is a good reason to change the default value, we highly recommend
     * that you not play around with this parameter.
     *
     * As of GATK 2.2 the genotyper can handle a very large number of events, so the default maximum has been increased to 6.
     */
    @Advanced
    @Argument(fullName = "max_alternate_alleles", shortName = "maxAltAlleles", doc = "Maximum number of alternate alleles to genotype", optional = true)
    public int MAX_ALTERNATE_ALLELES = 6;

    /**
     * By default, the prior specified with the argument --heterozygosity/-hets is used for variant discovery at a particular locus, using an infinite sites model,
     * see e.g. Waterson (1975) or Tajima (1996).
     * This model asserts that the probability of having a population of k variant sites in N chromosomes is proportional to theta/k, for 1=1:N
     *
     * There are instances where using this prior might not be desireable, e.g. for population studies where prior might not be appropriate,
     * as for example when the ancestral status of the reference allele is not known.
     * By using this argument, user can manually specify priors to be used for calling as a vector for doubles, with the following restriciotns:
     * a) User must specify 2N values, where N is the number of samples.
     * b) Only diploid calls supported.
     * c) Probability values are specified in double format, in linear space.
     * d) No negative values allowed.
     * e) Values will be added and Pr(AC=0) will be 1-sum, so that they sum up to one.
     * f) If user-defined values add to more than one, an error will be produced.
     *
     * If user wants completely flat priors, then user should specify the same value (=1/(2*N+1)) 2*N times,e.g.
     *   -inputPrior 0.33 -inputPrior 0.33
     * for the single-sample diploid case.
     */
    @Advanced
    @Argument(fullName = "input_prior", shortName = "inputPrior", doc = "Input prior for calls", optional = true)
    public List<Double> inputPrior = Collections.emptyList();

    /**
     *   Sample ploidy - equivalent to number of chromosomes per pool. In pooled experiments this should be = # of samples in pool * individual sample ploidy
     */
    @Argument(shortName="ploidy", fullName="sample_ploidy", doc="Ploidy (number of chromosomes) per sample. For pooled data, set to (Number of samples in each pool * Sample Ploidy).", optional=true)
    public int samplePloidy = HomoSapiensConstants.DEFAULT_PLOIDY;
}

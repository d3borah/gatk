package org.broadinstitute.hellbender.metrics;

import htsjdk.samtools.SamPairUtil;

import htsjdk.samtools.util.Histogram;
import htsjdk.samtools.metrics.MetricsFile;

import java.io.Serializable;
import java.util.*;

/** A Collector for individual InsertSizeMetrics for a given SAMPLE or SAMPLE/LIBRARY or SAMPLE/LIBRARY/READ_GROUP (depending on aggregation levels) */
public final class PerUnitInsertSizeMetricsCollector
        implements PerUnitMetricCollector<InsertSizeMetrics, Integer, InsertSizeMetricsCollectorArgs>,
                    Serializable
{
    private static final long serialVersionUID = 1L;

    // TODO: work around kryo serialization bug
    //final EnumMap<SamPairUtil.PairOrientation, Histogram<Integer>> Histograms = new EnumMap<>(SamPairUtil.PairOrientation.class);
    final Map<SamPairUtil.PairOrientation, Histogram<Integer>> Histograms = new HashMap<>();

    final String sample;
    final String library;
    final String readGroup;

    // When generating the Histogram, discard any data categories (out of FR, TANDEM, RF) that have fewer than this
    // percentage of overall reads. (Range: 0 to 1)
    private final double minimumPct;

    // Generate mean, sd and plots by trimming the data down to MEDIAN + DEVIATIONS*MEDIAN_ABSOLUTE_DEVIATION.
    // This is done because insert size data typically includes enough anomolous values from chimeras and other
    // artifacts to make the mean and sd grossly misleading regarding the real distribution.
    private final double deviations;

    //Explicitly sets the Histogram width, overriding automatic truncation of Histogram tail.
    //Also, when calculating mean and stdev, only bins <= HISTOGRAM_WIDTH will be included.
    private Integer HistogramWidth;


    private double totalInserts = 0;

    public PerUnitInsertSizeMetricsCollector(
            final String sample,
            final String library,
            final String readGroup,
            final double minimumPct,
            final double deviations,
            final Integer HistogramWidth) {

        super();

        this.sample = sample;
        this.library = library;
        this.readGroup = readGroup;
        this.minimumPct = minimumPct;
        this.HistogramWidth = HistogramWidth;
        this.deviations = deviations;

        String prefix = null;
        if (this.readGroup != null) {
            prefix = this.readGroup + ".";
        }
        else if (this.library != null) {
            prefix = this.library + ".";
        }
        else if (this.sample != null) {
            prefix = this.sample + ".";
        }
        else {
            prefix = "All_Reads.";
        }

        Histograms.put(SamPairUtil.PairOrientation.FR,     new Histogram<Integer>("insert_size", prefix + "fr_count"));
        Histograms.put(SamPairUtil.PairOrientation.TANDEM, new Histogram<Integer>("insert_size", prefix + "tandem_count"));
        Histograms.put(SamPairUtil.PairOrientation.RF,     new Histogram<Integer>("insert_size", prefix + "rf_count"));
    }

    public void acceptRecord(final InsertSizeMetricsCollectorArgs args) {
        Histograms.get(args.getPairOrientation()).increment(args.getInsertSize());
    }

    public void finish() { }

    public double getTotalInserts() {
        return totalInserts;
    }

    @SuppressWarnings("unchecked")
    public void addMetricsToFile(final MetricsFile<InsertSizeMetrics, Integer> file) {

        for(final Map.Entry<SamPairUtil.PairOrientation, Histogram<Integer>> entry : Histograms.entrySet()) {
            totalInserts += entry.getValue().getCount();
        }

        for(final Map.Entry<SamPairUtil.PairOrientation, Histogram<Integer>> entry : Histograms.entrySet()) {

            final SamPairUtil.PairOrientation pairOrientation = entry.getKey();
            final Histogram<Integer> Histogram = entry.getValue();

            final double total = Histogram.getCount();

            // Only include a category if it has a sufficient percentage of the data in it
            if( total > totalInserts * minimumPct ) {
                final InsertSizeMetrics metrics = new InsertSizeMetrics();
                metrics.SAMPLE             = this.sample;
                metrics.LIBRARY            = this.library;
                metrics.READ_GROUP         = this.readGroup;
                metrics.PAIR_ORIENTATION   = pairOrientation;
                metrics.READ_PAIRS         = (long) total;
                metrics.MAX_INSERT_SIZE    = (int) Histogram.getMax();
                metrics.MIN_INSERT_SIZE    = (int) Histogram.getMin();
                metrics.MEDIAN_INSERT_SIZE = Histogram.getMedian();
                metrics.MEDIAN_ABSOLUTE_DEVIATION = Histogram.getMedianAbsoluteDeviation();

                final double median  = Histogram.getMedian();
                double covered = 0;
                double low  = median;
                double high = median;

                while (low >= Histogram.getMin() || high <= Histogram.getMax()) {
                    final Histogram.Bin<Integer> lowBin = Histogram.get((int) low);
                    if (lowBin != null) covered += lowBin.getValue();

                    if (low != high) {
                        final Histogram.Bin<Integer> highBin = Histogram.get((int) high);
                        if (highBin != null) covered += highBin.getValue();
                    }

                    final double percentCovered = covered / total;
                    final int distance = (int) (high - low) + 1;
                    if (percentCovered >= 0.1  && metrics.WIDTH_OF_10_PERCENT == 0) metrics.WIDTH_OF_10_PERCENT = distance;
                    if (percentCovered >= 0.2  && metrics.WIDTH_OF_20_PERCENT == 0) metrics.WIDTH_OF_20_PERCENT = distance;
                    if (percentCovered >= 0.3  && metrics.WIDTH_OF_30_PERCENT == 0) metrics.WIDTH_OF_30_PERCENT = distance;
                    if (percentCovered >= 0.4  && metrics.WIDTH_OF_40_PERCENT == 0) metrics.WIDTH_OF_40_PERCENT = distance;
                    if (percentCovered >= 0.5  && metrics.WIDTH_OF_50_PERCENT == 0) metrics.WIDTH_OF_50_PERCENT = distance;
                    if (percentCovered >= 0.6  && metrics.WIDTH_OF_60_PERCENT == 0) metrics.WIDTH_OF_60_PERCENT = distance;
                    if (percentCovered >= 0.7  && metrics.WIDTH_OF_70_PERCENT == 0) metrics.WIDTH_OF_70_PERCENT = distance;
                    if (percentCovered >= 0.8  && metrics.WIDTH_OF_80_PERCENT == 0) metrics.WIDTH_OF_80_PERCENT = distance;
                    if (percentCovered >= 0.9  && metrics.WIDTH_OF_90_PERCENT == 0) metrics.WIDTH_OF_90_PERCENT = distance;
                    if (percentCovered >= 0.99 && metrics.WIDTH_OF_99_PERCENT == 0) metrics.WIDTH_OF_99_PERCENT = distance;

                    --low;
                    ++high;
                }

                // Trim the Histogram down to get rid of outliers that would make the chart useless.
                final Histogram<Integer> trimmedHisto = Histogram; //alias it
                if (HistogramWidth == null) {
                    HistogramWidth = (int) (metrics.MEDIAN_INSERT_SIZE + (deviations * metrics.MEDIAN_ABSOLUTE_DEVIATION));
                }

                trimmedHisto.trimByWidth(HistogramWidth);

                metrics.MEAN_INSERT_SIZE = trimmedHisto.getMean();
                metrics.STANDARD_DEVIATION = trimmedHisto.getStandardDeviation();

                file.addHistogram(trimmedHisto);
                file.addMetric(metrics);
            }
        }
    }

    /**
     * Combine this InsertSizeMetricsCollector with another InsertSizeMetricsCollector
     * @param sourceCollector InsertSizeMetricsCollector to combine
     * @return InsertSizeMetricsCollector representing the combination of the source collector with this collector
     */
    public PerUnitInsertSizeMetricsCollector combine(PerUnitInsertSizeMetricsCollector sourceCollector) {

        String validationMessage = "Internal error combining collectors";
        validateEquals(this.sample, sourceCollector.sample, validationMessage);
        validateEquals(this.library, sourceCollector.library, validationMessage);
        validateEquals(this.readGroup, sourceCollector.readGroup, validationMessage);

        PerUnitInsertSizeMetricsCollector combinedCollector = new PerUnitInsertSizeMetricsCollector(
                this.sample,
                this.library,
                this.readGroup,
                this.minimumPct,
                this.deviations,
                this.HistogramWidth);
        combinedCollector.totalInserts = this.totalInserts + sourceCollector.totalInserts;

        //TODO: the initialization code above adds an entry for each Pair orientation
        this.Histograms.forEach( //TODO: do we need to clone these histograms ?
                (po, h) -> {
                    // TODO: assert bin and values are the same ??
                    Histogram<Integer> targetHist = new Histogram<>(h.getBinLabel(), h.getValueLabel());
                    targetHist.addHistogram(sourceCollector.Histograms.get(po));
                    targetHist.addHistogram(this.Histograms.get(po));
                    combinedCollector.Histograms.put(po, targetHist);
                }
        );
        return combinedCollector;
    }

    //TODO Move these to a utility class
    private void validateEquals(final String source, final String target, final String message) {
        if (!safeEqualStrings(source, target)) {
            throw new IllegalStateException(
                    String.format("%s (%s : %s)",
                            message,
                            source == null ? "null" : source,
                            target == null ? "null" : target)
            );
        }
    }

    private boolean safeEqualStrings(final String source, final String target) {
        return source == null ?
                target == null ? true : false :
                source.equals(target);
    }

}

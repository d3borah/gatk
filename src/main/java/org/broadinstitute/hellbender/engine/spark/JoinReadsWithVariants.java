package org.broadinstitute.hellbender.engine.spark;

import com.beust.jcommander.internal.Lists;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.broadinstitute.hellbender.engine.dataflow.datasources.VariantShard;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.variant.Variant;
import scala.Tuple2;

import java.util.List;

/**
 * PairReadsAndVariants takes two RDDs (GATKRead and Variant) and returns an RDD with a
 * (GATKRead,Iterable<Variant>) for every read and all the variants that overlap. We do this by first
 * Making KVs of GATKRead,Variant, which means that a read or variant may be present
 * multiple times in the output. Also, there may be duplicate (GATKRead,Variant) pairs in the output.
 * Currently, all reads must be mapped. We then dedup (distinct) then group by key.
 *
 * The function works by creating a single RDD of KV where GATKRead is the key and Variant is the value.
 * We do this join by sharding both collections by "variant shard" and checking for overlap on each "shard."
 *
 * step 1: key reads and variants by shard
 * |---- shard 0 -----|---- shard 1 -----|---- shard 2 -----|---- shard 3 -----|---- shard 4 -----|
 *     |---------- read a ---------|               |----- read b ------|
 *   |- variant 1 -|    |- variant 2 -|               |- variant 3 -|
 *
 * step 2: shard read and variant by variant shard
 *                      |---- shard 0 -----|
 *                          |---------- read a ---------|
 *                        |- variant 1 -|
 *
 *
 *                      |---- shard 1 -----|
 *       |---------- read a ---------|
 *                        |- variant 2 -|
 *
 *
 *                      |---- shard 2 -----|
 *                                |----- read b ------|
 *                                   |- variant 3 -|
 *
 *                      |---- shard 3 -----|
 *             |----- read b ------|
 *                |- variant 3 -|
 *
 * step 3: pair reads and variants
 * KV<read a, variant 1> // from shard 0
 * KV<read a, variant 2> // from shard 1
 * KV<read b, variant 3> // from shard 2
 * KV<read b, variant 3> // from shard 3
 *
 * step 4: distinct and group by key
 * KV<read a, <variant 1, variant2>>
 * KV<read b, <variant 3>>
 */
public class JoinReadsWithVariants {
    public static JavaPairRDD<GATKRead, Iterable<Variant>> join(
            final JavaRDD<GATKRead> reads, final JavaRDD<Variant> variants) {

        JavaPairRDD<VariantShard, GATKRead> readsWShards = reads.flatMapToPair(gatkRead -> {
            List<VariantShard> shards = VariantShard.getVariantShardsFromInterval(gatkRead);
            List<Tuple2<VariantShard, GATKRead>> out = Lists.newArrayList();
            for (VariantShard shard : shards) {
                out.add(new Tuple2<>(shard, gatkRead));
            }
            return out;
        });

        JavaPairRDD<VariantShard, Variant> variantsWShards = variants.flatMapToPair(variant -> {
            List<VariantShard> shards = VariantShard.getVariantShardsFromInterval(variant);
            List<Tuple2<VariantShard, Variant>> out = Lists.newArrayList();
            for (VariantShard shard : shards) {
                out.add(new Tuple2<>(shard, variant));
            }
            return out;
        });

        JavaPairRDD<VariantShard, Tuple2<Iterable<GATKRead>, Iterable<Variant>>> cogroup = readsWShards.cogroup(variantsWShards);

        JavaPairRDD<GATKRead, Variant> allPairs = cogroup.flatMapToPair(cogroupValue -> {
            Iterable<GATKRead> iReads = cogroupValue._2()._1();
            Iterable<Variant> iVariants = cogroupValue._2()._2();

            List<Tuple2<GATKRead, Variant>> out = Lists.newArrayList();
            // For every read, find every overlapping variant.
            for (GATKRead r : iReads) {
                boolean foundVariants = false;
                SimpleInterval interval = new SimpleInterval(r);
                for (Variant v : iVariants) {
                    if (interval.overlaps(v)) {
                        foundVariants = true;
                        out.add(new Tuple2<>(r, v));
                    }
                }
                // If no variants are found, we still want to output the read.
                if (!foundVariants) {
                    out.add(new Tuple2<>(r, null));
                }
            }
            return out;
        });
        return allPairs.distinct().groupByKey();
    }
}
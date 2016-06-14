package org.broadinstitute.hellbender.tools.spark.utils;

import java.util.Collection;
import java.util.Map;
import java.util.Objects;
import java.util.function.BiPredicate;
import java.util.function.Function;

/**
 * A uniquely keyed map with O(1) operations.
 * Sadly, it's not a java.util.Map, but it does behave like a java.util.Map's entrySet.
 */
public final class HopscotchMap<K, V, T extends Map.Entry<K, V>>  extends HopscotchSet<T> {
    public HopscotchMap() {}
    public HopscotchMap( final int capacity ) { super(capacity); }
    public HopscotchMap( final Collection<? extends T> collection ) { super(collection); }

    @Override
    protected BiPredicate<T, T> entryCollides() { return (t1, t2) -> Objects.equals(t1.getKey(), t2.getKey()); }

    @Override
    protected Function<T, Object> toKey() { return Map.Entry::getKey; }
}

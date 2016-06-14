package org.broadinstitute.hellbender.tools.spark.utils;

import java.util.Collection;
import java.util.Map;
import java.util.function.Function;

/**
 * A map that can contain multiple values for a given key.
 */
public final class HopscotchMultiMap<K, V, T extends Map.Entry<K, V>>  extends HopscotchCollection<T> {
    public HopscotchMultiMap() {}
    public HopscotchMultiMap( final int capacity ) { super(capacity); }
    public HopscotchMultiMap( final Collection<? extends T> collection ) { super(collection); }

    protected Function<T, Object> toKey() { return Map.Entry::getKey; }
}

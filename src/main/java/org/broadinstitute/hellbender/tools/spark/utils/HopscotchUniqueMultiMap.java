package org.broadinstitute.hellbender.tools.spark.utils;

import java.util.Collection;
import java.util.Map;
import java.util.function.BiPredicate;
import java.util.function.Function;

/**
 * A map that can contain multiple values for a given key, but distinct entries.
 */
public final class HopscotchUniqueMultiMap<K, V, T extends Map.Entry<K, V>>  extends HopscotchMultiMap<K, V, T> {
    public HopscotchUniqueMultiMap() {}
    public HopscotchUniqueMultiMap( final int capacity ) { super(capacity); }
    public HopscotchUniqueMultiMap( final Collection<? extends T> collection ) { super(collection); }

    protected BiPredicate<T, T> entryCollides() { return Object::equals; }
}

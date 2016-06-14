package org.broadinstitute.hellbender.tools.spark.utils;

import java.util.Collection;
import java.util.Objects;
import java.util.Set;
import java.util.function.BiPredicate;

/**
 * Implements Set by imposing a unique-element constraint on HopscotchCollection.
 * Also implements equals and hashCode to be consistent with the documented requirements of the java Set interface.
 */
public class HopscotchSet<T> extends HopscotchCollection<T> implements Set<T> {
    public HopscotchSet() {}
    public HopscotchSet( final int capacity ) { super(capacity); }
    public HopscotchSet( final Collection<? extends T> collection ) { super(collection); }

    @Override
    public boolean equals( final Object obj ) {
        if ( this == obj ) return true;
        if ( !(obj instanceof Set) ) return false;
        final Set that = (Set)obj;
        return this.size() == that.size() && this.containsAll(that);
    }

    @Override
    public int hashCode() { return stream().mapToInt(Objects::hashCode).sum(); }

    @Override
    protected BiPredicate<T, T> entryCollides() { return Objects::equals; }
}

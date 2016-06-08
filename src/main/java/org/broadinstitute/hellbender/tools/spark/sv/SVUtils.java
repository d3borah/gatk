package org.broadinstitute.hellbender.tools.spark.sv;

import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * Useful scraps of this and that.
 */
public final class SVUtils {

    /** return a good initialCapacity for a HashMap that will hold a given number of elements */
    public static int hashMapCapacity( final int nElements )
    {
        return (int)((nElements*4L)/3) + 1;
    }

    /** count the number of items available from an iterator */
    public static <T> int iteratorSize( final Iterator<T> itr ) {
        int result = 0;
        while ( itr.hasNext() ) { result += 1; itr.next(); }
        return result;
    }
}

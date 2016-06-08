package org.broadinstitute.hellbender.tools.walkers.genotyper.afcalc;

import org.broadinstitute.hellbender.tools.walkers.genotyper.GenotpyeIterator;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.util.Arrays;
import java.util.List;

/**
 * Created by davidben on 6/9/16.
 */
public class GenotpyeIteratorUnitTest extends BaseTest {

    @Test(dataProvider = "expectedIterationSequences")
    public void testIteration(final int numAlleles, final int ploidy, final List<int[]> expectedIteration) {
        final GenotpyeIterator iterator = new GenotpyeIterator(numAlleles, ploidy);
        int linearIndex = 0;
        for (final int[] expectedAlleleCounts : expectedIteration) {
            int j = 30;
            Assert.assertTrue(iterator.hasNext());
            Assert.assertTrue(Arrays.equals(iterator.getAlleleCounts(), expectedAlleleCounts));
            Assert.assertEquals(GenotpyeIterator.getLinearIndex(iterator.getAlleleCounts(), numAlleles, ploidy), linearIndex);
            Assert.assertEquals(iterator.getLinearIndex(), linearIndex);
            iterator.next();
            linearIndex++;
        }
        Assert.assertFalse(iterator.hasNext());
    }


    // data is (int (numAlleles), int (ploidy), List<int[]> (sequence of allele count vectors))
    @DataProvider(name="expectedIterationSequences")
    public Object[][] correctRunData() {
        return new Object[][] {
                {1, 0, Arrays.asList(new int[] {0})},
                {1, 1, Arrays.asList(new int[] {1})},
                {1, 2, Arrays.asList(new int[] {2})},
                {2, 1, Arrays.asList(new int[] {1,0}, new int[] {0,1})},
                {2, 2, Arrays.asList(new int[] {2,0}, new int[] {1,1}, new int[] {0,2})},
                {3, 3, Arrays.asList(new int[] {3,0,0}, new int[] {2,1,0}, new int[] {1,2,0}, new int[] {0,3,0},
                        new int[] {2,0,1}, new int[] {1,1,1}, new int[] {0,2,1}, new int[] {1,0,2}, new int[] {0,1,2},
                        new int[] {0,0,3})}
        };
    }

}
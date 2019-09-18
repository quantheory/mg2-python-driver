#!/usr/bin/env python3

from functools import reduce
import unittest

import numpy as np

class Reduction:

    """Operation that reduces an array to a scalar.

    Methods:
    __init__
    __call__
    has_ltr
    """

    def __init__(self, array_op=None, ltr_op=None, ltr_start=None,
                 ltr_finalize=None):
        """Initialize a reduction from input functions.

        Arguments:
        array_op - A function that does a reduction given an array. If this
                   function is specified, it defines the reduction completely.
                   Otherwise, ltr_op must be specified.
        ltr_op - A function that takes two arguments, and produces an output
                 representing their combination. This function is used to
                 reduce the array if array_op is not present.
        ltr_start - If present, ltr_start is the first argument of ltr_op the
                    first time it is called, while the first element of the
                    array is the second argument. Otherwise, the first two
                    elements of the array are used as the arguments.
        ltr_finalize - If present, this function is called on the result of the
                       final ltr_op call. Otherwise, the result is returned
                       unchanged.
        """
        self.ltr_op = ltr_op
        self.ltr_start = ltr_start
        self.ltr_finalize = ltr_finalize
        if array_op is None:
            assert ltr_op is not None, \
                "Cannot create a Reduction with no input function."
            if ltr_finalize is None:
                ltr_finalize = lambda x: x
            if ltr_start is None:
                self.reduce = lambda arr: ltr_finalize(reduce(ltr_op, arr))
            else:
                self.reduce = lambda arr: ltr_finalize(reduce(ltr_op, arr, ltr_start))
        else:
            self.reduce = array_op

    def __call__(self, array):
        """Reduce an array."""
        return self.reduce(array)

    def has_ltr(self):
        """Whether or not this reduction uses an ltr_op."""
        return self.ltr_op is not None

# Save builtin sum for testing purposes.
_sum = sum

sum = Reduction(ltr_op=lambda x,y: x+y,
                ltr_start=0.)

mean = Reduction(ltr_op=lambda x,y:(x[0]+y,x[1]+1),
                 ltr_start=(0.,0),
                 ltr_finalize=lambda x: x[0]/x[1])
"""Reduction that takes the mean of an array."""

def _median_func(array):
    length = len(array)
    sorted_array = sorted(array)
    if length % 2 == 0:
        return 0.5 * (sorted_array[length//2 - 1] + sorted_array[length//2])
    else:
        return sorted_array[length//2]

median = Reduction(array_op=_median_func)
"""Reduction that takes the median of an array."""

max = Reduction(ltr_op=max)
"""Reduction that takes the max of an array."""

min = Reduction(ltr_op=min)
"""Reduction that takes the min of an array."""

def percentile(percent):
    """Function that produces a Reduction finding a given percentile.

    Note that percentile(0.) is equivalent to min, percentile(100.) is
    equivalent to max, and percentile(50.) is equivalent to median, up to
    rounding error.
    """
    if percent == 100.:
        return max
    frac = percent*0.01
    def array_op(array):
        ngaps = len(array)-1
        weight, index = np.modf(frac*ngaps)
        index = int(index)
        sort_array = sorted(array)
        return sort_array[index] + weight*(sort_array[index+1]-sort_array[index])
    return Reduction(array_op=array_op)


class Accumulator:

    """Object that accumulates a set of values being reduced.

    Methods:
    __init__
    push
    output
    merge
    merge_series
    """

    def __init__(self, reductions):
        """Create an Accumulator from a set of reductions."""
        self._reductions = reductions
        self._series = []

    def push(self, value):
        """Add a new value to the series being accumulated."""
        self._series.append(value)

    def output(self):
        """Output the accumulated value."""
        return [red(self._series) for red in self._reductions]

    def merge(self, other):
        """Merge another accumulator into this one."""
        self._series += other._series

    def merge_series(self, series):
        """Add a whole to this accumulator."""
        self._series += series


class TestReductions(unittest.TestCase):

    @staticmethod
    def closest_to_mean(array):
        mean = _sum(array) / len(array)
        closest = array[0]
        for x in array[1:]:
            if abs(x - mean) < abs(closest - mean):
                closest = x
        return closest

    def test_array_op(self):
        close_reduce = Reduction(array_op=self.closest_to_mean)
        self.assertEqual(close_reduce.reduce([0,1,2.1,3,4]), 2.1)

    def test_call_interface(self):
        close_reduce = Reduction(array_op=self.closest_to_mean)
        self.assertEqual(close_reduce([0,4,2.1,3,1]), 2.1)

    def test_ltr(self):
        sum_reduce = Reduction(ltr_op=lambda x,y: x+y)
        self.assertEqual(sum_reduce.ltr_op(1.1, 2.5), 3.6)

    def test_ltr_reduce(self):
        sum_reduce = Reduction(ltr_op=lambda x,y: x+y)
        self.assertEqual(sum_reduce([1.1, 2.5, 6.2]), 9.8)

    def test_no_argument_error(self):
        with self.assertRaises(AssertionError):
            Reduction()

    def test_ltr_start(self):
        sum_len_reduce = Reduction(ltr_op=lambda x,y:(x[0]+y,x[1]+1),
                                   ltr_start=(0,0))
        total, length = sum_len_reduce.ltr_op(sum_len_reduce.ltr_start, 4.)
        self.assertEqual(total, 4.)
        self.assertEqual(length, 1)

    def test_ltr_start_reduce(self):
        sum_len_reduce = Reduction(ltr_op=lambda x,y:(x[0]+y,x[1]+1),
                                   ltr_start=(0.,0))
        self.assertEqual(sum_len_reduce([1.1, 2.5, 6.2]), (9.8, 3))

    def test_ltr_finalize(self):
        sum_squared_reduce = Reduction(ltr_op=lambda x,y: x+y,
                                       ltr_finalize=lambda x: x*x)
        x = sum_squared_reduce.ltr_finalize(sum_squared_reduce.ltr_op(2., 3.))
        self.assertEqual(x, 25.)

    def test_ltr_finalize_reduce(self):
        sum_squared_reduce = Reduction(ltr_op=lambda x,y: x+y,
                                       ltr_finalize=lambda x: x*x)
        x = sum_squared_reduce([2., 3., 4.])
        self.assertEqual(x, 81.)

    def test_ltr_start_finalize_reduce(self):
        mean_reduce = Reduction(ltr_op=lambda x,y:(x[0]+y,x[1]+1),
                                ltr_start=(0.,0),
                                ltr_finalize=lambda x: x[0]/x[1])
        x = mean_reduce([2., 3., 4., 6., 12., 0., 1.])
        self.assertEqual(x, 4.)

    def test_has_ltr(self):
        close_reduce = Reduction(array_op=self.closest_to_mean)
        sum_reduce = Reduction(ltr_op=lambda x,y: x+y)
        self.assertFalse(close_reduce.has_ltr())
        self.assertTrue(sum_reduce.has_ltr())

    def test_mean(self):
        x = mean([2., 3., 4., 6., 12., 0., 1.])
        self.assertTrue(mean.has_ltr())
        self.assertEqual(x, 4.)

    def test_median(self):
        array = [2., 3., 4., 6., 12., 0., 1.]
        x = median(array)
        x2 = median(array[1:])
        self.assertEqual(x, 3.)
        self.assertEqual(x2, 3.5)

    def test_max(self):
        array = [2., 3., 4., 6., 12., 0., 1.]
        self.assertTrue(max.has_ltr())
        x = max(array)
        self.assertEqual(x, 12.)

    def test_min(self):
        array = [2., 3., 4., 6., 12., 0., 1.]
        self.assertTrue(min.has_ltr())
        x = min(array)
        self.assertEqual(x, 0.)

    def test_sum(self):
        array = [2., 3., 4., 6., 12., 0., 1.]
        self.assertTrue(min.has_ltr())
        x = sum(array)
        self.assertEqual(x, 28.)

    def test_sum_empty(self):
        x = sum([])
        self.assertEqual(x, 0.)

    def test_percentile(self):
        array = [20., 3., 4., 6., 12., 40., 25.]
        test_percent = 80.
        self.assertAlmostEqual(percentile(test_percent)(array), 24.)

    def test_percentile_0(self):
        array = [20., 3., 4., 6., 12., 40., 25.]
        test_percent = 0.
        self.assertEqual(percentile(test_percent)(array), 3.)

    def test_percentile_100(self):
        array = [20., 3., 4., 6., 12., 40., 25.]
        test_percent = 100.
        self.assertEqual(percentile(test_percent)(array), 40.)


class TestAccumulator(unittest.TestCase):

    def test_accumulator_array(self):
        median_acc = Accumulator([median])
        array = [2., 3., 4., 6., 12., 0., 1.]
        for num in array:
            median_acc.push(num)
        self.assertEqual(median_acc.output()[0], 3.)

    def test_accumulator_multi_array(self):
        median_acc = Accumulator([median, percentile(100./3.)])
        array = [2., 3., 4., 6., 12., 0., 1.]
        for num in array:
            median_acc.push(num)
        self.assertEqual(median_acc.output()[0], 3.)
        self.assertAlmostEqual(median_acc.output()[1], 2.)

    def test_accumulator_merge(self):
        median_acc = Accumulator([median, percentile(100./3.)])
        median2_acc = Accumulator([median, percentile(100./3.)])
        array = [2., 3., 4., 6., 12., 0., 1.]
        for i in range(len(array)):
            if i < 3:
                median_acc.push(array[i])
            else:
                median2_acc.push(array[i])
        median_acc.merge(median2_acc)
        self.assertEqual(median_acc.output()[0], 3.)
        self.assertAlmostEqual(median_acc.output()[1], 2.)

    def test_accumulator_merge_series(self):
        median_acc = Accumulator([median, percentile(100./3.)])
        series = []
        array = [2., 3., 4., 6., 12., 0., 1.]
        for i in range(len(array)):
            if i < 3:
                median_acc.push(array[i])
            else:
                series.append(array[i])
        median_acc.merge_series(series)
        self.assertEqual(median_acc.output()[0], 3.)
        self.assertAlmostEqual(median_acc.output()[1], 2.)

from __future__ import print_function
import numpy
from frigus.utils import linear_2d_index


def test_that_linear_2d_index_has_no_clashes():
    n = 1000000
    i, j = numpy.unique(
        numpy.vstack(
            (numpy.random.randint(0, 100, n),
             numpy.random.randint(0, 100, n))
        ), axis=1
    )

    linear_inds = linear_2d_index(i, j)
    assert linear_inds.size == numpy.unique(linear_inds).size

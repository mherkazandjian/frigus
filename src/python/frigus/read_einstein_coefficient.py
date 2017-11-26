# -*- coding: utf-8 -*-
"""
read the Einstein coefficients
"""

import numpy
from numpy import zeros

from astropy import units as u


def read_einstein_coppola():
    """read the data provided by Coppola for transitions (v', j') -> (v'', j'')
     with the limitation that delta j i.e |j'' - j'| = 1

    :return: A 4D matrix holding the A coefficients. A[v', j', v'', j'']

    .. code-block:: python

        .. todo:: add example
    """

    # lists that store the read data. These are the nonzero entries of A
    # that is also returned
    vp_nnz, jp_nnz, vpp_nnz, jpp_nnz, A_nnz = [], [], [], [], []

    data = numpy.loadtxt(
        '../../../data/read/lipovka/hd_einstein_coeffs.dat',
        skiprows=2)

    vp_nnz = data[:, 0].astype(numpy.int32)
    jp_nnz = data[:, 1].astype(numpy.int32)
    vpp_nnz = data[:, 2].astype(numpy.int32)
    jpp_nnz = data[:, 3].astype(numpy.int32)
    A_nnz = data[:, 5].astype(numpy.float64)

    jmax = jp_nnz.max()
    vmax = vp_nnz.max()

    # define the A matrix (vmax and jmax assume zero indexing that is how they
    # they are provided in the data files)
    A = zeros((vmax + 1, jmax + 1, vmax + 1, jmax + 1), 'f8')

    A[vp_nnz, jp_nnz, vpp_nnz, jpp_nnz] = A_nnz

    # .. todo:: add these tests
    # testing.assert_approx_equal(A.sum(), 9.3724e-4, significant=4)
    # testing.assert_approx_equal(A[7, 2, 4, 2], 4.133e-7, significant=4)

    # setting units to the Einstein coefficients that will be returned
    A_with_units = A / u.second
    A_nnz_with_units = A_nnz / u.second

    return A_with_units, (
        vp_nnz,
        jp_nnz,
        vpp_nnz,
        jpp_nnz,
        A_nnz_with_units)

from __future__ import print_function
import numpy
from numpy.testing import assert_allclose
from astropy import units as u
import py.test
from frigus.population import fit_lipovka


def test_that_lipovka_cooling_function_fit_is_computed_correctly():

    n_hd = 1e6 * u.m**-3
    T = numpy.linspace(100.0, 2000.0, 100) * u.K

    cooling_function = fit_lipovka(T, n_hd)

    expected_values = [2.218094e-25, 7.523573e-24, 1.984604e-23]
    assert_allclose(
        cooling_function[[0, 50, -1]].cgs.value,
        expected_values,
        rtol=1e-6,
        atol=0.0
    )


from __future__ import print_function
import numpy
from numpy.testing import assert_allclose
from astropy import units as u
import py.test
from population import fit_lipovka, fit_glover, fit_lipovka_low_density


def test_that_lipovka_cooling_function_fit_is_computed_correctly_at_low_densities():

    n_hd = 1e6 * u.m**-3
    T = numpy.linspace(100.0, 2000.0, 100) * u.K

    cooling_function = fit_lipovka(T, n_hd)

    expected_values = [2.218094e-25, 7.523573e-24, 1.984604e-23]
    assert_allclose(
        cooling_function[[0, 50, -1]].cgs.value,
        expected_values,
        rtol=1e-6,
        atol=0.0)


def test_that_lipovka_cooling_function_fit_is_computed_correctly_at_high_denisties():

    py.test.skip('not implemented')

def test_that_glover_cooling_function_fit_is_computed_correctly():
    expected_values = [3.6341638e-29, 1.5221534e-26, 4.8841725e-25]
    assert_allclose(
        u.Quantity([
            fit_glover(100.0),
            fit_glover(500.0),
            fit_glover(1000.0)]).cgs.value,
        expected_values,
        rtol=1e-6,
        atol=0.0)

from __future__ import print_function
import numpy
from numpy.testing import assert_allclose
from astropy import units as u
from frigus.cooling_function.fits import (fit_lipovka,
                                          fit_glover,
                                          fit_lipovka_low_density)


def test_that_glover_cooling_function_fit_is_computed_correctly():

    computed_values = u.Quantity(
        [
            fit_glover(100.0),
            fit_glover(500.0),
            fit_glover(1000.0)
        ]).cgs.value

    expected_values = [3.6341638e-29, 1.5221534e-26, 4.8841725e-25]

    assert_allclose(
        computed_values,
        expected_values,
        rtol=1e-6,
        atol=0.0
    )

def test_that_lipovka_cooling_function_fit_is_computed_correctly():

    t_kin = numpy.linspace(100.0, 2000.0, 3) * u.K
    n_hd = numpy.linspace(1.0e6, 1.0e14, 3) * u.m**-3

    computed_values = fit_lipovka(t_kin, n_hd)

    expected_values = [2.218094e-25, 6.03779267e-19, 2.81066527e-18]

    assert_allclose(
        computed_values.cgs.value,
        expected_values,
        rtol=1e-6,
        atol=0.0
    )


def test_that_lipovka_cooling_low_density_function_fit_is_computed_correctly():

    t_kin = numpy.linspace(100.0, 20000.0, 3) * u.K

    computed_values = fit_lipovka_low_density(t_kin)

    expected_values = [2.38935126e-25, 2.71299166e-22, 7.11558046e-22]

    assert_allclose(
        computed_values.cgs.value,
        expected_values,
        rtol=1e-6,
        atol=0.0
    )

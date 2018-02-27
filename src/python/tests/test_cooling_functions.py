from __future__ import print_function
from numpy.testing import assert_allclose
from astropy import units as u

from frigus.population import cooling_rate_at_steady_state
from frigus.readers.dataset import DataLoader


def test_that_the_lipovka_cooling_function_is_computed_correctly():

    t_rad = 0.0 * u.K
    nc_h = 1e6 * u.m ** -3

    species_data = DataLoader().load('HD_lipovka')

    t_array = u.Quantity([100.0, 500.0, 1000.0, 1500.0, 2000.0]) * u.K
    cooling_rate_expected = u.Quantity(
        [
            2.239044e-25,    # T=100K
            2.861033e-24,    # T=500K
            6.783145e-24,    # T=1000K
            1.108309e-23,    # T=1500K
            1.570218e-23     # T=2000K
        ]
    ) * u.erg / u.s

    cooling_rate_for_t_array = u.Quantity(
        [
            cooling_rate_at_steady_state(
                species_data,
                t,
                t_rad,
                nc_h
            )
            for t in t_array
        ]
    )

    assert_allclose(
        cooling_rate_for_t_array.cgs.value, cooling_rate_expected.cgs.value,
        rtol=1e-6, atol=0.0
    )

def test_that_the_lique_cooling_function_is_computed_correctly():

    t_rad = 0.0 * u.K
    nc_h = 1e6 * u.m ** -3

    species_data = DataLoader().load('H2_lique')

    t_array = u.Quantity([100.0, 500.0, 1000.0, 2000.0, 5000.0]) * u.K
    cooling_rate_expected = u.Quantity(
        [
            1.091851e-28,  # T=100K
            2.674649e-26,  # T=500K
            6.645044e-25,  # T=1000K
            8.351277e-24,  # T=2000K
            1.016337e-22,  # T=5000K
        ]
    ) * u.erg / u.s

    cooling_rate_for_t_array = u.Quantity(
        [
            cooling_rate_at_steady_state(
                species_data,
                t,
                t_rad,
                nc_h
            )
            for t in t_array
        ]
    )

    assert_allclose(
        cooling_rate_for_t_array.cgs.value, cooling_rate_expected.cgs.value,
        rtol=1e-6, atol=0.0
    )

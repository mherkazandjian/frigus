from __future__ import print_function
import pytest
import numpy
from numpy.testing import assert_allclose

from astropy import units as u

from frigus.readers.dataset import DataLoader
from frigus.population import population_density_at_steady_state
from frigus.solvers import analytic



def test_that_two_level_analytic_pop_dens_ratios_are_computed_correctly():

    #
    # load the two level species data
    #
    species_data = DataLoader().load('two_level_1')

    n_c = 1e6 * u.meter ** -3
    t_rad = 1000.0 * u.Kelvin

    #
    # for each value in the temperature range
    #  - compute the population density ration numerically
    #  - compute the population density analytically
    #
    # plot the comaprision
    #
    ratio_numeric_vs_t_kin = []
    ratio_analytic_vs_t_kin = []
    t_kin_range = numpy.logspace(1.0, 10.0, 50) * u.Kelvin
    for t_kin in t_kin_range:

        pop_dens_eq_numeric = population_density_at_steady_state(
            species_data,
            t_kin,
            t_rad,
            n_c)
        pop_dens_eq_ratio_numeric = pop_dens_eq_numeric[1] / pop_dens_eq_numeric[0]
        ratio_numeric_vs_t_kin.append(pop_dens_eq_ratio_numeric)

        pop_dens_eq_ratio_analytic = analytic.population_density_ratio_two_level(
            species_data.energy_levels.data['g'],
            species_data.energy_levels.data['E'],
            species_data.k_dex_matrix_interpolator(t_kin)[1, 0],
            species_data.a_matrix[1, 0],
            n_c,
            t_kin,
            t_rad
        )
        ratio_analytic_vs_t_kin.append(pop_dens_eq_ratio_analytic)
        print(t_kin)

    ratio_numeric_vs_t_kin = numpy.array(ratio_numeric_vs_t_kin).flatten()
    ratio_analytic_vs_t_kin = numpy.array(ratio_analytic_vs_t_kin).flatten()

    relative_error = numpy.abs(
        1.0 - numpy.array(ratio_numeric_vs_t_kin) / numpy.array(
            ratio_analytic_vs_t_kin)
    )

    assert relative_error.max() < 1.2e-15
    assert relative_error.std() < 5.0e-16

def test_that_three_level_analytic_pop_dens_ratios_are_computed_correctly():

    #
    # load the two level species data
    #
    species_data = DataLoader().load('three_level_1')

    n_c = 1e6 * u.meter ** -3
    t_rad = 1000.0 * u.Kelvin

    #
    # for each value in the temperature range
    #  - compute the population density ration numerically
    #  - compute the population density analytically
    #
    # plot the comaprision
    #
    ratio_numeric_vs_t_kin = []
    ratio_analytic_vs_t_kin = []
    t_kin_range = numpy.logspace(1.0, 10.0, 50) * u.Kelvin
    for t_kin in t_kin_range:
        # population densities using matrix inversion
        pop_dens_eq_numeric = population_density_at_steady_state(
            species_data,
            t_kin,
            t_rad,
            n_c
        )

        pop_dens_eq_ratio_numeric = numpy.array(
            [
                pop_dens_eq_numeric[1] / pop_dens_eq_numeric[0],
                pop_dens_eq_numeric[2] / pop_dens_eq_numeric[0]
            ]).flatten()
        ratio_numeric_vs_t_kin.append(pop_dens_eq_ratio_numeric)

        # popultation densities using analytic expressions
        pop_dens_eq_ratio_analytic = analytic.population_density_ratio_three_level(
            species_data.energy_levels.data['g'],
            species_data.energy_levels.data['E'],
            species_data.k_dex_matrix_interpolator(t_kin)[1, 0],
            species_data.k_dex_matrix_interpolator(t_kin)[2, 0],
            species_data.k_dex_matrix_interpolator(t_kin)[2, 1],
            species_data.a_matrix[1, 0],
            species_data.a_matrix[2, 0],
            species_data.a_matrix[2, 1],
            n_c,
            t_kin,
            t_rad
        )

        ratio_analytic_vs_t_kin.append(numpy.array(pop_dens_eq_ratio_analytic))

    ratio_numeric_vs_t_kin = numpy.array(ratio_numeric_vs_t_kin)
    ratio_analytic_vs_t_kin = numpy.array(ratio_analytic_vs_t_kin)

    relative_error = numpy.abs(
        1.0 - numpy.array(ratio_numeric_vs_t_kin) / numpy.array(
            ratio_analytic_vs_t_kin)
    )

    assert relative_error.max() < 5.0e-15
    assert relative_error.std() < 6.0e-16

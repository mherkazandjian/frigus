# -*- coding: utf-8 -*-
"""
<keywords>
frigus, time, dependent, equilibrium, steady, state, compare, check
</keywords>
<description>
  - evolve the population density of the levels of a species using an integrator
  - compare the solution to the equilibrium solution using matrix inversion.
</description>
<seealso>
</seealso>
"""
from __future__ import print_function
import numpy

from astropy import units as u
import matplotlib.pyplot as plt
from scipy.integrate import ode

from frigus.population import (
    population_density_at_steady_state,
    compute_transition_rate_matrix
)
from frigus.readers.dataset import DataLoader

# load the species data
species_data_lique = DataLoader().load('H2_lique')
species_data_wrathmall = DataLoader().load('H2_wrathmall')

species_data = species_data_lique

# time range and timestep
t_0 = 0.0
t_f = 1e14
dt = 1e7

# the environment parameters
nc_h = 1e6 * u.meter ** -3
t_kin_array = numpy.arange(100.0, 5000., 100.)
t_kin_array = t_kin_array * u.Kelvin
t_rad = 0. * u.Kelvin

# the rates matrix (that is fixed if the environment params above are fixed)
# so it is computed once
#m_matrix = compute_transition_rate_matrix(species_data, t_kin, t_rad, nc_h)

y_td = []

for t_kin in t_kin_array:
    m_matrix = \
        compute_transition_rate_matrix(species_data, t_kin, t_rad, nc_h)


    def ode_rhs(_, y):
        """
        define the function that computes the rhs of dy/dt
        """
        return numpy.dot(m_matrix, y)


    def jac(t, y):
        return m_matrix

    def time_dependent(t_rad, t_kin, nc_h):
        # set the initial abundances for t = 0

        n_levels = len(species_data.energy_levels.data)

    #initial_fractional_abundances_wrathmall = numpy.ones(n_levels_wrathmall) \
    #                                      / n_levels_wrathmall

        initial_fractional_abundances = numpy.random.rand(n_levels)

        #initial_fractional_abundances[0] = 0.8
        #initial_fractional_abundances[1] = 1.e-3

        initial_fractional_abundances = \
            initial_fractional_abundances/initial_fractional_abundances.sum()

        y_0 = initial_fractional_abundances

        solver = ode(
            ode_rhs, jac=jac
        ).set_integrator(
            'lsoda',
            #method='bdf',
            with_jacobian = True,
            rtol=1e-14, atol=1e-14
        )

        solver.set_initial_value(y_0, t_0)

        step_counter = 0
        t_all = []
        x_all = []
        while solver.t < t_f:

            solver.integrate(solver.t + dt)

            t_all.append(solver.t)
            x_all.append(solver.y)
            rel_diff = None
            if len(x_all) >= 2:
                y_2 = x_all[-2]
                y_1 = x_all[-1]
                rel_diff = numpy.abs((y_2 - y_1)/y_1)
                if (rel_diff < 1e-8).all():
                    break

            if step_counter % 100 == 0:
                print('{} t={:e} 1-x.sum()={:e}'.format(
                    step_counter, solver.t, 1.0 - solver.y.sum()
                ))
                if rel_diff is not None:
                    print(numpy.where((rel_diff < 1e-6).all())[0].size, rel_diff.mean())
            step_counter += 1


        t_all, x_all = numpy.array(t_all), numpy.vstack(x_all)

        pop_dens_eq = population_density_at_steady_state(
            species_data,
            t_kin,
            t_rad,
            nc_h)

        return t_all, x_all, pop_dens_eq



    t_all, x_all, pop_dens_eq = time_dependent(t_rad, t_kin, nc_h)

    x = t_all
    y = x_all

    y_td.append(y[-1, :])
#
# plot the solution
#
fig, axs = plt.subplots()

colors = {
    'r': '-',
    'g': '--',
    'b': '-.',
    'k': '-..',
    'y': '+',
    'r--': '+',
    'g--': '+',
    'b--': '+',
    'c--': '+',
    'k--': '+'
}

for level_index in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]:

    x = t_kin_array.value
    y = y_td

    axs.loglog(
        x,
        y_td[:],
        # pop_dens_eq[level_index],
        'o', list(colors.keys())[level_index] + list(colors.values())[
         level_index],
         label='level {}'.format(level_index)
    )

    # axs.loglog(
    #        x[-1],
    #        y_all[:, level_index],
    #        #pop_dens_eq[level_index],
    #        'o', list(colors.keys())[level_index] + list(colors.values())[level_index],
    #        label='level {}'.format(level_index)
    #    )


    # axs.loglog(
    #     x,
    #     y,
    #     list(colors.keys())[level_index] + list(colors.values())[level_index],
    #     alpha=0.5,
    #     label='level {}'.format(level_index)
    # )


axs.set_xlabel('T [K]')
axs.set_ylabel('abundance')
#plt.legend()

plt.ion()
plt.show()

print('done')

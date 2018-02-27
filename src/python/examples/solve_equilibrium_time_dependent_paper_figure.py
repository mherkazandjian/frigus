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
species_data = DataLoader().load('HD_lipovka')

# time range and timestep
t_0 = 0.0
t_f = 1e10
dt = 1e5

# the environment parameters
nc_h = 1e8 * u.meter ** -3
t_kin = 2000.0 * u.Kelvin
t_rad = 2.73 * u.Kelvin

# the rates matrix (that is fixed if the environemnt params above are fixed)
# so it is computed once
m_matrix = compute_transition_rate_matrix(species_data, t_kin, t_rad, nc_h)


def ode_rhs(_, y):
    """
    define the function that computes the rhs of dy/dt
    """
    return numpy.dot(m_matrix, y)


# set the initial abundances for t = 0
n_levels = len(species_data.energy_levels.data)
initial_fractional_abundances = numpy.ones(n_levels) / n_levels
y_0 = initial_fractional_abundances


# define the solver
solver = ode(
    ode_rhs,
    jac=None
).set_integrator(
    'vode',
    method='bdf',
    # with_jacobian = False,
    rtol=1e-6
)

solver.set_initial_value(y_0, t_0)

step_counter = 0
t_all = []
x_all = []
while solver.t < t_f:
    solver.integrate(solver.t + dt)
    t_all.append(solver.t)
    x_all.append(solver.y)

    if step_counter % 100 == 0:
        print('{} t={:e} 1-x.sum()={:e}'.format(
            step_counter, solver.t, 1.0 - solver.y.sum()
        ))
    step_counter += 1


t_all, x_all = numpy.array(t_all), numpy.vstack(x_all)

#
# plot the solution
#
fig, axs = plt.subplots()

colors = {
    'k': '-',
    'k': '--',
    'k': '-.',
    'k': '-..',
    'k': '+',
    'r--': '+',
    'g--': '+',
    'b--': '+',
    'c--': '+',
    'k--': '+'
}

linestyles = [
    '-',
    '--',
    '-.',
    '-..',
]

# get the equilibrium solution and plot them as crosses at t = t_f
pop_dens_eq = population_density_at_steady_state(
    species_data,
    t_kin,
    t_rad,
    nc_h)

for level_index in [0, 1, 2]:

    x = t_all
    y = numpy.abs(1.0 - x_all[:, level_index] / pop_dens_eq[level_index])

    axs.loglog(
        list(x[0:200]) + list(x[200::100]),
        list(y[0:200]) + list(y[200::100]),
        'k' + linestyles[level_index],
        alpha=0.5,
        label='level {}'.format(level_index)
    )

axs.set_xlabel('t [s]')
axs.set_ylabel('relative difference')
plt.legend()

plt.ion()
plt.show()

print('done')

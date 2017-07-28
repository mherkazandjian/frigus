# -*- coding: utf-8 -*-
"""
<keywords>
frigus, time, dependent, equilibrium, steady, state, compare, check
</keywords>
<description>
evolve the population density of the levels of a species using an integrator
and compare the solution to the equilibrium solution using matrix inversion.
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
from frigus.readers import DataLoader

#
# load the species data
#
species_data = DataLoader().load('HD_lipovka')

# time range and timestep
t_0 = 0.0
t_f = 1e10
dt = 1e6

# the environment parameters
nc_H = 1e6 * u.meter ** -3
T_kin = 100.0 * u.Kelvin
T_rad = 2.73 * u.Kelvin

# the rates matrix (that is fixed if the environemnt params above are fixed)
# so it is computed once
M_matrix = compute_transition_rate_matrix(species_data, T_kin, T_rad, nc_H)


def ode_rhs(_, y, __):
    """
    define the function that computes the rhs of dy/dt
    """
    return numpy.dot(M_matrix, y)


# set the initial abundances for t = 0
n_levels = len(species_data.energy_levels.data)
initial_fractional_abundances = numpy.ones(n_levels) / n_levels
y_0 = initial_fractional_abundances


# define the solver
solver = ode(ode_rhs, jac=None).set_integrator('vode',
                                               method='bdf',
                                               # with_jacobian = False,
                                               rtol=1e-3)

solver.set_initial_value(y_0, t_0).set_f_params(1.0)

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

axs.loglog(t_all, x_all[:, 0], 'r')
axs.loglog(t_all, x_all[:, 1], 'g')
axs.loglog(t_all, x_all[:, 2], 'b')
axs.loglog(t_all, x_all[:, 3], 'k')

# get the equilibrium solution and plot them as crosses at t = t_f
pop_dens_eq = population_density_at_steady_state(
    species_data,
    T_kin,
    T_rad,
    nc_H)

axs.loglog(t_all[-1], pop_dens_eq[0], 'r+')
axs.loglog(t_all[-1], pop_dens_eq[1], 'g+')
axs.loglog(t_all[-1], pop_dens_eq[2], 'b+')
axs.loglog(t_all[-1], pop_dens_eq[3], 'k+')

plt.ion()
plt.show()

print('done')
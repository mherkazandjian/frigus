# -*- coding: utf-8 -*-
"""
calculate equilibrium population density ratio of a two level system and
plot it as a function of kinetic temperature
"""

import pylab
from numpy import fabs, logspace, array

from astropy import units as u

from frigus.solvers.analytic import population_denisty_ratio_two_level_no_radiation
from frigus.population import population_density_at_steady_state
from frigus.readers.dataset import DataLoader


# density of the colliding species, in m^3
nc = 1e14 * u.meter**-3

# the kinetic temperature of the gas
t_kin = 3000.0 * u.Kelvin
t_rad = 30.0 * u.Kelvin

species_data = DataLoader().load('two_level_1')


t_kin_range = logspace(2.0, 6.0, 100) * u.K

x_exact = array(
    [
        population_denisty_ratio_two_level_no_radiation(species_data, T_kin, nc)
        for T_kin in t_kin_range
    ], 'f8'
)

x_numerical = array(
    [
        population_density_at_steady_state(species_data, T_kin, 0.0, nc)
        for T_kin in t_kin_range
    ], 'f8'
)[:, :, 0]

# plot the relative population densities as a function of temperatures
pylab.figure()
pylab.loglog(t_kin_range, x_exact[:, 0], 'r-', label='x_0 exact')
pylab.loglog(t_kin_range, x_numerical[:, 0], 'ro', label='x_0 numerical')
pylab.loglog(t_kin_range, x_exact[:, 1], 'b-', label='x_1 exact')
pylab.loglog(t_kin_range, x_numerical[:, 1], 'bo', label='x_1 numerical')
pylab.legend()
pylab.title('fractional population densities of levels 0 and level 1')
pylab.show()

rel_errors = fabs(1.0 - x_numerical / x_exact)

assert rel_errors.max() < 1e-13

# plot the relative errors in the relative population densities as a
# function of temperatures
pylab.figure()
pylab.loglog(t_kin_range, rel_errors[:, 0], 'r-', label='x_0')
pylab.loglog(t_kin_range, rel_errors[:, 1], 'k-', label='x_1')
pylab.legend()
pylab.title('relative errors between numerical and exact solution')
pylab.show()

print('done')

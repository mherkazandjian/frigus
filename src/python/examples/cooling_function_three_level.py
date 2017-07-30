# -*- coding: utf-8 -*-
"""
<keywords>
frigus, cooling, funtion, two, level, analytic, test
</keywords>
<description>
compute the cooling function of a test two level system and compare it to
the analytical result
</description>
<seealso>
</seealso>
"""
from __future__ import print_function
import numpy
from astropy import units as u
import matplotlib.pyplot as plt

from frigus.readers import DataLoader
from frigus.population import (
    population_density_at_steady_state,
    population_density_ratio_analytic_three_level_system
)

#
# load the two level species data
#
species_data = DataLoader().load('three_level_1')

n_c = 1e6 * u.meter ** -3
T_kin = 100000.0 * u.Kelvin
T_rad = 100.0 * u.Kelvin

#
# for each value in the temperature range
#  - compute the population density ration numerically
#  - compute the population density analytically
#
# plot the comaprision
#
ratio_numeric_vs_T_kin = []
ratio_analytic_vs_T_kin = []
T_kin_range = numpy.logspace(1.0, 10.0, 50) * u.Kelvin
for T_kin in T_kin_range:

    # population densities using matrix inversion
    pop_dens_eq_numeric = population_density_at_steady_state(
        species_data,
        T_kin,
        T_rad,
        n_c)

    pop_dens_eq_ratio_numeric = numpy.array(
        [pop_dens_eq_numeric[1] / pop_dens_eq_numeric[0],
         pop_dens_eq_numeric[2] / pop_dens_eq_numeric[0]]).flatten()
    ratio_numeric_vs_T_kin.append(pop_dens_eq_ratio_numeric)

    # popultation densities using analytic expressions
    # pop_dens_eq_ratio_analytic = population_density_ratio_analytic_no_radiation(
    pop_dens_eq_ratio_analytic = population_density_ratio_analytic_three_level_system(
        species_data.energy_levels.data['g'],
        species_data.energy_levels.data['E'],
        species_data.K_dex_matrix_interpolator(T_kin)[1, 0],
        species_data.K_dex_matrix_interpolator(T_kin)[2, 0],
        species_data.K_dex_matrix_interpolator(T_kin)[2, 1],
        species_data.A_matrix[1, 0],
        species_data.A_matrix[2, 0],
        species_data.A_matrix[2, 1],
        n_c,
        T_kin,
        T_rad
    )

    ratio_analytic_vs_T_kin.append(numpy.array(pop_dens_eq_ratio_analytic))
    print(T_kin)

ratio_numeric_vs_T_kin = numpy.array(ratio_numeric_vs_T_kin)
ratio_analytic_vs_T_kin = numpy.array(ratio_analytic_vs_T_kin)

#
# plot the population density computed numericall and analyticall and also
# plot the relative difference between them
#
fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2)
ax0.loglog(T_kin_range, ratio_numeric_vs_T_kin[:, 0], 'r+', label = 'x_1/x_0')
ax0.loglog(T_kin_range, ratio_analytic_vs_T_kin[:, 0], 'k--')

ax0.loglog(T_kin_range, ratio_numeric_vs_T_kin[:, 1], 'b+', label = 'x_2/x_0')
ax0.loglog(T_kin_range, ratio_analytic_vs_T_kin[:, 1], 'g--')

ax0.set_xlabel('T_kin')
ax0.legend()


relative_error = []
relative_error.append(numpy.abs(
     1.0 - numpy.array(ratio_numeric_vs_T_kin[:,0]) /
     numpy.array(ratio_analytic_vs_T_kin[:,0])
))
relative_error.append(numpy.abs(
     1.0 - numpy.array(ratio_numeric_vs_T_kin[:,1]) /
     numpy.array(ratio_analytic_vs_T_kin[:,1])
))
relative_error = numpy.array(relative_error)


ax1.loglog(T_kin_range, relative_error[0,:],'b--', label = 'x_1/x_0')
ax1.loglog(T_kin_range, relative_error[1,:],'g--', label = 'x_2/x_0')
ax1.set_xlabel('T_kin')
ax1.set_title('relative error')
ax1.legend()

plt.show()

print('done')

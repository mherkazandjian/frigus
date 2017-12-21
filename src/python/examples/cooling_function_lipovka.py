# -*- coding: utf-8 -*-
"""
Dalculate equilibrium population of species and the cooling function using
the Lipovka data
"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit
from frigus.multivariate_fit import func, fit_lambda


from astropy import units as u

from frigus.population import (fit_lipovka,
                               cooling_rate_at_steady_state,
                               population_density_at_steady_state)

from frigus.readers import DataLoader

species_data = DataLoader().load('HD_lipovka')


x_fit = []
y_fit = []
data_to_fit = []

plt.ion()


# Calculate the population density and the cooling rate per particle
# for one value of kinetic temperature and gas density of H
if False:

    # density of the colliding species, in m^3
    nc_H = 1e6 * u.meter ** -3
    T_kin = 100.0 * u.Kelvin
    T_rad = 0.0 * u.Kelvin

    pop_dens_equ = population_density_at_steady_state(
        species_data,
        T_kin,
        T_rad,
        nc_H)

    cooling_rate = cooling_rate_at_steady_state(
        species_data,
        T_kin,
        T_rad,
        nc_H)

if True:

    # density of the colliding species, in m^3
    # nc_H = 1e6 * u.meter ** -3
    nc_H_rng = [1.e6, 1.e7, 1.e8, 1.e9, 1.e10, 1.e11, 1.e12, 1.e13, 1.e14] *\
               u.meter ** -3
    T_rad = 0.0 * u.Kelvin

    T_rng = np.logspace(2, 3.2, 10) * u.Kelvin
    # T_rng = species_data.raw_data.collision_rates_T_range
    for nc_H in nc_H_rng:
        lambda_vs_T_kin = []
        pop_dens_vs_T_kin = []
        for T_kin in T_rng:
            print(T_kin, nc_H)


            lambda_vs_T_kin = cooling_rate_at_steady_state(
                    species_data,
                    T_kin,
                    T_rad,
                    nc_H
            )

            x_fit.append(T_kin.value)
            y_fit.append(nc_H.cgs.value)
            lambda_vs_T_kin = u.Quantity(lambda_vs_T_kin)
            data_to_fit.append(lambda_vs_T_kin.cgs.value)



fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_zscale('log')

ax.scatter(x_fit, y_fit, data_to_fit)
ax.set_xlabel('Temperature [K]')
ax.set_ylabel('Density')
ax.set_zlabel('cooling function')

plt.show()

popt, pcov = fit_lambda(x_fit, y_fit, data_to_fit)

#lambda_flatten = [item for sublist in data_to_fit for item in sublist]

for i in np.arange(len(x_fit)):
    X_new = x_fit[i], y_fit[i]
    new_func = func(np.log10(X_new), *popt)
    new_lambda = 10 ** new_func
    lambda_vs_T_kin_lipovka = fit_lipovka(x_fit[i] * u.Kelvin,
                                          y_fit[i] * u.cm**-3)
    print(new_lambda, data_to_fit[i], lambda_vs_T_kin_lipovka.si.value)

print('done')

"""
compute the cooling function as a function of kinetic temperature and gas
density and plot them in 3D
"""
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt

from astropy import units as u

from frigus.population import cooling_rate_at_steady_state

from frigus.readers.dataset import DataLoader
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.gca(projection='3d')

species_data = DataLoader().load('HD_lipovka')

plt.ion()


T_rad = 0.0 * u.Kelvin

t_grid, density_grid = np.meshgrid(
    np.logspace(2, 3.2, 20) * u.Kelvin,
    np.logspace(6, 14, 20) * u.m ** -3
)

cooling_rate = np.zeros_like(t_grid.value).flatten()


for i, (t_kin, nc_H) in enumerate(zip(t_grid.flat, density_grid.flat)):
    cooling_rate[i] = cooling_rate_at_steady_state(
                        species_data,
                        t_kin,
                        T_rad,
                        nc_H).cgs.value

cooling_rate_grid = cooling_rate.reshape(t_grid.shape)


ax.plot_surface(
    np.log10(t_grid.value),
    np.log10(density_grid.value),
    np.log10(cooling_rate.reshape(t_grid.shape))
)

plt.setp(
    ax,
    'xlabel', 'log10(T_kin [K])',
    'ylabel', 'log10(density n [cm-3])'
)
ax.set_zlabel('log10(cooling function [cgs])')

plt.show()

print('done')

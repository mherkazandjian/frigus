"""
compute the cooling function as a function of kinetic temperature and gas
density and plot them in 2D, comparing with existing fit, if available
"""
from __future__ import print_function
import numpy as np
from astropy import units as u
from frigus.readers.dataset import DataLoader
from frigus.cooling_function.grid import CoolingFunctionGrid


grid = CoolingFunctionGrid()

grid.set_species(DataLoader().load('HD_lipovka'))

grid.set_density(np.logspace(6, 14, 10) * u.m ** -3)
grid.set_t_kin(np.logspace(2, 3.2, 10) * u.Kelvin)
grid.set_t_rad(0.0 * u.Kelvin)

grid.plot_2d(x='n', y='t_kin', show=True)

print('done')

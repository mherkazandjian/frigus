"""
compute the cooling function as a function of kinetic temperature and gas
density and plot them in 3D
"""
from __future__ import print_function
import numpy as np
from astropy import units as u
from frigus.cooling_function.grid import CoolingFunctionGrid
from frigus.readers.dataset import DataLoader

grid = CoolingFunctionGrid()

grid.set_species(DataLoader().load('HD_lipovka'))

grid.set_density(np.logspace(6, 14, 10) * u.m ** -3)
grid.set_t_kin(np.logspace(2, 3.2, 10) * u.Kelvin)
grid.set_t_rad(0.0 * u.Kelvin)

grid.plot(x='n', y='t_kin')

print('done')

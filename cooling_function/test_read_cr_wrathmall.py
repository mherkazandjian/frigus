from __future__ import print_function
import os

import read_cr_wrathmall
from read_collision_coefficients import read_collision_coefficients

fname = os.path.expanduser('~/dynamics/cooling_function/fortran/Read/wrathmall/Rates_H_H2_flower_new.dat')

x = read_collision_coefficients(fname)

print('done')

from __future__ import print_function
import os

import read_cr_wrathmall
from read_collision_coefficients import read_collision_coefficients

fname = os.path.expanduser('~/dynamics/cooling_function/fortran/Read/wrathmall/Rates_H_H2_flower_new.dat')

x = read_collision_coefficients(fname)

print('done')
reader = read_cr_wrathmall.Reader('Read/wrathmall/H.ortho_H2')

'''
# definig the reader object
reader = uga('/home/mher/tmp/foo/H2XBCphotodiss.cs')

# showing info
reader.info()

# get the data of a certain transition
en, cs = reader.data[4,1,2]
'''

print('done')
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 21 23:45:39 2015

@author: carla
"""

import read_ei
import read_cr
import read_levels
from scipy.constants import c,h,Boltzmann

en_H2 = read_levels.read_levels("Read/H2Xvjlevels.cs")
cr, T, ini, fin, vj_unique = read_cr.read_coeff("Read/Rates_H_H2.dat")
A = read_ei.read_einstein()

kB = Boltzmann

print A,c,h,kB,en_H2

print vj_unique
print 'done reading!'


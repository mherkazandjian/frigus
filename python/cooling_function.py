# -*- coding: utf-8 -*-
"""
Created on Sat Nov 21 23:45:39 2015

@author: carla
"""

import read_ei
import read_cr
import read_levels
#import population
import matplotlib.pyplot as plt

en_H2 = read_levels.read_levels("Read/H2Xvjlevels.cs")
cr, T, v, j, vp, jp, ini, fin, vj_unique = read_cr.read_coeff("Read/Rates_H_H2.dat")
A = read_ei.read_einstein()

#matrix = population.computeRateMatrix()
#nvj = population.solveEquilibrium()

#print A,c,h,kB,en_H2

print vj_unique
print 'done reading!'

plt.yscale('log')
plt.xscale('log')
plt.plot(T, cr[ini[:,1][0],ini[:,1][1],fin[:,1][0],fin[:,1][1]])
#plt.show()
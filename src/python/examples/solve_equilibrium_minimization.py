# -*- coding: utf-8 -*-
"""
calculate equilibrium population of species and the cooling function.
"""
from __future__ import print_function
import pylab
import numpy

from astropy import units as u
import matplotlib.pyplot as plt

from population import population_density_at_steady_state
import utils

from readers import DataLoader

species_data = DataLoader().load('H2_lique')
# species_data = DataLoader().load('HD_lipovka')

# solving using minimization
full2 = computeRateMatrix(pNH3, Tkin, nc)

# generating the constraint equation and appending it to the Matrix
cons= ones(pNH3.nlevels)
full2 = vstack((full2, cons))

# generating the RHS
rhs = zeros(pNH3.nlevels + 1)
rhs[-1] = 1

# solving
sol = linalg.lstsq(full2, rhs)
f2 = sol[0]

print('done')
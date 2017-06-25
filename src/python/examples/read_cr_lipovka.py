"""
<keywords>
test, python, collision, rates, lipovka, wrathmal, roueff
</keywords>
<description>
read the collision data by flower and roueff that are used by lipovka.

.. todo:: the raw data can be read and parsed but the reduction into matrix
 form needs to be implemented.
</description>
<seealso>
</seealso>
"""
from __future__ import print_function
import pylab
import numpy
import pdb

from astropy import units as u
import matplotlib.pyplot as plt

import population
from population import (fit_glover,
                        cooling_rate_at_steady_state,
                        population_density_at_steady_state)
import utils

from readers import DataLoader

species_data = DataLoader().load('HD_lipovka')
# species_data = DataLoader().load('H2_lique')

# from __future__ import print_function
# from read_cr_lipovka import Reader
#
# reader = Reader('../../../data/read/lipovka/flower_roueff_data.dat')
#
# reader.read_data()
print('done')

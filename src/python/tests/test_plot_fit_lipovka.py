from __future__ import print_function
import numpy
from astropy import units as u
from population import fit_lipovka
import pylab
n_hd = 1e6 * u.m**-3
T = numpy.linspace(100.0, 2000.0, 100) * u.K

cooling_function = fit_lipovka(T, n_hd)

pylab.loglog(T, cooling_function, '-o')
pylab.xlabel('T [{}]'.format(str(T.unit)))
pylab.ylabel('Cooling function [{}]'.format(str(cooling_function.unit)))
pylab.show()

print('done')
from __future__ import print_function
import os
import numpy

fpath = os.path.expanduser('~/tmp/sandbox/to_be_fitted.dat')

data = numpy.loadtxt(fpath).T

_, x, y, f_xy = data

from matplotlib import pyplot
from scipy.interpolate import griddata

x = numpy.log10(x)
y = numpy.log10(y)
f_xy = numpy.log10(f_xy)

# define the data to be fit to be in a form that griddata expects
points = numpy.vstack((x, y)).T
values = f_xy

# defining the points in the uniform 2D grid (regular grid where the interpolation will be done)
grid_x, grid_y = numpy.meshgrid(numpy.unique(x), numpy.unique(y))

# interpolating with different methods
grid_z0 = griddata(points, values, (grid_x, grid_y), method='nearest')
grid_z1 = griddata(points, values, (grid_x, grid_y), method='linear')
grid_z2 = griddata(points, values, (grid_x, grid_y), method='cubic')

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.plot_surface(grid_x, grid_y, grid_z2,
                cmap=cm.coolwarm,
                linewidth=0, antialiased=True)

"""
minimize

      W( f*(x, y, {p}), f(x,y)| )

where f*({p}) is the analytical function with parameters {p}. where

      W = (1/(2*N) * sqrt( \sum_i,j (f*(x_i,y_j, {p})^2 - f(x_i, y_j))^2 )

or W could have another form

and The analytic function f*(x, y, {p}) is

      f*{p} = \sum_l,m=0^l,m=4 p_lm * x^l * y^m

so the minimizer would explore values of {p} each time evaluating f*{p} at
the locationx x, y of the data and taking its and minimize W

"""
from scipy.optimize import minimize

pyplot.show()

# pyplot.scatter(numpy.log10(x) ,numpy.log10(y))
print('done')
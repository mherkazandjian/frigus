import numpy
from numpy import isscalar

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from frigus.population import cooling_rate_at_steady_state


class CoolingFunctionGrid(object):
    def __init__(self):
        self.species = None
        self.n = None
        self.t_kin = None
        self.t_rad = None
        self.n_grid = None
        self.t_kin_grid = None
        self.t_rad_grid = None
        self.cooling_function = None

    def set_species(self, species):
        self.species = species

    def set_density(self, density):
        self.n = density

    def set_t_kin(self, t_kin):
        self.t_kin = t_kin

    def set_t_rad(self, t_rad):
        self.t_rad = t_rad


    def _compute_mesh(self):

        n_grid, t_kin_grid, t_rad_grid = numpy.meshgrid(
            self.n,
            self.t_kin,
            self.t_rad
        )

        # check the shape, one of the three quantities should be a scalar
        self.n_grid = n_grid.squeeze()
        self.t_kin_grid = t_kin_grid.squeeze()
        self.t_rad_grid = t_rad_grid.squeeze()

        msg = (
            'one of the three quantities should have either a dimension 1 or\n'
            'should be a scalar. The provided shapes are:\n'
            'n      = {n_shape}\n'
            't_kin  = {t_kin_shape}\n'
            't_rad  = {t_rad_shape}\n'
        ).format(
            n_shape=1 if isscalar(self.n) else self.n.shape,
            t_kin_shape=1 if isscalar(self.t_kin.shape) else self.t_kin.shape,
            t_rad_shape=1 if isscalar(self.t_rad.shape) else self.t_rad.shape
        )

        assert len(self.n_grid.shape) == 2, msg
        assert len(self.t_kin_grid.shape) == 2, msg
        assert len(self.t_rad_grid.shape) == 2, msg

    def compute(self):

        self._compute_mesh()

        cooling_rate = numpy.zeros_like(self.n_grid.value).flatten()
        for i, (n, t_kin, t_rad) in enumerate(
                zip(self.n_grid.flat,
                    self.t_kin_grid.flat,
                    self.t_rad_grid.flat)):

            cooling_rate[i] = cooling_rate_at_steady_state(
                self.species,
                t_kin,
                t_rad,
                n).cgs.value

        cooling_rate_grid = cooling_rate.reshape(self.n_grid.shape)
        self.cooling_function = cooling_rate_grid
        return cooling_rate_grid

    def _determine_x_y_quantities(self, x, y):
        xval = getattr(self, '{}_grid'.format(x))
        yval = getattr(self, '{}_grid'.format(y))
        return xval, yval

    def plot(self, x='n', y='t_kin', show=True):

        if self.cooling_function is None:
            self.compute()

        xval, yval = self._determine_x_y_quantities(x, y)
        assert xval.shape == yval.shape, 'requested quantities are not meshes'

        zval = self.cooling_function

        fig = plt.figure()
        ax = fig.gca(projection='3d')

        ax.plot_surface(
            numpy.log10(xval.value),
            numpy.log10(yval.value),
            numpy.log10(zval)
        )

        plt.setp(
            ax,
            'xlabel', 'log10(T_kin [K])',
            'ylabel', 'log10(density n [cm-3])'
        )
        ax.set_zlabel('log10(cooling function [cgs])')

        if show is True:
            plt.show()
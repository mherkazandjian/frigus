import numpy
from numpy import isscalar

from astropy import units as u
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from frigus.population import cooling_rate_at_steady_state
from frigus.cooling_function.fits import fit_lipovka


class CoolingFunctionGrid(object):
    """
    Compute the cooling function over a grid of densities and temperatures
    """
    def __init__(self):
        """
        Constructor
        """

        self.species = None
        """the specieis with the radiative and collisional data"""

        self.n = None
        """the gas density mesh points"""

        self.t_kin = None
        """the kinetic temperature mesh points"""

        self.t_rad = None
        """the radiation temperature mesh point"""

        self.n_grid = None
        """the gas density grid"""

        self.t_kin_grid = None
        """the kinetic temperature grid"""

        self.t_rad_grid = None
        """The radiation temperature grid"""

        self.cooling_function = None
        """The cooling function grid"""

    def set_species(self, species):
        """setter for the speicies object"""
        self.species = species

    def set_density(self, density):
        """setter for the gas density"""
        self.n = density

    def set_t_kin(self, t_kin):
        """setter for the kinetic temperature"""
        self.t_kin = t_kin

    def set_t_rad(self, t_rad):
        """setter for the radiation temperature"""
        self.t_rad = t_rad

    def _compute_mesh(self):
        """
        Compute the 3D grid where the cooling function will be computed
        """
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
        """
        Compute the cooling function for the specified grid

        :return: ndarray
        """

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
        """
        Get the attribute values by specifiying the quantity names

        :param str x: the x attribute
        :param str y: the y attribute
        :return: tuple: the grid value of x and y
        """
        xval = getattr(self, '{}_grid'.format(x))
        yval = getattr(self, '{}_grid'.format(y))
        return xval, yval

    def plot_3d(self, x='n', y='t_kin', show=True):
        """
        Plot the cooling function in 3d

        :param str x: the quantity of the x axis
        :param str y: the quantity of the y axis
        :param bool show: display the plot if True
        """

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
            'xlabel', 'log$_{10}$(density n [m$^{-3}$])',
            'ylabel', 'log$_{10}$(T$_{kin}$ [K])'
        )
        ax.set_zlabel('log$_{10}$(cooling function [cgs])')

        if show is True:
            plt.show()

    def plot_2d(self, x='n', y='t_kin', show=True):
        """
        Plot the cooling function in 2d

        :param str x: the quantity of the x axis
        :param str y: the quantity of the y axis
        :param bool show: display the plot if True
        """
        if self.cooling_function is None:
            self.compute()

        xval, yval = self._determine_x_y_quantities(x, y)

        plt.ion()
        fig, axs = plt.subplots(figsize=(8, 8))

        plot_markers = [(2 + i // 2, 1 + i % 2, 0) for i in range(16)]

        lambda_vs_t_kin = u.Quantity(self.cooling_function)

        lambda_vs_t_kin_lipovka = fit_lipovka(yval[:], xval[:])

        for nc_index, nc_h in enumerate(x):

            axs.loglog(
                yval.value, lambda_vs_t_kin.cgs.value,
                '-x', color='black', marker=plot_markers[nc_index], label=nc_h
            )

            axs.loglog(
                yval.value, lambda_vs_t_kin_lipovka.cgs.value,
                'r--', color='black', label=''
            )

        axs.set_xlabel('T$_\mathrm{kin}$ [K]')
        axs.set_ylabel('cooling function [erg s$^{-1}$]')

        if show is True:
            plt.show()

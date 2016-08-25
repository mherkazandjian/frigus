"""
module that provides functions to load various data in a standard format.
"""

import read_energy_levels
import read_einstien_coefficient
from read_collision_coefficients import read_collision_coefficients
import population


class DataSetRawBase(object):
    """

    """
    def __init__(self):
        """
        constructor
        """

        self.energy_levels = None
        """The energy levels"""

        self.A = None
        """The A transition rates"""

        self.A_info_nnz = None
        """The non zero info of self.A"""

        self.collision_rates = None
        """The collision rates"""

        self.collision_rates_T_range = None
        """The range over which the collision rates are provided"""

        self.collision_rates_info_nnz = None
        """The non zero info of the collision rates"""



class DataSetBase(object):
    """
    An abstract class than can be sub-classed to populate the data attributes.
    This data set is the one that would be used to do actual computations of
    e.g. time evolution on equilibrium solutions...etc..
    """
    def __init__(self):
        """
        constructor
        """
        self.energy_levels = None
        """The energy levels object"""

        self.A_matrix = None
        """The matrix of the Einstein coefficients"""

        self.K_dex_matrix_interpolator = None
        """An interpolation function that takes T_kin as an argument and
        returns an array of the same shape as self.A_matrix"""

        self.raw_data = DataSetRawBase()
        """the raw data from which the 2D matrices are computed"""

    def read_raw_data(self):
        """populate the self.raw_data object"""
        pass

    def reduce_raw_data(self):
        """use the raw data to produce the 2D matrices"""


class DataSetH2Lique(DataSetBase):
    """
    Data of H2 colliding with H using collisional data by F. Lique.

      - energy levels of H2 (vibrational and rotational)
      - collisional coefficients of H2 with H (K_ij)
      - radiative coefficients (A_ij, B_ij, B_ji)

    Limitations

      - The smallest data set of (A, B, K) determines the number of states to
       be inserted in the model.
    """
    def __init__(self):
        """
        constructor
        """
        super(DataSetH2Lique, self).__init__()
        self.read_raw_data()
        self.reduce_raw_data()

    def read_raw_data(self):
        """read the raw H2 data"""

        #
        # read the energy levels (v, j, energy)
        #
        energy_levels = read_energy_levels.read_levels_lique(
                                         'Read/H2Xvjlevels_francois_mod.cs')

        # en_H2 = read_levels.read_levels("Read/H2Xvjlevels.cs")
        # print('{:3}{:3}{:10}'.format('v', 'j', 'E(eV)'))
        # for level in levels:
        #     print('{:<3}{:<3}{:<10}'.format(level['v'], level['j'], level['E']))
        self.raw_data.energy_levels = energy_levels

        #
        # read the einstein coefficients for the H2 transitions
        #
        A, A_info_nnz = read_einstien_coefficient.read_einstein()
        self.raw_data.A = A
        self.raw_data.A_info_nnz = A_info_nnz

        # read the collisional rates for H2 with H
        collision_rates, T_rng, collision_rates_info_nnz =\
            read_collision_coefficients("Read/Rates_H_H2.dat")
        self.raw_data.collision_rates = collision_rates
        self.raw_data.collision_rates_T_range = T_rng
        self.raw_data.collision_rates_info_nnz = collision_rates_info_nnz

    def reduce_raw_data(self):
        """
        """

        # find the maximum v and j from the Einstein and collisional rates data
        # sets and adjust the labels of the energy levels according to that
        v_max_data, j_max_data = population.find_v_max_j_max_from_data(
            self.raw_data.A_info_nnz,
            self.raw_data.collision_rates_info_nnz)

        self.energy_levels = self.raw_data.energy_levels
        self.energy_levels.set_labels(v_max=v_max_data + 1)

        #
        # reduce the Einstein coefficients to a 2D matrix (construct the A
        # matrix) [n_levels, n_levels]
        # A_reduced_slow = population.reduce_einstein_coefficients_slow(
        #                                 self.raw_data.A_info_nnz,
        #                                 self.energy_levels)
        A_matrix = population.reduce_einstein_coefficients(
                          self.raw_data.A,
                          self.energy_levels)
        self.A_matrix = A_matrix

        # getting the collisional de-excitation matrix (K_dex) (for all
        # tabulated values)  [n_level, n_level, n_T_kin_values]
        K_dex_matrix = population.reduce_collisional_coefficients_slow(
                             self.raw_data.collision_rates_info_nnz,
                             self.energy_levels)

        # compute the interpolator that produces K_dex at a certain temperature
        K_dex_matrix_interpolator = population.compute_K_dex_matrix_interpolator(
            K_dex_matrix, self.raw_data.collision_rates_T_range)
        self.K_dex_matrix_interpolator = K_dex_matrix_interpolator


class DataLoader(object):
    """
    Load various data sets.
    """
    def __init__(self):
        """

        """
        pass

    def load(self, name):
        """
        :param string name: The name of the data set to be loaded
        :return:
        """
        if name == 'H2_lique':
            return DataSetH2Lique()

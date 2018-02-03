"""
module that provides functions to load various data in a standard format.
No parsing is done in this module. The parsing is done in other modules
that implement the specific readers and parsers for e.g. einstien coefficient,
energy levels and collision reaction rates.
"""
import os
import numpy
from astropy import units as u

from frigus import utils
from frigus import read_energy_levels, read_einstein_coefficient
from frigus.read_collision_coefficients import (
    read_collision_coefficients_lique_and_wrathmall,
    read_collision_coefficients_lipovka
)

from frigus import population

DATADIR = utils.datadir_path()


class DataSetRawBase(object):
    """
    Abstract class for storing raw data of species. 
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
        raise NotImplementedError("this method should be implemented by "
                                  "the subclass")

    def reduce_raw_data(self):
        """use the raw data to produce the 2D matrices"""
        raise NotImplementedError("this method should be implemented by "
                                  "the subclass")


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
            os.path.join(DATADIR, 'H2Xvjlevels_francois_mod.cs')
        )

        self.raw_data.energy_levels = energy_levels

        #
        # read the einstein coefficients for the H2 transitions
        #
        A, A_info_nnz = read_einstein_coefficient.read_einstein_simbotin()
        self.raw_data.A = A
        self.raw_data.A_info_nnz = A_info_nnz

        #
        # read the collisional rates for H2 with H
        #
        collision_rates, T_rng, collision_rates_info_nnz = \
            read_collision_coefficients_lique_and_wrathmall(
                os.path.join(DATADIR, "Rates_H_H2.dat")
            )
        self.raw_data.collision_rates = collision_rates
        self.raw_data.collision_rates_T_range = T_rng
        self.raw_data.collision_rates_info_nnz = collision_rates_info_nnz

    def reduce_raw_data(self):
        """
        use the raw data in self.raw_data to populate the A and K_dex matrices
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
            self.energy_levels,
            reduced_data_is_upper_to_lower_only=True
        )

        # compute the interpolator that produces K_dex at a certain temperature
        K_dex_matrix_interpolator = population.compute_K_dex_matrix_interpolator(
            K_dex_matrix, self.raw_data.collision_rates_T_range)
        self.K_dex_matrix_interpolator = K_dex_matrix_interpolator


class DataSetH2Wrathmall(DataSetBase):
    """
    Data of H2 colliding with H using collisional data by Wrathmall and Flower

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
        super(DataSetH2Wrathmall, self).__init__()
        self.read_raw_data()
        self.reduce_raw_data()

    def read_raw_data(self):
        """read the raw H2 data"""

        #
        # read the energy levels (v, j, energy)
        #
        energy_levels = read_energy_levels.read_levels_wrathmall_and_flower(
            os.path.join(DATADIR, 'H2Xvjlevels_flower.cs')
        )

        self.raw_data.energy_levels = energy_levels

        #
        # read the einstein coefficients for the H2 transitions
        #
        A, A_info_nnz = read_einstein_coefficient.read_einstein_simbotin()
        self.raw_data.A = A
        self.raw_data.A_info_nnz = A_info_nnz

        #
        # read the collisional rates for H2 with H
        #
        collision_rates, T_rng, collision_rates_info_nnz = \
            read_collision_coefficients_lique_and_wrathmall(
                os.path.join(
                    DATADIR,
                    'wrathmall', 'Rates_H_H2_flower_frigus_downwards.dat'
                )
            )

        self.raw_data.collision_rates = collision_rates
        self.raw_data.collision_rates_T_range = T_rng
        self.raw_data.collision_rates_info_nnz = collision_rates_info_nnz

    def reduce_raw_data(self):
        """
        use the raw data in self.raw_data to populate the A and K_dex matrices
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
            self.energy_levels,
            reduced_data_is_upper_to_lower_only=True
        )

        # compute the interpolator that produces K_dex at a certain temperature
        K_dex_matrix_interpolator = population.compute_K_dex_matrix_interpolator(
            K_dex_matrix, self.raw_data.collision_rates_T_range)
        self.K_dex_matrix_interpolator = K_dex_matrix_interpolator


class DataSetH2Glover(DataSetBase):
    """
    Data of H2 colliding with H using collisional data by Wrathmall and Flower,
    including only the low energy levels of H2 (to compare with the fit by Glover)

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
        super(DataSetH2Glover, self).__init__()
        self.read_raw_data()
        self.reduce_raw_data()

    def read_raw_data(self):
        """read the raw H2 data"""

        #
        # read the energy levels (v, j, energy)
        #
        energy_levels = read_energy_levels.read_levels_wrathmall_and_flower(
            os.path.join(DATADIR, 'H2Xvjlevels_low_energies.dat')
        )

        self.raw_data.energy_levels = energy_levels

        #
        # read the einstein coefficients for the H2 transitions
        #
        A, A_info_nnz = read_einstein_coefficient.read_einstein_simbotin()
        self.raw_data.A = A
        self.raw_data.A_info_nnz = A_info_nnz

        #
        # read the collisional rates for H2 with H
        #
        collision_rates, T_rng, collision_rates_info_nnz = \
            read_collision_coefficients_lique_and_wrathmall(
                os.path.join(
                    DATADIR,
                    "wrathmall",
                    "Rates_H_H2_flower_low_energy_frigus_downwards.dat"
                )
            )

        self.raw_data.collision_rates = collision_rates
        self.raw_data.collision_rates_T_range = T_rng
        self.raw_data.collision_rates_info_nnz = collision_rates_info_nnz

    def reduce_raw_data(self):
        """
        use the raw data in self.raw_data to populate the A and K_dex matrices
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
            self.energy_levels,
            reduced_data_is_upper_to_lower_only=True
        )

        # compute the interpolator that produces K_dex at a certain temperature
        K_dex_matrix_interpolator = population.compute_K_dex_matrix_interpolator(
            K_dex_matrix, self.raw_data.collision_rates_T_range)
        self.K_dex_matrix_interpolator = K_dex_matrix_interpolator


class DataSetTwoLevel_1(DataSetBase):
    """
    Data of a synthetic two level system
    """
    def __init__(self):
        """
        constructor. The raw data is already reduced and in usable form without
        the need to reducing it.
        """
        super(DataSetTwoLevel_1, self).__init__()
        self.read_raw_data()

    def read_raw_data(self):
        """read the raw two level system data"""

        #
        # read the energy levels (v, j, energy)
        #
        energy_levels = read_energy_levels.read_levels_two_levels_test_1(
            os.path.join(
                DATADIR, 'two_levels_1', 'energy_levels.txt')
        )

        self.energy_levels = energy_levels
        self.raw_data.energy_levels = energy_levels

        #
        # read the einstein coefficients
        #
        A_matrix = numpy.loadtxt(
            os.path.join(DATADIR, 'two_levels_1', 'A_matrix.txt')
        )

        A_matrix = A_matrix / u.second
        self.A_matrix = A_matrix
        self.raw_data.A = self.A_matrix / u.second

        #
        # read the collisional rates
        #
        collision_rates = numpy.loadtxt(
            os.path.join(DATADIR, 'two_levels_1', 'K_dex_matrix.txt')
        )

        collision_rates = collision_rates * u.m**3 / u.second
        T_rng = 0.0 * u.K, 1e6 * u.K

        self.K_dex_matrix_interpolator = lambda T_kin: collision_rates

        self.raw_data.collision_rates = collision_rates
        self.raw_data.collision_rates_T_range = T_rng


class DataSetThreeLevel_1(DataSetBase):
    """
    Data of a synthetic three level system
    """
    def __init__(self):
        """
        constructor. The raw data is already reduced and in usable form without
        the need to reducing it.
        """
        super(DataSetThreeLevel_1, self).__init__()
        self.read_raw_data()

    def read_raw_data(self):
        """read the raw three level system data"""

        #
        # read the energy levels (v, j, energy)
        #
        energy_levels = read_energy_levels.read_levels_two_levels_test_1(
            os.path.join(DATADIR, 'three_levels_1', 'energy_levels.txt')
        )

        self.energy_levels = energy_levels
        self.raw_data.energy_levels = energy_levels

        #
        # read the einstein coefficients
        #
        A_matrix = numpy.loadtxt(
            os.path.join(DATADIR, 'three_levels_1', 'A_matrix.txt')
        )

        A_matrix = A_matrix / u.second
        self.A_matrix = A_matrix
        self.raw_data.A = self.A_matrix / u.second

        #
        # read the collisional rates
        #
        collision_rates = numpy.loadtxt(
            os.path.join(DATADIR, 'three_levels_1', 'K_dex_matrix.txt')
        )

        collision_rates = collision_rates * u.m**3 / u.second
        T_rng = 0.0 * u.K, 1e6 * u.K

        self.K_dex_matrix_interpolator = lambda T_kin: collision_rates

        self.raw_data.collision_rates = collision_rates
        self.raw_data.collision_rates_T_range = T_rng


class DataSetHDLipovka(DataSetBase):
    """
    Data of H2 colliding with H using collisional data use by Lipovka.

      - energy levels of HD (only rotational for v = 0)
      - collisional coefficients of HD with H (K_ij)
      - radiative coefficients (A_ij, B_ij, B_ji)

    Limitations

      - The smallest data set of (A, B, K) determines the number of states to
       be inserted in the model.
    """
    def __init__(self):
        """
        constructor
        """
        super(DataSetHDLipovka, self).__init__()
        self.read_raw_data()
        self.reduce_raw_data()

    def read_raw_data(self):
        """read the raw HD data"""

        #
        # read the energy levels (v, j, energy)
        #
        energy_levels = read_energy_levels.read_levels_lipovka(
            os.path.join(DATADIR, 'lipovka', 'flower_roueff_data.dat')
        )

        self.raw_data.energy_levels = energy_levels

        #
        # read the einstein coefficients for the HD transitions
        #
        A, A_info_nnz = read_einstein_coefficient.read_einstein_coppola()
        self.raw_data.A = A
        self.raw_data.A_info_nnz = A_info_nnz

        #
        # read the collisional rates for HD with H
        #
        collision_rates, T_rng, collision_rates_info_nnz = \
            read_collision_coefficients_lipovka(
                os.path.join(
                    DATADIR, 'lipovka', 'flower_roueff_data.dat'
                )
            )

        self.raw_data.collision_rates = collision_rates
        self.raw_data.collision_rates_T_range = T_rng
        self.raw_data.collision_rates_info_nnz = collision_rates_info_nnz

    def reduce_raw_data(self):
        """
        use the raw data in self.raw_data to populate the A and K_dex matrices
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
            self.energy_levels,
            set_inelastic_coefficient_to_zero=True,
            set_excitation_coefficients_to_zero=True,
            reduced_data_is_upper_to_lower_only=False,
        )

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
        elif name == 'HD_lipovka':
            return DataSetHDLipovka()
        elif name == 'H2_wrathmall':
            return DataSetH2Wrathmall()
        elif name == 'H2_low_energy_levels':
            return DataSetH2Glover()
        elif name == 'two_level_1':
            return DataSetTwoLevel_1()
        elif name == 'three_level_1':
            return DataSetThreeLevel_1()
        else:
            msg = '''not data loader defined for {}'''.format(name)
            raise ValueError(msg)

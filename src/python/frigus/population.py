# -*- coding: utf-8 -*-

#    population.py is part of Frigus.

#    Frigus: software to compure the energy exchange in a multi-level system
#    Copyright (C) 2016-2018 Mher V. Kazandjian and Carla Maria Coppola

#    Frigus is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, version 3 of the License.    
#
#    Frigus is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with Frigus.  If not, see <http://www.gnu.org/licenses/>.


from __future__ import print_function

import numpy
from numpy import zeros, fabs, exp, where, in1d, pi

import scipy
from scipy import interpolate


from astropy import units as u
from astropy.constants import c as c_light
from astropy.constants import h as h_planck
from astropy.constants import k_B as kb
from astropy.modeling.blackbody import blackbody_nu as B_nu

from frigus.utils import linear_2d_index, find_matching_indices
from frigus.solvers.linear import solve_equilibrium


def find_v_max_j_max_from_data(a_einstein_nnz, cr_coefficients_nnz):
    """
    Find the v_max and j_max from the Einstien and collisional coefficients data

    The coeffients are for a two quantum number system (v,j).

    :param tuple a_einstein_nnz: The Einstien coefficients for the transitions
     provided as numpy arrays. (See DataSetRaw.A_info_nnz)
    :param tuple cr_coefficients_nnz: The collisional coeffieints for the
     transitions. (See DataSetRaw.collision_rates_info_nnz)
    :return: tuple: v_max, j_max
    """

    # find the non zeros elements and their corresponding indices of the
    # Einstein coefficients
    v_nnz, j_nnz, vp_nnz, jp_nnz, a_nnz = a_einstein_nnz
    v_max_a = max(v_nnz.max(), vp_nnz.max())
    j_max_a = max(j_nnz.max(), jp_nnz.max())
    # print('DEBUG: (v_max_A, j_max_A) = ', v_max_A, j_max_A)

    # find the non zeros elements and their corresponding indices of the
    # collisional coefficients
    (v_nnz, j_nnz), (vp_nnz, jp_nnz), unique_nnz, cr_nnz = cr_coefficients_nnz
    v_max_cr = max(v_nnz.max(), vp_nnz.max())
    j_max_cr = max(j_nnz.max(), jp_nnz.max())
    # print('DEBUG: (v_max_cr, j_max_cr) = ', v_max_cr, j_max_cr)

    return max(v_max_a, v_max_cr), max(j_max_a, j_max_cr)


def check_self_transitions_in_einstien_nnz_data(a_info_nnz):
    """
    Raise an error if there are self transitions in the Einstein coefficient

    :param DataSetRawBase.A_info_nnz like a_info_nnz: a tuple of elements v_nnz,
     j_nnz, vp_nnz, jp_nnz, A_nnz
    :return: True if everything ok, else error is raised
    """
    v_nnz, j_nnz, vp_nnz, jp_nnz, a_nnz = a_info_nnz

    for i, A in enumerate(a_nnz):
        v, j, vp, jp = v_nnz[i], j_nnz[i], vp_nnz[i], jp_nnz[i]
        if v == vp and j == jp:
            raise ValueError('v -> v, j -> j transition. This is not possible')
    else:
        return True


def reduce_einstein_coefficients_slow(a_info_nnz, energy_levels):
    """
    Construct the A_matrix from the nnz info of the Einstein coefficients

    Use the nonzero entries of the Einstein coefficients to construct the
    A matrix for the levels found in the energy_levels parameter.

    :param DataSetRawBase.A_info_nnz like a_info_nnz: a tuple of elements v_nnz,
     j_nnz, vp_nnz, jp_nnz, A_nnz
    :param energy_levels: numpy array array of N rows, the v and j levels and
     their corresponding energy. The elements of the array are stored in
     increasing order in the energy. This array can be obtained e.g. from
     read_levels.read_levels_lique()
    :return: A 2D matrix of shape (energy_levels.size, energy_levels.size) that
     is a lower triangular matrix containing the Einstein coefficients.
    """

    check_self_transitions_in_einstien_nnz_data(a_info_nnz)

    v_nnz, j_nnz, vp_nnz, jp_nnz, a_nnz = a_info_nnz

    levels = energy_levels

    # get the unique label for the (v,j) pairs
    labels_ini = linear_2d_index(v_nnz, j_nnz, n_i=levels.v_max_allowed)
    labels_fin = linear_2d_index(vp_nnz, jp_nnz, n_i=levels.v_max_allowed)

    a_reduced = zeros((levels.size, levels.size), 'f8')

    for i, A_i in enumerate(a_nnz):

        # print('{:4}/{:4}'.format(i+1, len(A_nnz)))

        # get the indices based on v,j, jp, jp comaprisons
        #     v, j, vp, jp = v_nnz[i], j_nnz[i], vp_nnz[i], jp_nnz[i]
        #     ind_ini = where((levels['v'] == v)*(levels['j'] == j))[0]
        #     ind_fin = where((levels['v'] == vp)*(levels['j'] == jp))[0]

        # get the indices based on label comparisons
        ind_ini = where(levels.data['label'] == labels_ini[i])[0]
        ind_fin = where(levels.data['label'] == labels_fin[i])[0]

        if ind_ini.size != 0 or ind_fin.size != 0:
            a_reduced[ind_ini, ind_fin] = A_i
        else:
            continue

    # DEBUG
    # A_reduced[A_reduced > 0.0] = numpy.log10(A_reduced[A_reduced > 0.0])
    # pylab.imshow(A_reduced, interpolation='none')
    # pylab.colorbar()
    # pylab.show()

    return a_reduced


def reduce_einstein_coefficients(a_mat, energy_levels):
    """
    Recude a 4D A transition ndarray into a 2D matrix

    .. todo:: rename this function to
     reduce_einstein_coeffients_two_quantum_numbers or a name that indicates a
     molecule whose levels are identified by two "quantum" numbers or labels.

    Given the array A that is indexed using four indices A[v, j, v', j']
    returns an array A_reduced that is indexed with two indices A_reduced[i, f]
    where i and f are the initial and final levels respectively.

    :param ndarray a_mat: A 4D matrix holding the A coefficients.
     A[v, j, v', j'] where (v,j) are the initial levels and (v',j') are the
     final levels.
    :param EnergyLevelsSpeciesBase energy_levels: numpy array array of N rows,
     the v and j levels and their corresponding energy. The elements of the
     array are stored in increasing order in the energy. This array can be
     obtained e.g. from read_levels.read_levels_lique() .. todo:: check lique!!
    :return: A_reduced that is a square matrix of shape
     (energy_levels.size, energy_levels.size) with the nonzero values of A
     mapped to the indices of the levels in the array energy_levels.
    """
    levels = energy_levels
    n_levels = len(energy_levels.data)
    labels = energy_levels.data['label']

    # find the non zeros elements and their corresponding indices in A
    (v_nnz, j_nnz, vp_nnz, jp_nnz), a_nnz = where(a_mat > 0), a_mat[a_mat > 0]

    # get the unique label for the (v,j) pairs
    labels_ini = linear_2d_index(v_nnz, j_nnz, n_i=levels.v_max_allowed)
    labels_fin = linear_2d_index(vp_nnz, jp_nnz, n_i=levels.v_max_allowed)

    # keep transitions whose initial levels labels and the final label of the
    # transition are found in energy_levels
    mask = in1d(labels_ini, labels)*in1d(labels_fin, labels)
    labels_ini, labels_fin = labels_ini[mask], labels_fin[mask]
    a_nnz = a_nnz[mask]

    # get the indices of the labels of the levels in the transitions (that are
    # now all a subset of the energy_levels)
    inds_ini = find_matching_indices(labels, labels_ini)
    inds_fin = find_matching_indices(labels, labels_fin)

    # define the reduced A matrix and fill it up using inds_ini and inds_fin
    a_reduced = zeros((n_levels, n_levels), 'f8') * a_nnz.unit

    a_reduced[inds_ini, inds_fin] = a_nnz

    # DEBUG
    # A_reduced[A_reduced > 0.0] = numpy.log10(A_reduced[A_reduced > 0.0])
    # pylab.imshow(A_reduced, interpolation='none')
    # pylab.colorbar()
    # pylab.show()

    return a_reduced


def compute_delta_energy_matrix(levels):
    """
    Compute the energy difference matrix

    Given the energy levels, returns the delta energy matrix \Delta E = E - E^T
    that is documented in the notebook doc/rate_equations_derivation.ipynb

    :param read_energy_levels.EnergyLevelsBase levels: The energy levels
     object or a subclass of it that has the energies defined in the attribute
     record levels.data['E'].
    :return: square matrix of shape n x n where n is the number of energy
     levels
    """
    n = len(levels.data)
    energies_as_column = levels.data['E'].reshape(1, n)

    # the energy matrix with identical columns
    E_matrix = numpy.repeat(energies_as_column, n, axis=0).T

    delta_E = E_matrix - E_matrix.T

    return delta_E


def compute_degeneracy_matrix(levels):
    """
    Compute the degeneracy matrix using an energy levels object as input

    Given the energy levels, returns the degeneracy matrix R that is strictly
    upper triangular that is documented in the notebook
    doc/rate_equations_derivation.ipynb

    :param EnergyLevelsSpeciesBase levels: The energy levels object
    :return: square strictly upper triangular matrix of shape n x n.
    """
    n = len(levels.data)
    degeneracies_as_column = levels.data['g'].reshape(1, n)

    G = numpy.array(numpy.repeat(degeneracies_as_column, n, axis=0).T)

    # a strict upper triangular matrix
    one_u_nn = numpy.tril(numpy.ones((n, n), 'f8'), -1).T

    # the R matrix (the ratios of the degeneracies)
    R = (G * (1.0 / G.T)).T * one_u_nn

    return R


def compute_b_j_nu_matrix_from_a_matrix(energy_levels,
                                        a_matrix,
                                        t_rad):
    """
    Compute the stimulated emission and absorption coefficients matrix

    These quantities are computed from the energy levels and the A_matrix.
    see the notebook doc/rate_equations_derivation.ipynb for more details

    https://en.wikipedia.org/wiki/Einstein_coefficients
    http://www.ifa.hawaii.edu/users/kud/teaching_12/3_Radiative_transfer.pdf

    .. todo:: replace the J_nu in the name of this function and in the body
    .. todo:: to something that represents energy density like u_nu

    :param EnergyLevelsBase energy_levels: The energy levels
     object or a subclass of it that has the energies defined in the attribute
     record levels.data['E'] (see the documentation of
     compute_delta_energy_matrix).
    :param astropy.units.quantity.Quantity a_matrix: The spontaneous emission
     coefficients matrix (A in the ipython notebook).
    :param Quantity t_rad: The radiation temperature.
    :return: The B matrix defined in the notebook multiplied by J_nu
    """
    delta_e = compute_delta_energy_matrix(energy_levels)

    nu_matrix = (fabs(delta_e) / h_planck).to(u.Hz)

    b_e_matrix = a_matrix / (8.0 * pi * h_planck * nu_matrix**3 / c_light**3)
    numpy.fill_diagonal(b_e_matrix, 0.0)

    r_matrix = compute_degeneracy_matrix(energy_levels)

    # B_nu is the planck function, when multiplied by 4pi/c we obtain the
    # spectral energy density usually called u_i and that has dimensions
    # of Energy / Length^3 / Hz
    j_nu_matrix = (4.0*pi*u.sr/c_light)*B_nu(nu_matrix, t_rad)
    numpy.fill_diagonal(j_nu_matrix, 0.0)

    b_a_matrix = b_e_matrix.T * r_matrix

    b_matrix = b_e_matrix + b_a_matrix

    b_j_nu_matrix = b_matrix * j_nu_matrix

    return b_j_nu_matrix


def reduce_collisional_coefficients_slow(
        cr_info_nnz,
        energy_levels,
        set_inelastic_coefficient_to_zero=False,
        set_excitation_coefficients_to_zero=False,
        reduced_data_is_upper_to_lower_only=True):
    """
    Construct the K_matrix from the sparese nnz collisional coefficeints data

    Given the data that is available for the collisional transitions as a
    function of initial and final v,j -> v',j' levels and for different values
    of temperature, returns an array of the data reduced as a matrix for each
    temperature value and mapped by the energy level labels instead of v and j.

    :param tuple cr_info_nnz: The information of the collisional data of the
     levels for which data is available. This parameter is usually returned
     by a read of raw collisional data e.g
     read_collision_coefficients_lique_and_wrathmall.

     The tuple has four elements:

         - (v,j): The first element is a 2D array of shape (2, n_transitions)
           that are the v and j of the initial level of the transitions (ini).
           The columns of this array (v, j = ini) are the initial v and j of
           the transitions for a certain T section. i.e. the number of non-zero
           elements in K for a certain temperature is equal to the number of
           elements in the v or j columns.
           v.size = j.size = where(K[..., 0] > 0)[0].size

         - (v',j') The second element is the same of the first element but for
          the final level of the transitions (v', j' = fin)

         - unique_nnz: The third element is a 2D array of shape
          (2, n_unique_transitions) that are the unique levels involved in all
          the transitions.

         - cr_nnz: The last element is an array of shape (T.size, n_transitions)
           which are the collisional coefficient rates with non-zero values
           for each value of temperature in the T array.

    :param EnergyLevelsMolecular energy_levels: The energy levels object whose
     evergy levels map to the collisional nnz data.
    :param bool set_inelastic_coefficient_to_zero: If True, then all the
     elements along the diagonal of the matrix for all the temperatures are set
     to zero. i.e the collision coefficient for the in-elastic collisions
     would be ignored.
    :param bool set_excitation_coefficients_to_zero: If True, then all the
     elements in the upper triangular part of the matrix for all the
     temperatures are set to zero. i.e the lower to upper transitions
     (excitation) are ignored. 
    :param bool reduced_data_is_upper_to_lower_only: If True, then it is assumed
     that the computed reduced matrix has only upper to lower
     (i.e de-excitation) coefficients. If there is any non-zero value along the
     diagonal (i.e in-elastic coefficients) or any non zero value in the upper
     triangular part (i.e excitation coeffients), then an exception is raised.
    :return: The reduced matrices of the collisional coefficients. One matrix
     for each temperature value. The shape of the matrix is
      (n_levels, n_levels, n_temperature_values)
    """
    # check_self_transitions_in_einstien_nnz_data(A_info_nnz)

    levels = energy_levels
    n_levels = len(energy_levels.data)
    labels = energy_levels.data['label']

    (v_nnz, j_nnz), (vp_nnz, jp_nnz), unique_nnz, cr_nnz = cr_info_nnz

    # get the unique label for the (v,j) pairs
    labels_ini = linear_2d_index(v_nnz, j_nnz, n_i=levels.v_max_allowed)
    labels_fin = linear_2d_index(vp_nnz, jp_nnz, n_i=levels.v_max_allowed)

    # number of temperature value for which collisional data is available
    n_T = cr_nnz.shape[0]

    k_dex_reduced = zeros((n_levels, n_levels, n_T), 'f8') * cr_nnz.unit

    for i, cr_i in enumerate(cr_nnz.T):

        # print('{:4}/{:4}'.format(i+1, len(A_nnz)))

        # get the indices based on v,j, jp, jp comparisons
        #     v, j, vp, jp = v_nnz[i], j_nnz[i], vp_nnz[i], jp_nnz[i]
        #     ind_ini = where((levels['v'] == v)*(levels['j'] == j))[0]
        #     ind_fin = where((levels['v'] == vp)*(levels['j'] == jp))[0]

        # get the indices based on label comparisons
        ind_ini = where(labels == labels_ini[i])[0]
        ind_fin = where(labels == labels_fin[i])[0]

        if ind_ini.size != 0 or ind_fin.size != 0:
            k_dex_reduced[ind_ini, ind_fin, :] = cr_i
        else:
            continue

    #
    # optionally zero out data above the diagonals
    #
    if set_inelastic_coefficient_to_zero:
        i_diag, j_diag = numpy.diag_indices(n_levels)
        k_dex_reduced[i_diag, j_diag, :] = 0.0
    if set_excitation_coefficients_to_zero:
        i_upper, j_upper = numpy.triu_indices(n_levels, 1)
        k_dex_reduced[i_upper, j_upper, :] = 0.0

    if reduced_data_is_upper_to_lower_only:
        # check that the upper triangular matrices for all the temperatures
        # including the diagonal are zero since this is the K_dex matrix
        # (see doc)
        assert numpy.triu(numpy.moveaxis(k_dex_reduced, -1, 0)).sum() == 0.0

    # DEBUG
    # K_dex_reduced[K_dex_reduced > 0.0] = numpy.log10(
    #     K_dex_reduced[K_dex_reduced > 0.0])
    # pylab.imshow(K_dex_reduced[:, :, 0], interpolation='none')
    # pylab.colorbar()
    # pylab.show()

    return k_dex_reduced


def reduce_collisional_coefficients(cr, energy_levels):
    """.. todo:: add doc

    :param cr: .. todo:: add doc
    :param energy_levels: .. todo:: add doc
    :return: .. todo:: add doc
    """
    raise NotImplementedError("not implemented yet")


def compute_k_dex_matrix_interpolator(k_dex_vs_tkin, t_range):
    """

    :param ndarray k_dex_vs_tkin: The K dexcitation matrix as a function of
     temperature. The last dimension k_dex_vs_tkin.shape[-1] and t_range.size
     should have the same value
    :param ndarray t_range: The values of the temperatures corresponding to the
     last dimension of k_dex_vs_tkin
    :return: callable: The interpolation function that returns a square matrix
     of the collisional coefficient given a temperature.
    """
    # get the linear interpolator of the upper to lower collision rates as a
    # function of temperature (the last axis). This function returns an array
    # that is the same shape of K_dex[..., 0]
    k_dex_interpolator = scipy.interpolate.interp1d(t_range, k_dex_vs_tkin)

    return lambda t_kin: k_dex_interpolator(t_kin) * k_dex_vs_tkin.unit


def compute_k_matrix_from_k_dex_matrix(energy_levels,
                                       k_dex_matrix_interpolator,
                                       t_kin):
    """
    Compute the K matrix from K de-excitation matrix

    The K_dex matrix is a lower traignular matrix, it is used to populte the
    the upper triangular matrix since they are related through the R matrix
    (see rate_quations_derivation.ipynb notebook for details)

    :param EnergyLevelsSpeciesBase energy_levels: The energy levels
    :param callable k_dex_matrix_interpolator:  a callable function that takes
     a kinetic temperature as an argument and returns k_dex matrix at that
     specfied temperature (through e.g. interpolation).
    :param float t_kin: The kinetic temperature
    :return: ndarray: The K matrix
    """
    delta_e_matrix = fabs(compute_delta_energy_matrix(energy_levels))

    r_matrix = compute_degeneracy_matrix(energy_levels)

    # R*K_{dex}^T(T) to be multiplied by the exp(-dE/kb*T) in the loop
    k_dex_t = k_dex_matrix_interpolator(t_kin)
    k_ex_t = r_matrix * k_dex_t.T * exp(-delta_e_matrix / (kb * t_kin))

    k_matrix = k_dex_t + k_ex_t

    return k_matrix


def compute_transition_rate_matrix(data_set,
                                   t_kin,
                                   t_rad,
                                   collider_density):
    """
    compute the matrix M that can be used to compute the dn/dt

    The rates dn/dt can be computed from M by multiplying it by the abundances
    :math:`dn/dt = M.n`

    :param frigus.readers.dataset.DataSetBase: The data of the species
    :param Quantity t_kin: The kinetic temperature at which the steady state
     computation  will be done.
    :param Quantity t_rad: The radiation temperature at which the steady state
     computation will be done.
    :param Quantity collider_density: The density of the collider species.
    :return: Quantity ndarray: The M matrix as a nxn ndarray
    """
    energy_levels = data_set.energy_levels
    a_matrix = data_set.a_matrix
    k_dex_matrix_interpolator = data_set.k_dex_matrix_interpolator

    # compute the stimulated emission and absorption coefficients matrix
    b_jnu_matrix = compute_b_j_nu_matrix_from_a_matrix(
        energy_levels,
        a_matrix,
        t_rad
    )

    # get the K matrix for a certain temperature in the tabulated range
    k_matrix = compute_k_matrix_from_k_dex_matrix(
        energy_levels,
        k_dex_matrix_interpolator,
        t_kin
    )

    # compute the M matrix that can be used to compute the equilibrium state of
    # the levels (see notebook)
    o_matrix = (a_matrix + b_jnu_matrix + k_matrix * collider_density).T

    d_matrix = -numpy.eye(o_matrix.shape[0]) * o_matrix.sum(axis=0)

    m_matrix = o_matrix + d_matrix

    return m_matrix


def population_density_at_steady_state(data_set,
                                       t_kin=None,
                                       t_rad=None,
                                       collider_density=None):
    """
    Compute the population density at steady state by solving the linear system

    compute the matrix M that can be used to compute the dn/dt

    The rates dn/dt can be computed from M by multiplying it by the abundances
    :math:`dn/dt = M.n`

    :param frigus.readers.dataset.DataSetBase: The data set of the species
    :param Quantity t_kin: The kinetic temperature at which the steady state
     computation  will be done.
    :param Quantity t_rad: The radiation temperature at which the steady state
     computation will be done.
    :param Quantity collider_density: The density of the collider species.
    :return: ndarray: The equilibrium population density as a column vector
    """

    m_matrix = compute_transition_rate_matrix(
        data_set,
        t_kin,
        t_rad,
        collider_density
    )

    x_equilibrium = solve_equilibrium(m_matrix.si.value)

    # assert bool((x_equilibrium < 0.0).any()) is False
    # assert numpy.fabs(1.0 - numpy.fabs(x_equilibrium.sum())) <= 1e-3

    return x_equilibrium


def cooling_rate_at_steady_state(data_set, t_kin, t_rad, collider_density):
    """
    Compute the cooling rate at steady state

    The cooling rate is computed for a certain:

       - kinetic temparature
       - radiation temperature
       - collider density

    :param DatasetBase data_set: The species dataset
    :param astropy.units.quantity.Quantity t_kin: the kinetic temperature
    :param astropy.units.quantity.Quantity t_rad: the radiation temperature
    :param astropy.units.quantity.Quantity collider_density: The density of the
     collider species
    :return: the cooling rate
    """

    x_equilibrium = population_density_at_steady_state(
        data_set,
        t_kin,
        t_rad,
        collider_density
    )

    # compute the cooling rate (per particle)
    return cooling_rate(
        x_equilibrium,
        data_set.energy_levels,
        data_set.a_matrix
    )


def cooling_rate(population_densities, energy_levels, a_matrix):
    """
    Compute the cooling rate due to the spontaneous transitions.

    :param array_like population_densities: A column vector of the population
     densities. This is a dimensionless vector of shape nx1, where n is the 
     number of energy levels.
    :param read_energy_levels.EnergyLevelsBase energy_levels: The energy levels
     object or a subclass of it that has the energies defined in the attribute
     record levels.data['E'].
    :param array_lik a_matrix: A square matrix that has a shape n x n where n
     is the number of energy levels in "energy_levels". The elements of the
     matrix are the spontaneous transition rates. An element A[upper, lower] 
     should be read (interpreted) as: The spontaneous transition probability
     per unit time from the level "upper" to the level "lower". A_matrix is
     assumed to be a strictly lower triangular matrix (this is not checked,
     thus it is the responsibility of the called to assure that).
    :return: scalar astropy.units.quantity.Quantity: The cooling rate due to
     all the transitions in units of A_matrix.unit * energy_levels['E'].units. 
    """
    energy_levels_unit = energy_levels.data['E'].unit
    a_matrix_unit = a_matrix.unit

    delta_e_matrix = fabs(compute_delta_energy_matrix(energy_levels)).value
    a_matrix = a_matrix.value

    retval = (a_matrix * delta_e_matrix * population_densities).sum()

    return retval * energy_levels_unit * a_matrix_unit

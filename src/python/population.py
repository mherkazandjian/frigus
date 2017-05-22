from __future__ import print_function

import os

import pylab

import numpy
from numpy import (zeros, fabs, arange, array_equal, exp, ones, log10,
                   linalg, eye, dot, where, intersect1d, setdiff1d, in1d, pi)
import scipy
from scipy import interpolate

from astropy import units as u
from astropy.constants import c as c_light
from astropy.constants import h as h_planck
from astropy.constants import k_B as kb
from astropy.analytic_functions import blackbody_nu as B_nu

from utils import linear_2d_index, find_matching_indices


dropbox_path = '~/dropbox/Dropbox/us/cooling_function/mher-1e6-100'


def find_v_max_j_max_from_data(A_einstein_nnz, cr_coefficients_nnz):
    """
    :param A_einstein_nnz: .. todo:: add doc
    :param cr_coefficients_nnz:  .. todo:: add doc
    :return:  .. todo:: add doc
    """

    # find the non zeros elements and their corresponding indices of the
    # Einstein coefficients
    v_nnz, j_nnz, vp_nnz, jp_nnz, A_nnz = A_einstein_nnz
    v_max_A = max(v_nnz.max(), vp_nnz.max())
    j_max_A = max(j_nnz.max(), jp_nnz.max())
    print('DEBUG: (v_max_A, j_max_A) = ', v_max_A, j_max_A)

    # find the non zeros elements and their corresponding indices of the
    # collisional coefficients
    (v_nnz, j_nnz), (vp_nnz, jp_nnz), unique_nnz, cr_nnz = cr_coefficients_nnz
    v_max_cr = max(v_nnz.max(), vp_nnz.max())
    j_max_cr = max(j_nnz.max(), jp_nnz.max())
    print('DEBUG: (v_max_cr, j_max_cr) = ', v_max_cr, j_max_cr)

    return max(v_max_A, v_max_cr), max(j_max_A, j_max_cr)


def check_self_transitions_in_Einstien_nnz_data(A_info_nnz):
    """raises an error if there are self transitions in the non-zero data of
    the Einstein coefficients.

    :param A_info_nnz: a tuple of elements v_nnz, j_nnz, vp_nnz, jp_nnz, A_nnz
    :return: True if everything ok, else error is raised
    """
    v_nnz, j_nnz, vp_nnz, jp_nnz, A_nnz = A_info_nnz

    for i, A in enumerate(A_nnz):
        v, j, vp, jp = v_nnz[i], j_nnz[i], vp_nnz[i], jp_nnz[i]
        if v == vp and j == jp:
            raise ValueError('v -> v, j -> j transition. This is'
                             'not possible')

    return True


def reduce_einstein_coefficients_slow(A_info_nnz, energy_levels):
    """Use the nonzero entries of the Einstein coefficients to construct the
    A matrix for the levels found in the energy_levels parameter.

    :param A_info_nnz: a tuple of elements v_nnz, j_nnz, vp_nnz, jp_nnz, A_nnz
    :param energy_levels: numpy array array of N rows, the v and j levels and
     their corresponding energy. The elements of the array are stored in
     increasing order in the energy. This array can be obtained e.g. from
     read_levels.read_levels_lique()
    :return: A 2D matrix of shape (energy_levels.size, energy_levels.size) that
     is a lower triangular matrix containing the Einstein coefficients.
    """

    check_self_transitions_in_Einstien_nnz_data(A_info_nnz)

    v_nnz, j_nnz, vp_nnz, jp_nnz, A_nnz = A_info_nnz

    levels = energy_levels

    # get the unique label for the (v,j) pairs
    labels_ini = linear_2d_index(v_nnz, j_nnz, n_i=levels.v_max_allowed)
    labels_fin = linear_2d_index(vp_nnz, jp_nnz, n_i=levels.v_max_allowed)

    A_reduced = zeros((levels.size, levels.size), 'f8')

    for i, A_i in enumerate(A_nnz):

        # print('{:4}/{:4}'.format(i+1, len(A_nnz)))

        # get the indices based on v,j, jp, jp comaprisons
        #     v, j, vp, jp = v_nnz[i], j_nnz[i], vp_nnz[i], jp_nnz[i]
        #     ind_ini = where((levels['v'] == v)*(levels['j'] == j))[0]
        #     ind_fin = where((levels['v'] == vp)*(levels['j'] == jp))[0]

        # get the indices based on label comparisons
        ind_ini = where(levels.data['label'] == labels_ini[i])[0]
        ind_fin = where(levels.data['label'] == labels_fin[i])[0]

        if ind_ini.size != 0 or ind_fin.size != 0:
            A_reduced[ind_ini, ind_fin] = A_i
        else:
            continue


    # DEBUG
    # A_reduced[A_reduced > 0.0] = numpy.log10(A_reduced[A_reduced > 0.0])
    # pylab.imshow(A_reduced, interpolation='none')
    # pylab.colorbar()
    # pylab.show()

    return A_reduced


def reduce_einstein_coefficients(A, energy_levels):
    """
    .. todo:: rename this function to 
     reduce_einstein_coeffients_two_quantum_numbers or a name that indicates a
     molecule whose levels are identified by two "quantum" numbers or labels.
    Given the array A that is indexed using four indices A[v, j, v', j']
    returns an array A_reduced that is indexed with two indices A_reduced[i, f]
    where i and f are the initial and final levels respectively.

    :param A: A 4D matrix holding the A coefficients. A[v, j, v', j'] where
     (v,j) are the initial levels and (v',j') are the final levels.
    :param energy_levels: numpy array array of N rows, the v and j levels and
     their corresponding energy. The elements of the array are stored in
     increasing order in the energy. This array can be obtained e.g. from
     read_levels.read_levels_lique()
    :return: A_reduced that is a square matrix of shape
     (energy_levels.size, energy_levels.size) with the nonzero values of A
     mapped to the indices of the levels in the array energy_levels.
    """

    levels = energy_levels
    n_levels = len(energy_levels.data)
    labels = energy_levels.data['label']

    # find the non zeros elements and their corresponding indices in A
    (v_nnz, j_nnz, vp_nnz, jp_nnz), A_nnz = where(A > 0.0), A[A > 0.0]

    # get the unique label for the (v,j) pairs
    labels_ini = linear_2d_index(v_nnz, j_nnz, n_i=levels.v_max_allowed)
    labels_fin = linear_2d_index(vp_nnz, jp_nnz, n_i=levels.v_max_allowed)

    # keep transitions whose initial levels labels and the final label of the
    # transition are found in energy_levels
    mask = in1d(labels_ini, labels)*in1d(labels_fin, labels)
    labels_ini, labels_fin = labels_ini[mask], labels_fin[mask]
    A_nnz = A_nnz[mask]

    # get the indices of the labels of the levels in the transitions (that are
    # now all a subset of the energy_levels)
    inds_ini = find_matching_indices(labels, labels_ini)
    inds_fin = find_matching_indices(labels, labels_fin)

    # define the reduced A matrix and fill it up using inds_ini and inds_fin
    A_reduced = zeros((n_levels, n_levels), 'f8') * A_nnz.unit

    A_reduced[inds_ini, inds_fin] = A_nnz

    # DEBUG
    # A_reduced[A_reduced > 0.0] = numpy.log10(A_reduced[A_reduced > 0.0])
    # pylab.imshow(A_reduced, interpolation='none')
    # pylab.colorbar()
    # pylab.show()

    return A_reduced


def compute_delta_energy_matrix(levels):
    """
    Compute the energy difference matrix (see notebook .. todo:: notebook ref)
    
    Given the energy levels, returns the delta energy matrix \Delta E = E - E^T
    that is documented in the notebook.
    :param read_energy_levels.EnergyLevelsBase levels: The energy levels
     object or a subclass of it that has the energies defined in the attribute
     record levels.data['E'].
    :return: square matrix of shape n x n where n is the number of energy
     levels (see notebook .. todo:: notebook ref).
    """
    n = len(levels.data)
    energies_as_column = levels.data['E'].reshape(1, n)

    # the energy matrix with identical columns
    E_matrix = numpy.repeat(energies_as_column, n, axis=0).T

    delta_E = E_matrix - E_matrix.T

    return delta_E


def compute_degeneracy_matrix(levels):
    """
    Given the energy levels, returns the degeneracy matrix R that is strictly
    upper triangular that is documented in the notebook.

    :param levels:  .. todo:: add doc
    :return: square strictly upper triangular matrix of shape n x n.
    """
    n = len(levels.data)
    degeneracies_as_column = levels.data['g'].reshape(1, n)

    G = numpy.array(numpy.repeat(degeneracies_as_column, n, axis=0).T)

    # a strict upper triangular matrix
    one_U_nn = numpy.tril(numpy.ones((n, n), 'f8'), -1).T

    # the R matrix (the ratios of the degeneracies)
    R = (G * (1.0 / G.T)).T * one_U_nn

    return R


def compute_B_J_nu_matrix_from_A_matrix(energy_levels,
                                        A_matrix,
                                        T_rad,
                                        debug=False):
    """
    Given the energy levels, returns the stimulated emission and absorption
    coefficients matrix.

    https://en.wikipedia.org/wiki/Einstein_coefficients
    http://www.ifa.hawaii.edu/users/kud/teaching_12/3_Radiative_transfer.pdf

    .. todo:: replace the J_nu in the name of this function and in the body
    .. todo:: to something that represents energy density like u_nu

    :param read_energy_levels.EnergyLevelsBase energy_levels: The energy levels
     object or a subclass of it that has the energies defined in the attribute
     record levels.data['E'] (see the documentation of
     compute_delta_energy_matrix).
    :param A_matrix: The spontaneous emission coefficients matrix (A in the
     ipython notebook).
    :param T_rad: The radiation temperature.
    :return: The B matrix defined in the notebook multiplied by J_nu
    """
    delta_e = compute_delta_energy_matrix(energy_levels)

    nu_matrix = (fabs(delta_e) / h_planck).to(u.Hz)

    B_e_matrix = A_matrix / (8.0*pi*h_planck*nu_matrix**3/c_light**3)
    numpy.fill_diagonal(B_e_matrix, 0.0)

    R_matrix = compute_degeneracy_matrix(energy_levels)

    # B_nu is the planck function, when multiplied by 4pi/c we obtain the
    # spectral energy density usually called u_i and that has dimensions
    # of Energy / Length^3 / Hz
    J_nu_matrix = (4.0*pi*u.sr/c_light)*B_nu(nu_matrix, T_rad)
    numpy.fill_diagonal(J_nu_matrix, 0.0)

    B_a_matrix = B_e_matrix.T * R_matrix

    B_matrix = B_e_matrix + B_a_matrix

    B_J_nu_matrix = B_matrix * J_nu_matrix

    if debug is True:
        numpy.savetxt(
            os.path.expanduser(
                os.path.join(
                    dropbox_path, 'A_matrix.txt'))
            ,
            A_matrix.si.value,
            fmt='%+-1.16e')

        numpy.savetxt(
            os.path.expanduser(
                os.path.join(
                    dropbox_path, 'B_matrix.txt'))
            ,
            B_matrix.si.value,
            fmt='%+-1.16e')

        numpy.savetxt(
            os.path.expanduser(
                os.path.join(
                    dropbox_path, 'J_nu_matrix.txt')),
            J_nu_matrix.to(u.N / u.m**2 * u.s).value,
            fmt='%+-1.16e')

        numpy.savetxt(
            os.path.expanduser(
                os.path.join(
                    dropbox_path, 'B_J_nu_matrix.txt')),
            B_J_nu_matrix.si.value,
            fmt='%+-1.16e')

    return B_J_nu_matrix


def reduce_collisional_coefficients_slow(
        cr_info_nnz,
        energy_levels,
        set_inelastic_coefficient_to_zero=False,
        set_excitation_coefficients_to_zero=False,
        reduced_data_is_upper_to_lower_only=True):
    """
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
    # check_self_transitions_in_Einstien_nnz_data(A_info_nnz)

    levels = energy_levels
    n_levels = len(energy_levels.data)
    labels = energy_levels.data['label']

    (v_nnz, j_nnz), (vp_nnz, jp_nnz), unique_nnz, cr_nnz = cr_info_nnz

    # get the unique label for the (v,j) pairs
    labels_ini = linear_2d_index(v_nnz, j_nnz, n_i=levels.v_max_allowed)
    labels_fin = linear_2d_index(vp_nnz, jp_nnz, n_i=levels.v_max_allowed)

    # number of temperature value for which collisional data is available
    n_T = cr_nnz.shape[0]

    K_dex_reduced = zeros((n_levels, n_levels, n_T), 'f8') * cr_nnz.unit

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
            K_dex_reduced[ind_ini, ind_fin, :] = cr_i
        else:
            continue

    #
    # optionally zero out data above the diagonals
    #
    if set_inelastic_coefficient_to_zero:
        i_diag, j_diag = numpy.diag_indices(n_levels)
        K_dex_reduced[i_diag, j_diag, :] = 0.0
    if set_excitation_coefficients_to_zero:
        i_upper, j_upper = numpy.triu_indices(n_levels, 1)
        K_dex_reduced[i_upper, j_upper, :] = 0.0

    if reduced_data_is_upper_to_lower_only:
        # check that the upper triangular matrices for all the temperatures
        # including the diagonal are zero since this is the K_dex matrix
        # (see doc)
        assert numpy.triu(numpy.moveaxis(K_dex_reduced, -1, 0)).sum() == 0.0

    # DEBUG
    # K_dex_reduced[K_dex_reduced > 0.0] = numpy.log10(
    #     K_dex_reduced[K_dex_reduced > 0.0])
    # pylab.imshow(K_dex_reduced[:, :, 0], interpolation='none')
    # pylab.colorbar()
    # pylab.show()

    return K_dex_reduced


def reduce_collisional_coefficients(cr, energy_levels):
    """.. todo:: add doc

    :param cr: .. todo:: add doc
    :param energy_levels: .. todo:: add doc
    :return: .. todo:: add doc
    """
    raise NotImplementedError("not implemented yet")


def compute_K_dex_matrix_interpolator(K_dex_vs_T, T_range):
    """

    :param K_dex_vs_T:
    :param T_range:
    :return:
    """
    # get the linear interpolator of the upper to lower collision rates as a
    # function of temperature (the last axis). This function returns an array
    # that is the same shape of K_dex[..., 0]
    K_dex_interpolator = scipy.interpolate.interp1d(T_range, K_dex_vs_T)

    return lambda T_kin: K_dex_interpolator(T_kin)*K_dex_vs_T.unit


def compute_K_matrix_from_K_dex_matrix(energy_levels,
                                       K_dex_matrix_interpolator,
                                       T_kin,
                                       debug=True):
    """ .. todo:: add doc

    :param energy_levels: .. todo:: add doc
    :param K_dex_matrix_interpolator:  .. todo:: add doc
    :param T_range: .. todo:: add doc
    :param T_kin: .. todo:: add doc
    :return: .. todo:: add doc
    """
    delta_e_matrix = fabs(compute_delta_energy_matrix(energy_levels))

    R_matrix = compute_degeneracy_matrix(energy_levels)

    # R*K_{dex}^T(T) to be multiplied by the exp(-dE/kb*T) in the loop
    K_dex_T = K_dex_matrix_interpolator(T_kin)
    K_ex_T = R_matrix * K_dex_T.T * exp(-delta_e_matrix/(kb * T_kin))

    K_matrix = K_dex_T + K_ex_T

    if debug is True:
        numpy.savetxt(
            os.path.expanduser(
                os.path.join(
                    dropbox_path, 'K_nc_matrix.txt')),
            K_matrix.si.value*1e14,
            fmt='%+-1.16e')

    return K_matrix


def solveEquilibrium(M_matrix, debug=False):
    """
    Solve for the equilibrium population densities given the right hand side of
    the linear system of the rate equations dx/dt as a matrix dx/dt = A.x
    where M_matrix = A in the case of this function. The first row of the A
    is replaces by the conservation equation.

    :param matrix_like M_matrix: The right hand side matrix of the rate
    equation as an n x n matrix.
    :return: The population densities as a column vector. 
    """

    sz = M_matrix.shape[0]

    # solving directly. replacing the first row with the conservation equation
    # i.e the sum of the independent variable is 1, i.e the sum of the
    # population levels is
    dxdt = zeros((sz, 1), 'f8')
    M_matrix[0, :], dxdt[0] = 1.0, 1.0

    # solving the system A.x = b
    # before solving, we will divide each row by the diagonal
    A, b = M_matrix, dxdt

    # scale the rows by normalizing w.r.t the diagonal element
    for i in arange(sz):
        A[i, :] = A[i, :] / A[i, i]

    numpy.savetxt(
        os.path.expanduser(
            os.path.join(
                dropbox_path, 'solver_matrix.txt')),
                A,
            fmt='%+-1.16e')

    x = linalg.solve(A, b)

    return x


#def Tg2Tr(Tg):
#    '''Given the array of kinetic temperatures, it returns the
#     array with the corresponding radiation temperatures (through
#     the redshift parameter z)
#     '''
#     return Tr


def population_density_at_steady_state(data_set,
                                       T_kin,
                                       T_rad,
                                       collider_density,
                                       debug=False):
    """
    Compute the population density at steady state by solving the linear
     system. 

    :param readers.DataSet data_set: The data of the species
    :param T_kin: The kinetic temperature at which the steady state computation
     will be done.
    :param T_rad: The radiation temperature at which the steady state computation
     will be done.
    :param collider_density: The density of the collider species.
    :return: The equilibrium population density
    """

    energy_levels = data_set.energy_levels
    A_matrix = data_set.A_matrix
    K_dex_matrix_interpolator = data_set.K_dex_matrix_interpolator

    # compute the stimulated emission and absorption coefficients matrix
    B_J_nu_matrix = compute_B_J_nu_matrix_from_A_matrix(
        energy_levels,
        A_matrix,
        T_rad,
        debug=debug)

    # get the K matrix for a certain temperature in the tabulated range
    K_matrix = compute_K_matrix_from_K_dex_matrix(
        energy_levels,
        K_dex_matrix_interpolator,
        T_kin,
        debug=debug)

    # compute the M matrix that can be used to compute the equilibrium state of
    # the levels (see notebook)
    O_matrix = (A_matrix + B_J_nu_matrix + K_matrix * collider_density).T
    if debug is True:
        numpy.savetxt(
            os.path.expanduser(
                os.path.join(
                    dropbox_path, 'O_matrix.txt')),
            O_matrix.si.value,
            fmt='%+-1.16e')

    D_matrix = -numpy.eye(O_matrix.shape[0]) * O_matrix.sum(axis=0)
    if debug is True:
        numpy.savetxt(
            os.path.expanduser(
                os.path.join(
                    dropbox_path, 'D_matrix.txt')),
            D_matrix.si.value,
            fmt='%+-1.16e')

    M_matrix = O_matrix + D_matrix
    if debug is True:
        numpy.savetxt(
            os.path.expanduser(
                os.path.join(
                    dropbox_path, 'M_matrix.txt')),
            M_matrix.si.value,
            fmt='%+-1.16e')

    # solve the equilibrium population densities
    x_equilibrium = solveEquilibrium(M_matrix.si.value, debug=debug)

    assert numpy.fabs(1.0 - numpy.fabs(x_equilibrium.sum())) <= 1e-3

    return x_equilibrium


def cooling_rate_at_steady_state(data_set, T_kin, T_rad, collider_density):
    """.. todo:: add doc

    :param data_set: .. todo:: add doc
    :param T_kin: .. todo:: add doc
    :param T_Rad: .. todo:: add doc
    :param collider_density: .. todo:: add doc
    :return: the cooling rate
    """

    x_equilibrium = population_density_at_steady_state(data_set,
                                                       T_kin,
                                                       T_rad,
                                                       collider_density)

    # compute the cooling rate (per particle)
    c_rate = cooling_rate(x_equilibrium,
                          data_set.energy_levels,
                          data_set.A_matrix)

    return c_rate


def cooling_rate(population_densities, energy_levels, A_matrix):
    """compute the cooling rate due to the spontaneous transitions

    :param array_like population_densities: A column vector of the population
     densities. This is a dimensionless vector of shape nx1, where n is the 
     number of energy levels.
    :param read_energy_levels.EnergyLevelsBase energy_levels: The energy levels
     object or a subclass of it that has the energies defined in the attribute
     record levels.data['E'].
    :param A_matrix: .. todo:: add doc
    :return: .. todo:: add doc
    """
    delta_e_matrix = fabs(compute_delta_energy_matrix(energy_levels)).si.value
    A_matrix = A_matrix.si.value

    retval = (A_matrix * delta_e_matrix * population_densities).sum()

    retval = retval * u.Joule * u.second**-1

    return retval


def fit_glover(T):
    """
    fit of the cooling rate of H2 as a function of temperature (in K) in units
    of erg x cm^3 x s^-1
    .. todo:: add ref

    :param T: .. todo:: add doc
    :return: .. todo:: add doc
    """
    if 100.0 <= T and T <= 1000:
      retval = 10**(-24.311209
               +3.5692468*log10(T/1000.)
               -11.332860*(log10(T/1000.))**2
               -27.850082*(log10(T/1000.))**3
               -21.328264*(log10(T/1000.))**4
               -4.2519023*(log10(T/1000.))**5)
    elif 1000 < T and T <=6000:
      retval = 10**(-24.311209
               +4.6450521*log10(T/1000.)
               -3.7209846*log10((T/1000.))**2
               +5.9369081*log10((T/1000.))**3
               -5.5108047*log10((T/1000.))**4
               +1.5538288*log10((T/1000.))**5)
    else:
        raise ValueError("""out of bound""")

    return retval * u.erg * u.s**-1 * u.cm**3


def fit_lipovka(T, n_hd):
    """
    fit of the cooling rate of HD as a function of temperature (in K) in units
    of erg / s.

       HD revised cooling function - Lipovka, Nunez, Avila Reese 2005 

    :param u.quantity.Quantity T: The temperature range as an astropy quantity
    :param u.quantity.Quantity n_hd: The density of HD
    :return: u.quantity.Quantity: The cooling function in erg / s for the
     inpute T values

    .. see-also:: test_plot_fit_lipovka.py
    """
    #
    # make sure the intput temperatures and density value are within the
    # validity range of the fit and are of the correct types
    #
    assert isinstance(T, u.quantity.Quantity)
    assert isinstance(n_hd, u.quantity.Quantity)

    assert T.min() >= 100.0 * u.K
    assert T.max() <= 2000.0 * u.K
    assert n_hd.min() >= 1e6 * u.m**-3
    assert n_hd.max() <= 1e14 * u.m**-3

    lt_kin = numpy.log10(T.value)
    ln_hd = numpy.log10(n_hd.cgs.value)

    retval = 10.**(- 42.57688 * lt_kin ** 0 * ln_hd ** 0
                   + 0.92433 * lt_kin ** 0 * ln_hd ** 1
                   + 0.54962 * lt_kin ** 0 * ln_hd ** 2
                   - 0.07676 * lt_kin ** 0 * ln_hd ** 3
                   + 0.00275 * lt_kin ** 0 * ln_hd ** 4

                   + 21.93385 * lt_kin ** 1 * ln_hd ** 0
                   + 0.77952 * lt_kin ** 1 * ln_hd ** 1
                   - 1.06447 * lt_kin ** 1 * ln_hd ** 2
                   + 0.11864 * lt_kin ** 1 * ln_hd ** 3
                   - 0.00366 * lt_kin ** 1 * ln_hd ** 4

                   - 10.19097 * lt_kin ** 2 * ln_hd ** 0
                   - 0.54263 * lt_kin ** 2 * ln_hd ** 1
                   + 0.62343 * lt_kin ** 2 * ln_hd ** 2
                   - 0.07366 * lt_kin ** 2 * ln_hd ** 3
                   + 0.002514 * lt_kin ** 2 * ln_hd ** 4

                   + 2.19906 * lt_kin ** 3 * ln_hd ** 0
                   + 0.11711 * lt_kin ** 3 * ln_hd ** 1
                   - 0.13768 * lt_kin ** 3 * ln_hd ** 2
                   + 0.01759 * lt_kin ** 3 * ln_hd ** 3
                   - 0.000666317 * lt_kin ** 3 * ln_hd ** 4

                   - 0.17334 * lt_kin ** 4 * ln_hd ** 0
                   - 0.00835 * lt_kin ** 4 * ln_hd ** 1
                   + 0.0106 * lt_kin ** 4 * ln_hd ** 2
                   - 0.001482 * lt_kin ** 4 * ln_hd ** 3
                   + 0.000061926 * lt_kin ** 4 * ln_hd ** 4)

    return retval * u.erg * u.s**-1


def fit_lipovka_low_density(T):
    """
    fit of the cooling rate of HD as a function of temperature (in K) in units
    of erg x s^-1 for low temperature

    .. todo:: add ref

    :param T: .. todo:: add doc
    :return: .. todo:: add doc
    """

    lt_kin = numpy.log10(T.value)

    if numpy.log10(100.0) <= lt_kin and lt_kin <= numpy.log10(20000.0):
        retval = 10.0**(- 42.45906
                        + 21.90083 * lt_kin
                        - 10.1954 * lt_kin**2
                        + 2.19788 * lt_kin**3
                        - 0.17286 * lt_kin**4)
    else:
        msg = """temperature value {} is out of bound""".format(T)
        raise ValueError(msg)

    return retval * u.erg * u.s**-1

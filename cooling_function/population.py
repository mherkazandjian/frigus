from __future__ import print_function
import pylab

import numpy
from numpy import (zeros, fabs, arange, array_equal, exp, ones, log10,
                   linalg, eye, dot, where, intersect1d, setdiff1d, in1d, pi)
import scipy

import pdb

from IPython.core.debugger import Tracer

import read_einstien_coefficient
import read_collision_coefficients
import read_levels
from scipy.constants import c as c_light
from scipy.constants import h as h_planck
from scipy.constants import electron_volt as ev
from scipy.constants import Boltzmann as kb
from scipy.constants import hbar as hbar_planck

from Leiden_ISM.ismUtils import planck_function as J_nu

from utils import linear_2d_index, find_matching_indices

kb_ev = kb / ev


def find_v_max_j_max_from_data(A_einstein_nnz, cr_coefficients_nnz):
    """"""

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
     is a lower triangular matrix containing the Einstien coefficients.
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
    """Given the array A that is indexed using four indices A[v, j, v', j']
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
    n_levels = energy_levels.data.size
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
    A_reduced = zeros((n_levels, n_levels), 'f8')

    A_reduced[inds_ini, inds_fin] = A_nnz

    # DEBUG
    # A_reduced[A_reduced > 0.0] = numpy.log10(A_reduced[A_reduced > 0.0])
    # pylab.imshow(A_reduced, interpolation='none')
    # pylab.colorbar()
    # pylab.show()

    return A_reduced


def compute_delta_energy_matrix(levels):
    """Given the energy levels, returns the delta energy matrix
     \Delta E = E - E^T that is documented in the notebook.
    :return: square matrix of shape n x n where n is the number of energy
     levels.
    """
    n = levels.data.size
    energies_as_column = levels.data['E'].reshape(1, n)

    # the energy matrix with identical columns
    E_matrix = numpy.repeat(energies_as_column, n, axis=0).T

    delta_E = E_matrix - E_matrix.T

    return delta_E


def compute_degeneracy_matrix(levels):
    """Given the energy levels, returns the degeneracy matrix R that is
     strictly upper triangular that is documented in the notebook.
    :return: square matrix of shape n x n.
    """
    n = levels.data.size
    degeneracies_as_column = levels.data['g'].reshape(1, n)

    G = numpy.repeat(degeneracies_as_column, n, axis=0).T

    # a strict upper triangular matrix
    one_U_nn = numpy.tril(numpy.ones((n, n), 'f8'), -1).T

    # the R matrix (the ratios of the degeneracies)
    R = (G * (1.0 / G.T)).T * one_U_nn

    return R


def compute_B_J_nu_matrix_from_A_matrix(energy_levels, A_matrix, T):
    """given the energy levels, returns the stimulated emission and absorption
    coefficients matrix.
    https://en.wikipedia.org/wiki/Einstein_coefficients

    :param A_matrix: The spontaneous emission coefficients matrix (A in the
     ipython notebook)
    :param energy_levels: The energy levels
    :return: The B matrix defined in the notebook
    """
    delta_e = compute_delta_energy_matrix(energy_levels)

    nu_matrix = fabs(delta_e*ev) / h_planck

    B_e_matrix = (8.0*pi*hbar_planck/c_light**3)*(nu_matrix**3)*A_matrix

    R_matrix = compute_degeneracy_matrix(energy_levels)

    J_nu_matrix = J_nu(h_planck, nu_matrix, kb, T, c_light)
    numpy.fill_diagonal(J_nu_matrix, 0.0)

    B_a_matrix = B_e_matrix.T * R_matrix

    B_matrix = B_e_matrix + B_a_matrix

    B_J_nu_matrix = B_matrix * J_nu_matrix

    return B_J_nu_matrix


def reduce_collisional_coefficients_slow(cr_info_nnz, energy_levels):
    """

    :return:
    """
    # check_self_transitions_in_Einstien_nnz_data(A_info_nnz)

    levels = energy_levels
    n_levels = energy_levels.data.size
    labels = energy_levels.data['label']

    (v_nnz, j_nnz), (vp_nnz, jp_nnz), unique_nnz, cr_nnz = cr_info_nnz

    n_T = cr_nnz.shape[0]

    # get the unique label for the (v,j) pairs
    labels_ini = linear_2d_index(v_nnz, j_nnz, n_i=levels.v_max_allowed)
    labels_fin = linear_2d_index(vp_nnz, jp_nnz, n_i=levels.v_max_allowed)

    K_ex_reduced = zeros((levels.data.size, levels.data.size, n_T), 'f8')

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
            K_ex_reduced[ind_ini, ind_fin, :] = cr_i
        else:
            continue


    # DEBUG
    # K_ex_reduced[K_ex_reduced > 0.0] = numpy.log10(
    #     K_ex_reduced[K_ex_reduced > 0.0])
    # pylab.imshow(K_ex_reduced[:, :, 0], interpolation='none')
    # pylab.colorbar()
    # pylab.show()

    return K_ex_reduced


def reduce_collisional_coefficients(cr, energy_levels):
    """

    :return:
    """
    raise NotImplementedError("not implemented yet")


def compute_K_matrix_from_K_dex_matrix(energy_levels, K_dex, T_range, T):
    """

    :return:
    """
    delta_e_matrix = fabs(compute_delta_energy_matrix(energy_levels))

    R_matrix = compute_degeneracy_matrix(energy_levels)

    # get the linear interpolator of the upper to lower collision rates as a
    # function of temperature (the last axis). This function returns an array
    # that is the same shape of K_dex[..., 0]
    K_dex_interpolator = scipy.interpolate.interp1d(T_range, K_dex)

    # R*K_{dex}^T(T) to be multiplied by the exp(-dE/kb*T) in the loop
    K_dex_T = K_dex_interpolator(T)
    K_ex_T = R_matrix * K_dex_T.T * exp(-delta_e_matrix/(kb_ev * T))

    K_matrix = K_dex_T + K_ex_T

    return K_matrix


def solveEquilibrium(M_matrix):
    """solve for the equilibrium population densities.

    i.e solving A.x = b
    where M_matrix = A
    """

    sz = M_matrix.shape[0]
    
    # solving directly
    # replacing the first row with the conservation equation
    dndt = zeros((sz, 1), 'f8')
    M_matrix[0, :], dndt[0] = 1.0, 1.0

    # solving the system A.x = b
    # before solving, we will divide each row by the diagonal
    A, b = M_matrix, dndt
    for i in arange(sz):
        A[i, :] = A[i, :] / A[i, i]

    x = linalg.solve(A, b)

    # the fractional population density
    f = x
    return f


#def Tg2Tr(Tg):
#    '''Given the array of kinetic temperatures, it returns the
#     array with the corresponding radiation temperatures (through
#     the redshift parameter z)
#     '''
#     return Tr


def cooling_rate_at_steady_state(A_matrix,
                                 energy_levels,
                                 K_dex_matrix,
                                 T_rng,
                                 T_kin,
                                 collider_density):
    """

    :return:
    """

    # compute the stimulated emission and absorption coefficients matrix
    B_J_nu_matrix = compute_B_J_nu_matrix_from_A_matrix(energy_levels,
                                                        A_matrix,
                                                        T_kin)

    # get the K matrix for a certain temperature in the tabulated range
    K_matrix = compute_K_matrix_from_K_dex_matrix(energy_levels,
                                                  K_dex_matrix,
                                                  T_rng,
                                                  T_kin)
    assert (numpy.fabs(1.0 - K_matrix.sum() / 1.8873371663e-08) < 1e-10,
            "asdadasd")

    # compute the M matrix that can be used to compute the equilibrium state of
    # the levels (see notebook)
    O_matrix = (A_matrix + B_J_nu_matrix + K_matrix * collider_density).T

    D_matrix = numpy.zeros(O_matrix.shape, 'f8')
    D_matrix[numpy.diag_indices(D_matrix.shape[0])] = -O_matrix.sum(axis=0)

    M_matrix = O_matrix + D_matrix

    # solve the equilibrium population densities
    n_vj_equilibrium = solveEquilibrium(M_matrix)

    # compute the cooling rate
    c_rate = cooling_rate(n_vj_equilibrium, energy_levels, A_matrix)

    return c_rate


def cooling_rate(population_densities, energy_levels, A_matrix):
    """compute the cooling rate due to the spontaneous transitions

    :param population_densities:
    :param energy_levels:
    :param A_matrix:
    :return:
    """
    delta_e_matrix = fabs(compute_delta_energy_matrix(energy_levels))

    retval = dot(dot(A_matrix, delta_e_matrix), population_densities).sum()

    return retval


def fit_glover(T):
    """

    :param T:
    :return:
    """
    if 100 < T and T <= 1000:
      return   10**(-24.311209
               +3.5692468*log10(T/1000.)
               -11.332860*(log10(T/1000.))**2
               -27.850082*(log10(T/1000.))**3
               -21.328264*(log10(T/1000.))**4
               -4.2519023*(log10(T/1000.))**5)
    elif 1000 < T and T <=6000:
      return 10**(-24.311209
               +4.6450521*log10(T/1000.)
               -3.7209846*log10((T/1000.))**2
               +5.9369081*log10((T/1000.))**3
               -5.5108047*log10((T/1000.))**4
               +1.5538288*log10((T/1000.))**5)
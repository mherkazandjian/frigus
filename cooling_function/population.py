from __future__ import print_function
import pylab

import numpy
from numpy import (zeros, fabs, arange, array_equal, exp, ones, log10,
                   linalg, eye, dot, where, intersect1d, setdiff1d, in1d, pi)
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

    K_ex_reduced = zeros((n_T, levels.data.size, levels.data.size), 'f8')

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
            K_ex_reduced[:, ind_ini, ind_fin] = cr_i.reshape(cr_i.size, 1)
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

    # find the temperature interval within which T is
    T_ind = where(T == T_range)

    # R*K^T_{dex}(T) to be multiplied by the exp(-dE/kb*T) in the loop
    K_dex_T = K_dex[T_ind][0]
    K_ex_T = R_matrix * K_dex_T.T * exp(-delta_e_matrix/(kb_ev * T))

    K_matrix = K_dex_T + K_ex_T

    return K_matrix


def reduce_vj_repr(en, a_eins, cr, T,  ini, fin, vj_unique,
                   debug=False):
    """.. todo:: add doc"""

    # number of v and j levels respectively
    nv, nj = en.shape

    # linear indices for all the energy levels
    lind = arange(en.size)

    # the linear indices indexed as a 2D array
    # lind2d[v,j] = the linear index
    lind2d = lind.copy().reshape(en.shape)

    # convert the energies from a v,j representation
    # to a linear representation
    en_l = en.flatten()
    if debug is True:
        assert en[:,:] - en_l[lind2d[:,:]]

    # convert the representation of the einstein
    # coefficients from a 4D representation v,j,v',j'
    # to a 2D representation i,j
    a_eins_l = a_eins.reshape(en.size, en.size)
    # added for checking the reshaping
    if debug is True:
        for v in arange(nv):
            for j in arange(nj):
                for vp in arange(nv):
                    for jp in arange(nj):
                        a_l = a_eins_l[lind2d[v,j], lind2d[vp,jp]]
                        a = a_eins[v,j,vp,jp]
                        assert a_l == a

    # initializing the linearized version of the
    # the collisional reaction rates and 
    # initial and final quantum numbers; the  
    # reshape is done according to the total 
    # number of rovibrational levels of H2 
    cr_l = zeros((en.size, en.size, T.size), 'f8')
    ini_l = zeros(ini.shape[1], 'i')
    fin_l = zeros(fin.shape[1], 'i')
    
    for i, ((v,j), (vp,jp)) in enumerate(zip(ini.T, fin.T)):
        cr_l[lind2d[v,j], lind2d[vp,jp], :] = cr[v, j, vp, jp,:]
        ini_l[i] = lind2d[v,j]
        fin_l[i] = lind2d[vp,jp]
    vj_unique_l = lind2d[vj_unique[0], vj_unique[1]]
    # added for checking the reshaping
    # .. todo:: think of a better test
    # .. todo:: see how ini_l and fin_l can be checked
    if debug is True:
        for (v,j), (vp,jp) in zip(ini.T, fin.T):
            T1 = cr_l[lind2d[v,j], lind2d[vp,jp], :]
            T2 = cr[v, j, vp, jp, :]
            assert array_equal(T1, T2)
    

    # computing the degeneracies
    g = vj_unique[1]*2 + 1

    return en_l, a_eins_l, cr_l, ini_l, fin_l, vj_unique_l, g


#def computeRateMatrix(pNH3, Tkin, nc):
def computeRateMatrix(en, a_eins, cr, ini, fin, unique_col,
                      g, tkin, tcmb, nc,tkin_ind,tcmb_ind):
    """compute the matrix of transition rates"""

    ###########################################################
    # constructing the matrix
    ###########################################################
    def fill_K_matrix(cr=None,
                      levels=None,
                      collider_density=None,
                      unique_col=None,
                      g=None,
                      ini=None,
                      fin=None, 
                      tkin=None,
                      tkin_ind=None):
        """fill the kij matrix from the collsion rates
        cr:  matrix containing the coolisional
             rates in the linearized form (58,58,50)
        levels: energies for all the rovibrational 
                state of H2 (480 from UCLA database)
        collider_density: density of the coolider for which
                          the collisional reaction rates are
                          given (H in Francois's case)
        unique_col: couples (v,j) for which the collisional  
                    reaction rates are provided
                    (58 in Francois's case)
        g: degeneracies of the levels for which the collisional 
           reaction rates are provided
        ini: couple of all possible initial (v,j) in the provided
             table of reaction rates
        fin: couple of all possible final (v',j') in the provided
             table of reaction rates
        tkin: array containing the temperature at which the reaction
              rates are provided
        .. todo:: improve documentation
        """

        # 1) get the rate coefficient for level i -> i' at
        # temperature T. The collision rate matrix should
        # look like this
        #
        #    k = rate_coeff_for_transition_cr[i, i', T]
        # 
        # 2) get the energy difference between i -> i'
        #
        #    e = fabs(E[i] - E[i']
        #
        # 3) fill the k matrix by repeating 1 for the direct transitions
        #    i.e. i -> i'. and use the detailed balance to
        #    get the rate coeff for the reverse transtion i' -> i

        n = unique_col.size
        K = zeros((n, n), 'f8')

        # mapping the level number to the entry column/row
        # in the K matrix
        ired = zeros(unique_col.max() + 1, 'i4')
        for i, iu in enumerate(unique_col):
            ired[iu] = i

        # filling the K matrix
        for i, ip in zip(ini, fin):

            # indicies of the levels in the matrix K
            ir, irp = ired[i], ired[ip]

            # the difference between the two energy levels
            dE = fabs(levels[i] - levels[ip])

            # the collision rate at tkin[0]. we are useing ir and ip in
            # indexing g since g has the same length as the unqiue levels
            criip = cr[i,ip, tkin_ind]
            criip_rev = cr[ip,i, tkin_ind] * exp( -dE / (kb*tkin[tkin_ind]) )\
                        * g[ir]/g[irp]
            
            K[ir,irp] = criip
            K[irp,ir] = criip_rev
        
        return K


    def fill_AP_matrix(a_eins=None,
                       levels=None,
                       unique_col=None,
                       ini=None,
                       fin=None,
                       tcmb=None,
                       tcmb_ind=None):
        """fill the (A prime)_ij matrix from the lambda radiative transitions
        """
        n = unique_col.size
        AP = zeros( (n, n), dtype = 'f8')

        # mapping the level number to the entry column/row
        # in the K matrix
        ired = zeros(unique_col.max() + 1, 'i4')
        for i, iu in enumerate(unique_col):
            ired[iu] = i

        # filling the AP matrix
        for i, ip in zip(ini, fin):

            # indicies of the levels in the matrix K
            ir, irp = ired[i], ired[ip]

            # the difference between the two energy levels
            dE = fabs(levels[i] - levels[ip])
            
            nu = dE / h # freq in Hz
 
            # the einstein coefficients. we are useing ir and ip in
            # indexing g since g has the same length as the unique levels   
            APiip = (1.0 + ng(h, nu, kb, tcmb[tcmb_ind]))*a_eins[i,ip]
            #APiip = (1.0)*a_eins[i,ip]

            AP[ir,irp] = APiip
        return AP
  
    def fill_ABS_matrix(a_eins=None,
                       levels=None,
                       unique_col=None,
                       g = None,
                       ini=None,
                       fin=None,
                       tcmb=None,
                       tcmb_ind=None):
        """fill the Aij matrix for absorbtion transitions from the lambda
        radiative transitions"""
        n = unique_col.size
        ABS = zeros( (n, n), dtype = 'f8')

        # mapping the level number to the entry column/row
        # in the K matrix
        ired = zeros(unique_col.max() + 1, 'i4')
        for i, iu in enumerate(unique_col):
            ired[iu] = i

        # filling the AP matrix
        for i, ip in zip(ini, fin):

            # indicies of the levels in the matrix K
            ir, irp = ired[i], ired[ip]

            # the difference between the two energy levels
            dE = fabs(levels[i] - levels[ip])
            nu = dE / h # freq in Hz
 
            # the einstein coefficients. we are useing ir and ip in
            # indexing g since g has the same length as the unique levels   
            ABSiip = (g[ir]/g[irp])*ng(h, nu, kb, tcmb[tcmb_ind])*a_eins[i,ip]
            #ABSiip = 0.0
            ABS[ir,irp] = ABSiip
        
        return ABS
    
    def fill_E_matrix(unique_col = None):
        """This is a matrix full of ones (it is a utility matrix"""
        n = unique_col.size
        return ones((n, n))


    n = unique_col.size

    k_mat = fill_K_matrix(cr=cr,
                          levels=en,
                          collider_density=nc,
                          unique_col=unique_col,
                          g=g,
                          ini=ini,
                          fin=fin,
                          tkin=tkin,
                          tkin_ind=tkin_ind)

    ap_mat = fill_AP_matrix(a_eins=a_eins,
                            levels=en,
                            unique_col=unique_col,
                            ini=ini,
                            fin=fin,
                            tcmb=tcmb,
                            tcmb_ind=tkin_ind)

    abs_mat = fill_ABS_matrix(a_eins=a_eins,
                              levels=en,
                              unique_col=unique_col,
                              g = g,
                              ini=ini,
                              fin=fin,
                              tcmb=tcmb,
                              tcmb_ind=tkin_ind)

    E = fill_E_matrix(unique_col = unique_col)

    F = nc * k_mat + ap_mat + abs_mat.T
    diag = eye(n)*dot(F, E)[:, 1]
    offdiag = F.T
    full = -diag + offdiag
    
    return full


def solveEquilibrium(full):
    """solve for the equlibrium population densities.
    
    i.e solving A.x = b
    where full = A
    """

    n = full.shape[0]
    
    # solving directly
    # replacing the first row with the conservation equation
    dndt = zeros((n, 1), 'f8')
    full[0,:] = 1.0
    dndt[0]   = 1.0

    # solving the system A.x = b
    # before solving, we will devide each row by the diagonal
    A, b = full, dndt
    for i in arange(n):
        A[i,:] = A[i,:]/A[i,i]
        x = linalg.solve(A, b)

    #print x.T

    # the fractional population density
    f = x
    return f

#def Tg2Tr(Tg):
#    '''Given the array of kinetic temperatures, it returns the
#     array with the corresponding radiation temperatures (through
#     the redshift parameter z)
#     '''
#     return Tr



def coolingFunction(x, en, eins, ini, fin, unique_col):
    cooling_rate = 0.0

    # mapping the level number to the entry column/row
    # in the K matrix
    ired = zeros(unique_col.max() + 1, 'i4')
    for i, iu in enumerate(unique_col):
        ired[iu] = i

    for i, ip in zip(ini, fin):
        # indicies of the levels in the matrix K
        ir, irp = ired[i], ired[ip]

        # the difference between the two energy levels:
        dE = fabs(en[i] - en[ip])

        # the einstein coefficient connecting i to ip:
        Ai_ip = eins[i,ip]

        #the population of the upper level:
        chi = x[ir]
        cooling_rate += Ai_ip*dE*chi

    return cooling_rate

f=0.
def fit_glover(T):
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
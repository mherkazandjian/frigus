# -*- coding: utf-8 -*-
"""
Read the collisional rates (collision coefficients)
"""

import numpy
from numpy import (loadtxt, arange, int32, zeros, unique, void,
                   ascontiguousarray, dtype, hstack, fabs, exp)
from scipy.constants import Boltzmann, elementary_charge

from scipy import interpolate
import pylab
from IPython.core.debugger import Tracer

from utils import linear_2d_index

def unique_level_pairs(vj):
    """from a list of levels find the list of unique levels and return them.
    :param iterable vj: the list of v levels where each item vj[x] is a level.
    The shape of vj should be (2,n) where n is the number of levels.
    :return: (ndarray) The unique levels. The shape of this array is
    (2,n_unique) where n_unique is the number of unique transitions.

    .. code-block: python

        vj = array([[0, 0, 0, 7, 4, 8, 4, 6, 9, 8, 9, 6, 2, 0, 5, 5, 9, 8, 1],
                    [4, 4, 3, 5, 3, 6, 3, 8, 2, 5, 7, 8, 8, 6, 6, 9, 1, 0, 9]])
        vj_unique = unique_level_pairs(vj)
    """

    assert vj.shape[0] == 2

    a = vj.T
    new_dtype = dtype((void, a.dtype.itemsize * 2))
    b = ascontiguousarray(a).view(new_dtype)
    _, idx = unique(b, return_index=True)
    unique_a = a[idx]
    return unique_a.T


def read_collision_coefficients(fname):
    """parse the collisional data by François. These are the coefficient rates
    K_ij where i > j (so these fill the lower triangular K matrix).

    Supplementary Information for manuscript:
    "Revisited study of the ro-vibrational excitation of H$_2$ by H: Towards a
    revision of the cooling of astrophysical media" by François LIQUE

    The table contains the H2-H collisional rate coefficients

    v, v' : initial and final vibrational state
    j, j' : initial and final rotational state

    v  j  v' j'    k(cm3 s-1) (T) , T= 100 to 5000K by steps of 100K

    0  1  0  0     0.3098E-21  0.1521E-18  0.1791E-16  0.3164E-15.....
    0  2  0  0     0.6561E-13  0.7861E-13  0.9070E-13  0.1266E-12.....

    :param string fname: The path to the ascii data.
    :return: a tuple of 5 elements.
      The first element is a 5D array that holds all the rate coefficients.
      K[v, j, v', j', T_index]

      The second element is a 1D temperature array (T) corresponding to the
      rate coefficients in the 5D array. This array has the same size as the
      last dimension of K (i.e K[0,0,0,0,:].size = T.size

      The third element is a tuple of 4 elements:

         - The first element is a 2D array of shape (2, n_transitions) (ini).
           The columns of this array (v, j = ini) are the initial v and j of
           the transitions for a certain T section. i.e. the number of non-zero
           elements in K for a certain temperature is equal to the number of
           elements in the v or j columns.
           v.size = j.size = where(K[..., 0] > 0)[0].size

         - The second element is the same of the first element but for the
           final transitions (v', j' = fin)

         - The third element is a 2D array of shape (2, n_unique_transitions)
           that are the unique levels involved in all the transitions.

         - The last element is an array of shape (T.size, n_transitions)
           which are the collisional coefficient rates with non-zero values
           for each value of temperature in the T array.
    """

    # the temperature in the data file is not provided explicitly. It
    # it provided as a range. So we refine the temperature and an array
    T = arange(100.0, 5000.1, 100.0)

    # read the data from the original ascii file
    data_read = loadtxt(fname, unpack=True, skiprows=10)
    (v, j, vp, jp), cr = int32(data_read[0:4]), data_read[4:]
    ini = zeros((2, v.size), 'i')
    fin = zeros((2, v.size), 'i')

    # declare the array where the data will be stored
    nv, nj, nvp, njp = v.max()+1, j.max()+1, vp.max()+1, jp.max()+1
    nv_max, nj_max = int(max(nv, nvp)), int(max(nj, njp))
    data = zeros((T.size, nv_max, nj_max, nv_max, nj_max), 'f8')

    # copy the read data into the container array
    for i, cri in enumerate(cr.T):
        data[:, v[i], j[i], vp[i], jp[i]] = cri
        ini[:, i] = v[i], j[i]
        fin[:, i] = vp[i], jp[i]

    # find the unique levels from from the transitions
    unique_levels = unique_level_pairs(hstack((unique_level_pairs(ini),
                                               unique_level_pairs(fin))))

    return data*1e-6, T, (ini, fin, unique_levels, cr)


def compute_lower_to_upper_collision_coefficients(cr,
                                                  ini,
                                                  fin,
                                                  T,
                                                  energy_levels):
    """
    compute the collision coefficients of the lower to upper transitions using
    the detailed balance.

    :param cr: is the collision rates 5D matrix returned by
      read_collision_coefficients.
    :param ini: the initial levels returned by read_collision_coefficients.
    :param fin: the final levels returned by read_collision_coefficients.
    :param T: The temperatures corresponding to the last index of the cr array.
     returned by read_collision_coefficients.
    :param energy_levels: The energy levels

    :return: a copy of the cr matrix with the updated lower to upper
    transitions coefficients.
    """

    retval = cr.copy()

    # convert Boltzmann to eV
    kB = Boltzmann/elementary_charge

    for T_index, T_i in enumerate(T):
        for transition_index in range(ini.shape[1]):
            vi, ji = ini[:, transition_index]
            vf, jf = fin[:, transition_index]

            # the initial and final and the difference of the energies
            # respectively
            energy_i = energy_levels[(energy_levels['v'] == vi)*\
                                     (energy_levels['j'] == ji)]['E'][0]
            energy_f = energy_levels[(energy_levels['v'] == vf)*\
                                     (energy_levels['j'] == jf)]['E'][0]

            dE = fabs(energy_f - energy_i)

            # the degeneracies of the initial and final levels respectively
            g_i, g_f = float(2*ji + 1), float(2*jf + 1)

            # compute the rate from lower to upper and update its value
            cr_upper_2_lower = cr[vi, ji, vf, jf, T_index]
            cr_lower_2_upper = (g_f / g_i)*cr_upper_2_lower*exp(-dE/(kB*T_i))

            retval[vf, jf, vi, ji, T_index] = cr_lower_2_upper

    return retval

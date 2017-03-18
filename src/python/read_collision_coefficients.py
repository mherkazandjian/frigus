# -*- coding: utf-8 -*-
"""
Read the collisional rates (collision coefficients)
"""
from StringIO import StringIO
import numpy
from numpy import (loadtxt, arange, int32, zeros, unique, void,
                   ascontiguousarray, dtype, hstack, fabs, exp)

from astropy import units as u


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
    :return: a tuple of 3 elements.
      The first element is a 5D array that holds all the rate coefficients.
      K[v, j, v', j', T_index] ( in m3/s)

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

    def find_temperature_array():
        """search the header of the data file and return the range of
         temperatures used. The temperature range is assumed to be on the
         10th line of the header"""
        with open(fname) as fobj:
            for line_num, line in enumerate(fobj):
                if line_num == 8:
                    T_values = loadtxt(StringIO(line), delimiter=',') * u.Kelvin
                    break
        return T_values

    T_values = find_temperature_array()

    # read the data from the original ascii file
    data_read = loadtxt(fname, unpack=True, skiprows=10)
    (v, j, vp, jp), cr = int32(data_read[0:4]), data_read[4:]
    ini = zeros((2, v.size), 'i')
    fin = zeros((2, v.size), 'i')

    # declare the array where the data will be stored
    nv, nj, nvp, njp = v.max()+1, j.max()+1, vp.max()+1, jp.max()+1
    nv_max, nj_max = int(max(nv, nvp)), int(max(nj, njp))
    data = zeros((T_values.size, nv_max, nj_max, nv_max, nj_max), 'f8')

    # copy the read data into the container array
    for i, cri in enumerate(cr.T):
        data[:, v[i], j[i], vp[i], jp[i]] = cri
        ini[:, i] = v[i], j[i]
        fin[:, i] = vp[i], jp[i]

    # find the unique levels from from the transitions
    unique_levels = unique_level_pairs(hstack((unique_level_pairs(ini),
                                               unique_level_pairs(fin))))

    # set the units of the data to be returned
    data_with_units = data * (u.cm**3 / u.second)
    cr_with_units = cr * (u.cm**3 / u.second)

    # convert the units to m^3/s
    data_with_units = data_with_units.to(u.m**3 / u.second)
    cr_with_units = cr_with_units.to(u.m**3 / u.second)

    return data_with_units, T_values, (ini, fin, unique_levels, cr_with_units)
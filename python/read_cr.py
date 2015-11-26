# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 01:51:39 2015

@author: carla
"""

import numpy
from numpy import loadtxt, arange, int32, zeros, unique, void
from numpy import ascontiguousarray, dtype, hstack
from numpy import log10
from scipy import interpolate
import pylab
from IPython.core.debugger import Tracer

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

def read_coeff(fname):
    '''parse the  data sent by François:

    Supplementary Information for manuscript:
    "Revisited study of the ro-vibrational excitation of H$_2$ by H: Towards a
    revision of the cooling of astrophysical media" by François LIQUE

    The table contains the H2-H collisional rate coefficients

    v, v' : initial and final vibrational state
    j, j' : initial and final rotational state

    v  j  v' j'    k(cm3 s-1) (T) , T= 100 to 5000K by steps of 100K

    0  1  0  0     0.3098E-21  0.1521E-18  0.1791E-16  0.3164E-15.....
    0  2  0  0     0.6561E-13  0.7861E-13  0.9070E-13  0.1266E-12....

    :param string fname: The path to the ascii data.
    :return: a tuple of two elemtns. 
      The first elemnt is a 5D array that holds all the rate coefficients
      The second elemnt is the temperature corresponding to the rate
      coefficients in the 5D array.

    .. todo:: update the return value
    '''

    # the tempareture in the data file is not provided explicity. It
    # it prived as a range. So we refine the temperature and an array
    T = arange(100.0, 5000.1, 100.0)

    # read the data from the original ascii file
    data_read = loadtxt( fname, unpack=True, skiprows=10)
    (v, j, vp, jp), cr = int32(data_read[0:4]), data_read[4:]
    ini = zeros((2, v.size), 'i')
    fin = zeros((2, v.size), 'i')

    # declare the array where the data will be stored
    nv, nj, nvp, njp = v.max()+1, j.max()+1, vp.max()+1, jp.max()+1
    data = zeros((int(nv), int(nj), int(nvp), int(njp), cr.shape[0]), 'f8')

    # copy the read data into the container array
    for i, cri in enumerate(cr.T):
        data[v[i], j[i], vp[i], jp[i], :] = cri
        ini[:,i] = v[i],j[i]
        fin[:,i] = vp[i],jp[i]

    # find the unique levels from from the transitions
    unique_levels = unique_level_pairs(hstack((unique_level_pairs(ini),
                                               unique_level_pairs(fin))))

    return data, T, v, j, vp, jp, ini, fin, unique_levels

def label(vj_unique):
    '''Given a series of different couples this function assign an integer
    to each of them'''
    couple_label = zeros(len(vj_unique[0])+1,'int')
    for icouple in range(len(vj_unique[0])+1):
        print icouple
        couple_label = [vj_unique[0][icouple],vj_unique[1][icouple]]
#    for i in range(len(v)):
#        vi, ji = v[i], j[i]
#        enh2[vi,ji] = energies[i]
    return couple_label



def interpol(T_grid, vi, ji, vf, jf, data):
    '''Given a grid of temperatures, it returns the interpolated reaction rates for
    a given set of quantum numbers. Moreover, if the provided temperatures are outside
    the range of available temperature in the data, a reaction rate of TINY will be
    assigned to them.

    T_grid must be provided in K
    EXAMPLE OF USAGE:
    T_grid   = numpy.logspace(2, 3, 100) #in K
    rc_francois =  interpol(T_grid ,0,0,1,0, rcinterp)

    def test_call_interpolav(rcu):
        en_grid   = numpy.logspace(10, 16, 3000) #in Hz
        cs_interp =  interpolav_uga(en_grid , 2, 1, 3, csu)
        return cs_interp
    test = test_call_interpolav_uga(cs_uga)
    '''

    data, T, ini, fin, vj_unique, interp_cr = read_cr.read_coeff("Read/Rates_H_H2.dat")

    cr_ret = numpy.zeros( T_grid.size, 'f8' )

    mask = ( T_grid >= T.min() )*( T_grid <= T.max() )
    inds_in  = numpy.where( mask )
    inds_out = numpy.where( numpy.invert(mask) )
                                         # the indexes of the points in the data range

    f_interp = interpolate.interp1d(T, log10(data))
                                         # constructing the interpolation function

    cr_ret[ inds_in  ] = 10.0**f_interp( T_grid[ inds_in ] )
                                         # interpolating the value of the points
                                         # within the range

    cr_ret[ inds_out ] = 1.e-50
                                         # setting the value of the points outside
                                         # the range to a tiny number
#
#    pylab.loglog( en*29979245852.398716, cs/1e20, '*' )
#    pylab.loglog( en_grid*29979245852.398716, cs_ret, 'o' )
#    pylab.show()

    return rc_ret

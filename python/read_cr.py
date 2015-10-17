# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 01:51:39 2015

@author: carla
"""

import numpy
from numpy import loadtxt, genfromtxt, arange, int32, zeros
import pylab
from IPython.core.debugger import Tracer

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
    '''

    # the tempareture in the data file is not provided explicity. It
    # it prived as a range. So we refine the temperature and an array
    T = arange(100.0, 5000.1, 100.0)

    # read the data from the original ascii file
    data_read = loadtxt( fname, unpack=True, skiprows=10)
    (v, j, vp, jp), cr = int32(data_read[0:4]), data_read[4:]

    # declare the array where the data will be stored
    nv, nj, nvp, njp = v.max()+1, j.max()+1, vp.max()+1, jp.max()+1
    data = zeros((int(nv), int(nj), int(nvp), int(njp), cr.shape[0]), 'f8')

    # copy the read data into the container array
    for i, cri in enumerate(cr.T):
        data[v[i], j[i], vp[i], jp[i], :] = cri

    return data, T
#

data, T = read_coeff('../../Rates_H_H2.dat')

print 'done'

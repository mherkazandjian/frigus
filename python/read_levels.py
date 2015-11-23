# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 01:51:39 2015

@author: carla
"""

import numpy 
from numpy import loadtxt, genfromtxt
import pylab

def read_levels( fname ):
    '''parse the  data of
     
    :rtype : object
            http://www.physast.uga.edu/ugamop/energies/H2Xvjlevels
            
    to numpy arrays and provide utilities to get information from the data.       
        # Rovibrational energy levels for H_2
      # v  J        -binding energy
      #           a.u.            cm-1   
        0   0  0.16456626E+00  0.36118114E+05
        0   1  0.16402653E+00  0.35999657E+05
        0   2  0.16295202E+00  0.35763830E+05
        0   3  0.16135249E+00  0.35412773E+05
        0   4  0.15924216E+00  0.34949609E+05
        0   5  0.15663934E+00  0.34378356E+05
        0   6  0.15356587E+00  0.33703808E+05
                   .
                   .
                   .
                   .
       13   7  0.23025618E-03  0.50535384E+02
       14   0  0.65673931E-03  0.14413760E+03
       14   1  0.57867630E-03  0.12700475E+03
       14   2  0.42998862E-03  0.94371580E+02
       14   3  0.22670691E-03  0.49756409E+02
    '''

    v, j, a , energies  = loadtxt( fname, unpack = True )  

    v, j = numpy.int32(v), numpy.int32(j)

    #conversion of the energies from cm-1 -> Joule
    energies  = energies*1.98630e-23

    nv_max, nj_max = v.max(), j.max()

    enh2 = numpy.zeros( (nv_max+1,nj_max+1), 'f8' )


#    bottom = 2179.0*1.98630e-23 # having v=0 as 0 for the energies, the bottom of the well
                                # is at -2179.0 cm-1 + convertion into Joule

    # filling the matrix enh2 with the values read from the text file
    for i in range(len(v)):
        vi, ji = v[i], j[i]
        enh2[vi,ji] = energies[i]

#        depth = enh2[0,0] + bottom
#        enh2[vi,ji]  =  enh2[vi,ji] - depth

    return enh2
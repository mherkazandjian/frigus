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

    return data*1e6, T, ini, fin, unique_levels



import numpy 
import pylab
import matplotlib.pyplot as plt

class uga(object):
    '''parse the  data by Wrathmall and Flower contained in the files
       H.ortho_H2 and H.para_H2
       to numpy arrays and provide utilities to get information from the data.
       
       A block of data is defined everything between :

       T = 100 K

       ( 0, 1) ( 0, 3) ( 0, 5) ( 0, 7) ( 1, 1) ( 1, 3) ( 0, 9) ( 1, 5) ( 0,11) ( 1, 7) ( 2, 1) ( 2, 3) ( 1, 9) ( 0,13) ( 2, 5)     
       ( 1,11) ( 2, 7) ( 3, 1) ( 0,15) ( 3, 3) ( 2, 9) ( 1,13) ( 3, 5) ( 3, 7)
       8.6D-12 7.8D-14 1.1D-15 6.0D-18 6.6D-18 1.7D-18 2.2D-19 1.9D-19 1.1D-20 1.8D-20 2.1D-19 1.2D-19 5.4D-20 2.0D-20 5.7D-20 ...
       ...................................................................................


       T = 300 K

       ( 0, 1) ( 0, 3) ( 0, 5) ( 0, 7) ( 1, 1) ( 1, 3) ( 0, 9) ( 1, 5) ( 0,11) ( 1, 7) ( 2, 1) ( 2, 3) ( 1, 9) ( 0,13) ( 2, 5)  
       ( 1,11) ( 2, 7) ( 3, 1) ( 0,15) ( 3, 3) ( 2, 9) ( 1,13) ( 3, 5) ( 3, 7)
       6.6D-10 1.4D-13 7.3D-15 5.2D-17 1.6D-16 5.6D-17 5.7D-18 1.4D-17 1.0D-18 2.9D-18 6.3D-18 3.8D-18 1.0D-18 2.2D-19 1.7D-18 ...
       ...................................................................................

    
            #!  nei=0, nef=2, vi= 0, ji= 0
            #!  et=118421.5
            #!  np=   2044
             118441.3,1.053E-06
             118455.3,9.454E-06
             118469.4,4.419E-05
                   .
                   .
                   .
                   .
             980392.2,4.963E-15
             990099.0,2.190E-15
            1000000.0,1.623E-16
    '''
    def __init__(self, fname, tiny=1e-70):

        self.fname = fname
        self.fd    = None
        self.tiny  = tiny
        self.eof   = False

        self.data  = None  #: all the crossection data as a function of energy
        self.datai = None  #: an ndarray of object of interpolation function which has the same shape of self.data
        self.v     = None  #: all the read v  levels
        self.j     = None  #: all the read j  levels

        self.read_data()

    def read_data(self):
        '''read all the data to temporary storage in self.data. The fromat 
        of the data can be transformed into more efficient storage by calling
        the method self.optmize_data
        '''
        
        # temporary lists to store the energies, crossection, v, j for all the block
        en_cs, v, j, ef = [], [], [], []
        
        self.fd = open(self.fname, 'r')

        block = 0
        while self.eof == False:

            print '.',

            header, data = self.read_block()

            if header != None:

                nei, nef, vi, ji, et, np = header

                en_cs.append( data )
                
                ef.append( nef )
                v.append( vi )
                j.append( ji )

            block += 1

            if header == None:
                break
        print '\nread data file\n\t%s' % self.fname
        
        self.optimize_data(en_cs, v, j, ef)
        
    def optimize_data(self, en_cs, v, j, ef):
        '''store the data in an efficient way. convert en_cs into a 2D matrix
        of numpy array objects, and v and j into numpy arrays.
        '''
        
        self.v   = numpy.array( v , 'i' )
        self.j   = numpy.array( j , 'i' )
        self.ef  = numpy.array( ef, 'i' )
        
        # create the 3D array of numpy array objects
        nv   = numpy.unique( self.v ).size
        nj   = numpy.unique( self.j ).size
        nef  = numpy.unique( self.ef).size
        
        self.data = numpy.ndarray( (self.ef.max() + 1,
                                    self.v.max()  + 1, 
                                    self.j.max()  + 1),
                                   'object' )

        self.data[:,:,:] = None

        # populating the data
        for i, item in enumerate(en_cs):
            self.data[ self.ef[i], self.v[i], self.j[i] ] = item 
            #print item[0][0],item[1][0]
            #conversion energy units cm-1 -> Hz
#            item[:] = [x*29979245852.398716 for x in item[0]]
            #conversion cross-sections units A**2 -> m**2
#            item[:] = [y*1.e-20 for y in item[1]]

            '''
            plt.plot(item[0][0],item[1][0], '-' )
            plt.yscale('log')
            plt.xscale('log')
            pylab.xlabel('energies [Hz]')
            pylab.ylabel('cross-sections [m**2]')
            pylab.show()
            '''            
   
#    def skip_header(self):
#        ''' skip reading the header at the top of the file'''
#        self.fd.next()
#        self.fd.next()
#        self.fd.next()
#        self.fd.next()
    
    def info(self):
        '''print a summary of the read info'''

        unq_ef = numpy.unique( self.ef ) 
        print 'unique final electronic states:'
        print unq_ef

        unq_v = numpy.unique( self.v ) 
        print 'unique v levels:'
        print unq_v

        unq_j = numpy.unique( self.j )
        print 'unique j levels:'
        print unq_j

        print 'number of unique (nef,v,j) levels:'
        print '\t %d %d %d' %  ( unq_ef.size, unq_v.size, unq_j.size )

        # computing the amount of memory consumed by the data
        nBytes = 0
        for i in numpy.arange(self.data.shape[0]):
            for j in numpy.arange(self.data.shape[1]):
                for k in numpy.arange(self.data.shape[2]):
                    if type(self.data[i,j,k]) != type(None):
                        nBytes += self.data[i,j,k].nbytes

        print 'total number of bytes in memory:%.2f MB' %\
                  ( numpy.float64(nBytes) / 1024.0**2.0 )

        # plot the available transitions irrespective of the fina electronic state)
#        pylab.plot( self.v, self.j, 'o' )
#        pylab.xlabel('v')
#        pylab.ylabel('j')
#        pylab.show()


    def read_block_header(self):
        '''read the header of a block
           T = float            
           and return
           T
           as
           f
        '''
        try:
            hd1 = self.fd.next().replace(' K','').replace('\n','')
        except:
            print 'end of file reached'
            return None

        ''' T = 300 K'''
        hd1 = self.fd.next().replace(' K','').replace('\n','')
        T =  numpy.float64(hd1.split('=')[1])

        '''  ( 0, 0) ( 0, 2) ( 0, 4) ( 0, 6) ( 0, 8) ( 1, 0) ( 1, 2) ( 1, 4) ( 0,10)'''
        hd2 = self.fd.next().replace(' K','').replace('\n','')
        T =  numpy.float64(hd1.split('=')[1])

        return T
        
    def read_block(self):
        '''reads one block of data. When an eof is encountered the attribute eof
        is changed to True. The header is returned as a list and the data of the
        block is returned as a n x n numpy float array.
        ''' 
        
        header = self.read_block_header()
        if header == None:
            self.eof == True
            return None, None 
        
        T  = header
        self.fd.next()
     
        data = numpy.zeros( (2,np), 'f8')
        
        for i in numpy.arange(np):
            
            try:
                en, cs = self.fd.next().replace('\n','').split(',')
            except:
                self.eof == True
                break
            
            data[0,i] = numpy.float64(en)
            
            if '*' in cs:
                data[1,i] = self.tiny
            else:
                data[1,i] = numpy.float64(cs)

        return header, data
    #
#

'''
# definig the reader object
reader = uga('/home/mher/tmp/foo/H2XBCphotodiss.cs')

# showing info
reader.info()

# get the data of a certain transition
en, cs = reader.data[4,1,2]
'''

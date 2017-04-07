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



class Reader(object):
    """parse the  data by Flower and Roueff contained in the file
       flower_roueff_data.dat downloaded from http://massey.dur.ac.uk/drf/HD_H

       A block of data (for a certain temperature) is defined as everything
       between:

       100
       (0,0) (0,1) (0,2) (0,3) (0,4) (0,5) (0,6) (0,7) (0,8) (0, 9)
       1.0D-09 1.5D-11 5.4D-13 4.0D-14 7.8D-15 6.1D-16 2.4D-17 8.7D-18 2.2D-18 2.3D-19
       ...
       ...

    .. code-block:: python

        reader = read_cr_lipovka.Reader('path_to_the_file')

        # print all the collision rates for the transition (v,j) -> (vp,jp) for
        # all the temperatures
        print reader.data[v, j, vp, jp, :]
    """
    def __init__(self, fname, tiny=1e-70):

        self.fname = fname
        self.fd    = None
        self.tiny  = tiny
        self.eof   = False

        self.data  = None  #: the collision rates for all the transitions for all the tempratures
        self.datai = None  #: an ndarray of object of interpolation function which has the same shape of self.data
        self.v     = None  #: all the read v  levels
        self.j     = None  #: all the read j  levels
        self.ini   = None  #: the initial v,j of all the transitions
        self.fin   = None  #: the final v,j of all the transitions
        self.tkin  = None  #: the temperatures at which the collisional data are given

        self.read_data()

    def read_data(self):
        '''read all the data to temporary storage in self.data. The format 
        of the data can be transformed into more efficient storage by calling
        the method self.optmize_data
        '''
        
        # temporary lists to store kinetic temperatures, levels and rates for all the block
        tkin, levels, rates = [], [], []

        with open(self.fname, 'r') as f:

            linesold = f.readlines()
            lines = []
            for line_no, line in enumerate(linesold):
                if line.isspace() is False and line_no >= 12:
                    lines.append(line)
            raw_data = ''.join(lines)

        cr, ini, fin, tkin = self.parse_data(raw_data)
        self.optimize_data(cr, ini, fin, tkin)

    def parse_data(self, raw_data):
        """parses the read data into blocks, one block for each temperature"""

        # split the raw data into blocks
        blocks = raw_data.split('T =')[1:]

        # get the needed info from each block
        blocks_parsed = []

        tkin = zeros(len(blocks), 'f8')
        cr_tkin = []
        ini_tkin = []
        fin_tkin = []

        for ib, block in enumerate(blocks):
            ini, fin, t, cr = self.parse_block(block)
            tkin[ib] = t
            cr_tkin.append(cr)
            ini_tkin.append(ini)
            fin_tkin.append(fin)

        return cr_tkin, ini_tkin, fin_tkin, tkin

    def parse_block(self, block):
        """parse a block of data and return the temperature, the levels
        and the collision rates"""
        lines = iter(block.split('\n'))
        # get the temperature
        T = lines.next()
        assert 'K' in T
        tkin = T.replace('K', '')
        print T, tkin

        # get the header (levels)
        levels = lines.next().replace(' ','').replace(')(',':')[1:-1].split(':')
        v , j = [] , []
        for level in levels:
            v.append(int32(level.split(',')[0]))
            j.append(int32(level.split(',')[1]))

        ini = numpy.repeat(numpy.vstack((v,j)), len(v), axis=1)
        fin = numpy.tile([v,j], len(v))

        cr = zeros((int(len(v)), int(len(j)), int(len(v)), int(len(j))), 'f8')
        for initial, line in enumerate(lines):
            if len(line.strip()) > 0:
               data = numpy.float64(line.replace('D','E').strip().split())    
               for final, (vp, jp) in enumerate(zip(v,j)):
                   cr[v[initial], j[initial], vp, jp] = data[final]
                     
               print data[0]

        return ini, fin, tkin, cr
    
    def optimize_data(self, cr, ini, fin, tkin):
        """store the data in an efficient way. convert en_cs into a 2D matrix
        of numpy array objects, and v and j into numpy arrays.

        .. todo:: update this docstring
        .. todo:: update also the parameter list
        """

        # set these attributes from the parsed data that is returned by
        # self.parse_data
        #
        # self.data
        # self.datai
        # self.v
        # self.j
        # self.ini
        # self.fin
        # self.tkin
        raise NotImplementedError('not implelented yet')

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



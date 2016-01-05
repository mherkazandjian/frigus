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
        '''read all the data to temporary storage in self.data. The format 
        of the data can be transformed into more efficient storage by calling
        the method self.optmize_data
        '''
        
        # temporary lists to store kinetic temperatures, levels and rates for all the block
        tkin, levels, rates = [], [], []

        with open(self.fname, 'r') as f:

            linesold = f.readlines()
            lines = []
            for line in linesold:
                if line.isspace() is False:
     		   lines.append(line)
            raw_data = ''.join(lines)

        self.parse_data(raw_data)

        # self.optimize_data(en_cs, v, j, ef)

    def parse_data(self, raw_data):
        """parses the read data into blocks, one block for each temperature"""

        # split the raw data into blocks
        blocks = raw_data.split('T =')[1:]

        # get the needed info from each block
        blocks_parsed = []
        for block in blocks:
            
            lines = block.split('\n')

            # get the temperature
            T = lines[0]
            assert 'K' in T
            print T

            # get the header (levels)
            levels = lines[1].replace(' ','').replace(')(',':')[1:-1].split(':')
            v , j = [] , []
            for level in levels:
		v.append(int32(level.split(',')[0]))
                j.append(int32(level.split(',')[1]))        
            
            Tracer()() 
	    #
            # get the data (put them in a 2D numpy array)
            #
            #            data = lines[2:24]
            #for ini in levels:
#		for fin in levels:
            #		    k(v[ini],j[ini],v[fin],j[fin]) = lines[2+ini]
            
            #
            # .. todo:: put the parsed data, header and T in blocks_parsed
            #           as a dictionary. For example, the extracted data
            #           should be accessed via:
            #
            #              T = block_parsed[0]['T']
            #              levels = block_parsed[0]['levels']
            #              cr = block_parsed[0]['cr']
	Tracer()()
        pass

    def parse_block(self, block):
        """parse a block of data and return the temperature, the levels
        and the collision rates"""

        #
        # put the stuff
        #   for block in blocks:
        #  
        # in the above method, into this method
        pass
    
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



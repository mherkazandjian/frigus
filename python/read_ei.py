# -*- coding: utf-8 -*-
"""
Created on Wed Jan 28 01:51:39 2015

@author: carla
"""

import numpy
from numpy import zeros
from IPython.core.debugger import Tracer

def read_einstein():
    '''parse the data by Simbotin:

     0   j''=j'
                  14              13              12              11              10               9         .....     1
     1    0   0.6094641D-13   0.1367911D-12   0.1860857D-12   0.1826816D-12   0.9029585D-13   0.8476571D-14   .....
    '''

    # defining the greatest rotational number for H2, starting from j = 0 (tot. rotational numbers: 32)
    jmax = 31

    # defining the greatest vibrational number for H2, starting from v = 0 (tot. vibrational numbers: 15)
    vmax = 14

    #initializing A. The indices are as follows:
    # A[vi, ji, vf, jf]
    A = zeros((vmax+1, jmax+1, vmax+1, jmax+1), 'f8')

    # opening the original ascii file j2jdown to read once for all the
    # maximum number of vibrational levelsfor each rotational level
    # for H2
    f = open('Read/j2jdown', 'r')
    read = f.readlines()
    f.close()
    ivmax = numpy.array(read[2].split(), dtype = 'i')


    def read_j2jdown(fname):

        #initializing Aj2j_down. The indices are as follows:
        # A[vi, ji, vf, jf]
        Aj2j_down = zeros((vmax+1, jmax+1, vmax+1, jmax+1), 'f8')

        # opening the original ascii file
        f = open(fname, 'r')
        read = f.readlines()
        f.close()

        # transition in j = k * |2| i.e  j'' = j' - 2
        k = read[0].split()[0]
        
        # start reading the blocks (after the 3rd line)
        # for each value of j one block (of v transitions) is read
        index = 3
        for j in range(2, jmax+1, 1):
            # initial vibrational level for the columns of the block
            vis = numpy.array(read[index].split(), dtype = 'i')

            # processing the block data one row at time
            for row_num, ivi in enumerate(range(0, ivmax[j]+1)):
                # the inital rotational level 
                ji = int(read[index+ivi+1].split()[0])

                # the final vibrational level
                vf = int(read[index+ivi+1].split()[1])

                # the einstein coefficients for the whole row whose
                # jf = ji -2
                # vi = vis
                es = numpy.array(read[index+ivi+1].replace('D', 'E').split()[2:])

                # A[vi, ji, vf, jf]
                for col_index, vi in enumerate(vis[0:(es.size-1)]):
                    Aj2j_down[vi, ji, vf, ji-2] = es[col_index] 
                #print ji,vf, a
                #Tracer()()

            # incrementing the index of the line to move it to the
            # next block
            index = index + ivmax[j]+1+1
        #
        
        # read the saved array that has been quality checked and compare
        # with what was read above
        A_checked = numpy.load('A_j2jdown_checked.npy')
        assert numpy.fabs(1.0 - Aj2j_down.sum() / A_checked.sum()) == 0.0

        return Aj2j_down
    #

    def read_j2j(fname):

        #initializing Aj2j. The indices are as follows:
        # A[vi, ji, vf, jf]
        Aj2j = zeros((vmax+1, jmax+1, vmax+1, jmax+1), 'f8')

    
        # opening the original ascii file
        f = open('Read/j2j', 'r')
        read= f.readlines()
        f.close

        # transition in j = k * |0| i.e  j'' = j'
        k = read[0].split()[0]

        # start reading the blocks (after the 1st line)
        # for each value of j one block (of v transitions) is read
        index = 1
        for j in range(1, jmax+1, 1):
            # initial vibrational level for the columns of the block
            vis = numpy.array(read[index].split(), dtype = 'i')

            # processing the block data one row at time
            for row_num, ivi in enumerate(range(0,ivmax[j])):
                # the inital rotational level
                ji = int(read[index+ivi+1].split()[0])

                # the final vibrational level
                vf = int(read[index+ivi+1].split()[1])

                # the einstein coefficients for the whole row whose
                # jf = ji -2
                # vi = vis
                es = numpy.array(read[index+ivi+1].replace('D', 'E').split()[2:])

                # A[vi, ji, vf, jf]
                for col_index, vi in enumerate(vis[0:(es.size-1)]):
                    Aj2j[vi, ji, vf, ji] = es[col_index]
                #print ji,vf, a
                #Tracer()()

            # incrementing the index of the line to move it to the
            # next block
            index = index + ivmax[j]+1

        # # read the saved array that has been quality checked and compare
        # # with what was read above
        # A_checked = numpy.load('A_j2j_checked.npy')
        # assert numpy.fabs(1.0 - Aj2j.sum() / A_checked.sum()) == 0.0

        return Aj2j

    def read_j2jup(fname):

        #initializing Aj2j_up. The indices are as follows:
        # A[vi, ji, vf, jf]
        Aj2j_up = zeros((vmax+1, jmax+1, vmax+1, jmax+1), 'f8')

        # opening the original ascii file
        f = open('Read/j2jup', 'r')
        read = f.readlines()
        f.close()

        # transition in j = k * |2| i.e  j'' = j'+2
        k = read[0].split()[0]

        # start reading the blocks (after the 1st line)
        # for each value of j one block (of v transitions) is read
        index = 1
        for j in range(0, jmax-1, 1):
            # initial vibrational level for the columns of the block
            vis = numpy.array(read[index].split(), dtype = 'i')

            # processing the block data one row at time
            for row_num, ivi in enumerate(range(0,min(ivmax[j+2],ivmax[j]-1)+1)):
                # the inital rotational level
                ji = int(read[index+ivi+1].split()[0])

                # the final vibrational level
                vf = int(read[index+ivi+1].split()[1])

                # the einstein coefficients for the whole row whose
                # jf = ji -2
                # vi = vis
                es = numpy.array(read[index+ivi+1].replace('D', 'E').split()[2:])

                # A[vi, ji, vf, jf]
                for col_index, vi in enumerate(vis[0:(es.size-1)]):
                    Aj2j_up[vi, ji, vf, ji] = es[col_index]
               #print ji,vf, a
               #Tracer()()

            # incrementing the index of the line to move it to the
            # next block
            index = index + (min(ivmax[j+2],ivmax[j]-1)+1)+1

        # # read the saved array that has been quality checked and compare
        # # with what was read above
        # A_checked = numpy.load('A_j2jup_checked.npy')
        # assert numpy.fabs(1.0 - Aj2jup.sum() / A_checked.sum()) == 0.0

        return Aj2j_up

    Aj2j_down = read_j2jdown('Read/j2jdown')
    Aj2j = read_j2j('Read/j2j')
    Aj2j_up = read_j2jup('Read/j2jup')
    A = Aj2j_down + Aj2j + Aj2j_up

    return A


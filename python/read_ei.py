# -*- coding: utf-8 -*-
"""
read the Einstein coefficients
"""

import numpy
from numpy import zeros
from IPython.core.debugger import Tracer

def read_einstein():
    """read the data provided by Simbotin from multiple files and returns the A
     matrix for transitions (v', j') -> (v'', j'') with the limitation that
     delta j i.e |j'' - j'| = 0 or 2

     0   j''=j'
                  14              13              12              11              10               9         .....     1
     1    0   0.6094641D-13   0.1367911D-12   0.1860857D-12   0.1826816D-12   0.9029585D-13   0.8476571D-14   .....

    The data is read from the following files:
          Read/j2jdown
          Read/j2j
          Read/j2jup

    :return: A 4D matrix holding the A coefficients. A[v', j', v'', j'']

    .. code-block:: python

        # get the data as a 4D matrix
        A = read_einstein()

        # get the Einstein coefficient for the transition (v'=3, j'=9) to
        # (v''=0, j''=18)
        print(A[3, 9, 0, 18])
    """

    # defining the greatest rotational number for H2, starting from j = 0
    # (tot. rotational numbers: 32)
    jmax = 31

    # defining the greatest vibrational number for H2, starting from v = 0
    # (tot. vibrational numbers: 15)
    vmax = 14

    # initializing A. The indices are as follows:
    # A[vi, ji, vf, jf]
    A = zeros((vmax+1, jmax+1, vmax+1, jmax+1), 'f8')

    # opening the original ascii file j2jdown to read once for all the
    # maximum number of vibrational levels for each rotational level
    # for H2
    f = open('Read/j2jdown', 'r')
    read = f.readlines()
    f.close()
    ivmax = numpy.array(read[2].split(), dtype='i')


    def read_j2jdown(fname):
        """read the Einstein coefficients from the file j2jdown.

        in the snippet below:
        - the first row is the final vibrational level v''
        - the first column is the initial rotational level j'
        - the second column is the initial vibrational level v'
        - the final rational level can be obtained from j' (the first column)
          subtracted by 2

           j'   v'| v''      14              13              12              11              10               9               8               7               6               5               4               3               2               1               0
           2    0 |     0.2237387D-13   0.7602406D-13   0.1977334D-12   0.5213964D-12   0.1451901D-11   0.4314544D-11   0.1378568D-10   0.4779910D-10   0.1820544D-09   0.7742611D-09   0.3770442D-08   0.2158490D-07   0.1274956D-06   0.2526477D-06   0.2941861D-10
           2    1 |     0.9136254D-12   0.2921039D-11   0.6942250D-11   0.1660788D-10   0.4217889D-10   0.1154684D-09   0.3432847D-09   0.1116306D-08   0.4011464D-08   0.1621802D-07   0.7435754D-07   0.3193196D-06   0.3682392D-06   0.2785271D-10
        .. todo:: add documentation

        :param fname: path to the file 'j2jdown'
        :return: A 4D matrix containing the Einstein coefficients

        .. code-block:: python

                .. todo:: add documentation
        """

        # initializing Aj2j_down. The indices are as follows:
        # A[vi, ji, vf, jf]
        Aj2j_down = zeros((vmax+1, jmax+1, vmax+1, jmax+1), 'f8')

        # opening the original ascii file
        with open(fname, 'r') as fobj:
            read = fobj.readlines()

        # transition in j = k * |2| i.e  j'' = j' - 2
        k = read[0].split()[0]
        
        # start reading the blocks (after the 3rd line)
        # for each value of j one block (of v transitions) is read
        index = 3
        for j in range(2, jmax+1, 1):
            # initial vibrational level for the columns of the block
            vis = numpy.array(read[index].split(), dtype='i')

            # processing the block data one row at time
            for row_num, ivi in enumerate(range(0, ivmax[j]+1)):
                # the initial rotational level
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
                # print ji,vf, a
                # Tracer()()

            # incrementing the index of the line to move it to the
            # next block
            index = index + ivmax[j] + 1 + 1

        # read the saved array that has been quality checked and compare
        # with what was read above
        # A_checked = numpy.load('A_j2jdown_checked.npy')
        # assert numpy.fabs(1.0 - Aj2j_down.sum() / A_checked.sum()) == 0.0

        return Aj2j_down
    #

    def read_j2j(fname):
        """
        .. todo:: add doc
        :param fname:  .. todo:: add doc
        :return:       .. todo:: add doc
        """

        # initializing Aj2j. The indices are as follows:
        # A[vi, ji, vf, jf]
        Aj2j = zeros((vmax+1, jmax+1, vmax+1, jmax+1), 'f8')

        # opening the original ascii file
        with open(fname, 'r') as fobj:
            read = fobj.readlines()

        # transition in j = k * |0| i.e  j'' = j'
        k = read[0].split()[0]

        # start reading the blocks (after the 1st line)
        # for each value of j one block (of v transitions) is read
        index = 1
        for j in range(1, jmax+1, 1):
            # initial vibrational level for the columns of the block
            vis = numpy.array(read[index].split(), dtype='i')

            # processing the block data one row at time
            for row_num, ivi in enumerate(range(0,ivmax[j])):
                # the initial rotational level
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
        """
        .. todo:: add doc
        :param fname:  .. todo:: add doc
        :return:       .. todo:: add doc
        """

        # initializing Aj2j_up. The indices are as follows:
        # A[vi, ji, vf, jf]
        Aj2j_up = zeros((vmax+1, jmax+1, vmax+1, jmax+1), 'f8')

        # opening the original ascii file
        with open(fname, 'r') as fobj:
            read = fobj.readlines()

        # transition in j = k * |2| i.e  j'' = j'+2
        k = read[0].split()[0]

        # start reading the blocks (after the 1st line)
        # for each value of j one block (of v transitions) is read
        index = 1
        for j in range(0, jmax-1, 1):
            # initial vibrational level for the columns of the block
            vis = numpy.array(read[index].split(), dtype='i')

            # processing the block data one row at time
            for row_num, ivi in enumerate(range(0, min(ivmax[j+2], ivmax[j]-1)+1)):
                # the initial rotational level
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


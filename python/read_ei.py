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

    # opening the original ascii file
    file = open('Read/j2jdown', 'r')
    read = file.readlines()

    # tansition in j = k * |2| i.e  j'' = j' - 2
    k = read[0].split()[0]

    # the maximum vibrational level for each rotational level
    ivmax = numpy.array(read[2].split(), dtype = 'i')

    # start reading the blocks (after the 3rd line)
    # for each value of j one block (of v transitions) is read
    index = 3
    for j in range(2, jmax+1, 1):
        # initial vibrational level for the columns of the block
        # .. todo:: rename this to vis
        ii = numpy.array(read[index].split(), dtype = 'i')

        # processing the block data one row at time
        for row_num, ivi in enumerate(range(0, ivmax[j]+1)):
            # the inital rotational level 
            ji = int(read[index+ivi+1].split()[0])

            # the final vibrational level
            vf = int(read[index+ivi+1].split()[1])

            # the einstein coefficients for the whole row whose
            # jf = ji -2
            # vi = ii
            es = numpy.array(read[index+ivi+1].replace('D', 'E').split()[2:])

            # A[vi, ji, vf, jf]
            for col_index, vi in enumerate(ii[0:(es.size-1)]):
                A[vi, ji, vf, ji-2] = es[col_index] 
            #print ji,vf, a
            #Tracer()()

        # incrementing the index of the line to move it to the
        # next block
        index = index + ivmax[j]+1+1
    file.close()
    Tracer()()
    
    # opening the original ascii file
    file = open('Read/j2j', 'r')
    read = file.readlines()

    k = read[0].split()[0]
    #ivmax = numpy.array(read[1].split(), dtype = 'i')
    #
    index = 1
    for j in range(1, jmax+1, 1):
         ii = numpy.array(read[index].split(), dtype = 'i')
         for ivi in range(0,ivmax[j]):
             ji = int(read[index+ivi+1].split()[0])
             vf = int(read[index+ivi+1].split()[1])
             a = numpy.array(read[index+ivi+1].split()[2:])
             print ji,vf, a
         index = index + ivmax[j]+1
    file.close()

    # opening the original ascii file
    file = open('Read/j2jup', 'r')
    read = file.readlines()

    k = read[0].split()[0]
    #ivmax = numpy.array(read[1].split(), dtype = 'i')
    #
    index = 1
    for j in range(0, jmax-1, 1):
         ii = numpy.array(read[index].split(), dtype = 'i')
         for ivi in range(0,min(ivmax[j+2],ivmax[j]-1)+1):
             print 'min(ivmax[j+2],ivmax[j]-1)+1',min(ivmax[j+2],ivmax[j]-1)+1
             ji = int(read[index+ivi+1].split()[0])
             vf = int(read[index+ivi+1].split()[1])
             a = numpy.array(read[index+ivi+1].split()[2:])
             print ji,vf, a
         index = index + (min(ivmax[j+2],ivmax[j]-1)+1)+1
         print index
    file.close()

#            A[k, vf, ivmax[j] - ivi, ji] = numpy.array(read[index+ivi+1].split()[2:])

#        for ivf in range(0, ivmax[j],1):
#            print 'j',j

#int32(data_read[0:2]), data_read[2:]

      # open(1,file='j2jdown')
      # read(1,*)k
      # read(1,*)(jj,j=0,jmax)
      # read(1,*)(ivmax(j),j=0,jmax)
      # do j=2,jmax
      #    read(1,*)(ii,i=ivmax(j),0,-1)
      #    do ivf=0,ivmax(j)
      #       read(1,*)jj,ii,(a(-1,ivf,i,j),i=ivmax(j),ivf,-1)
      #    enddo
      # enddo
      # close(1)
      # open(1,file='j2j')
      # read(1,*)k
      # do j=1,jmax-1
      #    read(1,*)(ii,i=ivmax(j),1,-1)
      #    do ivf=0,ivmax(j)-1
      #       read(1,*)jj,ii,(a(0,ivf,i,j),i=ivmax(j),ivf+1,-1)
      #    enddo
      # enddo
      # close(1)
      # open(1,file='j2jup')
      # read(1,*)k
      # do j=0,jmax-2
      #    read(1,*)(ii,i=ivmax(j),1,-1)
      #    do ivf=0,min0(ivmax(j+2),ivmax(j)-1)
      #       read(1,*)jj,ii,(a(1,ivf,i,j),i=ivmax(j),ivf+1,-1)
      #       do ivi=ivmax(j),ivf+1,-1
      #          if(a(1,ivf,ivi,j).lt.0.d0)then
      #             a(-1,ivi,ivf,j+2)=-a(1,ivf,ivi,j)
      #             a(1,ivf,ivi,j)=0.d0
      #          endif
      #       enddo
      #    enddo
      # enddo
      # close(1)

#    (v, j), A = int32(data_read[0:1]), data_read[1:]

#    # declare the array where the data will be stored
#    nv, nj, nvp, njp = v.max()+1, j.max()+1, vp.max()+1, jp.max()+1
#    data = zeros((int(nv), int(nj), int(nvp), int(njp), cr.shape[0]), 'f8')

#    # copy the read data into the container array
#    for i, cri in enumerate(cr.T):
#        data[v[i], j[i], vp[i], jp[i], :] = cri

#    return data, T


#data = read_einstein('Read/j2j')
read_einstein()
print 'done'


      # open(1,file='j2j')
      # read(1,*)k
      # do j=1,jmax-1
      #    read(1,*)(ii,i=ivmax(j),1,-1)
      #    do ivf=0,ivmax(j)-1
      #       read(1,*)jj,ii,(a(0,ivf,i,j),i=ivmax(j),ivf+1,-1)
      #    enddo
      # enddo
      # close(1)

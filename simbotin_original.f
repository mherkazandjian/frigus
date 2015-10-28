      program toy
c
      implicit real*8(a-h,o-z)
      parameter(jmax=31)
      dimension a(-1:1,0:14,0:14,0:jmax),ivmax(0:jmax)
c
c     This is a small program that uses the subroutine 'table'
c
      call table(a,ivmax)
      write(*,202)
      write(*,*)
      write(*,*)'  enter UPPER state rotational quantum number J ='
      read(*,*)j
      write(*,*)'  enter UPPER state vibrational quantum number v ='
      read(*,*)ivi
      write(*,*)'  ENTER lower STATE VIBRATIONAL QUANTUM NUMBER v ='
      read(*,*)ivf
      write(*,*)
      write(*,202)
 202  format(77('_'))
      if(j.lt.0.or.ivi.lt.0.or.ivf.lt.0
     $     .or.j.gt.jmax.or.ivi.gt.14.or.ivf.gt.14)then
         write(*,*)'wrong input -- program is beeing stopped'
         write(*,202)
         stop
      endif
c      write(*,201)ivi,j,ivf
      write(*,*)
      write(*,*)a(-1,ivf,ivi,j)
      write(*,*)a(0,ivf,ivi,j)
      write(*,*)a(1,ivf,ivi,j)
      write(*,202)
      write(*,*)
c 201  format(9x,'The probabilities A for the 
c     ^           quadrupole transitions')
c     $ ('v=',ivi,i3),('J=',j,i3), '--->', ('v=',ivf,i3),'J')
 101  format(22x,'J = J-2,',5x,'A =',e10.3)
 102  format(22x,'J = J,  ',5x,'A =',e10.3)
 103  format(22x,'J = J+2,',5x,'A =',e10.3)
      end
c==========================================
      subroutine table(a,ivmax)
c
      implicit real*8(a-h,o-z)
      parameter(jmax=31)
      dimension a(-1:1,0:14,0:14,0:jmax),ivmax(0:jmax)
c*****************************************************************************
c     The subroutine reads the quadrupole transition probabilities 'a' from  *
c     3 data-files: 'j2jdown,' 'j2j,' and 'j2jup.' Thus, a(k,ivf,ivi,j) is   *
c     the probability for the transition from an initial (upper) state with  *
c     rotational quantum number 'j' and vibrational quantum number 'ivi' to  *
c     a final (lower) state with vibrational quantum number 'ivf' and        *
c     rotational quantum number equal to j-2 for k=-1, j for k=0, and        *
c     j+2 for k=1.                                                           *
c                                                                            *
c     The array elements of 'a' for which the subscripts do NOT correspond   *
c     to a quadrupole transition are put equal to ZERO.                      *
c
c     The subroutine also gives the maximum of the vibrational number for    *
c     each rotational quantum number j. Thus, for a fixed value j, the       *
c     vibrational quantum number can be v=0,1,2,...,ivmax(j). The maximum    *
c     value for j is jmax=31, in which case there is only one possible       *
c     vibrational state, i.e. ivmax(31)=0.                                   *
c*****************************************************************************
      do j=0,jmax
         do i1=0,14
            do i2=0,14
               do k=-1,1
                  a(k,i2,i1,j)=0.d0
               enddo
            enddo
         enddo
      enddo
      open(1,file='j2jdown')
      read(1,*)k
      read(1,*)(jj,j=0,jmax)
      read(1,*)(ivmax(j),j=0,jmax)
      do j=2,jmax
         read(1,*)(ii,i=ivmax(j),0,-1)
         do ivf=0,ivmax(j)
            read(1,*)jj,ii,(a(-1,ivf,i,j),i=ivmax(j),ivf,-1)
         enddo
      enddo
      close(1)
      open(1,file='j2j')
      read(1,*)k
      do j=1,jmax-1
         read(1,*)(ii,i=ivmax(j),1,-1)
         do ivf=0,ivmax(j)-1
            read(1,*)jj,ii,(a(0,ivf,i,j),i=ivmax(j),ivf+1,-1)
         enddo
      enddo
      close(1)
      open(1,file='j2jup')
      read(1,*)k
      do j=0,jmax-2
         read(1,*)(ii,i=ivmax(j),1,-1)
         do ivf=0,min0(ivmax(j+2),ivmax(j)-1)
            read(1,*)jj,ii,(a(1,ivf,i,j),i=ivmax(j),ivf+1,-1)
            do ivi=ivmax(j),ivf+1,-1
               if(a(1,ivf,ivi,j).lt.0.d0)then
                  a(-1,ivi,ivf,j+2)=-a(1,ivf,ivi,j)
                  a(1,ivf,ivi,j)=0.d0
               endif
            enddo
         enddo
      enddo
      close(1)
 401  format(11x,15(6x,i3,7x))
 402  format(2(2x,i3),15(2x,e14.7))
      return
      end


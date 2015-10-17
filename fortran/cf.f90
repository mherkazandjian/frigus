integer                    :: i
integer                    :: v,j,vp,jp
integer, parameter         :: vimax=3,jimax=18,vfmax=2,jfmax=17
integer, parameter                          :: ntemp=50,ntrans=1653
real*8, dimension(1:ntemp)                  :: tg !gas temperature
real*8, allocatable, dimension(:,:,:,:,:)   :: rr !collisional reaction rates

open (12, file='/home/carla/Science/Francois/H+H2/Rates_H_H2.dat', status = 'unknown')

allocate(rr(0:vimax,0:jimax,0:vfmax,0:jfmax,1:ntemp))

tg=0.d0
rr=0.d0

tg =(/ (i, i=100,5000,100) /)

do i=1,10
   read(12,*) 
enddo
do i=1,ntrans
   read(12,*) v,j,vp,jp, (rr(v,j,vp,jp,it),it=1,ntemp)
enddo

end

subroutine tableh2(a,ivmax)

implicit real*8(a-h,o-z)
parameter(jmax=31)
integer ivmax
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
   do ii1=0,14
      do ii2=0,14
         do k=-1,1
            a(k,ii2,ii1,j)=0.d0
         enddo
      enddo
   enddo
enddo      

open(16,file='../Read/j2jdown')
read(16,*)k
read(16,*) (jj,j=0,jmax)
(read(16,*) (ivmax(j),j=0,jmax))
do j=2,jmax
   read(16,*)(ii,i=ivmax(j),1,-1)
   do ivf=0,ivmax(j)
      read(16,*)jj,ii,(a(-1,ivf,i,j),i=ivmax(j),ivf,-1)
   enddo
enddo
      close(16)
      open(17,file='../Read/j2j')
      read(17,*)k
      do j=1,jmax-1
         read(17,*)(ii,i=ivmax(j),1,-1)
         do ivf=0,ivmax(j)-1
            read(17,*)jj,ii,(a(0,ivf,i,j),i=ivmax(j),ivf+1,-1)
         enddo
      enddo
      close(17)
      open(18,file='../Read/j2jup')
      read(18,*)k
      do j=0,jmax-2
         read(18,*)(ii,i=ivmax(j),1,-1)
         do ivf=0,min0(ivmax(j+2),ivmax(j)-1)
            read(18,*)jj,ii,(a(1,ivf,i,j),i=ivmax(j),ivf+1,-1)
            do ivi=ivmax(j),ivf+1,-1
               if(a(1,ivf,ivi,j).lt.0.d0)then
                  a(-1,ivi,ivf,j+2)=-a(1,ivf,ivi,j)
                  a(1,ivf,ivi,j)=0.d0
               endif
            enddo
         enddo
      enddo
      close(18)
 401  format(11x,15(6x,i3,7x))
 402  format(2(2x,i3),15(2x,e14.7))
      return
      end subroutine tableh2
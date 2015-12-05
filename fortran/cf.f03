!GENERAL PURPOSE INDEXES DEFITION 
integer                             :: i
integer, allocatable, dimension (:) :: v,j,vp,jp
integer                             :: vi,ji,vf,jf

!COLLISIONAL TRANSITION
integer, parameter         :: vimax=3,jimax=18,vfmax=2,jfmax=17 !collisional transition
integer, parameter                          :: ntemp=50,ntrans=1653
real*8, dimension(1:ntemp)                  :: tg      !gas temperature
real*8, allocatable, dimension(:,:,:,:,:)   :: rr      !collisional reaction rates
real*8, allocatable, dimension(:,:,:,:,:)   :: K       !collisional reaction rates complete (with detailed balance + for all the transitions)

!RADIATIVE TRANSITION
integer, parameter                          :: jmax=31 !radiative transition
integer, parameter                          :: vmax=14 !radiative transition
integer, dimension(0:jmax)                  :: ivmax   !radiative transition
real*8, allocatable, dimension(:,:,:,:)     :: a       !radiative transition

!LEVEL POPULATION
integer, parameter                          :: nlev = 301 ! number of rovibrational levels according Stancil
real*8, allocatable, dimension(:,:)         :: nvj     !level population

!ENERGY LEVELS
real*8, dimension(:,:), allocatable         :: en

open (20, file='/home/carla/Science/Francois/H+H2/Rates_H_H2.dat', status = 'unknown')
open (21, file='Read/H2Xvjlevels.cs', status = 'unknown')


allocate(rr(0:vimax,0:jimax,0:vfmax,0:jfmax,1:ntemp))
allocate(K(0:vmax,0:jmax,0:vmax,0:jmax,1:ntemp))
allocate(a(-1:1,0:14,0:14,0:jmax))
allocate(en(0:vmax,0:jmax))

allocate(v(1:ntrans),j(1:ntrans),vp(1:ntrans),jp(1:ntrans))
allocate(nvj(0:vimax,0:jimax))
!initializing variables
rr = 0.d0
K = 0.d0
a = 0.d0
en = 0.d0
nvj = 0.d0
tg = (/ (i, i=100,5000,100) /)



!reading data 
!-energy levels
do i=1,10
   read(21,*)
enddo
do i=1,nlev
   read(21,*) vi,ji, en(vi,ji), b
   write(6,'(4(e10.4))') vi, ji, en(vi,ji)
enddo


!-collisional
do i=1,10 
   read(20,*) 
enddo
do i=1,ntrans
   read(20,*) v(i),j(i),vp(i),jp(i),  &
             (rr(v(i),j(i),vp(i),jp(i),it),it=1,ntemp)
   vi=v(i)
   ji=j(i)
   vf=vp(i)
   jf=jp(i)
!   dE=
   do it=1,ntemp
!      K(vi,ji,vf,jf,it) = dexp(-dE/(kb*temp(it))) * rr(vi,ji,vf,jf,it)
      K(vf,jf,vi,ji,it) = rr(vi,ji,vf,jf,it)
   enddo
enddo
!-radiative
call tableh2(a,ivmax)
!write(6,'(a24,2x,e10.4)') 'a(-1:1,0:14,0:14,0:jmax)',a(1,0,1,0)

!filling the matrix



end







subroutine tableh2(a,ivmax)

implicit real*8(a-h,o-z)
parameter(jmax=31)
integer ivmax
dimension a(-1:1,0:14,0:14,0:jmax),ivmax(0:jmax)
!*****************************************************************************
!     The subroutine reads the quadrupole transition probabilities 'a' from  *
!     3 data-files: 'j2jdown,' 'j2j,' and 'j2jup.' Thus, a(k,ivf,ivi,j) is   *
!     the probability for the transition from an initial (upper) state with  *
!     a final (lower) state with vibrational quantum number 'ivf' and        *
!     rotational quantum number equal to j-2 for k=-1, j for k=0, and        *
!     j+2 for k=1.                                                           *
!                                                                            *
!     The array elements of 'a' for which the subscripts do NOT correspond   *
!     to a quadrupole transition are put equal to ZERO.                      *
!
!     The subroutine also gives the maximum of the vibrational number for    *
!     each rotational quantum number j. Thus, for a fixed value j, the       *
!     value for j is jmax=31, in which case there is only one possible       *
!     vibrational state, i.e. ivmax(31)=0.                                   *
!*****************************************************************************

a=0.d0

open(16,file='Read/j2jdown')
read(16,*)k
read(16,*) (jj,j=0,jmax)
read(16,*) (ivmax(j),j=0,jmax)
do j=2,jmax
   read(16,*)(ii,i=ivmax(j),1,-1)
   do ivf=0,ivmax(j)
      read(16,*)jj,ii,(a(-1,ivf,i,j),i=ivmax(j),ivf,-1)
   enddo
enddo
      close(16)
      open(17,file='Read/j2j')
      read(17,*)k
      do j=1,jmax-1
         read(17,*)(ii,i=ivmax(j),1,-1)
         do ivf=0,ivmax(j)-1
            read(17,*)jj,ii,(a(0,ivf,i,j),i=ivmax(j),ivf+1,-1)
         enddo
      enddo
      close(17)
      open(18,file='Read/j2jup')
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

!PHYSICAL PARAMETERS
real*8, parameter                           :: kb = 1.38064852e-23 !J/K
real*8, parameter                           :: hp = 6.62607004e-34 !J s

!GENERAL PURPOSE INDEXES DEFITION 
integer                                     :: i,l, i0,i1,i2,i3,upper,lower
integer, allocatable, dimension(:)          :: v,j,vp,jp              ! labels for the collisional transitions
integer, allocatable, dimension(:)          :: vrad,jrad,vprad,jprad  ! labels for the radiative transitions with non-zero  with non-zero A
integer                                     :: vi,ji,vf,jf            ! integers for identifying the collisional transitions

!COLLISIONAL TRANSITION
integer, parameter                          :: vimax=3,jimax=18,vfmax=2,jfmax=17 !collisional transition
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

!ENERGY LEVELS
real*8, dimension(:,:), allocatable         :: en
real*8, dimension(:), allocatable           :: ene
integer, allocatable, dimension(:)          :: vl, jl     ! labels for the levels (for the final ordering)
real*8                                      :: dE

!LABELLING TRANSITIONS
integer, dimension(:), allocatable           :: couple1, couple2 !they contain the indexes of the couples in the collisional file
integer                                      :: nradtrans        !number of radiative transitions not equal to 0
integer, dimension(:), allocatable           :: coupler1, coupler2 !they contain the indexes of the couples in the radiative file

!MATRICES FOR THE SYSTEM (nlev * nlev)
real*8, dimension(:,:,:), allocatable        :: KK
real*8, dimension(:,:), allocatable          :: AA
real*8, dimension(:)  , allocatable          :: nvj        !level population

open (20, file='Read/Rates_H_H2.dat', status = 'unknown')
open (21, file='Read/H2Xvjlevels.cs', status = 'unknown')
open (22, file='Read/rad_trans_unique', status = 'unknown')
open (23, file='Read/lev_labels', status = 'unknown')


allocate(rr(0:vimax,0:jimax,0:vfmax,0:jfmax,1:ntemp))
allocate(K(0:vmax,0:jmax,0:vmax,0:jmax,1:ntemp))
allocate(a(-1:1,0:14,0:14,0:jmax))
allocate(en(0:vmax,0:jmax))
allocate(ene(1:nlev))

allocate(v(1:ntrans),j(1:ntrans),vp(1:ntrans),jp(1:ntrans))
allocate(vl(1:nlev),jl(1:nlev))
allocate(couple1(1:ntrans),couple2(1:ntrans))
allocate(nvj(1:nlev))
allocate(KK(1:nlev,1:nlev,1:ntemp),AA(1:nlev,1:nlev))

!initializing variables
rr  = 0.d0
K   = 0.d0
a   = 0.d0
en  = 0.d0
nvj = 0.d0
KK  = 0.d0
nradtrans = 0
AA  = 0.d0
tg = (/ (i, i=100,5000,100) /)



!READING DATA
!-ENERGY LEVELS
do i=1,10
   read(21,*)
enddo
do i=1,nlev
   read(21,*) vl(i),jl(i), en(vl(i),jl(i)), b
   ene(i) = en(vl(i),jl(i))
   write(23,'(3(i3,2x))') i, vl(i),jl(i)
!   write(6,'(i3,2x,i2,2x,i2,2x,e10.4)') i, vl(i), jl(i), en(vl(i),jl(i))
enddo
!conversion cm-1 -->  Joule
en = en*1.98630e-23

!-COLLISIONAL
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
   do l=1,nlev
      if(vi.eq.vl(l)) then
        if(ji.eq.jl(l)) then
          couple1(i) = l
        endif
      endif
      if(vf.eq.vl(l)) then
        if(jf.eq.jl(l)) then
          couple2(i) = l
        endif
      endif
   enddo
enddo
!do i=1,ntrans
!   print*, couple1(i),couple2(i)
!enddo
do i=1,ntrans
   vi=v(i)
   ji=j(i)
   vf=vp(i)
   jf=jp(i)
   dE=abs(en(vi,ji)-en(vf,jf))
!   print*, vi,ji,vf,jf,dE
   do it=1,ntemp
      K(vi,ji,vf,jf,it) = rr(vi,ji,vf,jf,it)
      K(vf,jf,vi,ji,it) = dexp(-dE/(kb*tg(it))) * rr(vi,ji,vf,jf,it)
      KK(couple1(i),couple2(i),it) = rr(vi,ji,vf,jf,it)
      KK(couple2(i),couple1(i),it) = dexp(-dE/(kb*tg(it))) * rr(vi,ji,vf,jf,it)
   enddo
enddo

!do i=1,ntrans
!   do l=1,ntemp 
!  if(K(v(i),j(i),vp(i),jp(i),l).ne.0.d0) print*, K(v(i),j(i),vp(i),jp(i),l)
!   enddo
!enddo


!write(6,'(i3,2x,i3)') couple1,couple2

!do i=1,ntrans
!   do l=1,ntemp 
!  if(KK(couple1(i),couple2(i),l).ne.0.d0) print*, i, KK(couple1(i),couple2(i),l)
!   enddo
!enddo




!-RADIATIVE
call tableh2(a,ivmax)
!write(6,'(a24,2x,e10.4)') 'a(-1:1,0:14,0:14,0:jmax)',a(1,0,1,0)
do i0=0,jmax
   do i1=0,14
      do i2=0,14
         do i3=-1,1
            if(a(i3,i2,i1,i0).ne.0.d0) then
            nradtrans = nradtrans + 1
            write(22,'(4(i2,2x),e10.4)') i1,i0,i2,i0+i3,a(i3,i2,i1,i0)
            endif
         enddo
      enddo
   enddo
enddo
rewind(22)
allocate(coupler1(1:nradtrans),coupler2(1:nradtrans))
allocate(vrad(1:nradtrans),jrad(1:nradtrans))
allocate(vprad(1:nradtrans),jprad(1:nradtrans))

!print*, nradtrans, (jmax+1)*15*15*3


do i=1,nradtrans
   read(22,*)  vrad(i),jrad(i),vprad(i),jprad(i), &
               a(vrad(i),jrad(i),vprad(i),jprad(i))
!   write(6,'(4(i2,2x),e10.4)')  vrad(i),jrad(i),vprad(i),jprad(i), &
!              b   
   vi=vrad(i)
   ji=jrad(i)
   vf=vprad(i)
   jf=jprad(i)
   do l=1,nlev
      if(vi.eq.vl(l)) then
        if(ji.eq.jl(l)) then
          coupler1(i) = l
        endif
      endif
      if(vf.eq.vl(l)) then
        if(jf.eq.jl(l)) then
          coupler2(i) = l
        endif
      endif
   enddo
enddo
!do i=1,nradtrans
!   print*, coupler1(i),coupler2(i)
!enddo
do i=1,nradtrans
   vi=vrad(i)
   ji=jrad(i)
   vf=vprad(i)
   jf=jprad(i)
   AA(coupler1(i),coupler2(i)) = a(vi,ji,vf,jf)
!   write(6,*) vi,ji,vf,jf,a(vi,ji,vf,jf)
enddo

!checking the transitions/coefficients
!for example,  1 -> 1 means from level 1 to level 1 in the list etc
! AA(in, fin) = a(vl(in),jl(in),vl(fin),jl(fin))
!upper  = 102 !in the nlev order
!lower  = 44 !in the nlev order
!write(6,'(a19,i2,a3,i2,a10,i2,a4,i2,a1)') 'transition from (v=',vl(upper),',j=',jl(upper),') to (vf=',vl(lower),',jf=',jl(lower),')'
!write(6,'(a4,e10.4,2x,a3,e10.4,2x,a14,e10.4)'),  'AA =', AA(upper,lower), 'KK=',KK(upper,lower,ntemp),'at gas temp T=',tg(ntemp)

! searching for the transitions that are in common between the collisional and the radiative ones
do i=1,nradtrans
   do l=1,ntrans
      if(couple1(l).eq.coupler1(i)) then
        if(couple2(l).eq.coupler2(i)) then
          write(6,'(i3,2x,i3,2x)') couple1(l),couple2(l)
        endif
      endif
   enddo      
enddo
 




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


do j=0,jmax
   do i1=0,14
      do i2=0,14
         do k=-1,1
            a(k,i2,i1,j)=0.d0
         enddo
      enddo
   enddo
enddo

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

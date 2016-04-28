!PHYSICAL PARAMETERS
real*8, parameter                           :: kb = 1.38064852e-23 !J/K
real*8, parameter                           :: hp = 6.62607004e-34 !J s

!GENERAL PURPOSE INDEXES DEFITION 
integer                                     :: i,l,m,n,i0,i1,i2,i3,upper,lower
integer, allocatable, dimension(:)          :: v,j,vp,jp              ! labels for the collisional transitions
integer                                     :: vi,ji,vf,jf            ! integers for identifying the collisional transitions

!COLLISIONAL TRANSITION
integer, parameter                          :: vimax=3,jimax=18,vfmax=2,jfmax=17 !collisional transition
integer, parameter                          :: ntemp=50,ntrans=1653
real*8, dimension(1:ntemp)                  :: tg      !gas temperature
real*8, allocatable, dimension(:,:,:,:,:)   :: rr      !collisional reaction rates
real*8, allocatable, dimension(:,:,:,:,:)   :: K       !collisional reaction rates complete (with detailed balance + for all the transitions)

!!RADIATIVE TRANSITION
!integer, parameter                          :: jmax=31 !radiative transition
!integer, parameter                          :: vmax=14 !radiative transition
!integer, dimension(0:jmax)                  :: ivmax   !radiative transition
!real*8, allocatable, dimension(:,:,:,:)     :: a       !radiative transition

!LEVEL POPULATION
!integer, parameter                          :: nlev = 301 ! number of rovibrational levels according to Stancil

!!ENERGY LEVELS
!real*8, dimension(:,:), allocatable         :: en
!real*8, dimension(:), allocatable           :: ene
!integer, allocatable, dimension(:)          :: vl, jl     ! labels for the levels (for the final ordering)
!real*8                                      :: dE

!LABELLING TRANSITIONS
integer, dimension(:), allocatable           :: couple1, couple2   !they contain the indexes of the couples in the collisional file
!integer                                      :: nradtrans          !number of radiative transitions not equal to 0
!integer, dimension(:), allocatable           :: coupler1, coupler2 !they contain the indexes of the couples in the radiative file

!MATRICES FOR THE SYSTEM (nlev * nlev)
real*8, dimension(:,:,:), allocatable        :: KK         !collisional transitions (assembled according to the couples vi,ji->vf,jf)
!real*8, dimension(:,:),   allocatable        :: AA         !einstein coefficients (radiative spontaneous transitions)
real*8, dimension(:,:,:), allocatable        :: MM          !containing both radiative and collisional transitions (notation Tinè1998)
real*8, dimension(:)  , allocatable          :: nvj        !level population

open (20, file='Read/Rates_H_H2.dat', status = 'unknown')
open (22, file='Read/rad_trans_unique', status = 'unknown')


allocate(rr(0:vimax,0:jimax,0:vfmax,0:jfmax,1:ntemp))
allocate(K(0:vmax,0:jmax,0:vmax,0:jmax,1:ntemp))
allocate(a(-1:1,0:14,0:14,0:jmax))
!allocate(en(0:vmax,0:jmax))
!allocate(ene(1:nlev))

allocate(v(1:ntrans),j(1:ntrans),vp(1:ntrans),jp(1:ntrans))
!allocate(vl(1:nlev),jl(1:nlev))
allocate(couple1(1:ntrans),couple2(1:ntrans))
allocate(nvj(1:nlev))
allocate(KK(1:nlev,1:nlev,1:ntemp),AA(1:nlev,1:nlev),MM(1:nlev,1:nlev,1:ntemp))

!initializing variables
rr  = 0.d0
K   = 0.d0
a   = 0.d0
!en  = 0.d0
nvj = 0.d0
KK  = 0.d0
nradtrans = 0
AA = 0.d0
MM = 0.d0
tg = (/ (i, i=100,5000,100) /)



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
!units conversion: cm3 s-1 -> m3 s-1
K=K*1.d-6

!check of the coefficients compared to the read data by François
!do i=1,ntrans
!   do l=1,ntemp 
!    if(K(v(i),j(i),vp(i),jp(i),l).ne.0.d0) &
!               write(6,'(4(i2,2x),e12.4)') &
!    v(i),j(i),vp(i),jp(i),K(v(i),j(i),vp(i),jp(i),l)
!   enddo
!enddo


!write(6,'(i3,2x,i3)') couple1,couple2

do i=1,ntrans
   do l=1,ntemp 
  if(KK(couple1(i),couple2(i),l).ne.0.d0) print*, i, KK(couple1(i),couple2(i),l)
   enddo
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
          write(6,'(i3,2x,i3,2x)')  couple1(l),couple2(l)
          write(6,'(a10,2x,i3,2x,i3,2x,i3,2x,i3,2x)') &
                    'radiation:',vrad(i),jrad(i),vprad(i),jprad(i)
          write(6,'(a10,2x,i3,2x,i3,2x,i3,2x,i3,2x)') &
                    'collision:',v(l),j(l),vp(l),jp(l)
        endif
      endif
   enddo      
enddo
 

!solution of the linear system of equations M \cdot n = 0 where n = nvj while the matrix M contains the 
!term of de/excitations due to collisions and radiative de-excitations (only spontaneous transitions included)
!filling M:
do l=1,ntemp
   do n=1, nlev
      do m=1,nlev
      if((ene(n)-ene(m)).gt.0.d0) then !de-excitation from n down to m (En - Em > 0)
            MM(m,n,l) = MM(m,n,l) + nvj(n)*AA(n,m)
          else                         !excitation from n up to m (En - Em < 0)
            MM(m,n,l) = 0.d0
         endif
      enddo
   enddo
enddo






end





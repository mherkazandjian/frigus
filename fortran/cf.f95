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






print*, size(rr)

end

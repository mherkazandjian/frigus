module types_and_parameters

    !PHYSICAL PARAMETERS
    real*8, parameter                           :: kb = 1.38064852e-23 !J/K
    real*8, parameter                           :: hp = 6.62607004e-34 !J s

    !GENERAL PURPOSE INDEXES DEFITION 
    integer                                     :: i,l,m,n,i0,i1,i2,i3,upper,lower
    integer, allocatable, dimension(:)          :: v,j,vp,jp   ! labels for the coll transitions
    integer                                     :: vi,ji,vf,jf ! integers for identifying coll transitions

    ! ENERGY LEVELS
    integer, parameter :: nlev = 301

    ! RADIATIVE TRANSITIONS
    integer, parameter :: jmax=31
    integer, parameter :: vmax=14

    ! COLLISIONAL TRANSITIONS
    integer, parameter :: vimax = 3, jimax = 18, vfmax = 2, jfmax = 17
    integer, parameter :: ntemp = 50, ntrans = 1653

    
    type :: energy_lev
        real*8, dimension(:,:), allocatable         :: en         ! per pair (v,j)
        real*8, dimension(:), allocatable           :: ene        ! labelled according to the order in the file
        integer, allocatable, dimension(:)          :: vl, jl     ! labels for the levels (for the final ordering)
        integer  :: jmax
        integer  :: vmax
    end type energy_lev
    

    type :: radiative_coeffs
        integer, dimension(0:jmax) :: ivmax
        real*8  :: reading(-1:1,0:14,0:14,0:jmax)
        real*8  :: M(1:nlev,1:nlev)
        integer :: ntransrad
        integer, allocatable, dimension(:) :: couple1r,couple2r 
        integer, allocatable, dimension(:) :: vir, jir, vfr, jfr
    end type radiative_coeffs


    type :: collisional_coeffs
        integer  :: ntemp   ! # of temperatures @ which collisional data are given
        real*8, dimension(1:ntemp)  :: temp    ! temperatures @ which collisional data are given
        integer, dimension(1:ntrans)  :: vi, ji, vf, jf ! initial and final rovibrational levels (transitions          
                                                        ! order)
        real*8  :: reading(1:nlev,1:nlev, 1:ntemp)
        real*8  :: detailedbalance(1:nlev,1:nlev, 1:ntemp)
        integer :: couple1(1:ntrans),couple2(1:ntrans)
    end type collisional_coeffs


end module
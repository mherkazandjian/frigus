module types_and_parameters

    !PHYSICAL PARAMETERS
    real*8, parameter                           :: kb = 1.38064852d-23 !J/K Boltzmann constant
    real*8, parameter                           :: hp = 6.62607004d-34 !J s Planck constant
    real*8, parameter                           :: c  = 2.99792458d8   !speed of light m/s

    !GENERAL PURPOSE INDEXES DEFITION 
    integer                                     :: i, l, m, n, i0, i1, i2, i3
    integer                                     :: lower, upper, ini, fin, it
    integer                                     :: vi, ji, vf, jf ! integers for identifying coll transitions

    ! GAS DENSITY AND RADIATION TEMPERATURE
!    integer, parameter                          :: ndensity = 5   ! dimension of the arra of density
    real*8, parameter                           :: Trad = 30.d0        ! radiation temperature in kelvin
!   real*8, parameter, dimension(1:ndensity)    :: nb =              ! baryon density
    
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
        real*8, dimension(:,:), allocatable          :: freq       ! frequencies, obtained as differences between 
                                                                  !              each couple of energy levels
        integer, allocatable, dimension(:)          :: vl, jl     ! labels for the levels (for the final ordering)
        integer  :: jmax
        integer  :: vmax
    end type energy_lev
    

    type :: radiative_coeffs
        integer, dimension(0:jmax) :: ivmax
        real*8  :: reading(-1:1,0:14,0:14,0:jmax)
        real*8  :: M(1:nlev,1:nlev)
        integer :: ntransrad
        integer, allocatable, dimension(:) :: couple1r, couple2r 
        integer, allocatable, dimension(:) :: vir, jir, vfr, jfr
    end type radiative_coeffs


    type :: collisional_coeffs
        integer  :: ntemp   ! # of temperatures @ which collisional data are given
        real*8, dimension(1:ntemp)  :: temp    ! temperatures @ which collisional data are given
        integer, dimension(1:ntrans)  :: vic, jic, vfc, jfc ! initial and final rovibrational levels (transitions          
                                                        ! order)
        real*8  :: reading(0:vimax, 0:jimax, 0:vfmax, 0:jfmax, 1:ntemp)
        integer :: couple1c(1:ntrans),couple2c(1:ntrans)
        real*8  :: matrix(1:nlev,1:nlev, 1:ntemp) ! n.b.: it includes the detailed balance of the read data
    end type collisional_coeffs

    type :: reaction_matrix
        real*8 :: A(1:nlev,1:nlev)        
    end type reaction_matrix

end module
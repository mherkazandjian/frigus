module types_and_parameters

    !PHYSICAL PARAMETERS
    real*8, parameter                           :: kb = 1.3806488d-23  !J/K Boltzmann constant
    real*8, parameter                           :: hp = 6.62606957d-34  !J s Planck constant
    real*8, parameter                           :: c  = 2.99792458d8    !speed of light m/s
    real*8, parameter                           :: q  = 1.602176565d-19  !charge of electron in C
    real*8, parameter                           :: pi = 3.14159265359

    !GENERAL PURPOSE INDEXES DEFITION 
    integer                                     :: i, l, m, n, i0, i1, i2, i3
    integer                                     :: lower, upper, ini, fin, row, col, it
    integer                                     :: vi, ji, vf, jf ! integers for identifying coll transitions
    integer                                     :: id_temp
    integer, parameter                          :: id_temp_test = 1
    real*8                                      :: xm, ym

    ! GAS DENSITY AND RADIATION TEMPERATURE
    integer, parameter                          :: ndensity = 9     ! dimension of the arra of density
    real*8,  parameter                          :: Trad = 0.d0   ! radiation temperature in kelvin
    real*8,  parameter, dimension(1:ndensity)   :: nc = [1.d6, 1.d7, 1.d8, 1.d9, 1.d10, 1.d11, 1.d12, 1.d13,1.d14]       ! baryon density [m-3]
    !real*8,  parameter                          :: nc = 1.d14       ! baryon density [m-3]
    
    ! ENERGY LEVELS
    integer, parameter :: nlev = 301
    integer, parameter :: nlev_lique = 108 ! 58 lique, 108 flower
    integer, parameter :: jmax_lique = 23  ! 18 lique,  23 flower
    integer, parameter :: vmax_lique = 6   ! 3  lique,   6 flower

    ! RADIATIVE TRANSITIONS
    integer, parameter :: jmax = 31
    integer, parameter :: vmax = 14

    ! COLLISIONAL TRANSITIONS
    integer, parameter :: vimax = 6    ! 3 lique,   6 flower
    integer, parameter :: jimax = 23   ! 18 lique, 23 flower
    integer, parameter :: vfmax = 6    !  3 lique,  6 flower
    integer, parameter :: jfmax = 23   ! 17 lique, 23 flower
    integer, parameter :: ntemp = 50   ! 50 lique, 50 flower
    integer, parameter :: ntrans = 5797 ! 1653 lique 5797 flower, tot 1305 24*24+27*27


    type :: energy_lev
        real*8, dimension(:,:), allocatable         :: en         ! per pair (v,j)
        real*8, dimension(:),   allocatable         :: ene        ! labelled according to the order in the file
        real*8, dimension(:,:), allocatable         :: freq       ! frequencies, obtained as differences between 
                                                                  !              each couple of energy levels
        integer, allocatable, dimension(:)          :: vl, jl     ! labels for the levels (for the final ordering)
        integer  :: jmax
        integer  :: vmax
        real*8, dimension(:,:), allocatable         :: en_lique         ! per pair (v,j)
        real*8, dimension(:),   allocatable         :: ene_lique        ! labelled according to the order in the file
        real*8, dimension(:,:), allocatable         :: freq_lique       ! frequencies, obtained as differences between 
                                                                  !              each couple of energy levels
        integer, allocatable, dimension(:)          :: vl_lique, jl_lique     ! labels for the levels (for the final ordering)
        integer, allocatable, dimension(:)          :: vl_flower, jl_flower
        integer  :: jmax_lique
        integer  :: vmax_lique        
    end type energy_lev
    

    type :: radiative_coeffs
        integer, dimension(0:jmax) :: ivmax
        real*8  :: reading(-1:1,0:14,0:14,0:jmax)
        real*8, allocatable, dimension(:,:,:,:)  :: arranging
        real*8  :: M(1:nlev,1:nlev)
        real*8  :: M_lique(1:nlev_lique,1:nlev_lique)        
        integer :: ntransrad
        integer, allocatable, dimension(:) :: couple1r, couple2r 
        integer, allocatable, dimension(:) :: vir, jir, vfr, jfr
    end type radiative_coeffs


    type :: parameter_collisional_coeff
        integer :: ntrans, ntemp, nlev, vimax, jimax, vfmax, jfmax
    end type parameter_collisional_coeff
 
 
!     type :: generic_collisional_coeffs
!         real*8, dimension(:), allocatable :: temp    ! temperatures @ which collisional data are given
!         integer, dimension(:), allocatable  :: vic, jic, vfc, jfc ! in and fin rovibrational levels (transitions order)
!         real*8, dimension(:,:,:,:,:)  :: reading
!         integer, dimension(:) :: couple1c,couple2c
!         real*8, dimension(:,:,:)  :: matrix           ! n.b.: it includes the detailed balance of the read data
!     end type generic_collisional_coeffs
    
    
    
    type :: collisional_coeffs
        integer  :: ntemp   ! # of temperatures @ which collisional data are given
        real*8, dimension(1:ntemp)  :: temp, temp_flower   ! temperatures @ which collisional data are given
        integer, dimension(1:ntrans)  :: vic, jic, vfc, jfc ! initial and final rovibrational levels (transitions          
                                                        ! order)
        real*8  :: reading(0:vimax, 0:jimax, 0:vfmax, 0:jfmax, 1:ntemp)
        integer :: couple1c(1:ntrans),couple2c(1:ntrans)
        !real*8  :: matrix(1:nlev,1:nlev, 1:ntemp) ! n.b.: it includes the detailed balance of the read data
        real*8  :: matrix_lique(1:nlev_lique,1:nlev_lique, 1:ntemp) ! n.b.: it includes the detailed balance of the read data        
    end type collisional_coeffs
    

    type :: reaction_matrix
        real*8 :: A(1:nlev, 1:nlev)
        real*8 :: M(1:nlev_lique, 1:nlev_lique)
        real*8 :: O(1:nlev_lique, 1:nlev_lique)
        real*8 :: D(1:nlev_lique, 1:nlev_lique)
    end type reaction_matrix
    
    
    type :: population
        real*8 :: pop(1:nlev)
        real*8 :: pop_t(1:nlev, 1:ntemp)
        real*8 :: pop_lique(1:nlev_lique)        
        real*8 :: pop_t_lique(1:nlev_lique, 1:ntemp)        
    end type population

end module
module linear_algebra

    use types_and_parameters, only: nlev
    ! LAPACK SOLVER VARIABLES
    integer, parameter :: ndim = nlev
    integer :: i, info
    integer, parameter :: lda = ndim, ldb = ndim ! leading dimension of a
                                                 ! leading dimension of b    
    integer, parameter :: nrhs = 1 ! number of right hand sides in b
    integer, dimension(1:ndim) :: ipiv


    contains
    
    subroutine sparsity_calc(b, sparsity)
        use types_and_parameters, only: reaction_matrix
        type(reaction_matrix) :: b
        integer               :: i, j
        real*8                :: sparsity       
        sparsity=0.d0
        !print*, nlev**2
        do i = 1, nlev
           do j = 1, nlev
              !print*, i, j, b%A(i,j)
              if(abs(b%A(i,j)).ne.0.d0) then
               sparsity = sparsity + 1
              endif
           enddo
        enddo
        sparsity = sparsity/nlev**2
        return
    end subroutine sparsity_calc

end module linear_algebra
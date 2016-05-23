module linear_algebra

    use types_and_parameters, only: nlev
    ! LAPACK SOLVER VARIABLES
    integer*4, parameter :: nsquared = 25
    integer*4, parameter :: ml = 1
    integer*4, parameter :: mu = 1
 
    integer*4, parameter :: lda = 2 * ml + mu + 1
 
    real*8, dimension(1:lda,1:nsquared) :: a
    real*8, dimension(1:nsquared)       :: b
    integer*4                           :: info
    integer*4, dimension(1:nsquared)    :: ipiv
    integer*4, parameter                :: bandwidth = ml + mu + 1
    integer*4                           :: sparsity

    contains
    
    subroutine sparsity_calc(b, sparsity)
        use types_and_parameters, only: reaction_matrix
        type(reaction_matrix) :: b
        integer               :: sparsity       
        sparsity=0
        print*, nlev
        do i = 1, nlev
           do j = 1, nlev
              print*, i, j
              if(b%A(i,j).ne.0.d0) sparsity = sparsity + 1
           enddo
        enddo
        sparsity = int(sparsity/nlev**2)
        return
    end subroutine sparsity_calc

end module linear_algebra
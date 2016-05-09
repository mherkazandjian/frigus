module matrix_construction

    ! this module evaluates the terms that go in the matrix that contains all the collisional 
    ! and radiative contribution for the formation/destruction of each level 
    use types_and_parameters, only: nlev,                                 &
                                    radiative_coeffs, collisional_coeffs, &
                                    reaction_matrix, Trad

    contains

    subroutine matrix_builder(a21, b21, r21, b12, r12, rr, coll_rad_matrix)

               type(radiative_coeffs) :: a21, b21, b12, r21, r12
               type(collisional_coeffs) :: rr
               type(reaction_matrix) :: coll_rad_matrix
               
               print*, 'nlev', nlev


    end subroutine matrix_builder

end module matrix_construction
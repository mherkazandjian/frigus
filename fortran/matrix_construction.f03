module matrix_construction

    ! this module evaluates the terms that go in the matrix that contains all the collisional 
    ! and radiative contribution for the formation/destruction of each level 
    use types_and_parameters, only: nlev, energy_lev,                     &
                                    radiative_coeffs, collisional_coeffs, &
                                    reaction_matrix, population,          &
                                    Trad,                                 &
                                    ini, fin, it, ntemp

    contains

    subroutine matrix_builder(energy, a21, b21, r21, b12, r12, rr, coll_rad_matrix)

               type(energy_lev)       :: energy
               type(radiative_coeffs) :: a21, b21, b12, r21, r12
               type(collisional_coeffs) :: rr
               type(reaction_matrix)    :: coll_rad_matrix
               type(population)         :: x
               
                coll_rad_matrix%A = 0.d0
                print*, 'nlev', nlev
                do it = 1, ntemp
                     do fin = 1, nlev
                         do ini = 1, nlev
                            if((energy%ene(ini)-energy%ene(fin)).lt.0.d0) then 
                                                                  ! according to the ordering in the  
                                                                  ! energy file downwards => E2 - E1 < 0                           
                               coll_rad_matrix%A(ini, fin) = coll_rad_matrix%A(ini, fin)  +    &
                                                             r21%M(ini,fin) * x%pop(ini)  +    &
                                                             r12%M(ini, fin) * x%pop(ini) +    &
                                                             nb * rr%matrix(ini, fin, itemp) * x%pop(ini)
                            else
                            
                            endif
                         enddo
                     enddo
                enddo

    end subroutine matrix_builder

!     subroutine initialize_level_population(x)
!         type(population) :: x
!         
!         do i = 1, nlev-1
!            x%pop(i) = 1.d-3
!         enddo
!         
!         x%pop(nlev) = 1.d0 - sum(x%pop, dim = 3)
!         
!     end subroutine initialize_level_population
    
end module matrix_construction
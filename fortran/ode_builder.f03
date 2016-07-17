module ode_contruction

    ! this module evaluates the rhs terms for the ODEs
    
    use types_and_parameters, only: nlev, energy_lev,                     &
                                    radiative_coeffs, collisional_coeffs, &
                                    reaction_matrix, population,          &
                                    Trad,                                 &
                                    it, ntemp

    contains

    subroutine ode_builder(energy, a21, b21, r21, b12, r12, rr, ydot)

               type(energy_lev)       :: energy
               type(radiative_coeffs) :: a21, b21, b12, r21, r12
               type(collisional_coeffs) :: rr
               type(reaction_matrix)    :: coll_rad_matrix
               type(population)         :: x
               
               real*8, dimension(1:nlev) :: ydot
               integer :: i, j
               
               ydot = 0.d0
               !print*, 'nlev', nlev
               do it = 1, ntemp
                  do j = 1, nlev
                     do i = 1, nlev
                        if((energy%ene(ini)-energy%ene(fin)).lt.0.d0) then  ! according to the ordering in the  
                                                                            ! energy file downwards => E2 - E1 < 0                           
                               ydot(i) = ydot(i) + 
                            else
                            
                            endif
                         enddo
                     enddo
                enddo

    end subroutine ode_builder

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
    
end module ode_construction
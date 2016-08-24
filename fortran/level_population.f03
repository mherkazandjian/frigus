module level_population

    use types_and_parameters, only: nlev, energy_lev, vmax, jmax,         &
                                    ntrans, ntemp, i,                      &
                                    radiative_coeffs, collisional_coeffs, &
                                    reaction_matrix,                      &
                                    population,                           &
                                    Trad, nb, ini, fin, it,               &
                                    nlev_lique, vi, vf, ji, jf
                                    
    use energy_levels, only: reading_data_energies
    use radiation,     only: reading_data_radiative,                      &
                             radiative_downwards, radiative_upwards
    use collisions,    only: reading_data_collisions

    use linear_algebra, only: sparsity_calc,                              &
                              ndim, info, lda, ldb, nrhs, ipiv
    
    use matrix_construction, only: matrix_builder, initialize_level_population
    

    contains
    
    subroutine lev_pop(energy, a21, b21, r21, b12, jnu, r12, rr, rr21, rr12, id_temp, coll_rad_matrix, x)
    
        type(energy_lev) :: energy
        type(radiative_coeffs) :: a21, b21, b12, jnu, r21, r12, rad
        type(collisional_coeffs) :: rr, rr21, rr12
        type(reaction_matrix)  :: coll_rad_matrix
        type(population) :: x, y

        call reading_data_energies(energy)

        call reading_data_radiative(energy, a21)
    
        call reading_data_collisions(energy, rr, rr21, rr12)
        
        call radiative_downwards(energy, Trad, a21, b21, r21)
            
        call radiative_upwards(energy, Trad, a21, b21, b12, jnu, r12)

!         ! building the total radiative matrix, including both stimulated and spontaneous
!         ! transitions; according to the convention adopted:
!         !     downward transitions -> upper triangular matrix
!         !     upward transitions   -> lower triangular matrix

         rad%M_lique = r21%M_lique + r12%M_lique


         call multiplication_by_nc(rr, rr21, rr12, nb(1))
         

 
!          
          call matrix_builder(rad, rr, id_temp, coll_rad_matrix)
!  
!          call initialize_level_population(y)
 
!         x%pop = y%pop
 
!         print*, 'before', x%pop
    
!         call dgesv(ndim, nrhs, coll_rad_matrix, lda, ipiv, x, ldb, info)
 
!         do i = 1, nlev_lique
!             write(6,'(3(i3), e14.7)') i, energy%vl(i), energy%jl(i), x%pop(i)
!         enddo
     
        return
    end subroutine lev_pop

    subroutine multiplication_by_nc(rr, rr21, rr12, nc) 
        real*8 nc
        type(collisional_coeffs) :: rr, rr21, rr12
        type(collisional_coeffs) :: rr_new, rr21_new, rr12_new
        print*, 'nc = ', nc
        rr%matrix_lique = rr%matrix_lique * nc
        rr21%matrix_lique = rr21%matrix_lique * nc
        rr12%matrix_lique = rr12%matrix_lique * nc
        return
    end subroutine multiplication_by_nc

end module level_population
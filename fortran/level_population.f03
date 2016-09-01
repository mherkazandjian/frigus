module level_population

    use types_and_parameters, only: nlev, energy_lev, vmax, jmax,         &
                                    ntrans, ntemp, i,                      &
                                    radiative_coeffs, collisional_coeffs, &
                                    reaction_matrix,                      &
                                    population,                           &
                                    Trad, nc, ini, fin, it,               &
                                    nlev_lique, vi, vf, ji, jf

    use linear_algebra, only: sparsity_calc,                              &
                              ndim, info, lda, ldb, nrhs, ipiv
    
    use matrix_construction, only: matrix_builder, initialize_level_population
    
    contains
    
    subroutine lev_pop(energy, a21, b21, r21, b12, jnu, r12, rr, rr21, rr12, id_temp, nc, coll_rad_matrix, x)
         type(energy_lev) :: energy
         type(radiative_coeffs) :: a21, b21, b12, jnu, r21, r12, rad
         type(collisional_coeffs) :: rr, rr21, rr12
         type(reaction_matrix)  :: coll_rad_matrix
         type(population) :: x
         real*8 :: nc
         

         rad%M_lique = r21%M_lique + r12%M_lique

         call matrix_builder(rad, rr, id_temp, coll_rad_matrix)
         
         call solve_steady_state(energy, coll_rad_matrix, x)

        return
    end subroutine lev_pop


    
    subroutine solve_steady_state(energy, coll_rad_matrix, x)
        type(energy_lev) :: energy
        type(reaction_matrix)  :: coll_rad_matrix
        type(population) :: x, y

          call initialize_level_population(y)

          ! normalizing the sum of the fractional abundances to 1
          coll_rad_matrix%M(nlev_lique, 1:nlev_lique) = 1.d0
 
          x%pop = y%pop
 
          !print*, 'before', x%pop
    
          call dgesv(ndim, nrhs, coll_rad_matrix%M, lda, ipiv, x, ldb, info)
 
          !do i = 1, nlev_lique
          !    write(6,'(3(i3), 2(ES23.15))') i, energy%vl(i), energy%jl(i), y%pop(i), x%pop(i)
          !enddo

        return
    end subroutine solve_steady_state

end module level_population
module level_population

    use types_and_parameters, only: nlev, energy_lev, vmax, jmax,         &
                                    ntrans, ntemp, i,                      &
                                    radiative_coeffs, collisional_coeffs, &
                                    reaction_matrix,                      &
                                    population,                           &
                                    Trad, nc, ini, fin, it,               &
                                    nlev_lique, vi, vf, ji, jf,           &
                                    ndensity, id_temp,                    &
                                    norm_first_row

    use linear_algebra, only: sparsity_calc,                              &
                              ndim, info, lda, ldb, nrhs, ipiv
    
    use matrix_construction, only: matrix_builder, initialize_level_population
    
    contains
    
    subroutine lev_pop(norm_first_row, energy, a21, b21, r21, b12, jnu, r12, rr, rr21, rr12, id_temp, nc, coll_rad_matrix, x)
         type(energy_lev) :: energy
         type(radiative_coeffs) :: a21, b21, b12, jnu, r21, r12, rad
         type(collisional_coeffs) :: rr, rr21, rr12
         type(reaction_matrix)  :: coll_rad_matrix
         type(population) :: x
         real*8, dimension(1:ndensity) :: nc
         

         rad%M_lique = r21%M_lique + r12%M_lique

         call matrix_builder(rad, rr, id_temp, coll_rad_matrix)
         
         call solve_steady_state(norm_first_row, energy, coll_rad_matrix, x)

        return
    end subroutine lev_pop


    
    subroutine solve_steady_state(norm_first_row, energy, coll_rad_matrix, x)
        type(energy_lev) :: energy
        type(reaction_matrix)  :: coll_rad_matrix, coll_rad_matrix_mher
        type(population) :: x, y

          call initialize_level_population(norm_first_row, y)

          
          if(norm_first_row.eq.1) then
          ! normalizing the sum of the fractional abundances to 1
            coll_rad_matrix%M(1, 1:nlev_lique) = 1.d0
          else
            coll_rad_matrix%M(nlev_lique, 1:nlev_lique) = 1.d0
          endif
          
          x%pop = y%pop


          call dgesv(ndim, nrhs, coll_rad_matrix%M, lda, ipiv, x, ldb, info)
          !dgesvx

           do i = 1, nlev_lique
               write(6,'(3(i3), 2(ES26.16E3))') i, energy%vl(i), energy%jl(i), y%pop(i), x%pop(i)
           enddo
           write(6,'(a26,es24.18)') 'sum fractional abundances:', sum(x%pop)
        return
    end subroutine solve_steady_state

end module level_population
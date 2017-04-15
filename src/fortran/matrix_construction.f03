module matrix_construction

    ! this module evaluates the terms that go in the matrix that contains all the collisional 
    ! and radiative contribution for the formation/destruction of each level 
    use types_and_parameters, only: nlev_lique,                           &
                                    radiative_coeffs, collisional_coeffs, &
                                    reaction_matrix, population,          &
                                    Trad, nc,                             &
                                    row, col, fin, ntemp, xm,             &
                                    id_temp, id_temp_test

    use testing_data, only: tests, writing_files


    contains

    subroutine matrix_builder(rad, rr, id_temp, coll_rad_matrix)

                 type(radiative_coeffs)   :: rad
                 type(collisional_coeffs) :: rr, rr21, rr12
                 type(reaction_matrix)    :: coll_rad_matrix
                 real*8, dimension(1:nlev_lique) :: temp

                 coll_rad_matrix%M = 0.d0 ! full matrix
                 coll_rad_matrix%O = 0.d0 ! off-diagonal terms
                 coll_rad_matrix%D = 0.d0 ! diagonal terms

                 do row = 1, nlev_lique
                     do col = 1, nlev_lique
                         if(row.ne.col) then 
                             coll_rad_matrix%O(row, col) = coll_rad_matrix%O(row, col)      &
                                                         + rad%M_lique(col, row)            &
                                                         + rr%matrix_lique(col, row, id_temp) 
                       endif
                     enddo
                 enddo

                !print*, shape(sum(coll_rad_matrix%O, dim = 1))
 
                 do row = 1, nlev_lique
                     do col = 1, nlev_lique
                         if(row.eq.col) then
                              temp = sum(coll_rad_matrix%O, dim = 1)
                              coll_rad_matrix%D(row, col) = -temp(row)
                         endif
                     enddo
                 enddo
                                  
                 coll_rad_matrix%M = coll_rad_matrix%O + coll_rad_matrix%D
                 open(40, file='carla_matrix_M', status = 'unknown')
                 if(id_temp.eq.1) then
                    do i = 1, nlev_lique
!                        do j = 1, nlev_lique
!                            write(40, '(3(i3), 1x, 59(es24.10,1x))')     &
                            write(40, '(59(es24.16e3, 1x))')     &
                              (coll_rad_matrix%M(i,j), j=1, nlev_lique)
!                        enddo
                    enddo
                endif


    end subroutine matrix_builder

    
    
    
     subroutine initialize_level_population(x)
         type(population) :: x
         
          do i = 1, nlev_lique-1 !i = 2, nlev_lique !
             x%pop(i) = 0.d0
          enddo

          x%pop(nlev_lique) = 1.d0
          !x%pop(1) = 1.d0          
         
     end subroutine initialize_level_population
    
end module matrix_construction
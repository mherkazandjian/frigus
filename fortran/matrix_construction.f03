module matrix_construction

    ! this module evaluates the terms that go in the matrix that contains all the collisional 
    ! and radiative contribution for the formation/destruction of each level 
    use types_and_parameters, only: nlev_lique,                           &
                                    radiative_coeffs, collisional_coeffs, &
                                    reaction_matrix, population,          &
                                    Trad, nb,                             &
                                    row, col, fin, it, ntemp, xm

    contains

    subroutine matrix_builder(rad, rr, id_t, coll_rad_matrix)

                 type(radiative_coeffs)   :: rad
                 type(collisional_coeffs) :: rr
                 type(reaction_matrix)    :: coll_rad_matrix
                 integer :: it
               
!              do it = 1, ntemp
                 coll_rad_matrix%M = 0.d0 ! full matrix
                 coll_rad_matrix%O = 0.d0 ! off-diagonal terms
                 coll_rad_matrix%D = 0.d0 ! diagonal terms
                 do row = 1, nlev_lique
                     do col = 1, nlev_lique
                         if(row.ne.col) then 
                             coll_rad_matrix%O(row, col, id_t) = coll_rad_matrix%O(row, col, id_t)   &
                                                             + rad%M_lique(row, col)             &
                                                             + rr%matrix_lique(row, col, id_t) 
                       endif
                     enddo
                 enddo
                
                

!                                   else
!                                  coll_rad_matrix%A(row,col) = coll_rad_matrix%A(row, col)   &
!                                                             + 0.d0
!                                 endif
!                              enddo
!                         elseif(row.eq.nlev) then
!                              do fin = 1, nlev
!                                 coll_rad_matrix%A(row,fin) = 1.d0    ! normalization equation
!                              enddo                        
!                         else   ! the upwards and downwards transitions are already 
!                                ! implemented in the proper way in the rr and rad matrix 
!                                ! according to the indexes
!                         !write(6,'(a2, i3, a1, i3, a9, i3, a1, i3, a5, i3, a1, i3, a1, i3, a4)')  &
!                         !'M(', row, ',', col, ') = +r12(', row, ',', col, ') +c(',row,',', col, ',', it, ')*nb'
!                          coll_rad_matrix%A(row,col) = rad%M(row, col) + rr%matrix(col, row, it)*nb(1)
!                         ! CHECK ORDER IN THE INDEXES in the previous line
! !                        elseif(row.lt.col) then  ! upper triangular matrix -> downward transitions
! !                                 coll_rad_matrix%A(row,col) = rad%M(row, col) + rr%matrix(row, col, it)*nb(1)
! !                        elseif(row.gt.col) then
! !                                 coll_rad_matrix%A(row,col) = rad%M(row, col) + rr%matrix(row, col, it)*nb(1)
!                         endif
!                      enddo
!                   enddo
!                 !enddo
    end subroutine matrix_builder

     subroutine initialize_level_population(x)
         type(population) :: x
         
         do i = 1, nlev-1
            x%pop(i) = 0.d0 !1.d-3
         enddo
         
         x%pop(nlev) = 1.d0 !1.d0 - sum(x%pop)
         
     end subroutine initialize_level_population
    
end module matrix_construction
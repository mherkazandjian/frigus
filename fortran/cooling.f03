program cooling
     use types_and_parameters, only: id_temp, energy_lev,                   &
                                     radiative_coeffs, collisional_coeffs,  &
                                     reaction_matrix, population,           &
                                     ntemp, nlev_lique, Trad

     use read_data, only: get_data                                     
                                     
     use level_population, only: lev_pop
         
     type(energy_lev) :: energy
     type(radiative_coeffs) :: a21, b21, r21, b12, r12, jnu
     type(collisional_coeffs) :: rr, rr21, rr12
     type(reaction_matrix)  :: coll_rad_matrix
     type(population) :: x
     real*8, dimension(1:ntemp) :: cooling_rate
    
     call get_data(Trad, id_temp, energy, a21, b21, r21, b12, jnu, r12, rr, rr21, rr12)
!     do itemp = 1, ntemp 
     itemp = 30
     call lev_pop(energy, a21, b21, r21, b12, jnu, r12, rr, rr21, rr12, itemp, coll_rad_matrix, x)     
        cooling_rate(itemp) = 0.d0     
        do i = 1, nlev_lique
            do j = 1, nlev_lique
                if(i.gt.j) then
                    cooling_rate(itemp) = cooling_rate(itemp) + a21%M_lique(i, j)                 &
                                        * abs(energy%ene_lique(i)-energy%ene_lique(j))    &
                                        * x%pop(i)
                endif
            enddo
        enddo
        write(6,'(ES23.15)') rr%temp(itemp), cooling_rate(itemp)
     !enddo
     
     
end program cooling
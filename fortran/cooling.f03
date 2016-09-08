program cooling
     use types_and_parameters, only: energy_lev,                            &
                                     radiative_coeffs, collisional_coeffs,  &
                                     reaction_matrix, population,           &
                                     ntemp, nlev_lique, Trad, nc,           &
                                     id_temp, id_temp_test, ndensity

     use read_data, only: get_data                                     
                                     
     use level_population, only: lev_pop
     
     use testing_data, only: tests, writing_files

         
     type(energy_lev) :: energy
     type(radiative_coeffs) :: a21, b21, r21, b12, r12, jnu
     type(collisional_coeffs) :: rr, rr21, rr12
     type(reaction_matrix)  :: coll_rad_matrix
     type(population) :: x
     real*8, dimension(1:ntemp) :: cooling_rate, glover
     character*7 :: char1
     character*19 :: filename

     call get_data(Trad, id_temp_test, energy, a21, b21, r21, b12, jnu, r12, rr, rr21, rr12)

     do idensity = 1, ndensity
        write(*,'(a5, ES23.15)') 'nc = ', nc(idensity)
        rr%matrix_lique = rr%matrix_lique * nc(idensity)
        rr21%matrix_lique = rr21%matrix_lique * nc(idensity)
        rr12%matrix_lique = rr12%matrix_lique * nc(idensity)

        do id_temp = 1, ntemp 
            call lev_pop(energy, a21, b21, r21, b12, jnu, r12, rr, rr21, rr12, id_temp, nc(idensity), coll_rad_matrix, x)    
            write(6, '(a17, ES23.15)') 'gas temperature: ', rr%temp(id_temp)
            cooling_rate(id_temp) = 0.d0
            glover(id_temp) = 0.d0
            do i = 1, nlev_lique
                do j = 1, nlev_lique
                    if(i.gt.j) then
                        cooling_rate(id_temp) = cooling_rate(id_temp) + a21%M_lique(i, j)                 &
                                            * abs(energy%ene_lique(i)-energy%ene_lique(j))    &
                                            * x%pop(i)
                    endif
                enddo
            enddo

            if(rr%temp(id_temp).ge.100.d0.and.rr%temp(id_temp).le.1000.d0) then
                    glover(id_temp) = 10**(-24.311209                                  &
                                        +3.5692468*log10(rr%temp(id_temp)/1000.)       &
                                        -11.332860*(log10(rr%temp(id_temp)/1000.))**2  &
                                        -27.850082*(log10(rr%temp(id_temp)/1000.))**3  &
                                        -21.328264*(log10(rr%temp(id_temp)/1000.))**4  &
                                        -4.2519023*(log10(rr%temp(id_temp)/1000.))**5)
            elseif(1000 < rr%temp(id_temp).and.rr%temp(id_temp)<=6000) then
                    glover(id_temp) = 10**(-24.311209                                  &
                                    +4.6450521*log10(rr%temp(id_temp)/1000.)            &
                                    -3.7209846*log10((rr%temp(id_temp)/1000.))**2       &  
                                    +5.9369081*log10((rr%temp(id_temp)/1000.))**3      &
                                    -5.5108047*log10((rr%temp(id_temp)/1000.))**4      &
                                    +1.5538288*log10((rr%temp(id_temp)/1000.))**5)
            endif
            glover(id_temp) = glover(id_temp)*1.d-7
            !call tests(energy, rr, rr21, rr12, a21, b21, r21, b12, jnu, r12, coll_rad_matrix, id_temp)
            !call writing_files(a21, b21, b12, jnu, id_temp, rr, coll_rad_matrix)
             write(char1, '(ES7.1)') nc(idensity)
             filename = 'cooling_nc=' // char1
             open(40, file = filename, status = 'unknown')
             write(40,'(4(ES23.15))') Trad, rr%temp(id_temp), cooling_rate(id_temp), glover(id_temp)            
        enddo ! loop on the kinetic temperatures
        rr%matrix_lique = rr%matrix_lique / nc(idensity)
        rr21%matrix_lique = rr21%matrix_lique / nc(idensity)
        rr12%matrix_lique = rr12%matrix_lique / nc(idensity)
     enddo ! loop on the density

end program cooling

     !subroutine multiplication_by_nc(rr, rr21, rr12, nc) 
     !   real*8 :: nc
     !   type(collisional_coeffs) :: rr, rr21, rr12
     !   write(*,'(a5, ES23.15)') 'nc = ', nc
     !   rr%matrix_lique = rr%matrix_lique * nc
     !   rr21%matrix_lique = rr21%matrix_lique * nc
     !   rr12%matrix_lique = rr12%matrix_lique * nc
     !   return
     !end subroutine multiplication_by_nc
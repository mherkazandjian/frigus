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
                        cooling_rate(id_temp) = cooling_rate(id_temp) + a21%M_lique(i, j)     &
                                            * abs(energy%ene_lique(i)-energy%ene_lique(j))    &
                                            * x%pop(i)
!                         cooling_rate(id_temp) = cooling_rate(id_temp) +                       & !lipovka expression
!                                             (a21%M_lique(i, j)                                &
!                                             * abs(energy%ene_lique(i)-energy%ene_lique(j))    &
!                                             * x%pop(i)-                                       &
!                                             a21%M_lique(j, i)                                 &
!                                             * abs(energy%ene_lique(j)-energy%ene_lique(i))    &
!                                             * x%pop(j))
                    endif
                enddo
            enddo

            if(rr%temp(id_temp).ge.100.d0.and.rr%temp(id_temp).le.1000.d0) then
                    glover(id_temp) = 0.25*10**(-24.216387                                    & ! para-H2
                                         +3.3237480*log10(rr%temp(id_temp)/1000.)        &
                                         -11.642384*(log10(rr%temp(id_temp)/1000.))**2   &
                                         -35.553366*(log10(rr%temp(id_temp)/1000.))**3   &
                                         -35.105689*(log10(rr%temp(id_temp)/1000.))**4   &
                                         -10.922078*(log10(rr%temp(id_temp)/1000.))**5)+ &
                                      0.75*10**(-24.330855                                    & ! ortho-H2
                                         +4.4404496*log10(rr%temp(id_temp)/1000.)        &
                                         -4.0460989*(log10(rr%temp(id_temp)/1000.))**2   &
                                         -1.1390725*(log10(rr%temp(id_temp)/1000.))**3   &
                                         +9.8094223*(log10(rr%temp(id_temp)/1000.))**4   &
                                         +8.6273872*(log10(rr%temp(id_temp)/1000.))**5)
             elseif(1000 < rr%temp(id_temp).and.rr%temp(id_temp)<=6000) then
                     glover(id_temp) = 0.25*10**(-24.216387                                &  ! para-H2
                                     +4.2046488*log10(rr%temp(id_temp)/1000.)         &
                                     -1.3155285*log10((rr%temp(id_temp)/1000.))**2    &  
                                     -1.6552763*log10((rr%temp(id_temp)/1000.))**3    &
                                     +4.1780102*log10((rr%temp(id_temp)/1000.))**4    &
                                     -0.56949697*log10((rr%temp(id_temp)/1000.))**5   &
                                     -3.3824407*log10((rr%temp(id_temp)/1000.))**6    &
                                     +1.0904027*log10((rr%temp(id_temp)/1000.))**7) + &
                                     0.75*10**(-24.329086                                  &  ! ortho-H2                                     
                                     +4.6105087*log10(rr%temp(id_temp)/1000.)         &
                                     -3.9505350*log10((rr%temp(id_temp)/1000.))**2    &  
                                     +12.363818*log10((rr%temp(id_temp)/1000.))**3    &
                                     -32.403165*log10((rr%temp(id_temp)/1000.))**4    &
                                     +48.853562*log10((rr%temp(id_temp)/1000.))**5    &
                                     -38.542008*log10((rr%temp(id_temp)/1000.))**6    &
                                     +12.066770*log10((rr%temp(id_temp)/1000.))**7)
             endif            
            
            
!            if(rr%temp(id_temp).ge.100.d0.and.rr%temp(id_temp).le.1000.d0) then
!                    glover(id_temp) = 10**(-24.311209                                  &
!                                         +3.5692468*log10(rr%temp(id_temp)/1000.)       &
!                                         -11.332860*(log10(rr%temp(id_temp)/1000.))**2  &
!                                         -27.850082*(log10(rr%temp(id_temp)/1000.))**3  &
!                                         -21.328264*(log10(rr%temp(id_temp)/1000.))**4  &
!                                         -4.2519023*(log10(rr%temp(id_temp)/1000.))**5)
!             elseif(1000 < rr%temp(id_temp).and.rr%temp(id_temp)<=6000) then
!                     glover(id_temp) = 10**(-24.311209                                  &
!                                     +4.6450521*log10(rr%temp(id_temp)/1000.)            &
!                                     -3.7209846*log10((rr%temp(id_temp)/1000.))**2       &  
!                                     +5.9369081*log10((rr%temp(id_temp)/1000.))**3      &
!                                     -5.5108047*log10((rr%temp(id_temp)/1000.))**4      &
!                                     +1.5538288*log10((rr%temp(id_temp)/1000.))**5)
!             endif
            glover(id_temp) = glover(id_temp)*1.d-7
            !call tests(energy, rr, rr21, rr12, a21, b21, r21, b12, jnu, r12, coll_rad_matrix, id_temp)
            !call writing_files(a21, b21, b12, jnu, id_temp, rr, coll_rad_matrix)
             write(char1, '(ES7.1)') nc(idensity)
             filename = 'cooling_nc=' // char1
             open(40, file = filename, status = 'unknown')
             open(41, file = 'to_be_fitted.dat', status = 'unknown')             
             write(40,'(4(ES23.15))') Trad, rr%temp(id_temp), cooling_rate(id_temp), glover(id_temp)
             write(41,'(4(ES23.15))') Trad, nc(idensity), rr%temp(id_temp), cooling_rate(id_temp)
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
program cooling
     use types_and_parameters, only: energy_lev,                            &
                                     radiative_coeffs, collisional_coeffs,  &
                                     reaction_matrix, population,           &
                                     ntemp, nlev_lique, Trad,               &
                                     id_temp, id_temp_test

     use read_data, only: get_data                                     
                                     
     use level_population, only: lev_pop
         
     type(energy_lev) :: energy
     type(radiative_coeffs) :: a21, b21, r21, b12, r12, jnu
     type(collisional_coeffs) :: rr, rr21, rr12
     type(reaction_matrix)  :: coll_rad_matrix
     type(population) :: x
     real*8, dimension(1:ntemp) :: cooling_rate, glover

     open(40, file = 'cooling.txt', status = 'unknown')
     
     call get_data(Trad, id_temp_test, energy, a21, b21, r21, b12, jnu, r12, rr, rr21, rr12)

!      do itemp = 1, ntemp 
!          call lev_pop(energy, a21, b21, r21, b12, jnu, r12, rr, rr21, rr12, itemp, coll_rad_matrix, x)     
!          write(6, '(a17, ES23.15)') 'gas temperature: ', rr%temp(itemp)
!          cooling_rate(itemp) = 0.d0
!          glover(itemp) = 0.d0
!          do i = 1, nlev_lique
!              do j = 1, nlev_lique
!                 if(i.gt.j) then
!                     cooling_rate(itemp) = cooling_rate(itemp) + a21%M_lique(i, j)                 &
!                                         * abs(energy%ene_lique(i)-energy%ene_lique(j))    &
!                                         * x%pop(i)
!                 endif
!              enddo
!          enddo
!          
!          if(rr%temp(itemp).ge.100.d0.and.rr%temp(itemp).le.1000.d0) then
!                 glover(itemp) = 10**(-24.311209                                  &
!                                     +3.5692468*log10(rr%temp(itemp)/1000.)       &
!                                     -11.332860*(log10(rr%temp(itemp)/1000.))**2  &
!                                     -27.850082*(log10(rr%temp(itemp)/1000.))**3  &
!                                     -21.328264*(log10(rr%temp(itemp)/1000.))**4  &
!                                     -4.2519023*(log10(rr%temp(itemp)/1000.))**5)
!          elseif(1000 < rr%temp(itemp).and.rr%temp(itemp)<=6000) then
!                 glover(itemp) = 10**(-24.311209                                  &
!                                +4.6450521*log10(rr%temp(itemp)/1000.)            &
!                                -3.7209846*log10((rr%temp(itemp)/1000.))**2       &  
!                                 +5.9369081*log10((rr%temp(itemp)/1000.))**3      &
!                                 -5.5108047*log10((rr%temp(itemp)/1000.))**4      &
!                                 +1.5538288*log10((rr%temp(itemp)/1000.))**5)
!         endif
!         glover(itemp) = glover(itemp)*1.d-7
!         write(40,'(4(ES23.15))') Trad, rr%temp(itemp), cooling_rate(itemp), glover(itemp)
!     enddo
     
     
end program cooling
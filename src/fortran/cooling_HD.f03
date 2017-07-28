program cooling
     use types_and_parameters, only: energy_lev,                                     &
                                     radiative_coeffs, collisional_coeffs,           &
                                     reaction_matrix, population,                    &
                                     ntemp, nlev_lique, Trad, nc,                    &
                                     id_temp, id_temp_test, ndensity

     use read_data, only: get_data

     use level_population, only: lev_pop

     use testing_data, only: tests, writing_files

     type(energy_lev) :: energy
     type(radiative_coeffs) :: a21, b21, r21, b12, r12, jnu
     type(collisional_coeffs) :: rr, rr21, rr12
     type(reaction_matrix)  :: coll_rad_matrix
     type(population) :: x
     real*8 :: ln_hd, lt_kin
     real*8, dimension(1:ntemp) :: cooling_rate, glover, lipovka
     character*7 :: char1
     character*19 :: filename

     call get_data(0, 0, 1, Trad,&
     id_temp_test, energy, a21, b21, r21, b12, jnu, r12, rr, rr21, rr12)
     


     do idensity = 1, ndensity
        write(*,'(a5, ES23.15)') 'nc = ', nc(idensity)
        rr%matrix_lique = rr%matrix_lique * nc(idensity)
        rr21%matrix_lique = rr21%matrix_lique * nc(idensity)
        rr12%matrix_lique = rr12%matrix_lique * nc(idensity)

        do id_temp = 1, ntemp 
            call lev_pop(norm_first_row, energy, a21, b21, r21, b12, jnu,&
            r12, rr, rr21, rr12, id_temp, nc(idensity), coll_rad_matrix, x)
            write(6, '(a17, ES23.15)') 'gas temperature: ', rr%temp(id_temp)
            cooling_rate(id_temp) = 0.d0
            glover(id_temp) = 0.d0
            lipovka(id_temp) = 0.d0
            do i = 1, nlev_lique
                do j = 1, nlev_lique
                    if(i.gt.j) then
                         cooling_rate(id_temp) = cooling_rate(id_temp) + a21%M_lique(i, j)                 &
                                             * abs(energy%ene_lique(i)-energy%ene_lique(j))    &
                                             * x%pop(i)
                    endif
                enddo
            enddo


            ln_hd = dlog10(nc(idensity)*1.d-6) ! the fit by lipovka is valid for densities in cm-3
            lt_kin = dlog10(rr%temp(id_temp))

            lipovka(id_temp) =                                          &
                      10**(-42.57688                                    &
                       + 0.92433d0 *                   (ln_hd ** 1)     &
                       + 0.54962d0 *                   (ln_hd ** 2)     &
                       - 0.07676d0 *                   (ln_hd ** 3)     &
                       + 0.00275d0 *                   (ln_hd ** 4)     &
                       + 21.93385d0 * (lt_kin ** 1)                     &
                       + 0.77952d0 *  (lt_kin ** 1)  * (ln_hd ** 1)     &
                       - 1.06447d0 *  (lt_kin ** 1)  * (ln_hd ** 2)     &
                       + 0.11864d0 *  (lt_kin ** 1)  * (ln_hd ** 3)     &
                       - 0.00366d0 *  (lt_kin ** 1)  * (ln_hd ** 4)     &
                       - 10.19097d0 * (lt_kin ** 2)                     &
                       - 0.54263d0 *  (lt_kin ** 2)  * (ln_hd ** 1)     &
                       + 0.62343d0 *  (lt_kin ** 2)  * (ln_hd ** 2)     &
                       - 0.07366d0 *  (lt_kin ** 2)  * (ln_hd ** 3)     &
                       + 0.002514d0 * (lt_kin ** 2)  * (ln_hd ** 4)     &
                       + 2.19906d0 *  (lt_kin ** 3)                     &
                       + 0.11711d0 *  (lt_kin ** 3)  * (ln_hd ** 1)     &
                       - 0.13768d0 *  (lt_kin ** 3)  * (ln_hd ** 2)     &
                       + 0.01759d0 *  (lt_kin ** 3)  * (ln_hd ** 3)     &
                    - 0.000666317d0 * (lt_kin ** 3)  * (ln_hd ** 4)     &
                       - 0.17334d0 *  (lt_kin ** 4)                     &
                       - 0.00835d0 *  (lt_kin ** 4)  * (ln_hd ** 1)     &
                       + 0.0106d0  *  (lt_kin ** 4)  * (ln_hd ** 2)     &
                       - 0.001482d0 * (lt_kin ** 4)  * (ln_hd ** 3)     &
                    + 0.000061926d0 * (lt_kin ** 4)  * (ln_hd ** 4))

            lipovka(id_temp) = lipovka(id_temp)*1.d-7

            call tests(energy, rr, rr21, rr12, a21, b21, r21, b12, jnu, r12, coll_rad_matrix, id_temp_test)


             write(char1, '(ES7.1)') nc(idensity)
             filename = 'cooling_nc=' // char1
             open(50, file = filename, status = 'unknown')
             open(51, file = 'to_be_fitted.dat', status = 'unknown')             
             write(50,'(5(ES23.15))') Trad, rr%temp(id_temp), cooling_rate(id_temp), lipovka(id_temp)
             write(51,'(4(ES23.15))') Trad, nc(idensity), rr%temp(id_temp), cooling_rate(id_temp)
        enddo ! loop on the kinetic temperatures
        rr%matrix_lique = rr%matrix_lique / nc(idensity)
        rr21%matrix_lique = rr21%matrix_lique / nc(idensity)
        rr12%matrix_lique = rr12%matrix_lique / nc(idensity)
     enddo ! loop on the density

end program cooling
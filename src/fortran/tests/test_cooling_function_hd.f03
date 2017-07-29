program test_cooling_function_hd
  ! test for the hd cooling function computed with frigus 
  ! in the fortran version (densities 10^6 m-3 up to 10^14 m-3);
  ! for n = 10^6 m-3 the numbers are the same reported in the 
  ! corresponding test in the frigus - python version
     use types_and_parameters, only: energy_lev,                                     &
                                     radiative_coeffs, collisional_coeffs,           &
                                     reaction_matrix, population,                    &
                                     nlev_lique, Trad, nc, ntemp,                    &
                                     id_temp, id_temp_test

     use read_data, only: get_data

     use level_population, only: lev_pop

     use testing_data, only: writing_files

     type(energy_lev) :: energy
     type(radiative_coeffs) :: a21, b21, r21, b12, r12, jnu
     type(collisional_coeffs) :: rr, rr21, rr12
     type(reaction_matrix)  :: coll_rad_matrix
     type(population) :: x
     real*8 :: ln_hd, lt_kin
     real*8, dimension(1:ntemp) :: cooling_rate, lipovka
     character*7 :: char1
     character*19 :: filename
     integer k


    real*8, parameter :: cooling_rate_expected_1e6m3(*) = &
					            [2.239044e-32, &     ! T = 100 K
					             2.861033e-31, &     ! T = 500 K
          					     6.783145e-31, &     ! T = 1000 K
					             1.108309e-30, &     ! T = 1500 K
					             1.570218e-30]       ! T = 2000 K

    real*8, parameter :: cooling_rate_expected_1e7m3(*) = &
					            [2.230320e-31, &     ! T = 100 K
                                                     2.833508e-30, &     ! T = 500 K
                                                     6.700052e-30, &     ! T = 1000 K
					             1.096455e-29, &     ! T = 1500 K
					             1.555448e-29]       ! T = 2000 K


    real*8, parameter :: cooling_rate_expected_1e8m3(*) = &
					            [2.147382e-30, &     ! T = 100 K
  					             2.626509e-29, &     ! T = 500 K
          					     6.189179e-29, &     ! T = 1000 K
					             1.029352e-28, &     ! T = 1500 K
					             1.475257e-28]       ! T = 2000 K

    real*8, parameter :: cooling_rate_expected_1e9m3(*) = &
					            [1.598687e-29, &     ! T = 100 K
  					             2.035301e-28, &     ! T = 500 K
          					     5.124798e-28, &     ! T = 1000 K
					             8.757443e-28, &     ! T = 1500 K
					             1.273190e-27]       ! T = 2000 K


    real*8, parameter :: cooling_rate_expected_1e10m3(*) = &
					            [5.986354e-29, &     ! T = 100 K
  					             1.203061e-27, &     ! T = 500 K
          					     3.385624e-27, &     ! T = 1000 K
					             5.984287e-27, &     ! T = 1500 K
					             8.893297e-27]       ! T = 2000 K

    real*8, parameter :: cooling_rate_expected_1e11m3(*) = &
					            [1.217292e-28, &     ! T = 100 K
  					             4.127840e-27, &     ! T = 500 K
          					     1.378790e-26, &     ! T = 1000 K
					             2.637886e-26, &     ! T = 1500 K
					             4.036089e-26]       ! T = 2000 K

    real*8, parameter :: cooling_rate_expected_1e12m3(*) = &
					            [1.455773e-28, &     ! T = 100 K
  					             7.216173e-27, &     ! T = 500 K
          					     2.882250e-26, &     ! T = 1000 K
					             5.670068e-26, &     ! T = 1500 K
					             8.261490e-26]       ! T = 2000 K

    real*8, parameter :: cooling_rate_expected_1e13m3(*) = &
					            [1.490467e-28, &     ! T = 100 K
  					             8.220277e-27, &     ! T = 500 K
          					     3.475609e-26, &     ! T = 1000 K
					             6.684426e-26, &     ! T = 1500 K
					             9.428664e-26]       ! T = 2000 K

    real*8, parameter :: cooling_rate_expected_1e14m3(*) = &
					            [1.494137e-28, &     ! T = 100 K
  					             8.359668e-27, &     ! T = 500 K
          					     3.558825e-26, &     ! T = 1000 K
					             6.813024e-26, &     ! T = 1500 K
					             9.567276e-26]       ! T = 2000 K

    real*8, parameter :: T_array(*) = [100.0, 500.0, 1000.0, 1500.0, 2000.0]

    integer, parameter :: id_temp_array(*) = [1, 21, 46, 71, 96]
    
    real*8, parameter :: n_array(*) = [1d6, 1d7, 1d8, 1d9, 1d10, 1d11, 1d12, 1d13, 1d14]
    

    call get_data(0, 0, 1, Trad, id_temp_test, energy, a21, b21, r21, b12, jnu, r12, rr, rr21, rr12)


     do idensity = 1, ndensity
        rr%matrix_lique = rr%matrix_lique * n_array(idensity)
        rr21%matrix_lique = rr21%matrix_lique * n_array(idensity)
        rr12%matrix_lique = rr12%matrix_lique * n_array(idensity)
        do k = 1, size(T_array)
            id_temp = id_temp_array(k)
            print*, id_temp, k
            call lev_pop(norm_first_row, energy, a21, b21, r21, b12, jnu,&
            r12, rr, rr21, rr12, id_temp, n_array(idensity), coll_rad_matrix, x)
            write(6, '(a17, ES23.15)') 'gas temperature: ', rr%temp(id_temp)
            cooling_rate(k) = 0.d0
            lipovka(k) = 0.d0
            do i = 1, nlev_lique
                do j = 1, nlev_lique
                    if(i.gt.j) then
                         cooling_rate(k) = cooling_rate(k) + a21%M_lique(i, j)                 &
                                             * abs(energy%ene_lique(i)-energy%ene_lique(j))    &
                                             * x%pop(i)
                    endif
                enddo
            enddo


             if(idensity.eq.1) then
                relative_diff = ((cooling_rate_expected_1e6m3(k)-cooling_rate(k))/(cooling_rate(k)))
                if(relative_diff.gt.1e-6) then 
                    print*, 'test failed for density:', n_array(idensity)
                else
                    write(6,'(a13, ES23.15)') 'test passing!', n_array(idensity)
                endif
            elseif(idensity.eq.2) then
                relative_diff = ((cooling_rate_expected_1e7m3(k)-cooling_rate(k))/(cooling_rate(k)))
                if(relative_diff.gt.1e-6) then 
                    print*, 'test failed for density:', n_array(idensity)
                else
                    write(6,'(a13, ES23.15)') 'test passing!', n_array(idensity)
                endif
            elseif(idensity.eq.3) then
                relative_diff = ((cooling_rate_expected_1e8m3(k)-cooling_rate(k))/(cooling_rate(k)))
                if(relative_diff.gt.1e-6) then 
                    print*, 'test failed for density:', n_array(idensity)
                else
                    write(6,'(a13, ES23.15)') 'test passing!', n_array(idensity)
                endif
            elseif(idensity.eq.4) then
                relative_diff = ((cooling_rate_expected_1e9m3(k)-cooling_rate(k))/(cooling_rate(k)))
                if(relative_diff.gt.1e-6) then 
                    print*, 'test failed for density:', n_array(idensity)
                else
                    write(6,'(a13, ES23.15)') 'test passing!', n_array(idensity)
                endif
            elseif(idensity.eq.5) then
                relative_diff = ((cooling_rate_expected_1e10m3(k)-cooling_rate(k))/(cooling_rate(k)))
                if(relative_diff.gt.1e-6) then 
                    print*, 'test failed for density:', n_array(idensity)
                else
                    write(6,'(a13, ES23.15)') 'test passing!', n_array(idensity)
                endif
            elseif(idensity.eq.6) then
                relative_diff = ((cooling_rate_expected_1e11m3(k)-cooling_rate(k))/(cooling_rate(k)))
                if(relative_diff.gt.1e-6) then 
                    print*, 'test failed for density:', n_array(idensity)
                else
                    write(6,'(a13, ES23.15)') 'test passing!', n_array(idensity)
                endif
            elseif(idensity.eq.7) then
                relative_diff = ((cooling_rate_expected_1e12m3(k)-cooling_rate(k))/(cooling_rate(k)))
                if(relative_diff.gt.1e-6) then 
                    print*, 'test failed for density:', n_array(idensity)
                else
                    write(6,'(a13, ES23.15)') 'test passing!', n_array(idensity)
                endif
            elseif(idensity.eq.8) then
                relative_diff = ((cooling_rate_expected_1e13m3(k)-cooling_rate(k))/(cooling_rate(k)))
                if(relative_diff.gt.1e-6) then 
                    print*, 'test failed for density:', n_array(idensity)
                else
                    write(6,'(a13, ES23.15)') 'test passing!', n_array(idensity)
                endif
            elseif(idensity.eq.9) then
                relative_diff = ((cooling_rate_expected_1e14m3(k)-cooling_rate(k))/(cooling_rate(k)))
                if(relative_diff.gt.1e-6) then 
                    print*, 'test failed for density:', n_array(idensity)
                else
                    write(6,'(a13, ES23.15)') 'test passing!', n_array(idensity)
                endif
            endif
            !print*, 'idensity = ', idensity
        enddo ! loop on the kinetic temperatures
        rr%matrix_lique = rr%matrix_lique / n_array(idensity)
        rr21%matrix_lique = rr21%matrix_lique / n_array(idensity)
        rr12%matrix_lique = rr12%matrix_lique / n_array(idensity)
     enddo ! loop on the density
 
end program test_cooling_function_hd

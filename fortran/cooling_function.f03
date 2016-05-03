program cooling_function

    use types_and_parameters, only: nlev, energy_lev, vmax, jmax, &
                                    ntrans, ntemp,               &
                                    radiative_coeffs, collisional_coeffs
    use energy_levels, only: reading_data_energies
    use radiation, only: reading_data_radiative    
    use collisions, only: reading_data_collisions
    
    type(energy_lev) :: e 
    type(radiative_coeffs) :: a21, b21, b12
    type(collisional_coeffs) :: rr

!    call reading_data_radiative(e, a21)
    
    call reading_data_collisions(e, rr)
    
   
! TEST READING ENERGY LEVELS
!    call reading_data(e)
!    do i = 0, vmax
!       do j = 0, jmax
!           write(6,'(2(i2,2x), e14.8)') i, j, e%en(i,j)
!       enddo
!    enddo

! TEST READING RADIATIVE COEFFICIENTS
!    call reading_data_radiative(a21)
!    print*, a21%ntransrad
!    print*, a21%reading

! TEST COUPLES FOR RADIATIVE TRANSITIONS
!    print*, a21%ntransrad, (jmax+1)*15*15*3
!    do l=1,a21%ntransrad
!        print*, a21%couple1r(l), a21%couple2r(l)
!    enddo

! TEST EINSTEIN COEFFICIENTS AND TRANSITIONS
!    do i = 1, nlev
!       do j = 1, nlev
!           write(6,'(2(i3,2x), e14.8)') i, j, a21%M(i,j)
!       enddo
!    enddo

! TEST COLLISIONAL COEFFICIENTS
! check of the coefficients compared to the read data by François + 
! detailed balance
!    do i=1, ntrans
!       do it=1, ntemp
!          write(6, '(6(i3,2x), 4(e14.5))') rr%couple1c(i),rr%couple2c(i),   &
!                                    rr%vic(i),rr%jic(i),rr%vfc(i),rr%jfc(i),&
!                                    rr%temp(it),                            &
!                     rr%matrix(rr%couple1c(i),rr%couple2c(i),it)*1.d6,      &
!                     rr%matrix(rr%couple2c(i),rr%couple1c(i),it)*1.d6,      &
!                          (rr%matrix(rr%couple1c(i),rr%couple2c(i),it)-     &
!                      rr%matrix(rr%couple2c(i),rr%couple1c(i),it))*1.d6/    &
!                      (rr%matrix(rr%couple1c(i),rr%couple2c(i),it)*1.d6)
!        enddo
!    enddo


end program cooling_function
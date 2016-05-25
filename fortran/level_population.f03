program level_population

    use types_and_parameters, only: nlev, energy_lev, vmax, jmax,         &
                                    ntrans, ntemp, i,                      &
                                    radiative_coeffs, collisional_coeffs, &
                                    reaction_matrix,                      &
                                    population,                           &
                                    Trad, nb, ini, fin
                                    
    use energy_levels, only: reading_data_energies
    use radiation,     only: reading_data_radiative,                      &
                             radiative_downwards, radiative_upwards
    use collisions,    only: reading_data_collisions

    use linear_algebra, only: sparsity_calc,                              &
                              ndim, info, lda, ldb, nrhs, ipiv
    
    use matrix_construction, only: matrix_builder, initialize_level_population
    
    type(energy_lev) :: energy
    type(radiative_coeffs) :: a21, b21, b12, r21, r12, rad
    type(collisional_coeffs) :: rr
    type(reaction_matrix)  :: coll_rad_matrix
    type(population) :: x, y

    call reading_data_energies(energy)
    
    call reading_data_radiative(energy, a21)
    
    call reading_data_collisions(energy, rr)
    
    call radiative_downwards(energy, Trad, a21, b21, r21)
    
    call radiative_upwards(energy, Trad, a21, b12, r12)

    ! building the total radiative matrix, including both stimulated and spontaneous
    ! transitions; according to the convention adopted:
    !     downward transitions -> upper triangular matrix
    !     upward transitions   -> lower triangular matrix

    rad%M = r21%M + r12%M

    call matrix_builder(rad, rr, coll_rad_matrix)

    call initialize_level_population(y)

    x%pop = y%pop

    print*, 'before', x%pop
    
    call dgesv(ndim, nrhs, coll_rad_matrix, lda, ipiv, x, ldb, info)

    do i = 1, nlev
       write(6,'(3(i3), e14.7)') i, energy%vl(i), energy%jl(i), x%pop(i)
    enddo
!    call dgbtrs('N', nlev, 0, 0, 1, coll_rad_matrix, 1, 1, y%pop, nlev, INFO)
    
    !print*, info
    
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
!    do ini = 1, nlev
!        do fin = 1, nlev
!           write(6,'(2(i3, 2x), e14.5)') ini, fin, a21%M(ini, fin)
!        enddo
!    enddo

! TEST COUPLES FOR RADIATIVE TRANSITIONS
!    print*, a21%ntransrad, (jmax+1)*15*15*3
!    do l=1,a21%ntransrad
!        print*, a21%couple1r(l), a21%couple2r(l)
!    enddo

! TEST EINSTEIN COEFFICIENTS AND TRANSITIONS
!  do i=1,a21%ntransrad
!     write(6,'(4(i2,2x),e14.5)') a21%vir(i), a21%jir(i), &
!                                 a21%vfr(i), a21%jfr(i), &
!                   a21%M(a21%couple1r(i),a21%couple2r(i))
!  enddo

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
!    enddo

! TEST RADIATIVE TRANSITIONS COEFFICIENTS
!  do ini = 1, nlev
!     do fin = 1, nlev
!        !if(a21%M(ini, fin).ne.0.d0)                                        &
!        write(6,'(2(i3, 2x), 10(e14.7))') ini, fin,                         &
!                                        energy%ene(ini), energy%ene(fin),   &
!                                        a21%M(ini, fin),                   &
!                                        b21%M(ini, fin),                   &
!                                        r21%M(ini, fin),                   &
!                                        b12%M(ini, fin),                   &
!                                        r12%M(ini, fin),                   &
!                                        ! to have them into the same line although for the reverse transition:
!                                        b12%M(fin, ini),                   &
!                                        r12%M(fin, ini),                   &
!                                        rad%M(ini, fin)
!     enddo
!  enddo

! TEST MATRIX LINEAR SYSTEM
!   do ini = 1, nlev
!        write(6,'(a11, i3, e14.7)') 'population:', ini, y%pop(ini)
!   enddo
!   do ini = 1, nlev  
!      do fin = 1, nlev
!         !if(a21%M(ini, fin).ne.0.d0)                                        &
!         write(6,'(2(i3, 2x), 4(e20.7))') ini, fin,                         &
!                                         energy%ene(ini), energy%ene(fin),   &
!                                         rad%M(ini, fin),                    &
!                                         coll_rad_matrix%A(ini, fin)
!      enddo
!   enddo

end program level_population
module level_population

    use types_and_parameters, only: nlev, energy_lev, vmax, jmax,         &
                                    ntrans, ntemp, i,                      &
                                    radiative_coeffs, collisional_coeffs, &
                                    reaction_matrix,                      &
                                    population,                           &
                                    Trad, nb, ini, fin, it,               &
                                    nlev_lique, vi, vf, ji, jf
                                    
    use energy_levels, only: reading_data_energies
    use radiation,     only: reading_data_radiative,                      &
                             radiative_downwards, radiative_upwards
    use collisions,    only: reading_data_collisions

    use linear_algebra, only: sparsity_calc,                              &
                              ndim, info, lda, ldb, nrhs, ipiv
    
    use matrix_construction, only: matrix_builder, initialize_level_population
    

    contains
    
    subroutine lev_pop(energy, a21, b21, r21, b12, jnu, r12, rr, rr21, rr12, x)
    
        type(energy_lev) :: energy
        type(radiative_coeffs) :: a21, b21, b12, jnu, r21, r12, rad
        type(collisional_coeffs) :: rr, rr21, rr12
        !type(reaction_matrix)  :: coll_rad_matrix
        type(population) :: x, y

        call reading_data_energies(energy)

        call reading_data_radiative(energy, a21)
    
        call reading_data_collisions(energy, rr, rr21, rr12)
        
        print*, 'T_Radiation', Trad
        
        call radiative_downwards(energy, Trad, a21, b21, r21)
        
        print*, 'T_Radiation', Trad
    
        call radiative_upwards(energy, Trad, a21, b21, b12, jnu, r12)

!         ! building the total radiative matrix, including both stimulated and spontaneous
!         ! transitions; according to the convention adopted:
!         !     downward transitions -> upper triangular matrix
!         !     upward transitions   -> lower triangular matrix
! 
!         rad%M = r21%M + r12%M
! 
!         it = 1
!         
!         call matrix_builder(rad, rr, it, coll_rad_matrix)
! 
!         call initialize_level_population(y)
! 
!         x%pop = y%pop
! 
!         print*, 'before', x%pop
!     
!         call dgesv(ndim, nrhs, coll_rad_matrix, lda, ipiv, x, ldb, info)
! 
!         do i = 1, nlev
!             write(6,'(3(i3), e14.7)') i, energy%vl(i), energy%jl(i), x%pop(i)
!         enddo
     
        return
    end subroutine lev_pop
    
    
    subroutine tests(energy, rr, rr21, rr12, a21, b21, r21, b12, jnu, r12)
        real*8 :: diagonal_a21, diagonal_b21, diagonal_r21
        real*8 :: diagonal_b12, diagonal_r12
        real*8 :: diagonal_rr21, diagonal_rr12        
        type(energy_lev) :: energy
        type(radiative_coeffs) :: a21, b21, r21, b12, r12, jnu!, rad
        type(collisional_coeffs) :: rr, rr21, rr12
        !type(reaction_matrix)  :: coll_rad_matrix
        !type(population) :: x, y    
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
    !  do i = 0, vmax-1
    !     do l = 0, jmax-1
    !        do m = 0, vmax-1
    !           do n = 0, jmax-1
    !           if(a21%arranging(i, l, m, n).ne.0.d0) then
    !     write(6,'(4(i2,2x),e20.7)') i, l, m, n, &
    !                   a21%arranging(i, l, m, n)
    !            endif
    !            enddo
    !        enddo
    !    enddo
    !  enddo
      !do i = 1, a21%ntransrad
      !   write(6,'(4(i2,2x),e14.5)') a21%vir(i), a21%jir(j), &
      !                               a21%vfr(l), a21%jfr(m), &
      !                 a21%M(a21%couple1r(i),a21%couple2r(i))
      !enddo      
      !print*, 'SUM(a21%M)', SUM(a21%M)
      !print*, 'max(a21%M)', maxval(a21%M)
      
      
    ! TEST COLLISIONAL COEFFICIENTS
    ! check of the coefficients compared to the read data by François + 
    ! detailed balance
    diagonal_rr21 = 0.d0
    diagonal_rr12 = 0.d0    
     do i=1, ntrans
        do it=1, ntemp
            if(rr%couple1c(i).ge.rr%couple2c(i))  then
                diagonal_rr21 = diagonal_rr21 + rr21%matrix_lique(rr%couple1c(i),rr%couple2c(i),it)
            else
                diagonal_rr12 = diagonal_rr12 + rr12%matrix_lique(rr%couple2c(i),rr%couple1c(i),it)
            endif        
    !          write(6, '(6(i3,2x), 4(e14.5))') rr%couple1c(i),rr%couple2c(i),   &
    !                                    rr%vic(i),rr%jic(i),rr%vfc(i),rr%jfc(i),&
    !                                    rr%temp(it),                            &
    !                     rr%matrix_lique(rr%couple1c(i),rr%couple2c(i),it)*1.d6,      &
    !                     rr%matrix_lique(rr%couple2c(i),rr%couple1c(i),it)*1.d6,      &
    !                          (rr%matrix_lique(rr%couple1c(i),rr%couple2c(i),it)-     &
    !                      rr%matrix_lique(rr%couple2c(i),rr%couple1c(i),it))*1.d6/    &
    !                      (rr%matrix_lique(rr%couple1c(i),rr%couple2c(i),it)*1.d6)
        enddo
    enddo
    print*, 'rr21', sum(rr21%matrix_lique), diagonal_rr21, 'max', maxval(rr21%matrix_lique)
    print*, 'rr12', sum(rr12%matrix_lique), diagonal_rr12, 'max', maxval(rr12%matrix_lique)
    
    print*, 'sum_collisional: ', sum(rr%matrix_lique)
    print*, 'max_collisional: ', maxval(rr%matrix_lique)
    print*, 'min_collisional: ', minval(rr%matrix_lique)    
    
    
    ! TEST RADIATIVE TRANSITIONS COEFFICIENTS
    diagonal_a21 = 0.d0
    diagonal_b21 = 0.d0
    diagonal_r21 = 0.d0
    diagonal_b12 = 0.d0
    diagonal_r12 = 0.d0
    diagonal_jnu = 0.d0    
      do ini = 1, nlev_lique
         do fin = 1, nlev_lique
    !        !if(a21%M(ini, fin).ne.0.d0)                                        &
    !        write(6,'(2(i3, 2x), 5(e14.7))')  ini, fin,                                     &
    !                                        energy%ene_lique(ini), energy%ene_lique(fin),   &
    !                                        a21%M_lique(ini, fin),                   &
    !                                        b21%M_lique(ini, fin),                   &
    !                                        r21%M_lique(ini, fin)!,                   &
    !                                        b12%M(ini, fin),                   &
    !                                        r12%M(ini, fin),                   &
    !                                        ! to have them into the same line although for the reverse transition:
    !                   b12%M(fin, ini),                   &
    !                                        r12%M(fin, ini),                   &
    !                                        rad%M(ini, fin)
            if(ini.ge.fin)  then
                diagonal_a21 = diagonal_a21 + a21%M_lique(ini, fin)   
                diagonal_b21 = diagonal_b21 + b21%M_lique(ini, fin)   
                diagonal_r21 = diagonal_r21 + r21%M_lique(ini, fin)
            else
                diagonal_b12 = diagonal_b12 + b12%M_lique(ini, fin)
                diagonal_r12 = diagonal_r12 + r12%M_lique(ini, fin)
            endif
            
            diagonal_jnu = diagonal_jnu + jnu%M_lique(ini, fin)
         enddo
      enddo
      print*, 'a21', sum(a21%M_lique), diagonal_a21
      write(6, '(a13,e14.7)') 'max frequency', maxval(energy%freq_lique)
      print*, 'b21', sum(b21%M_lique), diagonal_b21, 'max', maxval(b21%M_lique)
      print*, 'jnu', sum(jnu%M_lique), diagonal_jnu, 'max', maxval(jnu%M_lique)
      print*, 'r21', sum(r21%M_lique), diagonal_r21, 'max', maxval(r21%M_lique)

      print*, 'b12', sum(b12%M_lique), diagonal_b12, 'max', maxval(b12%M_lique)
      print*, 'r12', sum(r12%M_lique), diagonal_r12, 'max', maxval(r12%M_lique)
      print*, 'r21+r12', sum(r21%M_lique+r12%M_lique), 'max', maxval(r21%M_lique+r12%M_lique)
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

    ! TEST NUMBER OF TRANSITIONS IN COMMON AMONG RADIATIVE AND COLLISIONAL ONES  
    
     do i = 1, nlev_lique
        do j = 1, nlev_lique
            if(a21%couple1r(i).eq.rr%couple1c(j)) then
               if(a21%couple2r(j).eq.rr%couple2c(j)) then            
           write(6,'(a10, 6(i3, x), a12, 6(i3, x) )')                                                     &
                         'radiative:',                                                                    &
                         a21%couple1r(i), a21%couple2r(j), a21%vir(i), a21%jir(i), a21%vfr(j),a21%jfr(j), &
                         'collisional:',                                                                  &
                         rr%couple1c(i), rr%couple2c(j), rr%vic(i), rr%jic(i), rr%vfc(j),rr%jfc(j)
!            if(a21%M_lique(i, j).ne.0.d0) then
!               if(rr%matrix_lique(i, j, 1).ne.0.d0) then
!            write(6,'(6(i3, 2x), 2(e14.7, 2x))')                       &
!                          a21%couple1r(i), a21%couple2r(j),              &
!                          a21%vir(i),a21%jir(i),                         &
!                          a21%vfr(j),a21%jfr(j),                         &
!                          a21%M_lique(i, j),                           &
!                          rr%matrix_lique(i, j, 1)
              endif
           endif
        enddo
    enddo
!                    do i=1,a21%ntransrad
!                        vi=a21%vir(i)
!                        ji=a21%jir(i)
!                        vf=a21%vfr(i)
!                        jf=a21%jfr(i)
!                        !!if(a21%couple1r(i).ne.0) then
!                        !! if(a21%couple2r(i).ne.0) then
!                        !! if(energy%en(vi,ji).lt.energy%en(vf,jf)) then
!                        write(6,'(6(i3, 2x), 2(e14.7, 2x))')                         &
!                                      a21%couple1r(i),a21%couple2r(i),               &
!                                      a21%vir(i),a21%jir(i),                         &
!                                      a21%vfr(i),a21%jfr(i),                         &
!                                      a21%M_lique(a21%couple1r(i),a21%couple2r(i)),  &
!                                      a21%arranging(vi,ji,vf,jf)
!                        ! !endif
!                        ! !endif
!                        !!endif
!                    enddo

    return
    end subroutine tests

end module level_population
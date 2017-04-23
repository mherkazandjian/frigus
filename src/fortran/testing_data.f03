! This program writes the matrices on files with name corresponding to the name of the matrices

module testing_data
    use types_and_parameters, only: energy_lev, vmax, jmax,                  &
                                    ntrans, ntemp, i,                        &
                                    radiative_coeffs, collisional_coeffs,    &
                                    reaction_matrix,                         &
                                    Trad, nc, ini, fin,                      &
                                    nlev_lique, vi, vf, ji, jf, id_temp_test

 contains
 
    subroutine writing_files(a21, b21, b12, jnu, id_temp_test, rr, coll_rad_matrix)
        type(radiative_coeffs) :: a21, b21, b12, jnu
        type(collisional_coeffs) :: rr
        type(reaction_matrix) :: coll_rad_matrix
        character*59 :: command
        command = 'rm -fvr /home/carla/Dropbox/us/cooling_function/carla/*.txt'

        call system(command)

        open(30, file = '/home/carla/Dropbox/us/cooling_function/carla/A_matrix.txt', status = 'unknown')
        open(31, file = '/home/carla/Dropbox/us/cooling_function/carla/B_matrix.txt', status = 'unknown')
        open(32, file = '/home/carla/Dropbox/us/cooling_function/carla/J_nu_matrix.txt', status = 'unknown')
        open(33, file = '/home/carla/Dropbox/us/cooling_function/carla/K_nc_matrix.txt', status = 'unknown')
        open(34, file = '/home/carla/Dropbox/us/cooling_function/carla/O_matrix.txt', status = 'unknown')
        open(35, file = '/home/carla/Dropbox/us/cooling_function/carla/D_matrix.txt', status = 'unknown')        
        open(36, file = '/home/carla/Dropbox/us/cooling_function/carla/M_matrix.txt', status = 'unknown')        
        
        do j = 1, nlev_lique
           write(30, '(58(ES25.15E3, 1x))') (a21%M_lique(i, j), i = 1, nlev_lique)
           write(31, '(58(ES25.15E3, 1x))') (b21%M_lique(i, j)+b12%M_lique(i, j), i = 1, nlev_lique)
           write(32, '(58(ES25.15E3, 1x))') (jnu%M_lique(i, j), i = 1, nlev_lique)
           write(33, '(58(ES25.15E3, 1x))') (rr%matrix_lique(i, j, id_temp_test), i = 1, nlev_lique)
           write(34, '(58(ES25.15E3, 1x))') (coll_rad_matrix%O(i, j), i = 1, nlev_lique)
           write(35, '(58(ES25.15E3, 1x))') (coll_rad_matrix%D(i, j), i = 1, nlev_lique)
           write(36, '(58(ES25.15E3, 1x))') (coll_rad_matrix%M(i, j), i = 1, nlev_lique)
        enddo
        
        return
    end subroutine writing_files


    subroutine tests(energy, rr, rr21, rr12, a21, b21, r21, b12, jnu, r12, coll_rad_matrix, id_temp_test)
        real*8 :: diagonal_a21, diagonal_b21, diagonal_r21
        real*8 :: diagonal_b12, diagonal_r12
        real*8 :: diagonal_rr21, diagonal_rr12
        real*8 :: diagonal_jnu
        real*8 :: diagonal_M, diagonal_O, diagonal_D
        type(energy_lev) :: energy
        type(radiative_coeffs) :: a21, b21, r21, b12, r12, jnu
        type(collisional_coeffs) :: rr, rr21, rr12
        type(reaction_matrix)  :: coll_rad_matrix
        !type(population) :: x, y    

!----------------------------------------------------------------------------------------------------      

    ! TEST READING ENERGY LEVELS
    !!   call reading_data(energy)
!         do i = 0, vmax
!            do j = 0, jmax
    !    do i = 1, nlev_lique
    !       !do j = 1, nlev_lique
    !           !write(6,'(2(i2,2x), e14.8)') i, j, energy%en(i,j)
    !           write(6,'((i2,2x), e14.8)') i, energy%ene_lique(i)
    !       !enddo
    !    enddo

!----------------------------------------------------------------------------------------------------      

    ! TEST READING RADIATIVE COEFFICIENTS
    !!    call reading_data_radiative(a21)
    !    print*, a21%ntransrad
    !    do ini = 1, nlev_lique
    !        do fin = 1, nlev_lique
    !           write(6,'(2(i3, 2x), es7.1)') ini, fin, a21%M_lique(ini, fin)
    !        enddo
    !    enddo

!----------------------------------------------------------------------------------------------------      
    
    ! TEST COUPLES FOR RADIATIVE TRANSITIONS
    !    print*, a21%ntransrad, (jmax+1)*15*15*3
    !    do l=1,a21%ntransrad
    !        print*, a21%couple1r(l), a21%couple2r(l)
    !    enddo
    
!----------------------------------------------------------------------------------------------------      

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
      
!----------------------------------------------------------------------------------------------------      

    ! TEST COLLISIONAL COEFFICIENTS
    ! test 1: checks for comparison with the Python code
    ! diagonal_rr21 = 0.d0
    ! diagonal_rr12 = 0.d0    
    ! do i = 1, ntrans
    !    do itt = 1, ntemp
    !        if(rr%couple1c(i).ge.rr%couple2c(i))  then
    !            diagonal_rr21 = diagonal_rr21 + rr21%matrix_lique(rr%couple1c(i),rr%couple2c(i),itt)
    !        else
    !            diagonal_rr12 = diagonal_rr12 + rr12%matrix_lique(rr%couple2c(i),rr%couple1c(i),itt)
    !        endif
    !    enddo
    !enddo
    !print*, 'rr21', sum(rr21%matrix_lique), diagonal_rr21, 'max', maxval(rr21%matrix_lique)
    !print*, 'rr12', sum(rr12%matrix_lique), diagonal_rr12, 'max', maxval(rr12%matrix_lique)
    !print*, 'sum_collisional: ', sum(rr%matrix_lique),                        &
    !        'sum(rr21+rr12):',sum(rr21%matrix_lique + rr12%matrix_lique)
    !print*, 'max_collisional: ', maxval(rr%matrix_lique)
    !print*, 'min_collisional: ', minval(rr%matrix_lique)    
    
    ! test 2: checks with the read values from the input file + detailed balance
     do i = 1, ntrans
        do it = 1, 1!ntemp    
              write(6, '(6(i3,2x), 4(e14.5))') rr%couple1c(i),rr%couple2c(i),   &
                                        rr%vic(i),rr%jic(i),rr%vfc(i),rr%jfc(i),&
                                        rr%temp(it),                            &
                         rr%matrix_lique(rr%couple1c(i),rr%couple2c(i),it)*1.d6,      &
                         rr%matrix_lique(rr%couple2c(i),rr%couple1c(i),it)*1.d6,      &
                              (rr%matrix_lique(rr%couple1c(i),rr%couple2c(i),it)-     &
                          rr%matrix_lique(rr%couple2c(i),rr%couple1c(i),it))*1.d6/    &
                          (rr%matrix_lique(rr%couple1c(i),rr%couple2c(i),it)*1.d6)
        enddo
    enddo    

    ! test 3: check for collisional coefficients and corresponding temperature
    !     do i = 1, ntrans
    !           do it = 1, ntemp
    !              write(6,'(2(i2,2x), 2(e14.7,2x) )') rr%couple1c(i),rr%couple2c(i),   &
    !                                                  rr%temp(it),                     &
    !                                rr%matrix_lique(rr%couple1c(i),rr%couple2c(i), it)
    !           enddo
    !     enddo    
    
    
!----------------------------------------------------------------------------------------------------

    ! TEST RADIATIVE TRANSITIONS COEFFICIENTS
    ! test 1: checks for comparison with the Python code
    !diagonal_a21 = 0.d0
    !diagonal_b21 = 0.d0
    !diagonal_r21 = 0.d0
    !diagonal_b12 = 0.d0
    !diagonal_r12 = 0.d0
    !diagonal_jnu = 0.d0    
    !  do ini = 1, nlev_lique
    !     do fin = 1, nlev_lique
    !        if(ini.ge.fin)  then
    !            diagonal_a21 = diagonal_a21 + a21%M_lique(ini, fin)   
    !            diagonal_b21 = diagonal_b21 + b21%M_lique(ini, fin)   
    !            diagonal_r21 = diagonal_r21 + r21%M_lique(ini, fin)
    !        else
    !            diagonal_b12 = diagonal_b12 + b12%M_lique(ini, fin)
    !            diagonal_r12 = diagonal_r12 + r12%M_lique(ini, fin)
    !        endif
    !       
    !        diagonal_jnu = diagonal_jnu + jnu%M_lique(ini, fin)
    !     enddo
    !  enddo
    !  print*, 'a21', sum(a21%M_lique), diagonal_a21
    !  write(6, '(a13,e14.7)') 'max frequency', maxval(energy%freq_lique)
    !  print*, 'sum b21:', sum(b21%M_lique), diagonal_b21, 'max b21:', maxval(b21%M_lique)
    !  print*, 'sum b12:', sum(b12%M_lique), diagonal_b12, 'max b12:', maxval(b12%M_lique)
    !  print*, 'sum b21+b12:', sum(b21%M_lique + b12%M_lique), 'max b:', maxval(b21%M_lique + b12%M_lique)
    !  print*, 'sum jnu', sum(jnu%M_lique), diagonal_jnu, 'max', maxval(jnu%M_lique)
    !  print*, 'sum r21', sum(r21%M_lique), diagonal_r21, 'max', maxval(r21%M_lique)
    !  print*, 'sum r12', sum(r12%M_lique), diagonal_r12, 'max', maxval(r12%M_lique)
    !  print*, 'sum r21+r12', sum(r21%M_lique+r12%M_lique), 'max', maxval(r21%M_lique+r12%M_lique)
    
    ! test 2: to check the read data and filled matrices
    !  do ini = 1, nlev_lique
    !     do fin = 1, nlev_lique
    !!        !if(a21%M(ini, fin).ne.0.d0)                                        &
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
    !                                        rad%M(ini, fin) !!! to be defined

!---------------------------------------------------------------------------------------------------

    ! TEST MATRIX LINEAR SYSTEM
    ! test 1: checks for the off-diagonal terms of M
    !diagonal_M = 0.d0
    !diagonal_O = 0.d0
    !diagonal_D = 0.d0
    !  do ini = 1, nlev_lique
    !     do fin = 1, nlev_lique
    !            diagonal_O = diagonal_O + coll_rad_matrix%O(ini, fin)
    !            diagonal_D = diagonal_D + coll_rad_matrix%D(ini, fin)                
    !            diagonal_M = diagonal_M + coll_rad_matrix%M(ini, fin)                
    !     enddo
    !  enddo
    ! print*, 'sum O: ', sum(coll_rad_matrix%O), diagonal_O, 'max: ', maxval(coll_rad_matrix%O)
    ! print*, 'sum D: ', sum(coll_rad_matrix%D), diagonal_D, 'max: ', maxval(coll_rad_matrix%D)
    ! print*, 'min D: ', minval(coll_rad_matrix%D)
    ! write(6, '(a7, 2(e14.7, 2x), 2(a5, (e14.7, 2x)), 2(a11, e14.7, 2x) )')                    &
    !          'sum M: ', sum(coll_rad_matrix%M), diagonal_M,             &
    !          'max: ', maxval(coll_rad_matrix%M),                        &
    !          'min: ', minval(coll_rad_matrix%M),                        &
    !          'sum O + D: ', sum(coll_rad_matrix%O+coll_rad_matrix%D),   &
    !          'max O + D: ', maxval(coll_rad_matrix%M)
    !          
    
    
    
    ! test 2: population test
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

!----------------------------------------------------------------------------------------------------

    ! TEST NUMBER OF TRANSITIONS IN COMMON AMONG RADIATIVE AND COLLISIONAL ONES  
    
    ! do i = 1, nlev_lique
    !    do j = 1, nlev_lique
    !        if(a21%couple1r(i).eq.rr%couple1c(j)) then
    !           if(a21%couple2r(j).eq.rr%couple2c(j)) then            
    !       write(6,'(a10, 6(i3, x), a12, 6(i3, x) )')                                                     &
    !                     'radiative:',                                                                    &
    !                     a21%couple1r(i), a21%couple2r(j), a21%vir(i), a21%jir(i), a21%vfr(j),a21%jfr(j), &
    !                     'collisional:',                                                                  &
    !                     rr%couple1c(i), rr%couple2c(j), rr%vic(i), rr%jic(i), rr%vfc(j),rr%jfc(j)
    !!            if(a21%M_lique(i, j).ne.0.d0) then
    !!               if(rr%matrix_lique(i, j, 1).ne.0.d0) then
    !!            write(6,'(6(i3, 2x), 2(e14.7, 2x))')                       &
    !!                          a21%couple1r(i), a21%couple2r(j),              &
    !!                          a21%vir(i),a21%jir(i),                         &
    !!                          a21%vfr(j),a21%jfr(j),                         &
    !!                          a21%M_lique(i, j),                           &
    !!                          rr%matrix_lique(i, j, 1)
    !          endif
    !       endif
    !    enddo
    ! enddo
    !!                    do i=1,a21%ntransrad
    !!                        vi=a21%vir(i)
    !!                        ji=a21%jir(i)
    !!                        vf=a21%vfr(i)
    !!                        jf=a21%jfr(i)
    !!                        !!if(a21%couple1r(i).ne.0) then
    !!                        !! if(a21%couple2r(i).ne.0) then
    !!                        !! if(energy%en(vi,ji).lt.energy%en(vf,jf)) then
    !!                        write(6,'(6(i3, 2x), 2(e14.7, 2x))')                         &
    !!                                      a21%couple1r(i),a21%couple2r(i),               &
    !!                                      a21%vir(i),a21%jir(i),                         &
    !!                                      a21%M_lique(a21%couple1r(i),a21%couple2r(i)),  &
    !!                                      a21%arranging(vi,ji,vf,jf)
    !!                        ! !endif
    !!                        ! !endif
    !!                        !!endif
    !!                    enddo

    return
    end subroutine tests
     
     
     
end module testing_data
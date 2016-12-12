module collisions
    use energy_levels
    use types_and_parameters, only: energy_lev, collisional_coeffs, &
                                    it, nlev,       &
                                    vimax, jimax,   &
                                    vfmax, jfmax,   &
                                    ntemp, ntrans,  &
                                    vi, ji, vf, jf    
    ! in this module the reading, indexes arrangement and calculations of
    ! stimulated coefficients are performed, starting from the data by
    ! Wolniewicz et al 1998

    type(collisional_coeffs) :: rr
    
    contains

      subroutine reading_data_collisions(e, rr, rr21, rr12)
                 use energy_levels, only: reading_data_energies
                 use types_and_parameters, only: it, nlev, nlev_lique,  &
                                                 vimax, jimax,          &
                                                 vfmax, jfmax,          &
                                                 ntemp, ntrans,         &
                                                 vi, ji, vf, jf,        &
                                                 kb, nc, id_temp
 
                 type(collisional_coeffs) :: rr, rr21, rr12       ! reaction rate
                 type(energy_lev) :: e

                 real*8 :: dE

                 ! initializing the arrays
                 rr%reading = 0.d0
                 rr%matrix_lique = 0.d0
                 rr%temp = [ (i, i=100,5000,100) ]
                 !do i = 1, ntemp
                 !   print*, 'temperature', i, rr%temp(i)
                 !enddo

!                 open (20, file='Read/Rates_H_H2.dat', status = 'unknown')
                 open (20, file='Read/Rates_H_H2.dat', status = 'unknown')

!                 do i=1,10 ! data by francois
!                    read(20,*) 
!                 enddo
                 do i=1,8   ! data by wrathmall
                    read(20,*) 
                 enddo
                 do i=1,ntrans
                    read(20,*) rr%vic(i),rr%jic(i),rr%vfc(i),rr%jfc(i),      &
                    (rr%reading(rr%vic(i),rr%jic(i),rr%vfc(i),rr%jfc(i),it), &
                     it=1,ntemp)
                    vi=rr%vic(i)
                    ji=rr%jic(i)
                    vf=rr%vfc(i)
                    jf=rr%jfc(i)
                     do l=1, nlev_lique
                         if(vi.eq.e%vl_lique(l)) then
                             if(ji.eq.e%jl_lique(l)) then
                                rr%couple1c(i) = l
                             endif
                         endif
                         if(vf.eq.e%vl_lique(l)) then
                             if(jf.eq.e%jl_lique(l)) then
                                rr%couple2c(i) = l
                             endif
                         endif
                     enddo
                 enddo

                !do i=1,ntrans
                !   write(6,'(6(i3,2x))') rr%vic(i),rr%jic(i),rr%vfc(i),rr%jfc(i),      &
                !           rr%couple1c(i),rr%couple2c(i)
                !enddo

                ! detailed balance implementation
                 do i=1,ntrans
                    vi=rr%vic(i)
                    ji=rr%jic(i)
                    vf=rr%vfc(i)
                    jf=rr%jfc(i)
                    dE = abs(e%en_lique(vi,ji)-e%en_lique(vf,jf))
                     do it = 1, ntemp
                        rr%matrix_lique(rr%couple1c(i),rr%couple2c(i),it) = &
                                       rr%reading(vi,ji,vf,jf,it)
                        rr21%matrix_lique(rr%couple1c(i),rr%couple2c(i),it) = &
                          rr%matrix_lique(rr%couple1c(i),rr%couple2c(i),it)
                          
                        rr%matrix_lique(rr%couple2c(i),rr%couple1c(i),it) =       &
                          dexp(-dE/(kb*rr%temp(it))) * rr%reading(vi,ji,vf,jf,it) &
                          * ((2.*ji+1.))/(2.*jf+1.)
                        rr12%matrix_lique(rr%couple2c(i),rr%couple1c(i),it) = &
                          rr%matrix_lique(rr%couple2c(i),rr%couple1c(i),it)
                     enddo
                 enddo

                 !units conversion: cm3 s-1 -> m3 s-1
                 rr%matrix_lique = rr%matrix_lique*1.d-6
                 rr21%matrix_lique = rr21%matrix_lique*1.d-6                                  
                 rr12%matrix_lique = rr12%matrix_lique*1.d-6


      end subroutine reading_data_collisions


    end module collisions


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

      subroutine reading_data_collisions(e, rr)
                 use energy_levels, only: reading_data_energies
                 use types_and_parameters, only: it, nlev,       &
                                                 vimax, jimax,   &
                                                 vfmax, jfmax,   &
                                                 ntemp, ntrans,  &
                                                 vi, ji, vf, jf, &
                                                 kb
 
                 type(collisional_coeffs) :: rr                   ! reaction rate
                 type(energy_lev) :: e

                 real*8 :: dE

                 ! initializing the arrays
                 rr%reading = 0.d0
                 rr%matrix = 0.d0
                 rr%temp = (/ (i, i=100,5000,100) /)
                 print*, rr%temp

                 
                 call reading_data_energies(e)
                 
                 open (20, file='Read/Rates_H_H2.dat', status = 'unknown')

                 do i=1,10 
                    read(20,*) 
                 enddo
                 do i=1,ntrans
                    read(20,*) rr%vic(i),rr%jic(i),rr%vfc(i),rr%jfc(i),      &
                    (rr%reading(rr%vic(i),rr%jic(i),rr%vfc(i),rr%jfc(i),it), it=1,ntemp)
                    vi=rr%vic(i)
                    ji=rr%jic(i)
                    vf=rr%vfc(i)
                    jf=rr%jfc(i)
                     do l=1,nlev
                         if(vi.eq.e%vl(l)) then
                             if(ji.eq.e%jl(l)) then
                                rr%couple1c(i) = l
                             endif
                         endif
                         if(vf.eq.e%vl(l)) then
                             if(jf.eq.e%jl(l)) then
                                rr%couple2c(i) = l
                             endif
                         endif
                     enddo
                 enddo

!                do i=1,ntrans
!                   print*, rr%couple1c(i),rr%couple2c(i)
!                enddo
                 do i=1,ntrans
                    vi=rr%vic(i)
                    ji=rr%jic(i)
                    vf=rr%vfc(i)
                    jf=rr%jfc(i)
                    dE = abs(e%en(vi,ji)-e%en(vf,jf))
!                    print*, vi,ji,vf,jf,dE
                     do it=1,ntemp
                        rr%matrix(rr%couple1c(i),rr%couple2c(i),it) = &
                                       rr%reading(vi,ji,vf,jf,it)
                        rr%matrix(rr%couple2c(i),rr%couple1c(i),it) = &
                          dexp(-dE/(kb*rr%temp(it))) * rr%reading(vi,ji,vf,jf,it)
                     enddo
                 enddo
!                 print*, kb

                 !check of the coefficients compared to the read data by François
                 do i=1, ntrans
                     do it=1, ntemp
                        write(6, '(2(i3,2x), 4(e14.5))') rr%couple1c(i),rr%couple2c(i),   &
                                                     rr%temp(it),                         &
                                        rr%matrix(rr%couple1c(i),rr%couple2c(i),it),      &
                                        rr%matrix(rr%couple2c(i),rr%couple1c(i),it),      &
                                        (rr%matrix(rr%couple1c(i),rr%couple2c(i),it)-     &
                                         rr%matrix(rr%couple2c(i),rr%couple1c(i),it))/    &
                                         rr%matrix(rr%couple1c(i),rr%couple2c(i),it)
                     enddo
                 enddo

                 !units conversion: cm3 s-1 -> m3 s-1
                 rr%matrix = rr%matrix*1.d-6
 
      end subroutine reading_data_collisions


    end module collisions


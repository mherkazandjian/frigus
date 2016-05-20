module radiation
    
    use energy_levels
    use types_and_parameters, only: jmax, vi, ji, vf, jf, &
                                    radiative_coeffs, energy_lev

    ! in this module the reading, indexes arrangement and calculations of
    ! stimulated coefficients are performed, starting from the data by
    ! Wolniewicz et al 1998

    type(radiative_coeffs) :: a21, b21, b12
    type(energy_lev) :: e
    
    contains

      subroutine reading_data_radiative(e, a21)
                 use types_and_parameters, only: jmax, nlev,       &
                                                 energy_lev
 
                 type(radiative_coeffs) :: a21
                 type(energy_lev) :: e 
                 integer, dimension(0:jmax)  :: ivmax


                 a21%reading = 0.d0
                 a21%M = 0.d0
                 
                 call tableh2(a21%reading, ivmax)

                 !print*, e%en
                 
                 a21%ntransrad = 0
                !write(6,'(a24,2x,e10.4)') 'a(-1:1,0:14,0:14,0:jmax)',a(1,0,1,0)
                 do i0=0,jmax
                    do i1=0,14
                        do i2=0,14
                            do i3=-1,1
                               if(a21%reading(i3,i2,i1,i0).ne.0.d0) then
                               a21%ntransrad = a21%ntransrad + 1
                                write(22,'(4(i2,2x),e10.4)') i1,i0,i2,i0+i3,a21%reading(i3,i2,i1,i0)
                               endif
                            enddo
                        enddo
                    enddo
                enddo
                rewind(22)
                allocate(a21%couple1r(1:a21%ntransrad),a21%couple2r(1:a21%ntransrad))
                allocate(a21%vir(1:a21%ntransrad),a21%jir(1:a21%ntransrad))
                allocate(a21%vfr(1:a21%ntransrad),a21%jfr(1:a21%ntransrad))
                allocate(a21%arranging(0:vmax, 0:jmax, 0:vmax, 0:jmax))
                ! print*, a21%ntransrad, (jmax+1)*15*15*3
! 
                 do i=1,a21%ntransrad
                     read(22,*)  a21%vir(i),a21%jir(i),a21%vfr(i),a21%jfr(i), &
                     a21%arranging(a21%vir(i),a21%jir(i),a21%vfr(i),a21%jfr(i))
                     !write(6,'(4(i3,2x),e10.4)')  a21%vir(i),a21%jir(i),a21%vfr(i),a21%jfr(i), &
                     !a21%arranging(a21%vir(i),a21%jir(i),a21%vfr(i),a21%jfr(i))
                     vi=a21%vir(i)
                     ji=a21%jir(i)
                     vf=a21%vfr(i)
                     jf=a21%jfr(i)
                     do l=1,nlev
                         if(vi.eq.e%vl(l)) then
                           if(ji.eq.e%jl(l)) then
                             a21%couple1r(i) = l
                           endif
                         endif
                         if(vf.eq.e%vl(l)) then
                           if(jf.eq.e%jl(l)) then
                             a21%couple2r(i) = l
                           endif
                         endif
                     enddo
                 enddo

                  do i=1,a21%ntransrad
                      vi=a21%vir(i)
                      ji=a21%jir(i)
                      vf=a21%vfr(i)
                      jf=a21%jfr(i)
                      !a21%M(a21%couple1r(i),a21%couple2r(i)) = a21%reading(jf,vf,vi,ji)
                      a21%M(a21%couple1r(i),a21%couple2r(i)) = a21%arranging(vi,ji,vf,jf)
                  enddo
                !print*, ivmax
                return
      end subroutine reading_data_radiative


      subroutine tableh2(a,ivmax)

        implicit real*8(a-h,o-z)
        parameter(jmax=31)
        integer ivmax
        dimension a(-1:1,0:14,0:14,0:jmax),ivmax(0:jmax)
        !*****************************************************************************
        !     The subroutine reads the quadrupole transition probabilities 'a' from  *
        !     3 data-files: 'j2jdown,' 'j2j,' and 'j2jup.' Thus, a(k,ivf,ivi,j) is   *
        !     the probability for the transition from an initial (upper) state with  *
        !     a final (lower) state with vibrational quantum number 'ivf' and        *
        !     rotational quantum number equal to j-2 for k=-1, j for k=0, and        *
        !     j+2 for k=1.                                                           *
        !                                                                            *
        !     The array elements of 'a' for which the subscripts do NOT correspond   *
        !     to a quadrupole transition are put equal to ZERO.                      *
        !
        !     The subroutine also gives the maximum of the vibrational number for    *
        !     each rotational quantum number j. Thus, for a fixed value j, the       *
        !     value for j is jmax=31, in which case there is only one possible       *
        !     vibrational state, i.e. ivmax(31)=0.                                   *
        !*****************************************************************************

        do j=0,jmax
            do i1=0,14
                do i2=0,14
                    do k=-1,1
                        a(k,i2,i1,j)=0.d0
                    enddo
                enddo
            enddo
        enddo

        open(16,file='Read/j2jdown')
        read(16,*)k
        read(16,*) (jj,j=0,jmax)
        read(16,*) (ivmax(j),j=0,jmax)
        do j=2,jmax
            read(16,*)(ii,i=ivmax(j),1,-1)
            do ivf=0,ivmax(j)
                read(16,*)jj,ii,(a(-1,ivf,i,j),i=ivmax(j),ivf,-1)
            enddo
        enddo
        close(16)
        open(17,file='Read/j2j')
        read(17,*)k
        do j=1,jmax-1
            read(17,*)(ii,i=ivmax(j),1,-1)
            do ivf=0,ivmax(j)-1
                read(17,*)jj,ii,(a(0,ivf,i,j),i=ivmax(j),ivf+1,-1)
            enddo
        enddo
        close(17)
        open(18,file='Read/j2jup')
        read(18,*)k
        do j=0,jmax-2
            read(18,*)(ii,i=ivmax(j),1,-1)
            do ivf=0,min0(ivmax(j+2),ivmax(j)-1)
                read(18,*)jj,ii,(a(1,ivf,i,j),i=ivmax(j),ivf+1,-1)
                do ivi=ivmax(j),ivf+1,-1
                    if(a(1,ivf,ivi,j).lt.0.d0)then
                        a(-1,ivi,ivf,j+2)=-a(1,ivf,ivi,j)
                        a(1,ivf,ivi,j)=0.d0
                    endif
                enddo
            enddo
        enddo
        close(18)
        return
      end subroutine tableh2

      subroutine radiative_downwards(energy, Tr, a21, b21, r21)
                 ! this subroutine returns the Einstein coefficients
                 ! for downward stimulated transitions
                 ! input: a21
                 ! output: b21 = c**2/(2*hp*nu)
                 use types_and_parameters, only: hp, c, nlev,  &
                                                 energy_lev,   &
                                                 ini, fin

                 type(radiative_coeffs) :: a21, b21, r21
                 type(energy_lev) :: energy

                 real*8 :: Tr


                 b21%M = 0.d0
                 r21%M = 0.d0
                 
                 do ini = 1, nlev
                    do fin = 1, nlev
                        if((energy%ene(ini)-energy%ene(fin)).lt.0.d0) then  ! according to the ordering in the  
                                                                            ! energy file downwards => E2 - E1 < 0
                           if(ini.ne.fin) then
                                b21%M(ini, fin) = c**2/(2.d0*hp*(energy%freq(ini, fin))**3)* &
                                                  a21%M(ini,fin)
                                r21%M(ini, fin) = a21%M(ini, fin) +    &
                                                  b21%M(ini, fin) *    &
                                             planck(energy%freq(ini, fin), Tr)
                                !write(6,'(a26, 2(i3,2x), (e16.10,2x), a8, 3(e16.10,2x))') 'from radiative downwards: ', &
                                !                                         ini, fin,                                      &
                                !                                         b21%M(ini, fin),                               &
                                !                                         'planck: ',                                    &
                                !                                         planck(energy%freq(ini, fin), Tr),             &
                                !                                         r21%M(ini, fin),                               &
                                !                                         a21%M(ini, fin)
                           endif
                        endif
                    enddo
                enddo
                return
      end subroutine radiative_downwards 

      subroutine radiative_upwards(energy, Tr, a21, b12, r12)
                 use types_and_parameters, only: hp, nlev,   &
                                                 energy_lev, &
                                                 ini, fin
                 real*8  :: Tr
                 real*8  :: g1, g2
                 type(radiative_coeffs) :: a21, b21, b12, r21, r12
                 type(energy_lev) :: energy
                 

                 b12%M = 0.d0
                 r12%M = 0.d0
                 
                 call radiative_downwards(energy, Tr, a21, b21, r21)
                 !print*, 'a21%M', a21%M
                 !print*, 'b21%M', b21%M
                 
                 do ini = 1, nlev
                     do fin = 1, nlev
                        if((energy%ene(ini)-energy%ene(fin)).gt.0.d0) then  ! according to the ordering in the  
                                                                            ! energy file downwards => E2 - E1 < 0                             
                                if(energy%jl(ini).ne.0) then
                                   g2 = 2.*energy%jl(fin) + 1.
                                   g1 = 2.*energy%jl(ini) + 1.
                                   
                                   b12%M(ini, fin) = (g2/g1)*     &
                                                        b21%M(fin, ini)
                                   r12%M(ini, fin) = b12%M(ini, fin)*                     &
                                                     planck(energy%freq(ini, fin), Tr)
                                endif
                                !write(6,'(a24, 2(i3,2x), (e10.4, 2x), a8, 2(e10.4,2x))')                 &
                                !                                       'from radiative upwards: ', &
                                !                                        ini, fin, b12%M(ini, fin), &
                                !                                        'planck: ',                &
                                !                                planck(energy%freq(ini, fin), Tr), &
                                !                                        r12%M(ini, fin)
                             endif
                     enddo
                 enddo
                 !print*, 'b12%M', b12%M
                 
      end subroutine radiative_upwards



      real*8 function planck(ni, T)
    !    computes the planck function:
    !    inputs: frequencies, Trad
         use types_and_parameters, only: c, kb, hp
         real*8 :: xx, yy, ni, T
         xx = hp * ni / (kb * T)
         yy = (2.d0*hp*ni**3)/c**2
         if(xx.lt.0.6739d3) then 
            planck = yy* 1.0d0 / (dexp(xx) - 1.0d0)
         else
            planck = 1.d-250
         endif
         !write(6,'(a13, 3(e10.4, 2x), a8, e10.4)') 'from planck: ', xx, yy, ni, 'planck  ', planck
         end      
      
      real*8 function planckOccupation(hp, nu, kb, T)
    !    computes the planck function for input parameters.
    !    keywords: nu, T
         !PHYSICAL PARAMETERS
         real*8 :: kb, hp, x, nu, T
         x = hp * nu / (kb * T)
         planckOccupation = 1.0d0 / (dexp(x) - 1.0d0)
         end

! 
!                 !do i=1,nradtrans
!                 !   vi=vrad(i)
!                 !   ji=jrad(i)
!                 !   vf=vprad(i)
!                 !   jf=jprad(i)
!                 !   if(en(vi,ji).gt.en(vf,jf)) then 
!                 !      AA(coupler1(i),coupler2(i)) = a(vi,ji,vf,jf)+
!                 !   else
!                 !      AA(coupler1(i),coupler2(i)) =
!                 !   endif
!                 !!   write(6,*) vi,ji,vf,jf,a(vi,ji,vf,jf)
!                 !enddo


    end module radiation


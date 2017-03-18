module energy_levels
                   
       ! module that reads, converts to the appropriate units and 
       ! and store the energy levels of H2 (rovibrational); 
       ! input data from Stancil
       use types_and_parameters, only: hp, kb, nlev, jmax, vmax, energy_lev, &
                                       ilique_flag
       use sorting, only: piksrt

    contains

        subroutine reading_data_energies(e)
                   use types_and_parameters, only: nlev, jmax, vmax,       &
                                                   energy_lev, hp, kb,     &
                                                   ini, fin, q,            &
                                                   jmax_lique, vmax_lique, &
                                                   nlev_lique, ilique_flag
                   type(energy_lev) :: e
                   
                   
                   open(10, file='Read/H2Xvjlevels.cs', status = 'unknown')
                   open(11, file='Read/H2Xvjlevels_francois_mod.cs', status = 'unknown') ! modified version of
!                                                                                           ! data from lique;
!                                                                                           ! keeping only the
!                                                                                           ! actual levels given in ! the reaction rates 
                   open(12, file='Read/H2Xvjlevels_flower.cs', status = 'unknown') ! modified version of
                                                                                          ! data from lique;
                                                                                          ! keeping only the
                                                                                          ! actual levels given in ! the reaction rates 
                                                                                          ! file
                                                                 ! file
                   open(13, file='Read/lev_labels', status = 'unknown')
                   open(14, file='Read/lev_labels_lique', status = 'unknown')
                   open(15, file='Read/frequencies_lique', status = 'unknown')
                   
                    allocate(e%en(0:vmax,0:jmax))
                    allocate(e%ene(1:nlev))
                    allocate(e%vl(1:nlev),e%jl(1:nlev))
                    allocate(e%freq(1:nlev,1:nlev))
                    
                    allocate(e%en_lique(0:vmax_lique,0:jmax_lique))
                    allocate(e%ene_lique(1:nlev_lique))
                    allocate(e%vl_lique(1:nlev_lique),e%jl_lique(1:nlev_lique))
                    allocate(e%freq_lique(1:nlev_lique,1:nlev_lique))

                    e%en  = 0.d0
                    e%ene = 0.d0
                    e%en_lique  = 0.d0
                    e%ene_lique = 0.d0
                    
                    ! reading stancil data
                    do i = 1, 10
                        read(10,*)
                    enddo
                    do i = 1, nlev
                        read(10,*) e%vl(i),e%jl(i), b, e%en(e%vl(i),e%jl(i))
                        e%ene(i) = e%en(e%vl(i),e%jl(i))
                        !   write(6,'(i3,2x,i2,2x,i2,2x,e10.4)') i, vl(i), jl(i), en(vl(i),jl(i))
                    enddo
                    !conversion cm-1 -->  Joule
                    e%en = e%en*1.98630d-23
                    e%ene = e%ene*1.98630d-23

                    ! evaluation of the frequencies
                    e%freq = 0.d0
                    do ini = 1, nlev
                        do fin = 1, nlev
                            e%freq(ini, fin) = dabs(e%ene(ini)-e%ene(fin))/hp
                        enddo
                    enddo

                    print*, 'ilique_flag in energy levels file:', ilique_flag
                   
                   
                     if(ilique_flag.eq.1) then
                      ! reading lique data
                      do i = 1, 3
                          read(11,*)
                      enddo
                      do i = 1, nlev_lique
                          read(11,*) b, b, e%vl_lique(i), e%jl_lique(i), b,    &
                                     e%en_lique(e%vl_lique(i), e%jl_lique(i))
                          e%ene_lique(i) = e%en_lique(e%vl_lique(i), e%jl_lique(i))
                          !   write(6,'(i3,2x,i2,2x,i2,2x,e10.4)') i, vl(i), jl(i), en(vl(i),jl(i))
                      enddo
                      !conversion eV -->  Joule
                      e%en_lique = e%en_lique*q
                      e%ene_lique = e%ene_lique*q
  
                      e%ene = -e%ene   ! to revert the actual way the energies are given
                                       ! Emax @ (v=0,j=0) --> Emin @ (v=0,j=0)
  
                                       
                      ! ordering the levels according to their energies 
                      ! (for subsequent construction of the matrix in the 
                      ! linear system of equations to be solved)
                       call piksrt(nlev, nlev_lique, e)                                     
  
                      ! evaluation of the frequencies
                      e%freq_lique = 0.d0
                      do ini = 1, nlev_lique
                          do fin = 1, nlev_lique
                              e%freq_lique(ini, fin) = dabs(e%ene_lique(ini)-e%ene_lique(fin))/hp
                          enddo
                      enddo                     
                     else
                      ! reading flower data
                      do i = 1, 1
                          read(12,*)
                      enddo
                      do i = 1, nlev_lique
                          read(12,*) b, e%vl_lique(i), e%jl_lique(i),     &
                                     e%en_lique(e%vl_lique(i), e%jl_lique(i))
                          e%ene_lique(i) = e%en_lique(e%vl_lique(i), e%jl_lique(i))
                          !   write(6,'(i3,2x,i2,2x,i2,2x,e10.4)') i, vl(i), jl(i), en(vl(i),jl(i))
                      enddo
                      !conversion K -->  Joule
                      e%en_lique = e%en_lique*kb
                      e%ene_lique = e%ene_lique*kb
  
                      e%ene = -e%ene   ! to revert the actual way the energies are given
                                       ! Emax @ (v=0,j=0) --> Emin @ (v=0,j=0)
  
                                       
                      ! ordering the levels according to their energies 
                      ! (for subsequent construction of the matrix in the 
                      ! linear system of equations to be solved)
                       call piksrt(nlev, nlev_lique, e)                                     
  
                      ! evaluation of the frequencies
                      e%freq_lique = 0.d0
                      do ini = 1, nlev_lique
                          do fin = 1, nlev_lique
                              e%freq_lique(ini, fin) = dabs(e%ene_lique(ini)-e%ene_lique(fin))/hp
                          enddo
                      enddo                     
                     endif
                    
                     do i=1,nlev                     
                        write(13,'(3(i3,2x))') i, e%vl(i), e%jl(i)
                     enddo

                     write(14,'(a5, 2x, 2(a1,4x))') 'label', 'v' , 'j'
                     do i=1,nlev_lique
                        write(14,'(3(i3,2x),e14.7)') i, e%vl_lique(i), e%jl_lique(i), e%ene_lique(i)/q
                     enddo
                     
                     write(15,'(a5, 2x, 4(a2,3x))') 'label', 'v' , 'j', 'vp' , 'jp'
                     do i=1,nlev_lique
                        do j=1,nlev_lique                     
                            write(15,'(6(i3,2x),e14.7)') i, j,                            &
                                                            e%vl_lique(i), e%jl_lique(i), &
                                                            e%vl_lique(j), e%jl_lique(j), & 
                                                            e%freq_lique(i,j)
                        enddo
                     enddo                     

                    return
        end subroutine reading_data_energies

end module energy_levels



!                     do i=1,nlev
!                        write(6,'(i3,2x,2(i2,2x),e10.4)')                        &
!                                                          i, e%vl(i), e%jl(i),   &
!                                                             e%ene(i)
!                        !do j=1,nlev
!                        ! write(6,'(e10.4)') e%freq(i, j)
!                        !enddo
!                     enddo
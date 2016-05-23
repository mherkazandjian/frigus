module energy_levels
                   
       ! module that reads, converts to the appropriate units and 
       ! and store the energy levels of H2 (rovibrational); 
       ! input data from Stancil
       use types_and_parameters, only: hp, nlev, jmax, vmax, energy_lev
       use sorting, only: piksrt

    contains

        subroutine reading_data_energies(e)
                   use types_and_parameters, only: nlev, jmax, vmax, & 
                                                   energy_lev, hp,   &
                                                   ini, fin
                   type(energy_lev) :: e
                   
                   open (21, file='Read/H2Xvjlevels.cs', status = 'unknown')
                   open (23, file='Read/lev_labels', status = 'unknown')

                    allocate(e%en(0:vmax,0:jmax))
                    allocate(e%ene(1:nlev))
                    allocate(e%vl(1:nlev),e%jl(1:nlev))
                    allocate(e%freq(1:nlev,1:nlev))

                    e%en  = 0.d0
                    e%ene = 0.d0

                    do i=1,10
                        read(21,*)
                    enddo
                    do i=1,nlev
                        read(21,*) e%vl(i),e%jl(i), b, e%en(e%vl(i),e%jl(i))
                        e%ene(i) = e%en(e%vl(i),e%jl(i))
                        write(23,'(3(i3,2x))') i, e%vl(i),e%jl(i)
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

                    ! ordering the levels according to their energies 
                    ! (for subsequent construction of the matrix in the 
                    ! linear system of equations to be solved)
                    call piksrt(nlev, e)
                    
                    e%ene = -e%ene   ! to revert the actual way the energies are given
                                     ! Emax @ (v=0,j=0) --> Emin @ (v=0,j=0)

                    !do i=1,nlev
                    !   write(6,'(i3,2x,2(i2,2x),e10.4)') i, e%vl(i), e%jl(i), &
                    !                                        e%ene(i)
                    !   !do j=1,nlev
                    !   ! write(6,'(e10.4)') e%freq(i, j)
                    !   !enddo
                    !enddo

                    
                    
                    return
        end subroutine reading_data_energies

end module energy_levels
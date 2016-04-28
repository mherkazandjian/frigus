module energy_levels

       ! module that reads, converts to the appropriate units and 
       ! and store the energy levels of H2 (rovibrational); 
       ! input data from Stancil

       real*8, dimension(:,:), allocatable         :: en
       real*8, dimension(:), allocatable           :: ene
       integer, allocatable, dimension(:)          :: vl, jl     ! labels for the levels (for the final ordering)
       real*8                                      :: dE
       integer, parameter                          :: nlev = 301 ! # of rovibrational levels according to Stancil
       integer, parameter                          :: jmax = 31
       integer, parameter                          :: vmax = 14

    contains

        subroutine reading_data(vl, jl, en, ene)
            real*8, dimension(:,:), allocatable         :: en
            real*8, dimension(:), allocatable           :: ene
            integer, allocatable, dimension(:)          :: vl, jl     ! labels for the levels (for the final ordering)
            open (21, file='Read/H2Xvjlevels.cs', status = 'unknown')
            open (23, file='Read/lev_labels', status = 'unknown')


            allocate(en(0:vmax,0:jmax))
            allocate(ene(1:nlev))
            allocate(vl(1:nlev),jl(1:nlev))

            en  = 0.d0
            ene = 0.d0
            
            do i=1,10
                read(21,*)
            enddo
            do i=1,nlev
                read(21,*) vl(i),jl(i), b, en(vl(i),jl(i))
                ene(i) = en(vl(i),jl(i))
                write(23,'(3(i3,2x))') i, vl(i),jl(i)
                !   write(6,'(i3,2x,i2,2x,i2,2x,e10.4)') i, vl(i), jl(i), en(vl(i),jl(i))
            enddo
            !conversion cm-1 -->  Joule
            en = en*1.98630d-23
            ene = ene*1.98630d-23
            return
            end subroutine reading_data

end module energy_levels
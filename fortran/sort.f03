module sorting

   use types_and_parameters, only: energy_lev

   CONTAINS


   SUBROUTINE piksrt(n1, n2, arr)
    !Sorts an array arr(1:n) into ascending numerical order, by straight insertion. n is input;
    !arr is replaced on output by its sorted rearrangement.
    ! n1, n2 etc according to how many energy levels sets we have (this case stancil + lique)

    type(energy_lev) :: arr
    INTEGER n1, n2
    INTEGER in, jn
    REAL*8  a, a_lique
    integer b, c, b_lique, c_lique
    !arr%vl = 0
    !arr%jl = 0
    do jn = 2, n1
         !Pick out each element in turn.
         a = arr%ene(jn)
         b = arr%vl(jn)
         c = arr%jl(jn)
         do in = jn-1, 1, -1
         !Look for the place to insert it.
            if(arr%ene(in).le.a) goto 10
            arr%ene(in+1) = arr%ene(in)
            arr%vl(in+1) = arr%vl(in)
            arr%jl(in+1) = arr%jl(in)
         enddo
         in = 0
10       arr%ene(in+1) = a          !Insert it.
         arr%vl(in+1) = b
         arr%jl(in+1) = c
    enddo


    do jn = 2, n2
         !Pick out each element in turn.
         a_lique = arr%ene_lique(jn)
         b_lique = arr%vl_lique(jn)
         c_lique = arr%jl_lique(jn)
         do in = jn-1, 1, -1
         !Look for the place to insert it.
            if(arr%ene_lique(in).le.a_lique) goto 11
            arr%ene_lique(in+1) = arr%ene_lique(in)
            arr%vl_lique(in+1) = arr%vl_lique(in)
            arr%jl_lique(in+1) = arr%jl_lique(in)
         enddo
         in = 0
11       arr%ene_lique(in+1) = a_lique          !Insert it.
         arr%vl_lique(in+1)  = b_lique
         arr%jl_lique(in+1)  = c_lique
    enddo

    return
    END



END module  sorting
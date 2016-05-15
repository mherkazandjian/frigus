module sorting

   use types_and_parameters, only: energy_lev

   CONTAINS


   SUBROUTINE piksrt(n, arr)
    !Sorts an array arr(1:n) into ascending numerical order, by straight insertion. n is input;
    !arr is replaced on output by its sorted rearrangement.

    INTEGER n
    type(energy_lev) :: arr
    INTEGER in, jn
    REAL*8 a
    integer b, c
    !arr%vl = 0
    !arr%jl = 0
    do jn = 2, n
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
            !print*, 'case 2', i+1, j
            !temp1(i+1) = arr%vl(i)
            !temp2(i+1) = arr%jl(i)
         enddo
         in = 0
10       arr%ene(in+1) = a          !Insert it.
         arr%vl(in+1) = b
         arr%jl(in+1) = c
         !print*, 'case 1', i+1, j
         !temp1(i+1) = arr%vl(j)
         !temp2(i+1) = arr%jl(j)
    enddo
    return
    END



END module  sorting
       program interp

       implicit real*8 (a-h,o-z)

c - Load grid for interpolation

       call initWH2

c - Ask for data set

c       print *, " Temperature?"
c       read *, t_u
c       tt = log10(t_u)
c       print *, " tt = ", tt

       print *, " Density?"
       read *, d_u
       dd = log10(d_u)
       print *, " dd = ", dd

       print *, " H/H2?"
       read *, h_u
       hh = log10(h_u)
       print *, " hh = ", hh

       print *, " O/P?"
       read *, r_u
       rr = r_u
       print *, " rr = ", rr

c - give interpolated value

       open(61,file='a.dat',status='unknown')
       do i=200,430
       tt = i/100.
       call coolH2(tt, dd, hh, rr, WH2)
       write(61,99) tt, WH2
c       print *, tt, WH2
c       print *, " H2 cooling (erg s-1 mol-1)", WH2, 10.0d0**WH2
       enddo

      close(61)

       stop
   99  format(f5.2, 1x, f7.2)
       end

c------------------------------------------------------------------------

       subroutine coolH2(tt, dd, hh, rr, WH2)

c - Compute interpolated value in 4D-grid

       implicit real*8 (a-h,o-z)

       parameter (npt=23)
       parameter (npd=49)
       parameter (nph=36)
       parameter (npr=6)

       real*8 TN(npt), DN(npd)
       real*8 RH(nph), ROP(npr)

       real*8 WG(npt,npd,nph,npr)

       common / grid / TN, DN, RH, ROP, WG

       data un, dix / 1.0d0, 10.0d0 /

c - Check that the point is inside our grid

cc       print *, tt, dd, hh, rr

       if (tt .lt. TN(1) .or.
     1     tt .gt. TN(npt)) then
         print *, " Value out of range (tt)!"
         print *, " T should be in the range:", dix**TN(1)
     1          , " - ", dix**TN(npt)
       endif

       if (dd .lt. DN(1) .or.
     1     dd .gt. DN(npd)) then
         print *, " Value out of range (dd)!"
         print *, " nH should be in the range:"
     1          , dix**DN(1), " - ", dix**DN(npd)
       endif

       if (hh .lt. RH(1) .or.
     1     hh .gt. RH(nph)) then
         print *, " Value out of range (hh)"
         print *, " H/H2 should be in the range:"
     1          , dix**RH(1), " - ", dix**RH(nph)
       endif

       if (rr .lt. ROP(1) .or.
     1     rr .gt. ROP(npr)) then
         print *, " Value out of range (rr)"
         print *, " O/P should be in the range:"
     1          , ROP(1), " - ", ROP(npr)
       endif

c  Find place in grid

       do 10 i=1,npt
         if (TN(i) .gt. tt) then
           it = i-1
           goto 11
         endif
   10  continue
       it = npt - 1
   11  continue

       do 20 i=1,npd
         if (DN(i) .gt. dd) then
           id = i-1
           goto 21
         endif
   20  continue
       id = npd - 1
   21  continue

       do 30 i=1,nph
         if (RH(i) .gt. hh) then
           ih = i-1
           goto 31
         endif
   30  continue
       ih = nph - 1
   31  continue

       do 40 i=1,npr
         if (ROP(i) .gt. rr) then
           ir = i-1
           goto 41
         endif
   40  continue
       ir = npr - 1
   41  continue

cc       print *, it, id, ih, ir

c  Do a simple linear interpolation

       y0000 = WG(it  ,id  ,ih  ,ir  )
       y1000 = WG(it+1,id  ,ih  ,ir  )
       y0100 = WG(it  ,id+1,ih  ,ir  )
       y1100 = WG(it+1,id+1,ih  ,ir  )
       y0010 = WG(it  ,id  ,ih+1,ir  )
       y1010 = WG(it+1,id  ,ih+1,ir  )
       y0110 = WG(it  ,id+1,ih+1,ir  )
       y1110 = WG(it+1,id+1,ih+1,ir  )
       y0001 = WG(it  ,id  ,ih  ,ir+1)
       y1001 = WG(it+1,id  ,ih  ,ir+1)
       y0101 = WG(it  ,id+1,ih  ,ir+1)
       y1101 = WG(it+1,id+1,ih  ,ir+1)
       y0011 = WG(it  ,id  ,ih+1,ir+1)
       y1011 = WG(it+1,id  ,ih+1,ir+1)
       y0111 = WG(it  ,id+1,ih+1,ir+1)
       y1111 = WG(it+1,id+1,ih+1,ir+1)

       t = (tt - TN(it)) / (TN(it+1) - TN(it))
       u = (dd - DN(id)) / (DN(id+1) - DN(id))
       v = (hh - RH(ih)) / (RH(ih+1) - RH(ih))
       w = (rr - ROP(ir)) / (ROP(ir+1) - ROP(ir))

       WH2 = (un-t) * (un-u) * (un-v) * (un-w) * y0000
     1     +    t   * (un-u) * (un-v) * (un-w) * y1000
     1     + (un-t) *    u   * (un-v) * (un-w) * y0100
     1     +    t   *    u   * (un-v) * (un-w) * y1100
     1     + (un-t) * (un-u) *    v   * (un-w) * y0010
     1     +    t   * (un-u) *    v   * (un-w) * y1010
     1     + (un-t) *    u   *    v   * (un-w) * y0110
     1     +    t   *    u   *    v   * (un-w) * y1110
     1     + (un-t) * (un-u) * (un-v) *    w   * y0001
     1     +    t   * (un-u) * (un-v) *    w   * y1001
     1     + (un-t) *    u   * (un-v) *    w   * y0101
     1     +    t   *    u   * (un-v) *    w   * y1101
     1     + (un-t) * (un-u) *    v   *    w   * y0011
     1     +    t   * (un-u) *    v   *    w   * y1011
     1     + (un-t) *    u   *    v   *    w   * y0111
     1     +    t   *    u   *    v   *    w   * y1111

       return
       end

c----------------------------------------------------------------------------------

       subroutine initWH2

       implicit real*8 (a-h,o-z)

c  H2 global cooling
c      W = W(T, nH, H/H2, O/P)

c  npt : Nb of points in Temperature grid (TN)
c  npd : Nb of points in Density grid (DN)
c  nph : Nb of points in H/H2 ratio grid (RH)
c  npr : Nb of points in Ortho/Para ratio grid (ROP)

       parameter (npt=23)
       parameter (npd=49)
       parameter (nph=36)
       parameter (npr=6)

       real*8 TN(npt), DN(npd)
       real*8 RH(nph), ROP(npr)

       real*8 WG(npt,npd,nph,npr)

       common / grid / TN, DN, RH, ROP, WG

       data TNmin, TNmax / 100.0d0, 1.0d+04 /
       data DNmin, DNmax / 1.0d0, 1.0d+08 /

c  Read Cooling file (binary)

       open (10, file="le_cube", status="old")

       read(10,*) RH
       read(10,*) ROP
       read(10,*) TN
       read(10,*) DN
       read(10,*) ((((WG(it,id,ih,io), id=1,npd)
     1          , it=1,npt), io=1,npr), ih=1,nph)

       close (10)

       return
       end
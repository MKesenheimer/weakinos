      subroutine fsirgx(fi,sc,gc,fmass,fwidth , fsi)
c     
c- by Yoshitaro Takaesu - 2010/12/28
c
      implicit none
      double complex fi(6),sc(6),fsi(6),gc(2),gcc(2)
      double precision fmass,fwidth

      gcc(1) = dconjg(gc(2))
      gcc(2) = dconjg(gc(1))

      call fsixxx(fi,sc,gcc,fmass,fwidth , fsi)

      return
      end

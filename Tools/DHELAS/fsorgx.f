      subroutine fsorgx(fo,sc,gc,fmass,fwidth , fso)
c     
c- by Yoshitaro Takaesu - 2010/10/31
c
      implicit none
      double complex fo(6),sc(6),fso(6),gc(2)
      double precision fmass,fwidth

      call fsoxxx(fo,sc,gc,fmass,fwidth , fso)

      return
      end

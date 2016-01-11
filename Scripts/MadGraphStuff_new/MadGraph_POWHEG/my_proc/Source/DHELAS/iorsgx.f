      subroutine iorsgx(fi,fo,sc,gc , vertex)
c     
c- by Yoshitaro Takaesu - 2010/10/31
c
      implicit none
      double complex fi(6),fo(6),sc(3),gc(2),vertex

      call  iosxxx(fi,fo,sc,gc , vertex)

      return
      end

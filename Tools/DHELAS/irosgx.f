      subroutine irosgx(ri,fo,sc,gc , vertex)
c     
c- by Yoshitaro Takaesu - 2010/12/28
c
      implicit none
      double complex ri(6),fo(6),sc(3),gc(2),gcc(2),vertex

      gcc(1) = dconjg(gc(2))
      gcc(2) = dconjg(gc(1))

      call  iosxxx(ri,fo,sc,gcc , vertex)

      return
      end

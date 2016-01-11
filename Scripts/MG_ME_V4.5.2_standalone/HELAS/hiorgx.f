      subroutine hiorgx(fi,fo,gc,smass,swidth , hio)
c     
c- by Yoshitaro Takaesu - 2010/10/31
c
      implicit none
      double complex fi(6),fo(6),hio(3),gc(2)
      double precision smass,swidth

      call hioxxx(fi,fo,gc,smass,swidth , hio)

      return
      end

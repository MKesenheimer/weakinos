      subroutine hirogx(fi,fo,gc,smass,swidth , hio)
c     
c- by Yoshitaro Takaesu - 2010/12/28
c
      implicit none
      double complex fi(6),fo(6),hio(3),gc(2),smass,swidth,gcc(2)
      
      gcc(1) = dconjg(gc(2))
      gcc(2) = dconjg(gc(1))

      call hioxxx(fi,fo,gcc,smass,swidth , hio)

      return
      end

      subroutine born_phsp(xborn)
      implicit none
      include 'nlegborn.h'
      real * 8 xborn(ndiminteg-3)
      end



      subroutine born_suppression(fact)
      implicit none
      real * 8 fact
      fact=1d0
      end


      subroutine set_fac_ren_scales(muf,mur)
      implicit none
      real * 8 muf,mur
      muf=100d0
      mur=100d0
      end

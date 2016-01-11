      subroutine setvirtual(p,vflav,virtual)
c Virtual needs to be provided by the user and put here
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      real * 8 p(0:3,nlegborn)
      integer vflav(nlegborn)
      real * 8 virtual

      virtual = 0d0

      end

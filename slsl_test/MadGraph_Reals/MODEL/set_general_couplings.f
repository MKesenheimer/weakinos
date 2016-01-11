c############### set_general_couplings.f ##################################
c last modified by MK, date

c############### subroutine set_general_couplings #########################
      subroutine set_general_couplings(fin1,fin2)
        implicit none
   
#include "finalstate.h"

        include 'coupl.inc'
        integer fin1,fin2

        final1 = fin1
        final2 = fin2

      end

c############### end subroutine set_general_couplings #####################

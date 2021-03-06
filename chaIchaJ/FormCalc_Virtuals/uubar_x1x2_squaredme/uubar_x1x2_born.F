#if 0
* uubar_x1x2_born.F
* form factors for {F[3, {1, Col1}], -F[3, {1, Col2}]} -> {F[12, {1}], -F[12, {2}]}
* generated by FormCalc 8.4 on 17-Jun-2015 12:36
#endif

#include "uubar_x1x2_vars.h"
#include "contains.h"


	subroutine uubar_x1x2_born
	implicit none

#include "uubar_x1x2_vars.h"


        Ctree(HelInd(1)) = 
     &   Ctree(HelInd(1)) + (2*Alfa*Pi)/(3.D0*CW2)*Sub3/
     &                  (MZ2 - S -cI*MZ*WZ)

#ifdef DEBUGQQ
	DEB("uubar_x1x2_born: Ctree(1) =", Ctree(HelInd(1)))
#endif

#include "contains.h"

	end

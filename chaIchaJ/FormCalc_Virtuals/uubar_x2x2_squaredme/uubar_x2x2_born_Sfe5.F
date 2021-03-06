#if 0
* uubar_x2x2_born_Sfe5.F
* form factors for {F[3, {1, Col1}], -F[3, {1, Col2}]} -> {F[12, {2}], -F[12, {2}]}
* generated by FormCalc 8.4 on 18-Jun-2015 10:15
#endif

#include "uubar_x2x2_vars.h"
#include "contains.h"


	subroutine uubar_x2x2_born_Sfe5
	implicit none

#include "uubar_x2x2_vars.h"


        Ctree(HelInd(1)) = 
     &   Ctree(HelInd(1)) - 
     &    (Alfa*Pi)/(MW2*SW2)*
     &     (F10*F9*Sub1(Sfe5)*Sub2(Sfe5))/
     &      (CB2*(U - MSf2(Sfe5,4,1)))

#ifdef DEBUGQQ
	DEB("uubar_x2x2_born_Sfe5: Ctree(1) =", Ctree(HelInd(1)))
#endif

#include "contains.h"

	end

#if 0
* ddbar_x1x1_box_Sfe5_Sfe6.F
* form factors for {F[4, {1, Col1}], -F[4, {1, Col2}]} -> {F[12, {1}], -F[12, {1}]}
* generated by FormCalc 8.4 on 12-Jun-2015 9:32
#endif

#include "ddbar_x1x1_vars.h"
#include "contains.h"


	subroutine ddbar_x1x1_box_Sfe5_Sfe6
	implicit none

#include "ddbar_x1x1_vars.h"


        Cloop(HelInd(1)) = 
     &   Cloop(HelInd(1)) + 
     &    Alfa/(MW2*SW2)*(Alfas*
     &        (2/3.D0*Sub122(HelInd(Sfe5),Sfe6) - 
     &          4/3.D0*(pave4(PVC(dd00),Sfe5,Sfe6)*
     &             Sub39(HelInd(Sfe5),Sfe6)) - 
     &          2/3.D0*(Sub123(HelInd(Sfe5),Sfe6)*
     &             Sub96(HelInd(Sfe5),Sfe6)))*Sub104(Sfe5)*
     &        Sub43(Sfe6))/CB2

#ifdef DEBUGQQ
	DEB("ddbar_x1x1_box_Sfe5_Sfe6: Cloop(1) =", Cloop(HelInd(1)))
#endif

#include "contains.h"

	end

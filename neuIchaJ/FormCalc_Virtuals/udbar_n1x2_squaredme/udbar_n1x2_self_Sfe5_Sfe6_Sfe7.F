#if 0
* udbar_n1x2_self_Sfe5_Sfe6_Sfe7.F
* form factors for {F[3, {1, Col1}], -F[4, {1, Col2}]} -> {F[11, {1}], -F[12, {2}]}
* generated by FormCalc 8.4 on 19-Jun-2015 11:51
#endif

#include "udbar_n1x2_vars.h"
#include "contains.h"


	subroutine udbar_n1x2_self_Sfe5_Sfe6_Sfe7
	implicit none

#include "udbar_n1x2_vars.h"


        Cloop(HelInd(1)) = 
     &   Cloop(HelInd(1)) + 
     &    1/(CW*MW2*SW2)*(1/432.D0*
     &        (F3*pave1(Sfe7)*Sub20(Sfe6)*Sub47(Sfe5,Sfe6,Sfe7)*
     &           (F5*Sub4(Sfe5) - F4*Sub6(Sfe5)))/
     &         (SB2*(T - MSf2(Sfe5,3,1))*(T - MSf2(Sfe6,3,1))) - 
     &       1/432.D0*(F6*pave2(Sfe7)*
     &           (F8*Sub29(Sfe6) + F7*Sub30(Sfe6))*
     &           Sub55(Sfe5,Sfe6,Sfe7)*Sub7(Sfe5))/
     &         (CB2*(U - MSf2(Sfe5,4,1))*(U - MSf2(Sfe6,4,1))))

#ifdef DEBUGQQ
	DEB("udbar_n1x2_self_Sfe5_Sfe6_Sfe7: Cloop(1) =", Cloop(HelInd(1)))
#endif

#include "contains.h"

	end

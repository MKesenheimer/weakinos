#if 0
* ddbar_n1n3_self_Sfe5_Sfe6.F
* form factors for {F[4, {1, Col1}], -F[4, {1, Col2}]} -> {F[11, {1}], F[11, {3}]}
* generated by FormCalc 8.4 on 7-Jul-2015 9:13
#endif

#include "ddbar_n1n3_vars.h"
#include "contains.h"


	subroutine ddbar_n1n3_self_Sfe5_Sfe6
	implicit none

#include "ddbar_n1n3_vars.h"


        Cloop(HelInd(1)) = 
     &   Cloop(HelInd(1)) + 
     &    Alfa/(CW2*MW2*SW2)*
     &     (1/27.D0*(Sub255(Sfe5,Sfe6)*
     &            (F7*Sub30(Sfe6) + F2*Sub31(Sfe6))*
     &            (F3*Sub4(Sfe5) + F8*Sub6(Sfe5)))/
     &          ((U - MSf2(Sfe5,4,1))*(U - MSf2(Sfe6,4,1))) - 
     &        1/27.D0*((F9*Sub21(Sfe6) + F11*Sub22(Sfe6))*
     &            Sub256(Sfe5,Sfe6)*
     &            (F10*Sub7(Sfe5) + F12*Sub9(Sfe5)))/
     &          ((T - MSf2(Sfe5,4,1))*(T - MSf2(Sfe6,4,1))))/CB2

#ifdef DEBUGQ
	DEB("ddbar_n1n3_self_Sfe5_Sfe6: Cloop(1) =", Cloop(HelInd(1)))
#endif

#include "contains.h"

	end

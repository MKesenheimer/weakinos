#if 0
* ddbar_n3n3_abbr0s.F
* abbreviations for {F[4, {1, Col1}], -F[4, {1, Col2}]} -> {F[11, {3}], F[11, {3}]}
* generated by FormCalc 8.4 on 7-Jul-2015 9:52
#endif

#include "ddbar_n3n3_vars.h"
#include "contains.h"


	subroutine ddbar_n3n3_abbr0s
	implicit none

#include "ddbar_n3n3_vars.h"


	Sub4 = -(3*CW*ZNeu(3,2)) + SW*ZNeu(3,1)

	Sub2 = -(3*CW*ZNeuC(3,2)) + SW*ZNeuC(3,1)

	Sub8 = ZNeu(3,3)*ZNeuC(3,3) - ZNeu(3,4)*ZNeuC(3,4)

	MatSUN(1,1) = 3

	LOOP(Sfe5, 1,2,1)

        Sub1(Sfe5) = 2*MW*(CB*SW*USf(Sfe5,2,4,1)*ZNeu(3,1)) + 
     &    3*CW*MD*(USf(Sfe5,1,4,1)*ZNeu(3,3))

        Sub5(Sfe5) = MW*(CB*Sub4*USfC(Sfe5,1,4,1)) + 
     &    3*CW*MD*(USfC(Sfe5,2,4,1)*ZNeu(3,3))

        Sub3(Sfe5) = MW*(CB*Sub2*USf(Sfe5,1,4,1)) + 
     &    3*CW*MD*(USf(Sfe5,2,4,1)*ZNeuC(3,3))

        Sub6(Sfe5) = 2*MW*
     &     (CB*SW*USfC(Sfe5,2,4,1)*ZNeuC(3,1)) + 
     &    3*CW*MD*(USfC(Sfe5,1,4,1)*ZNeuC(3,3))

	ENDLOOP(Sfe5)

#include "contains.h"

	end

#if 0
* ddbar_n1n4_abbr1s.F
* abbreviations for {F[4, {1, Col1}], -F[4, {1, Col2}]} -> {F[11, {1}], F[11, {4}]}
* generated by FormCalc 8.4 on 7-Jul-2015 9:33
#endif

#include "ddbar_n1n4_vars.h"
#include "contains.h"


	subroutine ddbar_n1n4_abbr1s
	implicit none

#include "ddbar_n1n4_vars.h"


	Sub84 = 2*dCB1 - CB*(2*dZe1 + dZfL1(4,1,1))

	Sub85 = CB*(dMWsq1*SW + 2*MW2*dSW1) + MW2*(Sub84*SW)

	Sub98 = 2*dCB1 - CB*(2*dZe1 + dZbarfR1(4,1,1))

	Sub99 = CB*(dMWsq1*SW + 2*MW2*dSW1) + MW2*(Sub98*SW)

	Sub86 = MD*Sub85 - 2*MW2*(CB*SW*dMf1(4,1))

	Sub68 = 2*dCB1 - CB*(2*dZe1 + dZbarfL1(4,1,1))

	Sub81 = 2*dSW1 + SW*(-2*dZe1 - dZfL1(4,1,1))

	Sub61 = -(3*CW*ZNeu(2,2)) + SW*ZNeu(2,1)

	Sub62 = -(3*CW*ZNeu(3,2)) + SW*ZNeu(3,1)

	Sub89 = -(3*CW*ZNeuC(2,2)) + SW*ZNeuC(2,1)

	Sub90 = -(3*CW*ZNeuC(3,2)) + SW*ZNeuC(3,1)

        Sub132 = Sub8*dZfL1(11,1,4) + Sub89*dZfL1(11,2,4) + 
     &    Sub90*dZfL1(11,3,4) + Sub5*dZfL1(11,4,4)

        Sub102 = dZfL1(11,1,4)*ZNeuC(1,1) + 
     &    dZfL1(11,2,4)*ZNeuC(2,1) + dZfL1(11,3,4)*ZNeuC(3,1) + 
     &    dZfL1(11,4,4)*ZNeuC(4,1)

        Sub103 = dZfL1(11,1,4)*ZNeuC(1,3) + 
     &    dZfL1(11,2,4)*ZNeuC(2,3) + dZfL1(11,3,4)*ZNeuC(3,3) + 
     &    dZfL1(11,4,4)*ZNeuC(4,3)

	Sub55 = 2*MW2*(SW*dMf1(4,1)) - MD*(dMWsq1*SW + 2*MW2*dSW1)

	Sub100 = MD*Sub99 - 2*MW2*(CB*SW*dMf1(4,1))

	Sub60 = 2*dSW1 + SW*(-2*dZe1 - dZbarfL1(4,1,1))

        Sub52 = 2*MW*(CB*SW*USf(1,2,4,1)*ZNeu(1,1)) + 
     &    3*CW*MD*(USf(1,1,4,1)*ZNeu(1,3))

        Sub53 = 2*MW*(CB*SW*USf(2,2,4,1)*ZNeu(1,1)) + 
     &    3*CW*MD*(USf(2,1,4,1)*ZNeu(1,3))

        Sub112 = 2*MW*(CB*SW*USf(1,2,4,1)*ZNeu(4,1)) + 
     &    3*CW*MD*(USf(1,1,4,1)*ZNeu(4,3))

        Sub113 = 2*MW*(CB*SW*USf(2,2,4,1)*ZNeu(4,1)) + 
     &    3*CW*MD*(USf(2,1,4,1)*ZNeu(4,3))

        Sub91 = Sub8*dZfL1(11,1,1) + Sub89*dZfL1(11,2,1) + 
     &    Sub90*dZfL1(11,3,1) + Sub5*dZfL1(11,4,1)

        Sub116 = Sub1*dZfR1(11,1,1) + Sub61*dZfR1(11,2,1) + 
     &    Sub62*dZfR1(11,3,1) + Sub10*dZfR1(11,4,1)

        Sub63 = Sub1*dZfR1(11,1,4) + Sub61*dZfR1(11,2,4) + 
     &    Sub62*dZfR1(11,3,4) + Sub10*dZfR1(11,4,4)

        Sub56 = dZfR1(11,1,1)*ZNeu(1,3) + 
     &    dZfR1(11,2,1)*ZNeu(2,3) + dZfR1(11,3,1)*ZNeu(3,3) + 
     &    dZfR1(11,4,1)*ZNeu(4,3)

        Sub69 = dZfR1(11,1,4)*ZNeu(1,3) + 
     &    dZfR1(11,2,4)*ZNeu(2,3) + dZfR1(11,3,4)*ZNeu(3,3) + 
     &    dZfR1(11,4,4)*ZNeu(4,3)

        Sub127 = dZfL1(11,1,1)*ZNeuC(1,1) + 
     &    dZfL1(11,2,1)*ZNeuC(2,1) + dZfL1(11,3,1)*ZNeuC(3,1) + 
     &    dZfL1(11,4,1)*ZNeuC(4,1)

        Sub92 = dZfL1(11,1,1)*ZNeuC(1,3) + 
     &    dZfL1(11,2,1)*ZNeuC(2,3) + dZfL1(11,3,1)*ZNeuC(3,3) + 
     &    dZfL1(11,4,1)*ZNeuC(4,3)

	Sub45 = 4*CW2**2*dZe1 + CW2*(1 + 2*CW2)*dZZZ1

	Sub59 = 2*dSW1*SW + CW2*(2*dZe1 + dZbarfL1(4,1,1))

	Sub96 = 2*dSW1*SW + CW2*(2*dZe1 + dZbarfR1(4,1,1))

	Sub79 = 2*dSW1*SW + CW2*(2*dZe1 + dZfL1(4,1,1))

	Sub48 = 2*dCB1 - CB*(2*dZe1 + dZfR1(4,1,1))

	Sub49 = 2*dSW1*SW + CW2*(2*dZe1 + dZfR1(4,1,1))

	Sub117 = Sub116*SW + 3*CW*(Sub60*ZNeu(1,2))

	Sub57 = MD*MW2*(Sub56*SW) + Sub55*ZNeu(1,3)

	Sub64 = Sub63*SW + 3*CW*(Sub60*ZNeu(4,2))

	Sub111 = MD*MW2*(Sub69*SW) + Sub55*ZNeu(4,3)

        Sub50 = dZfR1(11,1,1)*ZNeu(1,1) + 
     &    dZfR1(11,2,1)*ZNeu(2,1) + dZfR1(11,3,1)*ZNeu(3,1) + 
     &    dZfR1(11,4,1)*ZNeu(4,1)

        Sub109 = dZfR1(11,1,4)*ZNeu(1,1) + 
     &    dZfR1(11,2,4)*ZNeu(2,1) + dZfR1(11,3,4)*ZNeu(3,1) + 
     &    dZfR1(11,4,4)*ZNeu(4,1)

	Sub46 = Sub45*SW + 2*SW2*dSW1 - CW2*(6*dSW1 - 2*dZe1*SW)

	Sub43 = 2*dZe1 + dZZZ1 + dZbarfR1(4,1,1) + dZfR1(4,1,1)

	Sub44 = CW2**2*dZAZ1 - CW*(2*dSW1 + CW2*(Sub43*SW))

        Sub47 = 2*CW2**2*SW2*dZAZ1 + 
     &    CW*(Sub46 + CW2*(1 + 2*CW2)*
     &        (SW*(dZbarfL1(4,1,1) + dZfL1(4,1,1))))

        Sub168 = 1/(CW*CW2)*
     &    (3*Pi*Sub44 + 2*CW*CW2*(Alfas*Finite*SW))/SW

        Sub169 = 1/(CW*CW2)*
     &    (3*Pi*Sub47 - CW*CW2*(2 + 4*CW2)*(Alfas*Finite*SW))/SW

	pave23 = A0(MGl2)

	call Bput(pave2(bb0),S,MD2,MD2)

	call Cput(pave3(cc0),MD2,S,MD2,0.D0,MD2,MD2)

	Sub119 = MD*Sub56 + 2*dMf1(4,1)*ZNeu(1,3)

	Sub70 = MD*Sub69 + 2*dMf1(4,1)*ZNeu(4,3)

	Sub120 = Sub119*SW - 2*MD*(dSW1*ZNeu(1,3))

	Sub71 = Sub70*SW - 2*MD*(dSW1*ZNeu(4,3))

	Sub121 = MW2*Sub120 - MD*(dMWsq1*SW*ZNeu(1,3))

	Sub72 = MW2*Sub71 - MD*(dMWsq1*SW*ZNeu(4,3))

	Sub122 = CB*Sub121 - MD*MW2*(Sub68*SW*ZNeu(1,3))

	Sub73 = CB*Sub72 - MD*MW2*(Sub68*SW*ZNeu(4,3))

	Sub275 = pave3(PVC(cc0)) + pave3(PVC(cc1)) + pave3(PVC(cc2))

	LOOP(Sfe5, 1,2,1)

        Sub83(Sfe5) = dZbarSf1(1,Sfe5,4,1)*USf(1,2,4,1) + 
     &    dZbarSf1(2,Sfe5,4,1)*USf(2,2,4,1)

        Sub133(Sfe5) = 
     &   MW*(CB*Sub132*USf(Sfe5,1,4,1)) + 
     &    3*CW*MD*(Sub103*USf(Sfe5,2,4,1))

        Sub87(Sfe5) = MD*MW2*(CB*SW*Sub83(Sfe5)) - 
     &    Sub86*USf(Sfe5,2,4,1)

        Sub65(Sfe5) = dZSf1(1,Sfe5,4,1)*USfC(1,1,4,1) + 
     &    dZSf1(2,Sfe5,4,1)*USfC(2,1,4,1)

        Sub104(Sfe5) = 
     &   3*CW*MD*(Sub103*USfC(Sfe5,1,4,1)) + 
     &    2*MW*(CB*Sub102*SW*USfC(Sfe5,2,4,1))

        Sub114(Sfe5) = 
     &   Sub112*dZbarSf1(1,Sfe5,4,1) + Sub113*dZbarSf1(2,Sfe5,4,1)

        Sub54(Sfe5) = Sub52*dZbarSf1(1,Sfe5,4,1) + 
     &    Sub53*dZbarSf1(2,Sfe5,4,1)

        Sub78(Sfe5) = dZbarSf1(1,Sfe5,4,1)*USf(1,1,4,1) + 
     &    dZbarSf1(2,Sfe5,4,1)*USf(2,1,4,1)

        Sub93(Sfe5) = MW*(CB*Sub91*USf(Sfe5,1,4,1)) + 
     &    3*CW*MD*(Sub92*USf(Sfe5,2,4,1))

        Sub66(Sfe5) = dZSf1(1,Sfe5,4,1)*USfC(1,2,4,1) + 
     &    dZSf1(2,Sfe5,4,1)*USfC(2,2,4,1)

        Sub101(Sfe5) = 
     &   MD*MW2*(CB*SW*Sub65(Sfe5)) - Sub100*USfC(Sfe5,1,4,1)

        Sub128(Sfe5) = 
     &   3*CW*MD*(Sub92*USfC(Sfe5,1,4,1)) + 
     &    2*MW*(CB*Sub127*SW*USfC(Sfe5,2,4,1))

        Sub105(Sfe5) = 
     &   MW2*(CB*SW*Sub104(Sfe5)) + 3*CW*(Sub101(Sfe5)*ZNeuC(4,3))

        Sub115(Sfe5) = 
     &   MW2*(SW*Sub114(Sfe5)) + 3*CW*(Sub111*USf(Sfe5,1,4,1))

        Sub58(Sfe5) = MW2*(SW*Sub54(Sfe5)) + 
     &    3*CW*(Sub57*USf(Sfe5,1,4,1))

	Sub80(Sfe5) = CW2*Sub78(Sfe5) + Sub79*USf(Sfe5,1,4,1)

	Sub97(Sfe5) = CW2*Sub66(Sfe5) + Sub96*USfC(Sfe5,2,4,1)

        Sub51(Sfe5) = 2*MW*SW2*
     &     (CB2*USf(Sfe5,2,4,1)*(CW2*Sub50 + Sub49*ZNeu(1,1))) - 
     &    3*CW*CW2*MD*(Sub48*SW*USf(Sfe5,1,4,1)*ZNeu(1,3))

        Sub110(Sfe5) = 
     &   2*MW*SW2*(CB2*USf(Sfe5,2,4,1)*
     &       (CW2*Sub109 + Sub49*ZNeu(4,1))) - 
     &    3*CW*CW2*MD*(Sub48*SW*USf(Sfe5,1,4,1)*ZNeu(4,3))

        Sub129(Sfe5) = 
     &   MW2*(CB*SW*Sub128(Sfe5)) + 3*CW*(Sub101(Sfe5)*ZNeuC(1,3))

        Sub106(Sfe5) = 
     &   CW2*Sub105(Sfe5) + 
     &    2*MW*MW2*SW2*(CB2*Sub97(Sfe5)*ZNeuC(4,1))

        Sub130(Sfe5) = 
     &   CW2*Sub129(Sfe5) + 
     &    2*MW*MW2*SW2*(CB2*Sub97(Sfe5)*ZNeuC(1,1))

        call Cput(pave12(cc0,Sfe5),S,MNeu2(1),MNeu2(4),MD2,MD2,
     &   MSf2(Sfe5,4,1))

        call Cput(pave13(cc0,Sfe5),S,MNeu2(4),MNeu2(1),MD2,MD2,
     &   MSf2(Sfe5,4,1))

	Sub82(Sfe5) = SW*Sub78(Sfe5) - Sub81*USf(Sfe5,1,4,1)

        Sub118(Sfe5) = 
     &   MW*(CB*Sub1*Sub65(Sfe5)) + 3*CW*MD*(Sub66(Sfe5)*ZNeu(1,3))

        Sub67(Sfe5) = MW*(CB*Sub10*Sub65(Sfe5)) + 
     &    3*CW*MD*(Sub66(Sfe5)*ZNeu(4,3))

        Sub131(Sfe5) = 
     &   MW*MW2*(CB2*Sub82(Sfe5)*ZNeuC(4,2)) - 
     &    Sub87(Sfe5)*ZNeuC(4,3)

	Sub134(Sfe5) = 3*CW*Sub131(Sfe5) - MW2*(CB*SW*Sub133(Sfe5))

        Sub123(Sfe5) = 
     &   MW2*(CB*SW*Sub118(Sfe5)) + 3*CW*(Sub122*USfC(Sfe5,2,4,1))

        Sub74(Sfe5) = MW2*(CB*SW*Sub67(Sfe5)) + 
     &    3*CW*(Sub73*USfC(Sfe5,2,4,1))

        Sub88(Sfe5) = MW*MW2*(CB2*Sub82(Sfe5)*ZNeuC(1,2)) - 
     &    Sub87(Sfe5)*ZNeuC(1,3)

	Sub94(Sfe5) = 3*CW*Sub88(Sfe5) - MW2*(CB*SW*Sub93(Sfe5))

        Sub124(Sfe5) = 
     &   CW2*Sub123(Sfe5) + 
     &    MW*MW2*(CB2*(CW2*Sub117 + SW2*(Sub59*ZNeu(1,1)))*
     &       USfC(Sfe5,1,4,1))

        Sub75(Sfe5) = CW2*Sub74(Sfe5) + 
     &    MW*MW2*(CB2*(CW2*Sub64 + SW2*(Sub59*ZNeu(4,1)))*
     &       USfC(Sfe5,1,4,1))

        Sub135(Sfe5) = 
     &   CW2*Sub134(Sfe5) - MW*MW2*SW2*(CB2*Sub80(Sfe5)*ZNeuC(4,1))

        Sub95(Sfe5) = CW2*Sub94(Sfe5) - 
     &    MW*MW2*SW2*(CB2*Sub80(Sfe5)*ZNeuC(1,1))

        Sub125(Sfe5) = 
     &   (MW2*Sub110(Sfe5) + CW2*(CB*Sub115(Sfe5)))*Sub2(Sfe5) + 
     &    Sub124(Sfe5)*Sub4(Sfe5)

        Sub76(Sfe5) = (MW2*Sub51(Sfe5) + CW2*(CB*Sub58(Sfe5)))*
     &     Sub11(Sfe5) + Sub7(Sfe5)*Sub75(Sfe5)

	ENDLOOP(Sfe5)

	LOOP(Sfe6, 1,2,1)

        Sub151(Sfe6) = 
     &   2*MW*(CB*SW*USf(Sfe6,2,4,1)*ZNeu(1,1)) + 
     &    3*CW*MD*(USf(Sfe6,1,4,1)*ZNeu(1,3))

        Sub143(Sfe6) = 
     &   2*MW*(CB*SW*USf(Sfe6,2,4,1)*ZNeu(4,1)) + 
     &    3*CW*MD*(USf(Sfe6,1,4,1)*ZNeu(4,3))

        Sub153(Sfe6) = 
     &   MW*(CB*Sub8*USf(Sfe6,1,4,1)) + 
     &    3*CW*MD*(USf(Sfe6,2,4,1)*ZNeuC(1,3))

        Sub144(Sfe6) = 
     &   MW*(CB*Sub5*USf(Sfe6,1,4,1)) + 
     &    3*CW*MD*(USf(Sfe6,2,4,1)*ZNeuC(4,3))

        Sub30(Sfe6) = MW*(CB*Sub1*USfC(Sfe6,1,4,1)) + 
     &    3*CW*MD*(USfC(Sfe6,2,4,1)*ZNeu(1,3))

        Sub21(Sfe6) = MW*(CB*Sub10*USfC(Sfe6,1,4,1)) + 
     &    3*CW*MD*(USfC(Sfe6,2,4,1)*ZNeu(4,3))

        Sub31(Sfe6) = 2*MW*
     &     (CB*SW*USfC(Sfe6,2,4,1)*ZNeuC(1,1)) + 
     &    3*CW*MD*(USfC(Sfe6,1,4,1)*ZNeuC(1,3))

        Sub22(Sfe6) = 2*MW*
     &     (CB*SW*USfC(Sfe6,2,4,1)*ZNeuC(4,1)) + 
     &    3*CW*MD*(USfC(Sfe6,1,4,1)*ZNeuC(4,3))

	ENDLOOP(Sfe6)

	LOOP(Sfe7, 1,2,1)

	pave1(Sfe7) = A0(MSf2(Sfe7,4,1))

	ENDLOOP(Sfe7)

	LOOP(Sfe6, 1,2,1)
	LOOP(Sfe5, 1,2,1)

        Sub16(Sfe5,Sfe6) = 
     &   dZbarSf1(Sfe5,Sfe6,4,1) + dZSf1(Sfe6,Sfe5,4,1)

        Sub17(Sfe5,Sfe6) = 
     &   2*dMSfsq1(Sfe6,Sfe5,4,1) + 
     &    dZbarSf1(Sfe5,Sfe6,4,1)*MSf2(Sfe5,4,1) + 
     &    dZSf1(Sfe6,Sfe5,4,1)*MSf2(Sfe6,4,1)

        Sub32(Sfe5,Sfe6) = 
     &   USf(Sfe6,1,4,1)*USfC(Sfe5,1,4,1) + 
     &    USf(Sfe6,2,4,1)*USfC(Sfe5,2,4,1)

        call Cput(pave20(cc0,Sfe5,Sfe6),MD2,S,MD2,MGl2,
     &   MSf2(Sfe5,4,1),MSf2(Sfe6,4,1))

        Sub138(Sfe5,Sfe6) = 
     &   -(2*SW2*(USf(Sfe6,2,4,1)*USfC(Sfe5,2,4,1))) + 
     &    USf(Sfe6,1,4,1)*(2*CW2*USfC(Sfe5,1,4,1) + 
     &       USfC(Sfe5,1,4,1))

        Sub152(Sfe5,Sfe6) = 
     &   Sub2(Sfe5)*Sub21(Sfe6)*USf(Sfe5,1,4,1)*USf(Sfe6,2,4,1) + 
     &    Sub151(Sfe6)*Sub4(Sfe5)*USfC(Sfe5,2,4,1)*USfC(Sfe6,1,4,1)

        Sub146(Sfe5,Sfe6) = 
     &   Sub11(Sfe5)*Sub30(Sfe6)*USf(Sfe5,1,4,1)*USf(Sfe6,2,4,1) + 
     &    Sub143(Sfe6)*Sub7(Sfe5)*USfC(Sfe5,2,4,1)*USfC(Sfe6,1,4,1)

        Sub147(Sfe5,Sfe6) = 
     &   Sub11(Sfe5)*Sub31(Sfe6)*USf(Sfe5,2,4,1)*USf(Sfe6,1,4,1) + 
     &    Sub143(Sfe6)*Sub9(Sfe5)*USfC(Sfe5,2,4,1)*USfC(Sfe6,1,4,1)

        Sub154(Sfe5,Sfe6) = 
     &   Sub21(Sfe6)*Sub3(Sfe5)*USf(Sfe5,1,4,1)*USf(Sfe6,2,4,1) + 
     &    Sub153(Sfe6)*Sub4(Sfe5)*USfC(Sfe5,1,4,1)*USfC(Sfe6,2,4,1)

	ENDLOOP(Sfe5)
	ENDLOOP(Sfe6)

	LOOP(Sfe7, 1,2,1)
	LOOP(Sfe5, 1,2,1)

        Sub35(Sfe5,Sfe7) = 
     &   USf(Sfe7,1,4,1)*USfC(Sfe5,1,4,1) - 
     &    USf(Sfe7,2,4,1)*USfC(Sfe5,2,4,1)

	ENDLOOP(Sfe5)
	ENDLOOP(Sfe7)

	LOOP(Sfe7, 1,2,1)
	LOOP(Sfe6, 1,2,1)

        Sub36(Sfe6,Sfe7) = 
     &   USf(Sfe6,1,4,1)*USfC(Sfe7,1,4,1) - 
     &    USf(Sfe6,2,4,1)*USfC(Sfe7,2,4,1)

	ENDLOOP(Sfe6)
	ENDLOOP(Sfe7)

	LOOP(Sfe7, 1,2,1)
	LOOP(Sfe6, 1,2,1)
	LOOP(Sfe5, 1,2,1)

        Sub42(Sfe5,Sfe6,Sfe7) = 
     &   48*Alfa*(Alfas*Sub35(Sfe5,Sfe7)*Sub36(Sfe6,Sfe7))

	ENDLOOP(Sfe5)
	ENDLOOP(Sfe6)
	ENDLOOP(Sfe7)

#include "contains.h"

	end

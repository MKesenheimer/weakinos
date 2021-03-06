#if 0
* dubar_n1x1_abbr1s.F
* abbreviations for {F[4, {1, Col1}], -F[3, {1, Col2}]} -> {F[11, {1}], F[12, {1}]}
* generated by FormCalc 8.4 on 19-Jun-2015 12:29
#endif

#include "dubar_n1x1_vars.h"
#include "contains.h"


	subroutine dubar_n1x1_abbr1s
	implicit none

#include "dubar_n1x1_vars.h"


	Sub104 = 2*dCB1 - CB*(2*dZe1 + dZfL1(4,1,1))

	Sub105 = CB*(dMWsq1*SW + 2*MW2*dSW1) + MW2*(Sub104*SW)

	Sub106 = MD*Sub105 - 2*MW2*(CB*SW*dMf1(4,1))

	Sub109 = -(3*CW*ZNeuC(2,2)) + SW*ZNeuC(2,1)

	Sub110 = -(3*CW*ZNeuC(3,2)) + SW*ZNeuC(3,1)

	Sub111 = -(3*CW*ZNeuC(4,2)) + SW*ZNeuC(4,1)

        Sub150 = dZfR1(11,1,1)*ZNeu(1,4) + 
     &    dZfR1(11,2,1)*ZNeu(2,4) + dZfR1(11,3,1)*ZNeu(3,4) + 
     &    dZfR1(11,4,1)*ZNeu(4,4)

	Sub75 = 2*dCB1 - CB*(2*dZe1 + dZbarfL1(3,1,1))

	Sub131 = dMWsq1 + MW2*(-2*dZe1 - dZfL1(4,1,1))

	Sub151 = MU*Sub150 + 2*dMf1(3,1)*ZNeu(1,4)

        Sub112 = Sub5*dZfL1(11,1,1) + Sub109*dZfL1(11,2,1) + 
     &    Sub110*dZfL1(11,3,1) + Sub111*dZfL1(11,4,1)

        Sub113 = dZfL1(11,1,1)*ZNeuC(1,3) + 
     &    dZfL1(11,2,1)*ZNeuC(2,3) + dZfL1(11,3,1)*ZNeuC(3,3) + 
     &    dZfL1(11,4,1)*ZNeuC(4,3)

	Sub76 = CB*(dMWsq1*SW + 2*MW2*dSW1) + MW2*(Sub75*SW)

	Sub132 = SB*Sub131*SW + 2*MW2*(dSW1*SB + dSB1*SW)

	Sub94 = 2*MW2*(SW*dMf1(4,1)) - MD*(dMWsq1*SW + 2*MW2*dSW1)

	Sub81 = dMWsq1 + MW2*(-2*dZe1 - dZbarfR1(3,1,1))

	Sub99 = 2*dSW1*SW + CW2*(2*dZe1 + dZfL1(4,1,1))

	Sub101 = 2*dSW1 + SW*(-2*dZe1 - dZfL1(4,1,1))

        Sub136 = dZbarfR1(12,1,1)*VChaC(1,1) + 
     &    dZbarfR1(12,1,2)*VChaC(2,1)

        Sub91 = 2*MW*(CB*SW*USf(1,2,4,1)*ZNeu(1,1)) + 
     &    3*CW*MD*(USf(1,1,4,1)*ZNeu(1,3))

        Sub92 = 2*MW*(CB*SW*USf(2,2,4,1)*ZNeu(1,1)) + 
     &    3*CW*MD*(USf(2,1,4,1)*ZNeu(1,3))

	Sub152 = Sub151*SW - 2*MU*(dSW1*ZNeu(1,4))

        Sub95 = dZfR1(11,1,1)*ZNeu(1,3) + 
     &    dZfR1(11,2,1)*ZNeu(2,3) + dZfR1(11,3,1)*ZNeu(3,3) + 
     &    dZfR1(11,4,1)*ZNeu(4,3)

	Sub82 = SB*Sub81*SW + 2*MW2*(dSW1*SB + dSB1*SW)

	Sub133 = MU*Sub132 - 2*MW2*(SB*SW*dMf1(3,1))

	Sub77 = MD*Sub76 - 2*MW2*(CB*SW*dMf1(4,1))

	Sub149 = dMWsq1 + MW2*(-2*dZe1 - dZbarfL1(3,1,1))

	Sub87 = 2*dCB1 - CB*(2*dZe1 + dZfR1(4,1,1))

	Sub88 = 2*dSW1*SW + CW2*(2*dZe1 + dZfR1(4,1,1))

        Sub73 = dZbarfL1(12,1,1)*UCha(1,1) + 
     &    dZbarfL1(12,1,2)*UCha(2,1)

        Sub78 = dZbarfL1(12,1,1)*UCha(1,2) + 
     &    dZbarfL1(12,1,2)*UCha(2,2)

	Sub96 = MD*MW2*(Sub95*SW) + Sub94*ZNeu(1,3)

	Sub153 = SB*Sub152 - 2*MU*(dSB1*SW*ZNeu(1,4))

	Sub141 = 3*CW*ZNeu(2,2) + SW*ZNeu(2,1)

	Sub142 = 3*CW*ZNeu(3,2) + SW*ZNeu(3,1)

	Sub143 = 3*CW*ZNeu(4,2) + SW*ZNeu(4,1)

        Sub89 = dZfR1(11,1,1)*ZNeu(1,1) + 
     &    dZfR1(11,2,1)*ZNeu(2,1) + dZfR1(11,3,1)*ZNeu(3,1) + 
     &    dZfR1(11,4,1)*ZNeu(4,1)

        Sub160 = dZfL1(11,1,1)*ZNeuC(1,1) + 
     &    dZfL1(11,2,1)*ZNeuC(2,1) + dZfL1(11,3,1)*ZNeuC(3,1) + 
     &    dZfL1(11,4,1)*ZNeuC(4,1)

        Sub161 = dZfL1(11,1,1)*ZNeuC(1,4) + 
     &    dZfL1(11,2,1)*ZNeuC(2,4) + dZfL1(11,3,1)*ZNeuC(3,4) + 
     &    dZfL1(11,4,1)*ZNeuC(4,4)

	Sub83 = MU*Sub82 - 2*MW2*(SB*SW*dMf1(3,1))

	Sub72 = 2*dSW1 + SW*(-2*dZe1 - dZbarfL1(3,1,1))

	Sub157 = 2*dSW1*SW + CW2*(2*dZe1 + dZbarfR1(3,1,1))

	Sub79 = MD*MW2*(CB*Sub78*SW) - Sub77*UCha(1,2)

	Sub154 = MW2*Sub153 - MU*(SB*Sub149*SW*ZNeu(1,4))

        Sub144 = Sub8*dZfR1(11,1,1) + Sub141*dZfR1(11,2,1) + 
     &    Sub142*dZfR1(11,3,1) + Sub143*dZfR1(11,4,1)

	Sub140 = 2*dSW1*SW + CW2*(2*dZe1 + dZbarfL1(3,1,1))

        Sub85 = dZbarfR1(12,1,1)*VChaC(1,2) + 
     &    dZbarfR1(12,1,2)*VChaC(2,2)

	Sub145 = Sub144*SW - 3*CW*(Sub72*ZNeu(1,2))

	Sub64 = -(3*CW*ZNeu(1,2)) + SW*ZNeu(1,1)

	Sub55 = 2*dZe1 + dZW1 + dZbarfL1(3,1,1) + dZfL1(4,1,1)

	pave23 = A0(MGl2)

        Sub175 = 1/Sqrt(2.D0)*
     &    (6*Pi*dSW1 + SW*(2*Alfas*Finite - 3*Pi*Sub55))/SW

	call Bput(pave3(bb0),S,MD2,MU2)

	call Cput(pave4(cc0),MD2,S,MU2,0.D0,MD2,MU2)

	Sub126 = CB*(dMWsq1*SW + 2*MW2*dSW1) + MW2*(Sub87*SW)

	Sub127 = MD*Sub126 - 2*MW2*(CB*SW*dMf1(4,1))

	Sub128 = MD*MW2*(CB*Sub78*SW) - Sub127*UCha(1,2)

	Sub234 = pave4(PVC(cc0)) + pave4(PVC(cc1)) + pave4(PVC(cc2))

	LOOP(Sfe5, 1,2,1)

        Sub103(Sfe5) = 
     &   dZbarSf1(1,Sfe5,4,1)*USf(1,2,4,1) + 
     &    dZbarSf1(2,Sfe5,4,1)*USf(2,2,4,1)

        Sub107(Sfe5) = 
     &   MD*MW2*(CB*SW*Sub103(Sfe5)) - Sub106*USf(Sfe5,2,4,1)

        Sub98(Sfe5) = dZbarSf1(1,Sfe5,4,1)*USf(1,1,4,1) + 
     &    dZbarSf1(2,Sfe5,4,1)*USf(2,1,4,1)

        Sub114(Sfe5) = 
     &   MW*(CB*Sub112*USf(Sfe5,1,4,1)) + 
     &    3*CW*MD*(Sub113*USf(Sfe5,2,4,1))

        Sub70(Sfe5) = dZSf1(1,Sfe5,4,1)*USfC(1,2,4,1) + 
     &    dZSf1(2,Sfe5,4,1)*USfC(2,2,4,1)

        Sub93(Sfe5) = Sub91*dZbarSf1(1,Sfe5,4,1) + 
     &    Sub92*dZbarSf1(2,Sfe5,4,1)

        Sub130(Sfe5) = 
     &   dZbarSf1(1,Sfe5,3,1)*USf(1,2,3,1) + 
     &    dZbarSf1(2,Sfe5,3,1)*USf(2,2,3,1)

	Sub100(Sfe5) = CW2*Sub98(Sfe5) + Sub99*USf(Sfe5,1,4,1)

        Sub146(Sfe5) = 
     &   dZSf1(1,Sfe5,3,1)*USfC(1,1,3,1) + 
     &    dZSf1(2,Sfe5,3,1)*USfC(2,1,3,1)

        Sub125(Sfe5) = 
     &   dZbarSf1(1,Sfe5,3,1)*USf(1,1,3,1) + 
     &    dZbarSf1(2,Sfe5,3,1)*USf(2,1,3,1)

        Sub97(Sfe5) = MW2*(SW*Sub93(Sfe5)) + 
     &    3*CW*(Sub96*USf(Sfe5,1,4,1))

        Sub134(Sfe5) = 
     &   MU*MW2*(SB*SW*Sub130(Sfe5)) - Sub133*USf(Sfe5,2,3,1)

        Sub69(Sfe5) = dZSf1(1,Sfe5,4,1)*USfC(1,1,4,1) + 
     &    dZSf1(2,Sfe5,4,1)*USfC(2,1,4,1)

        Sub147(Sfe5) = 
     &   dZSf1(1,Sfe5,3,1)*USfC(1,2,3,1) + 
     &    dZSf1(2,Sfe5,3,1)*USfC(2,2,3,1)

        Sub159(Sfe5) = 
     &   MU*MW2*(SB*SW*Sub146(Sfe5)) - Sub83*USfC(Sfe5,1,3,1)

        Sub162(Sfe5) = 
     &   3*CW*MU*(Sub161*USfC(Sfe5,1,3,1)) - 
     &    4*MW*(SB*Sub160*SW*USfC(Sfe5,2,3,1))

        Sub90(Sfe5) = 2*MW*SW2*
     &     (CB2*USf(Sfe5,2,4,1)*(CW2*Sub89 + Sub88*ZNeu(1,1))) - 
     &    3*CW*CW2*MD*(Sub87*SW*USf(Sfe5,1,4,1)*ZNeu(1,3))

        Sub148(Sfe5) = 
     &   MW*(SB*Sub8*Sub146(Sfe5)) + 
     &    3*CW*MU*(Sub147(Sfe5)*ZNeu(1,4))

        Sub84(Sfe5) = MU*MW2*(SB*SW*Sub69(Sfe5)) - 
     &    Sub83*USfC(Sfe5,1,4,1)

        Sub155(Sfe5) = 
     &   MW2*(SB*SW*Sub148(Sfe5)) + 3*CW*(Sub154*USfC(Sfe5,2,3,1))

	Sub158(Sfe5) = CW2*Sub147(Sfe5) + Sub157*USfC(Sfe5,2,3,1)

        Sub163(Sfe5) = 
     &   MW2*(SB*SW*Sub162(Sfe5)) + 3*CW*(Sub159(Sfe5)*ZNeuC(1,4))

        call Cput(pave10(cc0,Sfe5),S,MNeu2(1),MCha2(1),MD2,MU2,
     &   MSf2(Sfe5,3,1))

        Sub86(Sfe5) = MU*MW2*(SB*Sub85*SW*USfC(Sfe5,1,4,1)) + 
     &    Sub84(Sfe5)*VChaC(1,2)

        Sub156(Sfe5) = 
     &   CW2*Sub155(Sfe5) + 
     &    MW*MW2*(SB2*(CW2*Sub145 + SW2*(Sub140*ZNeu(1,1)))*
     &       USfC(Sfe5,1,3,1))

        Sub164(Sfe5) = 
     &   CW2*Sub163(Sfe5) - 
     &    4*MW*MW2*SW2*(SB2*Sub158(Sfe5)*ZNeuC(1,1))

        Sub211(Sfe5) = 
     &   MW*(CB*Sub64*USfC(Sfe5,1,4,1)) + 
     &    3*CW*MD*(USfC(Sfe5,2,4,1)*ZNeu(1,3))

	Sub102(Sfe5) = SW*Sub98(Sfe5) - Sub101*USf(Sfe5,1,4,1)

        Sub108(Sfe5) = 
     &   MW*MW2*(CB2*Sub102(Sfe5)*ZNeuC(1,2)) - 
     &    Sub107(Sfe5)*ZNeuC(1,3)

	Sub115(Sfe5) = 3*CW*Sub108(Sfe5) - MW2*(CB*SW*Sub114(Sfe5))

        Sub71(Sfe5) = -(2*MW*(CB*Sub69(Sfe5)*UCha(1,1))) + 
     &    Sqrt(2.D0)*MD*(Sub70(Sfe5)*UCha(1,2))

	Sub135(Sfe5) = SW*Sub125(Sfe5) - Sub101*USf(Sfe5,1,3,1)

        Sub137(Sfe5) = 
     &   -(2*MW*(SB*Sub136*USf(Sfe5,1,3,1))) + 
     &    Sqrt(2.D0)*MU*(Sub85*USf(Sfe5,2,3,1))

        Sub74(Sfe5) = -(2*MW*
     &       (CB2*(Sub73*SW - Sub72*UCha(1,1))*USfC(Sfe5,1,4,1)))+
     &      CB*SW*Sub71(Sfe5)

        Sub138(Sfe5) = 
     &   -(2*MW*(SB2*Sub135(Sfe5)*VChaC(1,1))) + SB*SW*Sub137(Sfe5)

        Sub116(Sfe5) = 
     &   CW2*Sub115(Sfe5) - 
     &    MW*MW2*SW2*(CB2*Sub100(Sfe5)*ZNeuC(1,1))

        Sub129(Sfe5) = 
     &   MD*MW2*(CB*SW*Sub125(Sfe5)*UCha(1,2)) + 
     &    Sub128*USf(Sfe5,1,3,1)

        Sub80(Sfe5) = MW2*Sub74(Sfe5) + 
     &    Sqrt(2.D0)*(Sub79*USfC(Sfe5,2,4,1))

        Sub139(Sfe5) = 
     &   MW2*Sub138(Sfe5) + Sqrt(2.D0)*(Sub134(Sfe5)*VChaC(1,2))

	ENDLOOP(Sfe5)

	LOOP(Sfe6, 1,2,1)

        Sub120(Sfe6) = 
     &   4*MW*(SB*SW*USf(Sfe6,2,3,1)*ZNeu(1,1)) - 
     &    3*CW*MU*(USf(Sfe6,1,3,1)*ZNeu(1,4))

        Sub66(Sfe6) = 2*MW*
     &     (CB*SW*USfC(Sfe6,2,4,1)*ZNeuC(1,1)) + 
     &    3*CW*MD*(USfC(Sfe6,1,4,1)*ZNeuC(1,3))

        Sub121(Sfe6) = 
     &   MW*(SB*USf(Sfe6,1,3,1)*ZNeuC(1,2)) + 
     &    MU*(USf(Sfe6,2,3,1)*ZNeuC(1,4))

        Sub15(Sfe6) = -(2*MW*
     &       (CB*UCha(1,1)*USfC(Sfe6,1,4,1))) + 
     &    Sqrt(2.D0)*MD*(UCha(1,2)*USfC(Sfe6,2,4,1))

        Sub63(Sfe6) = -(2*MW*
     &       (SB*USf(Sfe6,1,3,1)*VChaC(1,1))) + 
     &    Sqrt(2.D0)*MU*(USf(Sfe6,2,3,1)*VChaC(1,2))

        Sub28(Sfe6) = MW*(SB*Sub8*USfC(Sfe6,1,3,1)) + 
     &    3*CW*MU*(USfC(Sfe6,2,3,1)*ZNeu(1,4))

        Sub122(Sfe6) = 
     &   3*CW*Sub121(Sfe6) + MW*(SB*SW*USf(Sfe6,1,3,1)*ZNeuC(1,1))

        Sub29(Sfe6) = 4*MW*
     &     (SB*SW*USfC(Sfe6,2,3,1)*ZNeuC(1,1)) - 
     &    3*CW*MU*(USfC(Sfe6,1,3,1)*ZNeuC(1,4))

        Sub65(Sfe6) = MW*(CB*Sub64*USfC(Sfe6,1,4,1)) + 
     &    3*CW*MD*(USfC(Sfe6,2,4,1)*ZNeu(1,3))

	ENDLOOP(Sfe6)

	LOOP(Sfe7, 1,2,1)

	pave1(Sfe7) = A0(MSf2(Sfe7,4,1))

	pave2(Sfe7) = A0(MSf2(Sfe7,3,1))

	ENDLOOP(Sfe7)

	LOOP(Sfe6, 1,2,1)
	LOOP(Sfe5, 1,2,1)

        Sub14(Sfe5,Sfe6) = 
     &   dZbarSf1(Sfe5,Sfe6,3,1) + dZSf1(Sfe6,Sfe5,3,1)

        Sub12(Sfe5,Sfe6) = 
     &   dZbarSf1(Sfe5,Sfe6,4,1) + dZSf1(Sfe6,Sfe5,4,1)

        Sub38(Sfe5,Sfe6) = 
     &   USf(Sfe6,1,3,1)*USfC(Sfe5,1,3,1) + 
     &    USf(Sfe6,2,3,1)*USfC(Sfe5,2,3,1)

        Sub32(Sfe5,Sfe6) = 
     &   USf(Sfe6,1,4,1)*USfC(Sfe5,1,4,1) + 
     &    USf(Sfe6,2,4,1)*USfC(Sfe5,2,4,1)

        Sub13(Sfe5,Sfe6) = 
     &   2*dMSfsq1(Sfe6,Sfe5,3,1) + 
     &    dZbarSf1(Sfe5,Sfe6,3,1)*MSf2(Sfe5,3,1) + 
     &    dZSf1(Sfe6,Sfe5,3,1)*MSf2(Sfe6,3,1)

        Sub11(Sfe5,Sfe6) = 
     &   2*dMSfsq1(Sfe6,Sfe5,4,1) + 
     &    dZbarSf1(Sfe5,Sfe6,4,1)*MSf2(Sfe5,4,1) + 
     &    dZSf1(Sfe6,Sfe5,4,1)*MSf2(Sfe6,4,1)

        call Cput(pave18(cc0,Sfe5,Sfe6),MD2,S,MU2,MGl2,
     &   MSf2(Sfe5,4,1),MSf2(Sfe6,3,1))

	ENDLOOP(Sfe5)
	ENDLOOP(Sfe6)

	LOOP(Sfe7, 1,2,1)
	LOOP(Sfe5, 1,2,1)

        Sub47(Sfe5,Sfe7) = 
     &   USf(Sfe7,1,3,1)*USfC(Sfe5,1,3,1) - 
     &    USf(Sfe7,2,3,1)*USfC(Sfe5,2,3,1)

        Sub39(Sfe5,Sfe7) = 
     &   USf(Sfe7,1,4,1)*USfC(Sfe5,1,4,1) - 
     &    USf(Sfe7,2,4,1)*USfC(Sfe5,2,4,1)

	ENDLOOP(Sfe5)
	ENDLOOP(Sfe7)

	LOOP(Sfe7, 1,2,1)
	LOOP(Sfe6, 1,2,1)

        Sub48(Sfe6,Sfe7) = 
     &   USf(Sfe6,1,3,1)*USfC(Sfe7,1,3,1) - 
     &    USf(Sfe6,2,3,1)*USfC(Sfe7,2,3,1)

        Sub40(Sfe6,Sfe7) = 
     &   USf(Sfe6,1,4,1)*USfC(Sfe7,1,4,1) - 
     &    USf(Sfe6,2,4,1)*USfC(Sfe7,2,4,1)

	ENDLOOP(Sfe6)
	ENDLOOP(Sfe7)

	LOOP(Sfe7, 1,2,1)
	LOOP(Sfe6, 1,2,1)
	LOOP(Sfe5, 1,2,1)

        Sub46(Sfe5,Sfe6,Sfe7) = 
     &   48*Sqrt(2.D0)*Alfa*(Alfas*Sub39(Sfe5,Sfe7)*Sub40(Sfe6,Sfe7))

        Sub54(Sfe5,Sfe6,Sfe7) = 
     &   48*Sqrt(2.D0)*Alfa*(Alfas*Sub47(Sfe5,Sfe7)*Sub48(Sfe6,Sfe7))

	ENDLOOP(Sfe5)
	ENDLOOP(Sfe6)
	ENDLOOP(Sfe7)

#include "contains.h"

	end

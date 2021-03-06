#if 0
* udbar_n2x2_abbr1s.F
* abbreviations for {F[3, {1, Col1}], -F[4, {1, Col2}]} -> {F[11, {2}], -F[12, {2}]}
* generated by FormCalc 8.4 on 19-Jun-2015 12:48
#endif

#include "udbar_n2x2_vars.h"
#include "contains.h"


	subroutine udbar_n2x2_abbr1s
	implicit none

#include "udbar_n2x2_vars.h"


        Sub154 = dZfR1(11,1,2)*ZNeu(1,3) + 
     &    dZfR1(11,2,2)*ZNeu(2,3) + dZfR1(11,3,2)*ZNeu(3,3) + 
     &    dZfR1(11,4,2)*ZNeu(4,3)

	Sub111 = dMWsq1 + MW2*(-2*dZe1 - dZfL1(3,1,1))

	Sub155 = MD*Sub154 + 2*dMf1(4,1)*ZNeu(2,3)

	Sub112 = SB*Sub111*SW + 2*MW2*(dSW1*SB + dSB1*SW)

	Sub82 = dMWsq1 + MW2*(-2*dZe1 - dZbarfL1(4,1,1))

	Sub77 = 2*dCB1 - CB*(2*dZe1 + dZbarfR1(4,1,1))

	Sub130 = 2*dCB1 - CB*(2*dZe1 + dZfL1(3,1,1))

	Sub156 = Sub155*SW - 2*MD*(dSW1*ZNeu(2,3))

	Sub131 = CB*(dMWsq1*SW + 2*MW2*dSW1) + MW2*(Sub130*SW)

	Sub78 = CB*(dMWsq1*SW + 2*MW2*dSW1) + MW2*(Sub77*SW)

	Sub83 = SB*Sub82*SW + 2*MW2*(dSW1*SB + dSB1*SW)

	Sub113 = MU*Sub112 - 2*MW2*(SB*SW*dMf1(3,1))

	Sub153 = 2*dCB1 - CB*(2*dZe1 + dZbarfL1(4,1,1))

	Sub108 = 2*dSW1 + SW*(-2*dZe1 - dZfL1(3,1,1))

	Sub136 = dZfL1(12,1,2)*UChaC(1,1) + dZfL1(12,2,2)*UChaC(2,1)

	Sub145 = -(3*CW*ZNeu(1,2)) + SW*ZNeu(1,1)

	Sub157 = MW2*Sub156 - MD*(dMWsq1*SW*ZNeu(2,3))

	Sub147 = -(3*CW*ZNeu(3,2)) + SW*ZNeu(3,1)

	Sub146 = -(3*CW*ZNeu(4,2)) + SW*ZNeu(4,1)

	Sub116 = 3*CW*ZNeuC(1,2) + SW*ZNeuC(1,1)

	Sub117 = 3*CW*ZNeuC(2,2) + SW*ZNeuC(2,1)

	Sub118 = 3*CW*ZNeuC(3,2) + SW*ZNeuC(3,1)

	Sub119 = 3*CW*ZNeuC(4,2) + SW*ZNeuC(4,1)

        Sub164 = dZfL1(11,1,2)*ZNeuC(1,1) + 
     &    dZfL1(11,2,2)*ZNeuC(2,1) + dZfL1(11,3,2)*ZNeuC(3,1) + 
     &    dZfL1(11,4,2)*ZNeuC(4,1)

        Sub165 = dZfL1(11,1,2)*ZNeuC(1,3) + 
     &    dZfL1(11,2,2)*ZNeuC(2,3) + dZfL1(11,3,2)*ZNeuC(3,3) + 
     &    dZfL1(11,4,2)*ZNeuC(4,3)

	Sub100 = -(SB*SW*dMf1(3,1)) + MU*(dSW1*SB + dSB1*SW)

	Sub84 = MU*Sub83 - 2*MW2*(SB*SW*dMf1(3,1))

	Sub132 = MD*Sub131 - 2*MW2*(CB*SW*dMf1(4,1))

	Sub79 = MD*Sub78 - 2*MW2*(CB*SW*dMf1(4,1))

	Sub87 = 2*dSW1 + SW*(-2*dZe1 - dZbarfL1(4,1,1))

	Sub161 = 2*dSW1*SW + CW2*(2*dZe1 + dZbarfR1(4,1,1))

	Sub88 = dZfR1(12,1,2)*VCha(1,1) + dZfR1(12,2,2)*VCha(2,1)

	Sub158 = CB*Sub157 - MD*MW2*(Sub153*SW*ZNeu(2,3))

        Sub93 = 4*MW*(SB*SW*USf(1,2,3,1)*ZNeu(2,1)) - 
     &    3*CW*MU*(USf(1,1,3,1)*ZNeu(2,4))

        Sub94 = 4*MW*(SB*SW*USf(2,2,3,1)*ZNeu(2,1)) - 
     &    3*CW*MU*(USf(2,1,3,1)*ZNeu(2,4))

        Sub120 = Sub116*dZfL1(11,1,2) + Sub117*dZfL1(11,2,2) + 
     &    Sub118*dZfL1(11,3,2) + Sub119*dZfL1(11,4,2)

        Sub148 = Sub145*dZfR1(11,1,2) + Sub8*dZfR1(11,2,2) + 
     &    Sub147*dZfR1(11,3,2) + Sub146*dZfR1(11,4,2)

        Sub101 = dZfR1(11,1,2)*ZNeu(1,4) + 
     &    dZfR1(11,2,2)*ZNeu(2,4) + dZfR1(11,3,2)*ZNeu(3,4) + 
     &    dZfR1(11,4,2)*ZNeu(4,4)

        Sub121 = dZfL1(11,1,2)*ZNeuC(1,4) + 
     &    dZfL1(11,2,2)*ZNeuC(2,4) + dZfL1(11,3,2)*ZNeuC(3,4) + 
     &    dZfL1(11,4,2)*ZNeuC(4,4)

	Sub144 = 2*dSW1*SW + CW2*(2*dZe1 + dZbarfL1(4,1,1))

	Sub106 = 2*dSW1*SW + CW2*(2*dZe1 + dZfL1(3,1,1))

	Sub96 = 2*dSW1*SW + CW2*dZfR1(3,1,1)

	Sub99 = dMWsq1 + MW2*(-2*dZe1 - dZfR1(3,1,1))

	Sub75 = dZfL1(12,1,2)*UChaC(1,2) + dZfL1(12,2,2)*UChaC(2,2)

	Sub85 = dZfR1(12,1,2)*VCha(1,2) + dZfR1(12,2,2)*VCha(2,2)

	Sub149 = Sub148*SW + 3*CW*(Sub87*ZNeu(2,2))

	Sub102 = MU*(SB*Sub101*SW) - 2*Sub100*ZNeu(2,4)

        Sub97 = dZfR1(11,1,2)*ZNeu(1,1) + 
     &    (2*dZe1 + dZfR1(11,2,2))*ZNeu(2,1) + 
     &    dZfR1(11,3,2)*ZNeu(3,1) + dZfR1(11,4,2)*ZNeu(4,1)

	Sub103 = MW2*Sub102 - MU*(SB*Sub99*SW*ZNeu(2,4))

	Sub68 = 3*CW*ZNeu(2,2) + SW*ZNeu(2,1)

	Sub65 = -(3*CW*ZNeuC(2,2)) + SW*ZNeuC(2,1)

	Sub56 = 2*dZe1 + dZW1 + dZbarfL1(4,1,1) + dZfL1(3,1,1)

	pave22 = A0(MGl2)

        Sub178 = 1/Sqrt(2.D0)*
     &    (6*Pi*dSW1 + SW*(2*Alfas*Finite - 3*Pi*Sub56))/SW

	call Bput(pave3(bb0),S,MD2,MU2)

	call Cput(pave4(cc0),MD2,S,MU2,0.D0,MD2,MU2)

	Sub140 = SB*Sub99*SW + 2*MW2*(dSW1*SB + dSB1*SW)

	Sub141 = MU*Sub140 - 2*MW2*(SB*SW*dMf1(3,1))

	Sub86 = MU*MW2*(SB*Sub85*SW) - Sub84*VCha(2,2)

	Sub142 = MU*MW2*(SB*Sub85*SW) - Sub141*VCha(2,2)

	Sub234 = pave4(PVC(cc0)) + pave4(PVC(cc1)) + pave4(PVC(cc2))

	LOOP(Sfe5, 1,2,1)

        Sub110(Sfe5) = 
     &   dZbarSf1(1,Sfe5,3,1)*USf(1,2,3,1) + 
     &    dZbarSf1(2,Sfe5,3,1)*USf(2,2,3,1)

        Sub150(Sfe5) = 
     &   dZSf1(1,Sfe5,4,1)*USfC(1,1,4,1) + 
     &    dZSf1(2,Sfe5,4,1)*USfC(2,1,4,1)

        Sub89(Sfe5) = dZSf1(1,Sfe5,3,1)*USfC(1,2,3,1) + 
     &    dZSf1(2,Sfe5,3,1)*USfC(2,2,3,1)

        Sub129(Sfe5) = 
     &   dZbarSf1(1,Sfe5,4,1)*USf(1,2,4,1) + 
     &    dZbarSf1(2,Sfe5,4,1)*USf(2,2,4,1)

        Sub114(Sfe5) = 
     &   MU*MW2*(SB*SW*Sub110(Sfe5)) - Sub113*USf(Sfe5,2,3,1)

        Sub76(Sfe5) = dZSf1(1,Sfe5,3,1)*USfC(1,1,3,1) + 
     &    dZSf1(2,Sfe5,3,1)*USfC(2,1,3,1)

        Sub151(Sfe5) = 
     &   dZSf1(1,Sfe5,4,1)*USfC(1,2,4,1) + 
     &    dZSf1(2,Sfe5,4,1)*USfC(2,2,4,1)

        Sub163(Sfe5) = 
     &   MD*MW2*(CB*SW*Sub150(Sfe5)) - Sub79*USfC(Sfe5,1,4,1)

        Sub166(Sfe5) = 
     &   3*CW*MD*(Sub165*USfC(Sfe5,1,4,1)) + 
     &    2*MW*(CB*Sub164*SW*USfC(Sfe5,2,4,1))

        Sub90(Sfe5) = -(2*MW*(SB*Sub76(Sfe5)*VCha(2,1))) + 
     &    Sqrt(2.D0)*MU*(Sub89(Sfe5)*VCha(2,2))

        Sub152(Sfe5) = 
     &   MW*(CB*Sub8*Sub150(Sfe5)) + 
     &    3*CW*MD*(Sub151(Sfe5)*ZNeu(2,3))

        Sub95(Sfe5) = Sub93*dZbarSf1(1,Sfe5,3,1) + 
     &    Sub94*dZbarSf1(2,Sfe5,3,1)

        Sub105(Sfe5) = 
     &   dZbarSf1(1,Sfe5,3,1)*USf(1,1,3,1) + 
     &    dZbarSf1(2,Sfe5,3,1)*USf(2,1,3,1)

        Sub122(Sfe5) = 
     &   MW*(SB*Sub120*USf(Sfe5,1,3,1)) + 
     &    3*CW*MU*(Sub121*USf(Sfe5,2,3,1))

        Sub133(Sfe5) = 
     &   MD*MW2*(CB*SW*Sub129(Sfe5)) - Sub132*USf(Sfe5,2,4,1)

        Sub80(Sfe5) = MD*MW2*(CB*SW*Sub76(Sfe5)) - 
     &    Sub79*USfC(Sfe5,1,3,1)

        Sub159(Sfe5) = 
     &   MW2*(CB*SW*Sub152(Sfe5)) + 3*CW*(Sub158*USfC(Sfe5,2,4,1))

	Sub162(Sfe5) = CW2*Sub151(Sfe5) + Sub161*USfC(Sfe5,2,4,1)

        Sub91(Sfe5) = -(2*MW*
     &       (SB2*USfC(Sfe5,1,3,1)*(Sub88*SW - Sub87*VCha(2,1))))+
     &      SB*SW*Sub90(Sfe5)

        Sub167(Sfe5) = 
     &   MW2*(CB*SW*Sub166(Sfe5)) + 3*CW*(Sub163(Sfe5)*ZNeuC(2,3))

        Sub134(Sfe5) = 
     &   dZbarSf1(1,Sfe5,4,1)*USf(1,1,4,1) + 
     &    dZbarSf1(2,Sfe5,4,1)*USf(2,1,4,1)

	Sub107(Sfe5) = CW2*Sub105(Sfe5) + Sub106*USf(Sfe5,1,3,1)

        Sub81(Sfe5) = MD*MW2*(CB*Sub75*SW*USfC(Sfe5,1,3,1)) + 
     &    Sub80(Sfe5)*UChaC(2,2)

        Sub98(Sfe5) = CW2*(SB*SW*Sub95(Sfe5)) + 
     &    4*MW*SW2*(SB2*USf(Sfe5,2,3,1)*
     &       (CW2*Sub97 + Sub96*ZNeu(2,1)))

        Sub160(Sfe5) = 
     &   CW2*Sub159(Sfe5) + 
     &    MW*MW2*(CB2*(CW2*Sub149 + SW2*(Sub144*ZNeu(2,1)))*
     &       USfC(Sfe5,1,4,1))

        Sub168(Sfe5) = 
     &   CW2*Sub167(Sfe5) + 
     &    2*MW*MW2*SW2*(CB2*Sub162(Sfe5)*ZNeuC(2,1))

        call Cput(pave12(cc0,Sfe5),S,MCha2(2),MNeu2(2),MD2,MU2,
     &   MSf2(Sfe5,4,1))

        call Cput(pave13(cc0,Sfe5),S,MNeu2(2),MCha2(2),MD2,MU2,
     &   MSf2(Sfe5,3,1))

        Sub104(Sfe5) = 
     &   MW2*Sub98(Sfe5) - 3*CW*CW2*(Sub103*USf(Sfe5,1,3,1))

        Sub213(Sfe5) = 
     &   MW*(SB*Sub68*USfC(Sfe5,1,3,1)) + 
     &    3*CW*MU*(USfC(Sfe5,2,3,1)*ZNeu(2,4))

	Sub109(Sfe5) = SW*Sub105(Sfe5) - Sub108*USf(Sfe5,1,3,1)

	Sub135(Sfe5) = SW*Sub134(Sfe5) - Sub108*USf(Sfe5,1,4,1)

        Sub137(Sfe5) = 
     &   -(2*MW*(CB*Sub136*USf(Sfe5,1,4,1))) + 
     &    Sqrt(2.D0)*MD*(Sub75*USf(Sfe5,2,4,1))

        Sub138(Sfe5) = 
     &   -(2*MW*(CB2*Sub135(Sfe5)*UChaC(2,1))) + CB*SW*Sub137(Sfe5)

        Sub115(Sfe5) = 
     &   MW*MW2*(SB2*Sub109(Sfe5)*ZNeuC(2,2)) + 
     &    Sub114(Sfe5)*ZNeuC(2,4)

	Sub123(Sfe5) = 3*CW*Sub115(Sfe5) + MW2*(SB*SW*Sub122(Sfe5))

        Sub139(Sfe5) = 
     &   MW2*Sub138(Sfe5) + Sqrt(2.D0)*(Sub133(Sfe5)*UChaC(2,2))

        Sub92(Sfe5) = MW2*Sub91(Sfe5) + 
     &    Sqrt(2.D0)*(Sub86*USfC(Sfe5,2,3,1))

        Sub143(Sfe5) = 
     &   MU*MW2*(SB*SW*Sub134(Sfe5)*VCha(2,2)) + 
     &    Sub142*USf(Sfe5,1,4,1)

        Sub124(Sfe5) = 
     &   CW2*Sub123(Sfe5) + 
     &    MW*MW2*SW2*(SB2*Sub107(Sfe5)*ZNeuC(2,1))

	ENDLOOP(Sfe5)

	LOOP(Sfe6, 1,2,1)

        Sub70(Sfe6) = 4*MW*
     &     (SB*SW*USfC(Sfe6,2,3,1)*ZNeuC(2,1)) - 
     &    3*CW*MU*(USfC(Sfe6,1,3,1)*ZNeuC(2,4))

        Sub64(Sfe6) = 2*MW*(CB*SW*USf(Sfe6,2,4,1)*ZNeu(2,1)) + 
     &    3*CW*MD*(USf(Sfe6,1,4,1)*ZNeu(2,3))

        Sub72(Sfe6) = -(2*MW*
     &       (CB*UChaC(2,1)*USf(Sfe6,1,4,1))) + 
     &    Sqrt(2.D0)*MD*(UChaC(2,2)*USf(Sfe6,2,4,1))

        Sub20(Sfe6) = -(2*MW*
     &       (SB*USfC(Sfe6,1,3,1)*VCha(2,1))) + 
     &    Sqrt(2.D0)*MU*(USfC(Sfe6,2,3,1)*VCha(2,2))

        Sub30(Sfe6) = MW*(CB*Sub8*USfC(Sfe6,1,4,1)) + 
     &    3*CW*MD*(USfC(Sfe6,2,4,1)*ZNeu(2,3))

        Sub66(Sfe6) = MW*(CB*Sub65*USf(Sfe6,1,4,1)) + 
     &    3*CW*MD*(USf(Sfe6,2,4,1)*ZNeuC(2,3))

        Sub31(Sfe6) = 2*MW*
     &     (CB*SW*USfC(Sfe6,2,4,1)*ZNeuC(2,1)) + 
     &    3*CW*MD*(USfC(Sfe6,1,4,1)*ZNeuC(2,3))

        Sub69(Sfe6) = MW*(SB*Sub68*USfC(Sfe6,1,3,1)) + 
     &    3*CW*MU*(USfC(Sfe6,2,3,1)*ZNeu(2,4))

	ENDLOOP(Sfe6)

	LOOP(Sfe7, 1,2,1)

	pave1(Sfe7) = A0(MSf2(Sfe7,3,1))

	pave2(Sfe7) = A0(MSf2(Sfe7,4,1))

	ENDLOOP(Sfe7)

	LOOP(Sfe6, 1,2,1)
	LOOP(Sfe5, 1,2,1)

        Sub13(Sfe5,Sfe6) = 
     &   dZbarSf1(Sfe5,Sfe6,3,1) + dZSf1(Sfe6,Sfe5,3,1)

        Sub14(Sfe5,Sfe6) = 
     &   2*dMSfsq1(Sfe6,Sfe5,3,1) + 
     &    dZbarSf1(Sfe5,Sfe6,3,1)*MSf2(Sfe5,3,1) + 
     &    dZSf1(Sfe6,Sfe5,3,1)*MSf2(Sfe6,3,1)

        Sub12(Sfe5,Sfe6) = 
     &   dZbarSf1(Sfe5,Sfe6,4,1) + dZSf1(Sfe6,Sfe5,4,1)

        Sub37(Sfe5,Sfe6) = 
     &   USf(Sfe6,1,3,1)*USfC(Sfe5,1,3,1) + 
     &    USf(Sfe6,2,3,1)*USfC(Sfe5,2,3,1)

        Sub38(Sfe5,Sfe6) = 
     &   USf(Sfe6,1,4,1)*USfC(Sfe5,1,4,1) + 
     &    USf(Sfe6,2,4,1)*USfC(Sfe5,2,4,1)

        Sub11(Sfe5,Sfe6) = 
     &   2*dMSfsq1(Sfe6,Sfe5,4,1) + 
     &    dZbarSf1(Sfe5,Sfe6,4,1)*MSf2(Sfe5,4,1) + 
     &    dZSf1(Sfe6,Sfe5,4,1)*MSf2(Sfe6,4,1)

        call Cput(pave19(cc0,Sfe5,Sfe6),MU2,S,MD2,MGl2,
     &   MSf2(Sfe5,3,1),MSf2(Sfe6,4,1))

	ENDLOOP(Sfe5)
	ENDLOOP(Sfe6)

	LOOP(Sfe7, 1,2,1)
	LOOP(Sfe5, 1,2,1)

        Sub40(Sfe5,Sfe7) = 
     &   USf(Sfe7,1,3,1)*USfC(Sfe5,1,3,1) - 
     &    USf(Sfe7,2,3,1)*USfC(Sfe5,2,3,1)

        Sub48(Sfe5,Sfe7) = 
     &   USf(Sfe7,1,4,1)*USfC(Sfe5,1,4,1) - 
     &    USf(Sfe7,2,4,1)*USfC(Sfe5,2,4,1)

	ENDLOOP(Sfe5)
	ENDLOOP(Sfe7)

	LOOP(Sfe7, 1,2,1)
	LOOP(Sfe6, 1,2,1)

        Sub41(Sfe6,Sfe7) = 
     &   USf(Sfe6,1,3,1)*USfC(Sfe7,1,3,1) - 
     &    USf(Sfe6,2,3,1)*USfC(Sfe7,2,3,1)

        Sub49(Sfe6,Sfe7) = 
     &   USf(Sfe6,1,4,1)*USfC(Sfe7,1,4,1) - 
     &    USf(Sfe6,2,4,1)*USfC(Sfe7,2,4,1)

	ENDLOOP(Sfe6)
	ENDLOOP(Sfe7)

	LOOP(Sfe7, 1,2,1)
	LOOP(Sfe6, 1,2,1)
	LOOP(Sfe5, 1,2,1)

        Sub47(Sfe5,Sfe6,Sfe7) = 
     &   48*Sqrt(2.D0)*Alfa*(Alfas*Sub40(Sfe5,Sfe7)*Sub41(Sfe6,Sfe7))

        Sub55(Sfe5,Sfe6,Sfe7) = 
     &   48*Sqrt(2.D0)*Alfa*(Alfas*Sub48(Sfe5,Sfe7)*Sub49(Sfe6,Sfe7))

	ENDLOOP(Sfe5)
	ENDLOOP(Sfe6)
	ENDLOOP(Sfe7)

#include "contains.h"

	end

#if 0
* uubar_n1n2_abbr1s.F
* abbreviations for {F[3, {1, Col1}], -F[3, {1, Col2}]} -> {F[11, {1}], F[11, {2}]}
* generated by FormCalc 8.4 on 11-Jun-2015 10:57
#endif

#include "uubar_n1n2_vars.h"
#include "contains.h"


	subroutine uubar_n1n2_abbr1s
	implicit none

#include "uubar_n1n2_vars.h"


	Sub84 = dMWsq1 + MW2*(-2*dZe1 - dZfL1(3,1,1))

	Sub85 = SB*Sub84*SW + 2*MW2*(dSW1*SB + dSB1*SW)

	Sub101 = dMWsq1 + MW2*(-2*dZe1 - dZbarfR1(3,1,1))

	Sub102 = SB*Sub101*SW + 2*MW2*(dSW1*SB + dSB1*SW)

	Sub86 = MU*Sub85 - 2*MW2*(SB*SW*dMf1(3,1))

	Sub81 = 2*dSW1 + SW*(-2*dZe1 - dZfL1(3,1,1))

        Sub111 = 4*MW*(SB*SW*USf(1,2,3,1)*ZNeu(2,1)) - 
     &    3*CW*MU*(USf(1,1,3,1)*ZNeu(2,4))

        Sub112 = 4*MW*(SB*SW*USf(2,2,3,1)*ZNeu(2,1)) - 
     &    3*CW*MU*(USf(2,1,3,1)*ZNeu(2,4))

	Sub89 = 3*CW*ZNeuC(1,2) + SW*ZNeuC(1,1)

	Sub90 = 3*CW*ZNeuC(2,2) + SW*ZNeuC(2,1)

	Sub91 = 3*CW*ZNeuC(3,2) + SW*ZNeuC(3,1)

	Sub92 = 3*CW*ZNeuC(4,2) + SW*ZNeuC(4,1)

        Sub135 = Sub89*dZfL1(11,1,2) + Sub90*dZfL1(11,2,2) + 
     &    Sub91*dZfL1(11,3,2) + Sub92*dZfL1(11,4,2)

        Sub70 = dZfR1(11,1,2)*ZNeu(1,4) + 
     &    dZfR1(11,2,2)*ZNeu(2,4) + dZfR1(11,3,2)*ZNeu(3,4) + 
     &    dZfR1(11,4,2)*ZNeu(4,4)

        Sub129 = dZfL1(11,1,1)*ZNeuC(1,1) + 
     &    dZfL1(11,2,1)*ZNeuC(2,1) + dZfL1(11,3,1)*ZNeuC(3,1) + 
     &    dZfL1(11,4,1)*ZNeuC(4,1)

	Sub55 = -(SB*SW*dMf1(3,1)) + MU*(dSW1*SB + dSB1*SW)

	Sub103 = MU*Sub102 - 2*MW2*(SB*SW*dMf1(3,1))

	Sub69 = dMWsq1 + MW2*(-2*dZe1 - dZbarfL1(3,1,1))

        Sub48 = 4*MW*(SB*SW*USf(1,2,3,1)*ZNeu(1,1)) - 
     &    3*CW*MU*(USf(1,1,3,1)*ZNeu(1,4))

        Sub49 = 4*MW*(SB*SW*USf(2,2,3,1)*ZNeu(1,1)) - 
     &    3*CW*MU*(USf(2,1,3,1)*ZNeu(1,4))

	Sub116 = MU*(SB*Sub70*SW) - 2*Sub55*ZNeu(2,4)

	Sub63 = 3*CW*ZNeu(3,2) + SW*ZNeu(3,1)

	Sub62 = 3*CW*ZNeu(4,2) + SW*ZNeu(4,1)

        Sub93 = Sub89*dZfL1(11,1,1) + Sub90*dZfL1(11,2,1) + 
     &    Sub91*dZfL1(11,3,1) + Sub92*dZfL1(11,4,1)

        Sub119 = Sub11*dZfR1(11,1,1) + Sub1*dZfR1(11,2,1) + 
     &    Sub63*dZfR1(11,3,1) + Sub62*dZfR1(11,4,1)

        Sub114 = dZfR1(11,1,2)*ZNeu(1,1) + 
     &    (2*dZe1 + dZfR1(11,2,2))*ZNeu(2,1) + 
     &    dZfR1(11,3,2)*ZNeu(3,1) + dZfR1(11,4,2)*ZNeu(4,1)

        Sub56 = dZfR1(11,1,1)*ZNeu(1,4) + 
     &    dZfR1(11,2,1)*ZNeu(2,4) + dZfR1(11,3,1)*ZNeu(3,4) + 
     &    dZfR1(11,4,1)*ZNeu(4,4)

        Sub105 = dZfL1(11,1,2)*ZNeuC(1,1) + 
     &    dZfL1(11,2,2)*ZNeuC(2,1) + dZfL1(11,3,2)*ZNeuC(3,1) + 
     &    dZfL1(11,4,2)*ZNeuC(4,1)

        Sub94 = dZfL1(11,1,1)*ZNeuC(1,4) + 
     &    dZfL1(11,2,1)*ZNeuC(2,4) + dZfL1(11,3,1)*ZNeuC(3,4) + 
     &    dZfL1(11,4,1)*ZNeuC(4,4)

        Sub106 = dZfL1(11,1,2)*ZNeuC(1,4) + 
     &    dZfL1(11,2,2)*ZNeuC(2,4) + dZfL1(11,3,2)*ZNeuC(3,4) + 
     &    dZfL1(11,4,2)*ZNeuC(4,4)

	Sub45 = 8*CW2**2*dZe1 - (1 - 4*CW2)*CW2*dZZZ1

	Sub61 = 2*dSW1 + SW*(-2*dZe1 - dZbarfL1(3,1,1))

	Sub99 = 2*dSW1*SW + CW2*(2*dZe1 + dZbarfR1(3,1,1))

	Sub79 = 2*dSW1*SW + CW2*(2*dZe1 + dZfL1(3,1,1))

	Sub51 = 2*dSW1*SW + CW2*dZfR1(3,1,1)

	Sub54 = dMWsq1 + MW2*(-2*dZe1 - dZfR1(3,1,1))

	Sub120 = Sub119*SW - 3*CW*(Sub61*ZNeu(1,2))

	Sub57 = MU*(SB*Sub56*SW) - 2*Sub55*ZNeu(1,4)

	Sub117 = MW2*Sub116 - MU*(SB*Sub54*SW*ZNeu(2,4))

        Sub64 = Sub11*dZfR1(11,1,2) + Sub1*dZfR1(11,2,2) + 
     &    Sub63*dZfR1(11,3,2) + Sub62*dZfR1(11,4,2)

        Sub52 = (2*dZe1 + dZfR1(11,1,1))*ZNeu(1,1) + 
     &    dZfR1(11,2,1)*ZNeu(2,1) + dZfR1(11,3,1)*ZNeu(3,1) + 
     &    dZfR1(11,4,1)*ZNeu(4,1)

	Sub60 = 2*dSW1*SW + CW2*(2*dZe1 + dZbarfL1(3,1,1))

	Sub58 = MW2*Sub57 - MU*(SB*Sub54*SW*ZNeu(1,4))

	Sub65 = Sub64*SW - 3*CW*(Sub61*ZNeu(2,2))

	Sub46 = Sub45*SW - 2*SW2*dSW1 - CW2*(6*dSW1 + 2*dZe1*SW)

	Sub43 = 2*dZe1 + dZZZ1 + dZbarfR1(3,1,1) + dZfR1(3,1,1)

	Sub44 = CW2**2*dZAZ1 - CW*(2*dSW1 + CW2*(Sub43*SW))

        Sub47 = 4*CW2**2*SW2*dZAZ1 + 
     &    CW*(Sub46 - (1 - 4*CW2)*CW2*
     &        (SW*(dZbarfL1(3,1,1) + dZfL1(3,1,1))))

        Sub180 = 1/(CW*CW2)*
     &    (3*Pi*Sub44 + 2*CW*CW2*(Alfas*Finite*SW))/SW

        Sub181 = 1/(CW*CW2)*
     &    (3*Pi*Sub47 + CW*(2 - 8*CW2)*CW2*(Alfas*Finite*SW))/SW

	pave22 = A0(MGl2)

	call Bput(pave2(bb0),S,MU2,MU2)

	call Cput(pave3(cc0),MU2,S,MU2,0.D0,MU2,MU2)

	Sub122 = MU*Sub56 + 2*dMf1(3,1)*ZNeu(1,4)

	Sub123 = Sub122*SW - 2*MU*(dSW1*ZNeu(1,4))

	Sub71 = MU*Sub70 + 2*dMf1(3,1)*ZNeu(2,4)

	Sub124 = SB*Sub123 - 2*MU*(dSB1*SW*ZNeu(1,4))

	Sub72 = Sub71*SW - 2*MU*(dSW1*ZNeu(2,4))

	Sub125 = MW2*Sub124 - MU*(SB*Sub69*SW*ZNeu(1,4))

	Sub73 = SB*Sub72 - 2*MU*(dSB1*SW*ZNeu(2,4))

	Sub74 = MW2*Sub73 - MU*(SB*Sub69*SW*ZNeu(2,4))

	Sub275 = pave3(PVC(cc0)) + pave3(PVC(cc1)) + pave3(PVC(cc2))

	LOOP(Sfe5, 1,2,1)

        Sub83(Sfe5) = dZbarSf1(1,Sfe5,3,1)*USf(1,2,3,1) + 
     &    dZbarSf1(2,Sfe5,3,1)*USf(2,2,3,1)

        Sub113(Sfe5) = 
     &   Sub111*dZbarSf1(1,Sfe5,3,1) + Sub112*dZbarSf1(2,Sfe5,3,1)

        Sub87(Sfe5) = MU*MW2*(SB*SW*Sub83(Sfe5)) - 
     &    Sub86*USf(Sfe5,2,3,1)

        Sub66(Sfe5) = dZSf1(1,Sfe5,3,1)*USfC(1,1,3,1) + 
     &    dZSf1(2,Sfe5,3,1)*USfC(2,1,3,1)

        Sub50(Sfe5) = Sub48*dZbarSf1(1,Sfe5,3,1) + 
     &    Sub49*dZbarSf1(2,Sfe5,3,1)

        Sub78(Sfe5) = dZbarSf1(1,Sfe5,3,1)*USf(1,1,3,1) + 
     &    dZbarSf1(2,Sfe5,3,1)*USf(2,1,3,1)

        Sub95(Sfe5) = MW*(SB*Sub93*USf(Sfe5,1,3,1)) + 
     &    3*CW*MU*(Sub94*USf(Sfe5,2,3,1))

        Sub67(Sfe5) = dZSf1(1,Sfe5,3,1)*USfC(1,2,3,1) + 
     &    dZSf1(2,Sfe5,3,1)*USfC(2,2,3,1)

        Sub104(Sfe5) = 
     &   MU*MW2*(SB*SW*Sub66(Sfe5)) - Sub103*USfC(Sfe5,1,3,1)

        Sub107(Sfe5) = 
     &   3*CW*MU*(Sub106*USfC(Sfe5,1,3,1)) - 
     &    4*MW*(SB*Sub105*SW*USfC(Sfe5,2,3,1))

        Sub115(Sfe5) = 
     &   CW2*(SB*SW*Sub113(Sfe5)) + 
     &    4*MW*SW2*(SB2*USf(Sfe5,2,3,1)*
     &       (CW2*Sub114 + Sub51*ZNeu(2,1)))

        Sub68(Sfe5) = MW*(SB*Sub1*Sub66(Sfe5)) + 
     &    3*CW*MU*(Sub67(Sfe5)*ZNeu(2,4))

        Sub118(Sfe5) = 
     &   MW2*Sub115(Sfe5) - 3*CW*CW2*(Sub117*USf(Sfe5,1,3,1))

	Sub80(Sfe5) = CW2*Sub78(Sfe5) + Sub79*USf(Sfe5,1,3,1)

	Sub100(Sfe5) = CW2*Sub67(Sfe5) + Sub99*USfC(Sfe5,2,3,1)

        Sub53(Sfe5) = CW2*(SB*SW*Sub50(Sfe5)) + 
     &    4*MW*SW2*(SB2*USf(Sfe5,2,3,1)*
     &       (CW2*Sub52 + Sub51*ZNeu(1,1)))

        Sub108(Sfe5) = 
     &   MW2*(SB*SW*Sub107(Sfe5)) + 3*CW*(Sub104(Sfe5)*ZNeuC(2,4))

        Sub59(Sfe5) = MW2*Sub53(Sfe5) - 
     &    3*CW*CW2*(Sub58*USf(Sfe5,1,3,1))

        Sub109(Sfe5) = 
     &   CW2*Sub108(Sfe5) - 
     &    4*MW*MW2*SW2*(SB2*Sub100(Sfe5)*ZNeuC(2,1))

        call Cput(pave14(cc0,Sfe5),S,MNeu2(2),MNeu2(1),MU2,MU2,
     &   MSf2(Sfe5,3,1))

        call Cput(pave23(cc0,Sfe5),S,MNeu2(1),MNeu2(2),MU2,MU2,
     &   MSf2(Sfe5,3,1))

	Sub82(Sfe5) = SW*Sub78(Sfe5) - Sub81*USf(Sfe5,1,3,1)

        Sub136(Sfe5) = 
     &   MW*(SB*Sub135*USf(Sfe5,1,3,1)) + 
     &    3*CW*MU*(Sub106*USf(Sfe5,2,3,1))

        Sub130(Sfe5) = 
     &   3*CW*MU*(Sub94*USfC(Sfe5,1,3,1)) - 
     &    4*MW*(SB*Sub129*SW*USfC(Sfe5,2,3,1))

        Sub121(Sfe5) = 
     &   MW*(SB*Sub11*Sub66(Sfe5)) + 
     &    3*CW*MU*(Sub67(Sfe5)*ZNeu(1,4))

        Sub134(Sfe5) = 
     &   MW*MW2*(SB2*Sub82(Sfe5)*ZNeuC(2,2)) + 
     &    Sub87(Sfe5)*ZNeuC(2,4)

	Sub137(Sfe5) = 3*CW*Sub134(Sfe5) + MW2*(SB*SW*Sub136(Sfe5))

        Sub126(Sfe5) = 
     &   MW2*(SB*SW*Sub121(Sfe5)) + 3*CW*(Sub125*USfC(Sfe5,2,3,1))

        Sub131(Sfe5) = 
     &   MW2*(SB*SW*Sub130(Sfe5)) + 3*CW*(Sub104(Sfe5)*ZNeuC(1,4))

        Sub88(Sfe5) = MW*MW2*(SB2*Sub82(Sfe5)*ZNeuC(1,2)) + 
     &    Sub87(Sfe5)*ZNeuC(1,4)

	Sub96(Sfe5) = 3*CW*Sub88(Sfe5) + MW2*(SB*SW*Sub95(Sfe5))

        Sub75(Sfe5) = MW2*(SB*SW*Sub68(Sfe5)) + 
     &    3*CW*(Sub74*USfC(Sfe5,2,3,1))

        Sub127(Sfe5) = 
     &   CW2*Sub126(Sfe5) + 
     &    MW*MW2*(SB2*(CW2*Sub120 + SW2*(Sub60*ZNeu(1,1)))*
     &       USfC(Sfe5,1,3,1))

        Sub132(Sfe5) = 
     &   CW2*Sub131(Sfe5) - 
     &    4*MW*MW2*SW2*(SB2*Sub100(Sfe5)*ZNeuC(1,1))

        Sub138(Sfe5) = 
     &   CW2*Sub137(Sfe5) + MW*MW2*SW2*(SB2*Sub80(Sfe5)*ZNeuC(2,1))

        Sub128(Sfe5) = 
     &   Sub118(Sfe5)*Sub12(Sfe5) + Sub10(Sfe5)*Sub127(Sfe5)

        Sub133(Sfe5) = 
     &   -(Sub118(Sfe5)*Sub13(Sfe5)) + Sub10(Sfe5)*Sub132(Sfe5)

        Sub76(Sfe5) = CW2*Sub75(Sfe5) + 
     &    MW*MW2*(SB2*(CW2*Sub65 + SW2*(Sub60*ZNeu(2,1)))*
     &       USfC(Sfe5,1,3,1))

        Sub97(Sfe5) = CW2*Sub96(Sfe5) + 
     &    MW*MW2*SW2*(SB2*Sub80(Sfe5)*ZNeuC(1,1))

	Sub77(Sfe5) = Sub2(Sfe5)*Sub59(Sfe5) + Sub3(Sfe5)*Sub76(Sfe5)

	ENDLOOP(Sfe5)

	LOOP(Sfe6, 1,2,1)

        Sub147(Sfe6) = 
     &   MW*(SB*USf(Sfe6,1,3,1)*ZNeuC(1,2)) + 
     &    MU*(USf(Sfe6,2,3,1)*ZNeuC(1,4))

        Sub159(Sfe6) = 
     &   MW*(SB*USf(Sfe6,1,3,1)*ZNeuC(2,2)) + 
     &    MU*(USf(Sfe6,2,3,1)*ZNeuC(2,4))

        Sub145(Sfe6) = 
     &   4*MW*(SB*SW*USf(Sfe6,2,3,1)*ZNeu(1,1)) - 
     &    3*CW*MU*(USf(Sfe6,1,3,1)*ZNeu(1,4))

        Sub154(Sfe6) = 
     &   4*MW*(SB*SW*USf(Sfe6,2,3,1)*ZNeu(2,1)) - 
     &    3*CW*MU*(USf(Sfe6,1,3,1)*ZNeu(2,4))

        Sub148(Sfe6) = 
     &   3*CW*Sub147(Sfe6) + MW*(SB*SW*USf(Sfe6,1,3,1)*ZNeuC(1,1))

        Sub160(Sfe6) = 
     &   3*CW*Sub159(Sfe6) + MW*(SB*SW*USf(Sfe6,1,3,1)*ZNeuC(2,1))

        Sub28(Sfe6) = MW*(SB*Sub11*USfC(Sfe6,1,3,1)) + 
     &    3*CW*MU*(USfC(Sfe6,2,3,1)*ZNeu(1,4))

        Sub23(Sfe6) = MW*(SB*Sub1*USfC(Sfe6,1,3,1)) + 
     &    3*CW*MU*(USfC(Sfe6,2,3,1)*ZNeu(2,4))

        Sub29(Sfe6) = 4*MW*
     &     (SB*SW*USfC(Sfe6,2,3,1)*ZNeuC(1,1)) - 
     &    3*CW*MU*(USfC(Sfe6,1,3,1)*ZNeuC(1,4))

        Sub24(Sfe6) = 4*MW*
     &     (SB*SW*USfC(Sfe6,2,3,1)*ZNeuC(2,1)) - 
     &    3*CW*MU*(USfC(Sfe6,1,3,1)*ZNeuC(2,4))

	ENDLOOP(Sfe6)

	LOOP(Sfe7, 1,2,1)

	pave1(Sfe7) = A0(MSf2(Sfe7,3,1))

	ENDLOOP(Sfe7)

	LOOP(Sfe6, 1,2,1)
	LOOP(Sfe5, 1,2,1)

        Sub16(Sfe5,Sfe6) = 
     &   dZbarSf1(Sfe5,Sfe6,3,1) + dZSf1(Sfe6,Sfe5,3,1)

        Sub17(Sfe5,Sfe6) = 
     &   2*dMSfsq1(Sfe6,Sfe5,3,1) + 
     &    dZbarSf1(Sfe5,Sfe6,3,1)*MSf2(Sfe5,3,1) + 
     &    dZSf1(Sfe6,Sfe5,3,1)*MSf2(Sfe6,3,1)

        Sub32(Sfe5,Sfe6) = 
     &   USf(Sfe6,1,3,1)*USfC(Sfe5,1,3,1) + 
     &    USf(Sfe6,2,3,1)*USfC(Sfe5,2,3,1)

        call Cput(pave17(cc0,Sfe5,Sfe6),MU2,S,MU2,MGl2,
     &   MSf2(Sfe5,3,1),MSf2(Sfe6,3,1))

        Sub141(Sfe5,Sfe6) = 
     &   4*SW2*(USf(Sfe6,2,3,1)*USfC(Sfe5,2,3,1)) + 
     &    USf(Sfe6,1,3,1)*(-(4*CW2*USfC(Sfe5,1,3,1)) + 
     &       USfC(Sfe5,1,3,1))

        Sub146(Sfe5,Sfe6) = 
     &   Sub12(Sfe5)*Sub23(Sfe6)*USf(Sfe5,1,3,1)*USf(Sfe6,2,3,1) + 
     &    Sub10(Sfe5)*Sub145(Sfe6)*USfC(Sfe5,2,3,1)*
     &     USfC(Sfe6,1,3,1)

        Sub155(Sfe5,Sfe6) = 
     &   Sub2(Sfe5)*Sub28(Sfe6)*USf(Sfe5,1,3,1)*USf(Sfe6,2,3,1) + 
     &    Sub154(Sfe6)*Sub3(Sfe5)*USfC(Sfe5,2,3,1)*USfC(Sfe6,1,3,1)

        Sub156(Sfe5,Sfe6) = 
     &   Sub2(Sfe5)*Sub29(Sfe6)*USf(Sfe5,2,3,1)*USf(Sfe6,1,3,1) + 
     &    Sub154(Sfe6)*Sub5(Sfe5)*USfC(Sfe5,2,3,1)*USfC(Sfe6,1,3,1)

        Sub149(Sfe5,Sfe6) = 
     &   Sub13(Sfe5)*Sub23(Sfe6)*USf(Sfe5,1,3,1)*USf(Sfe6,2,3,1) + 
     &    Sub10(Sfe5)*Sub148(Sfe6)*USfC(Sfe5,1,3,1)*
     &     USfC(Sfe6,2,3,1)

	ENDLOOP(Sfe5)
	ENDLOOP(Sfe6)

	LOOP(Sfe7, 1,2,1)
	LOOP(Sfe5, 1,2,1)

        Sub35(Sfe5,Sfe7) = 
     &   USf(Sfe7,1,3,1)*USfC(Sfe5,1,3,1) - 
     &    USf(Sfe7,2,3,1)*USfC(Sfe5,2,3,1)

	ENDLOOP(Sfe5)
	ENDLOOP(Sfe7)

	LOOP(Sfe7, 1,2,1)
	LOOP(Sfe6, 1,2,1)

        Sub36(Sfe6,Sfe7) = 
     &   USf(Sfe6,1,3,1)*USfC(Sfe7,1,3,1) - 
     &    USf(Sfe6,2,3,1)*USfC(Sfe7,2,3,1)

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

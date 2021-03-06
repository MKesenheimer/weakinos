#if 0
* ddbar_n1n1_abbr1h.F
* abbreviations for {F[4, {1, Col1}], -F[4, {1, Col2}]} -> {F[11, {1}], F[11, {1}]}
* generated by FormCalc 8.4 on 5-Mar-2015 14:01
#endif

#include "ddbar_n1n1_vars.h"
#include "contains.h"


	subroutine ddbar_n1n1_abbr1h
	implicit none

#include "ddbar_n1n1_vars.h"


	F35 = ChainV1(Spinor(2,-1,2),0,k(1),1,Spinor(4,1,2))

	F36 = -ChainB1(Spinor(2,-1,2),0,k(1),1,Spinor(4,1,2))

	F29 = -ChainB1(Spinor(3,1,2),0,k(1),0,Spinor(4,-1,1))

	F20 = -ChainV0(Spinor(3,1,2),0,0,Spinor(4,-1,1))

	F19 = ChainB0(Spinor(3,1,2),0,0,Spinor(4,-1,1))

	F39 = ChainV1(Spinor(2,-1,2),0,k(1),0,Spinor(3,-1,1))

	F41 = ChainV1(Spinor(2,-1,2),0,k(1),1,Spinor(3,1,2))

	F37 = ChainV1(Spinor(2,-1,2),0,k(1),0,Spinor(4,-1,1))

	F38 = ChainB1(Spinor(2,-1,2),0,k(1),0,Spinor(3,-1,1))

	F42 = -ChainB1(Spinor(2,-1,2),0,k(1),1,Spinor(3,1,2))

	F40 = ChainB1(Spinor(2,-1,2),0,k(1),0,Spinor(4,-1,1))

	F32 = ChainV1(Spinor(3,1,2),0,k(1),0,Spinor(4,-1,1))

	F33 = ChainV1(Spinor(3,1,2),0,k(2),0,Spinor(1,1,1))

	F34 = ChainB1(Spinor(3,1,2),0,k(2),0,Spinor(1,1,1))

	F30 = ChainV1(Spinor(4,1,2),0,k(2),0,Spinor(1,1,1))

	F31 = ChainB1(Spinor(4,1,2),0,k(2),0,Spinor(1,1,1))

	F16 = ChainV0(Spinor(2,-1,2),0,0,Spinor(1,1,1))

	F21 = ChainV0(Spinor(2,-1,2),0,1,Spinor(3,1,2))

	F13 = -ChainB0(Spinor(2,-1,2),0,0,Spinor(1,1,1))

	F24 = -ChainB0(Spinor(2,-1,2),0,1,Spinor(3,1,2))

	F27 = -ChainV0(Spinor(3,-1,1),1,0,Spinor(1,1,1))

	F23 = ChainB0(Spinor(3,-1,1),1,0,Spinor(1,1,1))

	F25 = -ChainV0(Spinor(4,-1,1),1,0,Spinor(1,1,1))

	F22 = ChainB0(Spinor(4,-1,1),1,0,Spinor(1,1,1))

	F14 = -ChainV0(Spinor(4,1,2),0,0,Spinor(3,-1,1))

	F15 = ChainB0(Spinor(4,1,2),0,0,Spinor(3,-1,1))

	F17 = ChainV1(Spinor(4,1,2),0,k(1),0,Spinor(3,-1,1))

	F18 = -ChainB1(Spinor(4,1,2),0,k(1),0,Spinor(3,-1,1))

	F26 = ChainV1(Spinor(2,-1,2),0,k(3),0,Spinor(1,1,1))

	F28 = -ChainB1(Spinor(2,-1,2),0,k(3),0,Spinor(1,1,1))

	Sub147 = F11*F12 + F13*(F19 - F20)

	Sub155 = F13*(F14 - F15) - F1*F5

	Sub152 = (F14 - F15)*F16 + F2*F3

	Sub144 = F16*(F19 - F20) - F10*F9

	Sub148 = (F13*F29 - F12*F33)*MNeu(1) + Sub147*MNeu2(1)

	Sub156 = -((F13*F18 - F1*F30)*MNeu(1)) + Sub155*MNeu2(1)

	Sub178 = F2*F38*MNeu(1) + F1*(F2*Sub177 - F31*MNeu(1))

	Sub180 = -(F3*F31*MNeu(1)) + F2*(F3*Sub177 + F39*MNeu(1))

	Sub188 = -(F12*F34*MNeu(1)) + F10*(F12*Sub177 + F40*MNeu(1))

	Sub189 = -(F12*F33*MNeu(1)) + F11*(F12*Sub177 + F40*MNeu(1))

	Sub179 = F38*F5*MNeu(1) + F1*(F5*Sub177 - F30*MNeu(1))

	Sub181 = F39*F5*MNeu(1) + F3*(F5*Sub177 - F30*MNeu(1))

	Sub191 = -(F33*F9*MNeu(1)) + F11*(F9*Sub177 + F37*MNeu(1))

	Sub190 = -(F34*F9*MNeu(1)) + F10*(F9*Sub177 + F37*MNeu(1))

	Sub131 = F17 - F18 + (F14 - F15)*MNeu(1)

	Sub132 = F17 - F18 + (2*F14 - 2*F15)*MNeu(1)

        Sub123 = -(1/SW2*
     &      ((F3*F5 - F6*F7)*Sub122 + 
     &        SW2*((-2*F1*F2 + 2*F4*F8)*Sub121)))

	Sub173 = F38*MNeu(1) - F1*(T - MNeu2(1))

	Sub157 = F31*MNeu(1) + F2*(T - MNeu2(1))

	Sub171 = F39*MNeu(1) - F3*(T - MNeu2(1))

	Sub170 = (F27*F35 - F2*F39)*MNeu(1) + F2*F3*(T - MNeu2(1))

	Sub174 = (F23*F36 - F38*F5)*MNeu(1) + F1*F5*(T - MNeu2(1))

	Sub149 = F34*MNeu(1) + F10*(U - MNeu2(1))

	Sub154 = F30*MNeu(1) + F5*(T - MNeu2(1))

	Sub146 = F33*MNeu(1) + F11*(U - MNeu2(1))

	Sub165 = F40*MNeu(1) - F12*(U - MNeu2(1))

	Sub166 = (F11*F40 - F22*F42)*MNeu(1) - F11*F12*(U - MNeu2(1))

	Sub168 = F37*MNeu(1) - F9*(U - MNeu2(1))

	Sub167 = (F10*F37 - F25*F41)*MNeu(1) - F10*F9*(U - MNeu2(1))

        Sub161 = -((F12*F33 - F11*F40)*MNeu(1)) - 
     &    F11*F12*(U - MNeu2(1)) + (F11*F12 - F13*F20)*MNeu2(1)

        Sub160 = -((F12*F34 - F10*F40)*MNeu(1)) - 
     &    F10*F12*(U - MNeu2(1)) + (F10*F12 - F24*F25)*MNeu2(1)

        Sub185 = (F3*F31 - F2*F39)*MNeu(1) + 
     &    F2*F3*(T - MNeu2(1)) + (F15*F16 - F2*F3)*MNeu2(1)

        Sub184 = (F1*F30 - F38*F5)*MNeu(1) + 
     &    F1*F5*(T - MNeu2(1)) + (F13*F14 - F1*F5)*MNeu2(1)

        Sub186 = (F3*F30 - F39*F5)*MNeu(1) + 
     &    F3*F5*(T - MNeu2(1)) - (F3*F5 - F23*F7)*MNeu2(1)

        Sub183 = (F1*F31 - F2*F38)*MNeu(1) + 
     &    F1*F2*(T - MNeu2(1)) - (F1*F2 - F27*F8)*MNeu2(1)

        Sub162 = (F10*F37 - F34*F9)*MNeu(1) - 
     &    F10*F9*(U - MNeu2(1)) - (F16*F19 - F10*F9)*MNeu2(1)

        Sub163 = (F11*F37 - F33*F9)*MNeu(1) - 
     &    F11*F9*(U - MNeu2(1)) - (F21*F22 - F11*F9)*MNeu2(1)

        Sub145 = F10*F9*U + (-(F16*F32) + F34*F9)*MNeu(1) + 
     &    Sub144*MNeu2(1)

        Sub153 = F2*F3*T + (-(F16*F17) + F3*F31)*MNeu(1) - 
     &    Sub152*MNeu2(1)

	LOOP(Sfe5, 1,2,1)

        Sub136(HelInd(Sfe5)) = 
     &   (F10*F12 - F1*F2)*Sub1(Sfe5) + 
     &    (F11*F12 - F1*F5)*Sub3(Sfe5)

        Sub137(HelInd(Sfe5)) = 
     &   (F2*F3 - F10*F9)*Sub1(Sfe5) + (F3*F5 - F11*F9)*Sub3(Sfe5)

        Sub68(HelInd(Sfe5)) = 
     &   F11*Sub3(Sfe5)*Sub66(Sfe5) + F10*Sub67(Sfe5)

        Sub100(HelInd(Sfe5)) = 
     &   F5*Sub3(Sfe5)*Sub66(Sfe5) + F2*Sub67(Sfe5)

        Sub98(HelInd(Sfe5)) = 
     &   F10*(MW2*Sub42(Sfe5) + CW2*(CB*Sub49(Sfe5)))*Sub6(Sfe5) + 
     &    (F10*Sub1(Sfe5) + F11*Sub3(Sfe5))*Sub97(Sfe5)

        Sub101(HelInd(Sfe5)) = 
     &   F2*(MW2*Sub42(Sfe5) + CW2*(CB*Sub49(Sfe5)))*Sub6(Sfe5) + 
     &    (F2*Sub1(Sfe5) + F5*Sub3(Sfe5))*Sub97(Sfe5)

        Sub133(HelInd(Sfe5)) = 
     &   F10*F12*Sub1(Sfe5) + 2*F22*F24*Sub3(Sfe5)

	Sub139(HelInd(Sfe5)) = 2*F27*F7*Sub1(Sfe5) + F3*F5*Sub3(Sfe5)

	Sub140(HelInd(Sfe5)) = F1*F2*Sub1(Sfe5) + 2*F23*F8*Sub3(Sfe5)

        Sub134(HelInd(Sfe5)) = 
     &   2*F21*F25*Sub1(Sfe5) + F11*F9*Sub3(Sfe5)

        Sub138(HelInd(Sfe5)) = 
     &   Sub137(HelInd(Sfe5))*Sub5(Sfe5) - 
     &    Sub136(HelInd(Sfe5))*Sub6(Sfe5)

        Sub142(HelInd(Sfe5)) = 
     &   F21*F25*Sub1(Sfe5)*Sub5(Sfe5) + 
     &    F22*F24*Sub3(Sfe5)*Sub6(Sfe5)

        Sub143(HelInd(Sfe5)) = 
     &   F27*F7*Sub1(Sfe5)*Sub5(Sfe5) + 
     &    F23*F8*Sub3(Sfe5)*Sub6(Sfe5)

        Sub182(HelInd(Sfe5)) = 
     &   (Sub180*Sub1(Sfe5) + Sub181*Sub3(Sfe5))*Sub5(Sfe5) + 
     &    (Sub178*Sub1(Sfe5) + Sub179*Sub3(Sfe5))*Sub6(Sfe5)

        Sub192(HelInd(Sfe5)) = 
     &   (Sub190*Sub1(Sfe5) + Sub191*Sub3(Sfe5))*Sub5(Sfe5) + 
     &    (Sub188*Sub1(Sfe5) + Sub189*Sub3(Sfe5))*Sub6(Sfe5)

        Sub102(HelInd(Sfe5)) = 
     &   F3*Sub100(HelInd(Sfe5)) + F1*Sub101(HelInd(Sfe5)) - 
     &    F5*(F3*Sub5(Sfe5) + F1*Sub6(Sfe5))*Sub87(Sfe5)

        Sub99(HelInd(Sfe5)) = 
     &   -(F9*Sub68(HelInd(Sfe5))) + 
     &    F11*(F9*Sub5(Sfe5) + F12*Sub6(Sfe5))*Sub87(Sfe5) - 
     &    F12*Sub98(HelInd(Sfe5))

        Sub135(HelInd(Sfe5)) = 
     &   Sub134(HelInd(Sfe5))*Sub5(Sfe5) + 
     &    Sub133(HelInd(Sfe5))*Sub6(Sfe5)

        Sub141(HelInd(Sfe5)) = 
     &   Sub139(HelInd(Sfe5))*Sub5(Sfe5) + 
     &    Sub140(HelInd(Sfe5))*Sub6(Sfe5)

        Sub225(HelInd(Sfe5)) = 
     &   3*Pi*Sub102(HelInd(Sfe5))/(CB*SW) + 
     &    4*CW2*MW2*(Alfas*(Sub223(Sfe5) + 2*Sub224(Sfe5))*
     &       (F2*Sub1(Sfe5) + F5*Sub3(Sfe5))*
     &       (F3*Sub5(Sfe5) + F1*Sub6(Sfe5)))

        Sub241(HelInd(Sfe5)) = 
     &   pave6(PVC(dd12),Sfe5)*Sub142(HelInd(Sfe5)) + 
     &    Sub240(Sfe5)*(F2*Sub1(Sfe5) + F5*Sub3(Sfe5))*
     &     (F3*Sub5(Sfe5) + F1*Sub6(Sfe5))

        Sub243(HelInd(Sfe5)) = 
     &   pave7(PVC(dd12),Sfe5)*Sub143(HelInd(Sfe5)) + 
     &    Sub242(Sfe5)*(F10*Sub1(Sfe5) + F11*Sub3(Sfe5))*
     &     (F9*Sub5(Sfe5) + F12*Sub6(Sfe5))

        Sub228(HelInd(Sfe5)) = 
     &   -(4*CW2*MW2*(Alfas*(Sub226(Sfe5) + 2*Sub227(Sfe5))*
     &         (F10*Sub1(Sfe5) + F11*Sub3(Sfe5))*
     &         (F9*Sub5(Sfe5) + F12*Sub6(Sfe5)))) + 
     &    3*Pi*Sub99(HelInd(Sfe5))/(CB*SW)

        Sub158(HelInd(Sfe5)) = 
     &   Sub156*Sub3(Sfe5) + 
     &    F1*(Sub157*Sub1(Sfe5) + F5*T*Sub3(Sfe5))

        Sub150(HelInd(Sfe5)) = 
     &   -(Sub148*Sub3(Sfe5)) + 
     &    F12*(Sub149*Sub1(Sfe5) + F11*U*Sub3(Sfe5))

        Sub172(HelInd(Sfe5)) = 
     &   Sub170*Sub1(Sfe5) - F5*Sub171*Sub3(Sfe5)

        Sub175(HelInd(Sfe5)) = 
     &   F2*Sub173*Sub1(Sfe5) - Sub174*Sub3(Sfe5)

        Sub151(HelInd(Sfe5)) = 
     &   (Sub145*Sub1(Sfe5) + F9*Sub146*Sub3(Sfe5))*Sub5(Sfe5) + 
     &    Sub150(HelInd(Sfe5))*Sub6(Sfe5)

        Sub159(HelInd(Sfe5)) = 
     &   (Sub153*Sub1(Sfe5) + F3*Sub154*Sub3(Sfe5))*Sub5(Sfe5) + 
     &    Sub158(HelInd(Sfe5))*Sub6(Sfe5)

        Sub176(HelInd(Sfe5)) = 
     &   Sub172(HelInd(Sfe5))*Sub5(Sfe5) - 
     &    Sub175(HelInd(Sfe5))*Sub6(Sfe5)

        Sub164(HelInd(Sfe5)) = 
     &   (Sub162*Sub1(Sfe5) + Sub163*Sub3(Sfe5))*Sub5(Sfe5) + 
     &    (Sub160*Sub1(Sfe5) + Sub161*Sub3(Sfe5))*Sub6(Sfe5)

        Sub169(HelInd(Sfe5)) = 
     &   (Sub167*Sub1(Sfe5) + F11*Sub168*Sub3(Sfe5))*Sub5(Sfe5) + 
     &    (F10*Sub165*Sub1(Sfe5) + Sub166*Sub3(Sfe5))*Sub6(Sfe5)

        Sub187(HelInd(Sfe5)) = 
     &   (Sub185*Sub1(Sfe5) + Sub186*Sub3(Sfe5))*Sub5(Sfe5) + 
     &    (Sub183*Sub1(Sfe5) + Sub184*Sub3(Sfe5))*Sub6(Sfe5)

        Sub238(HelInd(Sfe5)) = 
     &   pave6(PVC(dd13),Sfe5)*Sub151(HelInd(Sfe5)) + 
     &    pave7(PVC(dd3),Sfe5)*Sub182(HelInd(Sfe5))

        Sub239(HelInd(Sfe5)) = 
     &   pave18(PVC(cc0),Sfe5)*Sub138(HelInd(Sfe5)) + 
     &    pave7(PVC(dd13),Sfe5)*Sub159(HelInd(Sfe5)) + 
     &    pave6(PVC(dd33),Sfe5)*Sub164(HelInd(Sfe5)) + 
     &    pave6(PVC(dd23),Sfe5)*Sub169(HelInd(Sfe5)) + 
     &    pave7(PVC(dd23),Sfe5)*Sub176(HelInd(Sfe5)) + 
     &    pave7(PVC(dd33),Sfe5)*Sub187(HelInd(Sfe5)) + 
     &    pave6(PVC(dd3),Sfe5)*Sub192(HelInd(Sfe5))

	ENDLOOP(Sfe5)

	LOOP(Sfe6, 1,2,1)
	LOOP(Sfe5, 1,2,1)

        Sub117(HelInd(Sfe5),Sfe6) = 
     &   F1*Sub116(Sfe5,Sfe6) + 
     &    F3*Sub14(Sfe6)*Sub5(Sfe5)*USf(Sfe5,1,4,1)*USf(Sfe6,1,4,1)

        Sub126(HelInd(Sfe5),Sfe6) = 
     &   F12*Sub116(Sfe5,Sfe6) + 
     &    F9*Sub14(Sfe6)*Sub5(Sfe5)*USf(Sfe5,1,4,1)*USf(Sfe6,1,4,1)

        Sub124(HelInd(Sfe5),Sfe6) = 
     &   F10*Sub15(Sfe6)*Sub5(Sfe5)*USf(Sfe5,2,4,1)*
     &     USf(Sfe6,2,4,1) + 
     &    Sub108(Sfe6)*(F10*Sub1(Sfe5) + F11*Sub3(Sfe5))*
     &     USfC(Sfe5,1,4,1)*USfC(Sfe6,1,4,1)

        Sub114(HelInd(Sfe5),Sfe6) = 
     &   F2*Sub15(Sfe6)*Sub5(Sfe5)*USf(Sfe5,2,4,1)*
     &     USf(Sfe6,2,4,1) + 
     &    Sub108(Sfe6)*(F2*Sub1(Sfe5) + F5*Sub3(Sfe5))*
     &     USfC(Sfe5,1,4,1)*USfC(Sfe6,1,4,1)

        Sub125(HelInd(Sfe5),Sfe6) = 
     &   F9*Sub124(HelInd(Sfe5),Sfe6) + 
     &    F10*F12*Sub15(Sfe6)*Sub6(Sfe5)*USf(Sfe5,2,4,1)*
     &     USf(Sfe6,2,4,1)

        Sub115(HelInd(Sfe5),Sfe6) = 
     &   F3*Sub114(HelInd(Sfe5),Sfe6) + 
     &    F1*F2*Sub15(Sfe6)*Sub6(Sfe5)*USf(Sfe5,2,4,1)*
     &     USf(Sfe6,2,4,1)

        Sub127(HelInd(Sfe5),Sfe6) = 
     &   F11*Sub126(HelInd(Sfe5),Sfe6) + 
     &    F10*F12*Sub1(Sfe5)*Sub105(Sfe6)*USfC(Sfe5,2,4,1)*
     &     USfC(Sfe6,2,4,1)

        Sub118(HelInd(Sfe5),Sfe6) = 
     &   F5*Sub117(HelInd(Sfe5),Sfe6) + 
     &    F1*F2*Sub1(Sfe5)*Sub105(Sfe6)*USfC(Sfe5,2,4,1)*
     &     USfC(Sfe6,2,4,1)

        Sub119(HelInd(Sfe5),Sfe6) = 
     &   Mino3*Sub115(HelInd(Sfe5),Sfe6) + 
     &    Mino3C*Sub118(HelInd(Sfe5),Sfe6)

        Sub128(HelInd(Sfe5),Sfe6) = 
     &   Mino3*Sub125(HelInd(Sfe5),Sfe6) + 
     &    Mino3C*Sub127(HelInd(Sfe5),Sfe6)

        Sub107(HelInd(Sfe5),Sfe6) = 
     &   F10*Sub106(Sfe5,Sfe6) + 
     &    F11*Sub105(Sfe6)*Sub3(Sfe5)*USfC(Sfe5,2,4,1)*
     &     USfC(Sfe6,1,4,1)

        Sub111(HelInd(Sfe5),Sfe6) = 
     &   F2*Sub106(Sfe5,Sfe6) + 
     &    F5*Sub105(Sfe6)*Sub3(Sfe5)*USfC(Sfe5,2,4,1)*
     &     USfC(Sfe6,1,4,1)

        Sub109(HelInd(Sfe5),Sfe6) = 
     &   F10*Sub14(Sfe6)*Sub6(Sfe5)*USf(Sfe5,1,4,1)*
     &     USf(Sfe6,2,4,1) + 
     &    Sub108(Sfe6)*(F10*Sub1(Sfe5) + F11*Sub3(Sfe5))*
     &     USfC(Sfe5,1,4,1)*USfC(Sfe6,2,4,1)

        Sub112(HelInd(Sfe5),Sfe6) = 
     &   F2*Sub14(Sfe6)*Sub6(Sfe5)*USf(Sfe5,1,4,1)*
     &     USf(Sfe6,2,4,1) + 
     &    Sub108(Sfe6)*(F2*Sub1(Sfe5) + F5*Sub3(Sfe5))*
     &     USfC(Sfe5,1,4,1)*USfC(Sfe6,2,4,1)

        Sub201(HelInd(Sfe5),Sfe6) = 
     &   (F17 - 2*F15*MNeu(1))*Sub108(Sfe6)*Sub5(Sfe5) + 
     &    (F18 - 2*F14*MNeu(1))*Sub105(Sfe6)*Sub6(Sfe5)

        Sub198(HelInd(Sfe5),Sfe6) = 
     &   (F18 + F15*MNeu(1))*Sub108(Sfe6)*Sub5(Sfe5) + 
     &    (F17 + F14*MNeu(1))*Sub105(Sfe6)*Sub6(Sfe5)

        Sub200(HelInd(Sfe5),Sfe6) = 
     &   (F18 + 2*F15*MNeu(1))*Sub108(Sfe6)*Sub5(Sfe5) + 
     &    (F17 + 2*F14*MNeu(1))*Sub105(Sfe6)*Sub6(Sfe5)

        Sub202(HelInd(Sfe5),Sfe6) = 
     &   (F17 - (F14 + F15)*MNeu(1))*Sub108(Sfe6)*Sub5(Sfe5) + 
     &    (F18 - (F14 + F15)*MNeu(1))*Sub105(Sfe6)*Sub6(Sfe5)

        Sub203(HelInd(Sfe5),Sfe6) = 
     &   (F18 + (F14 + F15)*MNeu(1))*Sub108(Sfe6)*Sub5(Sfe5) + 
     &    (F17 + (F14 + F15)*MNeu(1))*Sub105(Sfe6)*Sub6(Sfe5)

        Sub113(HelInd(Sfe5),Sfe6) = 
     &   F3*Sub111(HelInd(Sfe5),Sfe6) + 
     &    F1*Sub112(HelInd(Sfe5),Sfe6) + 
     &    F5*Sub15(Sfe6)*(F3*Sub5(Sfe5) + F1*Sub6(Sfe5))*
     &     USf(Sfe5,2,4,1)*USf(Sfe6,1,4,1)

        Sub110(HelInd(Sfe5),Sfe6) = 
     &   F9*Sub107(HelInd(Sfe5),Sfe6) + 
     &    F12*Sub109(HelInd(Sfe5),Sfe6) + 
     &    F11*Sub15(Sfe6)*(F9*Sub5(Sfe5) + F12*Sub6(Sfe5))*
     &     USf(Sfe5,2,4,1)*USf(Sfe6,1,4,1)

        Sub204(HelInd(Sfe5),Sfe6) = 
     &   F15*Sub108(Sfe6)*Sub5(Sfe5) + F14*Sub105(Sfe6)*Sub6(Sfe5)

        Sub206(HelInd(Sfe5),Sfe6) = 
     &   F14*Sub108(Sfe6)*Sub5(Sfe5) + F15*Sub105(Sfe6)*Sub6(Sfe5)

        Sub199(HelInd(Sfe5),Sfe6) = 
     &   (F17 - F15*MNeu(1))*Sub108(Sfe6)*Sub5(Sfe5) + 
     &    (F18 - F14*MNeu(1))*Sub105(Sfe6)*Sub6(Sfe5)

        Sub130(HelInd(Sfe5),Sfe6) = 
     &   F16*Mino3*USf(Sfe5,2,4,1)*USfC(Sfe6,1,4,1) + 
     &    F13*Mino3C*USf(Sfe5,1,4,1)*USfC(Sfe6,2,4,1)

        Sub194(HelInd(Sfe5),Sfe6) = 
     &   F6*F7*USf(Sfe5,1,4,1)*USfC(Sfe6,1,4,1) + 
     &    F1*F2*USf(Sfe5,2,4,1)*USfC(Sfe6,2,4,1)

        Sub193(HelInd(Sfe5),Sfe6) = 
     &   F3*F5*USf(Sfe5,1,4,1)*USfC(Sfe6,1,4,1) + 
     &    F4*F8*USf(Sfe5,2,4,1)*USfC(Sfe6,2,4,1)

        Sub104(HelInd(Sfe5),Sfe6) = 
     &   (F3*F5 - F6*F7)*USf(Sfe5,1,4,1)*USfC(Sfe6,1,4,1) + 
     &    (-(F1*F2) + F4*F8)*USf(Sfe5,2,4,1)*USfC(Sfe6,2,4,1)

        Sub220(HelInd(Sfe5),Sfe6) = 
     &   Sub131*pave16(PVC(cc0),Sfe5,Sfe6) + 
     &    (F17 - F18)*pave16(PVC(cc1),Sfe5,Sfe6) + 
     &    Sub132*pave16(PVC(cc2),Sfe5,Sfe6)

        Sub232(HelInd(Sfe5),Sfe6) = 
     &   pave5(PVC(dd22),Sfe5,Sfe6)*Sub204(HelInd(Sfe5),Sfe6) - 
     &    (F14 + F15)*Sub205(Sfe5,Sfe6)*Sub230(Sfe5,Sfe6) - 
     &    Sub206(HelInd(Sfe5),Sfe6)*Sub231(Sfe5,Sfe6)

        Sub229(HelInd(Sfe5),Sfe6) = 
     &   pave4(PVC(dd11),Sfe5,Sfe6)*Sub198(HelInd(Sfe5),Sfe6) + 
     &    pave5(PVC(dd11),Sfe5,Sfe6)*Sub199(HelInd(Sfe5),Sfe6) + 
     &    pave4(PVC(dd12),Sfe5,Sfe6)*Sub200(HelInd(Sfe5),Sfe6) + 
     &    pave5(PVC(dd12),Sfe5,Sfe6)*Sub201(HelInd(Sfe5),Sfe6) + 
     &    pave5(PVC(dd13),Sfe5,Sfe6)*Sub202(HelInd(Sfe5),Sfe6) + 
     &    pave4(PVC(dd13),Sfe5,Sfe6)*Sub203(HelInd(Sfe5),Sfe6) + 
     &    MNeu(1)*pave4(PVC(dd22),Sfe5,Sfe6)*
     &     Sub204(HelInd(Sfe5),Sfe6)

        Sub234(HelInd(Sfe5),Sfe6) = 
     &   pave5(PVC(dd2),Sfe5,Sfe6)*Sub204(HelInd(Sfe5),Sfe6) + 
     &    pave5(PVC(dd3),Sfe5,Sfe6)*Sub206(HelInd(Sfe5),Sfe6)

        Sub221(HelInd(Sfe5),Sfe6) = 
     &   2*pave16(PVC(cc00),Sfe5,Sfe6)*Sub104(HelInd(Sfe5),Sfe6) - 
     &    Sub130(HelInd(Sfe5),Sfe6)*Sub220(HelInd(Sfe5),Sfe6)

        Sub233(HelInd(Sfe5),Sfe6) = 
     &   Sub229(HelInd(Sfe5),Sfe6) - 
     &    MNeu(1)*Sub232(HelInd(Sfe5),Sfe6)

        Sub196(HelInd(Sfe5),Sfe6) = 
     &   Sub108(Sfe6)*Sub194(HelInd(Sfe5),Sfe6)*Sub5(Sfe5) + 
     &    Sub105(Sfe6)*Sub193(HelInd(Sfe5),Sfe6)*Sub6(Sfe5)

        Sub195(HelInd(Sfe5),Sfe6) = 
     &   Sub108(Sfe6)*Sub193(HelInd(Sfe5),Sfe6)*Sub5(Sfe5) + 
     &    Sub105(Sfe6)*Sub194(HelInd(Sfe5),Sfe6)*Sub6(Sfe5)

        Sub197(HelInd(Sfe5),Sfe6) = 
     &   F26*USf(Sfe5,1,4,1)*USfC(Sfe6,1,4,1) + 
     &    F28*USf(Sfe5,2,4,1)*USfC(Sfe6,2,4,1)

        Sub208(HelInd(Sfe5),Sfe6) = 
     &   -Sub130(HelInd(Sfe5),Sfe6) + Sub197(HelInd(Sfe5),Sfe6)

        Sub120(HelInd(Sfe5),Sfe6) = 
     &   MNeu2(1)*Sub113(HelInd(Sfe5),Sfe6) - 
     &    MNeu(1)*Sub119(HelInd(Sfe5),Sfe6)

        Sub129(HelInd(Sfe5),Sfe6) = 
     &   MNeu2(1)*Sub110(HelInd(Sfe5),Sfe6) - 
     &    MNeu(1)*Sub128(HelInd(Sfe5),Sfe6)

        Sub235(HelInd(Sfe5),Sfe6) = 
     &   pave4(PVC(dd2),Sfe5,Sfe6)*Sub204(HelInd(Sfe5),Sfe6) + 
     &    pave4(PVC(dd3),Sfe5,Sfe6)*Sub206(HelInd(Sfe5),Sfe6)

        Sub215(HelInd(Sfe5),Sfe6) = 
     &   Sub113(HelInd(Sfe5),Sfe6)*Sub213(Sfe6) - 
     &    Sub120(HelInd(Sfe5),Sfe6)*Sub214(Sfe6)

        Sub218(HelInd(Sfe5),Sfe6) = 
     &   Sub110(HelInd(Sfe5),Sfe6)*Sub216(Sfe6) - 
     &    Sub129(HelInd(Sfe5),Sfe6)*Sub217(Sfe6)

        Sub236(HelInd(Sfe5),Sfe6) = 
     &   pave4(PVC(dd1),Sfe5,Sfe6)*Sub198(HelInd(Sfe5),Sfe6) + 
     &    MNeu(1)*Sub235(HelInd(Sfe5),Sfe6)

        Sub207(HelInd(Sfe5),Sfe6) = 
     &   Sub130(HelInd(Sfe5),Sfe6) + Sub197(HelInd(Sfe5),Sfe6)

        Sub237(HelInd(Sfe5),Sfe6) = 
     &   pave5(PVC(dd1),Sfe5,Sfe6)*Sub199(HelInd(Sfe5),Sfe6)*
     &     Sub208(HelInd(Sfe5),Sfe6) + 
     &    Sub207(HelInd(Sfe5),Sfe6)*Sub236(HelInd(Sfe5),Sfe6)

        Sub219(HelInd(Sfe5),Sfe6) = 
     &   Sub215(HelInd(Sfe5),Sfe6)/(U - MSf2(Sfe5,4,1)) - 
     &    Sub218(HelInd(Sfe5),Sfe6)/(T - MSf2(Sfe5,4,1))

	ENDLOOP(Sfe5)
	ENDLOOP(Sfe6)

#include "contains.h"

	end

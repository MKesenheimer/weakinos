#if 0
* ddbar_n1n3_abbr1h.F
* abbreviations for {F[4, {1, Col1}], -F[4, {1, Col2}]} -> {F[11, {1}], F[11, {3}]}
* generated by FormCalc 8.4 on 7-Jul-2015 9:13
#endif

#include "ddbar_n1n3_vars.h"
#include "contains.h"


	subroutine ddbar_n1n3_abbr1h
	implicit none

#include "ddbar_n1n3_vars.h"


	F41 = ChainV1(Spinor(3,1,2),0,k(1),0,Spinor(4,-1,1))

	F36 = -ChainB1(Spinor(3,1,2),0,k(1),0,Spinor(4,-1,1))

	F23 = -ChainV0(Spinor(3,1,2),0,0,Spinor(4,-1,1))

	F19 = ChainB0(Spinor(3,1,2),0,0,Spinor(4,-1,1))

	F34 = ChainV1(Spinor(2,-1,2),0,k(1),0,Spinor(3,-1,1))

	F32 = ChainV1(Spinor(2,-1,2),0,k(1),1,Spinor(3,1,2))

	F40 = ChainV1(Spinor(2,-1,2),0,k(1),0,Spinor(4,-1,1))

	F38 = ChainV1(Spinor(2,-1,2),0,k(1),1,Spinor(4,1,2))

	F29 = ChainB1(Spinor(2,-1,2),0,k(1),0,Spinor(3,-1,1))

	F33 = -ChainB1(Spinor(2,-1,2),0,k(1),1,Spinor(3,1,2))

	F42 = ChainB1(Spinor(2,-1,2),0,k(1),0,Spinor(4,-1,1))

	F39 = -ChainB1(Spinor(2,-1,2),0,k(1),1,Spinor(4,1,2))

	F30 = ChainV1(Spinor(3,1,2),0,k(2),0,Spinor(1,1,1))

	F31 = ChainB1(Spinor(3,1,2),0,k(2),0,Spinor(1,1,1))

	F37 = ChainV1(Spinor(4,1,2),0,k(2),0,Spinor(1,1,1))

	F35 = ChainB1(Spinor(4,1,2),0,k(2),0,Spinor(1,1,1))

	F16 = ChainV0(Spinor(2,-1,2),0,0,Spinor(1,1,1))

	F21 = ChainV0(Spinor(2,-1,2),0,1,Spinor(3,1,2))

	F13 = -ChainB0(Spinor(2,-1,2),0,0,Spinor(1,1,1))

	F24 = -ChainB0(Spinor(2,-1,2),0,1,Spinor(3,1,2))

	F27 = -ChainV0(Spinor(3,-1,1),1,0,Spinor(1,1,1))

	F20 = ChainB0(Spinor(3,-1,1),1,0,Spinor(1,1,1))

	F25 = -ChainV0(Spinor(4,-1,1),1,0,Spinor(1,1,1))

	F22 = ChainB0(Spinor(4,-1,1),1,0,Spinor(1,1,1))

	F14 = -ChainV0(Spinor(4,1,2),0,0,Spinor(3,-1,1))

	F15 = ChainB0(Spinor(4,1,2),0,0,Spinor(3,-1,1))

	F17 = ChainV1(Spinor(4,1,2),0,k(1),0,Spinor(3,-1,1))

	F18 = -ChainB1(Spinor(4,1,2),0,k(1),0,Spinor(3,-1,1))

	F26 = ChainV1(Spinor(2,-1,2),0,k(3),0,Spinor(1,1,1))

	F28 = -ChainB1(Spinor(2,-1,2),0,k(3),0,Spinor(1,1,1))

	Sub193 = (F36 - F23*MNeu(1))*MNeu(3) + F19*MNeu2(3)

	Sub190 = (F41 - F19*MNeu(1))*MNeu(3) + F23*MNeu2(3)

	Sub210 = F35*MNeu(3) - F3*(Sub209 + MNeu2(3))

	Sub212 = F37*MNeu(3) - F8*(Sub209 + MNeu2(3))

	Sub244 = F18 - F15*MNeu(1) - F14*MNeu(3)

	Sub170 = F18 + F15*MNeu(1) - F14*MNeu(3)

	Sub243 = F17 - F14*MNeu(1) - F15*MNeu(3)

	Sub169 = F17 + F14*MNeu(1) - F15*MNeu(3)

	Sub204 = F31*F9*MNeu(1) + F10*(2*F9*Sub197 - F40*MNeu(3))

	Sub205 = F30*F9*MNeu(1) + F12*(2*F9*Sub197 - F40*MNeu(3))

	Sub206 = F11*F31*MNeu(1) + F10*(2*F11*Sub197 - F42*MNeu(3))

	Sub207 = F11*F30*MNeu(1) + F12*(2*F11*Sub197 - F42*MNeu(3))

	Sub173 = 2*F18 + F15*MNeu(1) - F14*MNeu(3)

	Sub172 = 2*F17 + F14*MNeu(1) - F15*MNeu(3)

        Sub166 = 1/SW2*
     &    ((F4*F5*Sub13 - F7*F8*Sub14)*Sub165 + 
     &      SW2*((2*F2*F3*Sub13 - 2*F1*F6*Sub14)*Sub164))

	Sub199 = F29*MNeu(1) - F2*(T - MNeu2(1))

	Sub198 = F34*MNeu(1) - F7*(T - MNeu2(1))

	Sub189 = F31*MNeu(1) + F10*(U - MNeu2(1))

	Sub191 = F30*MNeu(1) + F12*(U - MNeu2(1))

	Sub183 = F35*MNeu(3) + F3*(T - MNeu2(3))

	Sub185 = F37*MNeu(3) + F8*(T - MNeu2(3))

	Sub224 = F42*MNeu(3) - F11*(U - MNeu2(3))

	Sub222 = F40*MNeu(3) - F9*(U - MNeu2(3))

	Sub227 = F31*MNeu(1) + F10*(U - MNeu2(1) - MNeu2(3))

	Sub229 = F30*MNeu(1) + F12*(U - MNeu2(1) - MNeu2(3))

	Sub217 = F29*MNeu(1) - F2*(T - MNeu2(1) - MNeu2(3))

	Sub215 = F34*MNeu(1) - F7*(T - MNeu2(1) - MNeu2(3))

        Sub182 = -((F35*F7 + F15*F16*MNeu(1))*MNeu(3)) + 
     &    F16*(F17*MNeu(1) + F14*MNeu2(1)) - F3*F7*(T - MNeu2(3))

        Sub186 = -((F2*F37 + F13*F14*MNeu(1))*MNeu(3)) + 
     &    F13*(F18*MNeu(1) + F15*MNeu2(1)) - F2*F8*(T - MNeu2(3))

	Sub201 = F20*F39*MNeu(3) + F2*F8*(T - MNeu2(1))

	Sub220 = F8*Sub217 - (F2*F37 + F13*F14*MNeu(1))*MNeu(3)

	Sub216 = F3*Sub215 - (F35*F7 + F15*F16*MNeu(1))*MNeu(3)

	Sub228 = F9*Sub227 - (F10*F40 - F16*F19*MNeu(1))*MNeu(3)

	Sub230 = F9*Sub229 - (F12*F40 - F21*F22*MNeu(1))*MNeu(3)

	Sub232 = F11*Sub229 - (F12*F42 - F13*F23*MNeu(1))*MNeu(3)

	Sub231 = F11*Sub227 - (F10*F42 - F24*F25*MNeu(1))*MNeu(3)

	Sub219 = F8*Sub215 - (F37*F7 + F20*F5*MNeu(1))*MNeu(3)

	Sub218 = F3*Sub217 - (F2*F35 + F27*F6*MNeu(1))*MNeu(3)

	LOOP(Sfe5, 1,2,1)

        Sub211(HelInd(Sfe5)) = 
     &   (F7*Sub210 - F3*F34*MNeu(1))*Sub2(Sfe5) + 
     &    (F2*Sub210 - F29*F3*MNeu(1))*Sub3(Sfe5)

        Sub213(HelInd(Sfe5)) = 
     &   (F7*Sub212 - F34*F8*MNeu(1))*Sub2(Sfe5) + 
     &    (F2*Sub212 - F29*F8*MNeu(1))*Sub3(Sfe5)

        Sub124(HelInd(Sfe5)) = 
     &   F7*Sub123(Sfe5) + 
     &    F2*(MW2*Sub108(Sfe5) + CW2*(CB*Sub113(Sfe5)))*Sub3(Sfe5)

        Sub134(HelInd(Sfe5)) = 
     &   -(Sub133(Sfe5)*(F7*Sub2(Sfe5) + F2*Sub3(Sfe5))) + 
     &    F7*Sub122(Sfe5)*Sub6(Sfe5)

        Sub75(HelInd(Sfe5)) = 
     &   F10*Sub74(Sfe5) + F12*Sub73(Sfe5)*Sub9(Sfe5)

        Sub105(HelInd(Sfe5)) = 
     &   F10*(MW2*Sub49(Sfe5) + CW2*(CB*Sub56(Sfe5)))*
     &     Sub12(Sfe5) + Sub104(Sfe5)*
     &     (F10*Sub7(Sfe5) + F12*Sub9(Sfe5))

	Sub177(HelInd(Sfe5)) = 2*F27*F5*Sub2(Sfe5) + F2*F3*Sub3(Sfe5)

	Sub178(HelInd(Sfe5)) = F7*F8*Sub2(Sfe5) + 2*F20*F6*Sub3(Sfe5)

        Sub214(HelInd(Sfe5)) = 
     &   Sub211(HelInd(Sfe5))*Sub4(Sfe5) + 
     &    Sub213(HelInd(Sfe5))*Sub6(Sfe5)

        Sub181(HelInd(Sfe5)) = 
     &   F27*F5*Sub2(Sfe5)*Sub4(Sfe5) + 
     &    F20*F6*Sub3(Sfe5)*Sub6(Sfe5)

        Sub175(HelInd(Sfe5)) = 
     &   F10*F11*Sub7(Sfe5) + 2*F22*F24*Sub9(Sfe5)

        Sub174(HelInd(Sfe5)) = 
     &   2*F21*F25*Sub7(Sfe5) + F12*F9*Sub9(Sfe5)

        Sub180(HelInd(Sfe5)) = 
     &   F21*F25*Sub11(Sfe5)*Sub7(Sfe5) + 
     &    F22*F24*Sub12(Sfe5)*Sub9(Sfe5)

        Sub208(HelInd(Sfe5)) = 
     &   Sub11(Sfe5)*(Sub204*Sub7(Sfe5) + Sub205*Sub9(Sfe5)) + 
     &    Sub12(Sfe5)*(Sub206*Sub7(Sfe5) + Sub207*Sub9(Sfe5))

        Sub135(HelInd(Sfe5)) = 
     &   F3*Sub124(HelInd(Sfe5)) + F8*Sub134(HelInd(Sfe5)) + 
     &    F2*Sub128(Sfe5)*(F3*Sub4(Sfe5) + F8*Sub6(Sfe5))

        Sub106(HelInd(Sfe5)) = 
     &   F11*Sub105(HelInd(Sfe5)) + F9*Sub75(HelInd(Sfe5)) - 
     &    F12*(F9*Sub11(Sfe5) + F11*Sub12(Sfe5))*Sub93(Sfe5)

        Sub176(HelInd(Sfe5)) = 
     &   Sub11(Sfe5)*Sub174(HelInd(Sfe5)) + 
     &    Sub12(Sfe5)*Sub175(HelInd(Sfe5))

        Sub289(HelInd(Sfe5)) = 
     &   pave6(PVC(dd12),Sfe5)*Sub180(HelInd(Sfe5)) - 
     &    pave7(PVC(dd12),Sfe5)*Sub181(HelInd(Sfe5))

        Sub179(HelInd(Sfe5)) = 
     &   Sub177(HelInd(Sfe5))*Sub4(Sfe5) + 
     &    Sub178(HelInd(Sfe5))*Sub6(Sfe5)

        Sub277(HelInd(Sfe5)) = 
     &   3*Pi*Sub135(HelInd(Sfe5))/(CB*SW) + 
     &    2*CW2*MW2*(Alfas*(Sub274(Sfe5) + 2*Sub276(Sfe5))*
     &       (F7*Sub2(Sfe5) + F2*Sub3(Sfe5))*
     &       (F3*Sub4(Sfe5) + F8*Sub6(Sfe5)))

        Sub281(HelInd(Sfe5)) = 
     &   3*Pi*Sub106(HelInd(Sfe5))/(CB*SW) + 
     &    2*CW2*MW2*(Alfas*(F9*Sub11(Sfe5) + F11*Sub12(Sfe5))*
     &       (Sub278(Sfe5) + 2*Sub280(Sfe5))*
     &       (F10*Sub7(Sfe5) + F12*Sub9(Sfe5)))

        Sub194(HelInd(Sfe5)) = 
     &   Sub189*Sub7(Sfe5) + (F12*U + F30*MNeu(1))*Sub9(Sfe5)

        Sub184(HelInd(Sfe5)) = 
     &   -(Sub182*Sub2(Sfe5)) + F2*Sub183*Sub3(Sfe5)

        Sub187(HelInd(Sfe5)) = 
     &   F7*Sub185*Sub2(Sfe5) - Sub186*Sub3(Sfe5)

        Sub200(HelInd(Sfe5)) = 
     &   (F3*Sub198 - F27*F38*MNeu(3))*Sub2(Sfe5) + 
     &    F3*Sub199*Sub3(Sfe5)

        Sub202(HelInd(Sfe5)) = 
     &   Sub198*Sub2(Sfe5) + F29*MNeu(1)*Sub3(Sfe5)

        Sub192(HelInd(Sfe5)) = 
     &   -(F16*Sub190*Sub7(Sfe5)) + 
     &    F9*(Sub189*Sub7(Sfe5) + Sub191*Sub9(Sfe5))

        Sub223(HelInd(Sfe5)) = 
     &   (F10*Sub222 - F25*F32*MNeu(1))*Sub7(Sfe5) + 
     &    F12*Sub222*Sub9(Sfe5)

        Sub225(HelInd(Sfe5)) = 
     &   F10*Sub224*Sub7(Sfe5) + 
     &    (F12*Sub224 - F22*F33*MNeu(1))*Sub9(Sfe5)

        Sub195(HelInd(Sfe5)) = 
     &   F11*Sub194(HelInd(Sfe5)) - 
     &    (F13*Sub193 + F11*F12*MNeu2(1))*Sub9(Sfe5)

        Sub196(HelInd(Sfe5)) = 
     &   Sub11(Sfe5)*Sub192(HelInd(Sfe5)) + 
     &    Sub12(Sfe5)*Sub195(HelInd(Sfe5))

        Sub226(HelInd(Sfe5)) = 
     &   Sub11(Sfe5)*Sub223(HelInd(Sfe5)) + 
     &    Sub12(Sfe5)*Sub225(HelInd(Sfe5))

        Sub188(HelInd(Sfe5)) = 
     &   Sub184(HelInd(Sfe5))*Sub4(Sfe5) + 
     &    Sub187(HelInd(Sfe5))*Sub6(Sfe5)

        Sub203(HelInd(Sfe5)) = 
     &   Sub200(HelInd(Sfe5))*Sub4(Sfe5) + 
     &    (F8*Sub202(HelInd(Sfe5)) - Sub201*Sub3(Sfe5))*Sub6(Sfe5)

        Sub221(HelInd(Sfe5)) = 
     &   (Sub216*Sub2(Sfe5) + Sub218*Sub3(Sfe5))*Sub4(Sfe5) + 
     &    (Sub219*Sub2(Sfe5) + Sub220*Sub3(Sfe5))*Sub6(Sfe5)

        Sub233(HelInd(Sfe5)) = 
     &   Sub11(Sfe5)*(Sub228*Sub7(Sfe5) + Sub230*Sub9(Sfe5)) + 
     &    Sub12(Sfe5)*(Sub231*Sub7(Sfe5) + Sub232*Sub9(Sfe5))

        Sub290(HelInd(Sfe5)) = 
     &   pave7(PVC(dd13),Sfe5)*Sub188(HelInd(Sfe5)) + 
     &    pave7(PVC(dd3),Sfe5)*Sub214(HelInd(Sfe5)) + 
     &    pave6(PVC(dd23),Sfe5)*Sub226(HelInd(Sfe5))

        Sub291(HelInd(Sfe5)) = 
     &   pave6(PVC(dd13),Sfe5)*Sub196(HelInd(Sfe5)) + 
     &    pave7(PVC(dd23),Sfe5)*Sub203(HelInd(Sfe5)) + 
     &    pave6(PVC(dd3),Sfe5)*Sub208(HelInd(Sfe5)) + 
     &    pave7(PVC(dd33),Sfe5)*Sub221(HelInd(Sfe5)) + 
     &    pave6(PVC(dd33),Sfe5)*Sub233(HelInd(Sfe5))

	ENDLOOP(Sfe5)

	LOOP(Sfe6, 1,2,1)
	LOOP(Sfe5, 1,2,1)

        Sub167(HelInd(Sfe5),Sfe6) = 
     &   F8*Mino3C*Sub21(Sfe6)*USf(Sfe5,1,4,1)*USf(Sfe6,1,4,1) + 
     &    F3*Mino3*Sub22(Sfe6)*USf(Sfe5,2,4,1)*USf(Sfe6,2,4,1)

        Sub160(HelInd(Sfe5),Sfe6) = 
     &   F12*Mino3C*Sub30(Sfe6)*USf(Sfe5,1,4,1)*USf(Sfe6,1,4,1) + 
     &    F10*Mino3*Sub31(Sfe6)*USf(Sfe5,2,4,1)*USf(Sfe6,2,4,1)

        Sub158(HelInd(Sfe5),Sfe6) = 
     &   F9*Mino3*Sub138(Sfe6)*USfC(Sfe5,1,4,1)*USfC(Sfe6,1,4,1) + 
     &    F11*Mino3C*Sub137(Sfe6)*USfC(Sfe5,2,4,1)*USfC(Sfe6,2,4,1)

        Sub162(HelInd(Sfe5),Sfe6) = 
     &   F7*Mino3*Sub147(Sfe6)*USfC(Sfe5,1,4,1)*USfC(Sfe6,1,4,1) + 
     &    F2*Mino3C*Sub145(Sfe6)*USfC(Sfe5,2,4,1)*USfC(Sfe6,2,4,1)

        Sub157(HelInd(Sfe5),Sfe6) = 
     &   F8*Sub22(Sfe6)*USf(Sfe5,2,4,1)*USf(Sfe6,1,4,1) + 
     &    F3*Sub21(Sfe6)*USf(Sfe5,1,4,1)*USf(Sfe6,2,4,1)

        Sub136(HelInd(Sfe5),Sfe6) = 
     &   F12*Sub31(Sfe6)*USf(Sfe5,2,4,1)*USf(Sfe6,1,4,1) + 
     &    F10*Sub30(Sfe6)*USf(Sfe5,1,4,1)*USf(Sfe6,2,4,1)

        Sub139(HelInd(Sfe5),Sfe6) = 
     &   F9*Sub137(Sfe6)*USfC(Sfe5,2,4,1)*USfC(Sfe6,1,4,1) + 
     &    F11*Sub138(Sfe6)*USfC(Sfe5,1,4,1)*USfC(Sfe6,2,4,1)

        Sub150(HelInd(Sfe5),Sfe6) = 
     &   F7*Sub145(Sfe6)*USfC(Sfe5,2,4,1)*USfC(Sfe6,1,4,1) + 
     &    F2*Sub147(Sfe6)*USfC(Sfe5,1,4,1)*USfC(Sfe6,2,4,1)

        Sub142(HelInd(Sfe5),Sfe6) = 
     &   F10*Sub140(Sfe5,Sfe6) + F12*Sub141(Sfe5,Sfe6)

        Sub149(HelInd(Sfe5),Sfe6) = 
     &   F7*Sub146(Sfe5,Sfe6) + F2*Sub148(Sfe5,Sfe6)

        Sub151(HelInd(Sfe5),Sfe6) = 
     &   Sub150(HelInd(Sfe5),Sfe6)*Sub6(Sfe5) + 
     &    Sub22(Sfe6)*(F7*Sub2(Sfe5) + F2*Sub3(Sfe5))*
     &     USf(Sfe5,2,4,1)*USf(Sfe6,1,4,1)

        Sub143(HelInd(Sfe5),Sfe6) = 
     &   Sub12(Sfe5)*Sub136(HelInd(Sfe5),Sfe6) + 
     &    Sub138(Sfe6)*(F10*Sub7(Sfe5) + F12*Sub9(Sfe5))*
     &     USfC(Sfe5,1,4,1)*USfC(Sfe6,2,4,1)

        Sub257(HelInd(Sfe5),Sfe6) = 
     &   (F17*Sub13 - F18*Sub14)*pave20(PVC(cc1),Sfe5,Sfe6) + 
     &    (Sub13*Sub169 - Sub14*Sub170)*pave20(PVC(cc2),Sfe5,Sfe6)

        Sub144(HelInd(Sfe5),Sfe6) = 
     &   F9*Sub142(HelInd(Sfe5),Sfe6) + 
     &    F11*Sub143(HelInd(Sfe5),Sfe6)

        Sub250(HelInd(Sfe5),Sfe6) = 
     &   F15*Sub12(Sfe5)*Sub145(Sfe6) + 
     &    F14*Sub11(Sfe5)*Sub147(Sfe6)

        Sub248(HelInd(Sfe5),Sfe6) = 
     &   F14*Sub12(Sfe5)*Sub145(Sfe6) + 
     &    F15*Sub11(Sfe5)*Sub147(Sfe6)

        Sub245(HelInd(Sfe5),Sfe6) = 
     &   Sub244*Sub12(Sfe5)*Sub145(Sfe6) + 
     &    Sub243*Sub11(Sfe5)*Sub147(Sfe6)

        Sub240(HelInd(Sfe5),Sfe6) = 
     &   (F18 - 2*F14*MNeu(3))*Sub12(Sfe5)*Sub145(Sfe6) + 
     &    (F17 - 2*F15*MNeu(3))*Sub11(Sfe5)*Sub147(Sfe6)

        Sub246(HelInd(Sfe5),Sfe6) = 
     &   (F15*MNeu(1) + F14*MNeu(3))*Sub12(Sfe5)*Sub145(Sfe6) + 
     &    (F14*MNeu(1) + F15*MNeu(3))*Sub11(Sfe5)*Sub147(Sfe6)

        Sub152(HelInd(Sfe5),Sfe6) = 
     &   F3*Sub149(HelInd(Sfe5),Sfe6) + 
     &    F8*Sub151(HelInd(Sfe5),Sfe6)

        Sub239(HelInd(Sfe5),Sfe6) = 
     &   (F18 + 2*F15*MNeu(1))*Sub138(Sfe6)*Sub2(Sfe5) + 
     &    (F17 + 2*F14*MNeu(1))*Sub137(Sfe6)*Sub3(Sfe5)

        Sub249(HelInd(Sfe5),Sfe6) = 
     &   (F15*MNeu(1) + F14*MNeu(3))*Sub138(Sfe6)*Sub2(Sfe5) + 
     &    (F14*MNeu(1) + F15*MNeu(3))*Sub137(Sfe6)*Sub3(Sfe5)

        Sub241(HelInd(Sfe5),Sfe6) = 
     &   (F18 + F15*MNeu(1) + F14*MNeu(3))*Sub138(Sfe6)*
     &     Sub2(Sfe5) + (F17 + F14*MNeu(1) + F15*MNeu(3))*
     &     Sub137(Sfe6)*Sub3(Sfe5)

        Sub237(HelInd(Sfe5),Sfe6) = 
     &   (F18 - F14*MNeu(3))*Sub12(Sfe5)*Sub145(Sfe6) + 
     &    (F17 - F15*MNeu(3))*Sub11(Sfe5)*Sub147(Sfe6)

        Sub286(HelInd(Sfe5),Sfe6) = 
     &   MNeu(3)*pave5(PVC(dd2),Sfe5,Sfe6)*
     &     Sub248(HelInd(Sfe5),Sfe6) + 
     &    MNeu(1)*pave5(PVC(dd3),Sfe5,Sfe6)*
     &     Sub250(HelInd(Sfe5),Sfe6)

        Sub284(HelInd(Sfe5),Sfe6) = 
     &   pave5(PVC(dd23),Sfe5,Sfe6)*Sub246(HelInd(Sfe5),Sfe6) + 
     &    MNeu(1)*pave5(PVC(dd33),Sfe5,Sfe6)*
     &     Sub250(HelInd(Sfe5),Sfe6)

        Sub258(HelInd(Sfe5),Sfe6) = 
     &   (Sub13*Sub172 - Sub14*Sub173)*
     &     pave20(PVC(cc0),Sfe5,Sfe6) + 2*Sub257(HelInd(Sfe5),Sfe6)

        Sub242(HelInd(Sfe5),Sfe6) = 
     &   F15*Sub138(Sfe6)*Sub2(Sfe5) + F14*Sub137(Sfe6)*Sub3(Sfe5)

        Sub247(HelInd(Sfe5),Sfe6) = 
     &   F14*Sub138(Sfe6)*Sub2(Sfe5) + F15*Sub137(Sfe6)*Sub3(Sfe5)

        Sub238(HelInd(Sfe5),Sfe6) = 
     &   (F18 + F15*MNeu(1))*Sub138(Sfe6)*Sub2(Sfe5) + 
     &    (F17 + F14*MNeu(1))*Sub137(Sfe6)*Sub3(Sfe5)

        Sub171(HelInd(Sfe5),Sfe6) = 
     &   F16*Mino3*USf(Sfe5,2,4,1)*USfC(Sfe6,1,4,1) + 
     &    F13*Mino3C*USf(Sfe5,1,4,1)*USfC(Sfe6,2,4,1)

        Sub154(HelInd(Sfe5),Sfe6) = 
     &   F4*F5*USf(Sfe5,1,4,1)*USfC(Sfe6,1,4,1) + 
     &    F2*F3*USf(Sfe5,2,4,1)*USfC(Sfe6,2,4,1)

        Sub155(HelInd(Sfe5),Sfe6) = 
     &   F7*F8*USf(Sfe5,1,4,1)*USfC(Sfe6,1,4,1) + 
     &    F1*F6*USf(Sfe5,2,4,1)*USfC(Sfe6,2,4,1)

        Sub283(HelInd(Sfe5),Sfe6) = 
     &   pave5(PVC(dd11),Sfe5,Sfe6)*Sub237(HelInd(Sfe5),Sfe6) + 
     &    pave4(PVC(dd11),Sfe5,Sfe6)*Sub238(HelInd(Sfe5),Sfe6) + 
     &    pave4(PVC(dd12),Sfe5,Sfe6)*Sub239(HelInd(Sfe5),Sfe6) + 
     &    pave5(PVC(dd12),Sfe5,Sfe6)*Sub240(HelInd(Sfe5),Sfe6) + 
     &    pave4(PVC(dd13),Sfe5,Sfe6)*Sub241(HelInd(Sfe5),Sfe6) + 
     &    MNeu(1)*pave4(PVC(dd22),Sfe5,Sfe6)*
     &     Sub242(HelInd(Sfe5),Sfe6) + 
     &    pave5(PVC(dd13),Sfe5,Sfe6)*Sub245(HelInd(Sfe5),Sfe6) + 
     &    pave4(PVC(dd23),Sfe5,Sfe6)*Sub249(HelInd(Sfe5),Sfe6)

        Sub235(HelInd(Sfe5),Sfe6) = 
     &   Sub11(Sfe5)*Sub147(Sfe6)*Sub154(HelInd(Sfe5),Sfe6) + 
     &    Sub12(Sfe5)*Sub145(Sfe6)*Sub155(HelInd(Sfe5),Sfe6)

        Sub287(HelInd(Sfe5),Sfe6) = 
     &   pave5(PVC(dd1),Sfe5,Sfe6)*Sub237(HelInd(Sfe5),Sfe6) - 
     &    Sub286(HelInd(Sfe5),Sfe6)

        Sub234(HelInd(Sfe5),Sfe6) = 
     &   Sub138(Sfe6)*Sub155(HelInd(Sfe5),Sfe6)*Sub2(Sfe5) + 
     &    Sub137(Sfe6)*Sub154(HelInd(Sfe5),Sfe6)*Sub3(Sfe5)

        Sub236(HelInd(Sfe5),Sfe6) = 
     &   F26*USf(Sfe5,1,4,1)*USfC(Sfe6,1,4,1) + 
     &    F28*USf(Sfe5,2,4,1)*USfC(Sfe6,2,4,1)

        Sub288(HelInd(Sfe5),Sfe6) = 
     &   pave4(PVC(dd1),Sfe5,Sfe6)*Sub238(HelInd(Sfe5),Sfe6) + 
     &    MNeu(1)*pave4(PVC(dd2),Sfe5,Sfe6)*
     &     Sub242(HelInd(Sfe5),Sfe6) + 
     &    MNeu(3)*pave4(PVC(dd3),Sfe5,Sfe6)*
     &     Sub247(HelInd(Sfe5),Sfe6)

        Sub251(HelInd(Sfe5),Sfe6) = 
     &   -Sub171(HelInd(Sfe5),Sfe6) + Sub236(HelInd(Sfe5),Sfe6)

        Sub252(HelInd(Sfe5),Sfe6) = 
     &   Sub171(HelInd(Sfe5),Sfe6) + Sub236(HelInd(Sfe5),Sfe6)

        Sub159(HelInd(Sfe5),Sfe6) = 
     &   MNeu2(3)*Sub139(HelInd(Sfe5),Sfe6) - 
     &    MNeu(3)*Sub158(HelInd(Sfe5),Sfe6)

        Sub161(HelInd(Sfe5),Sfe6) = 
     &   MNeu2(1)*Sub136(HelInd(Sfe5),Sfe6) - 
     &    MNeu(1)*Sub160(HelInd(Sfe5),Sfe6)

        Sub163(HelInd(Sfe5),Sfe6) = 
     &   MNeu2(1)*Sub150(HelInd(Sfe5),Sfe6) - 
     &    MNeu(1)*Sub162(HelInd(Sfe5),Sfe6)

        Sub168(HelInd(Sfe5),Sfe6) = 
     &   MNeu2(3)*Sub157(HelInd(Sfe5),Sfe6) - 
     &    MNeu(3)*Sub167(HelInd(Sfe5),Sfe6)

        Sub261(HelInd(Sfe5),Sfe6) = 
     &   (MNeu2(3) - MSf2(Sfe6,4,1))*pave19(PVC(cc0),Sfe6)*
     &     Sub157(HelInd(Sfe5),Sfe6) + 
     &    Sub168(HelInd(Sfe5),Sfe6)*Sub260(Sfe6)

        Sub263(HelInd(Sfe5),Sfe6) = 
     &   (MNeu2(1) - MSf2(Sfe6,4,1))*pave18(PVC(cc0),Sfe6)*
     &     Sub150(HelInd(Sfe5),Sfe6) + 
     &    Sub163(HelInd(Sfe5),Sfe6)*Sub262(Sfe6)

        Sub267(HelInd(Sfe5),Sfe6) = 
     &   (MNeu2(1) - MSf2(Sfe6,4,1))*pave16(PVC(cc0),Sfe6)*
     &     Sub136(HelInd(Sfe5),Sfe6) + 
     &    Sub161(HelInd(Sfe5),Sfe6)*Sub266(Sfe6)

        Sub269(HelInd(Sfe5),Sfe6) = 
     &   (MNeu2(3) - MSf2(Sfe6,4,1))*pave17(PVC(cc0),Sfe6)*
     &     Sub139(HelInd(Sfe5),Sfe6) + 
     &    Sub159(HelInd(Sfe5),Sfe6)*Sub268(Sfe6)

        Sub264(HelInd(Sfe5),Sfe6) = 
     &   Sub261(HelInd(Sfe5),Sfe6)*
     &     (F7*Sub2(Sfe5) + F2*Sub3(Sfe5)) + 
     &    Sub263(HelInd(Sfe5),Sfe6)*(F3*Sub4(Sfe5) + F8*Sub6(Sfe5))

        Sub270(HelInd(Sfe5),Sfe6) = 
     &   (F9*Sub11(Sfe5) + F11*Sub12(Sfe5))*
     &     Sub267(HelInd(Sfe5),Sfe6) + 
     &    Sub269(HelInd(Sfe5),Sfe6)*
     &     (F10*Sub7(Sfe5) + F12*Sub9(Sfe5))

        Sub156(HelInd(Sfe5),Sfe6) = 
     &   Sub13*Sub154(HelInd(Sfe5),Sfe6) - 
     &    Sub14*Sub155(HelInd(Sfe5),Sfe6)

        Sub282(HelInd(Sfe5),Sfe6) = 
     &   -(pave4(PVC(dd33),Sfe5,Sfe6)*Sub247(HelInd(Sfe5),Sfe6)) + 
     &    pave5(PVC(dd22),Sfe5,Sfe6)*Sub248(HelInd(Sfe5),Sfe6)

        Sub265(HelInd(Sfe5),Sfe6) = 
     &   pave21(PVC(bb0))*Sub152(HelInd(Sfe5),Sfe6) - 
     &    Sub264(HelInd(Sfe5),Sfe6)

        Sub271(HelInd(Sfe5),Sfe6) = 
     &   pave22(PVC(bb0))*Sub144(HelInd(Sfe5),Sfe6) - 
     &    Sub270(HelInd(Sfe5),Sfe6)

        Sub259(HelInd(Sfe5),Sfe6) = 
     &   4*pave20(PVC(cc00),Sfe5,Sfe6)*Sub156(HelInd(Sfe5),Sfe6) + 
     &    Sub171(HelInd(Sfe5),Sfe6)*Sub258(HelInd(Sfe5),Sfe6)

        Sub285(HelInd(Sfe5),Sfe6) = 
     &   MNeu(3)*Sub282(HelInd(Sfe5),Sfe6) - 
     &    Sub283(HelInd(Sfe5),Sfe6) + Sub284(HelInd(Sfe5),Sfe6)

        Sub272(HelInd(Sfe5),Sfe6) = 
     &   Sub265(HelInd(Sfe5),Sfe6)/(U - MSf2(Sfe5,4,1)) - 
     &    Sub271(HelInd(Sfe5),Sfe6)/(T - MSf2(Sfe5,4,1))

	ENDLOOP(Sfe5)
	ENDLOOP(Sfe6)

#include "contains.h"

	end

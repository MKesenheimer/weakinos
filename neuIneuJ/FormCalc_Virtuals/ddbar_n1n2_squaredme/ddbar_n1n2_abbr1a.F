#if 0
* ddbar_n1n2_abbr1a.F
* abbreviations for {F[4, {1, Col1}], -F[4, {1, Col2}]} -> {F[11, {1}], F[11, {2}]}
* generated by FormCalc 8.4 on 10-Jun-2015 12:08
#endif

#include "ddbar_n1n2_vars.h"
#include "contains.h"


	subroutine ddbar_n1n2_abbr1a
	implicit none

#include "ddbar_n1n2_vars.h"


	Sub209 = -2*(T + U) + 2*MNeu2(1) + MNeu2(2)

	Sub189 = T + U - MNeu2(1) - MNeu2(2)

	call Bput(pave20(bb0),U,MD2,MGl2)

	call Bput(pave21(bb0),T,MD2,MGl2)

	LOOP(Sfe5, 1,2,1)

	call Cput(pave8(cc0,Sfe5),MD2,MNeu2(1),T,0.D0,MD2,MSf2(Sfe5,4,1))

	call Cput(pave9(cc0,Sfe5),MD2,MNeu2(2),T,0.D0,MD2,MSf2(Sfe5,4,1))

	call Cput(pave10(cc0,Sfe5),MD2,MNeu2(1),U,0.D0,MD2,MSf2(Sfe5,4,1))

	call Cput(pave11(cc0,Sfe5),MD2,MNeu2(2),U,0.D0,MD2,MSf2(Sfe5,4,1))

        Sub278(Sfe5) = 
     &   (MD2 + 2*T - 2*MNeu2(1))*pave8(PVC(cc0),Sfe5) + 
     &    (MD2 + 2*T - 2*MNeu2(2))*pave9(PVC(cc0),Sfe5)

        Sub274(Sfe5) = 
     &   (MD2 + 2*U - 2*MNeu2(1))*pave10(PVC(cc0),Sfe5) + 
     &    (MD2 + 2*U - 2*MNeu2(2))*pave11(PVC(cc0),Sfe5)

	call Bput(pave13(bb0,Sfe5),T,0.D0,MSf2(Sfe5,4,1))

	call Bput(pave14(bb0,Sfe5),U,0.D0,MSf2(Sfe5,4,1))

        Sub254(Sfe5) = 
     &   MSf2(Sfe5,4,1)*pave13(PVC(bb0),Sfe5) + 
     &    T*(3*pave13(PVC(bb0),Sfe5) + 2*pave13(PVC(bb1),Sfe5))

        Sub253(Sfe5) = 
     &   MSf2(Sfe5,4,1)*pave14(PVC(bb0),Sfe5) + 
     &    U*(3*pave14(PVC(bb0),Sfe5) + 2*pave14(PVC(bb1),Sfe5))

        call Dput(pave6(dd0,Sfe5),MD2,S,MNeu2(1),T,MD2,MNeu2(2),
     &   0.D0,MD2,MD2,MSf2(Sfe5,4,1))

        call Dput(pave7(dd0,Sfe5),MD2,S,MNeu2(2),U,MD2,MNeu2(1),
     &   0.D0,MD2,MD2,MSf2(Sfe5,4,1))

	Sub279(Sfe5) = pave8(PVC(cc2),Sfe5) + pave9(PVC(cc2),Sfe5)

	Sub275(Sfe5) = pave10(PVC(cc2),Sfe5) + pave11(PVC(cc2),Sfe5)

	Sub291(Sfe5) = pave6(PVC(dd0),Sfe5) + pave6(PVC(dd1),Sfe5)

	Sub287(Sfe5) = pave7(PVC(dd0),Sfe5) + pave7(PVC(dd1),Sfe5)

        Sub288(Sfe5) = 
     &   pave12(PVC(cc0),Sfe5) + Sub189*pave7(PVC(dd2),Sfe5)

	Sub292(Sfe5) = -(Sub189*pave6(PVC(dd2),Sfe5)) + S*Sub291(Sfe5)

        Sub276(Sfe5) = 
     &   (U - MNeu2(1))*pave10(PVC(cc1),Sfe5) + 
     &    (U - MNeu2(2))*pave11(PVC(cc1),Sfe5) + 
     &    pave14(PVC(bb0),Sfe5) + U*Sub275(Sfe5)

        Sub280(Sfe5) = 
     &   pave13(PVC(bb0),Sfe5) + 
     &    (T - MNeu2(1))*pave8(PVC(cc1),Sfe5) + 
     &    (T - MNeu2(2))*pave9(PVC(cc1),Sfe5) + T*Sub279(Sfe5)

	ENDLOOP(Sfe5)

	LOOP(Sfe6, 1,2,1)

        call Cput(pave15(cc0,Sfe6),T,MD2,MNeu2(1),MD2,MGl2,
     &   MSf2(Sfe6,4,1))

        call Cput(pave16(cc0,Sfe6),T,MD2,MNeu2(2),MD2,MGl2,
     &   MSf2(Sfe6,4,1))

        call Cput(pave17(cc0,Sfe6),U,MD2,MNeu2(1),MD2,MGl2,
     &   MSf2(Sfe6,4,1))

        call Cput(pave18(cc0,Sfe6),U,MD2,MNeu2(2),MD2,MGl2,
     &   MSf2(Sfe6,4,1))

	Sub268(Sfe6) = pave15(PVC(cc1),Sfe6) + pave15(PVC(cc2),Sfe6)

	Sub266(Sfe6) = pave16(PVC(cc1),Sfe6) + pave16(PVC(cc2),Sfe6)

	Sub262(Sfe6) = pave17(PVC(cc1),Sfe6) + pave17(PVC(cc2),Sfe6)

	Sub260(Sfe6) = pave18(PVC(cc1),Sfe6) + pave18(PVC(cc2),Sfe6)

	ENDLOOP(Sfe6)

	LOOP(Sfe6, 1,2,1)
	LOOP(Sfe5, 1,2,1)

        Sub255(Sfe5,Sfe6) = 
     &   3*Pi*(Sub16(Sfe5,Sfe6) - T*Sub17(Sfe5,Sfe6)) + 
     &    8*Alfas*(pave22 + T*pave21(PVC(bb1)))*Sub32(Sfe5,Sfe6)

        Sub256(Sfe5,Sfe6) = 
     &   3*Pi*(Sub16(Sfe5,Sfe6) - U*Sub17(Sfe5,Sfe6)) + 
     &    8*Alfas*(pave22 + U*pave20(PVC(bb1)))*Sub32(Sfe5,Sfe6)

        call Dput(pave4(dd0,Sfe5,Sfe6),T,MD2,S,MNeu2(2),MNeu2(1),
     &   MD2,MD2,MGl2,MSf2(Sfe5,4,1),MSf2(Sfe6,4,1))

        call Dput(pave5(dd0,Sfe5,Sfe6),U,MD2,S,MNeu2(1),MNeu2(2),
     &   MD2,MD2,MGl2,MSf2(Sfe5,4,1),MSf2(Sfe6,4,1))

	ENDLOOP(Sfe5)
	ENDLOOP(Sfe6)

#include "contains.h"

	end

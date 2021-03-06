#if 0
* uubar_n1n1_abbr1a.F
* abbreviations for {F[3, {1, Col1}], -F[3, {1, Col2}]} -> {F[11, {1}], F[11, {1}]}
* generated by FormCalc 8.4 on 5-Mar-2015 13:58
#endif

#include "uubar_n1n1_vars.h"
#include "contains.h"


	subroutine uubar_n1n1_abbr1a
	implicit none

#include "uubar_n1n1_vars.h"


	Sub172 = S - T - U + 2*MNeu2(1)

	call Bput(pave6(bb0),T,MGl2,MU2)

	call Bput(pave7(bb0),U,MGl2,MU2)

	Sub216 = pave6(PVC(bb0)) + pave6(PVC(bb1))

	Sub218 = pave7(PVC(bb0)) + pave7(PVC(bb1))

	LOOP(Sfe5, 1,2,1)

	call Cput(pave9(cc0,Sfe5),MU2,MNeu2(1),T,0.D0,MU2,MSf2(Sfe5,3,1))

	call Cput(pave11(cc0,Sfe5),MU2,MNeu2(1),U,0.D0,MU2,MSf2(Sfe5,3,1))

        Sub227(Sfe5) = 
     &   (T - MNeu2(1))*pave9(PVC(cc1),Sfe5) + 
     &    T*pave9(PVC(cc2),Sfe5)

        Sub230(Sfe5) = 
     &   (U - MNeu2(1))*pave11(PVC(cc1),Sfe5) + 
     &    U*pave11(PVC(cc2),Sfe5)

	call Bput(pave8(bb0,Sfe5),T,0.D0,MSf2(Sfe5,3,1))

	call Bput(pave10(bb0,Sfe5),U,0.D0,MSf2(Sfe5,3,1))

        call Dput(pave12(dd0,Sfe5),MU2,S,MNeu2(1),T,MU2,MNeu2(1),
     &   0.D0,MU2,MU2,MSf2(Sfe5,3,1))

        call Dput(pave13(dd0,Sfe5),MU2,S,MNeu2(1),U,MU2,MNeu2(1),
     &   0.D0,MU2,MU2,MSf2(Sfe5,3,1))

        Sub214(Sfe5) = 
     &   MSf2(Sfe5,3,1)*pave8(PVC(bb0),Sfe5) + 
     &    T*(3*pave8(PVC(bb0),Sfe5) + 2*pave8(PVC(bb1),Sfe5))

        Sub215(Sfe5) = 
     &   MSf2(Sfe5,3,1)*pave10(PVC(bb0),Sfe5) + 
     &    U*(3*pave10(PVC(bb0),Sfe5) + 2*pave10(PVC(bb1),Sfe5))

        Sub240(Sfe5) = 
     &   pave12(PVC(dd0),Sfe5) + pave12(PVC(dd1),Sfe5) + 
     &    pave12(PVC(dd2),Sfe5)

        Sub243(Sfe5) = 
     &   pave13(PVC(dd0),Sfe5) + pave13(PVC(dd1),Sfe5) + 
     &    pave13(PVC(dd2),Sfe5)

        Sub226(Sfe5) = 
     &   pave8(PVC(bb0),Sfe5) + 
     &    (MU2 + 2*T - 2*MNeu2(1))*pave9(PVC(cc0),Sfe5)

        Sub229(Sfe5) = 
     &   pave10(PVC(bb0),Sfe5) + 
     &    (MU2 + 2*U - 2*MNeu2(1))*pave11(PVC(cc0),Sfe5)

	ENDLOOP(Sfe5)

	LOOP(Sfe6, 1,2,1)

        call Cput(pave16(cc0,Sfe6),T,MNeu2(1),MU2,MGl2,MU2,
     &   MSf2(Sfe6,3,1))

        call Cput(pave17(cc0,Sfe6),U,MNeu2(1),MU2,MGl2,MU2,
     &   MSf2(Sfe6,3,1))

	ENDLOOP(Sfe6)

	LOOP(Sfe6, 1,2,1)
	LOOP(Sfe5, 1,2,1)

        call Dput(pave4(dd0,Sfe5,Sfe6),U,MNeu2(1),S,MU2,MU2,
     &   MNeu2(1),MGl2,MU2,MSf2(Sfe5,3,1),MSf2(Sfe6,3,1))

        call Dput(pave5(dd0,Sfe5,Sfe6),T,MNeu2(1),S,MU2,MU2,
     &   MNeu2(1),MGl2,MU2,MSf2(Sfe5,3,1),MSf2(Sfe6,3,1))

        Sub217(Sfe5,Sfe6) = 
     &   3*Pi*(Sub10(Sfe5,Sfe6) - T*Sub9(Sfe5,Sfe6)) + 
     &    8*Alfas*(pave15 - Sub216*T)*Sub23(Sfe5,Sfe6)

        Sub219(Sfe5,Sfe6) = 
     &   3*Pi*(Sub10(Sfe5,Sfe6) - U*Sub9(Sfe5,Sfe6)) + 
     &    8*Alfas*(pave15 - Sub218*U)*Sub23(Sfe5,Sfe6)

	ENDLOOP(Sfe5)
	ENDLOOP(Sfe6)

#include "contains.h"

	end

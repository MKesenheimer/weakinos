#if 0
* CalcRenConst.F
* renormalization constants
* generated by FormCalc 8.4 on 20-Jan-2015 10:09
#endif

#include "decl.h"
#include "contains.h"


        subroutine CalcRenConst
        implicit none

#include "decl.h"

        integer Sfe3
        double precision ltlambda
        common/debug/ltlambda

        ! ssbar dup's
        ComplexType dup1, dup10, dup11, dup12, dup13, dup14, dup15
        ComplexType dup16, dup17, dup18, dup19, dup2, dup20, dup21
        ComplexType dup22, dup23, dup24, dup25, dup26, dup27, dup28
        ComplexType dup3, dup4, dup5, dup6, dup7, dup8, dup9

        ComplexType CNull
        ! MK: FormCalc creates an additional Finite Term if the
        !     quark masses are zero.
        !     Setting Finite = 0 solves the problem for SM-RenConsts.
        !     This problem occurs only when there are UV- and IR-
        !     divergences at the same time.
        RealType xm0,xz0
        
        CNull = (0d0, 0d0)
        xz0 = 0d0 ! set Finite part in field renormalization constant equal to zero
        xm0 = 0d0 ! set Finite part in mass renormalization constant equal to zero
        
c######################## common renormalization constants #############
        
        dMChaOS1(1) = 0
        dMChaOS1(2) = 0
        dMNeuOS1(1) = 0
        dZH1        = 0
        dZH2        = 0
        dMZsq1      = 0
        dTB1        = 0
        dZAA1       = 0
        dZZA1       = 0
        dMWsq1      = 0
        dSB1        = 0
        dCB1        = 0
        dSW1        = 0
        dZAZ1       = 0
        dZe1        = 0
        dZZZ1       = 0
        dZNeu1(1)   = 0
        dMino11     = 0
        dMino21     = 0
        dMUE1       = 0
        dAf1(4,1,1) = 0
        dAf1(4,2,2) = 0
        dAf1(4,3,3) = 0
        dMNeu1(1,1) = 0
        dMNeu1(1,2) = 0
        dMNeu1(1,3) = 0
        dMNeu1(1,4) = 0
        dMNeu1(2,1) = 0
        dMNeu1(3,1) = 0
        dMNeu1(4,1) = 0
        dZfL1(11,1,1) = 0
        dZfL1(11,2,1) = 0
        dZfL1(11,3,1) = 0
        dZfL1(11,4,1) = 0
        dZfR1(11,1,1) = 0
        dZfR1(11,2,1) = 0
        dZfR1(11,3,1) = 0
        dZfR1(11,4,1) = 0
        dMsq12Sf1(4,1) = 0
        dZW1 = 0

#ifdef DEBUGQ
        print*,"Lambda = ", ltlambda
        print*,"xz0 = ", xz0
        print*,"xm0 = ", xm0
#endif

c######################## uubar renormalization constants ##############
        ! reset dup's for uubar CalcRenConst
        
        dup1  = CNull
        dup10 = CNull
        dup11 = CNull
        dup12 = CNull
        dup13 = CNull
        dup14 = CNull
        dup15 = CNull
        dup16 = CNull
        dup17 = CNull
        dup18 = CNull
        dup19 = CNull
        dup2  = CNull
        dup20 = CNull
        dup21 = CNull
        dup22 = CNull
        dup23 = CNull
        dup24 = CNull
        dup25 = CNull
        dup26 = CNull
        dup27 = CNull
        dup28 = CNull
        dup3  = CNull
        dup4  = CNull
        dup5  = CNull
        dup6  = CNull
        dup7  = CNull
        dup8  = CNull
        dup9  = CNull
        
        dMf1(3,1) = MU/Pi*
     &    (Alfas*(1/3.D0*Finite*xm0 - 2/3.D0*Re(B0i(bb0,MU2,0.D0,MU2))+ 
     &        2/3.D0*Re(B0i(bb1,MU2,0.D0,MU2))))

        LOOP(Sfe3, 1,2,1)

        dMf1(3,1) = dMf1(3,1) - 
     &    1/(3.D0*Pi)*(Alfas*
     &       (MU*(Re(B0i(bb1,MU2,MGl2,MSf2(Sfe3,3,1)))*
     &            (USf(Sfe3,1,3,1)*USfC(Sfe3,1,3,1) + 
     &              USf(Sfe3,2,3,1)*USfC(Sfe3,2,3,1))) + 
     &         Re(B0i(bb0,MU2,MGl2,MSf2(Sfe3,3,1)))*
     &          (Mino3*USf(Sfe3,2,3,1)*USfC(Sfe3,1,3,1) + 
     &            Mino3C*USf(Sfe3,1,3,1)*USfC(Sfe3,2,3,1))))

        ENDLOOP(Sfe3)

#ifdef DEBUGQ
        DEB("dMf1(3,1) =", dMf1(3,1))
#endif

        dup1 = MU*(Mino3*USf(1,2,3,1)) - 
     &    (MU2 - MSf2(1,3,1))*USf(1,1,3,1)

        dup2 = MU*(Mino3C*USf(1,1,3,1)) + 
     &    (-MU2 + MSf2(1,3,1))*USf(1,2,3,1)

        dMSfsq1(1,1,3,1) = 
     &   -(1/Pi*(Alfas*(-(4/3.D0*
     &             (Re(B0i(bb0,MSf2(1,3,1),MGl2,MU2))*
     &               (dup1*USfC(1,1,3,1) + dup2*USfC(1,2,3,1)))) + 
     &          4/3.D0*((Re(A0(MGl2)) - 
     &               MSf2(1,3,1)*Re(B0i(bb1,MSf2(1,3,1),MGl2,MU2)))
     &              *(USf(1,1,3,1)*USfC(1,1,3,1) + 
     &               USf(1,2,3,1)*USfC(1,2,3,1))) + 
     &          (4/3.D0*Re(B0i(bb0,MSf2(1,3,1),0.D0,MSf2(1,3,1))) + 
     &             2/3.D0*Re(B0i(bb1,MSf2(1,3,1),0.D0,MSf2(1,3,1))))*
     &           MSf2(1,3,1))))

        LOOP(Sfe3, 1,2,1)

        dMSfsq1(1,1,3,1) = 
     &   dMSfsq1(1,1,3,1) + 
     &    1/(3.D0*Pi)*(Alfas*Re(A0(MSf2(Sfe3,3,1)))*
     &       (USf(Sfe3,1,3,1)*USfC(1,1,3,1) - 
     &         USf(Sfe3,2,3,1)*USfC(1,2,3,1))*
     &       (USf(1,1,3,1)*USfC(Sfe3,1,3,1) - 
     &         USf(1,2,3,1)*USfC(Sfe3,2,3,1)))

        ENDLOOP(Sfe3)

#ifdef DEBUGQ
        DEB("dMSfsq1(1,1,3,1) =", dMSfsq1(1,1,3,1))
#endif

        dup3 = USf(1,1,3,1)*USfC(2,1,3,1) + USf(1,2,3,1)*USfC(2,2,3,1)

        dup4 = MU*(Mino3C*USf(1,1,3,1)) + 
     &    (-MU2 + MSf2(2,3,1))*USf(1,2,3,1)

        dMSfsq1(1,2,3,1) = 
     &   -(1/Pi*(Alfas*(dup3*
     &           (4/3.D0*Re(A0(MGl2)) - 
     &             2/3.D0*(MSf2(1,3,1)*
     &                 Re(B0i(bb1,MSf2(1,3,1),MGl2,MU2)) + 
     &                MSf2(2,3,1)*Re(B0i(bb1,MSf2(2,3,1),MGl2,MU2))
     &                )) - 
     &          2/3.D0*(Re(B0i(bb0,MSf2(1,3,1),MGl2,MU2))*
     &             (dup1*USfC(2,1,3,1) + dup2*USfC(2,2,3,1))) + 
     &          2/3.D0*(Re(B0i(bb0,MSf2(2,3,1),MGl2,MU2))*
     &             ((-(MU*(Mino3*USf(1,2,3,1))) + 
     &                  (MU2 - MSf2(2,3,1))*USf(1,1,3,1))*
     &                USfC(2,1,3,1) - dup4*USfC(2,2,3,1))))))

        LOOP(Sfe3, 1,2,1)

        dMSfsq1(1,2,3,1) = 
     &   dMSfsq1(1,2,3,1) + 
     &    1/(3.D0*Pi)*(Alfas*Re(A0(MSf2(Sfe3,3,1)))*
     &       (USf(Sfe3,1,3,1)*USfC(2,1,3,1) - 
     &         USf(Sfe3,2,3,1)*USfC(2,2,3,1))*
     &       (USf(1,1,3,1)*USfC(Sfe3,1,3,1) - 
     &         USf(1,2,3,1)*USfC(Sfe3,2,3,1)))

        ENDLOOP(Sfe3)

#ifdef DEBUGQ
        DEB("dMSfsq1(1,2,3,1) =", dMSfsq1(1,2,3,1))
#endif

        dup5 = USf(2,1,3,1)*USfC(1,1,3,1) + USf(2,2,3,1)*USfC(1,2,3,1)

        dup6 = MU*(Mino3*USf(2,2,3,1)) - 
     &    (MU2 - MSf2(2,3,1))*USf(2,1,3,1)

        dup7 = MU*(Mino3C*USf(2,1,3,1)) + 
     &    (-MU2 + MSf2(1,3,1))*USf(2,2,3,1)

        dup8 = MU*(Mino3C*USf(2,1,3,1)) + 
     &    (-MU2 + MSf2(2,3,1))*USf(2,2,3,1)

        dMSfsq1(2,1,3,1) = 
     &   -(1/Pi*(Alfas*(dup5*
     &           (4/3.D0*Re(A0(MGl2)) - 
     &             2/3.D0*(MSf2(1,3,1)*
     &                 Re(B0i(bb1,MSf2(1,3,1),MGl2,MU2)) + 
     &                MSf2(2,3,1)*Re(B0i(bb1,MSf2(2,3,1),MGl2,MU2))
     &                )) + 
     &          2/3.D0*(Re(B0i(bb0,MSf2(1,3,1),MGl2,MU2))*
     &             ((-(MU*(Mino3*USf(2,2,3,1))) + 
     &                  (MU2 - MSf2(1,3,1))*USf(2,1,3,1))*
     &                USfC(1,1,3,1) - dup7*USfC(1,2,3,1))) - 
     &          2/3.D0*(Re(B0i(bb0,MSf2(2,3,1),MGl2,MU2))*
     &             (dup6*USfC(1,1,3,1) + dup8*USfC(1,2,3,1))))))

        LOOP(Sfe3, 1,2,1)

        dMSfsq1(2,1,3,1) = 
     &   dMSfsq1(2,1,3,1) + 
     &    1/(3.D0*Pi)*(Alfas*Re(A0(MSf2(Sfe3,3,1)))*
     &       (USf(Sfe3,1,3,1)*USfC(1,1,3,1) - 
     &         USf(Sfe3,2,3,1)*USfC(1,2,3,1))*
     &       (USf(2,1,3,1)*USfC(Sfe3,1,3,1) - 
     &         USf(2,2,3,1)*USfC(Sfe3,2,3,1)))

        ENDLOOP(Sfe3)

#ifdef DEBUGQ
        DEB("dMSfsq1(2,1,3,1) =", dMSfsq1(2,1,3,1))
#endif

        dMSfsq1(2,2,3,1) = 
     &   -(1/Pi*(Alfas*(-(4/3.D0*
     &             (Re(B0i(bb0,MSf2(2,3,1),MGl2,MU2))*
     &               (dup6*USfC(2,1,3,1) + dup8*USfC(2,2,3,1)))) + 
     &          4/3.D0*((Re(A0(MGl2)) - 
     &               MSf2(2,3,1)*Re(B0i(bb1,MSf2(2,3,1),MGl2,MU2)))
     &              *(USf(2,1,3,1)*USfC(2,1,3,1) + 
     &               USf(2,2,3,1)*USfC(2,2,3,1))) + 
     &          (4/3.D0*Re(B0i(bb0,MSf2(2,3,1),0.D0,MSf2(2,3,1))) + 
     &             2/3.D0*Re(B0i(bb1,MSf2(2,3,1),0.D0,MSf2(2,3,1))))*
     &           MSf2(2,3,1))))

        LOOP(Sfe3, 1,2,1)

        dMSfsq1(2,2,3,1) = 
     &   dMSfsq1(2,2,3,1) + 
     &    1/(3.D0*Pi)*(Alfas*Re(A0(MSf2(Sfe3,3,1)))*
     &       (USf(Sfe3,1,3,1)*USfC(2,1,3,1) - 
     &         USf(Sfe3,2,3,1)*USfC(2,2,3,1))*
     &       (USf(2,1,3,1)*USfC(Sfe3,1,3,1) - 
     &         USf(2,2,3,1)*USfC(Sfe3,2,3,1)))

        ENDLOOP(Sfe3)

#ifdef DEBUGQ
        DEB("dMSfsq1(2,2,3,1) =", dMSfsq1(2,2,3,1))
#endif

        dup9 = 1/Pi*(Alfas*
     &      (1/3.D0*Finite*xz0 - 
     &        2/3.D0*(B0i(bb0,MU2,0.D0,MU2) + B0i(bb1,MU2,0.D0,MU2)) + 
     &        MU2*(4/3.D0*B0i(dbb0,MU2,0.D0,MU2) - 
     &           4/3.D0*B0i(dbb1,MU2,0.D0,MU2))))

        dZbarfL1(3,1,1) = dup9

        LOOP(Sfe3, 1,2,1)

        dZbarfL1(3,1,1) = 
     &   dZbarfL1(3,1,1) + 
     &    1/Pi*(Alfas*(2/3.D0*
     &          ((B0i(bb1,MU2,MGl2,MSf2(Sfe3,3,1)) + 
     &              MU2*B0i(dbb1,MU2,MGl2,MSf2(Sfe3,3,1)))*
     &            USf(Sfe3,1,3,1)*USfC(Sfe3,1,3,1)) + 
     &         1/MU*(1/3.D0*((2*MU*MU2*
     &                  (B0i(dbb1,MU2,MGl2,MSf2(Sfe3,3,1))*
     &                    USf(Sfe3,2,3,1)) + 
     &                 Mino3C*
     &                  (B0i(bb0,MU2,MGl2,MSf2(Sfe3,3,1)) + 
     &                    2*MU2*B0i(dbb0,MU2,MGl2,MSf2(Sfe3,3,1)))*
     &                  USf(Sfe3,1,3,1))*USfC(Sfe3,2,3,1)) + 
     &            Mino3*(-(1/3.D0*
     &                  B0i(bb0,MU2,MGl2,MSf2(Sfe3,3,1))) + 
     &               (2*MU2)/3.D0*B0i(dbb0,MU2,MGl2,MSf2(Sfe3,3,1)))*
     &             USf(Sfe3,2,3,1)*USfC(Sfe3,1,3,1))))

        ENDLOOP(Sfe3)

#ifdef DEBUGQ
        DEB("dZbarfL1(3,1,1) =", dZbarfL1(3,1,1))
#endif

        dZbarfR1(3,1,1) = dup9

        LOOP(Sfe3, 1,2,1)

        dZbarfR1(3,1,1) = 
     &   dZbarfR1(3,1,1) + 
     &    1/Pi*(Alfas*(-(1/(3.D0*MU)*
     &            ((-(2*MU*
     &                   ((B0i(bb1,MU2,MGl2,MSf2(Sfe3,3,1)) + 
     &                       MU2*B0i(dbb1,MU2,MGl2,MSf2(Sfe3,3,1)))
     &                      *USf(Sfe3,2,3,1))) + 
     &                Mino3C*
     &                 (B0i(bb0,MU2,MGl2,MSf2(Sfe3,3,1)) - 
     &                   2*MU2*B0i(dbb0,MU2,MGl2,MSf2(Sfe3,3,1)))*
     &                 USf(Sfe3,1,3,1))*USfC(Sfe3,2,3,1))) + 
     &         ((2*MU2)/3.D0*
     &             (B0i(dbb1,MU2,MGl2,MSf2(Sfe3,3,1))*
     &               USf(Sfe3,1,3,1)) + 
     &            1/MU*(Mino3*
     &               (1/3.D0*B0i(bb0,MU2,MGl2,MSf2(Sfe3,3,1)) + 
     &                 (2*MU2)/3.D0*B0i(dbb0,MU2,MGl2,MSf2(Sfe3,3,1))
     &                 )*USf(Sfe3,2,3,1)))*USfC(Sfe3,1,3,1)))

        ENDLOOP(Sfe3)

#ifdef DEBUGQ
        DEB("dZbarfR1(3,1,1) =", dZbarfR1(3,1,1))
#endif

        dup10 = -MSf2(1,3,1) + MSf2(2,3,1)

        dup11 = -A0(MGl2) + B0i(bb1,MSf2(2,3,1),MGl2,MU2)*MSf2(2,3,1)

        dup12 = -A0(MGl2) - 
     &    B0i(bb0,MSf2(2,3,1),MGl2,MU2)*(MU2 - MSf2(2,3,1)) + 
     &    B0i(bb1,MSf2(2,3,1),MGl2,MU2)*MSf2(2,3,1)

        dZbarSf1(1,2,3,1) = 
     &   (2*dMSfsq1(2,1,3,1) - 
     &      8/(3.D0*Pi)*(Alfas*
     &         ((MU*(Mino3*B0i(bb0,MSf2(2,3,1),MGl2,MU2)*
     &                 USf(2,2,3,1)) + dup12*USf(2,1,3,1))*
     &            USfC(1,1,3,1) + 
     &           (dup8*B0i(bb0,MSf2(2,3,1),MGl2,MU2) + 
     &              dup11*USf(2,2,3,1))*USfC(1,2,3,1))))/dup10

        LOOP(Sfe3, 1,2,1)

        dZbarSf1(1,2,3,1) = 
     &   dZbarSf1(1,2,3,1) - 
     &    2/(3.D0*Pi)*(Alfas*A0(MSf2(Sfe3,3,1))*
     &        (USf(Sfe3,1,3,1)*USfC(1,1,3,1) - 
     &          USf(Sfe3,2,3,1)*USfC(1,2,3,1))*
     &        (USf(2,1,3,1)*USfC(Sfe3,1,3,1) - 
     &          USf(2,2,3,1)*USfC(Sfe3,2,3,1)))/dup10

        ENDLOOP(Sfe3)

#ifdef DEBUGQ
        DEB("dZbarSf1(1,2,3,1) =", dZbarSf1(1,2,3,1))
#endif

        dup13 = MSf2(1,3,1) - MSf2(2,3,1)

        dup14 = -A0(MGl2) + B0i(bb1,MSf2(1,3,1),MGl2,MU2)*MSf2(1,3,1)

        dup15 = -A0(MGl2) - 
     &    B0i(bb0,MSf2(1,3,1),MGl2,MU2)*(MU2 - MSf2(1,3,1)) + 
     &    B0i(bb1,MSf2(1,3,1),MGl2,MU2)*MSf2(1,3,1)

        dZbarSf1(2,1,3,1) = 
     &   (2*dMSfsq1(1,2,3,1) - 
     &      8/(3.D0*Pi)*(Alfas*
     &         ((MU*(Mino3*B0i(bb0,MSf2(1,3,1),MGl2,MU2)*
     &                 USf(1,2,3,1)) + dup15*USf(1,1,3,1))*
     &            USfC(2,1,3,1) + 
     &           (dup2*B0i(bb0,MSf2(1,3,1),MGl2,MU2) + 
     &              dup14*USf(1,2,3,1))*USfC(2,2,3,1))))/dup13

        LOOP(Sfe3, 1,2,1)

        dZbarSf1(2,1,3,1) = 
     &   dZbarSf1(2,1,3,1) - 
     &    2/(3.D0*Pi)*(Alfas*A0(MSf2(Sfe3,3,1))*
     &        (USf(Sfe3,1,3,1)*USfC(2,1,3,1) - 
     &          USf(Sfe3,2,3,1)*USfC(2,2,3,1))*
     &        (USf(1,1,3,1)*USfC(Sfe3,1,3,1) - 
     &          USf(1,2,3,1)*USfC(Sfe3,2,3,1)))/dup13

        ENDLOOP(Sfe3)

#ifdef DEBUGQ
        DEB("dZbarSf1(2,1,3,1) =", dZbarSf1(2,1,3,1))
#endif

        dZfL1(3,1,1) = dup9

        LOOP(Sfe3, 1,2,1)

        dZfL1(3,1,1) = 
     &   dZfL1(3,1,1) + 1/Pi*
     &     (Alfas*(2/3.D0*((B0i(bb1,MU2,MGl2,MSf2(Sfe3,3,1)) + 
     &              MU2*B0i(dbb1,MU2,MGl2,MSf2(Sfe3,3,1)))*
     &            USf(Sfe3,1,3,1)*USfC(Sfe3,1,3,1)) + 
     &         1/MU*(-(1/3.D0*
     &               ((-(2*MU*MU2*
     &                      (B0i(dbb1,MU2,MGl2,MSf2(Sfe3,3,1))*
     &                       USf(Sfe3,2,3,1))) + 
     &                   Mino3C*
     &                    (B0i(bb0,MU2,MGl2,MSf2(Sfe3,3,1)) - 
     &                      2*MU2*B0i(dbb0,MU2,MGl2,MSf2(Sfe3,3,1))
     &                      )*USf(Sfe3,1,3,1))*USfC(Sfe3,2,3,1)))+
     &              Mino3*
     &             (1/3.D0*B0i(bb0,MU2,MGl2,MSf2(Sfe3,3,1)) + 
     &               (2*MU2)/3.D0*B0i(dbb0,MU2,MGl2,MSf2(Sfe3,3,1)))*
     &             USf(Sfe3,2,3,1)*USfC(Sfe3,1,3,1))))

        ENDLOOP(Sfe3)

#ifdef DEBUGQ
        DEB("dZfL1(3,1,1) =", dZfL1(3,1,1))
#endif

        dZfR1(3,1,1) = dup9

        LOOP(Sfe3, 1,2,1)

        dZfR1(3,1,1) = 
     &   dZfR1(3,1,1) + 1/Pi*
     &     (Alfas*(1/(3.D0*MU)*
     &          ((2*MU*((B0i(bb1,MU2,MGl2,MSf2(Sfe3,3,1)) + 
     &                   MU2*B0i(dbb1,MU2,MGl2,MSf2(Sfe3,3,1)))*
     &                 USf(Sfe3,2,3,1)) + 
     &              Mino3C*
     &               (B0i(bb0,MU2,MGl2,MSf2(Sfe3,3,1)) + 
     &                 2*MU2*B0i(dbb0,MU2,MGl2,MSf2(Sfe3,3,1)))*
     &               USf(Sfe3,1,3,1))*USfC(Sfe3,2,3,1)) + 
     &         ((2*MU2)/3.D0*
     &             (B0i(dbb1,MU2,MGl2,MSf2(Sfe3,3,1))*
     &               USf(Sfe3,1,3,1)) + 
     &            1/MU*(Mino3*
     &               (-(1/3.D0*B0i(bb0,MU2,MGl2,MSf2(Sfe3,3,1))) + 
     &                 (2*MU2)/3.D0*B0i(dbb0,MU2,MGl2,MSf2(Sfe3,3,1))
     &                 )*USf(Sfe3,2,3,1)))*USfC(Sfe3,1,3,1)))

        ENDLOOP(Sfe3)

#ifdef DEBUGQ
        DEB("dZfR1(3,1,1) =", dZfR1(3,1,1))
#endif

        dZSf1(1,1,3,1) = 
     &   1/(3.D0*Pi)*(Alfas*(3*B0i(bb0,MSf2(1,3,1),0.D0,MSf2(1,3,1)) + 
     &        2*(B0i(bb1,MSf2(1,3,1),0.D0,MSf2(1,3,1)) + 
     &           (2*B0i(dbb0,MSf2(1,3,1),0.D0,MSf2(1,3,1)) + 
     &              B0i(dbb1,MSf2(1,3,1),0.D0,MSf2(1,3,1)))*
     &            MSf2(1,3,1) - 
     &           2*((MU*(Mino3*B0i(dbb0,MSf2(1,3,1),MGl2,MU2)*
     &                    USf(1,2,3,1)) + 
     &                 (B0i(bb0,MSf2(1,3,1),MGl2,MU2) + 
     &                    B0i(bb1,MSf2(1,3,1),MGl2,MU2) + 
     &                    B0i(dbb1,MSf2(1,3,1),MGl2,MU2)*
     &                     MSf2(1,3,1) + 
     &                    B0i(dbb0,MSf2(1,3,1),MGl2,MU2)*
     &                     (-MU2 + MSf2(1,3,1)))*USf(1,1,3,1))*
     &               USfC(1,1,3,1) + 
     &              (dup2*B0i(dbb0,MSf2(1,3,1),MGl2,MU2) + 
     &                 (B0i(bb0,MSf2(1,3,1),MGl2,MU2) + 
     &                    B0i(bb1,MSf2(1,3,1),MGl2,MU2) + 
     &                    B0i(dbb1,MSf2(1,3,1),MGl2,MU2)*
     &                     MSf2(1,3,1))*USf(1,2,3,1))*USfC(1,2,3,1)
     &              ))))

#ifdef DEBUGQ
        DEB("dZSf1(1,1,3,1) =", dZSf1(1,1,3,1))
#endif

        dZSf1(1,2,3,1) = 
     &   -((2*dMSfsq1(1,2,3,1) - 
     &        8/(3.D0*Pi)*(Alfas*
     &           ((MU*(Mino3*B0i(bb0,MSf2(2,3,1),MGl2,MU2)*
     &                   USf(1,2,3,1)) + dup12*USf(1,1,3,1))*
     &              USfC(2,1,3,1) + 
     &             (dup4*B0i(bb0,MSf2(2,3,1),MGl2,MU2) + 
     &                dup11*USf(1,2,3,1))*USfC(2,2,3,1))))/dup13)

        LOOP(Sfe3, 1,2,1)

        dZSf1(1,2,3,1) = 
     &   dZSf1(1,2,3,1) + 2/(3.D0*Pi)*
     &     (Alfas*A0(MSf2(Sfe3,3,1))*
     &        (USf(Sfe3,1,3,1)*USfC(2,1,3,1) - 
     &          USf(Sfe3,2,3,1)*USfC(2,2,3,1))*
     &        (USf(1,1,3,1)*USfC(Sfe3,1,3,1) - 
     &          USf(1,2,3,1)*USfC(Sfe3,2,3,1)))/dup13

        ENDLOOP(Sfe3)

#ifdef DEBUGQ
        DEB("dZSf1(1,2,3,1) =", dZSf1(1,2,3,1))
#endif

        dZSf1(2,1,3,1) = 
     &   -((2*dMSfsq1(2,1,3,1) - 
     &        8/(3.D0*Pi)*(Alfas*
     &           ((MU*(Mino3*B0i(bb0,MSf2(1,3,1),MGl2,MU2)*
     &                   USf(2,2,3,1)) + dup15*USf(2,1,3,1))*
     &              USfC(1,1,3,1) + 
     &             (dup7*B0i(bb0,MSf2(1,3,1),MGl2,MU2) + 
     &                dup14*USf(2,2,3,1))*USfC(1,2,3,1))))/dup10)

        LOOP(Sfe3, 1,2,1)

        dZSf1(2,1,3,1) = 
     &   dZSf1(2,1,3,1) + 2/(3.D0*Pi)*
     &     (Alfas*A0(MSf2(Sfe3,3,1))*
     &        (USf(Sfe3,1,3,1)*USfC(1,1,3,1) - 
     &          USf(Sfe3,2,3,1)*USfC(1,2,3,1))*
     &        (USf(2,1,3,1)*USfC(Sfe3,1,3,1) - 
     &          USf(2,2,3,1)*USfC(Sfe3,2,3,1)))/dup10

        ENDLOOP(Sfe3)

#ifdef DEBUGQ
        DEB("dZSf1(2,1,3,1) =", dZSf1(2,1,3,1))
#endif

        dZSf1(2,2,3,1) = 
     &   1/(3.D0*Pi)*(Alfas*(3*B0i(bb0,MSf2(2,3,1),0.D0,MSf2(2,3,1)) + 
     &        2*(B0i(bb1,MSf2(2,3,1),0.D0,MSf2(2,3,1)) + 
     &           (2*B0i(dbb0,MSf2(2,3,1),0.D0,MSf2(2,3,1)) + 
     &              B0i(dbb1,MSf2(2,3,1),0.D0,MSf2(2,3,1)))*
     &            MSf2(2,3,1) - 
     &           2*((MU*(Mino3*B0i(dbb0,MSf2(2,3,1),MGl2,MU2)*
     &                    USf(2,2,3,1)) + 
     &                 (B0i(bb0,MSf2(2,3,1),MGl2,MU2) + 
     &                    B0i(bb1,MSf2(2,3,1),MGl2,MU2) + 
     &                    B0i(dbb1,MSf2(2,3,1),MGl2,MU2)*
     &                     MSf2(2,3,1) + 
     &                    B0i(dbb0,MSf2(2,3,1),MGl2,MU2)*
     &                     (-MU2 + MSf2(2,3,1)))*USf(2,1,3,1))*
     &               USfC(2,1,3,1) + 
     &              (dup8*B0i(dbb0,MSf2(2,3,1),MGl2,MU2) + 
     &                 (B0i(bb0,MSf2(2,3,1),MGl2,MU2) + 
     &                    B0i(bb1,MSf2(2,3,1),MGl2,MU2) + 
     &                    B0i(dbb1,MSf2(2,3,1),MGl2,MU2)*
     &                     MSf2(2,3,1))*USf(2,2,3,1))*USfC(2,2,3,1)
     &              ))))

#ifdef DEBUGQ
        DEB("dZSf1(2,2,3,1) =", dZSf1(2,2,3,1))
#endif

        dZbarSf1(1,1,3,1) = dZSf1(1,1,3,1)

#ifdef DEBUGQ
        DEB("dZbarSf1(1,1,3,1) =", dZbarSf1(1,1,3,1))
#endif

        dZbarSf1(2,2,3,1) = dZSf1(2,2,3,1)

#ifdef DEBUGQ
        DEB("dZbarSf1(2,2,3,1) =", dZbarSf1(2,2,3,1))
#endif

c######################## ddbar renormalization constants ##############
        ! reset dup's for ddbar CalcRenConst
        dup1  = CNull
        dup10 = CNull
        dup11 = CNull
        dup12 = CNull
        dup13 = CNull
        dup14 = CNull
        dup15 = CNull
        dup16 = CNull
        dup17 = CNull
        dup18 = CNull
        dup19 = CNull
        dup2  = CNull
        dup20 = CNull
        dup21 = CNull
        dup22 = CNull
        dup23 = CNull
        dup24 = CNull
        dup25 = CNull
        dup26 = CNull
        dup27 = CNull
        dup28 = CNull
        dup3  = CNull
        dup4  = CNull
        dup5  = CNull
        dup6  = CNull
        dup7  = CNull
        dup8  = CNull
        dup9  = CNull

        dMf1(4,1) = MD/Pi*
     &    (Alfas*(1/3.D0*Finite*xm0 - 2/3.D0*Re(B0i(bb0,MD2,0.D0,MD2)) + 
     &        2/3.D0*Re(B0i(bb1,MD2,0.D0,MD2))))

        LOOP(Sfe3, 1,2,1)

        dMf1(4,1) = dMf1(4,1) - 
     &    1/(3.D0*Pi)*(Alfas*
     &       (MD*(Re(B0i(bb1,MD2,MGl2,MSf2(Sfe3,4,1)))*
     &            (USf(Sfe3,1,4,1)*USfC(Sfe3,1,4,1) + 
     &              USf(Sfe3,2,4,1)*USfC(Sfe3,2,4,1))) + 
     &         Re(B0i(bb0,MD2,MGl2,MSf2(Sfe3,4,1)))*
     &          (Mino3*USf(Sfe3,2,4,1)*USfC(Sfe3,1,4,1) + 
     &            Mino3C*USf(Sfe3,1,4,1)*USfC(Sfe3,2,4,1))))

        ENDLOOP(Sfe3)

#ifdef DEBUGQ
        DEB("dMf1(4,1) =", dMf1(4,1))
#endif

        dup4 = -(MD*(Mino3C*USf(2,1,4,1))) + MD2*USf(2,2,4,1)

        dup5 = MD2*USf(2,1,4,1) - MD*(Mino3*USf(2,2,4,1))

        dup6 = USf(2,1,4,1)*USfC(2,1,4,1) + USf(2,2,4,1)*USfC(2,2,4,1)

        dMSfsq1(2,2,4,1) = 
     &   -(2/(3.D0*Pi)*(Alfas*
     &        (MSf2(2,4,1)*
     &           (2*Re(B0i(bb0,MSf2(2,4,1),0.D0,MSf2(2,4,1))) + 
     &             Re(B0i(bb1,MSf2(2,4,1),0.D0,MSf2(2,4,1)))) + 
     &          2*(dup6*(Re(A0(MGl2)) + 
     &                MSf2(2,4,1)*Re(B0i(bb1,MSf2(2,4,1),MD2,MGl2))
     &                ) + Re(B0i(bb0,MSf2(2,4,1),MD2,MGl2))*
     &              (dup5*USfC(2,1,4,1) + dup4*USfC(2,2,4,1))))))

        LOOP(Sfe3, 1,2,1)

        dMSfsq1(2,2,4,1) = 
     &   dMSfsq1(2,2,4,1) + 
     &    1/(3.D0*Pi)*(Alfas*Re(A0(MSf2(Sfe3,4,1)))*
     &       (USf(Sfe3,1,4,1)*USfC(2,1,4,1) - 
     &         USf(Sfe3,2,4,1)*USfC(2,2,4,1))*
     &       (USf(2,1,4,1)*USfC(Sfe3,1,4,1) - 
     &         USf(2,2,4,1)*USfC(Sfe3,2,4,1)))

        ENDLOOP(Sfe3)

#ifdef DEBUGQ
        DEB("dMSfsq1(2,2,4,1) =", dMSfsq1(2,2,4,1))
#endif

        dup7 = 1/Pi*(Alfas*
     &      (1/3.D0*Finite*xz0 - 
     &        2/3.D0*(B0i(bb0,MD2,0.D0,MD2) + B0i(bb1,MD2,0.D0,MD2)) + 
     &        MD2*(4/3.D0*B0i(dbb0,MD2,0.D0,MD2) - 
     &           4/3.D0*B0i(dbb1,MD2,0.D0,MD2))))

        dZbarfL1(4,1,1) = dup7

        LOOP(Sfe3, 1,2,1)

        dZbarfL1(4,1,1) = 
     &   dZbarfL1(4,1,1) + 
     &    1/Pi*(Alfas*(2/3.D0*
     &          ((B0i(bb1,MD2,MGl2,MSf2(Sfe3,4,1)) + 
     &              MD2*B0i(dbb1,MD2,MGl2,MSf2(Sfe3,4,1)))*
     &            USf(Sfe3,1,4,1)*USfC(Sfe3,1,4,1)) + 
     &         1/MD*(1/3.D0*((2*MD*MD2*
     &                  (B0i(dbb1,MD2,MGl2,MSf2(Sfe3,4,1))*
     &                    USf(Sfe3,2,4,1)) + 
     &                 Mino3C*
     &                  (B0i(bb0,MD2,MGl2,MSf2(Sfe3,4,1)) + 
     &                    2*MD2*B0i(dbb0,MD2,MGl2,MSf2(Sfe3,4,1)))*
     &                  USf(Sfe3,1,4,1))*USfC(Sfe3,2,4,1)) + 
     &            Mino3*(-(1/3.D0*
     &                  B0i(bb0,MD2,MGl2,MSf2(Sfe3,4,1))) + 
     &               (2*MD2)/3.D0*B0i(dbb0,MD2,MGl2,MSf2(Sfe3,4,1)))*
     &             USf(Sfe3,2,4,1)*USfC(Sfe3,1,4,1))))

        ENDLOOP(Sfe3)

#ifdef DEBUGQ
        DEB("dZbarfL1(4,1,1) =", dZbarfL1(4,1,1))
#endif

        dZbarfR1(4,1,1) = dup7

        LOOP(Sfe3, 1,2,1)

        dZbarfR1(4,1,1) = 
     &   dZbarfR1(4,1,1) + 
     &    1/Pi*(Alfas*(-(1/(3.D0*MD)*
     &            ((-(2*MD*
     &                   ((B0i(bb1,MD2,MGl2,MSf2(Sfe3,4,1)) + 
     &                       MD2*B0i(dbb1,MD2,MGl2,MSf2(Sfe3,4,1)))
     &                      *USf(Sfe3,2,4,1))) + 
     &                Mino3C*
     &                 (B0i(bb0,MD2,MGl2,MSf2(Sfe3,4,1)) - 
     &                   2*MD2*B0i(dbb0,MD2,MGl2,MSf2(Sfe3,4,1)))*
     &                 USf(Sfe3,1,4,1))*USfC(Sfe3,2,4,1))) + 
     &         ((2*MD2)/3.D0*
     &             (B0i(dbb1,MD2,MGl2,MSf2(Sfe3,4,1))*
     &               USf(Sfe3,1,4,1)) + 
     &            1/MD*(Mino3*
     &               (1/3.D0*B0i(bb0,MD2,MGl2,MSf2(Sfe3,4,1)) + 
     &                 (2*MD2)/3.D0*B0i(dbb0,MD2,MGl2,MSf2(Sfe3,4,1))
     &                 )*USf(Sfe3,2,4,1)))*USfC(Sfe3,1,4,1)))

        ENDLOOP(Sfe3)

#ifdef DEBUGQ
        DEB("dZbarfR1(4,1,1) =", dZbarfR1(4,1,1))
#endif

        dZfL1(4,1,1) = dup7

        LOOP(Sfe3, 1,2,1)

        dZfL1(4,1,1) = 
     &   dZfL1(4,1,1) + 1/Pi*
     &     (Alfas*(2/3.D0*((B0i(bb1,MD2,MGl2,MSf2(Sfe3,4,1)) + 
     &              MD2*B0i(dbb1,MD2,MGl2,MSf2(Sfe3,4,1)))*
     &            USf(Sfe3,1,4,1)*USfC(Sfe3,1,4,1)) + 
     &         1/MD*(-(1/3.D0*
     &               ((-(2*MD*MD2*
     &                      (B0i(dbb1,MD2,MGl2,MSf2(Sfe3,4,1))*
     &                       USf(Sfe3,2,4,1))) + 
     &                   Mino3C*
     &                    (B0i(bb0,MD2,MGl2,MSf2(Sfe3,4,1)) - 
     &                      2*MD2*B0i(dbb0,MD2,MGl2,MSf2(Sfe3,4,1))
     &                      )*USf(Sfe3,1,4,1))*USfC(Sfe3,2,4,1)))+
     &              Mino3*
     &             (1/3.D0*B0i(bb0,MD2,MGl2,MSf2(Sfe3,4,1)) + 
     &               (2*MD2)/3.D0*B0i(dbb0,MD2,MGl2,MSf2(Sfe3,4,1)))*
     &             USf(Sfe3,2,4,1)*USfC(Sfe3,1,4,1))))

        ENDLOOP(Sfe3)

#ifdef DEBUGQ
        DEB("dZfL1(4,1,1) =", dZfL1(4,1,1))
#endif

        dZfR1(4,1,1) = dup7

        LOOP(Sfe3, 1,2,1)

        dZfR1(4,1,1) = 
     &   dZfR1(4,1,1) + 1/Pi*
     &     (Alfas*(1/(3.D0*MD)*
     &          ((2*MD*((B0i(bb1,MD2,MGl2,MSf2(Sfe3,4,1)) + 
     &                   MD2*B0i(dbb1,MD2,MGl2,MSf2(Sfe3,4,1)))*
     &                 USf(Sfe3,2,4,1)) + 
     &              Mino3C*
     &               (B0i(bb0,MD2,MGl2,MSf2(Sfe3,4,1)) + 
     &                 2*MD2*B0i(dbb0,MD2,MGl2,MSf2(Sfe3,4,1)))*
     &               USf(Sfe3,1,4,1))*USfC(Sfe3,2,4,1)) + 
     &         ((2*MD2)/3.D0*
     &             (B0i(dbb1,MD2,MGl2,MSf2(Sfe3,4,1))*
     &               USf(Sfe3,1,4,1)) + 
     &            1/MD*(Mino3*
     &               (-(1/3.D0*B0i(bb0,MD2,MGl2,MSf2(Sfe3,4,1))) + 
     &                 (2*MD2)/3.D0*B0i(dbb0,MD2,MGl2,MSf2(Sfe3,4,1))
     &                 )*USf(Sfe3,2,4,1)))*USfC(Sfe3,1,4,1)))

        ENDLOOP(Sfe3)

#ifdef DEBUGQ
        DEB("dZfR1(4,1,1) =", dZfR1(4,1,1))
#endif

        dup8 = -(MD*(Mino3C*USf(1,1,4,1))) + MD2*USf(1,2,4,1)

        dup9 = MD2*USf(1,1,4,1) - MD*(Mino3*USf(1,2,4,1))

        dup10 = USf(1,1,4,1)*USfC(1,1,4,1) + 
     &    USf(1,2,4,1)*USfC(1,2,4,1)

        dZSf1(1,1,4,1) = 
     &   1/(3.D0*Pi)*(Alfas*(3*B0i(bb0,MSf2(1,4,1),0.D0,MSf2(1,4,1)) + 
     &        2*B0i(bb1,MSf2(1,4,1),0.D0,MSf2(1,4,1)) + 
     &        (4*B0i(dbb0,MSf2(1,4,1),0.D0,MSf2(1,4,1)) + 
     &           2*B0i(dbb1,MSf2(1,4,1),0.D0,MSf2(1,4,1)))*
     &         MSf2(1,4,1) + 
     &        4*(dup10*(B0i(bb1,MSf2(1,4,1),MD2,MGl2) + 
     &              B0i(dbb1,MSf2(1,4,1),MD2,MGl2)*MSf2(1,4,1)) + 
     &           B0i(dbb0,MSf2(1,4,1),MD2,MGl2)*
     &            (dup9*USfC(1,1,4,1) + dup8*USfC(1,2,4,1)))))

#ifdef DEBUGQ
        DEB("dZSf1(1,1,4,1) =", dZSf1(1,1,4,1))
#endif

        dZSf1(2,2,4,1) = 
     &   1/(3.D0*Pi)*(Alfas*(3*B0i(bb0,MSf2(2,4,1),0.D0,MSf2(2,4,1)) + 
     &        2*B0i(bb1,MSf2(2,4,1),0.D0,MSf2(2,4,1)) + 
     &        (4*B0i(dbb0,MSf2(2,4,1),0.D0,MSf2(2,4,1)) + 
     &           2*B0i(dbb1,MSf2(2,4,1),0.D0,MSf2(2,4,1)))*
     &         MSf2(2,4,1) + 
     &        4*(dup6*(B0i(bb1,MSf2(2,4,1),MD2,MGl2) + 
     &              B0i(dbb1,MSf2(2,4,1),MD2,MGl2)*MSf2(2,4,1)) + 
     &           B0i(dbb0,MSf2(2,4,1),MD2,MGl2)*
     &            (dup5*USfC(2,1,4,1) + dup4*USfC(2,2,4,1)))))

#ifdef DEBUGQ
        DEB("dZSf1(2,2,4,1) =", dZSf1(2,2,4,1))
#endif

        dup13 = MSf2(1,4,1) - MSf2(2,4,1)

        dMsq11Sf1(4,1) = 
     &   -(C2B*dMWsq1) + Abs(USf(1,1,3,1))**2*dMSfsq1(1,1,3,1) + 
     &    Abs(USf(1,2,3,1))**2*dMSfsq1(2,2,3,1) + 
     &    2*(MW2*(CB2*dTB1*S2B) + MD*dMf1(4,1)) - 
     &    2*(MU*dMf1(3,1) + 
     &       Re(dMSfsq1(1,2,3,1)*USf(2,2,3,1)*USfC(1,2,3,1)))

#ifdef DEBUGQ
        DEB("dMsq11Sf1(4,1) =", dMsq11Sf1(4,1))
#endif

        dMsq12Sf1(4,1) = 
     &   -((MUE*TB - AfC(4,1,1))*dMf1(4,1)) - 
     &    MD*(dTB1*MUE + dMUE1*TB - Conjugate(dAf1(4,1,1)))

#ifdef DEBUGQ
        print*,"dMf1",dMf1(4,1)
        DEB("dMsq12Sf1(4,1) =", dMsq12Sf1(4,1))
#endif

        dMSfsq1(1,1,4,1) = 
     &   (Abs(USf(1,2,4,1))**2*
     &       (dMSfsq1(2,2,4,1) - dMsq11Sf1(4,1)) + 
     &      Abs(USf(1,1,4,1))**2*dMsq11Sf1(4,1) + 
     &      2*Re(dMsq12Sf1(4,1)*USf(1,1,4,1)*USfC(1,2,4,1)))/
     &    Abs(USf(1,1,4,1))**2

#ifdef DEBUGQ
        DEB("dMSfsq1(1,1,4,1) =", dMSfsq1(1,1,4,1))
#endif

        dMSfsq1(1,2,4,1) = 
     &   -((Conjugate(dMsq12Sf1(4,1))*USf(1,2,4,1)*USfC(2,1,4,1) - 
     &        USf(1,1,4,1)*
     &         ((dMSfsq1(1,1,4,1) - dMSfsq1(2,2,4,1))*
     &            USfC(2,1,4,1) + dMsq12Sf1(4,1)*USfC(2,2,4,1)))/
     &      (Abs(USf(1,1,4,1))**2 - Abs(USf(1,2,4,1))**2))

#ifdef DEBUGQ
        DEB("dMSfsq1(1,2,4,1) =", dMSfsq1(1,2,4,1))
#endif

        dMSfsq1(2,1,4,1) = Conjugate(dMSfsq1(1,2,4,1))

#ifdef DEBUGQ
        DEB("dMSfsq1(2,1,4,1) =", dMSfsq1(2,1,4,1))
#endif

        dZbarSf1(1,1,4,1) = dZSf1(1,1,4,1)

#ifdef DEBUGQ
        DEB("dZbarSf1(1,1,4,1) =", dZbarSf1(1,1,4,1))
#endif

        dup15 = -MSf2(1,4,1) + MSf2(2,4,1)

        dup16 = USf(2,1,4,1)*USfC(1,1,4,1) + 
     &    USf(2,2,4,1)*USfC(1,2,4,1)

        dup17 = A0(MGl2) + B0i(bb1,MSf2(2,4,1),MD2,MGl2)*MSf2(2,4,1)

        dZbarSf1(1,2,4,1) = 
     &   (2*dMSfsq1(2,1,4,1) + 
     &      8/(3.D0*Pi)*(Alfas*
     &         (dup16*dup17 + 
     &           B0i(bb0,MSf2(2,4,1),MD2,MGl2)*
     &            (dup5*USfC(1,1,4,1) + dup4*USfC(1,2,4,1)))))/
     &    dup15

        LOOP(Sfe3, 1,2,1)

        dZbarSf1(1,2,4,1) = 
     &   dZbarSf1(1,2,4,1) - 
     &    2/(3.D0*Pi)*(Alfas*A0(MSf2(Sfe3,4,1))*
     &        (USf(Sfe3,1,4,1)*USfC(1,1,4,1) - 
     &          USf(Sfe3,2,4,1)*USfC(1,2,4,1))*
     &        (USf(2,1,4,1)*USfC(Sfe3,1,4,1) - 
     &          USf(2,2,4,1)*USfC(Sfe3,2,4,1)))/dup15

        ENDLOOP(Sfe3)

#ifdef DEBUGQ
        DEB("dZbarSf1(1,2,4,1) =", dZbarSf1(1,2,4,1))
#endif

        dup18 = USf(1,1,4,1)*USfC(2,1,4,1) + 
     &    USf(1,2,4,1)*USfC(2,2,4,1)

        dup19 = A0(MGl2) + B0i(bb1,MSf2(1,4,1),MD2,MGl2)*MSf2(1,4,1)

        dZbarSf1(2,1,4,1) = 
     &   (2*dMSfsq1(1,2,4,1) + 
     &      8/(3.D0*Pi)*(Alfas*
     &         (dup18*dup19 + 
     &           B0i(bb0,MSf2(1,4,1),MD2,MGl2)*
     &            (dup9*USfC(2,1,4,1) + dup8*USfC(2,2,4,1)))))/
     &    dup13

        LOOP(Sfe3, 1,2,1)

        dZbarSf1(2,1,4,1) = 
     &   dZbarSf1(2,1,4,1) - 
     &    2/(3.D0*Pi)*(Alfas*A0(MSf2(Sfe3,4,1))*
     &        (USf(Sfe3,1,4,1)*USfC(2,1,4,1) - 
     &          USf(Sfe3,2,4,1)*USfC(2,2,4,1))*
     &        (USf(1,1,4,1)*USfC(Sfe3,1,4,1) - 
     &          USf(1,2,4,1)*USfC(Sfe3,2,4,1)))/dup13

        ENDLOOP(Sfe3)

#ifdef DEBUGQ
        DEB("dZbarSf1(2,1,4,1) =", dZbarSf1(2,1,4,1))
#endif

        dZbarSf1(2,2,4,1) = dZSf1(2,2,4,1)

#ifdef DEBUGQ
        DEB("dZbarSf1(2,2,4,1) =", dZbarSf1(2,2,4,1))
#endif

        dZSf1(1,2,4,1) = 
     &   -((2*dMSfsq1(1,2,4,1) + 
     &        8/(3.D0*Pi)*(Alfas*
     &           (dup17*dup18 + 
     &             B0i(bb0,MSf2(2,4,1),MD2,MGl2)*
     &              (dup9*USfC(2,1,4,1) + dup8*USfC(2,2,4,1)))))/
     &      dup13)

        LOOP(Sfe3, 1,2,1)

        dZSf1(1,2,4,1) = 
     &   dZSf1(1,2,4,1) + 2/(3.D0*Pi)*
     &     (Alfas*A0(MSf2(Sfe3,4,1))*
     &        (USf(Sfe3,1,4,1)*USfC(2,1,4,1) - 
     &          USf(Sfe3,2,4,1)*USfC(2,2,4,1))*
     &        (USf(1,1,4,1)*USfC(Sfe3,1,4,1) - 
     &          USf(1,2,4,1)*USfC(Sfe3,2,4,1)))/dup13

        ENDLOOP(Sfe3)

#ifdef DEBUGQ
        DEB("dZSf1(1,2,4,1) =", dZSf1(1,2,4,1))
#endif

        dZSf1(2,1,4,1) = 
     &   -((2*dMSfsq1(2,1,4,1) + 
     &        8/(3.D0*Pi)*(Alfas*
     &           (dup16*dup19 + 
     &             B0i(bb0,MSf2(1,4,1),MD2,MGl2)*
     &              (dup5*USfC(1,1,4,1) + dup4*USfC(1,2,4,1)))))/
     &      dup15)

        LOOP(Sfe3, 1,2,1)

        dZSf1(2,1,4,1) = 
     &   dZSf1(2,1,4,1) + 2/(3.D0*Pi)*
     &     (Alfas*A0(MSf2(Sfe3,4,1))*
     &        (USf(Sfe3,1,4,1)*USfC(1,1,4,1) - 
     &          USf(Sfe3,2,4,1)*USfC(1,2,4,1))*
     &        (USf(2,1,4,1)*USfC(Sfe3,1,4,1) - 
     &          USf(2,2,4,1)*USfC(Sfe3,2,4,1)))/dup15

        ENDLOOP(Sfe3)

#ifdef DEBUGQ
        DEB("dZSf1(2,1,4,1) =", dZSf1(2,1,4,1))
#endif

#include "contains.h"

        end

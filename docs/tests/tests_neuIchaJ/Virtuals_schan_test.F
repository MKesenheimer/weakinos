c############### Virtuals.f ############################################
c last modified by MK, 10.02.2015
c neutralino pair production
c virtual corrections at NLO SQCD:
c parton parton -> neutralino neutralino

c############### subroutine setvirtual #################################
c returns (M_B * M_V)/(as/(2pi)) as real number,
c where M_B is the Born amplitude and
c M_V is the finite part of the virtual amplitude
c The as/(2pi) factor is attached at a later point
c A factor (4 Pi)^ep/Gamma(1-ep) with ep=(4-d)/2 is implied.
      subroutine setvirtual(p,vflav,ampvirt)
        implicit none

#include "PhysPars.h"
#include "pwhg_st.h"
#include "pwhg_math.h"
#include "nlegborn.h"
#include "looptools.h"
#include "Flags.h"

        double precision p(0:3,nlegborn)
        integer vflav(nlegborn)
        ! copy of vflav and p
        double precision pc(0:3,nlegborn)
        integer flavc(nlegborn)
        ! results. Note: bornjk and bmunu is only used to be able to
        ! call the MadGraph Born amplitude, when using fake virtuals.
        double precision born, bornjk(nlegborn,nlegborn)
        double precision bmunu(0:3,0:3,nlegborn)
        double precision virtual, ampvirt
        ! formcalc results
        double precision fc_result(2)
        ! symmetry factor
        double precision symfac
        ! function to calculate the squared momentum sum
        double precision momsum2sq
        ! kinematics
        double precision s, t, u
        double precision MassQ1, MassQ2, MassF1, MassF2
        ! status
        logical verbose
        logical init_virtual
        data init_virtual/.true./
        ! divergences
        double precision ltlambda
        common/debug/ltlambda
        ! functions
        double precision getAnalyticAmp

#ifdef DEBUG
#define DEBUG1
        ! debugging
        integer counter
        data counter/0/
        save counter
#endif

        if (init_virtual) then
          print*
          print*, "=== Init virtual... ==="
          verbose = .true.
          init_virtual = .false.
        else
          verbose = .false.
        endif
        
        ! backup flavor list and momenta, so we can operate on it without
        ! affecting the other program
        flavc = vflav
        pc = p
        
#ifdef DEBUGQ
        ! overwrite when precompiler flag DEBUGQ is set
        verbose = .true.
#endif

        if (.not.( (abs(vflav(2)).eq.abs(vflav(1)))
     &      .or.   (abs(vflav(2)+vflav(1)).eq.1) )) then
          print*, "Error in fct. setvirt: invalid initial state"
          print*, "vflav(1) = ", vflav(1)
          print*, "vflav(2) = ", vflav(2)
          stop
        endif
        
        ! check four momentum conservation
        call check_4conservation(p, nlegborn)

        ! calculate kinematics
        s = momsum2sq(p(0:3,1), p(0:3,2))
        t = momsum2sq(p(0:3,1),-p(0:3,3))
        u = momsum2sq(p(0:3,2),-p(0:3,3))

        ! calculate the couplings for this PS point
        call set_ebe_couplings

        ! load process dependent masses
        call set_process(vflav(1:4), MassQ1, MassQ2, MassF1, MassF2)

        !call setversionkey(4)
        !call ltini
        ltlambda = 0d0
        call setlambda(ltlambda) ! Only finite part
        call setuvdiv(0d0)  ! no 1/eps_UV terms
        call setdelta(0d0)  ! set divergences to zero
        call setmudim(st_muren2) ! renormalization scale
        
#define CHECKUV
#ifdef CHECKUV
        ! test for finiteness
        ! print*, "checking for UV finiteness"
        call setuvdiv(1d0)
        call setDelta(1D7)    ! Test UV finiteness
#endif
        
        ! crossing (FormCalc provides only amplitudes for qqbar 
        ! and not for qbarq)
        if(vflav(1).gt.0) then
          p(:,1) = pc(:,1)
          p(:,2) = pc(:,2)
        else
          p(:,1) = pc(:,2)
          p(:,2) = pc(:,1)
        endif

        ! symmetry factor
        if(vflav(3).eq.vflav(4)) then
          ! 72 = 2*2*3*3*2 (spin factor 2, colour factor 3, equal particles in the final state 2)
          symfac = 72d0 ! take spin and colour factors into account
        else
          ! 36 = 2*2*3*3 (spin factor 2, colour factor 3)
          symfac = 36d0 ! take spin and colour factors into account
        endif

        ! reset the amplitudes
        fc_result(1) = 0
        fc_result(2) = 0

        if(((vflav(1).eq.-2 .and. vflav(2).eq.1) .or.
     &           (vflav(1).eq.1 .and. vflav(2).eq.-2)) .and. 
     &           (vflav(3).eq.1000022 .and. vflav(4).eq.-1000024) ) then
          fc_result(1) = getAnalyticAmp(s,t,u,1)
          fc_result(2) = getAnalyticAmp(s,t,u,2)
        else if(((vflav(1).eq.-4 .and. vflav(2).eq.3) .or.
     &           (vflav(1).eq.3 .and. vflav(2).eq.-4)) .and.
     &           (vflav(3).eq.1000022 .and. vflav(4).eq.-1000024) ) then
          fc_result(1) = getAnalyticAmp(s,t,u,1)
          fc_result(2) = getAnalyticAmp(s,t,u,2)
        else
          print*,"error: process not implemented yet",vflav
          stop
        endif


        ! helicity factors (for every external particle factor of 2)
        fc_result(1) = 16D0*fc_result(1)
        fc_result(2) = 16D0*fc_result(2)
        ! symmetry factors
        born    = fc_result(1)/symfac
        virtual = fc_result(2)/symfac ! = 2 x Re[M_tree* x M_1-loop], FC8Guid.pdf p.48
        ! need to extract a factor N=(4Pi)^eps/Ga(1-eps)
        ! -> LoopTools does already extract this factor
        ampvirt = virtual/(Alfas/(2d0*pi))
        !ampvirt = virtual*2d0*pi ! if FormCalc Alfas = 1


        if(isnan(virtual)) then
          print*,"got strange value for virtual amplitude...", virtual
          stop
        endif
        
        ! Print informations
        if (verbose) then
          print*, "process: ", vflav
          print*, "s     = ", s
          print*, "t     = ", t
          print*, "u     = ", u
          print*, "AlphaS         = ", Alfas
          print*, "muren2         = ", st_muren2
          print*, "born           = ", born
          print*, "2*Born*Virtual = ", ampvirt
          print*
        endif

        ! restore flavor list
        vflav = flavc
        p = pc
        
#ifdef DEBUG1
        if(counter < 5) then
          print*, "process: ", vflav
          print*, "s     = ", s
          print*, "t     = ", t
          print*, "u     = ", u
          print*,"virt = ", ampvirt
          counter = counter + 1
        endif
#endif

      end

c############### end subroutine setvirtual #############################

c############### function getAnalyticAmp ###############################
c returns analytic born and virtual amplitude calculated by FormCalc
      double precision function getAnalyticAmp(s, t, u, index)
        implicit none

#include "PhysPars.h"
#include "pwhg_st.h"
#include "pwhg_math.h"
#include "nlegborn.h"
#include "looptools.h"
#include "Flags.h"

        ! input
        double precision s, t, u
        integer index

        ! renConstants
        double complex dMWsq1, dMZsq1, dZAA1, dZZA1, dSW1, dZe1, dZW1
        double complex dZbarfL1(4,1,1),dZfL1(3,1,1)
        double precision Finite
        double precision cI
        parameter (cI = (0D0,1D0))

        Finite = 1

        dMWsq1 = dcmplx(0)
        dMZsq1 = dcmplx(0)
        dZAA1  = dcmplx(0)
        dZZA1  = dcmplx(0)
        dSW1   = dcmplx(0)
        dZe1   = dcmplx(0)
        dZW1   = dcmplx(0)

        MU2 = 0D0
        MD2 = 0D0

        dZbarfL1(4,1,1) = dcmplx((Alfas*
     -      (Finite - 2*B0i(bb0,MD2,0,MD2) - 2*B0i(bb1,MD2,0,MD2) +
     -        4*MD2*B0i(dbb0,MD2,0,MD2) - 4*MD2*B0i(dbb1,MD2,0,MD2)))/
     -    (3.*Pi))
        dZfL1(3,1,1) = dcmplx((Alfas*(Finite - 2*B0i(bb0,MU2,0,MU2) -
     -        2*B0i(bb1,MU2,0,MU2) + 4*MU2*B0i(dbb0,MU2,0,MU2) -
     -        4*MU2*B0i(dbb1,MU2,0,MU2)))/(3.*Pi))


        ! born amplitude for process p p > n1 x1-, only s-channels
        if(index.eq.1) then
          getAnalyticAmp = (-3*Alfa2*Pi**2*(-((T**2 + U**2
     -          + 2*MCha2(1)*MNeu2(1) -
     -           T*(MCha2(1) + MNeu2(1)) - U*(MCha2(1) + MNeu2(1)))*
     -         (Sqrt(2.)*UCha(1,1)*UCha(1,2)*ZNeu(1,2)*ZNeu(1,3) + 
     -           UCha(1,2)**2*ZNeu(1,3)**2 + 
     -           UCha(1,1)*ZNeu(1,2)*
     -            (2*UCha(1,1)*ZNeu(1,2) + 
     -              Sqrt(2.)*UCha(1,2)*ZNeu(1,3)) - 
     -           Sqrt(2.)*VCha(1,1)*VCha(1,2)*ZNeu(1,2)*ZNeu(1,4) + 
     -           VCha(1,2)**2*ZNeu(1,4)**2 + 
     -           VCha(1,1)*ZNeu(1,2)*
     -            (2*VCha(1,1)*ZNeu(1,2) - 
     -              Sqrt(2.)*VCha(1,2)*ZNeu(1,4)))) + 
     -      S**2*(Sqrt(2.)*UCha(1,1)*UCha(1,2)*ZNeu(1,2)*ZNeu(1,3) + 
     -         UCha(1,2)**2*ZNeu(1,3)**2 + 
     -         UCha(1,1)*ZNeu(1,2)*
     -          (2*UCha(1,1)*ZNeu(1,2) + Sqrt(2.)*UCha(1,2)*ZNeu(1,3))
     -           + Sqrt(2.)*VCha(1,1)*VCha(1,2)*ZNeu(1,2)*ZNeu(1,4) - 
     -         VCha(1,2)**2*ZNeu(1,4)**2 + 
     -         VCha(1,1)*ZNeu(1,2)*
     -          (-2*VCha(1,1)*ZNeu(1,2) + 
     -            Sqrt(2.)*VCha(1,2)*ZNeu(1,4))) - 
     -      S*(2*MCha(1)*MNeu(1)*
     -          (2*UCha(1,1)*VCha(1,1)*ZNeu(1,2)**2 + 
     -            Sqrt(2.)*UCha(1,2)*VCha(1,1)*ZNeu(1,2)*ZNeu(1,3) - 
     -            Sqrt(2.)*UCha(1,1)*VCha(1,2)*ZNeu(1,2)*ZNeu(1,4) - 
     -            UCha(1,2)*VCha(1,2)*ZNeu(1,3)*ZNeu(1,4) + 
     -            UCha(1,2)*ZNeu(1,3)*
     -             (Sqrt(2.)*VCha(1,1)*ZNeu(1,2) - 
     -               VCha(1,2)*ZNeu(1,4)) + 
     -            UCha(1,1)*ZNeu(1,2)*
     -             (2*VCha(1,1)*ZNeu(1,2) - 
     -               Sqrt(2.)*VCha(1,2)*ZNeu(1,4))) - 
     -         (2*U - MCha2(1) - MNeu2(1))*
     -          (Sqrt(2.)*UCha(1,1)*UCha(1,2)*ZNeu(1,2)*ZNeu(1,3) + 
     -            UCha(1,2)**2*ZNeu(1,3)**2 + 
     -            UCha(1,1)*ZNeu(1,2)*
     -             (2*UCha(1,1)*ZNeu(1,2) + 
     -               Sqrt(2.)*UCha(1,2)*ZNeu(1,3)) + 
     -            Sqrt(2.)*VCha(1,1)*VCha(1,2)*ZNeu(1,2)*ZNeu(1,4) - 
     -            VCha(1,2)**2*ZNeu(1,4)**2 + 
     -            VCha(1,1)*ZNeu(1,2)*
     -             (-2*VCha(1,1)*ZNeu(1,2) + 
     -               Sqrt(2.)*VCha(1,2)*ZNeu(1,4))))))/
     -  (2.*(-MW2 + S + cI*MW*WW)**2*SW2**2)

        ! virt amplitude for process p p > n1 x1-, only s-channels
        ! the SUSY counter term is only needed for t-channel diagrams
        ! (they have yukawa couplings)
        else
          getAnalyticAmp =
     -       (Alfa2*Pi*(6*dSW1*Pi + 4*Alfas*Finite*SW - 6*dZe1*Pi*SW -
     -      3*dZW1*Pi*SW - 4*Alfas*SW*B0i(bb0,S,MD2,MU2) + 
     -      4*Alfas*S*SW*C0i(cc0,MD2,S,MU2,0,MD2,MU2) + 
     -      8*Alfas*SW*C0i(cc00,MD2,S,MU2,0,MD2,MU2) + 
     -      4*Alfas*S*SW*C0i(cc1,MD2,S,MU2,0,MD2,MU2) + 
     -      4*Alfas*S*SW*C0i(cc2,MD2,S,MU2,0,MD2,MU2) + 
     -      6*Pi*Conjugate(dSW1) - 6*Pi*SW*Conjugate(dZe1) - 
     -      3*Pi*SW*Conjugate(dZW1) - 
     -      4*Alfas*SW*Conjugate(B0i(bb0,S,MD2,MU2)) + 
     -      4*Alfas*S*SW*Conjugate(C0i(cc0,MD2,S,MU2,0,MD2,MU2)) + 
     -      8*Alfas*SW*Conjugate(C0i(cc00,MD2,S,MU2,0,MD2,MU2)) + 
     -      4*Alfas*S*SW*Conjugate(C0i(cc1,MD2,S,MU2,0,MD2,MU2)) + 
     -      4*Alfas*S*SW*Conjugate(C0i(cc2,MD2,S,MU2,0,MD2,MU2)) - 
     -      3*Pi*SW*Conjugate(dZbarfL1(4,1,1)) - 
     -      3*Pi*SW*Conjugate(dZfL1(3,1,1)) - 
     -      3*Pi*SW*dZbarfL1(4,1,1) - 3*Pi*SW*dZfL1(3,1,1))*
     -    (-((T**2 + U**2 + 2*MCha2(1)*MNeu2(1) - 
     -           T*(MCha2(1) + MNeu2(1)) - U*(MCha2(1) + MNeu2(1)))*
     -         (Sqrt(2.)*UCha(1,1)*UCha(1,2)*ZNeu(1,2)*ZNeu(1,3) + 
     -           UCha(1,2)**2*ZNeu(1,3)**2 + 
     -           UCha(1,1)*ZNeu(1,2)*
     -            (2*UCha(1,1)*ZNeu(1,2) + 
     -              Sqrt(2.)*UCha(1,2)*ZNeu(1,3)) - 
     -           Sqrt(2.)*VCha(1,1)*VCha(1,2)*ZNeu(1,2)*ZNeu(1,4) + 
     -           VCha(1,2)**2*ZNeu(1,4)**2 + 
     -           VCha(1,1)*ZNeu(1,2)*
     -            (2*VCha(1,1)*ZNeu(1,2) - 
     -              Sqrt(2.)*VCha(1,2)*ZNeu(1,4)))) + 
     -      S**2*(Sqrt(2.)*UCha(1,1)*UCha(1,2)*ZNeu(1,2)*ZNeu(1,3) + 
     -         UCha(1,2)**2*ZNeu(1,3)**2 + 
     -         UCha(1,1)*ZNeu(1,2)*
     -          (2*UCha(1,1)*ZNeu(1,2) + Sqrt(2.)*UCha(1,2)*ZNeu(1,3))
     -           + Sqrt(2.)*VCha(1,1)*VCha(1,2)*ZNeu(1,2)*ZNeu(1,4) - 
     -         VCha(1,2)**2*ZNeu(1,4)**2 + 
     -         VCha(1,1)*ZNeu(1,2)*
     -          (-2*VCha(1,1)*ZNeu(1,2) + 
     -            Sqrt(2.)*VCha(1,2)*ZNeu(1,4))) - 
     -      S*(2*MCha(1)*MNeu(1)*
     -          (2*UCha(1,1)*VCha(1,1)*ZNeu(1,2)**2 + 
     -            Sqrt(2.)*UCha(1,2)*VCha(1,1)*ZNeu(1,2)*ZNeu(1,3) - 
     -            Sqrt(2.)*UCha(1,1)*VCha(1,2)*ZNeu(1,2)*ZNeu(1,4) - 
     -            UCha(1,2)*VCha(1,2)*ZNeu(1,3)*ZNeu(1,4) + 
     -            UCha(1,2)*ZNeu(1,3)*
     -             (Sqrt(2.)*VCha(1,1)*ZNeu(1,2) - 
     -               VCha(1,2)*ZNeu(1,4)) + 
     -            UCha(1,1)*ZNeu(1,2)*
     -             (2*VCha(1,1)*ZNeu(1,2) - 
     -               Sqrt(2.)*VCha(1,2)*ZNeu(1,4))) - 
     -         (2*U - MCha2(1) - MNeu2(1))*
     -          (Sqrt(2.)*UCha(1,1)*UCha(1,2)*ZNeu(1,2)*ZNeu(1,3) + 
     -            UCha(1,2)**2*ZNeu(1,3)**2 + 
     -            UCha(1,1)*ZNeu(1,2)*
     -             (2*UCha(1,1)*ZNeu(1,2) + 
     -               Sqrt(2.)*UCha(1,2)*ZNeu(1,3)) + 
     -            Sqrt(2.)*VCha(1,1)*VCha(1,2)*ZNeu(1,2)*ZNeu(1,4) - 
     -            VCha(1,2)**2*ZNeu(1,4)**2 + 
     -            VCha(1,1)*ZNeu(1,2)*
     -             (-2*VCha(1,1)*ZNeu(1,2) + 
     -               Sqrt(2.)*VCha(1,2)*ZNeu(1,4))))))/
     -  (4.*(-MW2 + S + cI*MW*WW)**2*SW*SW2**2)
        endif


      end

c############### end subroutine setvirtual #############################

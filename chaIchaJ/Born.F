c############### Born.f ################################################
c last modified by MK, 22.01.2015
c -6  -5  -4  -3  -2  -1  0  1  2  3  4  5  6
c t~  b~  c~  s~  u~  d~  g  d  u  s  c  b  t

c############### setborn subroutine ####################################
c sets Born cross section
c Born contribution:
c parton parton -> weakino weakino

      subroutine setborn(p,bflav,born,bornjk,bmunu)
        implicit none

#include "pwhg_math.h"
#include "nlegborn.h"
#include "PhysPars.h"

        integer bflav(nlegborn)
        double precision p(0:3,nlegborn), bornjk(nlegborn,nlegborn)
        double precision bmunu(0:3,0:3,nlegborn), born, colcf
        ! copy of bflav and p
        double precision pc(0:3,nlegborn)
        integer  flavc(nlegborn)
        ! formcalc results
        double precision fc_result(2)
        ! symmetry factor
        double precision symfac
        ! helicity and flags parameter for formcalc       
        integer*8 helicities
        integer flags
        ! indices
        integer j,k,nmu,nnu
        ! kinematics: mandelstam variables, quark masses, weakino masses
        double precision s, t, u, MassQ1, MassQ2, MassF1, MassF2
        ! check 4-momentum conservation
        logical lresult
        ! phase space
        double precision ex, ey, ez, psNorm
        ! function to calculate the squared momentum sum
        double precision momsq, momsum3sq, momsum2sq
        external momsq, momsum2sq, momsum3sq
        ! calculate the relative error between formcalc and madgraph borns
        double precision fborn, mborn
        ! initialization variables
        logical verbose
        logical init_born
        data init_born/.true./
#ifdef DEBUGQ
#define DEBUG1
        ! debugging
        integer counter
        data counter/0/
        save counter
#endif

        if (init_born) then
           print*
           print*, "=== Init born... ==="
           verbose = .true.
           init_born = .false.
        else
           verbose = .false.
        endif

        if (.not.( (abs(bflav(2)).eq.abs(bflav(1)))
     &      .or.   (abs(bflav(2)+bflav(1)).eq.1) )) then
          print*, "Error in fct. setborn: invalid initial state"
          print*, "bflav(1) = ", bflav(1)
          print*, "bflav(2) = ", bflav(2)
          stop
        endif

        ! check four momentum conservation
        call check_4conservation(p, nlegborn, 1, lresult)
        
        !if(verbose) then
          ! calculate kinematics only for debug purpose
          s = momsum2sq(p(0:3,1), p(0:3,2))
          t = momsum2sq(p(0:3,1),-p(0:3,3))
          u = momsum2sq(p(0:3,2),-p(0:3,3))
        !endif
               
#if !defined(FORM_BORN) && !defined(MAD_BORN)
        print*,"Error: Preprocessor flag FORM_BORN or"
     &       //" MAD_BORN not set in Makefile."
        stop
#endif
        
#if defined(CHECK_FORM_MAD) || defined(FORM_BORN)
c ---------------------------- Formcalc Born ---------------------------

        ! store flavor list, so we can operate on it without
        ! affecting the other program
        flavc = bflav
        pc = p

        ! calculate the couplings for this PS point
        call set_ebe_couplings

        ! load process dependent masses
        call set_process(bflav(1:4), MassQ1, MassQ2, MassF1, MassF2)
        
        ! crossing (FormCalc provides only amplitudes for qqbar 
        ! and not for qbarq)
        if(bflav(1).gt.0) then
          p(:,1) = pc(:,2)
          p(:,2) = pc(:,1)
        else
          p(:,1) = pc(:,1)
          p(:,2) = pc(:,2)
        endif
        
        ! init phase space
        ! calculate unit vectors
        psNorm=dsqrt(p(1,1)**2+p(2,1)**2+p(3,1)**2)
        ex=p(1,1)/psNorm
        ey=p(2,1)/psNorm
        ez=p(3,1)/psNorm
        call VecSet(1, MassQ1, psNorm, ex,ey,ez)
        psNorm=dsqrt(p(1,2)**2+p(2,2)**2+p(3,2)**2)
        ex=p(1,2)/psNorm
        ey=p(2,2)/psNorm
        ez=p(3,2)/psNorm
        call VecSet(2, MassQ2, psNorm, ex,ey,ez)
        psNorm=dsqrt(p(1,3)**2+p(2,3)**2+p(3,3)**2)
        ex=p(1,3)/psNorm
        ey=p(2,3)/psNorm
        ez=p(3,3)/psNorm
        call VecSet(3, MassF1, psNorm, ex,ey,ez)
        psNorm=dsqrt(p(1,4)**2+p(2,4)**2+p(3,4)**2)
        ex=p(1,4)/psNorm
        ey=p(2,4)/psNorm
        ez=p(3,4)/psNorm
        call VecSet(4, MassF2, psNorm, ex,ey,ez)

        ! unpolarized particles: B01010 01010 01010 01010 = D338250, FC8Guid.pdf p.49
        helicities = 338250
        ! flags: Bit0 (reset) = 1, Bit1 (loop) = 0 -> B01 = D1
        flags = 1
        
        !call clearcache
        fc_result(1) = 0.D0    ! born
        fc_result(2) = 0.D0    ! virtual (not needed here)
        
        ! symmetry factor
        if(bflav(3).eq.bflav(4)) then
          ! 72 = 2*2*3*3*2 (spin factor 2, colour factor 3, equal particles in the final state 2)
          symfac = 72d0 ! take spin and colour factors into account
        else
          ! 36 = 2*2*3*3 (spin factor 2, colour factor 3)
          symfac = 36d0 ! take spin and colour factors into account
        endif
        
#ifndef NEGLECTBMASS
          print*,"ERROR: FormCalc massive b-quark amplitudes not yet "
                 //"implemented."
          stop
#endif

        if(     abs(bflav(1)).eq.2 .and. abs(bflav(2)).eq.2 .and. 
     &          (bflav(3).eq.1000024 .and. bflav(4).eq.-1000024)  ) then ! u + ubar -> x1 + x1
          call uubar_x1x1_squaredME(fc_result, helicities, flags)
          goto 20
        else if(abs(bflav(1)).eq.1 .and. abs(bflav(2)).eq.1 .and. 
     &          (bflav(3).eq.1000024 .and. bflav(4).eq.-1000024)  ) then ! d + dbar -> x1 + x1
          call ddbar_x1x1_squaredME(fc_result, helicities, flags)
          goto 20
        else if(abs(bflav(1)).eq.4 .and. abs(bflav(2)).eq.4 .and. 
     &          (bflav(3).eq.1000024 .and. bflav(4).eq.-1000024)  ) then ! c + cbar -> x1 + x1
          call uubar_x1x1_squaredME(fc_result, helicities, flags)
          goto 20
        else if(abs(bflav(1)).eq.3 .and. abs(bflav(2)).eq.3 .and. 
     &          (bflav(3).eq.1000024 .and. bflav(4).eq.-1000024)  ) then ! s + sbar -> x1 + x1
          call ddbar_x1x1_squaredME(fc_result, helicities, flags)
          goto 20
        else if(abs(bflav(1)).eq.5 .and. abs(bflav(2)).eq.5 .and. 
     &          (bflav(3).eq.1000024 .and. bflav(4).eq.-1000024)  ) then ! b + bbar -> x1 + x1
          call ddbar_x1x1_squaredME(fc_result, helicities, flags)
          goto 20
        else if(abs(bflav(1)).eq.2 .and. abs(bflav(2)).eq.2 .and. 
     &          (bflav(3).eq.1000024 .and. bflav(4).eq.-1000037)  ) then ! u + ubar -> x1 + x2
          call uubar_x1x2_squaredME(fc_result, helicities, flags)
          goto 20
        else if(abs(bflav(1)).eq.1 .and. abs(bflav(2)).eq.1 .and. 
     &          (bflav(3).eq.1000024 .and. bflav(4).eq.-1000037)  ) then ! d + dbar -> x1 + x2
          call ddbar_x1x2_squaredME(fc_result, helicities, flags)
          goto 20
        else if(abs(bflav(1)).eq.4 .and. abs(bflav(2)).eq.4 .and. 
     &          (bflav(3).eq.1000024 .and. bflav(4).eq.-1000037)  ) then ! c + cbar -> x1 + x2
          call uubar_x1x2_squaredME(fc_result, helicities, flags)
          goto 20
        else if(abs(bflav(1)).eq.3 .and. abs(bflav(2)).eq.3 .and. 
     &          (bflav(3).eq.1000024 .and. bflav(4).eq.-1000037)  ) then ! s + sbar -> x1 + x2
          call ddbar_x1x2_squaredME(fc_result, helicities, flags)
          goto 20
        else if(abs(bflav(1)).eq.5 .and. abs(bflav(2)).eq.5 .and. 
     &          (bflav(3).eq.1000024 .and. bflav(4).eq.-1000037)  ) then ! b + bbar -> x1 + x2
          call ddbar_x1x2_squaredME(fc_result, helicities, flags)
          goto 20
        else if(abs(bflav(1)).eq.2 .and. abs(bflav(2)).eq.2 .and. 
     &          (bflav(3).eq.1000037 .and. bflav(4).eq.-1000037)  ) then ! u + ubar -> x2 + x2
          call uubar_x2x2_squaredME(fc_result, helicities, flags)
          goto 20
        else if(abs(bflav(1)).eq.1 .and. abs(bflav(2)).eq.1 .and. 
     &          (bflav(3).eq.1000037 .and. bflav(4).eq.-1000037)  ) then ! d + dbar -> x2 + x2
          call ddbar_x2x2_squaredME(fc_result, helicities, flags)
          goto 20
        else if(abs(bflav(1)).eq.4 .and. abs(bflav(2)).eq.4 .and. 
     &          (bflav(3).eq.1000037 .and. bflav(4).eq.-1000037)  ) then ! c + cbar -> x2 + x2
          call uubar_x2x2_squaredME(fc_result, helicities, flags)
          goto 20
        else if(abs(bflav(1)).eq.3 .and. abs(bflav(2)).eq.3 .and. 
     &          (bflav(3).eq.1000037 .and. bflav(4).eq.-1000037)  ) then ! s + sbar -> x2 + x2
          call ddbar_x2x2_squaredME(fc_result, helicities, flags)
          goto 20
        else if(abs(bflav(1)).eq.5 .and. abs(bflav(2)).eq.5 .and. 
     &          (bflav(3).eq.1000037 .and. bflav(4).eq.-1000037)  ) then ! b + bbar -> x2 + x2
          call ddbar_x2x2_squaredME(fc_result, helicities, flags)
          goto 20
        else
          print*, "ERROR: Wrong flavor structure in subroutine setborn."
          print*, "bflav = ", bflav
          stop
        endif

 20   continue

        born = fc_result(1)/symfac ! take spin and colour factors into account
        fborn = born

        ! spin and color correlations
        do j=1,nlegborn
          if(abs(bflav(j)).le.6) then
          ! gluon spin correlations as in (eq. 2.8) of 1002.2581v1
          ! (A general framework for implementing NLO calculations in [...] the POWHEG BOX)
          ! only different from zero for processes involving gluons in 
          ! the initial or final state.: Bmunu=0
            if(bflav(j).eq.0) then
              do nmu=0,3
                do nnu=0,3
                  bmunu(nmu,nnu,j) = 0D0
                enddo
              enddo
            endif  
          ! color correlations as in (eq. 2.6) of 1002.2581v1
          ! (A general framework for implementing NLO calculations in [...] the POWHEG BOX)
          ! Bij = Cf B for i != j incoming q/qbar
            do k=j+1,nlegborn
              if(abs(bflav(k)).le.6) then
                colcf = cf !4/3
              else
                colcf = 0D0
              endif
                bornjk(j,k) = born*colcf
                bornjk(k,j) = bornjk(j,k)
            enddo
          endif
        enddo
        
        if (verbose) then
          print*, "process: ", bflav
          print*, "s  = ", s
          print*, "t  = ", t
          print*, "u  = ", u
#ifdef DEBUGQ
          print*, "MN = ", par_Fin1mass
          print*, "MN = ", par_Fin2mass
          print*, "MZ = ", par_Zmass
#endif
          print*, "born (formcalc) = ", fborn
#ifdef DEBUGQ 
          do j=1,nlegborn
            print*, "bornjk = ", bornjk(j,:)
          enddo
#endif
          print*
        endif
        
        ! restore flavor list
        bflav = flavc
        p = pc
        
c ---------------------------- End FormCalc Born -----------------------
c endif preprocessor condition "if defined(CHECK_FORM_MAD) || defined(FORM_BORN)"
#endif


#if defined(CHECK_FORM_MAD) || defined(MAD_BORN)
! ---------------------------- MadGraph Born ---------------------------
        
        ! store flavor list, so we can operate on it without
        ! affecting the other program
        flavc = bflav
        pc = p
        
        ! calculate the couplings for this PS point
        call set_ebe_couplings
        call sborn_proc(p,bflav,born,bornjk,bmunu)
        mborn = born
        
        if (verbose) then
          print*, "process: ", bflav
          print*, "s  = ", s
          print*, "t  = ", t
          print*, "u  = ", u
#ifdef DEBUGQ
          print*, "MN = ", par_Fin1mass
          print*, "MN = ", par_Fin2mass
          print*, "MZ = ", par_Zmass
#endif
          print*, "born (madgraph) = ", mborn
#ifdef DEBUGQ
          do j=1,nlegborn
            print*, "bornjk = ", bornjk(j,:)
          enddo
#endif
          print*
        endif
        
        ! restore flavor list
        bflav = flavc
        p = pc

c ---------------------------- End Madgraph Born -----------------------
c endif preprocessor condition "if defined(CHECK_FORM_MAD) || defined(MAD_BORN)"
#endif


#ifdef CHECK_FORM_MAD
        ! compare the formcalc and madgraph amplitudes
        if( (abs((fborn - mborn)/(fborn + mborn)) > 1d-09) ) then ! the relative error is probably larger for smaller amplitudes
          print*, "Error: formcalc and madgraph borns are too different"
          print*, "in process: ", bflav
          print*, "born (formcalc) = ", fborn
          print*, "born (madgraph) = ", mborn
          print*, "rel. err. = ", abs((fborn - mborn)/(fborn + mborn))
          print*
          print*, "Check input parameters or spin and colour factors in"
     &          //" formcalc amplitudes."
          stop
        endif
#endif
 

! If the program calculates both madgraph and formcalc borns (preprocessor flag CHECK_FORM_MAD)
! the user can choose which born amplitude powheg should use (formcalc or madgraph)
#if defined(CHECK_FORM_MAD) && defined(MAD_BORN)
        born = mborn ! (default)
#endif
#if defined(CHECK_FORM_MAD) && defined(FORM_BORN)
        born = fborn
#endif

        if(isnan(born)) then
          print*,"got strange value for born amplitude...",born
          stop
        endif

#ifdef DEBUG1
        if(counter < 5) then
          print*, "process: ", bflav
          print*, "s  = ", s
          print*, "t  = ", t
          print*, "u  = ", u
          print*,"born = ", born
          counter = counter + 1
        endif
#endif

      end

c############### end setborn subroutine ################################

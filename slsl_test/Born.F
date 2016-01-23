c############### Born.f ################################################
c last modified by MK, 22.01.2015
c -6  -5  -4  -3  -2  -1  0  1  2  3  4  5  6
c t~  b~  c~  s~  u~  d~  g  d  u  s  c  b  t

c############### setborn subroutine ####################################
c sets Born cross section
c Born contribution:
c parton parton -> slepton slepton

      subroutine setborn(p,bflav,born,bornjk,bmunu)
        implicit none

#include "pwhg_math.h"
#include "nlegborn.h"
#include "PhysPars.h"

        integer bflav(nlegborn)
        double precision p(0:3,nlegborn), bornjk(nlegborn,nlegborn)
        double precision bmunu(0:3,0:3,nlegborn), born, colcf
        ! formcalc results
        double precision fc_result(2)
        ! helicity and flags parameter for formcalc       
        integer*8 helicities
        integer flags
        ! indices
        integer j,k,nmu,nnu
        ! kinematics: mandelstam variables, quark masses, slepton masses
        double precision s, t, u, MassQ1, MassQ2, MassS1, MassS2
        ! phase space
        double precision ex, ey, ez, psNorm
        ! function to calculate the squared momentum sum
        double precision momsq, momsum3sq, momsum2sq
        ! calculate the relative error between formcalc and madgraph borns
        double precision fborn, mborn
        ! check 4-momentum conservation
        logical lresult
        ! initialization variables 
        logical verbose
        logical init_born
        data init_born/.true./

        if (init_born) then
           print*
           print*, "=== Init born... ==="
           verbose = .true.
           init_born = .false.
        else
           verbose = .false.
        endif
        
#ifdef DEBUGQ
        verbose = .true.
        bflav(1) =  5
        bflav(2) = -5
        bflav(3) =  2000011
        bflav(4) = -2000011
        p(0,1) = 174.87247233432291d0        
        p(1,1) =  0d0        
        p(2,1) =  0d0        
        p(3,1) =  174.87247233432291d0     
        p(0,2) =  174.87247233432291d0        
        p(1,2) =  0d0        
        p(2,2) =  0d0       
        p(3,2) = -174.87247233432291d0     
        p(0,3) =  174.87247233432291d0        
        p(1,3) =  72.645879273957476d0        
        p(2,3) =  0d0       
        p(3,3) = -67.359788636824902d0     
        p(0,4) =  174.87247233432291d0      
        p(1,4) = -72.645879273957476d0       
        p(2,4) =  0d0        
        p(3,4) =  67.359788636824902d0
#endif
         
        if ( bflav(2).ne.-bflav(1) ) then
           print*,"the current Born subprocess is ",bflav
           print*,"check your settings in init_processes"
           stop
        endif
        
        ! check four momentum conservation
        call check_4conservation(p,nlegborn,1,lresult)

#ifdef DEBUGQ
        verbose = .true.
#endif

        if(verbose) then
          ! calculate kinematics only for debug purpose
          s = momsum2sq(p(0:3,1), p(0:3,2))
          t = momsum2sq(p(0:3,1),-p(0:3,3))
          u = momsum2sq(p(0:3,2),-p(0:3,3))
        endif
               
#if !defined(FORM_BORN) && !defined(MAD_BORN)
        print*,"Error: Preprocessor flag FORM_BORN or"
     &       //" MAD_BORN not set in Makefile."
        stop
#endif
        
#if defined(CHECK_FORM_MAD) || defined(FORM_BORN)
        ! -------------------- Formcalc Born ---------------------------
        ! calculate the couplings for this PS point
        call set_ebe_couplings
        ! load process dependent masses and reload Formcalc parameters
        ! (this resets the Formcalc masses to their original values)
        call set_process(bflav(1:4), MassQ1, MassQ2, MassS1, MassS2)
        
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
        call VecSet(3, MassS1, psNorm, ex,ey,ez)
        psNorm=dsqrt(p(1,4)**2+p(2,4)**2+p(3,4)**2)
        ex=p(1,4)/psNorm
        ey=p(2,4)/psNorm
        ez=p(3,4)/psNorm
        call VecSet(4, MassS2, psNorm, ex,ey,ez)

        ! unpolarized particles: B01010 01010 00100 00100 = D338250, FC8Guid.pdf p.49
        helicities = 338052
        ! flags: Bit0 (reset) = 1, Bit1 (loop) = 0 -> B01 = D1
        flags = 1
        
        !call clearcache
        fc_result(1) = 0.D0    ! born
        fc_result(2) = 0.D0    ! virtual (not needed here)
        
        if(abs(bflav(1)).eq.2) then ! u + ubar -> slept + slept
          call uubar_squaredME(fc_result, helicities, flags)
        else if(abs(bflav(1)).eq.1) then ! d + dbar -> slept + slept
          call ddbar_squaredME(fc_result, helicities, flags)
        else if(abs(bflav(1)).eq.4) then ! c + cbar -> slept + slept
          call uubar_squaredME(fc_result, helicities, flags)
        else if(abs(bflav(1)).eq.3) then ! s + sbar -> slept + slept
          call ddbar_squaredME(fc_result, helicities, flags)
        else if(abs(bflav(1)).eq.5) then ! b + bbar -> slept + slept
          call ddbar_squaredME(fc_result, helicities, flags)
        else
          print*, "ERROR: Wrong flavor structure in subroutine setborn."
          print*, "bflav = ", bflav
          stop
        endif

        born = fc_result(1)/36d0 ! take spin and colour factors into account
        fborn = born

        ! spin and color correlations
        do j=1,nlegborn
        
          if(abs(bflav(j)).le.6) then
          ! gluon spin correlations as in (2.8) of 1002.2581v1
          ! (A general framework for implementing NLO calculations in [...] the POWHEG BOX)
          ! only different from zero for processes involving gluons in 
          ! the intial or final state.: Bmunu=0
            if(bflav(j).eq.0) then
              do nmu=0,3
                do nnu=0,3
                  bmunu(nmu,nnu,j) = 0D0
                enddo
              enddo
            endif  
          ! color correlations as in (2.6) of 1002.2581v1
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
          print*, "born (formcalc) = ", fborn
          print*
        endif

        ! -------------------- End Formcalc Born -----------------------
c endif preprocessor condition "if defined(CHECK_FORM_MAD) || defined(FORM_BORN)"
#endif


#if defined(CHECK_FORM_MAD) || defined(MAD_BORN)
        ! -------------------- Madgraph Born ---------------------------
        ! calculate the couplings for this PS point
        call set_ebe_couplings
        call sborn_proc(p,bflav,born,bornjk,bmunu)
        mborn = born

        if (verbose) then
          print*, "process: ", bflav
          print*, "s  = ", s
          print*, "t  = ", t
          print*, "u  = ", u 
          print*, "born (madgraph) = ", mborn
          print*
        endif

        ! -------------------- End Madgraph Born -----------------------
c endif preprocessor condition "if defined(CHECK_FORM_MAD) || defined(MAD_BORN)"
#endif


#ifdef CHECK_FORM_MAD
        ! compare the formcalc and madgraph amplitudes
        ! Note: This only makes sense if the flags neglectlightM and
        !       neglectBMass are set.
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
        
      end

c############### end setborn subroutine ################################

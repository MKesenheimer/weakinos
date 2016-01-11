c############### Virtuals.f ############################################
c last modified by MK, 10.02.2015
c Slepton pair production
c virtual corrections at NLO SQCD:
c parton parton -> slepton slepton

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

        double precision p(0:3,nlegborn)
        integer vflav(nlegborn)
        double precision born, virtual, ampvirt
        ! formcalc results
        double precision fc_result(2)
        ! helicity and flags parameter for formcalc       
        integer*8 helicities
        integer flags
        ! phase space
        double precision ex, ey, ez, psNorm
        ! function to calculate the squared momentum sum
        double precision momsum2sq
        ! kinematics
        double precision s, t, u
        double precision MassQ1, MassQ2, MassS1, MassS2
        ! Loops
        integer i,j,k
        ! status
        logical verbose
        logical init_virtual
        data init_virtual/.true./

        if (init_virtual) then
           print*
           print*, "=== Init virtual... ==="
           verbose = .true.
           init_virtual = .false.
        else
           verbose = .false.
        endif
        
        
#ifdef DEBUGQ
        verbose = .true.
#endif
#ifdef DEBUGQ
#define DEBUG1
        verbose = .true.
        st_alpha =  9.6934150432177207D-2
        vflav(1) =  1
        vflav(2) = -1
        vflav(3) =  2000011
        vflav(4) = -2000011
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
        
        if (vflav(1).ne.(-vflav(2))) then
          print*, "Error in fct. compvirtual: invalid initial state"
          print*, "vflav(1) = ", vflav(1)
          print*, "vflav(2) = ", vflav(2)
          stop
        endif
        
        ! check four momentum conservation
        call check_4conservation(p, nlegborn)
        call set_ebe_couplings
        ! load process dependent masses and reload Formcalc parameters
        ! (this resets the Formcalc masses to their original values)
        call set_process(vflav(1:4), MassQ1, MassQ2, MassS1, MassS2)
        
        ! calculate kinematics
        s = momsum2sq(p(0:3,1), p(0:3,2))
        t = momsum2sq(p(0:3,1),-p(0:3,3))
        u = momsum2sq(p(0:3,2),-p(0:3,3))     
       
        ! unpolarized particles: B01010 01010 00100 00100 = D338052, FC8Guid.pdf p.49
        helicities = 338052
        ! flags: Bit0 (reset) = 1, Bit1 (loop) = 1 -> B11 = D3
        flags = 3
        
        ! used by Formcalc
#ifdef DEBUGQ
        call setdebugkey(2)  ! DEBUG for LoopTools
#endif
        !call setversionkey(4)
        call setlambda(0d0) ! Only finite part
        call setuvdiv(0d0)  ! no 1/eps_UV terms 
        call setdelta(0d0)  ! set divergences to zero
        call setmudim(st_muren2) ! renormalization scale
        !print*,"DEBUG: minMass = ", getminmass()
        !call setminmass(10d0)
        
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

c#define CHECKUV
#ifdef CHECKUV
        ! test for finiteness
        ! call setMuDim(1D100)
        call setuvdiv(1d0)
        call setDelta(1D7)    ! Test UV finiteness
#endif

        call calcRenConst
        fc_result(1) = 0.D0    ! born
        fc_result(2) = 0.D0    ! virtual amplitude

        if(abs(vflav(1)).eq.2) then ! u + ubar -> slept + slept
          call uubar_squaredME(fc_result, helicities, flags)
        else if(abs(vflav(1)).eq.1) then ! d + dbar -> slept + slept
          call ddbar_squaredME(fc_result, helicities, flags)
        else if(abs(vflav(1)).eq.4) then ! c + cbar -> slept + slept
          call uubar_squaredME(fc_result, helicities, flags)
        else if(abs(vflav(1)).eq.3) then ! s + sbar -> slept + slept
          call ddbar_squaredME(fc_result, helicities, flags)
        else if(abs(vflav(1)).eq.5) then ! b + bbar -> slept + slept
          call ddbar_squaredME(fc_result, helicities, flags)
        else
          print*, "ERROR: Wrong flavor structure in subr. setvirtual."
          print*, "vflav(1) = ", vflav(1)
          stop
        endif

        born    = fc_result(1)/36d0
        virtual = fc_result(2)/36d0 ! = 2 x Re[M_tree* x M_1-loop], FC8Guid.pdf p.48
        
        ampvirt = virtual/(st_alpha/(2d0*pi))
        !ampvirt = virtual*2d0*pi ! FormCalc Alfas = 1


        if(isnan(virtual)) then
          print*,"got strange value for virtual amplitude...", virtual
          stop
        endif
        
        if (verbose) then
          print*, "process: ", vflav
          print*, "s     = ", s
          print*, "t     = ", t
          print*, "u     = ", u
          print*, "muren2         = ", st_muren2
          print*, "born           = ", born
          print*, "2*Born*Virtual = ", ampvirt
          print*
        end if
        
#ifdef DEBUG1
        stop
#endif
        
      end

c############### end subroutine setvirtual #############################

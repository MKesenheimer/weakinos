c############### Reals.f ###############################################
c last modified by MK, 18.12.2014
c adapted from dislepton
c neutralino pair production
c real emission contributions at NLO SQCD:
c parton parton -> neutralino neutralino + parton

c############### subroutine setreal ####################################

c neutral current neutralino pairs + jet
c q(1) qbar(2) -> Z -> neutralino(3) neutralino(4) gluon(5)
c average over incoming helicites
      subroutine setreal(p,rflav,ampreal)
        implicit none

#include "PhysPars.h"
#include "nlegborn.h"
#include "pwhg_flst.h"
#include "pwhg_math.h"
#include "pwhg_st.h"
#include "Flags.h"

        double precision p(0:3,nlegreal),p_in(0:3,nlegreal)
        integer rflav(nlegreal), dflav(nlegreal)
        ! the real amplitudes
        double precision ampreal
        double precision amp
        ! function to calculate the squared momentum sum
        double precision momsum3sq, momsum2sq
        ! kinematics
        double precision s125, s12, s25, s15
        ! initialization variables
        logical verbose
        logical init_real
        data init_real/.true./
#ifdef DEBUG
#define DEBUG1
        ! debugging
        integer counter
        data counter/0/
        save counter
#endif

        if (init_real) then
           print*
           print*, "=== Init real... ==="
           verbose = .true.
           init_real = .false.
        else
           verbose = .false.
        endif
        
        ! store flavor list and momenta
        dflav = rflav
        p_in(:,:) = p(:,:)
        
        if (.not.( (abs(rflav(2)).eq.abs(rflav(1)))
     &      .or.   (abs(rflav(2)+rflav(1)).eq.1)
     &      .or.   (rflav(2).ne.0)
     &      .or.   (rflav(1).ne.0)   )) then
           print*,"the current Real subprocess is ",rflav
           print*,"check your settings in init_processes"
           stop
        endif
        
        ! check four momentum conservation
        call check_4conservation(p, nlegborn+1)
        
        if(verbose) then
          ! calculate kinematics only for debugging purpose
          s125 = momsum3sq(p(0:3,1),p(0:3,2),p(0:3,5)) ! invariant mass of particle 1 + 2 + 5
          s12  = momsum2sq(p(0:3,1),p(0:3,2)) ! invariant mass of particle 1 + 2
          s25  = momsum2sq(p(0:3,2),p(0:3,5)) ! invariant mass ofparticle 2 + 5
          s15  = momsum2sq(p(0:3,1),p(0:3,5)) ! invariant mass of particle 1 + 5
        endif
        
        call set_ebe_couplings
        call sreal_proc(p,rflav,amp)

        ! Cancel as/(2pi) associated with amp2. It will be put back by real_ampsq
        ampreal = amp/(st_alpha/(2d0*pi))

        if (verbose) then
          print*, "process = ", rflav
          print*, "s125 = ", s125
          print*, "s12  = ", s12
          print*, "s25  = ", s25 
          print*, "s15  = ", s15
          print*, "Real amplitude = ", ampreal
          print*
        end if
        
        ! restore flavor list
        rflav = dflav
        p(:,:) = p_in(:,:)

#ifdef DEBUG1
        if(counter < 5) then
          print*,"real = ", ampreal
          counter = counter + 1
        endif
#endif
        
      end

c############### end subroutine setreal ################################
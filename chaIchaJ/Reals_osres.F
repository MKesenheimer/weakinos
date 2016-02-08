c############### Real_osres.f ##########################################
c last modified by MK, 25.11.2015
c weakino pair production
c -6  -5  -4  -3  -2  -1  0  1  2  3  4  5  6
c t~  b~  c~  s~  u~  d~  g  d  u  s  c  b  t

c############### subroutine setosresreal ###############################
c this is the routine which gives back |M_R|**2 - BW*corrfac*|M_R^OS|**2
c (the real contributions with a resonant intermediate squark subtracted 
c by on-shell counter terms) evaluated at a phase space point where 
c sij = mij**2 holds.
c the corrfac factor accounts for the restriction of the pahse-space and
c modifies the jacobian.
c
c In principle, we can have a resonance between particle i=3 and j=5
c or i=4 and j=5. This can be distinguished by the flag "chan":
c osresID 1 = dl35: left handed up-type squark gets resonant with particles 3 & 5
c osresID 2 = dl45: left handed up-type squark gets resonant with particles 4 & 5
c osresID 3 = ul35: left handed down-type squark gets resonant with particles 3 & 5
c osresID 4 = ul45: left handed down-type squark gets resonant with particles 4 & 5
c (only required for neutralino pair production:)
c osresID 5 = dr35: right handed up-type squark gets resonant with particles 3 & 5
c osresID 6 = dr45: right handed up-type squark gets resonant with particles 4 & 5
c osresID 7 = ur35: right handed down-type squark gets resonant with particles 3 & 5
c osresID 8 = ur45: right handed down-type squark gets resonant with particles 4 & 5
      subroutine  setosresreal(p,rflav,chan,ampreal)
        implicit none

#include "PhysPars.h"
#include "nlegborn.h"
#include "pwhg_flst.h"
#include "pwhg_math.h"
#include "pwhg_st.h"
#include "Flags.h"
#include "osres.h"

        double precision p(0:3,nlegreal),p_in(0:3,nlegreal)
        double precision p_os(0:3,nlegreal)
        integer rflav(nlegreal), dflav(nlegreal)
        ! the real amplitudes
        double precision ampreal
        ! variable to determine the channel
        character*4 chan
        ! local variables
        double precision ampCT, amp
        ! function to calculate the squared momentum sum
        double precision momsum3sq, momsum2sq
        external momsum2sq, momsum3sq
        ! jacobian for the on-shell terms
        double precision corrfac
        external corrfac
        ! kinematics
        double precision s12, sij
        ! breit-wigner ratios
        double precision bwratio
        ! masses
        double precision MassNi, MassNj, MassXi, MassXj
        double precision mij, mi, mj, mk
        ! indices
        integer i,j,k
        ! check 4-momentum conservation
        logical lresult
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
        
#if defined(DSUB_II) || defined(DSUB_II_TEST)

        if (init_real) then
           print*
           print*, "=== Init resonant real... ==="
           print*, "Computing the integral of the "//
     &             "remnant cross section consisting of "
           print*, "remaining terms from on-shell resonant "//
     &             "reals, incl. on-shell-subtraction "
           print*, "to set up the adaptive grid"
           verbose = .true.
           init_real = .false.
        else
           verbose = .false.
        endif
        
        ! store flavor list and momenta
        dflav = rflav
        p_in(:,:) = p(:,:)
        
        ! tests
        if (.not.( (abs(rflav(2)).eq.abs(rflav(1)))
     &      .or.   (abs(rflav(2)+rflav(1)).eq.1)
     &      .or.   (rflav(2).ne.0)
     &      .or.   (rflav(1).ne.0)   )) then
           print*,"the current Real subprocess is ",rflav
           print*,"check your settings in init_processes"
           stop
        endif
        
        if( chan.ne."ul35" .and. chan.ne."ul45" .and. 
     &      chan.ne."ur35" .and. chan.ne."ur45" .and.
     &      chan.ne."dl35" .and. chan.ne."dl45" .and. 
     &      chan.ne."dr35" .and. chan.ne."dr45" .and.
     &      chan.ne."allr ") then
          print*,"wrong channel number: ", chan
          print*,"choose ul35, ul45, ur35, ur45, dl35, dl45, dr35, "//
     &           "dr45 or allr."     
          stop
        endif
        
#ifndef DSUB_II_TEST
        if(flg_btilde) then
          print*, "the realresonant-routine should not cald for Btilde."
          stop
        endif
#endif

        ! check four momentum conservation
        call check_4conservation(p, nlegborn+1, 1, lresult)
        
        ! reset the breit-wigner ratios
        bwratio = 1d0
        amp   = 0d0
        ampCT = 0d0
        
        ! calculate the couplings for this PS point
        call set_ebe_couplings
        
        ! set the channel-related indices i,j,k and masses mi,mj,mk,mij
        call set_channel(chan,i,j,k,mi,mj,mk,mij) ! OK, but better check again!

        ! calculate the invariants
        s12  = momsum2sq(p(0:3,1),p(0:3,2)) ! invariant mass of particle 1 + 2
        sij  = momsum2sq(p(0:3,i),p(0:3,j)) ! invariant mass of particle i + j

        ! this is |M_R(res)|**2 (still with on-shell resonances)
        call sreal_proc_res(p,rflav,chan,amp)

        ! calculate the on-shell counter terms BW*corrfac*|M_R^OS(res)|**2
        ! (don't subtract anything if the squark is off-shell, but if 
        ! the quark is on-shell subtract the divergent part) 
        ! -> Prospino scheme 1211.0286
        if( (s12.ge.(mij+mk)**2) .and. (mij.ge.(mi+mj)) ) then
          call off_to_on(p,chan,p_os)                             ! off_to_on the momenta p to on-shell momenta p_os        
          bwratio = (mij*wreg)**2/((sij-mij**2)**2+(mij*wreg)**2) ! calculate the ratio of the breit wigner functions
          call sreal_proc_res(p_os,rflav,chan,ampCT)              ! new routine to calculate only the resonant matrix elements
          ampCT = bwratio*ampCT                                   ! multiply with the breit-wigner ratio
          ampCT = corrfac(s12,mi,mj,mk,sij,mij)*ampCT             ! rescale the counter term to fit into the new ps
        else
          ampCT = 0D0
        endif  

        ! calculate the total resonant amplitude and substract the
        ! on-shell counter terms -> no on-shell resonances anymore
        ampreal = (+amp-ampCT)
 
        if (verbose) then
          print*, "process = ", rflav
          print*, "chan = ", chan
          print*, "s12  = ", s12
          print*, "on-shell resonant Real amplitude = ", ampreal
          print*
        end if

        ! restore flavor list
        rflav = dflav
        p(:,:) = p_in(:,:)
        
#ifdef DEBUG1
        if(counter .lt. 5) then
          print*,"real_osres = ", ampreal
          counter = counter + 1
        endif
#endif
        
#else
        ampreal = 0d0
#endif

      end
      
c############### end subroutine setosresreal ###########################
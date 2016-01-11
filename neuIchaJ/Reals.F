c############### Reals.f ###############################################
c last modified by MK, 02.12.2015
c adapted from dislepton
c weakino pair production
c real emission contributions at NLO SQCD:
c parton parton -> weakino + weakino + parton

c############### subroutine setreal ####################################
c neutral current weakino pairs + jet
c q(1) qbar(2) -> Z/W -> weakino(3) weakino(4) gluon(5)
c average over incoming helicites
c
c NOTE: with flag DSUB_II defined, this routine only calculates the 
c real contributions without the resonant matrix element squared:
c |M_NR|**2 + 2 M_R* M_NR.
c The remaining terms and the on-shell counter terms are added back
c in subroutine setosresreal:
c |M_R|**2 - BW*corrfac*|M_R^OS|**2.

      subroutine setreal(p,rflav,ampreal)
        implicit none

#include "PhysPars.h"
#include "nlegborn.h"
#include "pwhg_flst.h"
#include "pwhg_math.h"
#include "pwhg_st.h"
#include "Flags.h"
#include "osres.h"

        double precision p(0:3,nlegreal),p_in(0:3,nlegreal)
        integer rflav(nlegreal), dflav(nlegreal)
        ! the real amplitudes
        double precision ampreal
#ifdef DSUB_II_TEST
        ! the real resonant amplitudes (testing)
        double precision amposres
        integer i
        double precision amprealres_ul35, amprealres_ul45
        double precision amprealres_dl35, amprealres_dl45
#endif
        double precision amp
        ! function to calculate the squared momentum sum
        double precision momsum3sq, momsum2sq
        external momsum2sq, momsum3sq
        ! masses
        double precision MassNi, MassXj
        ! kinematics
        double precision s125, s12, s25, s15, s35, s45
        ! cut on the onvariant mass sij
        logical doradcut
        double precision theta_dl35, theta_ul35, theta_dl45, theta_ul45
        double precision delta_cut
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
        
        ! calculate the invariants
        s125 = momsum3sq(p(0:3,1),p(0:3,2),p(0:3,5)) ! invariant mass of particle 1 + 2 + 5
        s12  = momsum2sq(p(0:3,1),p(0:3,2)) ! invariant mass of particle 1 + 2
        s25  = momsum2sq(p(0:3,2),p(0:3,5)) ! invariant mass ofparticle 2 + 5
        s15  = momsum2sq(p(0:3,1),p(0:3,5)) ! invariant mass of particle 1 + 5
        s35  = momsum2sq(p(0:3,3),p(0:3,5))
        s45  = momsum2sq(p(0:3,4),p(0:3,5))
        
        ! calculate the couplings for this PS point
        call set_ebe_couplings
        
        ! get the real amplitude from MadGraph
        call sreal_proc(p,rflav,amp)

        ! Cancel as/(2pi) associated with amp2. It will be put back by real_ampsq
        ampreal = amp/(st_alpha/(2d0*pi))
        
#ifdef DSUB_II_TEST
        ! test the DSUB_II scheme (here without the correct jacobian)
        ! this gives the same result as with flag DSUB_I set
        ! if the real amplitude is zero, do not subtract anything
        !if(.not.(ampreal.eq.0d0)) then
          do i=1,nosres
            amposres = 0D0
            call setosresreal(p,rflav,osresID(i),amposres)
            ! Note: amposres has the factor as/(2pi)
            ! included, so cancel as/(2pi) here
            ampreal = ampreal + amposres/(st_alpha/(2d0*pi))
          enddo
        !endif
#endif

        ! ==- cut on invariant mass when in radiation-phase -==
        ! We investigated that there are many negative events generated
        ! close to the resonant region. We introduced a cut on the 
        ! invariant mass Sqrt(sij) of the intermediate squark to avoid
        ! negative amplitudes squared. We don't have to do this
        ! in the DR-I case.
        doradcut = flg_rad ! cut only in radiation-phase
#ifdef DR_I
        doradcut = .false.
#endif

        if(doradcut) then
          ! reset the temp. variables
          theta_dl35 = 0D0
          theta_ul35 = 0D0
          theta_dl45 = 0D0
          theta_ul45 = 0D0
          delta_cut = 1D0
          
          ! load process dependent masses      
          MassNi = par_Fin1mass
          MassXj = par_Fin2mass
          
          call increasecnt('total calls to real-routine in rad-gen')

          ! we need to cut only in case of on-shell-resonances!
          ! u-g-channel
          if((((abs(rflav(1)).eq.2).and.(abs(rflav(2)).eq.0)) .or. 
     &       ((abs(rflav(2)).eq.2).and.(abs(rflav(1)).eq.0)) .or. 
     &       ((abs(rflav(1)).eq.4).and.(abs(rflav(2)).eq.0)) .or. 
     &       ((abs(rflav(2)).eq.4).and.(abs(rflav(1)).eq.0))) ) then
            ! s35 = mdl^2, m3 = mni, m5 = 0D0
            if((s12.ge.(mdl+MassXj)**2).and.(mdl.ge.dabs(MassNi))) then
              ! calculate the argument of the theta function
              theta_dl35 = dabs(mdl-dsqrt(s35))-delta_cut
            endif
            ! s45 = mul^2, m4 = MassXj, m5 = 0d0
            if((s12.ge.(mul+dabs(MassNi))**2).and.(mul.ge.MassXj)) then
              ! calculate the argument of the theta function
              theta_ul45 = dabs(mul-dsqrt(s45))-delta_cut
            endif
          elseif((((abs(rflav(1)).eq.1).and.(abs(rflav(2)).eq.0)) .or. 
     &       ((abs(rflav(2)).eq.1).and.(abs(rflav(1)).eq.0)) .or. 
     &       ((abs(rflav(1)).eq.3).and.(abs(rflav(2)).eq.0)) .or. 
     &       ((abs(rflav(2)).eq.3).and.(abs(rflav(1)).eq.0))) ) then
            ! s35 = mul^2, m3 = mni, m5 = 0D0
            if((s12.ge.(mul+MassXj)**2).and.(mul.ge.dabs(MassNi))) then
              ! calculate the argument of the theta function
              theta_ul35 = dabs(mul-dsqrt(s35))-delta_cut
            endif
            ! s45 = mdl^2, m4 = MassXj, m5 = 0d0
            if((s12.ge.(mdl+dabs(MassNi))**2).and.(mdl.ge.MassXj)) then
              ! calculate the argument of the theta function
              theta_dl45 = dabs(mdl-dsqrt(s45))-delta_cut
            endif
          endif
        
          ! theta functions
          if((theta_dl35.le.0D0) .and. (ampreal.le.0.d0)) then
            ampreal = 1D-20
            call increasecnt('radiation in on-shell-res. region vetoed')
          endif
          if((theta_ul35.le.0D0) .and. (ampreal.le.0.d0)) then
            ampreal = 1D-20
            call increasecnt('radiation in on-shell-res. region vetoed')
          endif
          if((theta_dl45.le.0D0) .and. (ampreal.le.0.d0)) then
            ampreal = 1D-20
            call increasecnt('radiation in on-shell-res. region vetoed')
          endif
          if((theta_ul45.le.0D0) .and. (ampreal.le.0.d0)) then
            ampreal = 1D-20
            call increasecnt('radiation in on-shell-res. region vetoed')
          endif
        endif
        
        
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
        if(counter .lt. 5) then
          print*,"real = ", ampreal
          counter = counter + 1
        endif
#endif

      end

c############### end subroutine setreal ################################
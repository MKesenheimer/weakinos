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
        
#ifdef DEBUG
#define RESCALE 1D-1
      ! set 1  
      rflav(1) = 1
      rflav(2) = 0
      rflav(3) = 1000022
      rflav(4) = -1000024
      rflav(5) = 2
      p(0,1) =    19.438476758172818D0*RESCALE     
      p(1,1) =    0.0000000000000000D0*RESCALE     
      p(2,1) =    0.0000000000000000D0*RESCALE     
      p(3,1) =    19.438476758172818D0*RESCALE     
      p(0,2) =    6785.2278051874255D0*RESCALE     
      p(1,2) =    0.0000000000000000D0*RESCALE     
      p(2,2) =    0.0000000000000000D0*RESCALE     
      p(3,2) =   -6785.2278051874255D0*RESCALE     
      p(0,3) =    527.72993002991370D0*RESCALE     
      p(1,3) =    127.18745113035945D0*RESCALE     
      p(2,3) =  -0.14061483021967258D0*RESCALE     
      p(3,3) =   -502.96483867384944D0*RESCALE     
      p(0,4) =    1750.9038308799525D0*RESCALE     
      p(1,4) =   -127.29569773903935D0*RESCALE     
      p(2,4) =  -0.18443033468507580D0*RESCALE     
      p(3,4) =   -1736.7919816859594D0*RESCALE     
      p(0,5) =    4526.0325210356950D0*RESCALE     
      p(1,5) =   0.10824660867989543D0*RESCALE     
      p(2,5) =   0.32504516490474800D0*RESCALE     
      p(3,5) =   -4526.0325080694065D0*RESCALE
      !rpowheg =   0.34201897787940250D0
#endif
        
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

#ifdef DEBUGQ
        print*
        print*,"     rflav",rflav
        print*,"     p(0,1) = ",p(0,1)
        print*,"     p(1,1) = ",p(1,1)
        print*,"     p(2,1) = ",p(2,1)
        print*,"     p(3,1) = ",p(3,1)
        print*,"     p(0,2) = ",p(0,2)
        print*,"     p(1,2) = ",p(1,2)
        print*,"     p(2,2) = ",p(2,2)
        print*,"     p(3,2) = ",p(3,2)
        print*,"     p(0,3) = ",p(0,3)
        print*,"     p(1,3) = ",p(1,3)
        print*,"     p(2,3) = ",p(2,3)
        print*,"     p(3,3) = ",p(3,3)
        print*,"     p(0,4) = ",p(0,4)
        print*,"     p(1,4) = ",p(1,4)
        print*,"     p(2,4) = ",p(2,4)
        print*,"     p(3,4) = ",p(3,4)
        print*,"     p(0,5) = ",p(0,5)
        print*,"     p(1,5) = ",p(1,5)
        print*,"     p(2,5) = ",p(2,5)
        print*,"     p(3,5) = ",p(3,5)
        print*,"     rpowheg = ",ampreal
        stop
        !verbose = .true.
#endif
        
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

        stop
        
      end

c############### end subroutine setreal ################################
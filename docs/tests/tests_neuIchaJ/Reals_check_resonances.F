c############### Reals.f ###############################################
c last modified by MK, 18.12.2014
c adapted from dislepton
c weakino pair production
c real emission contributions at NLO SQCD:
c parton parton -> weakino + weakino + parton

c############### subroutine setreal ####################################

c neutral current weakino pairs + jet
c q(1) qbar(2) -> Z -> weakino(3) weakino(4) gluon(5)
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
        double precision ampreal, amprealres135, amprealres145
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
        
#ifdef DEBUGQ
      rflav(1) = 1
      rflav(2) = 0
      rflav(3) = 1000022
      rflav(4) = -1000024
      rflav(5) = 2
#endif
        
#ifdef DEBUGQ
!       ! set 1 (Resonanz in s45 und keine in s35)
!       p(0,1) =    19.438476758172818D0     
!       p(1,1) =    0.0000000000000000D0     
!       p(2,1) =    0.0000000000000000D0     
!       p(3,1) =    19.438476758172818D0     
!       p(0,2) =    6785.2278051874255D0     
!       p(1,2) =    0.0000000000000000D0     
!       p(2,2) =    0.0000000000000000D0     
!       p(3,2) =   -6785.2278051874255D0     
!       p(0,3) =    527.72993002991370D0     
!       p(1,3) =    127.18745113035945D0     
!       p(2,3) =  -0.14061483021967258D0     
!       p(3,3) =   -502.96483867384944D0     
!       p(0,4) =    1750.9038308799525D0     
!       p(1,4) =   -127.29569773903935D0     
!       p(2,4) =  -0.18443033468507580D0     
!       p(3,4) =   -1736.7919816859594D0     
!       p(0,5) =    4526.0325210356950D0     
!       p(1,5) =   0.10824660867989543D0     
!       p(2,5) =   0.32504516490474800D0     
!       p(3,5) =   -4526.0325080694065D0
      
!       ! set 2 (keine Resonanz) 
!       p(0,1) =   5693.7127698464756D0  
!       p(1,1) =   0.0000000000000000D0
!       p(2,1) =   0.0000000000000000D0
!       p(3,1) =   5693.7127698464756D0
!       p(0,2) =   4.2972373564317135D0
!       p(1,2) =   0.0000000000000000D0
!       p(2,2) =   0.0000000000000000D0
!       p(3,2) =  -4.2972373564317135D0
!       p(0,3) =   2505.3760117060488D0
!       p(1,3) =   2.8472481258550761D0
!       p(2,3) =   5.4717115253711510D0
!       p(3,3) =   2503.5020133751846D0
!       p(0,4) =   2809.4949033458938D0
!       p(1,4) =  -21.073561284857632D0
!       p(2,4) =   9.8925281006224033D0
!       p(3,4) =   2803.5167283547771D0
!       p(0,5) =   383.13909215077155D0
!       p(1,5) =   18.226313159002505D0
!       p(2,5) =  -15.364239625993513D0
!       p(3,5) =   382.39679075988914D0

!       ! set 3 (nahe der Resonanz von s35)
!       p(0,1) =   669.30634556759014D0
!       p(1,1) =   0.0000000000000000D0
!       p(2,1) =   0.0000000000000000D0
!       p(3,1) =   669.30634556759014D0
!       p(0,2) =   1420.4614051200483D0
!       p(1,2) =   0.0000000000000000D0
!       p(2,2) =   0.0000000000000000D0
!       p(3,2) =  -1420.4614051200483D0
!       p(0,3) =   152.30166358523340D0
!       p(1,3) =  -12.194819086243115D0
!       p(2,3) =  -53.305913651724374D0
!       p(3,3) =  -104.19683334170966D0
!       p(0,4) =   1061.6243902695826D0
!       p(1,4) =  -827.39814691512720D0
!       p(2,4) =  -193.30143137038658D0
!       p(3,4) =  -609.98327883869388D0
!       p(0,5) =   875.84169683282278D0
!       p(1,5) =   839.59296600137009D0
!       p(2,5) =   246.60734502211093D0
!       p(3,5) =  -36.974947372054885D0

      ! set 4 (nahe der Resonanz von s45)  
      p(0,1) =   206.34511479462293D0
      p(1,1) =   0.0000000000000000D0
      p(2,1) =   0.0000000000000000D0
      p(3,1) =   206.34511479462293D0
      p(0,2) =   1247.0216746997478D0
      p(1,2) =   0.0000000000000000D0
      p(2,2) =   0.0000000000000000D0
      p(3,2) =  -1247.0216746997478D0
      p(0,3) =   537.26597598975775D0
      p(1,3) =   338.15989307456852D0
      p(2,3) =   53.127339285862099D0
      p(3,3) =  -402.65558344921777D0
      p(0,4) =   787.58623384230759D0
      p(1,4) =  -394.55062587246186D0
      p(2,4) =   61.228417262226635D0
      p(3,4) =  -654.10958685612366D0
      p(0,5) =   128.51457966230558D0
      p(1,5) =   56.390732797893314D0
      p(2,5) =  -114.35575654808872D0
      p(3,5) =   16.088610400216368D0
#endif


#ifdef DEBUGQ
        call setresonantreal(p,rflav,amprealres135,135)
        call setresonantreal(p,rflav,amprealres145,145)
        !print*,amprealres145,rflav
        if(.true.!abs(rflav(1)).eq.2.and.abs(rflav(2)).eq.0
     &     .and.amprealres145.gt.1d-5) then
        print*
        print*, "process = ", rflav
        print*,"     p(0,1) =",p(0,1)
        print*,"     p(1,1) =",p(1,1)
        print*,"     p(2,1) =",p(2,1)
        print*,"     p(3,1) =",p(3,1)
        print*,"     p(0,2) =",p(0,2)
        print*,"     p(1,2) =",p(1,2)
        print*,"     p(2,2) =",p(2,2)
        print*,"     p(3,2) =",p(3,2)
        print*,"     p(0,3) =",p(0,3)
        print*,"     p(1,3) =",p(1,3)
        print*,"     p(2,3) =",p(2,3)
        print*,"     p(3,3) =",p(3,3)
        print*,"     p(0,4) =",p(0,4)
        print*,"     p(1,4) =",p(1,4)
        print*,"     p(2,4) =",p(2,4)
        print*,"     p(3,4) =",p(3,4)
        print*,"     p(0,5) =",p(0,5)
        print*,"     p(1,5) =",p(1,5)
        print*,"     p(2,5) =",p(2,5)
        print*,"     p(3,5) =",p(3,5)
        print*,"amprealres_135 = ", amprealres135
        print*,"amprealres_145 = ", amprealres145
        endif
#endif
        
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
        
        ! calculate the couplings for this PS point
        call set_ebe_couplings
        
        ! get the real amplitude from MadGraph
        call sreal_proc(p,rflav,amp,135)

        ! Cancel as/(2pi) associated with amp2. It will be put back by real_ampsq
        ampreal = amp/(st_alpha/(2d0*pi))
        
#ifdef DSUB_II
        call setresonantreal(p,rflav,amprealres135,135)
        call setresonantreal(p,rflav,amprealres145,145)
        ampreal = ampreal - amprealres135 - amprealres145
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
        
#ifdef DEBUGQ
        print*,"amprealres_135 = ", amprealres135
        print*,"amprealres_145 = ", amprealres145
        stop
#endif
        
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
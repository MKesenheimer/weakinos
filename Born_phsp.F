c############### Born.f ################################################
c last modified by MK, 22.06.2015
c weakino pair production

c############### subroutine born_phsp ##################################
c fills the four-momenta of the Born-process particles,
c both in the laboratory and in the center-of-mass frame,
c the Bjorken x of the two incoming partons,
c the partonic Mandelstam variable s,
c the masses for Born and real emission,
c and the Jacobi factor for the phase space volume
      subroutine born_phsp(xborn)
        implicit none

#include "PhysPars.h"
#include "nlegborn.h"
#include "pwhg_kn.h"

        ! input:
        double precision xborn(ndiminteg-3)
        ! output, local variables:
        double precision vec(3),beta
        double precision MassQ1, MassQ2, MassF1, MassF2
        integer nu
        
        kn_jacborn = 1d0

#ifdef DEBUGQ
        print*,"xborn(1) = ", xborn(1)
        print*,"xborn(2) = ", xborn(2)
        print*,"xborn(3) = ", xborn(3)
#define DEBUG1
#endif

        ! set masses for Born and real emission
        kn_masses(1) = 0d0
        kn_masses(2) = 0d0
        kn_masses(3) = dabs(par_Fin1mass)
        kn_masses(4) = dabs(par_Fin2mass)
        kn_masses(5) = 0d0

        ! lower bound on the mass of the final state
        kn_minmass = kn_masses(3) + kn_masses(4)


#ifdef DEBUG1
        print*,"minmass = ", kn_minmass
        print*,"sbeams  = ", kn_sbeams
#endif

        ! calculate the parton momentum fractions kn_xb1, kn_xb2
        ! the squared CM energy of the Born process sborn
        ! and the jacobi factor of the PS volume
        call x1x2phspace(kn_sbeams,kn_minmass,xborn(1),kn_xb1,kn_xb2,
     &                   kn_sborn,kn_jacborn) ! OK!

        ! center of mass momenta
        kn_cmpborn(0,1) =  sqrt(kn_sborn)/2d0
        kn_cmpborn(1,1) =  0d0
        kn_cmpborn(2,1) =  0d0
        kn_cmpborn(3,1) =  kn_cmpborn(0,1)
        kn_cmpborn(0,2) =  kn_cmpborn(0,1)
        kn_cmpborn(1,2) =  0d0
        kn_cmpborn(2,2) =  0d0
        kn_cmpborn(3,2) = -kn_cmpborn(0,1)

        ! generate four-momenta in the cms frame
        ! for the particles produced in the Born process,
        ! calculates the Jacobi factor for the phase space volume
        call twobody0(xborn(3),kn_sborn,kn_masses(3),kn_masses(4),
     &                kn_cmpborn(:,3),kn_cmpborn(:,4),kn_jacborn) ! OK!

        ! velocity of subsystem in laboratory frame
        beta = (kn_xb1-kn_xb2)/(kn_xb1+kn_xb2)

        ! set up boost vector
        vec(1) = 0
        vec(2) = 0
        vec(3) = 1

        ! calculate momenta kn_pborn in the laboratory frame
        ! boost kn_cmpborn (CM) -> kn_pborn (lab. frame)
        call mboost(nlegborn-2,vec,beta,kn_cmpborn(0,3),kn_pborn(0,3))

        do nu=0,3 ! OK!
           ! kn_beams(0,1) =  dsqrt(kn_sbeams)/2
           ! kn_beams(3,1) =  dsqrt(kn_sbeams)/2
           ! kn_beams(0,2) =  dsqrt(kn_sbeams)/2
           ! kn_beams(3,2) =  -dsqrt(kn_sbeams)/2
           kn_pborn(nu,1) = kn_xb1*kn_beams(nu,1)
           kn_pborn(nu,2) = kn_xb2*kn_beams(nu,2)
        enddo

        ! check four momentum conservation
        call check_4conservation(kn_cmpborn, nlegborn)
        ! check four momentum conservation
        call check_4conservation(kn_pborn, nlegborn)

#ifdef DEBUG1
        print*,"bjorken x:"
        print*,"kn_xb1 = ", kn_xb1
        print*,"kn_xb2 = ", kn_xb2
        print*,"energy and momenta in lab frame:"
        print*,"E1, pz1  = ", kn_xb1*dsqrt(kn_sbeams)/2
        print*,"E2, -pz2 = ", kn_xb2*dsqrt(kn_sbeams)/2
        print*,"invariant masses"
        print*,"sborn = ", kn_sborn
        print*,"sborn = ", kn_xb1*kn_xb2*kn_sbeams
        print*,"laboratory frame:"
        print*,"kn_pborn(1) = ", kn_pborn(:,1)
        print*,"kn_pborn(2) = ", kn_pborn(:,2)
        print*,"kn_pborn(3) = ", kn_pborn(:,3)
        print*,"kn_pborn(4) = ", kn_pborn(:,4)
        print*,"center of mass frame:"
        print*,"kn_cmpborn(1) = ", kn_cmpborn(:,1)
        print*,"kn_cmpborn(2) = ", kn_cmpborn(:,2)
        print*,"kn_cmpborn(3) = ", kn_cmpborn(:,3)
        print*,"kn_cmpborn(4) = ", kn_cmpborn(:,4)
#endif
      end

c############### end subroutine born_phsp ##############################

c############### subroutine x1x2phspace ################################
c uses the two random numbers provided in xx
c to generate the fraction of the beam momentum
c for the two partons entering the Born process,
c calculates the partonic Mandelstam variable s
c and the Jacobi factor for the phase space volume
c Parameter to select phase space importance sampling (always flat in y):
c psgen=0:     flat in 1/tau
c psgen=1:     flat in tau
c psgen=2:     flat in log tau
c psgen=3:     flat in log tau (second choice)
      subroutine x1x2phspace(sbeams,minmass,xx,x1,x2,s,jac)
        implicit none

        ! input:
        double precision sbeams,minmass,xx(2)
        ! output, local variables:
        double precision taumin,taumax,tau,y,x1,x2,s,jac,tmp
        integer psgen
        
        ! select phase space importance sampling
        psgen = 3

        ! min and max values of tau
        taumin = minmass**2/sbeams
        taumax = 1d0
        
        ! map xx(1) to tau = x1*x2
        ! with condition:
        ! (m3+m4)**2 <= sborn <= sbeams
        if(psgen.eq.0)then
          ! Sampling flat in 1/tau
          tmp  = 1.0d0/taumax+xx(1)*(1d0/taumin-1d0/taumax)
          tau  = 1.0d0/tmp
          jac  = jac*tau**2*(1d0/taumin-1d0/taumax)
        elseif(psgen.eq.1) then
          ! Sampling flat in tau
          tau  = taumin + xx(1)*(taumax-taumin)
          jac  = jac*(taumax-taumin)
        elseif(psgen.eq.2) then
          ! Flat in log(tau)
          tau = taumin*dexp(xx(1)*dlog(taumax/taumin))
          jac = jac*tau*dabs(dlog(taumax/taumin))
        elseif(psgen.eq.3) then
          ! Flat in log(tau) (second choice, default for dislepton)
          tau = dexp(dlog(taumin)*(1-xx(1)**2))
          jac = jac*tau*dabs(dlog(taumin))*2*xx(1)
        else
         print*, 'Wrong psgen in Born_phsp.F'
         stop
        endif        

#ifdef DEBUG1
        print*,"tau = ", tau
        print*,"jac = ", jac
#endif

        ! map xx(2) to rapidity y
        ! with condition:
        ! 1/2*log(tau) <= y <= -1/2*log(tau)
        y   = -(1-2*xx(2))*dlog(tau)/2
        jac = jac*dlog(tau) ! changed sign in respect to dislepton here

        ! calculate parton momentum fractions
        ! and partonic s
        s  = sbeams*tau
        x1 = dsqrt(tau)*dexp(y)
        x2 = tau/x1

#ifdef DEBUG1
        print*,"y   = ", y
        print*,"jac = ", jac
#endif
      
      end

c############### end subroutine x1x2phspace ############################

c############### subroutine twobody0 ###################################
c uses the random number provided in xx
c to generate four-momenta in the cms frame
c for the particles produced in the Born process,
c calculates the Jacobi factor for the phase space volume
c sampling flat in cos(theta)
      subroutine twobody0(xx3,sborn,m3,m4,p3,p4,jac)
        implicit none

#include "pwhg_math.h"

        ! input: random variable, sborn, masses of outgoing particles
        double precision xx3,sborn,m3,m4
        ! output: momenta of outgoing particles, jacobian
        double precision p3(0:3),p4(0:3),jac
        ! local variables: absoulute momentum, invariant mass, 
        ! cosine, sine theta, kaellen-function, energy of outgoing particles
        double precision k,m,cosQ,sinQ,kaellen,e3,e4

        ! invariant mass
        m = dsqrt(sborn)

        ! sort out bad phase space points
        ! should never be reached since taumin = (m3+m4)**2/sbeams
        if(m.le.m3+m4) then
           jac=0
           return
        endif

        ! calculate absolute of 3-momenta and energies
        ! in center of mass frame of outgoing particles
        kaellen=m**4+m3**4+m4**4-2d0*(m**2*m3**2+m**2*m4**2+m3**2*m4**2)
        e3 = (m**2+m3**2-m4**2)/(2d0*m)
        e4 = (m**2+m4**2-m3**2)/(2d0*m)
        k  = dsqrt(kaellen)/(2d0*m)

        ! map random variable xx3 to cos and sin
        cosQ = 1d0-2d0*xx3
        sinQ = dsqrt(dabs(1-cosQ**2))
        jac  = -jac*2d0 ! changed sign in respect to dislepton here (*)

        ! set up 4-momenta in center of mass frame of outgoing particles
        p3(0) = e3
        p3(1) = sinQ*k
        p3(2) = 0
        p3(3) = cosQ*k

        p4(0) = e4
        p4(1) = -p3(1)
        p4(2) = -p3(2)
        p4(3) = -p3(3)

        ! Achtung: in dislepton stand hier 2*p3(0) statt m=sqrt(s)
        ! das macht für m3=m4 keinen Unterschied, allerdings schon für
        ! Teilchen die unterschiedlich schwer sind!
        jac = jac*k/(8d0*pi*m)

        ! Außerdem wurde hier ein Vorzeichen vergessen (*)
        ! was sich allerdings mit einem falschen Vorzeichen in
        ! x1x2phsp weghob.

      end

c############### subroutine twobody0 ###################################

c############### subroutine born_suppression ###########################
      subroutine born_suppression(fact)
        implicit none

#include "nlegborn.h"

        double precision fact
        double precision powheginput
        external powheginput
        logical init
        data init/.true./
        save init
        
        if (init) then
           if(powheginput("#bornsuppfact").gt.0) then
              fact = powheginput("bornsuppfact")
              print*, "Using Born suppression factor ", fact
           else if(powheginput("#ptsupp").gt.0) then
              fact = powheginput("ptsupp")
              print*, "WARNING: ptsupp is deprecated"
              print*, "Replace it with bornsuppfact"
           else
              fact = 1d0  
           endif
           init=.false.
        endif
        
        if(fact.le.0d0) fact = 1d0
        
#ifdef DEBUGQ
        print*, "bornsuppfact = ", fact
        !stop
#endif
      end

c############### end subroutine born_suppression #######################

c############### subroutine set_fac_ren_scales #########################
      subroutine set_fac_ren_scales(muf1,mur1)
        implicit none

#include "PhysPars.h"
#include "nlegborn.h"
#include "pwhg_kn.h"
#include "pwhg_st.h"

        ! renamed to prevent conflicts with madgraph parameter mur
        double precision muf1,mur1
        double precision muref
        double precision kn34(0:3)
        logical init
        data init/.true./
        logical runningscale
        save runningscale
        double precision powheginput
        external powheginput

        if(init) then
           if(powheginput("#runningscale").gt.0) then
              if (powheginput("runningscale").ne.1) then
                 print*, "you set runningscale to ", runningscale
                 print*, "this runningscale value is not allowed"
                 stop
              endif
              runningscale = .true.
           else
              runningscale = .false.
           endif
        endif
        
        if (runningscale) then
           kn34(:) = kn_pborn(:,3)+kn_pborn(:,4)
           muref = dsqrt(dabs(kn34(0)**2-kn34(1)**2
     &                      -kn34(2)**2-kn34(3)**2))
        else
           muref = dabs(par_Fin1mass) + dabs(par_Fin2mass)
        endif

        if (init) then
          print*
          print*, "===================================================="
          print*, "Factorization and renormalization scales:"
          if (runningscale) then
            print*, "using running scale and"
            print*, "invariant mass of neutralinos for muref"
          else
            print*, "using fixed scale and"
            print*, "sum of neutralino masses for muref"
            print*, "muref = ", muref
          endif
          ! NOTE: scales will be multiplied (outside of this function) by
          ! an additional factor if renscfact or facscfact is set > 0
          ! in powheg.input
          print*, "muren = muref * ", st_renfact
          print*, "mufac = muref * ", st_facfact
          print*, "===================================================="
          init = .false.
        endif

        muref = max(muref,dsqrt(2d0))  ! avoid too small values of the scales

        muf1 = muref
        mur1 = muref

      end

c############### end subroutine set_fac_ren_scales #####################

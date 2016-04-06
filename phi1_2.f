c the definitions used here are according to Byckling & Kajantie,
c ISBN 0 471 12885 6 (the best book for kinematics and phase space)
c
c you can use the subroutines defined here to generate any phase space
c you want:
c
c
c 
c P_b \                  / P_n         / P_n-1           / P_2
c      \           ____ /        ___  /            ____ /
c       \         |    |        |    |            |    |
c        \________| R2 |________| R2 |________..._| R2 |_____ P_1
c        /  P=K_n |____|  K_n-1 |____|  K_n-2     |____|
c       /
c      /
c P_a /
c
c
c
c
c
c              /(M_n-m_n)^2
c             |
c R_n(Mn^2) = | d(M_n-1)^2 R_2(K_n; (K_n-1)^2, P_n^2) R_n-1((M_n-1)^2)
c             |
c             / (mu_n-1)^2
c
c 
c where:
c                                                                           /
c                                                                          |
c R_2(K_n;(K_n-1)^2, P_n^2) = kaellenSqrt(K_n^2,(k_n-1)^2,P_n^2)/(8 K_n^2) | dOmega_n 
c                                                                          |
c                                                                          /  
c K_i  = P_1 + P_2 + ... + P_i
c mu_i = m_1 + m_2 + ... + m_i
c (M_n-1)^2 = (P - P_n), P = P_a + P_b
c 
c and:
c
c sigma_n = 1/F I_n(s)
c F = 2 kaellenSqrt(m_a^2,m_b^2,s) (2*Pi)^(3n-4)
c I_n = 1/(2 S_za + 1) 1/(2 S_zb + 1) R_n(s)  Sum_pol |M|^2 (for unpolarized particles)
c I_n = R_n(s) |M|^2 (for polarized particles)
c
c Notes: The flux factor F and the spin factors should be supplied in
c your phase space routines Born_phsp or Real_osres_phsp if not already 
c done by the POWHEG-BOX.



c############### subroutine x1x2phspace ################################
c uses the two random numbers provided in xx
c to generate the fraction of the beam momentum
c for the two partons entering the Born process,
c calculates the partonic Mandelstam variable s
c and the Jacobi factor for the phase space volume
c Parameter to select phase space importance sampling (always flat in y):
c psgen=0:     flat in 1/tau
c psgen=1:     flat in tau
c psgen=2:     flat in log tau with arbitrary exponent
c psgen=3:     flat in tan tau with arbitrary exponent 
      subroutine x1x2phspace(sbeams,minmass,xx,x1,x2,s,jac)
        implicit none

        ! input:
        double precision sbeams,minmass,xx(2)
        ! output, local variables:
        double precision taumin,taumax,tau,y,x1,x2,s,jac,tmp
        double precision tanmax,tanmin
        integer psgen
        ! choose the exponent for logarithmic sampling
        double precision xexp

        ! reset jacobian
        jac = 1D0

        ! min and max values of tau
        taumin = minmass**2/sbeams
        taumax = 1D0

#ifdef DEBUGQ
#define DEBUG1
        xx(1) = 1D0
#endif        

        ! select phase space importance sampling
        psgen = 2

        ! map xx(1) to tau = x1*x2
        ! with condition:
        ! (m3+m4)**2 <= sborn <= sbeams
        if(psgen.eq.0)then
          ! Sampling flat in 1/tau with arbitrary exponent
          xexp = 0.5D0
          tmp  = 1D0/taumin+xx(1)**xexp*(1D0/taumax-1D0/taumin)
          tau  = 1D0/tmp
          jac  = -jac*xexp*tau**2*(1D0/taumax-1D0/taumin)*
     &                xx(1)**(xexp-1)
        elseif(psgen.eq.1) then
          ! Sampling flat in tau with arbitrary exponent
          xexp = 4D0
          tau  = taumin + xx(1)**xexp*(taumax-taumin)
          jac  = jac*xexp*xx(1)**(xexp-1)*(taumax-taumin)
        elseif(psgen.eq.2) then
          ! Flat in log(tau) with arbitrary exponent
          xexp = 2D0
          tau = taumax*dexp(dlog(taumin/taumax)*(1-xx(1)**xexp))
          jac = -jac*tau*xexp*xx(1)**(xexp-1)*dlog(taumin/taumax)
        else
         print*, 'Wrong psgen in Born_phsp.F'
         stop
        endif        

#ifdef DEBUG1
        print*,"taumin", taumin
        print*,"taumax", taumax
        print*,"tau = ", tau
        print*,"jac = ", jac
        stop
#endif

        ! map xx(2) to rapidity y
        ! with condition:
        ! 1/2*log(tau) <= y <= -1/2*log(tau)
        y   = (1D0-2D0*xx(2))*dlog(tau)/2D0
        jac = -jac*dlog(tau)

        ! calculate parton momentum fractions
        ! and partonic s
        s  = sbeams*tau
        x1 = dsqrt(tau)*dexp(y)
        x2 = tau/x1

#ifdef DEBUG1
        print*,"y   = ", y
        print*,"jac = ", jac
#endif
      
        ! check if NaN occured
        if(isnan(x1) .or. isnan(x2) .or. isnan(s) .or. isnan(jac)) then
          print*,"warning in phi1_2:142: NaN occured"
          jac = 0D0
          return
        endif
      end
c############### end subroutine x1x2phspace ############################

c############### subroutine R2phsp #####################################
c massive particle p0 in rest frame decaying into p1 with mass m1 
c and p2 with mass m2. Vectors returned p1 and p2 are in the frame in
c which p0 is supplied.
c Result is R2(s) = 1/4 * |p|/sqrts * domega
c Expression evaluated is 
c R2(s) = d^3p1/(2 E1) d^3p2/(2 E2) delta^4(p0 - p1 - p2)
c if you don't want to integrate over the azimuthal degree of freedom
c just set xphi to zero - the jacobian will still be correct.
      subroutine R2phsp(xth,xphi,m1,m2,p0,p1,p2,jac)
        implicit none
        ! masses of decaying particles
        double precision m1,m2
        ! integration variables
        double precision xth,xphi
        double precision xexp
        ! momenta of incoming particle
        double precision p0(0:3)
        ! momenta of outgoing particles
        double precision p1(0:3),p2(0:3)
        ! momenta in CMS
        double precision p1_cms(0:3),p2_cms(0:3)
        ! angles
        double precision phi,cosTh,sinTh
        ! energies abs. momenta and invariants
        double precision E1,E2,Pabs,s,sqrts
        ! jacobian
        double precision jac
        ! variables for boosting into rest frame
        double precision beta, vec(1:3),norm
        ! indices
        integer i
        ! functions
        double precision kaellenSqrt,dotp
        external kaellenSqrt,dotp
        ! constants
        double precision m_pi
        parameter (m_pi = 4.D0*datan(1.D0))
        double precision tiny, tiny2
        parameter (tiny = 1d-6)
        parameter (tiny2 = 1d-8)

        ! reset the jacobian
        jac = 1D0

        s = dotp(p0,p0)
        ! sort out bad phase space points
        if ( s .le. ((1D0-tiny2)*(m1+m2)**2) ) then
          print*, "warning: s is less than the sum of provided masses "//
     &            "m1 and m2"
          print*, "s, (m1+m2)**2 =",s,(m1+m2)**2
          s   = 0D0
          jac = 0D0
          return
        ! no warning if not so bad phase space point
        ! but sort out to avoid NaNs
        else if ( (s .le. (m1+m2)**2) .and.
     &            (s .ge. ((1D0-tiny2)*(m1+m2)**2)) ) then
          s   = 0D0
          jac = 0D0
          return
        endif
        sqrts = dsqrt(s)

        ! choose here sampling exponent for theta integration
        xexp  = 1D0
        cosTh = 2D0*xth**xexp-1D0
        sinTh = dsqrt(dabs(1D0-cosTh**2))
        jac   = jac*2D0*xexp*xth**(xexp-1D0)

        phi   = 2D0*m_pi*xphi
        jac   = jac*2D0*m_pi

        E1 = (s+m1**2-m2**2)/(2D0*sqrts)
        E2 = (s+m2**2-m1**2)/(2D0*sqrts)
        Pabs = kaellenSqrt(s,m1**2,m2**2)/(2D0*sqrts)
        
        p1_cms(0) = E1
        p1_cms(1) = Pabs*sinTh*dcos(phi)
        p1_cms(2) = Pabs*sinTh*dsin(phi)
        p1_cms(3) = Pabs*cosTh

        p2_cms(0) = E2
        p2_cms(1) = -p1_cms(1)
        p2_cms(2) = -p1_cms(2)
        p2_cms(3) = -p1_cms(3)
        
        ! call boost(sqrts,p0,p1_cms,p1)
        ! new: use existing POWHEG-routines
        norm = dsqrt(p0(1)**2 + p0(2)**2 + p0(3)**2)
        if(norm.ne.0D0) then
          do i = 1,3
            vec(i) = p0(i)/norm
          enddo
          beta = norm/p0(0)
          call mboost(1,vec,beta,p1_cms(0:3),p1(0:3))
          call mboost(1,vec,beta,p2_cms(0:3),p2(0:3))
        else
          do i = 0,3
            p2(i) = p2_cms(i)
            p1(i) = p1_cms(i)
          enddo
        endif

        ! physical phase space jacobian
        jac = jac*Pabs/(4D0*sqrts)

        ! filter unphysical momenta configurations
        if ( p2(0) .lt. 0D0 .and. p2(0) .gt. -tiny ) jac = 0D0  
        if ( p1(0) .lt. 0D0 .and. p1(0) .gt. -tiny ) jac = 0D0 
        
        ! tests
        if (((p0(0).lt.0D0).or.(p1(0).lt.0D0).or.(p2(0).lt.0D0))
     &       .and.jac.ne.0D0 ) then
          print*,"warning: one of E1,E2,E3 is less than zero"
          print*,"p0",p0(0),p0(0)**2-p0(1)**2-p0(2)**2-p0(3)**2,s
          print*,"p1",p1(0),p1(0)**2-p1(1)**2-p1(2)**2-p1(3)**2
          print*,"p2",p2(0),p2(0)**2-p2(1)**2-p2(2)**2-p2(3)**2
          jac = 0D0
          return
        endif
        
        ! check if NaN occured
        do i=0,3
          if(isnan(p1(i)) .or. isnan(p2(i)) .or. isnan(jac)) then
            print*,"warning in phi1_2:271: NaN occured"
#ifdef DEBUGQ
            print*,"xth",xth
            print*,"xphi",xphi
            print*,"m1",m1
            print*,"m2",m2
            print*,"p0",p0
            print*,"s",s
            print*,"cosTh",cosTh
            print*,"sinTh",sinTh
            print*,"phi",phi
            print*,"E1",E1
            print*,"E2",E2
            print*,"Pabs",Pabs
            print*,"norm",norm
            print*,"beta",beta
            print*,"vec",vec
            print*,"p1",p1
            print*,"p2",p2
            print*,"jac",jac
#endif
            jac = 0D0
            return
          endif
        enddo
      end
c############### end subroutine R2phsp #################################

c############### subroutine R2phsp_s2 ##################################
c This routine differs from R2phsp in that s2 is integrated over.
c Massive particle p0 decaying into p1 mass m1 and p2 mass-squared s2.
c with invariant mass of particle three s2 integrated over.
c s2min is the minimum value of s2.
c Vectors returned p1 and p2 are in the same frame as p0 is supplied.
c Result is R2(s) = 1/4 * |p|/sqrts * domega
c Expression evaluated is 
c R2(s) = ds2 d^3p1/(2 E1) d^3p2/(2 E2) delta^4(p0 - p1 - p2) 
c delta(p2^2-s2)
c Parameter to select phase space importance sampling:
c psgen=0:     flat in s2 (bwmass and bwwidth is not required)
c psgen=1:     breit wigner in s2
c psgen=2:     breit wigner in s2 and flat below resonance
      subroutine R2phsp_s2(psgen,x2,xth,xphi,s2min,m1,bwmass,bwwidth,
     &                     p0,p1,p2,jac)
        implicit none
        ! parameter to select the PS sampling for s2
        integer psgen
        ! masses of decaying particles
        double precision m1,m2
        ! integration variables
        double precision x2,xth,xphi
        double precision xexp
        ! borders of s2-integration
        double precision s2max,s2min
        ! momenta of incoming particle
        double precision p0(0:3)
        ! momenta of outgoing particles
        double precision p1(0:3),p2(0:3)
        ! momenta in CMS
        double precision p2_cms(0:3),p1_cms(0:3)
        ! angles
        double precision phi,cosTh,sinTh
        ! energies abs. momenta and invariants
        double precision E1,E2,Pabs,s,sqrts,s1,s2
        ! jacobians
        double precision jac,jc1
        ! variables for boosting into rest frame
        double precision beta,vec(1:3),norm
        ! breit wigner
        double precision bwmass,bwwidth
        ! indices
        integer i
        ! functions
        double precision kaellenSqrt,dotp
        external kaellenSqrt,dotp
        ! constants
        double precision m_pi
        parameter (m_pi = 4.D0*datan(1.D0))
        double precision tiny
        parameter (tiny = 1d-6)

        ! reset the jacobian
        jac = 1D0
        jc1 = 1D0

        s = dotp(p0,p0)
        ! sort out bad phase space points
        if (s .le. 0D0) then
         print*, "warning: s is less than zero"
         print*, "s =",s
         print*, " => set s to 0 with jacobian 0"
         s   = 0D0
         jac = 0D0
         return
        endif
        sqrts = dsqrt(s)

        ! choose here sampling exponent for theta integration
        xexp  = 1D0
        cosTh = 2D0*xth**xexp-1D0
        sinTh = dsqrt(dabs(1D0-cosTh**2))
        jac   = jac*2D0*xexp*xth**(xexp-1D0)

        phi   = 2D0*m_pi*xphi
        jac   = jac*2D0*m_pi

        ! integration borders
        s1    = m1**2
        s2max = (sqrts-m1)**2

        ! sort out bad phase space points
        if (s2min .gt. s2max) then 
          print*,"warning: s2min should be smaller than s2max"
          s2min = s2max 
          jac   = 0D0
          return
        endif
        
        ! 0 = flat
        ! 1 = breit wigner
        ! 2 = breit wigner and flat below resonance
        if( psgen .lt. 0 .or. psgen .gt. 2) then
          print*,"error in R2phsp_s2: unknown psgen: ", psgen
          stop
        endif
        if( (psgen.ne.0) .and. ((psgen.eq.1) .or. 
     &       (s.ge.(bwmass+m1)**2.and.psgen.eq.2)) ) then
          call breitw(x2,s2min,s2max,bwmass,bwwidth,s2,jc1)
        else
          s2  = (s2max-s2min)*x2+s2min
          jc1 = (s2max-s2min)
        endif
        m2  = dsqrt(s2)
        jac = jac*jc1

        E1 = (s+m1**2-m2**2)/(2D0*sqrts)
        E2 = (s+m2**2-m1**2)/(2D0*sqrts)
        Pabs = kaellenSqrt(s,m1**2,m2**2)/(2D0*sqrts)
        
        p1_cms(0) = E1
        p1_cms(1) = Pabs*sinTh*dcos(phi)
        p1_cms(2) = Pabs*sinTh*dsin(phi)
        p1_cms(3) = Pabs*cosTh

        p2_cms(0) = E2
        p2_cms(1) = -p1_cms(1)
        p2_cms(2) = -p1_cms(2)
        p2_cms(3) = -p1_cms(3)
                
        ! call boost(sqrts,p0,p1_cms,p1)
        ! new: use existing POWHEG-routines
        norm = dsqrt(p0(1)**2 + p0(2)**2 + p0(3)**2)
        if(norm.ne.0D0) then
          do i = 1,3
            vec(i) = p0(i)/norm
          enddo
          beta = norm/p0(0)
          call mboost(1,vec,beta,p1_cms(0:3),p1(0:3))
          call mboost(1,vec,beta,p2_cms(0:3),p2(0:3))
        else
          do i = 0,3
            p2(i) = p2_cms(i)
            p1(i) = p1_cms(i)
          enddo
        endif

        ! physical phase space jacobian
        jac = jac*Pabs/(4D0*sqrts)

        ! filter unphysical momenta configurations
        if ( p2(0) .lt. 0D0 .and. p2(0) .gt. -tiny ) jac = 0D0  
        if ( p1(0) .lt. 0D0 .and. p1(0) .gt. -tiny ) jac = 0D0 
        
        ! tests
        if (((p0(0).lt.0D0).or.(p1(0).lt.0D0).or.(p2(0).lt.0D0))
     &       .and.jac.ne.0D0 ) then
          print*,"warning: one of E1,E2,E3 is less than zero"
          print*,"p0",p0(0),p0(0)**2-p0(1)**2-p0(2)**2-p0(3)**2,s
          print*,"p1",p1(0),p1(0)**2-p1(1)**2-p1(2)**2-p1(3)**2
          print*,"p2",p2(0),p2(0)**2-p2(1)**2-p2(2)**2-p2(3)**2
          jac = 0D0
          return
        endif
        
        ! check if NaN occured
        do i=0,3
          if(isnan(p1(i)) .or. isnan(p2(i)) .or. isnan(jac)) then
            print*,"warning in phi1_2:438: NaN occured"
            jac = 0D0
            return
          endif
        enddo  
      end
c############### end subroutine R2phsp_s2 ##############################

c############### subroutine breitw #####################################
c Given a number 0<x<1 generate a mass-squared msq and a jacobian jac 
c such that mminsq<msq<mmaxsq
c points are generated around resonance position rmass, but 
c breit-wigner should still be included in the matrix element
c jac is the jacobian between integration in msq and integration in x1
      subroutine breitw(x1,mminsq,mmaxsq,rmass,rwidth,msq,jac)       
        implicit none
        double precision x1,mminsq,mmaxsq,rmass,rwidth,msq,jac
        double precision xmin,xmax,x,tanx
        ! constants
        double precision m_pi
        parameter (m_pi = 4.D0*datan(1.D0))

        ! reset jacobian
        jac = 1D0
        
        ! in case the maximum msq is very small
        ! just generate linearly for safety
        if(mmaxsq .lt. rmass*1d-3) then
          msq = mminsq + x1*(mmaxsq - mminsq)
          jac = mmaxsq - mminsq
          return
        endif

        if(rwidth.eq.0D0) then
          tanx = 0D0
          xmax = +m_pi/2D0 
          xmin = -m_pi/2D0 
        else
          xmin = datan((mminsq-rmass**2)/rmass/rwidth)
          xmax = datan((mmaxsq-rmass**2)/rmass/rwidth)
          x    = (xmax-xmin)*x1+xmin
          tanx = dtan(x)
        endif

        msq = rmass**2+rmass*rwidth*tanx
        !bw = (1D0+tanx**2)*rmass**2*rwidth**2
        jac = (xmax-xmin)*rmass*rwidth*(1D0+tanx**2) ! OK
      end
c############### end subroutine breitw #################################
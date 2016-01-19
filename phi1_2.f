c copied and modified from POWHEG-BOX-V2/VBF_Z_Z

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

        ! reset jacobian
        jac = 1D0

        ! min and max values of tau
        taumin = minmass**2/sbeams
        taumax = 1d0
        
        ! select phase space importance sampling
        psgen = 3

        ! map xx(1) to tau = x1*x2
        ! with condition:
        ! (m3+m4)**2 <= sborn <= sbeams
        if(psgen.eq.0)then
          ! Sampling flat in 1/tau
          tmp  = 1d0/taumax+xx(1)*(1d0/taumin-1d0/taumax)
          tau  = 1d0/tmp
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
        y   = -(1D0-2D0*xx(2))*dlog(tau)/2D0
        jac = jac*dlog(tau)

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

c############### subroutine phi1_2m ####################################
c massive particle p0 in rest frame decaying into p1 with mass m1 
c and p2 with mass m2. Vectors returned p1 and p2 are in the frame in
c which p0 is supplied.
c Result is 1/8/pi * 2|p|/sqrts * domega/(4*pi)
c factor of (2*pi)^4 included in definition of phase space.
c Expression evaluated is 
c d^4 p1 d^4 p2 (2 pi)^4 delta(p0-p1-p2)/(2 pi)^6
c delta(p2^2) delta(p3^2)
c if you don't want to integrate over the azimuthal degree of freedom
c just set xphi to zero - the jacobian will be still correct.
      subroutine phi1_2m(xth,xphi,m1,m2,p0,p1,p2,jac)
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
        ! functions
        double precision kaellenSqrt,dotp
        external kaellenSqrt,dotp
        ! constants
        double precision m_pi
        parameter (m_pi = 4.D0*datan(1.D0))
        ! indices
        integer i

        ! reset the jacobian
        jac = 1D0

        s = dotp(p0,p0)
        ! sort out bad phase space points
        if (s .le. (m1+m2)**2) then
         print*, "warning: s is less than the sum of provided masses "//
     &           "m1 and m2"
         print*, "s, (m1+m2)**2 =",s,(m1+m2)**2
         print*, " => set s to 0 with jacobian 0"
         s   = 0D0
         jac = 0D0
         return
        endif
        sqrts = dsqrt(s)

        ! choose here sampling exponent for theta integration
        xexp  = 1D0
        cosTh = 1D0-2D0*xth**xexp
        sinTh = dsqrt(1D0-cosTh**2)
        jac   = -jac*2D0*xexp*xth**(xexp-1D0)

        phi   = 2D0*m_pi*xphi
        !jac   = jac*2D0*m_pi

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

        jac = jac*Pabs/(8d0*m_pi*sqrts)

        ! filter unphysical momenta configurations
        if ( p2(0) .lt. 0D0 .and. p2(0) .gt. -1D-10 ) jac = 0D0  
        if ( p1(0) .lt. 0D0 .and. p1(0) .gt. -1D-10 ) jac = 0D0 
        
        ! tests
        if (((p0(0).lt.0D0).or.(p1(0).lt.0D0).or.(p2(0).lt.0D0))
     &       .and.jac.ne.0D0 ) then
          print*,"error in phi3m"
          print*,"p0",p0(0),p0(0)**2-p0(1)**2-p0(2)**2-p0(3)**2,s
          print*,"p1",p1(0),p1(0)**2-p1(1)**2-p1(2)**2-p1(3)**2
          print*,"p2",p2(0),p2(0)**2-p2(1)**2-p2(2)**2-p2(3)**2
          stop
        endif
      end
c############### end subroutine phi1_2m ################################

c############### subroutine phi1_2m_bw #################################
c This routine differs from phi1_2m in that s2 is always generated
c according to a Breit-Wigner described by the additional
c arguments bwmass and bwwidth.
c massive particle p0 decaying into p1 mass m1 and p2 mass-squared s2.
c with invariant mass of particle three s2 integrated over.
c s2min is the minimum value of s2.
c Vectors returned p1 and p2 are in the same frame as p0 is supplied.
c Expression evaluated is 
c ds2 d^4 p1 d^4 p2 (2 pi)^4 delta(p0-p1-p2)/(2 pi)^6
c delta(p1^2-m1) delta(p2^2-s2)
      subroutine phi1_2m_bw(x2,xth,xphi,s2min,m1,bwmass,bwwidth,
     & p0,p1,p2,jac)
        implicit none
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
        cosTh = 1D0-2D0*xth**xexp
        sinTh = dsqrt(1D0-cosTh**2)
        jac   = -jac*2D0*xexp*xth**(xexp-1D0)

        phi   = 2D0*m_pi*xphi
        !jac   = jac*2D0*m_pi

        ! integration borders
        s1    = m1**2
        s2max = (m1-sqrts)**2

        ! sort out bad phase space points
        if (s2min .gt. s2max) then 
          print*,"warning: s2min should be smaller than s2max"
          s2min = s2max 
          jac   = 0D0
          return
        endif

        call breitw(x2,s2min,s2max,bwmass,bwwidth,s2,jc1)
        m2 = dsqrt(s2)

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

        jac = jac*jc1*Pabs/(8d0*m_pi*sqrts)

        ! filter unphysical momenta configurations
        if ( p2(0) .lt. 0D0 .and. p2(0) .gt. -1D-10 ) jac = 0D0  
        if ( p1(0) .lt. 0D0 .and. p1(0) .gt. -1D-10 ) jac = 0D0 
        
        ! tests
        if (((p0(0).lt.0D0).or.(p1(0).lt.0D0).or.(p2(0).lt.0D0))
     &       .and.jac.ne.0D0 ) then
          print*,"error in phi3m"
          print*,"p0",p0(0),p0(0)**2-p0(1)**2-p0(2)**2-p0(3)**2,s
          print*,"p1",p1(0),p1(0)**2-p1(1)**2-p1(2)**2-p1(3)**2
          print*,"p2",p2(0),p2(0)**2-p2(1)**2-p2(2)**2-p2(3)**2
          stop
        endif
      end
c############### end subroutine phi1_2m_bw #############################

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
        jac = (xmax-xmin)*rmass*rwidth*(1D0+tanx**2)
        return
      end
c############### end subroutine breitw #################################
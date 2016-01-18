c copied and modified from POWHEG-BOX-V2/VBF_Z_Z

c############### subroutine phi1_2m ####################################
c massive particle p0 in rest frame decaying into p1 with mass m1 
c and p2 with mass m2. Vectors returned p1 and p2 are in the frame in
c which p0 is supplied.
c Result is 1/8/pi * 2|p|/sqrts * domega/(4*pi)
c factor of (2*pi)^4 included in definition of phase space.
c Expression evaluated is 
c d^4 p1 d^4 p2 (2 pi)^4 delta(p0-p1-p2)/(2 pi)^6
c delta(p2^2) delta(p3^2)
      subroutine phi1_2m(xth,xphi,m1,m2,p0,p1,p2,jac)
        implicit none
        ! masses of decaying particles
        double precision m1,m2
        ! integration variables
        double precision xth,xphi
        ! momenta of incoming particle
        double precision p0(0:3)
        ! momenta of outgoing particles
        double precision p1(0:3),p2(0:3)
        ! momenta in CMS
        double precision p1_cms(0:3),p2_cms(0:3)
        ! angles
        double precision phi,costh,sinth
        ! energies abs. momenta and invariants
        double precision E1,E2,P1abs,P2abs,s,sqrts
        ! jacobian
        double precision jac0,jac
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
        parameter(jac0=1D0/8D0/m_pi)
        jac=0D0

        s = dotp(p0,p0)
        if (s .lt. 0D0) then
         print*, "warning in phi1_2m: s is less than zero, s =",s
         print*, " => s set to 0 with jacobian 0"
         s   = 0D0
         jac = 0D0
        endif

        sqrts = dsqrt(s)
        costh = 2D0*xth-1D0    
        sinth = dsqrt(1D0-costh**2)
        phi   = 2D0*m_pi*xphi

        jac = jac0

        E1 = (s+m1**2-m2**2)/(2D0*sqrts)
        E2 = (s+m2**2-m1**2)/(2D0*sqrts)
        P1abs = kaellenSqrt(s,m1**2,m2**2)/(2D0*sqrts)
        P2abs = kaellenSqrt(s,m2**2,m1**2)/(2D0*sqrts)
        
        p1_cms(0) = E1
        p1_cms(1) = P1abs*sinth*dsin(phi)
        p1_cms(2) = P1abs*sinth*dcos(phi)
        p1_cms(3) = P1abs*costh

        p2_cms(0) = E2
        p2_cms(1) = -P2abs*sinth*dsin(phi)
        p2_cms(2) = -P2abs*sinth*dcos(phi)
        p2_cms(3) = -P2abs*costh

        
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
          p2(:) = p2_cms(:)
          p1(:) = p1_cms(:)
        endif 

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
        double precision costh,sinth,phi,cphi,sphi
        ! energies abs. momenta and invariants
        double precision E1,E2,P1abs,P2abs,s,sqrts,m1,m2,s1,s2
        ! jacobians
        double precision jac,jac0,jac1,xjac,rtxth
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
        parameter(jac0=1D0/8D0/m_pi)
        jac   = 0D0
        xjac  = 1D0
        rtxth = 1D0
        
        s    = dotp(p0,p0)
        if (s .lt. 0D0) then 
          print*,"s <0" 
          s = 0D0 
          jac = 0D0 
        endif
        
        sqrts = dsqrt(s)
        s1    = m1**2
        s2max = (m1-sqrts)**2
        
        if (s2min .gt. s2max) then 
          print*,"warning: s2min should be smaller s2max"
          s2min = s2max 
          jac = 0D0 
        endif

        call breitw(x2,s2min,s2max,bwmass,bwwidth,s2,jac1) 
        m2 = dsqrt(s2)

        ! choose here sampling exponent
        xexp  = 1D0 
        rtxth = xth**xexp 
        xjac  = 1D0/(xexp*xth**(xexp-1D0))

        costh = 2D0*rtxth-1D0        
        phi   = 2D0*m_pi*xphi
        sinth = dsqrt(1D0-costh**2)
        cphi  = dcos(phi)
        sphi  = dsin(phi)

        E1 = (s+m1**2-m2**2)/(2D0*sqrts)
        E2 = (s+m2**2-m1**2)/(2D0*sqrts)
        P1abs = kaellenSqrt(s,m1**2,m2**2)/(2D0*sqrts)
        P2abs = kaellenSqrt(s,m2**2,m1**2)/(2D0*sqrts)
        
        p2_cms(0) = E2
        p2_cms(1) = P1abs*sinth*sphi
        p2_cms(2) = P1abs*sinth*cphi
        p2_cms(3) = P1abs*costh

        p1_cms(0) = E1
        p1_cms(1) = -P2abs*sinth*sphi
        p1_cms(2) = -P2abs*sinth*cphi
        p1_cms(3) = -P2abs*costh
        
        jac = jac0*jac1*2D0*P1abs/sqrts/xjac
                
        !call boost(sqrts,p0,p2_cms,p2)
        ! new: use existing POWHEG-routines
        norm = dsqrt(p0(1)**2 + p0(2)**2 + p0(3)**2)
        if(norm .ne. 0D0) then
          do i = 1,3
            vec(i) = p0(i)/norm
          enddo
          beta = norm/p0(0)
          call mboost(1,vec,beta,p2_cms(0:3),p2(0:3))
          call mboost(1,vec,beta,p1_cms(0:3),p1(0:3))
        else
          p2(:) = p2_cms(:)
          p1(:) = p1_cms(:)
        endif  

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
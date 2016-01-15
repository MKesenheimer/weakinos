c############### real_osres_phsp.f #####################################
c last modified by MK, 05.12.2015
c weakino pair production
c copied and modified version of disquark/Tools/Procindependent/
c Real_phsp.f

c############### subroutine real_osres_phsp_pwg ########################
c the following subroutine creates a phase-space for 2 -> 3 with
c a tan-mapping for the invariant mass between particles i and j is used
c to calculate the resonant contributions with intermediate resonant squarks.
c The numer "ichan" is an additional channel identifier "ichan" which
c selects the resonance specified in init_processes.F.
c Only perform a tan-mapping if there is a (possible) resonance, 
c otherwise: flat sampling in sij
c this subroutine uses only POWHEG routines.
c
c if the phase-space is to be restricted: this effects only sik 
c -> do so by modifying the jacobian kn_jacborn.
      subroutine real_osres_phsp_pwg(xx,ichan)
        implicit none
        
#include "PhysPars.h"
#include "nlegborn.h"
#include "pwhg_kn.h"
#include "osres.h"

        !include integration over sij -> 6 integration-variables
        double precision xx(ndiminteg)
        ! x1-x2 phase space
        integer psgen
        double precision tau,taumin,taumax
        double precision y,tmp
        ! mandelstams
        double precision s12,sij,sik,sjk
        double precision sjkmax,sjkmin
        double precision sijmin,sijmax
        ! map sij to tan if sij is resonant
        double precision arctan,arctanmin,arctanmax
        ! jacobian
        double precision jac
        ! energies and momenta
        double precision Ei,Ej,Ek,Pi,Pj,Pk
        ! angles
        double precision cosTh,sinTh
        double precision cosXi,sinXi,eta
        ! boost
        double precision vec(3),beta
        double precision mi,mj,mk,mij
        ! indices
        integer ichan,i,j,k
        integer iter, nu
        ! external functions
        double precision kaellenSqrt
        external kaellenSqrt
        ! constants
        double precision m_pi
        parameter (m_pi = 4.D0*datan(1.D0))
        double precision tiny
        parameter (tiny = 1d-12)
        
        ! variables for pwg routines phi1_2, phi1_2m_bw, phi1_2m_nobw...
        double precision ttmp,jc1,jc2,xth,pij(4),p12(4)
        
        ! check the channel identifier
        if( ichan.gt.nosres .or. ichan.lt.0) then
          print*,"wrong channel number: ", ichan
          stop
        endif

        ! set the channel-related indices i,j,k and masses mi,mj,mk,mij
        call set_channel(osresID(ichan),i,j,k,mi,mj,mk,mij)
        
#ifdef DEBUGQ
        print*,osresID(ichan),i,j,k,mi,mj,mk,mij
#endif
        
        ! reset the jacobian
        jac = 1D0

        ! do not restrict the integration over IS when performing the 
        ! subtraction: simply attach a theta-function like 
        ! sqrt(shat)-par_Fin1mass-mij when performing the subtraction
        taumin = (mi+mj+mk)**2/kn_sbeams ! OK
        taumax = 1D0
        
        ! select phase space importance sampling
        psgen = 3

        ! map xx(1) to tau = x1*x2
        ! with condition:
        ! (m3+m4+m5)**2 <= sborn <= sbeams
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
         print*, "Wrong psgen in Real_osres_phsp.F"
         stop
        endif        

#ifdef DEBUGQ
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
        kn_x1 = dsqrt(tau)*dexp(y)
        kn_x2 = tau/kn_x1
        kn_sreal  = kn_sbeams*tau
        s12 = kn_sreal ! for convenience

#ifdef DEBUGQ
        print*,"y   = ", y
        print*,"jac = ", jac
#endif

        ! kinematics of IS-particles already fixed in partonic CMS:
        kn_cmpreal(0,1) = dsqrt(s12)/2D0
        kn_cmpreal(1,1) = 0D0
        kn_cmpreal(2,1) = 0D0
        kn_cmpreal(3,1) = kn_cmpreal(0,1)

        kn_cmpreal(0,2) = kn_cmpreal(0,1)
        kn_cmpreal(1,2) = 0D0
        kn_cmpreal(2,2) = 0D0
        kn_cmpreal(3,2) = -kn_cmpreal(0,1)
  
        ! phi1_2(x1,x2,x3,x4,p1,p2,p3,n2,n3,mass2,mass3,width2,width3,wt)
        ! phi1_2m_bw(m2,x3,xth,xphi,s3min,p1,p2,p3,bwmass,bwwidth,wt)
        ! This routine differs from phi1_2m in that s3 is always generated
        ! according to a Breit-Wigner described by the additional
        ! arguments bwmass and bwwidth
        ! massive particle p1 decaying into p2 mass m2 and p3 mass-squared s3.
        ! with invariant mass of particle three s3 integrated over.
        ! s3min is the minimum value of s3.
        ! Vectors returned p2 and p3 are in the same frame as p1 is supplied.
        ! wt is the jacobian between integration in msq and integration in x1
        p12 = kn_cmpreal(:,1)+kn_cmpreal(:,2)
        sijmin = (mi+mj)**2
        !call phi1_2m_bw(mk,xx(3),xx(4),xx(5),sijmin,p12,kn_cmpreal(:,k),pij,mij,WREG,jc1)
        !call phi1_2m_nobw(mi,xx(6),xx(7),xx(8),mi**2,pij,kn_cmpreal(:,i),kn_cmpreal(:,j),jc2)
        call phi1_2(xx(3),xx(4),xx(5),xx(6),p12,kn_cmpreal(:,k),pij,1,1,mk,mij,0D0,WREG,jc1)
        call phi1_2(xx(7),xx(8),xx(9),xx(10),pij,kn_cmpreal(:,i),kn_cmpreal(:,j),1,1,mi,mj,0D0,0D0,jc2)
        jac = jac*jc1*jc2
     
     
        ! check four-momentum conservation
        call check_4conservation(kn_cmpreal,nlegreal)

        ! physical phase space jacobian
        jac = jac*1D0

        ! flux factors
        jac = jac/(2D0*s12)
        
        ! do not introduce a new jacobian-variable
        ! -> should be unique
        kn_jacborn = jac

        ! Boost to the LAB frame:
        beta = (kn_x1-kn_x2)/(kn_x1+kn_x2)
        vec(1) = 0D0
        vec(2) = 0D0
        vec(3) = 1D0
        call mboost(nlegreal,vec,beta,kn_cmpreal(0,1),kn_preal(0,1))

        call check_4conservation(kn_preal,nlegreal)

      end
c############### end subroutine real_osres_phsp_pwg ####################

c############### subroutine real_osres_phsp ############################
c the following subroutine creates a phase-space for 2 -> 3 with
c a tan-mapping for the invariant mass between particles i and j is used
c to calculate the resonant contributions with intermediate resonant squarks.
c The numer "ichan" is an additional channel identifier "ichan" which
c selects the resonance specified in init_processes.F.
c Only perform a tan-mapping if there is a (possible) resonance, 
c otherwise: flat sampling in sij
c
c if the phase-space is to be restricted: this effects only sik 
c -> do so by modifying the jacobian kn_jacborn.
      subroutine real_osres_phsp(xx,ichan)
        implicit none
        
#include "PhysPars.h"
#include "nlegborn.h"
#include "pwhg_kn.h"
#include "osres.h"

        !include integration over sij -> 6 integration-variables
        double precision xx(ndiminteg)
        ! x1-x2 phase space
        integer psgen
        double precision tau,taumin,taumax
        double precision y,tmp
        ! mandelstams
        double precision s12,sij,sik,sjk
        double precision sjkmax,sjkmin
        double precision sijmin,sijmax
        ! map sij to tan if sij is resonant
        double precision arctan,arctanmin,arctanmax
        ! jacobian
        double precision jac
        ! energies and momenta
        double precision Ei,Ej,Ek,Pi,Pj,Pk
        ! angles
        double precision cosTh,sinTh
        double precision cosXi,sinXi,eta
        ! boost
        double precision vec(3),beta
        double precision mi,mj,mk,mij
        ! indices
        integer ichan,i,j,k
        integer iter, nu
        ! external functions
        double precision kaellenSqrt
        external kaellenSqrt
        ! constants
        double precision m_pi
        parameter (m_pi = 4.D0*datan(1.D0))
        double precision tiny
        parameter (tiny = 1d-12)
        
        ! check the channel identifier
        if( ichan.gt.nosres .or. ichan.lt.0) then
          print*,"wrong channel number: ", ichan
          stop
        endif

        ! set the channel-related indices i,j,k and masses mi,mj,mk,mij
        call set_channel(osresID(ichan),i,j,k,mi,mj,mk,mij)
        
#ifdef DEBUGQ
        print*,osresID(ichan),i,j,k,mi,mj,mk,mij
#endif
        
        ! reset the jacobian
        jac = 1D0

        ! do not restrict the integration over IS when performing the 
        ! subtraction: simply attach a theta-function like 
        ! sqrt(shat)-par_Fin1mass-mij when performing the subtraction
        taumin = (mi+mj+mk)**2/kn_sbeams ! OK
        taumax = 1D0
        
        ! select phase space importance sampling
        psgen = 3

        ! map xx(1) to tau = x1*x2
        ! with condition:
        ! (m3+m4+m5)**2 <= sborn <= sbeams
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
         print*, "Wrong psgen in Real_osres_phsp.F"
         stop
        endif        

#ifdef DEBUGQ
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
        kn_x1 = dsqrt(tau)*dexp(y)
        kn_x2 = tau/kn_x1
        kn_sreal  = kn_sbeams*tau
        s12 = kn_sreal ! for convenience

#ifdef DEBUGQ
        print*,"y   = ", y
        print*,"jac = ", jac
#endif

        ! kinematics of IS-particles already fixed in partonic CMS:
        kn_cmpreal(0,1) = dsqrt(s12)/2D0
        kn_cmpreal(1,1) = 0D0
        kn_cmpreal(2,1) = 0D0
        kn_cmpreal(3,1) = kn_cmpreal(0,1)

        kn_cmpreal(0,2) = kn_cmpreal(0,1)
        kn_cmpreal(1,2) = 0D0
        kn_cmpreal(2,2) = 0D0
        kn_cmpreal(3,2) = -kn_cmpreal(0,1)

        ! perfrom sij-integration via tan-mapping if intermediate
        ! particle is on-shell, if not sample flat
        sijmax = (dsqrt(s12)-mk)**2
        sijmin = (mi+mj)**2       
        if( (s12.ge.(mij+mk)**2) .and. (mij.ge.(mi+mj)) ) then
          arctanmax = datan((sijmax-mij**2)/(mij*WREG))
          arctanmin = datan((sijmin-mij**2)/(mij*WREG))
          arctan = (arctanmax-arctanmin)*xx(3)+arctanmin
          sij = mij**2+mij*WREG*dtan(arctan)
          ! include extra-jacobian from tan-mapping
          jac = jac*(arctanmax-arctanmin)*mij*WREG/dcos(arctan)**2
        else
          sij = (sijmax-sijmin)*xx(3)+sijmin
          jac = jac*(sijmax-sijmin)
        endif

        ! avoid unstable sij-values...
        if(sijmin+tiny.gt.sij .or. sijmax-tiny.lt.sij) then
          jac = 0D0
          return
        endif

        ! integration over sik
        sjkmax = mj**2+mk**2-0.5D0*((sij-s12+mk**2)
     &           *(sij+mj**2-mi**2)-kaellenSqrt(sij,s12,mk**2)
     &           *kaellenSqrt(sij,mj**2,mi**2))/sij
        sjkmin = mj**2+mk**2-0.5D0*((sij-s12+mk**2)
     &           *(sij+mj**2-mi**2)+kaellenSqrt(sij,s12,mk**2)
     &           *kaellenSqrt(sij,mj**2,mi**2))/sij
        sjk = (sjkmax-sjkmin)*xx(4)+sjkmin
        jac = jac*(sjkmax-sjkmin)

        ! use mandelstam relations to calculate sjk
        sik = s12-sjk-sij+mi**2+mj**2+mk**2
        
        ! now set the energies and the momenta of outgoing particles in CMS
        Ei = (s12+mi**2-sjk)/(2D0*dsqrt(s12))
        Ej = (s12+mj**2-sik)/(2D0*dsqrt(s12))
        Ek = (s12+mk**2-sij)/(2D0*dsqrt(s12))
        
        Pi = kaellenSqrt(s12,mi**2,sjk)/(2D0*dsqrt(s12))
        Pj = kaellenSqrt(s12,mj**2,sik)/(2D0*dsqrt(s12))
        Pk = kaellenSqrt(s12,mk**2,sij)/(2D0*dsqrt(s12))
        
        ! relative angle between particles k & i in CMS
        cosXi = (-mj**2-Pk**2-Pi**2+(Ek+Ei-dsqrt(s12))**2)/(2D0*Pk*Pi) !OK
        sinXi = dsqrt(1D0-cosXi**2)
        
        ! reject numerically unstable values
        if(abs(cosXi).gt.(1D0-tiny)) then
          kn_jacborn = 0D0
          return
        endif
        
        ! integration over cosTh
        cosTh = 2D0*xx(5)-1D0
        sinTh = dsqrt(1D0-cosTh**2)
        jac = 2D0*jac

        ! integration over eta
        eta = 2D0*m_pi*xx(6)
        jac = jac*2D0*m_pi

        ! phi is already integrated out
        jac = jac*2D0*m_pi
        
        ! build the 4-momenta in CMS of incoming particles
        kn_cmpreal(0,k) = Ek
        kn_cmpreal(1,k) = Pk*sinTh
        kn_cmpreal(2,k) = 0D0
        kn_cmpreal(3,k) = Pk*cosTh
            
        kn_cmpreal(0,i) = Ei
        kn_cmpreal(1,i) = Pi*(sinXi*sinTh+sinXi*dsin(eta)*cosTh)
        kn_cmpreal(2,i) = Pi*(sinXi*dcos(eta))
        kn_cmpreal(3,i) = Pi*(cosXi*cosTh-sinXi*dsin(eta)*sinTh)

        kn_cmpreal(0,j) = Ej
        kn_cmpreal(1,j) = -(kn_cmpreal(1,i)+kn_cmpreal(1,k))
        kn_cmpreal(2,j) = -(kn_cmpreal(2,i)+kn_cmpreal(2,k))
        kn_cmpreal(3,j) = -(kn_cmpreal(3,i)+kn_cmpreal(3,k))
        
#ifdef DEBUGQ
        print*,"pi",kn_cmpreal(:,i)
        print*,"pj",kn_cmpreal(:,j)
        print*,"pk",kn_cmpreal(:,k)
        print*,"sij   ",sij
        print*,"mij**2",mij**2
        print*
#endif

        ! check four-momentum conservation
        call check_4conservation(kn_cmpreal,nlegreal)

        ! physical phase space jacobian
        jac = jac*(1/8D0)*(1/(2D0*m_pi)**5)*(1/(4D0*s12))

        ! flux factors
        jac = jac/(2D0*s12)
        
        ! do not introduce a new jacobian-variable
        ! -> should be unique
        kn_jacborn = jac

        ! Boost to the LAB frame:
        beta = (kn_x1-kn_x2)/(kn_x1+kn_x2)
        vec(1) = 0D0
        vec(2) = 0D0
        vec(3) = 1D0
        call mboost(nlegreal,vec,beta,kn_cmpreal(0,1),kn_preal(0,1))

        call check_4conservation(kn_preal,nlegreal)

      end
c############### end subroutine real_osres_phsp ########################

c############### function pt2max_osres #################################
c this function defines the scale for the radiation from the osres-reals
c per default simply take the value of mur,muf (which are set to a 
c fixed-values if flg_osres is set!)
      double precision function pt2max_osres()
        implicit none

#include "pwhg_st.h"

        ! per default: use the scale of the hard process 
        ! (so far a "running" scale is not supported here!)
        pt2max_osres = st_mufact2
        
#ifdef DEBUG
        print*
        print*,"DEBUG, real_osres_phsp.F"
        print*,"pt2max_osres",pt2max_osres
        stop
#endif

        return
      end
      
c############### end function pt2max_osres #############################
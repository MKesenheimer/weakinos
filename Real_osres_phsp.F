c############### real_osres_phsp.f #####################################
c last modified by MK, 05.12.2015
c weakino pair production
c copied and modified version of disquark/Tools/Procindependent/
c Real_phsp.f

c############### subroutine real_osres_phsp ############################
c the following subroutine creates a phase-space for 2 -> 3 with
c a tan-mapping for the invariant mass between particles i and j is used
c to calculate the resonant contributions with intermediate resonant squarks.
c The numer "ichan" is an additional channel identifier "ichan" which
c selects the resonance specified in init_processes.F.
c Only perform a tan-mapping if there is a (possible) resonance, 
c otherwise: flat sampling in sij
c
c if the phase-space is to be restricted: this effects only tik 
c -> do so by modifying the jacobian kn_jacborn.
      subroutine real_osres_phsp(xrealosres,ichan)
        implicit none
        
#include "PhysPars.h"
#include "nlegborn.h"
#include "pwhg_math.h"
#include "pwhg_kn.h"
#include "osres.h"

        !include integration over sij -> 6 integration-variables
        double precision xrealosres(ndiminteg)
        integer index, nu
        double precision tau,y
        double precision tau_min,tau_max,y_min,y_max
        double precision s12,sij,sijmin,sijmax,x,xmin,xmax
        double precision jaco
        double precision prest(0:3,1:2),beta_rest(3)
        double precision pij(0:3)
        double precision absp_cms,absp_rest
        double precision cosThij,phib,tik
        double precision tikmax,tikmin
        double precision cosTh1
        double precision vec(3),beta,norm
        double precision mi,mj,mk,mij
        integer ichan,i,j,k
        double precision kaellenSqrt
        external kaellenSqrt
        double precision tiny
        parameter (tiny = 1d-12)
        
        ! check the channel identifier
        if( ichan.gt.nosres ) then
          print*,"wrong channel number: ", ichan
          stop
        endif

        ! set the channel-related indices i,j,k and masses mi,mj,mk,mij
        call set_channel(osresID(ichan),i,j,k,mi,mj,mk,mij)

        ! reset the jacobian
        jaco = 1D0

        ! do not restrict the integration over IS when performing the 
        ! subtraction: simply attach a theta-function like 
        ! sqrt(shat)-par_Fin1mass-mij when performing the subtraction
        tau_min = (mi+mj+mk)**2/kn_sbeams ! OK
        tau_max = 1D0
        tau  = tau_min+xrealosres(1)*(tau_max-tau_min)
        jaco = jaco*(tau_max-tau_min)

        y_min = dlog(tau)/2D0
        y_max = -dlog(tau)/2D0
        y = (y_max-y_min)*xrealosres(2) + y_min
        jaco = jaco*(y_max-y_min)

        ! proton momenta fractions
        ! use kn_x1 here, as this routine is only used for REAL amplitudes
        kn_x1 = dsqrt(tau)*exp(y)
        kn_x2 = dsqrt(tau)*exp(-y)
        kn_sreal = kn_sbeams*tau
        s12 = kn_sreal ! for convenience

        ! kinematics of IS-particles already fixed in partonic CMS:
        kn_cmpreal(0,1) = dsqrt(s12)/2D0
        kn_cmpreal(1,1) = 0D0
        kn_cmpreal(2,1) = 0D0
        kn_cmpreal(3,1) = dsqrt(s12)/2D0

        kn_cmpreal(0,2) = dsqrt(s12)/2D0
        kn_cmpreal(1,2) = 0D0
        kn_cmpreal(2,2) = 0D0
        kn_cmpreal(3,2) = -dsqrt(s12)/2D0

        ! perfrom "sij"-integration via tan-mapping
        ! old:
        sijmax = (dsqrt(s12)-par_Fin1mass)**2
        sijmin = par_Fin2mass**2
        ! new:
        !sijmax = (dsqrt(s12)-mk)**2
        !sijmin = (mi+mj)**2
        ! old from disquark
        if(par_Fin2mass.le.mij) then
        ! new
        !if( (s12.ge.(mij+mk)**2) .and. (mij.ge.(mi+mj)) ) then
          xmax = atan((sijmax-mij**2)/(mij*WREG))
          xmin = atan((sijmin-mij**2)/(mij*WREG))
          x = (xmax-xmin)*xrealosres(3)+xmin
          sij = mij**2+mij*WREG*tan(x)
          !include extra-jacobian from tan-mapping
          jaco = jaco*(xmax-xmin)*mij*WREG/dcos(x)**2
        else
          sij = (sijmax-sijmin)*xrealosres(3)+sijmin
          jaco = jaco*(sijmax-sijmin)
        endif

        ! avoid unstable sij-values...
        if(sijmin+tiny.gt.sij .or. sijmax-tiny.lt.sij) then
          jaco = 0D0
          return
        endif

        ! get tik:
        ! old:
        tikmax = par_Fin1mass**2-0.5D0*((s12-sij+par_Fin1mass**2)
     &        -kaellenSqrt(s12,sij,par_Fin1mass**2))
        tikmin = par_Fin1mass**2-0.5D0*((s12-sij+par_Fin1mass**2)
     &        +kaellenSqrt(s12,sij,par_Fin1mass**2))

        ! new:
!        tikmax = mi**2+mk**2-0.5D0*((sij-s12+mk**2)
!     &           *(sij+mi**2-mj**2)-kaellenSqrt(sij,s12,mk**2)
!     &           *kaellenSqrt(sij,mi**2,mj**2))/sij
!        tikmin = mi**2+mk**2-0.5D0*((sij-s12+mk**2)
!     &           *(sij+mi**2-mj**2)+kaellenSqrt(sij,s12,mk**2)
!     &           *kaellenSqrt(sij,mi**2,mj**2))/sij


        ! hier weiter:
        ! check old: 1.0790077776021418E-002 +-   9.8366725364809160E-004 (mit ncall = 2000, itmx = 1)

        tik = (tikmax-tikmin)*xrealosres(4)+tikmin
        jaco = jaco*(tikmax-tikmin)

        ! get cosThij_b3
        cosThij = 2D0*xrealosres(5)-1D0
        jaco = 2D0*jaco

        ! get phib:
        phib = 2D0*pi*xrealosres(6)
        jaco = jaco*2D0*pi
        
        ! p1 ( =  spectator) in the CMS:
        absp_cms = 0.5D0/dsqrt(s12)*
     &             dsqrt(s12**2-2D0*s12*(sij+par_Fin1mass**2)+
     &             (sij-par_Fin1mass**2)**2)
        ! get cosTh1 from tik:
        cosTh1 = (tik-par_Fin1mass**2+dsqrt(s12)*
     &           dsqrt(absp_cms**2+par_Fin1mass**2))/
     &           dsqrt(s12)/absp_cms

        ! reject numerically unstable values
        if(abs(cosTh1).gt.1D0) then
          kn_jacborn = 0D0
          return
        endif
        ! phi1 = 0D0 ! by convention-> *2pi in the end

        ! now we have p1 (or p2) in the CMS.
        ! depending on which particles form the resonant pair,
        ! either set p3 or p4 as spectator
        kn_cmpreal(0,k) = dsqrt(par_Fin1mass**2+absp_cms**2)
        kn_cmpreal(1,k) = absp_cms*dsqrt(1D0-cosTh1**2)!*dcos(phi1)
        ! kn_cmpreal(2,k) = absp_cms*dsqrt(1D0-cosTh1**2)*dsin(phi1)
        kn_cmpreal(2,k) = 0D0
        kn_cmpreal(3,k) = absp_cms*cosTh1
        ! p2/p3 in the rest-system of the squark:
        absp_rest = abs(sij-par_Fin2mass**2)/2D0/dsqrt(sij) 

        ! the weakino from squark-decay: in the squark-rest-system
        prest(0,1) = dsqrt(absp_rest**2+par_Fin2mass**2)
        prest(1,1) = absp_rest*dsqrt(1D0-cosThij**2)*dcos(phib)
        prest(2,1) = absp_rest*dsqrt(1D0-cosThij**2)*dsin(phib)
        prest(3,1) = absp_rest*cosThij

        ! the q(-bar)
        prest(0,2) =  absp_rest
        prest(1,2) = -prest(1,1)
        prest(2,2) = -prest(2,1)
        prest(3,2) = -prest(3,1)

        ! now boost these particles to the CMS:
        ! new: use existing POWHEG-routines
        do index = 1,3
          vec(index) = -kn_cmpreal(index,k)/dsqrt(sij+absp_cms**2)
        enddo
        beta   = dsqrt(vec(1)**2+vec(2)**2+vec(3)**2)
        vec(1) = vec(1)/beta
        vec(2) = vec(2)/beta
        vec(3) = vec(3)/beta
        call mboost(1,vec,beta,prest(0:3,1),kn_cmpreal(0:3,i))
        call mboost(1,vec,beta,prest(0:3,2),kn_cmpreal(0:3,nlegreal))
        
        
        ! the further phase-space-factors (for mij>par_Fin2mass)
        ! if we restrict the phase-space: this jacobian is modified
        ! in the main-routine
        jaco = jaco/16D0/s12/2D0/s12/2D0/sij*
     &         abs(sij-par_Fin2mass**2)/16D0/pi**4
        
        ! do not introduce a new jacobian-variable
        ! -> should be unique
        kn_jacborn = jaco

        ! Boost to the LAB frame:
        beta = (kn_x1-kn_x2)/(kn_x1+kn_x2)
        vec(1) = 0D0
        vec(2) = 0D0
        vec(3) = 1D0
        call mboost(nlegreal,vec,beta,kn_cmpreal(0,1),kn_preal(0,1))

        ! check four-momentum conservation
        call check_4conservation(kn_cmpreal, nlegborn+1)
        call check_4conservation(kn_preal, nlegborn+1)

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
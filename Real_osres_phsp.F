c############### real_osres_phsp.f #####################################
c last modified by MK, 05.12.2015
c weakino pair production
c copied and modified version of disquark/Tools/Procindependent/
c Real_phsp.f

c############### subroutine real_osres_phsp ############################
c the following subroutine creates a phase-space for 2 -> 3 with
c a tan-mapping for the invariant mass between particles 4 and 5 
c or 3 and 5 is used to calculate the resonant contributions with 
c intermediate resonant squark-contribution.
c Pass an additional channel identifier "chan" which selects the
c resonance specified in init_processes.F.
c in all casses the masses have ALWAYS the meaning 
c par_Fin1mass -> spectator, par_Fin2mass -> weakino from squark-decays.
c Only perform a tan-mapping if there is a (possible) resonance, 
c otherwise: flat sampling in s4
c
c if the phase-space is to be restricted: this effects only t1 
c -> do so by modifying the jacobian kn_jacborn.
      subroutine real_osres_phsp(xrealosres,chan)
        implicit none
        
#include "PhysPars.h"
#include "nlegborn.h"
#include "pwhg_math.h"
#include "pwhg_kn.h"
#include "osres.h"

        !include integration over s4 -> 6 integration-variables
        double precision xrealosres(ndiminteg)
        integer i,j, nu
        integer specnum,resnum
        double precision tau,y,ran1,ran2
        double precision tau_min,tau_max,y_min,y_max
        double precision s4,s4min,s4max,x,xmin,xmax
        double precision jaco
        double precision prest(0:3,1:2),beta_rest(3)
        double precision pij(0:3)
        double precision absp_cms,absp_rest
        double precision cosTh23,phib,t1
        double precision t1max,t1min
        double precision cosTh1
        double precision mij
        double precision vec(3),beta,norm
        integer chan
        double precision kaellenSqrt
        external kaellenSqrt
        
        ! check the channel identifier
        if( chan.gt.nosres ) then
          print*,"wrong channel number: ", chan
          stop
        endif

        ! set the mass of the intermediate particle, selected by chan
        mij = osresM(chan)
        
        ! reset the jacobian
        jaco = 1D0
        ! create the IS-kinematics:
        ran1 = xrealosres(1)
        ran2 = xrealosres(2)

        ! do not restrict the integration over IS when performing the 
        ! subtraction: simply attach a theta-function like 
        ! sqrt(shat)-par_Fin1mass-mij when performing the subtraction
        tau_min = (par_Fin1mass+par_Fin2mass)**2/kn_sbeams

        tau_max = 1D0
        tau  = tau_min+ran1*(tau_max-tau_min)
        jaco = jaco*(tau_max-tau_min)

        y_min = dlog(tau)/2D0
        y_max = -dlog(tau)/2D0
        y = (y_max-y_min)*ran2 + y_min
        jaco = jaco*(y_max-y_min)

        ! proton momenta fractions
        ! use kn_x1 here, as this routine is only used for REAL amplitudes
        kn_x1 = dsqrt(tau)*exp(y)
        kn_x2 = dsqrt(tau)*exp(-y)
        kn_sreal = kn_sbeams*tau

        ! kinematics of IS-particles already fixed in partonic CMS!!:
        kn_cmpreal(0,1) = dsqrt(kn_sreal)/2D0
        kn_cmpreal(1,1) = 0D0
        kn_cmpreal(2,1) = 0D0
        kn_cmpreal(3,1) = dsqrt(kn_sreal)/2D0

        kn_cmpreal(0,2) = dsqrt(kn_sreal)/2D0
        kn_cmpreal(1,2) = 0D0
        kn_cmpreal(2,2) = 0D0
        kn_cmpreal(3,2) = -dsqrt(kn_sreal)/2D0

        ! perfrom "s4"-integration via tan-mapping:
        ! see comment above: masses are always 
        ! par_Fin1mass-> spectator, par_Fin2mass-> from squark-decay
        ! only necessary if par_Fin2mass <= mij
        ! otherwise: set ga to 0, sample s4 flat.
        s4max = (dsqrt(kn_sreal)-par_Fin1mass)**2
        s4min = par_Fin2mass**2
        if(par_Fin2mass.le.mij) then
          xmax = atan((s4max-mij**2)/(mij*WREG))
          xmin = atan((s4min-mij**2)/(mij*WREG))
          x = (xmax-xmin)*xrealosres(3)+xmin
          s4 = mij**2+mij*WREG*tan(x)
          !include extra-jacobian from tan-mapping
          jaco = jaco*(xmax-xmin)*mij*WREG/dcos(x)**2
        else
          s4 = (s4max-s4min)*xrealosres(3)+s4min
          jaco = jaco*(s4max-s4min)
        endif

        ! avoid unstable s4-values...
        if(s4min+1d-12.gt.s4 .or. s4max-1d-12.lt.s4) then
          jaco = 0D0
          return
        endif

        ! get t1:
        t1max = par_Fin1mass**2-0.5D0*((kn_sreal-s4+par_Fin1mass**2)
     &        -kaellenSqrt(kn_sreal,s4,par_Fin1mass**2))
        t1min = par_Fin1mass**2-0.5D0*((kn_sreal-s4+par_Fin1mass**2)
     &        +kaellenSqrt(kn_sreal,s4,par_Fin1mass**2))
         
        t1 = (t1max-t1min)*xrealosres(4)+t1min
        jaco = jaco*(t1max-t1min)

        ! get cosTh23_b3
        cosTh23 = 2D0*xrealosres(5)-1D0
        jaco = 2D0*jaco

        ! get phib:
        phib = 2D0*pi*xrealosres(6)
        jaco = jaco*2D0*pi
        
        ! p1 ( =  spectator) in the CMS:
        absp_cms = 0.5D0/dsqrt(kn_sreal)*
     &             dsqrt(kn_sreal**2-2D0*kn_sreal*(s4+par_Fin1mass**2)+
     &             (s4-par_Fin1mass**2)**2)
        ! get cosTh1 from t1:
        cosTh1 = (t1-par_Fin1mass**2+dsqrt(kn_sreal)*
     &           dsqrt(absp_cms**2+par_Fin1mass**2))/
     &           dsqrt(kn_sreal)/absp_cms
        ! reject numerically unstable values
        if(abs(cosTh1).gt.1D0) then
          kn_jacborn = 0D0
          return
        endif
        !  phi1 = 0D0 ! by convention-> *2pi in the end

        ! now we have p1 (or p2) in the CMS:
        ! depending on which particles form the resonant pair,
        ! either set p3 or p4 (the spectator)
        ! new
        specnum = osreslegs(chan,3)
        resnum  = osreslegs(chan,1)
        ! old
        !if(osreslegs(chan).eq.45) then
        !  specnum = 3
        !  resnum  = 4
        !elseif(osreslegs(chan).eq.35) then
        !  specnum = 4
        !  resnum  = 3
        !else
        !   print*,"wrong channel number: ",chan
        !   print*,"or problem with variable osreslegs: ", osreslegs
        !   stop
        !endif
        
        kn_cmpreal(0,specnum) = dsqrt(par_Fin1mass**2+absp_cms**2)
        kn_cmpreal(1,specnum) = absp_cms*dsqrt(1D0-cosTh1**2)!*dcos(phi1)
        ! kn_cmpreal(2,specnum) = absp_cms*dsqrt(1D0-cosTh1**2)*dsin(phi1)
        kn_cmpreal(2,specnum) = 0D0
        kn_cmpreal(3,specnum) = absp_cms*cosTh1
        ! p2/p3 in the rest-system of the squark:
        absp_rest = abs(s4-par_Fin2mass**2)/2D0/dsqrt(s4) 

        ! the sq from go-decay: in the go-rest-system
        prest(0,1) = dsqrt(absp_rest**2+par_Fin2mass**2)
        prest(1,1) = absp_rest*dsqrt(1D0-cosTh23**2)*dcos(phib)
        prest(2,1) = absp_rest*dsqrt(1D0-cosTh23**2)*dsin(phib)
        prest(3,1) = absp_rest*cosTh23

        ! the qbar
        prest(0,2) = absp_rest
        prest(1,2) = -prest(1,1)
        prest(2,2) = -prest(2,1)
        prest(3,2) = -prest(3,1)

        ! old from disquark:
        ! now boost these particles to the CMS:
        !do i = 1,3
        !   beta_rest(i) = -kn_cmpreal(i,specnum)/dsqrt(s4+absp_cms**2)
        !enddo
        ! resnum is either p4 (for res. in 4_5) or p3 (for res in 3_5)
        !call boost_dis(beta_rest,prest(0:3,1),kn_cmpreal(0:3,resnum))
        !call boost_dis(beta_rest,prest(0:3,2),kn_cmpreal(0:3,nlegreal))
        ! new: use existing POWHEG-routines
        do i = 1,3
           vec(i) = -kn_cmpreal(i,specnum)/dsqrt(s4+absp_cms**2)
        enddo
        beta   = dsqrt(vec(1)**2+vec(2)**2+vec(3)**2)
        vec(1) = vec(1)/beta
        vec(2) = vec(2)/beta
        vec(3) = vec(3)/beta
        call mboost(1,vec,beta,prest(0:3,1),kn_cmpreal(0:3,resnum))
        call mboost(1,vec,beta,prest(0:3,2),kn_cmpreal(0:3,nlegreal))
        
        
        ! the further phase-space-factors (for mij>par_Fin2mass!!!)
        ! if we restrict the phase-space: this jacobian is modified
        ! in the main-routine
        jaco = jaco/16D0/kn_sreal/2D0/kn_sreal/2D0/s4*
     &         abs(s4-par_Fin2mass**2)/16D0/pi**4
        
        ! do not introduce a new jacobian-variable
        ! -> should be unambiguous
        kn_jacborn = jaco

        ! Boost to the LAB frame:
        beta = (kn_x1-kn_x2)/(kn_x1+kn_x2)
        vec(1) = 0D0
        vec(2) = 0D0
        vec(3) = 1D0
        call mboost(nlegreal,vec,beta,kn_cmpreal(0,1),kn_preal(0,1))
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
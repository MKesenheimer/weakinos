c
c this routine is adapted from the POWHEG-BOX/ttJ implementation
c
      subroutine real_phsp(xreal)
      implicit none
#include "PhysPars.h"
#include "nlegborn.h"
#include "pwhg_flst.h"
#include "pwhg_kn.h"
#include "pwhg_math.h"
      double precision xreal(ndiminteg)
      double precision m2,m2_min,m2_max,tau_min,tau,tau_max,y,beta,
     $     vec(3),cth,cth_max,cthdec,phidec,s,z,sl_mass2,sl_mass,
     $     rho,dir(3),kzed
      integer nu,k,j
      logical ini,debug,verbose
      parameter(debug=.false.,verbose=.false.)
      data ini/.true./
      save ini
      double precision tiny
      parameter (tiny=1d-5)

c     - Parameter to select phase space importance sampling
c     - psgen=0: fixed ss masses, sampling flat in m2ss pair mass and in 1/tau
c     - psgen=1: fixed ss masses, sampling flat in 1/m2ss and in 1/tau
c     - psgen=2: fixed ss masses, sampling flat in log(m2ss) and in 1/tau
c     - psgen=4: fixed ss masses, sampling flat in 1/m2ss and in log(tau)
c     - psgen=5: fixed ss masses, sampling flat in log(m2ss) and in log(tau)
c     - psgen=6: slepton masses along a BW (not yet implemented)
      integer psgen
      save psgen

      double precision powheginput
      external powheginput
c     Phase space:
c     d m2/(2 pi) d omega /(8*(2 pi)^2)  (s-m2)/s *
c     d omegadec/(8*(2 pi)^2) * sqrt(1- 4m/m2)

c     m2:  heavy quark pair invariant mass
c     m: heavy quark mass      
c     omega: 3d angle in CM system of jet ss pair
c     omegadec: 3d angle in CM system of ss pair

      if(ini) then
         write(*,*) '************************************'
         write(*,*) '  Minimum kT of jet in underlying  '
         write(*,*) '  real is limited to',kn_ktmin,' GeV'
         write(*,*) '************************************'
         write(*,*)
         psgen=powheginput('#psgen')
         if (psgen.lt.0) psgen=1
      if(verbose) then   
      write(*,*) '****************************************************'
      write(*,*) 'Select phase space importance sampling psgen=',psgen
      write(*,*) '0: sampling flat in m2ss mass, 1/tau'
      write(*,*) '1: sampling flat in 1/m2ss, 1/tau'
      write(*,*) '2: sampling flat in log(m2ss), 1/tau'
      write(*,*) '4: sampling flat in 1/m2ss, log(tau)'
      write(*,*) '5: sampling flat in log(m2ss), log(tau)'
      write(*,*) '>5: slepton masses along BW (not yet implemented)'
      write(*,*) '****************************************************'
      endif
         ini=.false.     
      endif
      
      kn_jacreal=1

      ! set masses for real and real emission:
      kn_masses(1)=0
      kn_masses(2)=0
      kn_masses(3)=dabs(par_Fin1mass) !slepton_mass
      kn_masses(4)=dabs(par_Fin2mass) !slepton_mass
      kn_masses(5)=0
      kn_minmass =  kn_masses(3)+kn_masses(4)

      if(psgen.lt.6) then
c     decide whether to generate the slepton mass along a BW
c     or to keep it fixed
         sl_mass = 0.5d0*(dabs(par_Fin1mass)+dabs(par_Fin2mass))
         sl_mass2 = sl_mass**2
         if(psgen.eq.3) stop "psgen =3 not implemented"
      else
         stop "psgen >5 not yet implemented"
      endif

c   Set minimum and maximum ss invariant mass values; then get:

      m2_min = 4d0*sl_mass2
c   threshold
      m2_max = sqrt(kn_sbeams)*(sqrt(kn_sbeams)-2*kn_ktmin)
c   maximum virtuality allowed for the ss pair
c   if the momentum is kn_ktmin
      if(psgen.eq.0)then
c        Sampling flat in m2
         m2 = m2_min + xreal(1)*(m2_max-m2_min)
         kn_jacreal = kn_jacreal*(m2_max-m2_min)
      elseif((psgen.eq.1).or.(psgen.eq.4)) then
c        Sampling flat in 1/m2
         m2 = 1.0d0/(1.0d0/ m2_max + xreal(1) * (1d0/m2_min-1d0/m2_max))
         kn_jacreal =  kn_jacreal * m2**2 * (1d0/m2_min-1d0/m2_max)
      elseif((psgen.eq.2).or.(psgen.eq.5)) then
c        Sampling flat in log(m2)
         m2 = m2_min * exp(xreal(1)*(1.-tiny) * log(m2_max/m2_min))
         kn_jacreal = kn_jacreal*m2*log(m2_max/m2_min)*(1.-tiny)
      elseif(psgen.gt.5) then
         stop "psgen >5 not yet implemented"
      endif
      
c     d m^2/(2pi) jacobian
      kn_jacreal=kn_jacreal/(2*pi)
c     d x1 d x2 = d tau d y
      tau_min=( dsqrt(m2+kn_ktmin**2) + kn_ktmin )**2/kn_sbeams
      if (tau_min.ge.1d0) then
c     Check if, due to truncation errors, tau_min >1 
         call increasecnt('tau_min > 1')
         print *, "tau_min,m2,kn_ktmin",tau_min,m2,kn_ktmin
         tau_min=1d0-tiny
      endif
      tau_max=1d0
      
      if(psgen.lt.3) then
c   Sampling flat in 1/tau
         tau  = 1.0d0/(1.0d0/ tau_max + xreal(2) * (1d0/tau_min-1d0
     $        /tau_max))
         kn_jacreal = kn_jacreal * tau**2 * (1d0/tau_min-1d0/tau_max)
      elseif ((psgen.ge.4).and.(psgen.lt.6)) then     
c    Flat in log(tau)
         tau = tau_min*exp(xreal(2) * log(tau_max/tau_min))
         kn_jacreal=kn_jacreal*tau*log(tau_max/tau_min)
c        alternatively
c        tau = tau_min*exp(1-xreal(2)**2)
c        kn_jacreal=kn_jacreal*tau*dabs(log(tau_min))*2*xreal(2)
      endif

      if ((tau.lt.tau_min).or.(tau.gt.tau_max)) then
c     Check if, due to truncation errors, tau is wrong 
         call increasecnt('tau not in tau_min < tau < tau_max')
         print *,'tau is not in tau_min < tau < tau_max'
         print *, "tau_min,tau,tau_max,xreal(2)",tau_min, tau,tau_max
     $        ,xreal(2)
         stop
      endif
      
c     Now we can work out \hat{s}
      kn_sreal = tau * kn_sbeams
c     local copy
      s=kn_sreal
c     Set minimum sqrt(\hat{s}) value
      kn_minmass = sqrt(tau_min*kn_sbeams)
c     Compute ss momentum in partonic cm
      kzed=(s-m2)/(2*sqrt(s))
      if(kzed.le.kn_ktmin) then
c     It may happen that due to truncation errors
c     kzed<kn_ktmin. In these cases, put kzed=kn_ktmin 
         kzed=kn_ktmin
      endif
c     Set minimum and maximum y values 
c     ymin,ymax= +/- |log(tau)|/2
c     and cast y uniformly in the allowed range
      y=-(1-2*xreal(3))*log(tau)/2
      kn_jacreal=-kn_jacreal*log(tau)

c     Feynman x's can be computed and saved now:
      kn_xb1 = sqrt(tau) * exp(y)
      kn_xb2 = sqrt(tau) * exp(-y)

c     Check if rounding errors can cause problems 
      if ((kn_xb1.le.0d0).or.(kn_xb1.ge.1d0)) then
         call increasecnt('kn_xb1 outside allowed range')
         print *,'kn_xb1 outside allowed range', kn_xb1
         stop
      endif
      if ((kn_xb2.le.0d0).or.(kn_xb2.ge.1d0)) then
         call increasecnt('kn_xb2 outside allowed range')
         print *,'kn_xb2 outside allowed range', kn_xb2
         stop
       endif
       
c     Set minimum and maximum cos(theta) values (theta is the 
c     angle of ss pair wrt +z axis).
c     abs to protect from tiny negative values
      cth_max=sqrt(abs(1-(kn_ktmin/kzed)**2))
      z=1-2*xreal(4)
      kn_jacreal=kn_jacreal*2d0
      cth=1.5d0*(z-z**3/3)
      kn_jacreal=kn_jacreal*1.5d0*(1-z**2)
      cth=cth*cth_max
      kn_jacreal=kn_jacreal*cth_max
c     Now we know the jet transverse momentum
      kn_born_pt2=(1-cth**2)*kzed**2

c     Supply 2 pi for azimuthal integration (not performed)
      kn_jacreal=kn_jacreal*2*pi
c     Supply physical phase space factor
      kn_jacreal=kn_jacreal*(s-m2)/s/(8*(2*pi)**2)

c     ttb pair decay products in its rest frame
      cthdec=1-2*xreal(5)
      kn_cthdec=cthdec
      kn_jacreal=kn_jacreal*2
      phidec=2*pi*xreal(6)
      kn_jacreal=kn_jacreal*2*pi
      kn_jacreal=kn_jacreal*dsqrt(1d0-(m2_min/m2))/(8*(2*pi)**2)
    
c   build their momenta
      kn_cmpreal(0,3)=dsqrt(m2)/2
      kn_cmpreal(0,4)=kn_cmpreal(0,3)
      kn_cmpreal(3,3)=cthdec
     $     *dsqrt(kn_cmpreal(0,3)**2-sl_mass2)
      kn_cmpreal(1,3)=dsqrt(1-cthdec**2)*dsin(phidec)
     $     *dsqrt(kn_cmpreal(0,3)**2-sl_mass2)
      kn_cmpreal(2,3)=dsqrt(1-cthdec**2)*dcos(phidec)
     $     *dsqrt(kn_cmpreal(0,3)**2-sl_mass2)
      kn_cmpreal(1,4)=-kn_cmpreal(1,3)
      kn_cmpreal(2,4)=-kn_cmpreal(2,3)
      kn_cmpreal(3,4)=-kn_cmpreal(3,3)

c     velocity of ss pair in partonic CM
      beta=(s-m2)/(s+m2)
      vec(1)=dsqrt(1-cth**2)
      vec(2)=0
      vec(3)=cth
      call mboost(2,vec,beta,kn_cmpreal(0,3),kn_cmpreal(0,3))
      kn_cmpreal(1,5)=-kn_cmpreal(1,3)-kn_cmpreal(1,4)
      kn_cmpreal(2,5)=-kn_cmpreal(2,3)-kn_cmpreal(2,4)
      kn_cmpreal(3,5)=-kn_cmpreal(3,3)-kn_cmpreal(3,4)
      kn_cmpreal(0,5)=dsqrt(s)-kn_cmpreal(0,3)-kn_cmpreal(0,4)
      kn_cmpreal(0,1)=dsqrt(s)/2
      kn_cmpreal(0,2)=kn_cmpreal(0,1)
      kn_cmpreal(3,1)=kn_cmpreal(0,1)
      kn_cmpreal(3,2)=-kn_cmpreal(0,2)
      kn_cmpreal(1,1)=0
      kn_cmpreal(1,2)=0
      kn_cmpreal(2,1)=0
      kn_cmpreal(2,2)=0      
c     now boost everything along 3rd axis 
      beta=(kn_xb1-kn_xb2)/(kn_xb1+kn_xb2)
      vec(1)=0
      vec(2)=0
      vec(3)=1
      call mboost(nlegreal-2,vec,beta,kn_cmpreal(0,3),kn_preal(0,3))
      do nu=0,3
         kn_preal(nu,1)=kn_xb1*kn_beams(nu,1)
         kn_preal(nu,2)=kn_xb2*kn_beams(nu,2)
      enddo
      if(debug) then
         call checkmomzeronew(nlegreal,kn_preal)
         call checkmassnew(1,kn_preal)
         call checkmassnew(2,kn_preal)
         call checkmassnew(3,kn_preal)
         call checkmassnew(4,kn_preal)
         call checkmassnew(5,kn_preal)
      endif
      end

      subroutine checkmassnew(n,p)
      implicit none
#include "nlegborn.h"
#include "pwhg_kn.h"
      integer n
      double precision p(0:3,n)
      double precision s(0:3)
      double precision tiny 
      parameter (tiny=1d-5)
      integer nu,j
      
      do nu=0,3
         s(nu)=p(nu,n)
      enddo
      if (kn_masses(n).eq.0d0) then
         if (abs(s(0)**2-s(1)**2-s(2)**2-s(3)**2).gt.tiny) then
            write(*,*) ' mass check not working', n, kn_masses(n)
     $           ,sqrt(abs(s(0)**2-s(1)**2-s(2)**2-s(3)**2))
            call exit(1)
         endif
      else
         if ((abs(sqrt(abs(s(0)**2-s(1)**2-s(2)**2-s(3)**2))
     $        -kn_masses(n))/kn_masses(n)).gt.sqrt(tiny)) then
            write(*,*) ' mass check not working', n, kn_masses(n)
     $           ,sqrt(abs(s(0)**2-s(1)**2-s(2)**2-s(3)**2))
            call exit(1)
         endif
      endif
      end
      
 
      subroutine checkmomzeronew(n,p)
      implicit none
      integer n
      double precision p(0:3,n)
      double precision s(0:3),r(0:3)
      double precision ep
      parameter (ep=1d-10)
      integer nu,j
      do nu=0,3
         s(nu)=0
         r(nu)=0
      enddo
      do j=1,2
         do nu=0,3
            s(nu)=s(nu)+p(nu,j)
            r(nu)=r(nu)+dabs(p(nu,j))
         enddo
      enddo
      do j=3,n
         do nu=0,3
            s(nu)=s(nu)-p(nu,j)
            r(nu)=r(nu)+dabs(p(nu,j))
         enddo
      enddo
      if((s(0)**2+s(1)**2+s(2)**2+s(3)**2).ne.0d0) then
         if((s(0)**2+s(1)**2+s(2)**2+s(3)**2)
     #        /(r(0)**2+r(1)**2+r(2)**2+r(3)**2).gt.ep) then
         write(*,*) ' momentum check not working',s
      endif
      endif
      end

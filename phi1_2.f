      subroutine phi1_2(x1,x2,x3,x4,p1,p2,p3,n2,n3,
     .     mass2,mass3,width2,width3,wt)
c     massive particle p1 decaying into p2 mass m2 and p3 mass m3.
c     with invariant mass 
c     of particle two s2 and particle three s3 integrated over.
c     vectors returned p2 and p3 are in the same frame as p1 is supplied
c     Expression evaluate is 
c     ds2 ds3 d^4 p2 d^4 p3 (2 pi)^4 delta(p1-p2-p3)/(2 pi)^6 
c     delta(p2^2-s2) delta(p3^2-s3)
      implicit none
#include "nlegborn.h"
#include "pwhg_flst.h"
#include "pwhg_kn.h"
#include "pwhg_math.h"
#include "PhysPars.h"
      double precision p1(4),p2(4),p3(4),p3cm(4)
      double precision x1,x2,x3,x4,costh,sinth,phi,cphi,sphi
      double precision wt,wt0,w2,w3
      double precision s2max,s2min,s3max,s3min,one,two
      double precision m1,m2,s1,s2,s3,lambda,mass2,width2,mass3,width3
      integer j,n2,n3
      parameter(wt0=1d0/8d0/pi)
      parameter(one=1d0)
      parameter(two=2d0)
      logical zerowidth
      data zerowidth/.false./
      
      wt=0d0
      s1=p1(4)**2-p1(1)**2-p1(2)**2-p1(3)**2  
      if (s1 .lt. 0d0) stop 'phi1_2: s1 < 0 ' 
      m1=dsqrt(s1)

c--- if both particles are produced on-shell, reject if m1 too small 
c     (should never happen) 
      if (zerowidth.and.(m1 .lt. mass2*dfloat(n2)+mass3*dfloat(n3))) then 
         write(*,*) 'p1', p1 
         write(*,*) 'm1,m2,m3', m1, mass2, mass3, n2, n3 
         stop 'phi1_2: m1 < m2+m3 ?'  
      endif

      s2min=1d-15
      s2max=s1
      if (s2min .gt. s2max) stop 'phi1_2: s2min < s2max' 

      if (n2 .eq. 0) then
C     when generating p3456 make sure invariant mass is enough 
c      to decay to ZZ (ZZ specific)  
         if (zerowidth) then
            s2min = (2d0*par_Zmass)**2
         else
            s2min = 1d-15
         endif
         w2=s2max-s2min
         s2=s2max*x1+s2min*(1d0-x1)
      elseif (n2 .eq. 1) then
         call breitw(x1,s2min,s2max,mass2,width2,s2,w2)       
      endif
      
      m2=dsqrt(s2)
      s3min=1d-15
      s3max=(m2-m1)**2
      if (n3 .eq. 0) then
         w3=s3max-s3min
         s3=s3max*x2+s3min*(1d0-x2)
      elseif (n3 .eq. 1) then
         call breitw(x2,s3min,s3max,mass3,width3,s3,w3)       
      endif

      costh=two*x3-one      
      phi=two*pi*x4
      sinth=dsqrt(one-costh**2)
      cphi=dcos(phi)
      sphi=dsin(phi)
      lambda=((s1-s2-s3)**2-4d0*s2*s3)
c      if (lambda .lt. 0d0) stop 'phi1_2: lambda < 0 ?' 
c BJ test only:
      if (lambda .lt. 0d0) then
      	print*,'phi1_2: lambda = ',lambda
      	lambda = 0d0
      endif	

      lambda=dsqrt(lambda)
      wt=wt0*w2*w3*lambda/s1

      p3cm(4)=m1/two*(s1+s3-s2)/s1
      p3cm(1)=m1/two*lambda/s1*sinth*sphi
      p3cm(2)=m1/two*lambda/s1*sinth*cphi
      p3cm(3)=m1/two*lambda/s1*costh
      call boost(m1,p1,p3cm,p3)
      do j=1,4
      p2(j)=p1(j)-p3(j)
      enddo
      if (  (p1(4) .lt. 0d0) 
     & .or. (p2(4) .lt. 0d0) 
     & .or. (p3(4) .lt. 0d0)) then 
         stop 'phi1_2: one of E1,E2,E3 < 0' 
      endif

      end

      subroutine boost(mass,p1,p_in,p_out)
c     take momenta p_in in frame in which one particle is at rest with mass "mass" 
c     and convert to frame in which the one particle has fourvector p1
      implicit none
      double precision mass,p1(4),p_in(4),p_out(4)
      double precision gam,beta(1:3),bdotp,one
      parameter(one=1d0)
      integer j,k
      gam=p1(4)/mass
      bdotp=0d0
      do j=1,3
      beta(j)=-p1(j)/p1(4)
      bdotp=bdotp+p_in(j)*beta(j)
      enddo
      p_out(4)=gam*(p_in(4)-bdotp)
      do k=1,3
      p_out(k)=p_in(k)+gam*beta(k)*(gam/(gam+one)*bdotp-p_in(4))
      enddo
      return
      end      

      subroutine breitw(x1,mminsq,mmaxsq,rmass,rwidth,msq,wt)       
      implicit none
c---- Given a number 0<x<1 generate a mass-squared msq and a weight wt 
c---- such that mminsq<msq<mmaxsq
c---- points are generated around resonance position rmass, but 
c---- breit-wigner should still be included in the matrix element
c     wt is the jacobian between integration in msq and integration in x1
      double precision x1,mminsq,mmaxsq,rmass,rwidth,msq,wt
      double precision almin,almax,al,tanal
      include 'pwhg_math.h'
      logical zerowidth
      data zerowidth/.false./
      
c--- in case the maximum msq is very small, just generate linearly for safety
      if (mmaxsq .lt. rmass*1d-3) then
        msq=mminsq+x1*(mmaxsq-mminsq)
        wt=mmaxsq-mminsq
        return
      endif

      if (zerowidth) then
          tanal=0d0
          almax=+pi/2d0 
          almin=-pi/2d0 
      else
          almin=datan((mminsq-rmass**2)/rmass/rwidth)
          almax=datan((mmaxsq-rmass**2)/rmass/rwidth)
          al=(almax-almin)*x1+almin
          tanal=dtan(al)
      endif

      msq=rmass**2+rmass*rwidth*tanal
c---- bw=(1d0+tanal**2)*rmass**2*rwidth**2
      wt=(almax-almin)*rmass*rwidth*(1d0+tanal**2)
      return
      end



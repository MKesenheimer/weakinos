c kopiert von Projekt POWHEG-BOX-V2/VBF_Z_Z
      subroutine phi1_2m_nobw(m2,x3,xth,xphi,s3min,p1,p2,p3,wt)
c     massive particle p1 decaying into p2 mass m2 and p3 mass-squared s3.
c     with invariant mass of particle three s3 integrated over.
c     vectors returned p2 and p3 are in the same frame as p1 is supplied
c     Expression evaluate is 
c     ds2 d^4 p2 d^4 p3 (2 pi)^4 delta(p1-p2-p3)/(2 pi)^6
c     delta(p2^2-s2) delta(p3^2-s3)
      implicit none
#include "pwhg_math.h"
      double precision p1(0:3),p2(0:3),p3(0:3),p3cm(0:3)
      double precision x3,xth,xphi,costh,sinth,phi,cphi,sphi
      double precision wt,wt0,w3
      double precision s3max,s3min,one,two 
      double precision m1,m2,m3,s1,s2,s3,lambda,
     . mass2,width2,mass3,width3
      integer j,n2,n3
      parameter(one=1d0)
      parameter(two=2d0)
      parameter(wt0=one/8d0/pi)

      wt=0d0
      s1=p1(0)**2-p1(1)**2-p1(2)**2-p1(3)**2  
      if (s1 .lt. 0d0) then 
         write(*,*) 's1 < 0', s1 
         wt = 0d0 
      endif
      m1=dsqrt(s1)
      s2=m2**2
      s3max=(m2-m1)**2
      if (s3min .gt. s3max) then 
         write(*,*) 'smin , smax', s3min, s3max 
         s3max = s3min 
         wt = 0d0 
      endif
      w3=s3max-s3min
      s3=s3max*x3+s3min*(1d0-x3)

      m3=dsqrt(s3)
      costh=two*xth-one      
      phi=2d0*pi*xphi
      sinth=dsqrt(one-costh**2)
      cphi=dcos(phi)
      sphi=dsin(phi)
      lambda=((s1-s2-s3)**2-4d0*s2*s3)
      if (lambda .lt. 0d0) then
      write(6,*) 'lambda in phi1_2m',lambda
      write(6,*) 's1 in phi1_2m',s1
      write(6,*) 's2 in phi1_2m',s2
      write(6,*) 's3 in phi1_2m',s3
      write(6,*) 'm1 in phi1_2m',m1
      write(6,*) 'm2 in phi1_2m',m2
      write(6,*) 'm3 in phi1_2m',m3
      write(6,*) 'x3 in phi1_2m',x3
      write(6,*) 'mass3 in phi1_2m',mass3
      wt = 0d0 
      endif
      lambda=dsqrt(lambda)
      wt=wt0*w3*lambda/s1

      p3cm(0)=m1/two*(s1+s3-s2)/s1
      p3cm(1)=m1/two*lambda/s1*sinth*sphi
      p3cm(2)=m1/two*lambda/s1*sinth*cphi
      p3cm(3)=m1/two*lambda/s1*costh
      call boost(m1,p1,p3cm,p3)
      do j=0,3
      p2(j)=p1(j)-p3(j)
      enddo
      if (  (p1(0) .lt. 0d0) 
     & .or. (p2(0) .lt. 0d0) 
     & .or. (p3(0) .lt. 0d0)) then  
      write(6,*) 'E1',p1(0)
      write(6,*) 'E2',p2(0)
      write(6,*) 'E3',p3(0)
      write(6,*) 'p1sq',p1(0)**2-p1(1)**2-p1(2)**2-p1(3)**2,s1
      write(6,*) 'p2sq',p2(0)**2-p2(1)**2-p2(2)**2-p2(3)**2,s2
      write(6,*) 'p3sq',p3(0)**2-p3(1)**2-p3(2)**2-p3(3)**2,s3
      write(6,*) 'in phi1_2m.f'
      endif
      return
      end




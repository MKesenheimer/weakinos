cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c KINEMATIC FUNCTIONS THAT CAN BE PLOTTED OR USED FOR CUTS
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c To plot a function give the name of the option in the td_card.dat
c
c To use a function to apply cuts on the events use the function name and
c 'min' or 'max' in the td_card.dat
c    Example:      ptmin 1 1 30d0   The first particle in the first class needs
c                                   at least 30 GeV of pt to pass the cuts. 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      double precision function pt(p)
c************************************************************************
c     Returns transverse momentum of particle
c************************************************************************
      IMPLICIT NONE
      double precision p(0:3)

      pt = dsqrt(p(1)**2+p(2)**2)

      end

      double precision function et(p)
c************************************************************************
c     Returns transverse energy of particle
c************************************************************************
      IMPLICIT NONE
      double precision p(0:3)
      real*8 pt

      pt = dsqrt(p(1)**2+p(2)**2)
      if (pt .gt. 0d0) then
         et = p(0)*pt/dsqrt(pt**2+p(3)**2)
      else
         et = 0d0
      endif
      end

      double precision function e(p)
c************************************************************************
c     Returns energy of particle in lab frame
c************************************************************************
      IMPLICIT NONE
      double precision p(0:3)

      e = p(0)

      end

      double precision function mom(p)
c************************************************************************
c     Returns absolute value of momentum of particle in lab frame
c************************************************************************
      IMPLICIT NONE
      double precision p(0:3)

      mom = dsqrt(p(1)**2+p(2)**2+p(3)**2)

      end

      DOUBLE PRECISION  FUNCTION y(p)
c************************************************************************
c     Returns rapidity of particle in the lab frame
c************************************************************************
      IMPLICIT NONE
      double precision  p(0:3)
      double precision pm

      pm = p(0)
      y = .5d0*dlog((pm+p(3))/(pm-p(3)))

      end

       DOUBLE PRECISION  FUNCTION eta(p)
c************************************************************************
c     Returns pseudo-rapidity of particle in the lab frame
c************************************************************************
      IMPLICIT NONE
      double precision  p(0:3)
      double precision costh

      if (abs(costh(p)).lt.0.9999999999d0)then
         eta = -dlog(dtan(0.5d0*dacos(costh(p))))
      elseif (costh(p).ge.0.9999999999d0)then
         eta =  10d0
      else
         eta = -10d0
      endif

      end

       DOUBLE PRECISION  FUNCTION costh(p)
c************************************************************************
c     Returns the cosine of the angle between the particle and the +z axis
c************************************************************************
      IMPLICIT NONE
      double precision  p(0:3)
      double precision pt,pm

      pm = dsqrt(p(3)**2+pt(p)**2)
      if (abs(pm-p(3)).ge.abs(1d-10*p(3)))then
         costh = p(3)/pm
      else
         costh = 1d0
      endif

      end

      double precision function phi(p)
c************************************************************************
c     Returns azimuthal angle phi of a particle
c************************************************************************
      implicit none
      double precision p(0:3)
      double precision denom,temp
      double precision pi
      parameter ( pi=3.1415926535897932385 )

      denom = dsqrt(p(1)**2 + p(2)**2)
      temp = max(-0.99999999d0, p(1) / denom)
      temp = min( 0.99999999d0, temp)
      if (p(2).ge.0d0)then
         phi =  dacos(temp)
      else
         phi = -dacos(temp) + 2d0*pi
      endif

      end

      double precision function dRij(P1,P2)
c************************************************************************
c     Distance in eta,phi between two particles.
c************************************************************************
      IMPLICIT NONE
      double precision p1(0:3),p2(0:3)
      double precision ETA,DELTA_PHI,R2
      external eta,delta_phi

      R2 = (DELTA_PHI(P1,P2))**2+(ETA(p1)-ETA(p2))**2
      dRij = dsqrt(R2)

      RETURN
      END

      double precision function mij(P1,P2)
c************************************************************************
c     Invarient mass of 2 particles
c************************************************************************
      IMPLICIT NONE
      double precision p1(0:3),p2(0:3)
      double precision sumdot
      external sumdot

      mij = dsqrt(Sumdot(p1,p2,1d0))

      RETURN
      END

      double precision function delta_phi(p1,p2)
c************************************************************************
c     Returns separation in phi of two particles
c************************************************************************
      implicit none
      double precision p1(0:3),p2(0:3)
      REAL*8 DENOM, TEMP

      DENOM = SQRT(P1(1)**2 + P1(2)**2) * SQRT(P2(1)**2 + P2(2)**2)
      TEMP = MAX(-0.99999999D0, (P1(1)*P2(1) + P1(2)*P2(2)) / DENOM)
      TEMP = MIN( 0.99999999D0, TEMP)
      DELTA_PHI = ACOS(TEMP)

      END

      double precision function delta_eta(p1,p2)
c************************************************************************
c     Returns separation in eta of two particles
c************************************************************************
      implicit none
      double precision p1(0:3),p2(0:3)
      double precision eta

      delta_eta=abs(eta(p1)-eta(p2))

      end

      double precision function kTij(p1,p2)
c************************************************************************
c     Returns kT distance of two particles (jet algorithm)
c************************************************************************
      implicit none
      double precision p1(0:3),p2(0:3)
      double precision pt,dRij

      kTij = max( pt(p1) , pt(p2) )*dRij(p1,p2)

      end





      double precision function X1(p1)
C****************************************************************************
C     User function for plotting and cuts
C****************************************************************************
      implicit none
      double precision p1(0:3)

      X1=0d0

      end

      double precision function X2(p1)
C****************************************************************************
C     User function for plotting and cuts
C****************************************************************************
      implicit none
      double precision p1(0:3)

      X2=0d0

      end

      double precision function X3(p1)
C****************************************************************************
C     User function for plotting and cuts
C****************************************************************************
      implicit none
      double precision p1(0:3)

      X3=0d0

      end

      double precision function XY1(p1,p2)
C****************************************************************************
C     User function for plotting and cuts
C****************************************************************************
      implicit none
      double precision p1(0:3),p2(0:3)

      XY1=0d0

      end

      double precision function XY2(p1,p2)
C****************************************************************************
C     User function for plotting and cuts
C****************************************************************************
      implicit none
      double precision p1(0:3),p2(0:3)

      XY2=0d0

      end

      double precision function XY3(p1,p2)
C****************************************************************************
C     User function for plotting and cuts
C****************************************************************************
      implicit none
      double precision p1(0:3),p2(0:3)

      XY3=0d0

      end

      double precision function XYZ1(p1,p2,p3)
C****************************************************************************
C     User function for plotting and cuts
C****************************************************************************
      implicit none
      double precision p1(0:3),p2(0:3),p3(0:3)

      XYZ1=0d0

      end

      double precision function XYZ2(p1,p2,p3)
C****************************************************************************
C     User function for plotting and cuts
C****************************************************************************
      implicit none
      double precision p1(0:3),p2(0:3),p3(0:3)

      XYZ2=0d0

      end

      double precision function XYZ3(p1,p2,p3)
C****************************************************************************
C     User function for plotting and cuts
C****************************************************************************
      implicit none
      double precision p1(0:3),p2(0:3),p3(0:3)

      XYZ3=0d0

      end

      double precision function XYZA1(p1,p2,p3,p4)
C****************************************************************************
C     User function for plotting and cuts
C****************************************************************************
      implicit none
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3)

      XYZA1=0d0

      end

      double precision function XYZA2(p1,p2,p3,p4)
C****************************************************************************
C     User function for plotting and cuts
C****************************************************************************
      implicit none
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3)

      XYZA2=0d0

      end

      double precision function XYZA3(p1,p2,p3,p4)
C****************************************************************************
C     User function for plotting and cuts
C****************************************************************************
      implicit none
      double precision p1(0:3),p2(0:3),p3(0:3),p4(0:3)

      XYZA3=0d0

      end






cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c AUXILIARY FUNCTIONS & SUBROUTINES
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function dot(p1,p2)
C****************************************************************************
C     4-Vector Dot product
C****************************************************************************
      implicit none
      double precision p1(0:3),p2(0:3)
      dot=p1(0)*p2(0)-p1(1)*p2(1)-p1(2)*p2(2)-p1(3)*p2(3)
      end

      double precision function SumDot(P1,P2,dsign)
c************************************************************************
c     Invarient mass of 2 particles
c************************************************************************
      IMPLICIT NONE
      double precision p1(0:3),p2(0:3),dsign
      integer i
      double precision ptot(0:3)
      double precision dot
      external dot

      do i=0,3
         ptot(i)=p1(i)+dsign*p2(i)
      enddo
      SumDot = dot(ptot,ptot)

      RETURN
      END

      subroutine boostx(p,q , pboost)
c
c This subroutine performs the Lorentz boost of a four-momentum.  The
c momentum p is assumed to be given in the rest frame of q.  pboost is
c the momentum p boosted to the frame in which q is given.  q must be a
c timelike momentum.
c
c input:
c       real    p(0:3)         : four-momentum p in the q rest  frame
c       real    q(0:3)         : four-momentum q in the boosted frame
c
c output:
c       real    pboost(0:3)    : four-momentum p in the boosted frame
c
      implicit none
      double precision p(0:3),q(0:3),pboost(0:3),pq,qq,m,lf

      double precision rZero
      parameter( rZero = 0.0d0 )

      qq = q(1)**2+q(2)**2+q(3)**2

      if ( qq.ne.rZero ) then
         pq = p(1)*q(1)+p(2)*q(2)+p(3)*q(3)
         m = sqrt(q(0)**2-qq)
         lf = ((q(0)-m)*pq/qq+p(0))/m
         pboost(0) = (p(0)*q(0)+pq)/m
         pboost(1) =  p(1)+q(1)*lf
         pboost(2) =  p(2)+q(2)*lf
         pboost(3) =  p(3)+q(3)*lf
      else
         pboost(0) = p(0)
         pboost(1) = p(1)
         pboost(2) = p(2)
         pboost(3) = p(3)
      endif
c
      return
      end



      double precision function ordering(p)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This function is used to order the particles within a class
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      double precision p(0:3)
      double precision pt
      external pt

      ordering = pt(p)
      
      end
      


      subroutine setplotrange
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Here you can set the range to plot and the bin size.
c     ...range(1)    is bin size
c     ...range(2)    is smallest value to plot
c     ...range(3)    is largest value to plot
c
c     Max number of bins is 100
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'info.inc'
      integer l1
      double precision binsize, minvalue,maxvalue
      character*10 plotname
      character*140 buffer

c-- pt
      ptrange(1)=4d0
      ptrange(2)=0d0
      ptrange(3)=200d0
c-- et
      etrange(1)=4d0
      etrange(2)=0d0
      etrange(3)=200d0
c-- e
      erange(1)=4d0
      erange(2)=0d0
      erange(3)=200d0
c-- mom
      momrange(1)=4d0
      momrange(2)=0d0
      momrange(3)=200d0
c-- etmiss
      etmissrange(1)=5d0
      etmissrange(2)=0d0
      etmissrange(3)=250d0
c-- eta
      etarange(1)=0.5d0
      etarange(2)=-5d0
      etarange(3)=5d0
c-- y
      yrange(1)=0.5d0
      yrange(2)=-5d0
      yrange(3)=5d0
c-- costh
      costhrange(1)=0.1d0
      costhrange(2)=-1d0
      costhrange(3)=1d0
c-- phi
      phirange(1)=0.3d0
      phirange(2)=0d0
      phirange(3)=6.3d0
c-- dRij
      dRijrange(1)=0.1d0
      dRijrange(2)=0d0
      dRijrange(3)=5d0
c-- mij
      mijrange(1)=5d0
      mijrange(2)=0d0
      mijrange(3)=500d0
c-- delta_phi
      delta_phirange(1)=0.1d0
      delta_phirange(2)=0d0
      delta_phirange(3)=3.1d0
c-- delta_eta
      delta_etarange(1)=0.1d0
      delta_etarange(2)=0d0
      delta_etarange(3)=7d0
c-- kTij
      kTijrange(1)=4d0
      kTijrange(2)=0d0
      kTijrange(3)=200d0
c-- User defined
c  one momentum
      X1range(1)=1d0
      X1range(2)=0d0
      X1range(3)=100d0

      X2range(1)=1d0
      X2range(2)=0d0
      X2range(3)=100d0

      X3range(1)=1d0
      X3range(2)=0d0
      X3range(3)=100d0
c  two momenta
      XY1range(1)=1d0
      XY1range(2)=0d0
      XY1range(3)=100d0

      XY2range(1)=1d0
      XY2range(2)=0d0
      XY2range(3)=100d0

      XY3range(1)=1d0
      XY3range(2)=0d0
      XY3range(3)=100d0
c  three momenta
      XYZ1range(1)=1d0
      XYZ1range(2)=0d0
      XYZ1range(3)=100d0

      XYZ2range(1)=1d0
      XYZ2range(2)=0d0
      XYZ2range(3)=100d0

      XYZ3range(1)=1d0
      XYZ3range(2)=0d0
      XYZ3range(3)=100d0
c  four momenta
      XYZA1range(1)=1d0
      XYZA1range(2)=0d0
      XYZA1range(3)=100d0

      XYZA2range(1)=1d0
      XYZA2range(2)=0d0
      XYZA2range(3)=100d0

      XYZA3range(1)=1d0
      XYZA3range(2)=0d0
      XYZA3range(3)=100d0

      rewind 97
 56   read(97,'(a)',end=85,err=85)buffer
      if(buffer(1:18).eq. "# Begin PlotRange ")then
 55      read(97,'(a)')buffer
            if (buffer(1:16).eq. "# End PlotRange ") goto 58 ! End of PlotOptions
            l1=40
            if(index(buffer,"#").ne.0) l1=index(buffer,"#")-1 ! ignore comments
            read(buffer(1:l1),*,end=57,err=57) plotname,binsize,minvalue,maxvalue
 57         continue
            if (plotname(1:3).eq.'pt ')then
               ptrange(1)=binsize
               ptrange(2)=minvalue
               ptrange(3)=maxvalue
            elseif (plotname(1:3).eq.'et ')then
               etrange(1)=binsize
               etrange(2)=minvalue
               etrange(3)=maxvalue
            elseif (plotname(1:2).eq.'e ')then
               erange(1)=binsize
               erange(2)=minvalue
               erange(3)=maxvalue
            elseif (plotname(1:4).eq.'mom ')then
               momrange(1)=binsize
               momrange(2)=minvalue
               momrange(3)=maxvalue
            elseif (plotname(1:7).eq.'etmiss ')then
               etmissrange(1)=binsize
               etmissrange(2)=minvalue
               etmissrange(3)=maxvalue
            elseif (plotname(1:4).eq.'eta ')then
               etarange(1)=binsize
               etarange(2)=minvalue
               etarange(3)=maxvalue
            elseif (plotname(1:2).eq.'y ')then
               yrange(1)=binsize
               yrange(2)=minvalue
               yrange(3)=maxvalue
            elseif (plotname(1:6).eq.'costh ')then
               costhrange(1)=binsize
               costhrange(2)=minvalue
               costhrange(3)=maxvalue
            elseif (plotname(1:4).eq.'phi ')then
               phirange(1)=binsize
               phirange(2)=minvalue
               phirange(3)=maxvalue
            elseif (plotname(1:3).eq.'X1 ')then
               X1range(1)=binsize
               X1range(2)=minvalue
               X1range(3)=maxvalue
            elseif (plotname(1:3).eq.'X2 ')then
               X2range(1)=binsize
               X2range(2)=minvalue
               X2range(3)=maxvalue
            elseif (plotname(1:3).eq.'X3 ')then
               X3range(1)=binsize
               X3range(2)=minvalue
               X3range(3)=maxvalue
            elseif (plotname(1:5).eq.'dRij ')then
               dRijrange(1)=binsize
               dRijrange(2)=minvalue
               dRijrange(3)=maxvalue
            elseif (plotname(1:4).eq.'mij ')then
               mijrange(1)=binsize
               mijrange(2)=minvalue
               mijrange(3)=maxvalue
            elseif (plotname(1:10).eq.'delta_phi ')then
               delta_phirange(1)=binsize
               delta_phirange(2)=minvalue
               delta_phirange(3)=maxvalue
            elseif (plotname(1:10).eq.'delta_eta ')then
               delta_etarange(1)=binsize
               delta_etarange(2)=minvalue
               delta_etarange(3)=maxvalue
            elseif (plotname(1:5).eq.'kTij ')then
               kTijrange(1)=binsize
               kTijrange(2)=minvalue
               kTijrange(3)=maxvalue
            elseif (plotname(1:4).eq.'XY1 ')then
               XY1range(1)=binsize
               XY1range(2)=minvalue
               XY1range(3)=maxvalue
            elseif (plotname(1:4).eq.'XY2 ')then
               XY2range(1)=binsize
               XY2range(2)=minvalue
               XY2range(3)=maxvalue
            elseif (plotname(1:4).eq.'XY3 ')then
               XY3range(1)=binsize
               XY3range(2)=minvalue
               XY3range(3)=maxvalue
            elseif (plotname(1:5).eq.'XYZ1 ')then
               XYZ1range(1)=binsize
               XYZ1range(2)=minvalue
               XYZ1range(3)=maxvalue
            elseif (plotname(1:5).eq.'XYZ2 ')then
               XYZ2range(1)=binsize
               XYZ2range(2)=minvalue
               XYZ2range(3)=maxvalue
            elseif (plotname(1:5).eq.'XYZ3 ')then
               XYZ3range(1)=binsize
               XYZ3range(2)=minvalue
               XYZ3range(3)=maxvalue
            elseif (plotname(1:6).eq.'XYZA1 ')then
               XYZA1range(1)=binsize
               XYZA1range(2)=minvalue
               XYZA1range(3)=maxvalue
            elseif (plotname(1:6).eq.'XYZA2 ')then
               XYZA2range(1)=binsize
               XYZA2range(2)=minvalue
               XYZA2range(3)=maxvalue
            elseif (plotname(1:6).eq.'XYZA3 ')then
               XYZA3range(1)=binsize
               XYZA3range(2)=minvalue
               XYZA3range(3)=maxvalue
            endif
         goto 55  !Not yet at end of PlotOptions -> Read next line
 58      continue               !Found all PlotOptions ->  continue here
      endif     
      goto 56   !Go back to beginning to read next line
 85   continue                  !End of file -> continue




      end





C
C Intermediate functions
C
      double precision function pt_t(p,ss1)
      implicit none
      include "info.inc"
      double precision p(0:4,maxparticles),pt
      integer a1,ss1
      a1=eventinfo(ss1/1000,mod(ss1,1000))
      pt_t=pt(p(0,a1))
      return
      end

      double precision function et_t(p,ss1)
      implicit none
      include "info.inc"
      double precision p(0:4,maxparticles),et
      integer a1,ss1
      a1=eventinfo(ss1/1000,mod(ss1,1000))
      et_t=et(p(0,a1))
      return
      end

      double precision function e_t(p,ss1)
      implicit none
      include "info.inc"
      double precision p(0:4,maxparticles),e
      integer a1,ss1
      a1=eventinfo(ss1/1000,mod(ss1,1000))
      e_t=e(p(0,a1))
      return
      end

      double precision function mom_t(p,ss1)
      implicit none
      include "info.inc"
      double precision p(0:4,maxparticles),mom
      integer a1,ss1
      a1=eventinfo(ss1/1000,mod(ss1,1000))
      mom_t=mom(p(0,a1))
      return
      end

      double precision function eta_t(p,ss1)
      implicit none
      include "info.inc"
      double precision p(0:4,maxparticles),eta
      integer a1,ss1
      a1=eventinfo(ss1/1000,mod(ss1,1000))
      eta_t=eta(p(0,a1))
      return
      end

      double precision function y_t(p,ss1)
      implicit none
      include "info.inc"
      double precision p(0:4,maxparticles),y
      integer a1,ss1
      a1=eventinfo(ss1/1000,mod(ss1,1000))
      y_t=y(p(0,a1))
      return
      end

      double precision function costh_t(p,ss1)
      implicit none
      include "info.inc"
      double precision p(0:4,maxparticles),costh
      integer a1,ss1
      a1=eventinfo(ss1/1000,mod(ss1,1000))
      costh_t=costh(p(0,a1))
      return
      end

      double precision function phi_t(p,ss1)
      implicit none
      include "info.inc"
      double precision p(0:4,maxparticles),phi
      integer a1,ss1
      a1=eventinfo(ss1/1000,mod(ss1,1000))
      phi_t=phi(p(0,a1))
      return
      end

      double precision function X1_t(p,ss1)
      implicit none
      include "info.inc"
      double precision p(0:4,maxparticles),X1
      integer a1,ss1
      a1=eventinfo(ss1/1000,mod(ss1,1000))
      X1_t=X1(p(0,a1))
      return
      end

      double precision function X2_t(p,ss1)
      implicit none
      include "info.inc"
      double precision p(0:4,maxparticles),X2
      integer a1,ss1
      a1=eventinfo(ss1/1000,mod(ss1,1000))
      X2_t=X2(p(0,a1))
      return
      end

      double precision function X3_t(p,ss1)
      implicit none
      include "info.inc"
      double precision p(0:4,maxparticles),X3
      integer a1,ss1
      a1=eventinfo(ss1/1000,mod(ss1,1000))
      X3_t=X3(p(0,a1))
      return
      end

      double precision function dRij_t(p,ss1,ss2)
      implicit none
      include "info.inc"
      double precision p(0:4,maxparticles),dRij
      integer a1,a2,ss1,ss2
      a1=eventinfo(ss1/1000,mod(ss1,1000))
      a2=eventinfo(ss2/1000,mod(ss2,1000))
      dRij_t=dRij(p(0,a1),p(0,a2))
      return
      end

      double precision function mij_t(p,ss1,ss2)
      implicit none
      include "info.inc"
      double precision p(0:4,maxparticles),mij
      integer a1,a2,ss1,ss2
      a1=eventinfo(ss1/1000,mod(ss1,1000))
      a2=eventinfo(ss2/1000,mod(ss2,1000))
      mij_t=mij(p(0,a1),p(0,a2))
      return
      end

      double precision function delta_phi_t(p,ss1,ss2)
      implicit none
      include "info.inc"
      double precision p(0:4,maxparticles),delta_phi
      integer a1,a2,ss1,ss2
      a1=eventinfo(ss1/1000,mod(ss1,1000))
      a2=eventinfo(ss2/1000,mod(ss2,1000))
      delta_phi_t=delta_phi(p(0,a1),p(0,a2))
      return
      end

      double precision function delta_eta_t(p,ss1,ss2)
      implicit none
      include "info.inc"
      double precision p(0:4,maxparticles),delta_eta
      integer a1,a2,ss1,ss2
      a1=eventinfo(ss1/1000,mod(ss1,1000))
      a2=eventinfo(ss2/1000,mod(ss2,1000))
      delta_eta_t=delta_eta(p(0,a1),p(0,a2))
      return
      end

      double precision function kTij_t(p,ss1,ss2)
      implicit none
      include "info.inc"
      double precision p(0:4,maxparticles),kTij
      integer a1,a2,ss1,ss2
      a1=eventinfo(ss1/1000,mod(ss1,1000))
      a2=eventinfo(ss2/1000,mod(ss2,1000))
      kTij_t=kTij(p(0,a1),p(0,a2))
      return
      end

      double precision function XY1_t(p,ss1,ss2)
      implicit none
      include "info.inc"
      double precision p(0:4,maxparticles),XY1
      integer a1,a2,ss1,ss2
      a1=eventinfo(ss1/1000,mod(ss1,1000))
      a2=eventinfo(ss2/1000,mod(ss2,1000))
      XY1_t=XY1(p(0,a1),p(0,a2))
      return
      end

      double precision function XY2_t(p,ss1,ss2)
      implicit none
      include "info.inc"
      double precision p(0:4,maxparticles),XY2
      integer a1,a2,ss1,ss2
      a1=eventinfo(ss1/1000,mod(ss1,1000))
      a2=eventinfo(ss2/1000,mod(ss2,1000))
      XY2_t=XY2(p(0,a1),p(0,a2))
      return
      end

      double precision function XY3_t(p,ss1,ss2)
      implicit none
      include "info.inc"
      double precision p(0:4,maxparticles),XY3
      integer a1,a2,ss1,ss2
      a1=eventinfo(ss1/1000,mod(ss1,1000))
      a2=eventinfo(ss2/1000,mod(ss2,1000))
      XY3_t=XY3(p(0,a1),p(0,a2))
      return
      end

      double precision function XYZ1_t(p,ss1,ss2,ss3)
      implicit none
      include "info.inc"
      double precision p(0:4,maxparticles),XYZ1
      integer a1,a2,a3,ss1,ss2,ss3
      a1=eventinfo(ss1/1000,mod(ss1,1000))
      a2=eventinfo(ss2/1000,mod(ss2,1000))
      a3=eventinfo(ss3/1000,mod(ss3,1000))
      XYZ1_t=XYZ1(p(0,a1),p(0,a2),p(0,a3))
      return
      end

      double precision function XYZ2_t(p,ss1,ss2,ss3)
      implicit none
      include "info.inc"
      double precision p(0:4,maxparticles),XYZ2
      integer a1,a2,a3,ss1,ss2,ss3
      a1=eventinfo(ss1/1000,mod(ss1,1000))
      a2=eventinfo(ss2/1000,mod(ss2,1000))
      a3=eventinfo(ss3/1000,mod(ss3,1000))
      XYZ2_t=XYZ2(p(0,a1),p(0,a2),p(0,a3))
      return
      end

      double precision function XYZ3_t(p,ss1,ss2,ss3)
      implicit none
      include "info.inc"
      double precision p(0:4,maxparticles),XYZ3
      integer a1,a2,a3,ss1,ss2,ss3
      a1=eventinfo(ss1/1000,mod(ss1,1000))
      a2=eventinfo(ss2/1000,mod(ss2,1000))
      a3=eventinfo(ss3/1000,mod(ss3,1000))
      XYZ3_t=XYZ3(p(0,a1),p(0,a2),p(0,a3))
      return
      end

      double precision function XYZA1_t(p,ss1,ss2,ss3,ss4)
      implicit none
      include "info.inc"
      double precision p(0:4,maxparticles),XYZA1
      integer a1,a2,a3,a4,ss1,ss2,ss3,ss4
      a1=eventinfo(ss1/1000,mod(ss1,1000))
      a2=eventinfo(ss2/1000,mod(ss2,1000))
      a3=eventinfo(ss3/1000,mod(ss3,1000))
      a4=eventinfo(ss4/1000,mod(ss4,1000))
      XYZA1_t=XYZA1(p(0,a1),p(0,a2),p(0,a3),p(0,a4))
      return
      end

      double precision function XYZA2_t(p,ss1,ss2,ss3,ss4)
      implicit none
      include "info.inc"
      double precision p(0:4,maxparticles),XYZA2
      integer a1,a2,a3,a4,ss1,ss2,ss3,ss4
      a1=eventinfo(ss1/1000,mod(ss1,1000))
      a2=eventinfo(ss2/1000,mod(ss2,1000))
      a3=eventinfo(ss3/1000,mod(ss3,1000))
      a4=eventinfo(ss4/1000,mod(ss4,1000))
      XYZA2_t=XYZA2(p(0,a1),p(0,a2),p(0,a3),p(0,a4))
      return
      end

      double precision function XYZA3_t(p,ss1,ss2,ss3,ss4)
      implicit none
      include "info.inc"
      double precision p(0:4,maxparticles),XYZA3
      integer a1,a2,a3,a4,ss1,ss2,ss3,ss4
      a1=eventinfo(ss1/1000,mod(ss1,1000))
      a2=eventinfo(ss2/1000,mod(ss2,1000))
      a3=eventinfo(ss3/1000,mod(ss3,1000))
      a4=eventinfo(ss4/1000,mod(ss4,1000))
      XYZA3_t=XYZA3(p(0,a1),p(0,a2),p(0,a3),p(0,a4))
      return
      end


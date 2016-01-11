
      SUBROUTINE SETCUTS
C**************************************************************************
C     SET THE CUTS 
c     Cuts are set for particles in classes. Particles in a class are
c     ordered in pt. 
C**************************************************************************
      IMPLICIT NONE
      include 'info.inc'
c
c     Constants
c
      double precision zero
      parameter       (ZERO = 0d0)
      real*8 Pi
      parameter( Pi = 3.14159265358979323846d0 )
c
c     LOCAL
c
      integer i,j,m,n,k,l,r,s
      character*140 buffer
      character*15 cutname
      integer l1
      integer cclass,cparticle,cclass2,cparticle2
      integer cclass3,cparticle3,cclass4,cparticle4
      double precision cvalue

C
C-----
C  BEGIN CODE
C-----

      write(*,*) 'Setting up acceptance cuts'


      call set_pointers

C
C     First set all inactive
C
      do i=1,max1
         etmin_active(i)  =.false.
         etmax_active(i)  =.false.
         ptmin_active(i)  =.false.
         ptmax_active(i)  =.false.
         emin_active(i)  =.false.
         emax_active(i)  =.false.
         mommin_active(i)  =.false.
         mommax_active(i)  =.false.
         etamax_active(i) =.false.
         etamin_active(i) =.false.
         ymax_active(i) =.false.
         ymin_active(i) =.false.
         costhmax_active(i) =.false.
         costhmin_active(i) =.false.
         phimax_active(i) =.false.
         phimin_active(i) =.false.
         X1min_active(i)  =.false.
         X1max_active(i)  =.false.
         X2min_active(i)  =.false.
         X2max_active(i)  =.false.
         X3min_active(i)  =.false.
         X3max_active(i)  =.false.
      enddo
      do i=1,max2
         dRijmin_active(i)      =.false.
         dRijmax_active(i)      =.false.
         mijmin_active(i)       =.false.
         mijmax_active(i)       =.false.
         delta_phimin_active(i) =.false.
         delta_phimax_active(i) =.false.
         delta_etamin_active(i) =.false.
         delta_etamax_active(i) =.false.
         kTijmin_active(i)      =.false.
         kTijmax_active(i)      =.false.
         XY1min_active(i)       =.false.
         XY1max_active(i)       =.false.
         XY2min_active(i)       =.false.
         XY2max_active(i)       =.false.
         XY3min_active(i)       =.false.
         XY3max_active(i)       =.false.
      enddo
      do i=1,max3
         XYZ1min_active(i) =.false.
         XYZ1max_active(i) =.false.
         XYZ2min_active(i) =.false.
         XYZ2max_active(i) =.false.
         XYZ3min_active(i) =.false.
         XYZ3max_active(i) =.false.
      enddo
      do i=1,max4
         XYZA1min_active(i) =.false.
         XYZA1max_active(i) =.false.
         XYZA2min_active(i) =.false.
         XYZA2max_active(i) =.false.
         XYZA3min_active(i) =.false.
         XYZA3max_active(i) =.false.
      enddo

c
c     missing ET
c
      etmissmin=0d0
      etmissmax=0d0


C
C     Start reading cuts from sa_card.dat
C
      rewind(97)
 21   read(97,'(a)',end=20,err=20) buffer
      if (buffer(1:13).eq. "# Begin Cuts ") then
 25      read(97,'(a)') buffer
         if (buffer(1:11).eq. "# End Cuts ") goto 23 ! End of Cuts
         l1=50
         if(index(buffer,"#").ne.0) l1=index(buffer,"#")-1 ! ignore comments
C Set cuts as read from file         
         if(buffer(1:6).eq.'ptmin ')then
            read(buffer,*) cutname, cclass, cparticle,cvalue
            call fill1(cclass,cparticle,cvalue,ptmin,ptmin_active)
         elseif(buffer(1:6).eq.'ptmax ')then
            read(buffer,*) cutname, cclass, cparticle,cvalue
            call fill1(cclass,cparticle,cvalue,ptmax,ptmax_active)
         elseif (buffer(1:6).eq.'etmin ') then
            read(buffer,*) cutname, cclass, cparticle,cvalue
            call fill1(cclass,cparticle,cvalue,etmin,etmin_active)
         elseif(buffer(1:6).eq.'etmax ')then
            read(buffer,*)  cutname, cclass, cparticle,cvalue
            call fill1(cclass,cparticle,cvalue,etmax,etmax_active)
         elseif (buffer(1:5).eq.'emin ') then
            read(buffer,*) cutname, cclass, cparticle,cvalue
            call fill1(cclass,cparticle,cvalue,emin,emin_active)
         elseif(buffer(1:5).eq.'emax ')then
            read(buffer,*)  cutname, cclass, cparticle,cvalue
            call fill1(cclass,cparticle,cvalue,emax,emax_active)
         elseif (buffer(1:7).eq.'mommin ') then
            read(buffer,*) cutname, cclass, cparticle,cvalue
            call fill1(cclass,cparticle,cvalue,mommin,mommin_active)
         elseif(buffer(1:7).eq.'mommax ')then
            read(buffer,*)  cutname, cclass, cparticle,cvalue
            call fill1(cclass,cparticle,cvalue,mommax,mommax_active)
         elseif(buffer(1:10).eq.'etmissmin ')then
            read(buffer,*) cutname, cvalue
            etmissmin=cvalue
         elseif(buffer(1:10).eq.'etmissmax ')then
            read(buffer,*) cutname, cvalue
            etmissmax=cvalue
         elseif(buffer(1:7).eq.'etamin ')then
            read(buffer,*) cutname, cclass, cparticle,cvalue
            call fill1(cclass,cparticle,cvalue,etamin,etamin_active)
         elseif(buffer(1:7).eq.'etamax ')then
            read(buffer,*) cutname, cclass, cparticle,cvalue
            call fill1(cclass,cparticle,cvalue,etamax,etamax_active)
         elseif(buffer(1:5).eq.'ymin ')then
            read(buffer,*) cutname, cclass, cparticle,cvalue
            call fill1(cclass,cparticle,cvalue,ymin,ymin_active)
         elseif(buffer(1:5).eq.'ymax ')then
            read(buffer,*) cutname, cclass, cparticle,cvalue
            call fill1(cclass,cparticle,cvalue,ymax,ymax_active)
         elseif(buffer(1:9).eq.'costhmin ')then
            read(buffer,*) cutname, cclass, cparticle,cvalue
            call fill1(cclass,cparticle,cvalue,costhmin,costhmin_active)
         elseif(buffer(1:9).eq.'costhmax ')then
            read(buffer,*) cutname, cclass, cparticle,cvalue
            call fill1(cclass,cparticle,cvalue,costhmax,costhmax_active)
         elseif(buffer(1:7).eq.'phimin ')then
            read(buffer,*) cutname, cclass, cparticle,cvalue
            call fill1(cclass,cparticle,cvalue,phimin,phimin_active)
         elseif(buffer(1:7).eq.'phimax ')then
            read(buffer,*) cutname, cclass, cparticle,cvalue
            call fill1(cclass,cparticle,cvalue,phimax,phimax_active)
         elseif(buffer(1:6).eq.'X1min ')then
            read(buffer,*)  cutname, cclass, cparticle,cvalue
            call fill1(cclass,cparticle,cvalue,X1min,X1min_active)
         elseif(buffer(1:6).eq.'X1max ')then
            read(buffer,*)  cutname, cclass, cparticle,cvalue
            call fill1(cclass,cparticle,cvalue,X1max,X1max_active)
         elseif(buffer(1:6).eq.'X2min ')then
            read(buffer,*)  cutname, cclass, cparticle,cvalue
            call fill1(cclass,cparticle,cvalue,X2min,X2min_active)
         elseif(buffer(1:6).eq.'X2max ')then
            read(buffer,*)  cutname, cclass, cparticle,cvalue
            call fill1(cclass,cparticle,cvalue,X2max,X2max_active)
         elseif(buffer(1:6).eq.'X3min ')then
            read(buffer,*)  cutname, cclass, cparticle,cvalue
            call fill1(cclass,cparticle,cvalue,X3min,X3min_active)
         elseif(buffer(1:6).eq.'X3max ')then
            read(buffer,*)  cutname, cclass, cparticle,cvalue
            call fill1(cclass,cparticle,cvalue,X3max,X3max_active)
         elseif(buffer(1:8).eq.'dRijmin ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cvalue
            call fill2(cclass,cparticle,cclass2,cparticle2,cvalue,dRijmin,dRijmin_active)
         elseif(buffer(1:8).eq.'dRijmax ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cvalue
            call fill2(cclass,cparticle,cclass2,cparticle2,cvalue,dRijmax,dRijmax_active)
         elseif(buffer(1:7).eq.'mijmin ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cvalue
            call fill2(cclass,cparticle,cclass2,cparticle2,cvalue,mijmin,mijmin_active)
         elseif(buffer(1:7).eq.'mijmax ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cvalue
            call fill2(cclass,cparticle,cclass2,cparticle2,cvalue,mijmax,mijmax_active)
         elseif(buffer(1:13).eq.'delta_phimin ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cvalue
            call fill2(cclass,cparticle,cclass2,cparticle2,cvalue,delta_phimin,delta_phimin_active)
         elseif(buffer(1:13).eq.'delta_phimax ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cvalue
            call fill2(cclass,cparticle,cclass2,cparticle2,cvalue,delta_phimax,delta_phimax_active)
         elseif(buffer(1:13).eq.'delta_etamin ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cvalue
            call fill2(cclass,cparticle,cclass2,cparticle2,cvalue,delta_etamin,delta_etamin_active)
         elseif(buffer(1:13).eq.'delta_etamax ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cvalue
            call fill2(cclass,cparticle,cclass2,cparticle2,cvalue,delta_etamax,delta_etamax_active)
         elseif(buffer(1:8).eq.'kTijmin ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cvalue
            call fill2(cclass,cparticle,cclass2,cparticle2,cvalue,kTijmin,kTijmin_active)
         elseif(buffer(1:8).eq.'kTijmax ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cvalue
            call fill2(cclass,cparticle,cclass2,cparticle2,cvalue,kTijmax,kTijmax_active)
         elseif(buffer(1:7).eq.'XY1min ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cvalue
            call fill2(cclass,cparticle,cclass2,cparticle2,cvalue,XY1min,XY1min_active)
         elseif(buffer(1:7).eq.'XY1max ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cvalue
            call fill2(cclass,cparticle,cclass2,cparticle2,cvalue,XY1max,XY1max_active)
         elseif(buffer(1:7).eq.'XY2min ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cvalue
            call fill2(cclass,cparticle,cclass2,cparticle2,cvalue,XY2min,XY2min_active)
         elseif(buffer(1:7).eq.'XY2max ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cvalue
            call fill2(cclass,cparticle,cclass2,cparticle2,cvalue,XY2max,XY2max_active)
         elseif(buffer(1:7).eq.'XY3min ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cvalue
            call fill2(cclass,cparticle,cclass2,cparticle2,cvalue,XY3min,XY3min_active)
         elseif(buffer(1:7).eq.'XY3max ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cvalue
            call fill2(cclass,cparticle,cclass2,cparticle2,cvalue,XY3max,XY3max_active)
         elseif(buffer(1:8).eq.'XYZ1min ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cclass3,cparticle3,cvalue
            call fill3(cclass,cparticle,cclass2,cparticle2,cclass3,cparticle3,cvalue,XYZ1min,XYZ1min_active)
         elseif(buffer(1:8).eq.'XYZ1max ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cclass3,cparticle3,cvalue
            call fill3(cclass,cparticle,cclass2,cparticle2,cclass3,cparticle3,cvalue,XYZ1max,XYZ1max_active)
         elseif(buffer(1:8).eq.'XYZ2min ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cclass3,cparticle3,cvalue
            call fill3(cclass,cparticle,cclass2,cparticle2,cclass3,cparticle3,cvalue,XYZ2min,XYZ2min_active)
         elseif(buffer(1:8).eq.'XYZ2max ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cclass3,cparticle3,cvalue
            call fill3(cclass,cparticle,cclass2,cparticle2,cclass3,cparticle3,cvalue,XYZ2max,XYZ2max_active)
         elseif(buffer(1:8).eq.'XYZ3min ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cclass3,cparticle3,cvalue
            call fill3(cclass,cparticle,cclass2,cparticle2,cclass3,cparticle3,cvalue,XYZ3min,XYZ3min_active)
         elseif(buffer(1:8).eq.'XYZ3max ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cclass3,cparticle3,cvalue
            call fill3(cclass,cparticle,cclass2,cparticle2,cclass3,cparticle3,cvalue,XYZ3max,XYZ3max_active)
         elseif(buffer(1:9).eq.'XYZA1min ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cclass3,cparticle3,cclass4,cparticle4,cvalue
            call fill4(cclass ,cparticle ,cclass2,cparticle2,cclass3,cparticle3,cclass4,cparticle4,
     &                 cvalue,XYZA1min,XYZA1min_active)
         elseif(buffer(1:9).eq.'XYZA1max ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cclass3,cparticle3,cclass4,cparticle4,cvalue
            call fill4(cclass ,cparticle ,cclass2,cparticle2,cclass3,cparticle3,cclass4,cparticle4,
     &                 cvalue,XYZA1max,XYZA1max_active)
         elseif(buffer(1:9).eq.'XYZA2min ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cclass3,cparticle3,cclass4,cparticle4,cvalue
            call fill4(cclass ,cparticle ,cclass2,cparticle2,cclass3,cparticle3,cclass4,cparticle4,
     &                 cvalue,XYZA2min,XYZA2min_active)
         elseif(buffer(1:9).eq.'XYZA2max ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cclass3,cparticle3,cclass4,cparticle4,cvalue
            call fill4(cclass ,cparticle ,cclass2,cparticle2,cclass3,cparticle3,cclass4,cparticle4,
     &                 cvalue,XYZA2max,XYZA2max_active)
         elseif(buffer(1:9).eq.'XYZA3min ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cclass3,cparticle3,cclass4,cparticle4,cvalue
            call fill4(cclass ,cparticle ,cclass2,cparticle2,cclass3,cparticle3,cclass4,cparticle4,
     &                 cvalue,XYZA3min,XYZA3min_active)
         elseif(buffer(1:9).eq.'XYZA3max ')then
            read(buffer,*) cutname, cclass, cparticle,cclass2,cparticle2,cclass3,cparticle3,cclass4,cparticle4,cvalue
            call fill4(cclass ,cparticle ,cclass2,cparticle2,cclass3,cparticle3,cclass4,cparticle4,
     &                 cvalue,XYZA3max,XYZA3max_active)
         endif
         goto 25
 23      continue               ! End of cuts
      endif
      goto 21
 20   continue

      call write_cuts

      write (*,*) 'Cuts set up'
      write (*,*) 'SEE plots.log FOR MORE INFO'

      RETURN
      END


      LOGICAL FUNCTION PASSCUTS(P,nexternal)
C**************************************************************************
C     INPUT:  MOMENTA
C     OUTPUT: TRUE IF EVENTS PASSES ALL CUTS LISTED
c
c     Use eventinfo(i,j) to specify particles.
c     i              class  (i=0 for missing Et class)
c     j              particle in class (ordered from large to small pt)
c     eventinfo(i,j) is equal to 0 if particle does not exists
C**************************************************************************
      IMPLICIT NONE
      include "info.inc"
C
C     ARGUMENTS
C
      REAL*8 P(0:4,Maxparticles)
      integer nexternal
C
C     LOCAL
C
      LOGICAL FIRSTTIME
      data firsttime/.true./
      integer i,j,m,n,k,l,r,s
      double precision pmiss(0:3)
C
C     EXTERNAL
C
      REAL*8 ET
      real*8 et_t,eta_t,pt_t,X1_t,X2_t,X3_t,dRij_t,mij_t,delta_phi_t,XY1_t,XY2_t,XY3_t
      real*8 e_t,mom_t,y_t,costh_t,phi_t,delta_eta_t,kTij_t
      real*8 XYZ1_t,XYZ2_t,XYZ3_t,XYZA1_t,XYZA2_t,XYZA3_t
      logical cut_bw
      external cut_bw
C
C     GLOBAL
C
      integer maxclass,maxevparticles,partoclass(-300:300,0:MaxParticles)
      integer maxevpart(0:maxclasses),maxpart(0:maxclasses)
      common/classes/maxclass,maxevparticles,partoclass,maxevpart
C-----
C  BEGIN CODE
C-----

       PASSCUTS=.TRUE.           !EVENT IS OK UNLESS OTHERWISE CHANGED
c       RETURN                   ! uncomment if you do not want to put cuts
c
c

c
c Cuts for one particle
c
       do i=1,max1
c     pt min cut 
          if (ptmin_active(i))then
             if (pt_t(p,s1(i)).lt.ptmin(i).or.eventinfo(s1(i)/1000,mod(s1(i),1000)).eq.0) then
                passcuts=.false.
                return
             endif
          endif
c     pt max cut
          if (ptmax_active(i).and.eventinfo(s1(i)/1000,mod(s1(i),1000)).ne.0)then
             if (pt_t(p,s1(i)) .gt. ptmax(i)) then
                passcuts=.false.
                return
             endif
          endif
c     Et min cut 
          if (etmin_active(i))then
             if (et_t(p,s1(i)).lt.etmin(i).or.eventinfo(s1(i)/1000,mod(s1(i),1000)).eq.0) then
                passcuts=.false.
                return
             endif
          endif
c     Et max
          if (etmax_active(i).and.eventinfo(s1(i)/1000,mod(s1(i),1000)).ne.0)then
             if (et_t(p,s1(i)) .gt. etmax(i)) then
                passcuts=.false.
                return
             endif
          endif
c     E min cut 
          if (emin_active(i))then
             if (e_t(p,s1(i)).lt.emin(i).or.eventinfo(s1(i)/1000,mod(s1(i),1000)).eq.0) then
                passcuts=.false.
                return
             endif
          endif
c     E max
          if (emax_active(i).and.eventinfo(s1(i)/1000,mod(s1(i),1000)).ne.0)then
             if (e_t(p,s1(i)) .gt. emax(i)) then
                passcuts=.false.
                return
             endif
          endif
c     mom min cut 
          if (mommin_active(i))then
             if (mom_t(p,s1(i)).lt.mommin(i).or.eventinfo(s1(i)/1000,mod(s1(i),1000)).eq.0) then
                passcuts=.false.
                return
             endif
          endif
c     mom max
          if (mommax_active(i).and.eventinfo(s1(i)/1000,mod(s1(i),1000)).ne.0)then
             if (mom_t(p,s1(i)) .gt. mommax(i)) then
                passcuts=.false.
                return
             endif
          endif
c     eta min cut
          if (etamin_active(i))then
             if (abs(eta_t(p,s1)).lt.etamin(i).or.eventinfo(s1(i)/1000,mod(s1(i),1000)).eq.0)then
                passcuts=.false.
                return
             endif
          endif
c     eta max cut
          if (etamax_active(i).and.eventinfo(s1(i)/1000,mod(s1(i),1000)).ne.0)then
             if (abs(eta_t(p,s1(i))) .gt. etamax(i)) then
                passcuts=.false.
                return
             endif
          endif
c     y min cut
          if (ymin_active(i))then
             if (abs(y_t(p,s1)).lt.ymin(i).or.eventinfo(s1(i)/1000,mod(s1(i),1000)).eq.0)then
                passcuts=.false.
                return
             endif
          endif
c     y max cut
          if (ymax_active(i).and.eventinfo(s1(i)/1000,mod(s1(i),1000)).ne.0)then
             if (abs(y_t(p,s1(i))) .gt. ymax(i)) then
                passcuts=.false.
                return
             endif
          endif
c     costh min cut
          if (costhmin_active(i))then
             if (costh_t(p,s1).lt.costhmin(i).or.eventinfo(s1(i)/1000,mod(s1(i),1000)).eq.0)then
                passcuts=.false.
                return
             endif
          endif
c     costh max cut
          if (costhmax_active(i).and.eventinfo(s1(i)/1000,mod(s1(i),1000)).ne.0)then
             if (costh_t(p,s1(i)) .gt. costhmax(i)) then
                passcuts=.false.
                return
             endif
          endif
c     phi min cut
          if (phimin_active(i))then
             if (phi_t(p,s1).lt.phimin(i).or.eventinfo(s1(i)/1000,mod(s1(i),1000)).eq.0)then
                passcuts=.false.
                return
             endif
          endif
c     phi max cut
          if (phimax_active(i).and.eventinfo(s1(i)/1000,mod(s1(i),1000)).ne.0)then
             if (phi_t(p,s1(i)) .gt. phimax(i)) then
                passcuts=.false.
                return
             endif
          endif
c     X1 min cut 
          if (X1min_active(i))then
             if (X1_t(p,s1(i)).lt.X1min(i).or.eventinfo(s1(i)/1000,mod(s1(i),1000)).eq.0) then
                passcuts=.false.
                return
             endif
          endif
c     X1 max cut
          if (X1max_active(i).and.eventinfo(s1(i)/1000,mod(s1(i),1000)).ne.0)then
             if (X1_t(p,s1(i)) .gt. X1max(i)) then
                passcuts=.false.
                return
             endif
          endif
c     X2 min cut 
          if (X2min_active(i))then
             if (X2_t(p,s1(i)).lt.X2min(i).or.eventinfo(s1(i)/1000,mod(s1(i),1000)).eq.0) then
                passcuts=.false.
                return
             endif
          endif
c     X2 max cut
          if (X2max_active(i).and.eventinfo(s1(i)/1000,mod(s1(i),1000)).ne.0)then
             if (X2_t(p,s1(i)) .gt. X2max(i)) then
                passcuts=.false.
                return
             endif
          endif
c     X3 min cut 
          if (X3min_active(i))then
             if (X3_t(p,s1(i)).lt.X3min(i).or.eventinfo(s1(i)/1000,mod(s1(i),1000)).eq.0) then
                passcuts=.false.
                return
             endif
          endif
c     X3 max cut
          if (X3max_active(i).and.eventinfo(s1(i)/1000,mod(s1(i),1000)).ne.0)then
             if (X3_t(p,s1(i)) .gt. X3max(i)) then
                passcuts=.false.
                return
             endif
          endif
       enddo
c
c Cuts between two particles
c
       do i=1,max2
c     dRij min cut
          if (dRijmin_active(i)
     &         .and.eventinfo(s2(1,i)/1000,mod(s2(1,i),1000)).ne.0
     &         .and.eventinfo(s2(2,i)/1000,mod(s2(2,i),1000)).ne.0 )then
             if (dRij_t(p,s2(1,i),s2(2,i)).lt. dRijmin(i)) then
                passcuts=.false.
                return
             endif
          endif
c     dRij max cut
          if(dRijmax_active(i))then
             if (dRij_t(p,s2(1,i),s2(2,i)).gt.dRijmax(i)
     &            .or.eventinfo(s2(1,i)/1000,mod(s2(1,i),1000)).eq.0
     $            .or.eventinfo(s2(2,i)/1000,mod(s2(2,i),1000)).eq.0)then 
                passcuts=.false.
                return
             endif
          endif
c     mij min cut
          if (mijmin_active(i))then
             if (mij_t(p,s2(1,i),s2(2,i)).lt. mijmin(i)
     &            .or.eventinfo(s2(1,i)/1000,mod(s2(1,i),1000)).eq.0
     &            .or.eventinfo(s2(2,i)/1000,mod(s2(2,i),1000)).eq.0) then
                passcuts=.false.
                return
             endif
          endif
c     mij max cut
          if (mijmax_active(i)
     &         .and.eventinfo(s2(1,i)/1000,mod(s2(1,i),1000)).ne.0
     &         .and.eventinfo(s2(2,i)/1000,mod(s2(2,i),1000)).ne.0 )then
             if (mij_t(p,s2(1,i),s2(2,i)).gt.mijmax(i) )then 
                passcuts=.false.
                return
             endif
          endif
c     delta_phi min cut
          if (delta_phimin_active(i))then
             if (delta_phi_t(p,s2(1,i),s2(2,i)).lt. delta_phimin(i)
     &            .or.eventinfo(s2(1,i)/1000,mod(s2(1,i),1000)).eq.0
     &            .or.eventinfo(s2(2,i)/1000,mod(s2(2,i),1000)).eq.0) then
                passcuts=.false.
                return
             endif
          endif
c     delta_phi max cut
          if (delta_phimax_active(i)
     &         .and.eventinfo(s2(1,i)/1000,mod(s2(1,i),1000)).ne.0
     &         .and.eventinfo(s2(2,i)/1000,mod(s2(2,i),1000)).ne.0 )then
             if (delta_phi_t(p,s2(1,i),s2(2,i)).gt.delta_phimax(i) )then 
                passcuts=.false.
                return
             endif
          endif
c     delta_eta min cut
          if (delta_etamin_active(i))then
             if (delta_eta_t(p,s2(1,i),s2(2,i)).lt. delta_etamin(i)
     &            .or.eventinfo(s2(1,i)/1000,mod(s2(1,i),1000)).eq.0
     &            .or.eventinfo(s2(2,i)/1000,mod(s2(2,i),1000)).eq.0) then
                passcuts=.false.
                return
             endif
          endif
c     delta_eta max cut
          if (delta_etamax_active(i)
     &         .and.eventinfo(s2(1,i)/1000,mod(s2(1,i),1000)).ne.0
     &         .and.eventinfo(s2(2,i)/1000,mod(s2(2,i),1000)).ne.0 )then
             if (delta_eta_t(p,s2(1,i),s2(2,i)).gt.delta_etamax(i) )then 
                passcuts=.false.
                return
             endif
          endif
c     kTij min cut
          if (kTijmin_active(i)
     &         .and.eventinfo(s2(1,i)/1000,mod(s2(1,i),1000)).ne.0
     &         .and.eventinfo(s2(2,i)/1000,mod(s2(2,i),1000)).ne.0 )then
             if (kTij_t(p,s2(1,i),s2(2,i)).lt. kTijmin(i)) then
                passcuts=.false.
                return
             endif
          endif
c     kTij max cut
          if(kTijmax_active(i))then
             if (kTij_t(p,s2(1,i),s2(2,i)).gt.kTijmax(i)
     &            .or.eventinfo(s2(1,i)/1000,mod(s2(1,i),1000)).eq.0
     $            .or.eventinfo(s2(2,i)/1000,mod(s2(2,i),1000)).eq.0)then 
                passcuts=.false.
                return
             endif
          endif
c     XY1 min cut
          if (XY1min_active(i))then
             if (XY1_t(p,s2(1,i),s2(2,i)).lt. XY1min(i)
     &            .or.eventinfo(s2(1,i)/1000,mod(s2(1,i),1000)).eq.0
     &            .or.eventinfo(s2(2,i)/1000,mod(s2(2,i),1000)).eq.0) then
                passcuts=.false.
                return
             endif
          endif
c     XY1 max cut
          if (XY1max_active(i)
     &         .and.eventinfo(s2(1,i)/1000,mod(s2(1,i),1000)).ne.0
     &         .and.eventinfo(s2(2,i)/1000,mod(s2(2,i),1000)).ne.0 )then
             if (XY1_t(p,s2(1,i),s2(2,i)).gt.XY1max(i) )then 
                passcuts=.false.
                return
             endif
          endif
c     XY2 min cut
          if (XY2min_active(i))then
             if (XY2_t(p,s2(1,i),s2(2,i)).lt. XY2min(i)
     &            .or.eventinfo(s2(1,i)/1000,mod(s2(1,i),1000)).eq.0
     &            .or.eventinfo(s2(2,i)/1000,mod(s2(2,i),1000)).eq.0) then
                passcuts=.false.
                return
             endif
          endif
c     XY2 max cut
          if (XY2max_active(i)
     &         .and.eventinfo(s2(1,i)/1000,mod(s2(1,i),1000)).ne.0
     &         .and.eventinfo(s2(2,i)/1000,mod(s2(2,i),1000)).ne.0 )then
             if (XY2_t(p,s2(1,i),s2(2,i)).gt.XY2max(i) )then 
                passcuts=.false.
                return
             endif
          endif
c     XY3 min cut
          if (XY3min_active(i))then
             if (XY3_t(p,s2(1,i),s2(2,i)).lt. XY3min(i)
     &            .or.eventinfo(s2(1,i)/1000,mod(s2(1,i),1000)).eq.0
     &            .or.eventinfo(s2(2,i)/1000,mod(s2(2,i),1000)).eq.0) then
                passcuts=.false.
                return
             endif
          endif
c     XY3 max cut
          if (XY3max_active(i)
     &         .and.eventinfo(s2(1,i)/1000,mod(s2(1,i),1000)).ne.0
     &         .and.eventinfo(s2(2,i)/1000,mod(s2(2,i),1000)).ne.0 )then
             if (XY3_t(p,s2(1,i),s2(2,i)).gt.XY3max(i) )then 
                passcuts=.false.
                return
             endif
          endif
       enddo
c
c Cuts between three particles
c
       do i=1,max3
c     XYZ1 min cut
          if (XYZ1min_active(i))then
             if (XYZ1_t(p,s3(1,i),s3(2,i),s3(3,i)).lt. XYZ1min(i)
     &            .or.eventinfo(s3(1,i)/1000,mod(s3(1,i),1000)).eq.0
     &            .or.eventinfo(s3(2,i)/1000,mod(s3(2,i),1000)).eq.0
     &            .or.eventinfo(s3(3,i)/1000,mod(s3(3,i),1000)).eq.0) then
                passcuts=.false.
                return
             endif
          endif
c     XYZ1 max cut
          if (XYZ1max_active(i)
     &         .and.eventinfo(s3(1,i)/1000,mod(s3(1,i),1000)).ne.0
     &         .and.eventinfo(s3(2,i)/1000,mod(s3(2,i),1000)).ne.0
     &         .and.eventinfo(s3(3,i)/1000,mod(s3(3,i),1000)).ne.0 )then
             if (XYZ1_t(p,s3(1,i),s3(2,i),s3(3,i)).gt.XYZ1max(i) )then 
                passcuts=.false.
                return
             endif
          endif
c     XYZ2 min cut
          if (XYZ2min_active(i))then
             if (XYZ2_t(p,s3(1,i),s3(2,i),s3(3,i)).lt. XYZ2min(i)
     &            .or.eventinfo(s3(1,i)/1000,mod(s3(1,i),1000)).eq.0
     &            .or.eventinfo(s3(2,i)/1000,mod(s3(2,i),1000)).eq.0
     &            .or.eventinfo(s3(3,i)/1000,mod(s3(3,i),1000)).eq.0) then
                passcuts=.false.
                return
             endif
          endif
c     XYZ2 max cut
          if (XYZ2max_active(i)
     &         .and.eventinfo(s3(1,i)/1000,mod(s3(1,i),1000)).ne.0
     &         .and.eventinfo(s3(2,i)/1000,mod(s3(2,i),1000)).ne.0
     &         .and.eventinfo(s3(3,i)/1000,mod(s3(3,i),1000)).ne.0 )then
             if (XYZ2_t(p,s3(1,i),s3(2,i),s3(3,i)).gt.XYZ2max(i) )then 
                passcuts=.false.
                return
             endif
          endif
c     XY31 min cut
          if (XYZ3min_active(i))then
             if (XYZ3_t(p,s3(1,i),s3(2,i),s3(3,i)).lt. XYZ3min(i)
     &            .or.eventinfo(s3(1,i)/1000,mod(s3(1,i),1000)).eq.0
     &            .or.eventinfo(s3(2,i)/1000,mod(s3(2,i),1000)).eq.0
     &            .or.eventinfo(s3(3,i)/1000,mod(s3(3,i),1000)).eq.0) then
                passcuts=.false.
                return
             endif
          endif
c     XYZ3 max cut
          if (XYZ3max_active(i)
     &         .and.eventinfo(s3(1,i)/1000,mod(s3(1,i),1000)).ne.0
     &         .and.eventinfo(s3(2,i)/1000,mod(s3(2,i),1000)).ne.0
     &         .and.eventinfo(s3(3,i)/1000,mod(s3(3,i),1000)).ne.0 )then
             if (XYZ3_t(p,s3(1,i),s3(2,i),s3(3,i)).gt.XYZ3max(i) )then 
                passcuts=.false.
                return
             endif
          endif
       enddo
c
c Cuts between four particles
c
       do i=1,max4
c     XYZA1 min cut
          if (XYZA1min_active(i))then
             if (XYZA1_t(p,s4(1,i),s4(2,i),s4(3,i),s4(4,i)).lt. XYZA1min(i)
     &            .or.eventinfo(s4(1,i)/1000,mod(s4(1,i),1000)).eq.0
     &            .or.eventinfo(s4(2,i)/1000,mod(s4(2,i),1000)).eq.0
     &            .or.eventinfo(s4(3,i)/1000,mod(s4(3,i),1000)).eq.0
     &            .or.eventinfo(s4(4,i)/1000,mod(s4(4,i),1000)).eq.0) then
                passcuts=.false.
                return
             endif
          endif
c     XYZA1 max cut
          if (XYZA1max_active(i)
     &         .and.eventinfo(s4(1,i)/1000,mod(s4(1,i),1000)).ne.0
     &         .and.eventinfo(s4(2,i)/1000,mod(s4(2,i),1000)).ne.0
     &         .and.eventinfo(s4(3,i)/1000,mod(s4(3,i),1000)).ne.0
     &         .and.eventinfo(s4(4,i)/1000,mod(s4(4,i),1000)).ne.0 )then
             if (XYZA1_t(p,s4(1,i),s4(2,i),s4(3,i),s4(4,i)).gt.XYZA1max(i) )then 
                passcuts=.false.
                return
             endif
          endif
c     XYZA2 min cut
          if (XYZA2min_active(i))then
             if (XYZA2_t(p,s4(1,i),s4(2,i),s4(3,i),s4(4,i)).lt. XYZA2min(i)
     &            .or.eventinfo(s4(1,i)/1000,mod(s4(1,i),1000)).eq.0
     &            .or.eventinfo(s4(2,i)/1000,mod(s4(2,i),1000)).eq.0
     &            .or.eventinfo(s4(3,i)/1000,mod(s4(3,i),1000)).eq.0
     &            .or.eventinfo(s4(4,i)/1000,mod(s4(4,i),1000)).eq.0) then
                passcuts=.false.
                return
             endif
          endif
c     XYZA2 max cut
          if (XYZA2max_active(i)
     &         .and.eventinfo(s4(1,i)/1000,mod(s4(1,i),1000)).ne.0
     &         .and.eventinfo(s4(2,i)/1000,mod(s4(2,i),1000)).ne.0
     &         .and.eventinfo(s4(3,i)/1000,mod(s4(3,i),1000)).ne.0
     &         .and.eventinfo(s4(4,i)/1000,mod(s4(4,i),1000)).ne.0 )then
             if (XYZA2_t(p,s4(1,i),s4(2,i),s4(3,i),s4(4,i)).gt.XYZA2max(i) )then 
                passcuts=.false.
                return
             endif
          endif
c     XYZA3 min cut
          if (XYZA3min_active(i))then
             if (XYZA3_t(p,s4(1,i),s4(2,i),s4(3,i),s4(4,i)).lt. XYZA3min(i)
     &            .or.eventinfo(s4(1,i)/1000,mod(s4(1,i),1000)).eq.0
     &            .or.eventinfo(s4(2,i)/1000,mod(s4(2,i),1000)).eq.0
     &            .or.eventinfo(s4(3,i)/1000,mod(s4(3,i),1000)).eq.0
     &            .or.eventinfo(s4(4,i)/1000,mod(s4(4,i),1000)).eq.0) then
                passcuts=.false.
                return
             endif
          endif
c     XYZA3 max cut
          if (XYZA3max_active(i)
     &         .and.eventinfo(s4(1,i)/1000,mod(s4(1,i),1000)).ne.0
     &         .and.eventinfo(s4(2,i)/1000,mod(s4(2,i),1000)).ne.0
     &         .and.eventinfo(s4(3,i)/1000,mod(s4(3,i),1000)).ne.0
     &         .and.eventinfo(s4(4,i)/1000,mod(s4(4,i),1000)).ne.0 )then
             if (XYZA3_t(p,s4(1,i),s4(2,i),s4(3,i),s4(4,i)).gt.XYZA3max(i) )then 
                passcuts=.false.
                return
             endif
          endif
       enddo

c
c     missing Et cut
c
      if (maxevpart(0).ge.1) then  ! If Missing Et class exist it is class 0, check if it is not empty
         do j=0,3
            pmiss(j)=0d0
         enddo
         do k=1,maxevpart(0)
            do j=0,3
               pmiss(j)=pmiss(j)+p(j,eventinfo(0,k))
            enddo
         enddo
c miss Et min cut
         if(et(pmiss).lt.etmissmin) then
            passcuts=.false.
            return
         endif
c miss Et max cut
         if (etmissmax.gt.0d0)then
            if(et(pmiss).gt.etmissmax) then
               passcuts=.false.
               return
            endif
         endif
      endif
c cut also if missing Et class is empty and etmissmin > 0
      if (maxevpart(0).lt.1.and.etmissmin.gt.0d0)then
         passcuts=.false.
         return
      endif


      RETURN
      END






      subroutine set_pointers
      implicit none
      include "info.inc"
      integer maxclass,maxevparticles,partoclass(-300:300,0:MaxParticles)
      integer maxevpart(0:maxclasses),maxpart(0:maxclasses)
      common/classes/maxclass,maxevparticles,partoclass,maxevpart
      integer i,j,k,l,m,totpart

      totpart=0
      do i=0,maxclasses
         if (maxevpart(i).ne.0)then
            do j=1,maxevpart(i)
               totpart=totpart+1
               s1(totpart)=i*1000+j
            enddo
         endif
      enddo
      i=1
      m=1
      do while (i.lt.totpart)
         do j=1,totpart-i
            s2(1,m)=s1(i)
            s2(2,m)=s1(i+j)
            m=m+1
         enddo
         i=i+1
      enddo
      i=1
      m=1
      do while (i.lt.totpart)
         do j=1,totpart-i
            do k=1,totpart-i-j
               s3(1,m)=s1(i)
               s3(2,m)=s1(i+j)
               s3(3,m)=s1(i+j+k)
               m=m+1
            enddo
         enddo
         i=i+1
      enddo
      i=1
      m=1
      do while (i.lt.totpart)
         do j=1,totpart-i
            do k=1,totpart-i-j
               do l=1,totpart-i-j-k
                  s4(1,m)=s1(i)
                  s4(2,m)=s1(i+j)
                  s4(3,m)=s1(i+j+k)
                  s4(4,m)=s1(i+j+k+l)
                  m=m+1
               enddo
            enddo
         enddo
         i=i+1
      enddo

      return
      end



      subroutine fill1(c1,p1,value,var,active)
      implicit none
      include "info.inc"
      integer c1,p1,ctemp,ptemp,i
      double precision value,var(max1)
      logical active(max1)

      do i=1, max1
         if ( (s1(i)/1000.eq.c1.and.mod(s1(i),1000).le.p1))then
            var(i)=value
            active(i)=.true.
         endif
      enddo

      return 
      end


      subroutine fill2(c1,p1,c2,p2,value,var,active)
      implicit none
      include "info.inc"
      integer c1,c2,p1,p2,ctemp,ptemp,i
      double precision value,var(max2)
      logical active(max2)

      if (c1.gt.c2.or.(c1.eq.c2.and.p1.gt.p2))then
         ctemp=c1
         c1=c2
         c2=ctemp
         ptemp=p1
         p1=p2
         p2=ptemp         
      endif

      do i=1, max2
         if ( (s2(1,i)/1000.eq.c1.and.mod(s2(1,i),1000).le.p1).and.
     &        (s2(2,i)/1000.eq.c2.and.mod(s2(2,i),1000).le.p2)     )then
            var(i)=value
            active(i)=.true.
         endif
      enddo

      return
      end


      subroutine fill3(c1,p1,c2,p2,c3,p3,value,var,active)
      implicit none
      include "info.inc"
      integer c1,c2,c3,p1,p2,p3,ctemp,ptemp,i
      double precision value,var(max3)
      logical active(max3)

      do i=1,2
      if (c1.gt.c2.or.(c1.eq.c2.and.p1.gt.p2))then
         ctemp=c1
         c1=c2
         c2=ctemp
         ptemp=p1
         p1=p2
         p2=ptemp         
      endif
      if (c2.gt.c3.or.(c2.eq.c3.and.p2.gt.p3))then
         ctemp=c2
         c2=c3
         c3=ctemp
         ptemp=p2
         p2=p3
         p3=ptemp         
      endif
      enddo

      do i=1, max3
         if ( (s3(1,i)/1000.eq.c1.and.mod(s3(1,i),1000).le.p1).and.
     &        (s3(2,i)/1000.eq.c2.and.mod(s3(2,i),1000).le.p2).and.
     &        (s3(3,i)/1000.eq.c3.and.mod(s3(3,i),1000).le.p3)     )then
            var(i)=value
            active(i)=.true.
         endif
      enddo

      return
      end


      subroutine fill4(c1,p1,c2,p2,c3,p3,c4,p4,value,var,active)
      implicit none
      include "info.inc"
      integer c1,c2,c3,c4,p1,p2,p3,p4,ctemp,ptemp,i
      double precision value,var(max4)
      logical active(max4)

      do i=1,3
      if (c1.gt.c2.or.(c1.eq.c2.and.p1.gt.p2))then
         ctemp=c1
         c1=c2
         c2=ctemp
         ptemp=p1
         p1=p2
         p2=ptemp         
      endif
      if (c2.gt.c3.or.(c2.eq.c3.and.p2.gt.p3))then
         ctemp=c2
         c2=c3
         c3=ctemp
         ptemp=p2
         p2=p3
         p3=ptemp         
      endif
      if (c3.gt.c4.or.(c3.eq.c4.and.p3.gt.p4))then
         ctemp=c3
         c3=c4
         c4=ctemp
         ptemp=p3
         p3=p4
         p4=ptemp         
      endif
      enddo

      do i=1, max4
         if ( (s4(1,i)/1000.eq.c1.and.mod(s4(1,i),1000).le.p1).and.
     &        (s4(2,i)/1000.eq.c2.and.mod(s4(2,i),1000).le.p2).and.
     &        (s4(3,i)/1000.eq.c3.and.mod(s4(3,i),1000).le.p3).and.
     &        (s4(4,i)/1000.eq.c4.and.mod(s4(4,i),1000).le.p4)     )then
            var(i)=value
            active(i)=.true.
         endif
      enddo


      return
      end



      subroutine write_cuts
      implicit none
      include "info.inc"
      integer i

      write (52,*) ' '
      write (52,*) 'applied cuts:'
      do i=1, max1
         if (ptmin_active(i))then
            call write_cutsnow('ptmin       ',ptmin(i),1,i)
         endif
      enddo
      do i=1, max1
         if (ptmax_active(i))then
            call write_cutsnow('ptmax       ',ptmax(i),1,i)
         endif
      enddo
      do i=1, max1
         if (etmin_active(i))then
            call write_cutsnow('etmin       ',etmin(i),1,i)
         endif
      enddo
      do i=1, max1
         if (etmax_active(i))then
            call write_cutsnow('etmax       ',etmax(i),1,i)
         endif
      enddo
      do i=1, max1
         if (emin_active(i))then
            call write_cutsnow('emin        ',emin(i),1,i)
         endif
      enddo
      do i=1, max1
         if (emax_active(i))then
            call write_cutsnow('emax        ',emax(i),1,i)
         endif
      enddo
      do i=1, max1
         if (mommin_active(i))then
            call write_cutsnow('mommin      ',mommin(i),1,i)
         endif
      enddo
      do i=1, max1
         if (mommax_active(i))then
            call write_cutsnow('mommax      ',mommax(i),1,i)
         endif
      enddo
      do i=1, max1
         if (etamin_active(i))then
            call write_cutsnow('etamin      ',etamin(i),1,i)
         endif
      enddo
      do i=1, max1
         if (etamax_active(i))then
            call write_cutsnow('etamax      ',etamax(i),1,i)
         endif
      enddo
      do i=1, max1
         if (ymin_active(i))then
            call write_cutsnow('ymin        ',ymin(i),1,i)
         endif
      enddo
      do i=1, max1
         if (ymax_active(i))then
            call write_cutsnow('ymax        ',ymax(i),1,i)
         endif
      enddo
      do i=1, max1
         if (costhmin_active(i))then
            call write_cutsnow('costhmin    ',costhmin(i),1,i)
         endif
      enddo
      do i=1, max1
         if (costhmax_active(i))then
            call write_cutsnow('costhmax    ',costhmax(i),1,i)
         endif
      enddo
      do i=1, max1
         if (phimin_active(i))then
            call write_cutsnow('phimin      ',phimin(i),1,i)
         endif
      enddo
      do i=1, max1
         if (phimax_active(i))then
            call write_cutsnow('phimax      ',phimax(i),1,i)
         endif
      enddo
      do i=1, max1
         if (X1min_active(i))then
            call write_cutsnow('X1min       ',X1min(i),1,i)
         endif
      enddo
      do i=1, max1
         if (X1max_active(i))then
            call write_cutsnow('X1max       ',X1max(i),1,i)
         endif
      enddo
      do i=1, max1
         if (X2min_active(i))then
            call write_cutsnow('X2min       ',X2min(i),1,i)
         endif
      enddo
      do i=1, max1
         if (X2max_active(i))then
            call write_cutsnow('X2max       ',X2max(i),1,i)
         endif
      enddo
      do i=1, max1
         if (X3min_active(i))then
            call write_cutsnow('X3min       ',X3min(i),1,i)
         endif
      enddo
      do i=1, max1
         if (X3max_active(i))then
            call write_cutsnow('X3max       ',X3max(i),1,i)
         endif
      enddo
      do i=1, max2
         if (mijmin_active(i))then
            call write_cutsnow('mijmin      ',mijmin(i),2,i)
         endif
      enddo
      do i=1, max2
         if (mijmax_active(i))then
            call write_cutsnow('mijmax      ',mijmax(i),2,i)
         endif
      enddo
      do i=1, max2
         if (dRijmin_active(i))then
            call write_cutsnow('dRijmin     ',dRijmin(i),2,i)
         endif
      enddo
      do i=1, max2
         if (dRijmax_active(i))then
            call write_cutsnow('dRijmax     ',dRijmax(i),2,i)
         endif
      enddo
      do i=1, max2
         if (delta_phimin_active(i))then
            call write_cutsnow('delta_phimin',delta_phimin(i),2,i)
         endif
      enddo
      do i=1, max2
         if (delta_phimax_active(i))then
            call write_cutsnow('delta_phimax',delta_phimax(i),2,i)
         endif
      enddo
      do i=1, max2
         if (delta_etamin_active(i))then
            call write_cutsnow('delta_etamin',delta_etamin(i),2,i)
         endif
      enddo
      do i=1, max2
         if (delta_etamax_active(i))then
            call write_cutsnow('delta_etamax',delta_etamax(i),2,i)
         endif
      enddo
      do i=1, max2
         if (kTijmin_active(i))then
            call write_cutsnow('kTijmin     ',kTijmin(i),2,i)
         endif
      enddo
      do i=1, max2
         if (kTijmax_active(i))then
            call write_cutsnow('kTijmax     ',kTijmax(i),2,i)
         endif
      enddo
      do i=1, max2
         if (XY1min_active(i))then
            call write_cutsnow('XY1min      ',XY1min(i),2,i)
         endif
      enddo
      do i=1, max2
         if (XY1max_active(i))then
            call write_cutsnow('XY1max      ',XY1max(i),2,i)
         endif
      enddo
      do i=1, max2
         if (XY2min_active(i))then
            call write_cutsnow('XY2min      ',XY2min(i),2,i)
         endif
      enddo
      do i=1, max2
         if (XY2max_active(i))then
            call write_cutsnow('XY2max      ',XY2max(i),2,i)
         endif
      enddo
      do i=1, max2
         if (XY3min_active(i))then
            call write_cutsnow('XY3min      ',XY3min(i),2,i)
         endif
      enddo
      do i=1, max2
         if (XY3max_active(i))then
            call write_cutsnow('XY3max      ',XY3max(i),2,i)
         endif
      enddo
      do i=1, max3
         if (XYZ1min_active(i))then
            call write_cutsnow('XYZ1min     ',XYZ1min(i),3,i)
         endif
      enddo
      do i=1, max3
         if (XYZ1max_active(i))then
            call write_cutsnow('XYZ1max     ',XYZ1max(i),3,i)
         endif
      enddo
      do i=1, max3
         if (XYZ2min_active(i))then
            call write_cutsnow('XYZ2min     ',XYZ2min(i),3,i)
         endif
      enddo
      do i=1, max3
         if (XYZ2max_active(i))then
            call write_cutsnow('XYZ2max     ',XYZ2max(i),3,i)
         endif
      enddo
      do i=1, max3
         if (XYZ3min_active(i))then
            call write_cutsnow('XYZ3min     ',XYZ3min(i),3,i)
         endif
      enddo
      do i=1, max3
         if (XYZ3max_active(i))then
            call write_cutsnow('XYZ3max     ',XYZ3max(i),3,i)
         endif
      enddo
      do i=1, max4
         if (XYZA1min_active(i))then
            call write_cutsnow('XYZA1min    ',XYZA1min(i),4,i)
         endif
      enddo
      do i=1, max4
         if (XYZA1max_active(i))then
            call write_cutsnow('XYZA1max    ',XYZA1max(i),4,i)
         endif
      enddo
      do i=1, max4
         if (XYZA2min_active(i))then
            call write_cutsnow('XYZA2min    ',XYZA2min(i),4,i)
         endif
      enddo
      do i=1, max4
         if (XYZA2max_active(i))then
            call write_cutsnow('XYZA2max    ',XYZA2max(i),4,i)
         endif
      enddo
      do i=1, max4
         if (XYZA3min_active(i))then
            call write_cutsnow('XYZA3min    ',XYZA3min(i),4,i)
         endif
      enddo
      do i=1, max4
         if (XYZA3max_active(i))then
            call write_cutsnow('XYZA3max    ',XYZA3max(i),4,i)
         endif
      enddo


      return
      end



      subroutine write_cutsnow(cutname,cutvalue,number,ii)
      implicit none
      include "info.inc"
      character*12 cutname
      double precision cutvalue
      integer number,ii
      integer c1(max1),p1(max1),c2(2,max2),p2(2,max2),i
      integer c3(3,max3),p3(3,max3),c4(4,max4),p4(4,max4)
      
      if (number.eq.1) then
         do i=1, max1
            c1(i)=    s1(i)/1000
            p1(i)=mod(s1(i),1000)
         enddo
         write (52,*) cutname,'('//name(c1(ii),p1(ii))(1:lname(c1(ii),p1(ii)))//') =',cutvalue
      elseif (number.eq.2)then
         do i=1, max2
            c2(1,i)=    s2(1,i)/1000
            p2(1,i)=mod(s2(1,i),1000)
            c2(2,i)=    s2(2,i)/1000
            p2(2,i)=mod(s2(2,i),1000)
         enddo
         write (52,*) cutname,'('//name(c2(1,ii),p2(1,ii))(1:lname(c2(1,ii),p2(1,ii)))//','
     &                           //name(c2(2,ii),p2(2,ii))(1:lname(c2(2,ii),p2(2,ii)))//') =',cutvalue
      elseif (number.eq.3)then
         do i=1, max3
            c3(1,i)=    s3(1,i)/1000
            p3(1,i)=mod(s3(1,i),1000)
            c3(2,i)=    s3(2,i)/1000
            p3(2,i)=mod(s3(2,i),1000)
            c3(3,i)=    s3(3,i)/1000
            p3(3,i)=mod(s3(3,i),1000)
         enddo
         write (52,*) cutname,'('//name(c3(1,ii),p3(1,ii))(1:lname(c3(1,ii),p3(1,ii)))//','
     &                           //name(c3(2,ii),p3(2,ii))(1:lname(c3(2,ii),p3(2,ii)))//','
     &                           //name(c3(3,ii),p3(3,ii))(1:lname(c3(3,ii),p3(3,ii)))//') =',cutvalue
      elseif (number.eq.4)then
         do i=1, max3
            c4(1,i)=    s4(1,i)/1000
            p4(1,i)=mod(s4(1,i),1000)
            c4(2,i)=    s4(2,i)/1000
            p4(2,i)=mod(s4(2,i),1000)
            c4(3,i)=    s4(3,i)/1000
            p4(3,i)=mod(s4(3,i),1000)
            c4(4,i)=    s4(4,i)/1000
            p4(4,i)=mod(s4(4,i),1000)
         enddo
         write (52,*) cutname,'('//name(c4(1,ii),p4(1,ii))(1:lname(c4(1,ii),p4(1,ii)))//','
     &                           //name(c4(2,ii),p4(2,ii))(1:lname(c4(2,ii),p4(2,ii)))//','
     &                           //name(c4(3,ii),p4(3,ii))(1:lname(c4(3,ii),p4(3,ii)))//','
     &                           //name(c4(4,ii),p4(4,ii))(1:lname(c4(4,ii),p4(4,ii)))//') =',cutvalue
      endif

      return
      end

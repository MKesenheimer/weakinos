      Subroutine storecolor(color,iflow)
!**********************************************
! Stores color information of color(...) in larger
! color matrix color_tot(2,1000,0:2maxlines,0:maxterms)
!**********************************************

      implicit none

      include 'params.inc'

c    Arguments
      integer color(0:2*maxlines,0:maxfactors,0:maxterms),iflow

c    Local
      integer ifactor,nterms,iterm,i

c    Global
      integer color_tot(2,200,0:2*maxlines,0:maxfactors,0:maxterms),
     .        maxigraphs,max_ifactor(maxterms)
      integer graphnr
      common /dipol/ color_tot,maxigraphs,max_ifactor,graphnr

      if(iflow .gt.graphnr) then
      graphnr=graphnr+1
      nterms=color(0,0,0)
      color_tot(1,graphnr,0,0,0)=color(0,0,0)

      do iterm=1,nterms
      max_ifactor(iterm)=color(0,0,iterm)
      color_tot(1,graphnr,1,1,iterm)=color(1,1,iterm)
      color_tot(1,graphnr,2,1,iterm)=color(2,1,iterm)

      do ifactor=0,max_ifactor(iterm)
       color_tot(1,graphnr,0,ifactor,iterm)=color(0,ifactor,iterm)
       color_tot(1,graphnr,1,ifactor,iterm)=color(1,ifactor,iterm)
       color_tot(1,graphnr,2,ifactor,iterm)=color(2,ifactor,iterm)
        if(abs(color(0,ifactor,iterm)) .gt.2) then
          do i=3,abs(color(0,ifactor,iterm))
            color_tot(1,graphnr,i,ifactor,iterm)=color(i,ifactor,iterm)
          enddo
        endif
      enddo
      enddo
      endif

      end


      Subroutine getflowinfo(li,info)

      implicit none

      include 'params.inc'

c    Global
      integer flowinfo(max_particles,2)
      common /flow_info/ flowinfo

      character*(max_string) iwave(max_particles),owave(max_particles)
      character*(8) str2(3,max_particles)
      integer info_p(5,max_particles),iposx(3,max_particles)
      common/to_external/iwave,owave,iposx,info_p,str2
      integer        iline(-maxlines:maxlines),idir(-maxlines:maxlines)
      integer        this_coup(max_coup) ,goal_coup(max_coup)
      common/to_proc/iline,idir,this_coup,goal_coup
      integer         nincoming,nincfirst
      common/to_proc2/nincoming,nincfirst

c    Arguments
      integer li(4),info(0:4)

c    Local
      integer i,nr,found,nexternal,cl,cu

      nexternal=iline(0)

      nr=0
      if(info(0).eq.4 .or. info(0).eq.1 .or.info(0).eq.2) then
      if(info(0).eq.4) then
        cu=3
      else
        cu=3
      endif
      do i=1,cu
       found=0
       if(li(i) .gt. 0 .and. li(i) .le. nexternal) then
        if(nr .gt. 0 .and. (info(0).eq.4)) then
          if(li(i) .le. nincoming) then
            flowinfo(li(i),1)=0
            flowinfo(li(i),2)=1
            found=found+1
          endif
          if(li(i).gt. nincoming) then
            flowinfo(li(i),1)=1
            flowinfo(li(i),2)=-1
            found=found+1
          endif
        endif

        if(nr.eq.0 .and.(info(0).eq.4)) then
         if(li(i) .le. nincoming) then
            flowinfo(li(i),1)=0
            flowinfo(li(i),2)=-1
            found=found+1
          endif
          if(li(i).gt. nincoming) then
            flowinfo(li(i),1)=1
            flowinfo(li(i),2)=1
            found=found+1
          endif
        endif

         if(nr .eq. 0.and. found .gt.0) nr=nr+1
       endif
      enddo

       if(info(0).eq.4.and.li(1).lt.0.and.li(2).lt.0) then
        if(li(3).gt. nincoming) then
          flowinfo(li(3),1)=1
          flowinfo(li(3),2)=-1
        endif
        if(li(3).lt.nincoming.and.li(3).gt.0) then
          flowinfo(li(3),1)=0
          flowinfo(li(3),2)=1
        endif
       endif

      if(info(0).eq.2) then
       do i=1,3
        if(li(i).gt.0) then
         flowinfo(li(i),1)=0
         flowinfo(li(i),2)=0
        endif
       enddo
      endif

      if(info(0).eq.1) then
        if(li(1).gt.0) then
         flowinfo(li(1),1)=0
         flowinfo(li(1),2)=0
        endif
        if(li(2).gt.0.and. li(2).le.nincoming) then
           flowinfo(li(2),1)=0
           flowinfo(li(2),2)=1
        endif
        if(li(2).gt.0.and. li(2).gt.nincoming) then
           flowinfo(li(2),1)=1
           flowinfo(li(2),2)=-1
        endif
        if(li(3).gt.0.and. li(3).le.nincoming) then
           flowinfo(li(3),1)=0
           flowinfo(li(3),2)=-1
        endif
        if(li(3).gt.0.and. li(3).gt.nincoming) then
           flowinfo(li(3),1)=1
           flowinfo(li(3),2)=1
        endif
      endif
      endif
c      print*, flowinfo(1,1),flowinfo(1,2)
c      print*, flowinfo(2,1),flowinfo(2,2)
c      print*, flowinfo(3,1),flowinfo(3,2)
c      print*, flowinfo(4,1),flowinfo(4,2)
c      print*, flowinfo(5,1),flowinfo(5,2)
      end


      Subroutine getflowinfo1()
!**********************************************
! Provides information whether the external particles
! are particles/antiparticles and in-/outgoing
!**********************************************

      implicit none
      include 'params.inc'

      integer        iline(-maxlines:maxlines),idir(-maxlines:maxlines)
      integer        this_coup(max_coup) ,goal_coup(max_coup)
      common/to_proc/iline,idir,this_coup,goal_coup

c  Global
      integer flowinfo(max_particles,2)
      common /flow_info/ flowinfo
      integer         nincoming,nincfirst
      common/to_proc2/nincoming,nincfirst

c  Local
      integer i,nexternal

      nexternal=iline(0)

      do i=1,nexternal
       if(i.le. nincoming) then
c  1. initial quark
        if(mod(iline(i),2).ne.0) then
          flowinfo(i,1)=0
          flowinfo(i,2)=1
        endif
c  2. initial anti-quark
        if(mod(iline(i),2).eq.0) then
         flowinfo(i,1)=0
         flowinfo(i,2)=-1
        endif
       endif
       if(i.gt.nincoming) then
c  3. final quark
        if(mod(iline(i),2).eq.0) then
          flowinfo(i,1)=1
          flowinfo(i,2)=1
        endif
c  4. final anti-quark
        if(mod(iline(i),2).ne.0) then
          flowinfo(i,1)=1
          flowinfo(i,2)=-1
        endif
       endif
c  5. gluon or no qcd particle
       if(iline(i) .gt. 25) then
         flowinfo(i,1)=0
         flowinfo(i,2)=0
       endif

      enddo
c      print*, iline(1),iline(2),iline(3),iline(4),iline(5)
c      print*, flowinfo(1,1),flowinfo(1,2)
c      print*, flowinfo(2,1),flowinfo(2,2)
c      print*, flowinfo(3,1),flowinfo(3,2)
c      print*, flowinfo(4,1),flowinfo(4,2)
c      print*, flowinfo(5,1),flowinfo(5,2)

      end


      Subroutine getchargeinfo()
! **********************************************
! Provides the electromagnetic charge of the external
! particles. This is needed for QED dipoles.
! **********************************************
      implicit none
      include 'params.inc'


c  Global
      integer        iline(-maxlines:maxlines),idir(-maxlines:maxlines)
      integer        this_coup(max_coup) ,goal_coup(max_coup)
      common/to_proc/iline,idir,this_coup,goal_coup

      integer flowinfo(max_particles,2)
      common /flow_info/ flowinfo
      integer         nincoming,nincfirst
      common/to_proc2/nincoming,nincfirst

      real*8 chargeinfo(max_particles)
      integer chargint(maxdipoles,2)
      common /charge_info/ chargeinfo,chargint


c  Local
      integer i,nexternal

      nexternal=iline(0)

      do i=1, nexternal
c  1. initial state particles
       if(i.le.nincoming) then
c  2. particles or antiparticles 
        if(mod(iline(i),2).ne.0) then
c     now it's a particle, but quark,lepton,etc.?
         if(iline(i).le.25) then
c      uptype quark
          if(iline(i).eq.7.or.iline(i).eq.15.or.iline(i).eq.23) then
             chargeinfo(i)=-2d0/3d0
c      downtype quark
          elseif(iline(i).eq.3.or.iline(i).eq.11.or.iline(i).eq.19) then
             chargeinfo(i)=1d0/3d0
          endif
         elseif(iline(i).gt.25.and.iline(i).le.37) then
          chargeinfo(i)=1d0
         elseif(iline(i).gt.37) then
          chargeinfo(i)=0d0
         endif
        elseif(mod(iline(i),2).eq.0) then
c     now it's an anti-particle
         if(iline(i).le.25) then
c      uptype quark
          if(iline(i).eq.6.or.iline(i).eq.14.or.iline(i).eq.22) then
             chargeinfo(i)=2d0/3d0
c      downtype quark
          elseif(iline(i).eq.2.or.iline(i).eq.10.or.iline(i).eq.18) then
             chargeinfo(i)=-1d0/3d0
          endif
         elseif(iline(i).gt.25.and.iline(i).le.37) then
          chargeinfo(i)=-1d0
         elseif(iline(i).gt.37) then
          chargeinfo(i)=0d0
         endif
        endif


c  3. final state particles
       elseif(i.gt.nincoming) then
c  4. particle or antiparticle
        if(mod(iline(i),2).ne.0) then
c     now it's an antiparticle
         if(iline(i).le.25) then
c      uptype quark
          if(iline(i).eq.7.or.iline(i).eq.15.or.iline(i).eq.23) then
             chargeinfo(i)=-2d0/3d0
c      downtype quark
          elseif(iline(i).eq.3.or.iline(i).eq.11.or.iline(i).eq.19) then
             chargeinfo(i)=1d0/3d0
          endif
         elseif(iline(i).gt.25.and.iline(i).le.37) then
          chargeinfo(i)=1d0
         elseif(iline(i).gt.37) then
          chargeinfo(i)=0d0
         endif
        elseif(mod(iline(i),2).eq.0) then
c     now it's a particle
         if(iline(i).le.25) then
c      uptype quark
          if(iline(i).eq.6.or.iline(i).eq.14.or.iline(i).eq.22) then
             chargeinfo(i)=2d0/3d0
c      downtype quark
          elseif(iline(i).eq.2.or.iline(i).eq.10.or.iline(i).eq.18) then
             chargeinfo(i)=-1d0/3d0
          endif
         elseif(iline(i).gt.25.and.iline(i).le.37) then
          chargeinfo(i)=-1d0
         elseif(iline(i).gt.37) then
          chargeinfo(i)=0d0
         endif
        endif


       endif

      enddo

      end



      Subroutine Insertdipol(color,emitter,spectator)
!******************************************
! Insert the color information of emitter and
! spectator in the color matrix color_tot
!*****************************************

      implicit none

      include 'params.inc'

c    Arguments
      integer color(0:2*maxlines,0:maxfactors,0:maxterms)
      integer emitter, spectator

c    Global
      integer color_tot(2,200,0:2*maxlines,0:maxfactors,0:maxterms),
     .        maxigraphs,max_ifactor(maxterms)
      integer graphnr
      common /dipol/ color_tot,maxigraphs,max_ifactor,graphnr

      character*(max_string) iwave(max_particles),owave(max_particles)
      character*(8) str2(3,max_particles)
      integer info_p(5,max_particles),iposx(3,max_particles)
      common/to_external/iwave,owave,iposx,info_p,str2

      integer flowinfo(max_particles,2)
      common /flow_info/ flowinfo

c    Local
      integer nterms,iterm,j

      do j=1,graphnr
      nterms=color_tot(1,j,0,0,0)
      do iterm=1,nterms
       max_ifactor(iterm)=color(0,0,iterm)!-1
c    quark as emitter
      if(flowinfo(emitter,2) .ne.0) then
        color_tot(1,j,0,max_ifactor(iterm)+1,iterm)=-3
c    final state quark or initial state anti-quark as emitter
        if((flowinfo(emitter,1) .eq. 1 .and. flowinfo(emitter,2).eq.1) .or.
     &      (flowinfo(emitter,1) .eq. 0 .and. flowinfo(emitter,2).eq.-1))then
          color_tot(1,j,1,max_ifactor(iterm)+1,iterm)=emitter+max_particles
          color_tot(1,j,2,max_ifactor(iterm)+1,iterm)=emitter
          color_tot(1,j,3,max_ifactor(iterm)+1,iterm)=-99
        endif
c    final state anti-quark or initial state quark as emitter
        if((flowinfo(emitter,1) .eq. 1 .and. flowinfo(emitter,2).eq.-1).or.
     &       (flowinfo(emitter,1) .eq. 0 .and. flowinfo(emitter,2).eq.1))then
          color_tot(1,j,1,max_ifactor(iterm)+1,iterm)=emitter
          color_tot(1,j,2,max_ifactor(iterm)+1,iterm)=emitter+max_particles
          color_tot(1,j,3,max_ifactor(iterm)+1,iterm)=-99
        endif
      endif

c    quark as spectator
       if(flowinfo(spectator,2).ne.0) then
        color_tot(1,j,0,max_ifactor(iterm)+2,iterm)=-3
c    final state quark or initial state anti-quark as spectator
        if((flowinfo(spectator,1) .eq. 1 .and. flowinfo(spectator,2).eq.1) .or.
     &      (flowinfo(spectator,1) .eq. 0 .and. flowinfo(spectator,2).eq.-1))then
          color_tot(1,j,1,max_ifactor(iterm)+2,iterm)=spectator+max_particles
          color_tot(1,j,2,max_ifactor(iterm)+2,iterm)=spectator
          color_tot(1,j,3,max_ifactor(iterm)+2,iterm)=-99
        endif
c    final state anti-quark or initial state quark as spectator
        if((flowinfo(spectator,1) .eq. 1 .and. flowinfo(spectator,2).eq.-1).or.
     &       (flowinfo(spectator,1) .eq. 0 .and. flowinfo(spectator,2).eq.1))then
          color_tot(1,j,1,max_ifactor(iterm)+2,iterm)=spectator
          color_tot(1,j,2,max_ifactor(iterm)+2,iterm)=spectator+max_particles
          color_tot(1,j,3,max_ifactor(iterm)+2,iterm)=-99
        endif
      endif
      

c    gluon as spectator
      if(flowinfo(spectator,2).eq.0) then
         color_tot(1,j,0,max_ifactor(iterm)+2,iterm)=3
         color_tot(1,j,1,max_ifactor(iterm)+2,iterm)=spectator+max_particles
         color_tot(1,j,2,max_ifactor(iterm)+2,iterm)=-99
         color_tot(1,j,3,max_ifactor(iterm)+2,iterm)=spectator
      endif

c    gluon as emitter
      if(flowinfo(emitter,2).eq.0) then
         color_tot(1,j,0,max_ifactor(iterm)+1,iterm)=3
         color_tot(1,j,1,max_ifactor(iterm)+1,iterm)=emitter+max_particles
         color_tot(1,j,2,max_ifactor(iterm)+1,iterm)=-99
         color_tot(1,j,3,max_ifactor(iterm)+1,iterm)=emitter
      endif

      enddo
      enddo
c      print*,'flowinfo'
c      print*, flowinfo(1,1),flowinfo(1,2)
c      print*, flowinfo(2,1),flowinfo(2,2)
c      print*, flowinfo(3,1),flowinfo(3,2)
c      print*, flowinfo(4,1),flowinfo(4,2)
c      print*, flowinfo(5,1),flowinfo(5,2)

      end



      Subroutine doublecolor_tot(color,emitter,spectator)
!***********************************************
! Calculates the conjugate color matrix of the
! amplitudes and stores it in color_tot(2,....)
!***********************************************

      implicit none

      include 'params.inc'

c    Arguments
      integer emitter,spectator
      integer color(0:2*maxlines,0:maxfactors,0:maxterms)


c    Global
      integer color_tot(2,200,0:2*maxlines,0:maxfactors,0:maxterms),
     .        maxigraphs,max_ifactor(maxterms),graphnr
      common /dipol/ color_tot,maxigraphs,max_ifactor,graphnr

      character*(max_string) iwave(max_particles),owave(max_particles)
      character*(8) str2(3,max_particles)
      integer info_p(5,max_particles),iposx(3,max_particles)
      common/to_external/iwave,owave,iposx,info_p,str2

      integer flowinfo(max_particles,2)
      common /flow_info/ flowinfo

c    Local
      integer i,j,k,ifactor,factor,nterms,iterm,c,c1,ntermsmax,ratio,l,ncp
      logical ne,ns

      ntermsmax=1
      c=0
      if(flowinfo(spectator,2).eq.0) then
       c=c+1
      endif
      if(flowinfo(emitter,2).eq.0) then
       c=c+1
      endif
      do j=1,graphnr
       if(color_tot(1,j,0,0,0) .gt. ntermsmax) then
         ntermsmax=color_tot(1,j,0,0,0)
       endif
      enddo

      do j=1,graphnr
      nterms=color_tot(1,j,0,0,0)

      if(c.ne.0)then
       do iterm=1,nterms
        max_ifactor(nterms+iterm)=max_ifactor(iterm)
        color(0,0,nterms+iterm)=color(0,0,iterm) 
        color(0,1,nterms+iterm)=2
        color(1,1,nterms+iterm)=color(1,1,iterm)
        color(2,1,nterms+iterm)=color(2,1,iterm)
        color_tot(1,j,1,1,nterms+iterm)=color(1,1,nterms+iterm)
        color_tot(1,j,2,1,nterms+iterm)=color(2,1,nterms+iterm)
        color_tot(1,j,0,0,nterms+iterm)=color(0,0,nterms+iterm) 
       enddo
       do i=1,4*c**2
       do iterm=1,nterms
        color(0,0,i*nterms+iterm)=color(0,0,(i-1)*nterms+iterm) 
        color(0,1,i*nterms+iterm)=2
        max_ifactor(i*nterms+iterm)=max_ifactor((i-1)*nterms+iterm)
        color(1,1,i*nterms+iterm)=color(1,1,(i-1)*nterms+iterm)
        color(2,1,i*nterms+iterm)=color(2,1,(i-1)*nterms+iterm)
        color_tot(1,j,1,1,i*nterms+iterm)=color(1,1,(i-1)*nterms+iterm)
        color_tot(1,j,2,1,i*nterms+iterm)=color(2,1,(i-1)*nterms+iterm)
        color_tot(1,j,0,0,i*nterms+iterm)=color(0,0,(i-1)*nterms+iterm)
       enddo
       enddo
      endif

       do iterm=1,nterms
        do i=1,nterms-1
         max_ifactor(i*nterms+iterm)=max_ifactor(iterm)
        enddo
       enddo

c  if emitter or spectator is a gluon, enlarge colormatrix
      if(flowinfo(spectator,2).eq.0) then
       do iterm=1,nterms
c    copy terms
        do i=0,max_ifactor(iterm)+2
          do k=0,abs(color_tot(1,j,0,i,iterm))
          color_tot(1,j,k,i,nterms+iterm)=color_tot(1,j,k,i,iterm)
          enddo
        enddo
c    replace spectator term by modified one
        color_tot(1,j,0,max_ifactor(iterm)+2,nterms+iterm)=3
        color_tot(1,j,1,max_ifactor(iterm)+2,nterms+iterm)=spectator
        color_tot(1,j,2,max_ifactor(iterm)+2,nterms+iterm)=-99
        color_tot(1,j,3,max_ifactor(iterm)+2,nterms+iterm)=spectator+max_particles
       enddo

      endif

      if(flowinfo(emitter,2).eq.0) then
       if(flowinfo(spectator,2).ne.0)then
        do iterm=1,nterms
c    copy terms
        do i=0,max_ifactor(iterm)+2
          do k=0,abs(color_tot(1,j,0,i,iterm))
          color_tot(1,j,k,i,nterms+iterm)=color_tot(1,j,k,i,iterm)
          enddo
        enddo
c    replace emitter term by modified one
        color_tot(1,j,0,max_ifactor(iterm)+1,nterms+iterm)=3
        color_tot(1,j,1,max_ifactor(iterm)+1,nterms+iterm)=emitter
        color_tot(1,j,2,max_ifactor(iterm)+1,nterms+iterm)=-99
        color_tot(1,j,3,max_ifactor(iterm)+1,nterms+iterm)=emitter+max_particles
        enddo
       endif
       if(flowinfo(spectator,2).eq.0)then
        do iterm=1,2*nterms
c     copy terms
         do i=0,max_ifactor(iterm)+2
          do k=0,abs(color_tot(1,j,0,i,iterm))
           color_tot(1,j,k,i,2*nterms+iterm)=color_tot(1,j,k,i,iterm)
c          print*, j,i,iterm,k,color_tot(1,j,k,i,iterm)
          enddo
         enddo
        color_tot(1,j,0,max_ifactor(iterm)+1,2*nterms+iterm)=3
        color_tot(1,j,1,max_ifactor(iterm)+1,2*nterms+iterm)=emitter
        color_tot(1,j,2,max_ifactor(iterm)+1,2*nterms+iterm)=-99
        color_tot(1,j,3,max_ifactor(iterm)+1,2*nterms+iterm)=emitter+max_particles
        enddo
       endif
      endif

      if(c.ne.0) then
       c1=2*c
      else
        c1=1
      endif
       print*,ntermsmax,max_particles

c    copying numerator and denominator
      if(ntermsmax .gt. nterms) then
        do iterm=1,c1*nterms
         color(1,1,c1*nterms+iterm)=color(1,1,iterm)
         color(2,1,c1*nterms+iterm)=color(2,1,iterm)
         color_tot(1,j,1,1,c1*nterms+iterm)=color(1,1,c1*nterms+iterm)
         color_tot(1,j,2,1,c1*nterms+iterm)=color(2,1,c1*nterms+iterm)
         do i=1,2*max_ifactor(iterm)+2
          do k=0,abs(color_tot(1,j,0,i,iterm))
          color_tot(1,j,k,i,c1*nterms+iterm)=color_tot(1,j,k,i,iterm)
          enddo
         enddo
        enddo
      endif


      do iterm=1,c1*nterms
       ne=.true.
       ns=.true.
       max_ifactor(iterm)=color_tot(1,j,0,0,iterm)
       do ifactor=1,max_ifactor(iterm)
        color_tot(2,j,0,ifactor,iterm)=color_tot(1,j,0,ifactor,iterm)
        do k=1,7
         if(color_tot(1,j,k,ifactor,iterm).eq.emitter) ne=.false.
         if(color_tot(1,j,k,ifactor,iterm).eq. spectator) ns=.false.
        enddo


c    Kronecker Delta
        if(color_tot(1,j,0,ifactor,iterm) .eq. 2 .or. color_tot(1,j,0,ifactor,iterm) .eq. -2) then
c         do i=2,3
          if(color_tot(1,j,1,ifactor,iterm) .ne. 0) then
           if(color_tot(1,j,2,ifactor,iterm).eq.emitter .or. 
     &         color_tot(1,j,2,ifactor,iterm).eq.spectator.or.
     &         color_tot(1,j,2,ifactor,iterm).lt.0) then
               factor=1
           else
               factor=0
           endif
           color_tot(2,j,1,ifactor,iterm)
     &          = color_tot(1,j,2,ifactor,iterm)+
     &            sign(max_particles,color_tot(1,j,2,ifactor,iterm))*factor
          endif
          if(color_tot(1,j,2,ifactor,iterm) .ne. 0) then
           if(color_tot(1,j,1,ifactor,iterm).eq.emitter .or. 
     &         color_tot(1,j,1,ifactor,iterm).eq.spectator.or.
     &          color_tot(1,j,1,ifactor,iterm).lt.0) then
               factor=1
           else
               factor=0
           endif
           color_tot(2,j,2,ifactor,iterm)
     &          = color_tot(1,j,1,ifactor,iterm)+
     &            sign(max_particles,color_tot(1,j,1,ifactor,iterm))*factor
          endif
c         enddo
        endif

c   n colored particles
c   T[a,i,j]T[b,j,k]T[c,k,l]...
        if(color_tot(1,j,0,ifactor,iterm).le. -3 ) then
         ncp=abs(color_tot(1,j,0,ifactor,iterm))
           if(color_tot(1,j,2,ifactor,iterm).eq.emitter .or. 
     &         color_tot(1,j,2,ifactor,iterm).eq.spectator.or.
     &          color_tot(1,j,2,ifactor,iterm).lt.0) then
               factor=1
           else
               factor=0
           endif
          color_tot(2,j,1,ifactor,iterm) 
     &          = color_tot(1,j,2,ifactor,iterm)+
     &            sign(max_particles,color_tot(1,j,2,ifactor,iterm))*factor
           if(color_tot(1,j,1,ifactor,iterm).eq.emitter .or. 
     &         color_tot(1,j,1,ifactor,iterm).eq.spectator.or.
     &          color_tot(1,j,1,ifactor,iterm).lt.0) then
               factor=1
           else
               factor=0
           endif
          color_tot(2,j,2,ifactor,iterm) 
     &          = color_tot(1,j,1,ifactor,iterm)+
     &            sign(max_particles,color_tot(1,j,1,ifactor,iterm))*factor
          do l=3,ncp
           if(color_tot(1,j,ncp-l+3,ifactor,iterm).eq.emitter .or. 
     &         color_tot(1,j,ncp-l+3,ifactor,iterm).eq.spectator.or.
     &          color_tot(1,j,ncp-l+3,ifactor,iterm).lt.0) then
               factor=1
           else
               factor=0
           endif
          color_tot(2,j,l,ifactor,iterm) 
     &          = color_tot(1,j,ncp-l+3,ifactor,iterm)
     &            +sign(max_particles,color_tot(1,j,ncp-l+3,ifactor,iterm))*factor
          enddo

        endif

c   f[a,b,c]f[c,d,e]f[e,f,g]...
        if(color_tot(1,j,0,ifactor,iterm).ge. 3 ) then
         ncp=color_tot(1,j,0,ifactor,iterm)
         do l=1,ncp
           if(color_tot(1,j,ncp-l+1,ifactor,iterm).eq.emitter .or. 
     &         color_tot(1,j,ncp-l+1,ifactor,iterm).eq.spectator.or.
     &          color_tot(1,j,ncp-l+1,ifactor,iterm).lt.0) then
               factor=1
           else
               factor=0
           endif
          color_tot(2,j,l,ifactor,iterm) 
     &          = color_tot(1,j,ncp-l+1,ifactor,iterm)
     &            +sign(max_particles,color_tot(1,j,ncp-l+1,ifactor,iterm))*factor
         enddo
        endif

       enddo
       enddo


      enddo


      end


      Subroutine squarecolor_tot(color,color_tot,sqcolor,j,i,max_ifactor,emitter,spectator)
!********************************************
! Puts together color_tot(1,.....) of one amplitude
! and color_tot(2,....) of another amplitude into one
! color matrix sqcolor.
!********************************************

      implicit none

      include 'params.inc'

c    Arguments
      integer color_tot(2,200,0:2*maxlines,0:maxfactors,0:maxterms)
      integer sqcolor(0:2*maxlines,0:maxfactors,0:maxterms)
      integer i,j,max_ifactor(maxterms)
      integer color(0:2*maxlines,0:maxfactors,0:maxterms)
      integer emitter,spectator

c    Global
      integer flowinfo(max_particles,2)
      common /flow_info/ flowinfo

c    Local
      integer c,k,ifactor,iterm,itermi,itermj,ntermsi,ntermsj,l,nterms,m

      c=0
      if(flowinfo(spectator,2).eq.0) then
       c=c+1
      endif
      if(flowinfo(emitter,2).eq.0) then
       c=c+1
      endif

      ntermsi=color_tot(1,i,0,0,0)
      ntermsj=color_tot(1,j,0,0,0)
      sqcolor(0,0,0)=ntermsi*ntermsj
      if(c.ne.0) then
       sqcolor(0,0,0)=sqcolor(0,0,0)*2*c
      endif
        nterms=sqcolor(0,0,0)

        if(flowinfo(spectator,2).eq.0) then
         if(color_tot(1,i,0,0,0).ne.color_tot(1,j,0,0,0)) then
          do iterm=1,min(color_tot(1,i,0,0,0),color_tot(1,j,0,0,0))
           color_tot(1,j,0,max_ifactor(iterm)+2,iterm)=3
           color_tot(1,j,1,max_ifactor(iterm)+2,iterm)=spectator+max_particles
           color_tot(1,j,2,max_ifactor(iterm)+2,iterm)=-99
           color_tot(1,j,3,max_ifactor(iterm)+2,iterm)=spectator
          enddo
          do iterm=min(color_tot(1,i,0,0,0),color_tot(1,j,0,0,0))+1,nterms
           color_tot(1,j,0,max_ifactor(iterm)+2,iterm)=3
           color_tot(1,j,1,max_ifactor(iterm)+2,iterm)=spectator
           color_tot(1,j,2,max_ifactor(iterm)+2,iterm)=-99
           color_tot(1,j,3,max_ifactor(iterm)+2,iterm)=spectator+max_particles
          enddo
         endif
        endif

        if(flowinfo(emitter,2).eq.0) then
          if(color_tot(1,i,0,0,0).ne.color_tot(1,j,0,0,0)) then
           do iterm=1,min(color_tot(1,i,0,0,0),color_tot(1,j,0,0,0))
            color_tot(1,j,0,max_ifactor(iterm)+1,iterm)=3
            color_tot(1,j,1,max_ifactor(iterm)+1,iterm)=emitter+max_particles
            color_tot(1,j,2,max_ifactor(iterm)+1,iterm)=-99
            color_tot(1,j,3,max_ifactor(iterm)+1,iterm)=emitter
           enddo
           do iterm=min(color_tot(1,i,0,0,0),color_tot(1,j,0,0,0))+1,nterms
            color_tot(1,j,0,max_ifactor(iterm)+1,iterm)=3
            color_tot(1,j,1,max_ifactor(iterm)+1,iterm)=emitter
            color_tot(1,j,2,max_ifactor(iterm)+1,iterm)=-99
            color_tot(1,j,3,max_ifactor(iterm)+1,iterm)=emitter+max_particles
           enddo
          endif
        endif

      l=1
      if(c.ne.0) then
       ntermsi=2*c*ntermsi
      endif
      do itermi=1,ntermsi
       do itermj=1,ntermsj
        do ifactor=0,max_ifactor(itermi)+2
         do k=0,abs(color_tot(1,i,0,ifactor,itermi))
          sqcolor(k,ifactor,l) =color_tot(1,i,k,ifactor,itermi)
         enddo
        enddo
        do ifactor=max_ifactor(itermi)+3,2*max_ifactor(itermi)+1
         do k=0,abs(color_tot(2,j,0,ifactor-max_ifactor(itermi)-1,itermj))
         sqcolor(k,ifactor,l)=color_tot(2,j,k,ifactor-max_ifactor(itermi)-1,itermj)
         enddo
        enddo
       l=l+1
       enddo
      enddo

      l=1
      do itermi=1,ntermsi
       do itermj=1,ntermsj
        sqcolor(1,1,l)=color_tot(1,i,1,1,itermi)*color_tot(1,j,1,1,itermj)
        sqcolor(2,1,l)=color_tot(1,i,2,1,itermi)*color_tot(1,j,2,1,itermj)
        l=l+1
       enddo
      enddo

      end



      Subroutine dipolcolormatrix(sqcolor,color,i,j,emitter,spectator)
!**************************************************
! Copies sqcolor into color an evaluates the color-
! factor and returns the output in the usual variable
! color(...).
!***************************************************

      implicit none
      include 'params.inc'

c    Arguments
      integer color(0:2*maxlines,0:maxfactors,0:maxterms)
      integer i,j,emitter,spectator,sqcolor(0:2*maxlines,0:maxfactors,0:maxterms)

c    Global
      integer color_tot(2,200,0:2*maxlines,0:maxfactors,0:maxterms),
     .        maxigraphs,max_ifactor(maxterms),graphnr
      common /dipol/ color_tot,maxigraphs,max_ifactor,graphnr

      integer flowinfo(max_particles,2)
      common /flow_info/ flowinfo

      character*(max_string) iwave(max_particles),owave(max_particles)
      character*(8) str2(3,max_particles)
      integer info_p(5,max_particles),iposx(3,max_particles)
      common/to_external/iwave,owave,iposx,info_p,str2
      integer        iline(-maxlines:maxlines),idir(-maxlines:maxlines)
      integer        this_coup(max_coup) ,goal_coup(max_coup)
      common/to_proc/iline,idir,this_coup,goal_coup

      integer         nincoming,nincfirst
      common/to_proc2/nincoming,nincfirst
c    Local
      integer k,l,m,n,listintern(4),m1,m2,iterm,nterms,c
      logical check

      c=0
      if(flowinfo(spectator,2).eq.0) then
       c=c+1
      endif
      if(flowinfo(emitter,2).eq.0) then
       c=c+1
      endif

       nterms=sqcolor(0,0,0)
      color(0,0,0)=nterms
      do iterm=1,nterms
      color(1,1,iterm)=sqcolor(1,1,iterm)
      color(2,1,iterm)=sqcolor(2,1,iterm)

      do k=0,2*max_ifactor(iterm)+2
       do l=0,abs(sqcolor(0,k,iterm))
       color(l,k,iterm)=sqcolor(l,k,iterm)
       enddo
      enddo
      color(0,0,iterm)=2*max_ifactor(iterm)+1
      enddo

c     modifying numerator and denominator
      if(c.eq.1) then
       do iterm=1, nterms/2
        color(1,1,nterms/2+iterm)=2*color(1,1,iterm)
        color(2,1,nterms/2+iterm)=color(2,1,iterm)
        color(1,1,iterm)=-2*color(1,1,iterm)
       enddo
      endif

      if(c.eq.2) then
       do iterm=1, nterms/2/c
        color(1,1,nterms/2/c+iterm)=2*color(1,1,iterm)
        color(1,1,iterm)=-2*color(1,1,iterm)
       enddo
       do iterm=1,nterms/2
        color(1,1,nterms/2+iterm)=color(1,1,iterm)
        color(2,1,nterms/2+iterm)=color(2,1,iterm)
       enddo
       do iterm=1,nterms/2
        color(1,1,iterm)=-2*color(1,1,iterm)
        color(1,1,nterms/2+iterm)=2*color(1,1,nterms/2+iterm)
       enddo
      endif


      do n=1,4
       listintern(n)=0
      enddo

      n=1
      check=.true.
      do iterm=1,nterms
      do m=2, 2*max_ifactor(iterm)+4
       do m1=1,3
        check=.true.
        if(color(m1,m,iterm) .lt.0 .and. color(m1,m,iterm) .ne. -99) then
         do m2=1,n
          if(color(m1,m,iterm).eq.listintern(m2)) check=.false.
         enddo
         if(check.eqv..true.) then
          listintern(n)=color(m1,m,iterm)
          n=n+1
         endif
        endif
       enddo
      enddo
      enddo


      n=1
      do while(listintern(n) .ne. 0)
        call sumcolors(color,listintern(n))
        call clearzeros(color)
        n=n+1
      enddo
c      call printdipolcolor(color)
      call sumcolors(color,-99)
      call orderF(color)
      call addterms(color)
c      call ordercolors(color)
      call sumcolors(color,emitter+max_particles)
      call clearzeros(color)
      call sumcolors(color,spectator+max_particles)
      call clearzeros(color)
      call simpcolors(color,1)
      call addterms(color)
      if((flowinfo(emitter,1) .eq. 1 .and. flowinfo(emitter,2).eq.-1).or.
     &       (flowinfo(emitter,1) .eq. 0 .and. flowinfo(emitter,2).eq.1))then
           color(1,1,1)=-color(1,1,1)
      endif
      if((flowinfo(spectator,1) .eq. 1 .and. flowinfo(spectator,2).eq.-1).or.
     &       (flowinfo(spectator,1) .eq. 0 .and. flowinfo(spectator,2).eq.1))then
           color(1,1,1)=-color(1,1,1)
      endif
      if(flowinfo(spectator,2).eq.0) color(1,1,1)=-color(1,1,1)
      if(flowinfo(emitter,2).eq.0) color(1,1,1)=-color(1,1,1)

      print*, 'numerator and denominator: ',color(1,1,1),color(2,1,1)

      if(c.ne.0) then
       do k=1,graphnr
        color_tot(1,k,0,1,1)=color_tot(2,k,0,1,1)
       enddo
      endif

      end



      Subroutine orderdipol()

      implicit none

      include 'params.inc'

c    Global
      integer dipollist(maxdipoles,2,2),nrdipols
      common /list/ dipollist,nrdipols

      integer        iline(-maxlines:maxlines),idir(-maxlines:maxlines)
      integer        this_coup(max_coup) ,goal_coup(max_coup)
      common/to_proc/iline,idir,this_coup,goal_coup


c    Local
      integer i,j,k
      k=1
c      print*,'orderdipol',iline(0)
      do i=1,iline(0)
       do j=1,iline(0)
         if(iline(j) .lt. 24) then
          dipollist(k,1,1)=j
          dipollist(k,2,1)=3
         endif
         if(iline(j) .eq.50) then
          dipollist(k,1,1)=j
          dipollist(k,2,1)=8
         endif
         if(iline(i) .lt. 24) then
          dipollist(k,1,2)=i
          dipollist(k,2,2)=3
         endif
         if(iline(i) .eq.50) then
          dipollist(k,1,2)=i
          dipollist(k,2,2)=8
         endif
         if((iline(j).lt.24 .or. iline(j).eq.50) .and.
     &      (iline(i) .lt.24 .or. iline(i) .eq. 50).and.
     &        i.ne. j) then
           nrdipols=k
           k=k+1
         endif
       enddo
      enddo

      end



      Subroutine write_dipols()
!***************************************
! writes all dipoles in the file dipolsum.f
! and includes the theta function of the
! alpha parameter
!***************************************

      implicit none

      include 'params.inc'

c     Global
      integer dipollist(maxdipoles,2,2),nrdipols
      common /list/ dipollist,nrdipols

      integer info_ijk(maxdipoles,2,5),nr
      common /ijk/ info_ijk,nr

      character*70       cversion
      common/to_version/ cversion

      character*60 proc
      integer iproc
      common/to_write/iproc,proc

      integer         nincoming,nincfirst
      common/to_proc2/nincoming,nincfirst

      character*(8) mass(max_particles)
      common /masses/ mass

      integer iline_copy(-maxlines:maxlines,maxdipoles)
      common /to_writedipole/iline_copy
      character*(4*max_particles) particle(4)
      integer                               charge_c(max_particles)
      integer iparticle(0:max_particles,0:4),inverse(max_particles)
      common/to_model/iparticle,  particle,  inverse, charge_c
      integer info_p(5,max_particles),iposx(3,max_particles)
      character*(max_string) iwave(max_particles),owave(max_particles)
      character*(8) str(3,max_particles)
      common/to_external/iwave,owave,iposx,info_p,str


      integer iline_orig(-maxlines:maxlines)
      common /ciline_orig/iline_orig

c     Local
      integer ln,i,ascii_emitter,ascii_spectator,ascii_nr
      integer nrdipol
      parameter (ln=81)

      nrdipol=nr
      if(nrdipol.eq.0) nrdipol=1

      open(unit=ln,file='dipolsum.f',status='unknown')

      write(ln,10) 'SUBROUTINE DIPOLSUM(P,X,PSWGT,SUBTRACT)'
      write(ln,5) ' '
      write(ln,5) cversion
      write(ln,5) 'RETURNS SUM OVER ALL DIPOLS'
      write(ln,5) 'FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)'
      write(ln,5) ' '
      write(ln,5) 'FOR PROCESS '//proc(1:iproc)
      write(ln,5) ' '

      write(ln,10) 'IMPLICIT NONE'
      write(ln,10) 'include "nexternal.inc" '
      write(ln,10) 'include "dipole.inc" '
      write(ln,10) 'include "coupl.inc"'

      write(ln,5) ' '
      write(ln,5) 'ARGUMENTS '
      write(ln,5) ' '
      write(ln,10) 'REAL*8 P(0:3,NEXTERNAL),X(2),PSWGT,SUBTRACT'
      write(ln,5) ' '
      write(ln,5) 'GLOBAL '
      write(ln,5) ' '
      write(ln,10) 'REAL*8 dot'
      write(ln,10) 'LOGICAL PASSCUTSDIP'
      write(ln,10) 'EXTERNAL dot, PASSCUTSDIP'
      write(ln,5) ' '
      write(ln,5) 'LOCAL'
      write(ln,5) ' '
      write(ln,10) 'REAL*8 sij,sik,sjk,q(0:3),qsq,yp,P1(0:3,NEXTERNAL),pdfwgt'
      write(ln,10) 'REAL*8 y_ijk,omx_ija,u_i,vtilde_i'
      write(ln,10) 'REAL*8 mass_i,mass_j,mass_k,mass_ij'
      write(ln,10) 'REAL*8 msq_i,msq_j,msq_k,msq_ij,mu_i,mu_j,mu_k'
      write(ln,10) 'REAL*8 one,two,three,four'
      write(ln,10) 'parameter (one=1d0,two=2d0,three=3d0,four=4d0)'
      write(ln,60) 'REAL*8 ANS(',nrdipol,')'
      write(ln,10) 'REAL*8 ZERO'
      write(ln,10) 'INTEGER NU'
      write(ln,10) 'PARAMETER (ZERO=0D0)'
      write(ln,5) ' '
      write(ln,5) 'BEGIN CODE'
      write(ln,5) ' '
      write(ln,10) 'INTDIPOLE=.FALSE.'
      write(ln,10) 'SUBTRACT=0d0'
      do i=1,nr
c     First check if unresolved is massless. If not do not print it here.
         if (mass(info_ijk(i,1,2)).eq.'ZERO')then
c       final-final theta funtion
         if(info_ijk(i,1,1).gt.nincoming .and. info_ijk(i,1,3) .gt.nincoming) then
            write(ln,10) ' '
            write(ln,70)'sij=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,2),'))'
            write(ln,70)'sjk=dot(p(0,',info_ijk(i,1,2)
     .           ,'),p(0,',info_ijk(i,1,3),'))'
            write(ln,70)'sik=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,3),'))'
            write(ln,10)'y_ijk=sij/(sij+sjk+sik)'
           if(mass(info_ijk(i,1,3)).ne.'ZERO') then
            write(ln,10)'DO NU=0,3'
            write(ln,73)' q(nu)=p(nu,',info_ijk(i,1,1),')+p(nu,',info_ijk(i,1,2),
     &                    ')+p(nu,',info_ijk(i,1,3),')'
            write(ln,10)'ENDDO'
            write(ln,10)'qsq=dot(q,q)'
            write(ln,71)'mass_i=',mass(info_ijk(i,1,1))
            write(ln,71)'mass_j=',mass(info_ijk(i,1,2))
            write(ln,71)'mass_k=',mass(info_ijk(i,1,3))
            write(ln,10)'mu_i=mass_i/SQRT(qsq)'
            write(ln,10)'mu_j=mass_j/SQRT(qsq)'
            write(ln,10)'mu_k=mass_k/SQRT(qsq)'
            write(ln,10)'yp=1d0-2d0*mu_k*(1d0-mu_k)/(1d0-mu_i**2-mu_j**2-mu_k**2)'
            write(ln,14)'CALL TRANSFORM(P,P1,',info_ijk(i,1,1),',',info_ijk(i,1,2),',',
     &                   info_ijk(i,1,3),',',info_ijk(i,1,4),',',info_ijk(i,1,5),',',
     &                   mass(info_ijk(i,1,1)),',',mass(info_ijk(i,1,2)),',',mass(info_ijk(i,1,3)),')'
            write(ln,10)'IF(PASSCUTSDIP(P1)) THEN'
            write(ln,10)'IF (alpha_ff*yp .gt.y_ijk) THEN'
           else
            write(ln,14)'CALL TRANSFORM(P,P1,',info_ijk(i,1,1),',',info_ijk(i,1,2),',',
     &                   info_ijk(i,1,3),',',info_ijk(i,1,4),',',info_ijk(i,1,5),',',
     &                   mass(info_ijk(i,1,1)),',',mass(info_ijk(i,1,2)),',',mass(info_ijk(i,1,3)),')'
            write(ln,10)'IF(PASSCUTSDIP(P1)) THEN'
            write(ln,10)'IF (alpha_ff .gt.y_ijk) THEN'
           endif
            write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &           info_p(5,inverse(iline_orig(1))),',',
     &           info_p(5,inverse(iline_orig(2))),',pdfwgt)'
            if (i.le.9) then
               write(ln,40)
     .          'CALL SDIPOL'//char(48)//char(48)//char(48+i)//'(P,ANS(',i,'))'
            elseif(i.le.99) then
               write(ln,40)'CALL SDIPOL'//char(48)//char(i/10+48)
     &              //char(mod(i,10)+48)//'(P,ANS(',i,'))'
            elseif(i.le.999) then
               write(ln,41)'CALL SDIPOL'//char(i/100+48)
     &              //char((i-(i/100)*100)/10+48)
     &              //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
            else
               write (*,*) 'ERROR: write_dipoles: too many dipoles'
            endif
            write(ln,51)'SUBTRACT = SUBTRACT - ANS(',i,')*PDFWGT'
            write(ln,51)'CALL WRITEHIST(P1,ANS(',i,')*PSWGT*PDFWGT)'
            write(ln,10)'ENDIF'
            write(ln,10)'ENDIF'
         endif
c       final-initial theta function
         if(info_ijk(i,1,1).gt.nincoming .and. info_ijk(i,1,3) .le.nincoming) then
            write(ln,10) ' '
            write(ln,70)'sij=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,2),'))'
            write(ln,70)'sjk=dot(p(0,',info_ijk(i,1,2)
     .           ,'),p(0,',info_ijk(i,1,3),'))'
            write(ln,70)'sik=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,3),'))'
            if (mass(info_ijk(i,1,1)).ne.'ZERO'.or.
     .          mass(info_ijk(i,1,2)).ne.'ZERO') then
               write(ln,71)'mass_i=',mass(info_ijk(i,1,1))
               write(ln,71)'mass_j=',mass(info_ijk(i,1,2))
c               write(ln,71)'mass_k=',mass(info_ijk(i,1,3))
               write(ln,10)'msq_i=mass_i**2'
               write(ln,10)'msq_j=mass_j**2'
c            write(ln,10)'msq_k=mass_k**2'
               write(ln,10)'if(mass_i .eq.mass_j) then'
               write(ln,10)' msq_ij=0d0'
               write(ln,10)'else'
               write(ln,10)' msq_ij=msq_i'
               write(ln,10)'endif'
               write(ln,10)'omx_ija='//
     &              '(sij-0.5d0*(msq_ij-msq_i-msq_j))/(sjk+sik)'
            else
               write(ln,10)'omx_ija=sij/(sjk+sik)'
            endif
            write(ln,14)'CALL TRANSFORM(P,P1,',info_ijk(i,1,1),',',info_ijk(i,1,2),',',
     &                   info_ijk(i,1,3),',',info_ijk(i,1,4),',',info_ijk(i,1,5),',',
     &                   mass(info_ijk(i,1,1)),',',mass(info_ijk(i,1,2)),',',mass(info_ijk(i,1,3)),')'
            write(ln,10)'IF(PASSCUTSDIP(P1)) THEN'
            write(ln,10)'IF (alpha_fi .gt.omx_ija) THEN'
            write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &           info_p(5,inverse(iline_orig(1))),',',
     &           info_p(5,inverse(iline_orig(2))),',pdfwgt)'
            if (i.le.9) then
               write(ln,40)'CALL SDIPOL'//char(48)//char(48)
     .              //char(48+i)//'(P,ANS(',i,'))'
            elseif(i.le.99) then
               write(ln,40)'CALL SDIPOL'//char(48)//char(i/10+48)
     &              //char(mod(i,10)+48)//'(P,ANS(',i,'))'
            elseif(i.le.999) then
               write(ln,41)'CALL SDIPOL'//char(i/100+48)
     &              //char((i-(i/100)*100)/10+48)
     &              //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
            else
               write (*,*) 'ERROR: write_dipoles: too many dipoles'
            endif
            write(ln,51)'SUBTRACT = SUBTRACT - ANS(',i,')*PDFWGT'
            write(ln,51)'CALL WRITEHIST(P1,ANS(',i,')*PSWGT*PDFWGT)'
            write(ln,10)'ENDIF'
            write(ln,10)'ENDIF'
         endif
c       initial-final theta function
         if(info_ijk(i,1,1).le.nincoming .and. info_ijk(i,1,3) .gt.nincoming) then
            write(ln,10) ' '
            write(ln,70)'sij=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,2),'))'
            write(ln,70)'sik=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,3),'))'
            write(ln,10)'u_i=sij/(sij+sik)'
            write(ln,14)'CALL TRANSFORM(P,P1,',info_ijk(i,1,1),',',info_ijk(i,1,2),',',
     &                   info_ijk(i,1,3),',',info_ijk(i,1,4),',',info_ijk(i,1,5),',',
     &                   mass(info_ijk(i,1,1)),',',mass(info_ijk(i,1,2)),',',mass(info_ijk(i,1,3)),')'
            write(ln,10)'IF(PASSCUTSDIP(P1)) THEN'
            write(ln,10)'IF (alpha_if .gt.u_i) THEN'
            write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &           info_p(5,inverse(iline_orig(1))),',',
     &           info_p(5,inverse(iline_orig(2))),',pdfwgt)'
            if (i.le.9) then
               write(ln,40)'CALL SDIPOL'//char(48)//char(48)
     .              //char(48+i)//'(P,ANS(',i,'))'
            elseif(i.le.99) then
               write(ln,40)'CALL SDIPOL'//char(48)//char(i/10+48)
     &              //char(mod(i,10)+48)//'(P,ANS(',i,'))'
            elseif(i.le.999) then
               write(ln,41)'CALL SDIPOL'//char(i/100+48)
     &              //char((i-(i/100)*100)/10+48)
     &              //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
            else
               write (*,*) 'ERROR: write_dipoles: too many dipoles'
            endif
            write(ln,51)'SUBTRACT = SUBTRACT - ANS(',i,')*PDFWGT'
            write(ln,51)'CALL WRITEHIST(P1,ANS(',i,')*PSWGT*PDFWGT)'
            write(ln,10)'ENDIF'
            write(ln,10)'ENDIF'
         endif
c     initial-initial theta function
         if(info_ijk(i,1,1).le.nincoming .and. info_ijk(i,1,3) .le.nincoming) then
            write(ln,10) ' '
            write(ln,70)'sij=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,2),'))'
            write(ln,70)'sik=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,3),'))'
            write(ln,10)'vtilde_i=sij/sik'
            write(ln,14)'CALL TRANSFORM(P,P1,',info_ijk(i,1,1),',',info_ijk(i,1,2),',',
     &                   info_ijk(i,1,3),',',info_ijk(i,1,4),',',info_ijk(i,1,5),',',
     &                   mass(info_ijk(i,1,1)),',',mass(info_ijk(i,1,2)),',',mass(info_ijk(i,1,3)),')'
            write(ln,10)'IF(PASSCUTSDIP(P1)) THEN'
            write(ln,10)'IF (alpha_ii .gt.vtilde_i) THEN'
            write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &           info_p(5,inverse(iline_orig(1))),',',
     &           info_p(5,inverse(iline_orig(2))),',pdfwgt)'
            if (i.le.9) then
               write(ln,40)'CALL SDIPOL'//char(48)//char(48)
     .              //char(48+i)//'(P,ANS(',i,'))'
            elseif(i.le.99) then
               write(ln,40)'CALL SDIPOL'//char(48)//char(i/10+48)
     &              //char(mod(i,10)+48)//'(P,ANS(',i,'))'
            elseif(i.le.999) then
               write(ln,41)'CALL SDIPOL'//char(i/100+48)
     &              //char((i-(i/100)*100)/10+48)
     &              //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
            else
               write (*,*) 'ERROR: write_dipoles: too many dipoles'
            endif
            write(ln,51)'SUBTRACT = SUBTRACT - ANS(',i,')*PDFWGT'
            write(ln,51)'CALL WRITEHIST(P1,ANS(',i,')*PSWGT*PDFWGT)'
            write(ln,10)'ENDIF'
            write(ln,10)'ENDIF'
         endif
      endif
      enddo

      write(ln,10) 'RETURN'
      write(ln,10) 'END'
      write(ln,5) ' '
      write(ln,5) ' '
      write(ln,5) ' '
      write(ln,5) ' '
      write(ln,5) ' '

      write(ln,10) 'SUBROUTINE DIPOLSUMFINITE(P,X,PSWGT,SUBTRACT)'
      write(ln,5) ' '
      write(ln,5) cversion
      write(ln,5) 'RETURNS SUM OVER "NON-DIVERGENT" DIPOLS,'
      write(ln,5) 'I.E. DIPOLES FOR WHICH THE UNRESOLVED IS'
      write(ln,5) 'MASSIVE (IN particles.dat MASS IS NOT'
      write(ln,5) 'EQUAL TO "ZERO")'
      write(ln,5) 'FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)'
      write(ln,5) ' '
      write(ln,5) 'FOR PROCESS '//proc(1:iproc)
      write(ln,5) ' '

      write(ln,10) 'IMPLICIT NONE'
      write(ln,10) 'include "nexternal.inc" '
      write(ln,10) 'include "dipole.inc" '
      write(ln,10) 'include "coupl.inc"'

      write(ln,5) ' '
      write(ln,5) 'ARGUMENTS '
      write(ln,5) ' '
      write(ln,10) 'REAL*8 P(0:3,NEXTERNAL),X(2),PSWGT,SUBTRACT'
      write(ln,5) ' '
      write(ln,5) 'GLOBAL '
      write(ln,5) ' '
      write(ln,10) 'REAL*8 dot'
      write(ln,10) 'LOGICAL PASSCUTSDIP'
      write(ln,10) 'EXTERNAL dot, PASSCUTSDIP'
      write(ln,5) ' '
      write(ln,5) 'LOCAL'
      write(ln,5) ' '
      write(ln,10) 'REAL*8 sij,sik,sjk,q(0:3),qsq,yp,P1(0:3,NEXTERNAL),pdfwgt'
      write(ln,10) 'REAL*8 y_ijk,omx_ija,u_i,vtilde_i'
      write(ln,10) 'REAL*8 mass_i,mass_j,mass_k,mass_ij'
      write(ln,10) 'REAL*8 msq_i,msq_j,msq_k,msq_ij,mu_i,mu_j,mu_k'
      write(ln,10) 'REAL*8 one,two,three,four'
      write(ln,10) 'parameter (one=1d0,two=2d0,three=3d0,four=4d0)'
      write(ln,60) 'REAL*8 ANS(',nrdipol,')'
      write(ln,10) 'REAL*8 ZERO'
      write(ln,10) 'INTEGER NU'
      write(ln,10) 'PARAMETER (ZERO=0D0)'
      write(ln,5) ' '
      write(ln,5) 'BEGIN CODE'
      write(ln,5) ' '
      write(ln,10) 'INTDIPOLE=.FALSE.'
      write(ln,10) 'SUBTRACT=0d0'
      do i=1,nr
         if (mass(info_ijk(i,1,2)).ne.'ZERO')then
c       final-final theta funtion
         if(info_ijk(i,1,1).gt.nincoming .and. info_ijk(i,1,3) .gt.nincoming) then
            write(ln,70)'sij=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,2),'))'
            write(ln,70)'sjk=dot(p(0,',info_ijk(i,1,2)
     .           ,'),p(0,',info_ijk(i,1,3),'))'
            write(ln,70)'sik=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,3),'))'
            write(ln,10)'y_ijk=sij/(sij+sjk+sik)'
           if(mass(info_ijk(i,1,3)).ne.'ZERO') then
            write(ln,10)'DO NU=0,3'
            write(ln,73)' q(nu)=p(nu,',info_ijk(i,1,1),')+p(nu,',info_ijk(i,1,2),
     &                    ')+p(nu,',info_ijk(i,1,3),')'
            write(ln,10)'ENDDO'
            write(ln,10)'qsq=dot(q,q)'
            write(ln,71)'mass_i=',mass(info_ijk(i,1,1))
            write(ln,71)'mass_j=',mass(info_ijk(i,1,2))
            write(ln,71)'mass_k=',mass(info_ijk(i,1,3))
            write(ln,10)'mu_i=mass_i/SQRT(qsq)'
            write(ln,10)'mu_j=mass_j/SQRT(qsq)'
            write(ln,10)'mu_k=mass_k/SQRT(qsq)'
            write(ln,10)'yp=1d0-2d0*mu_k*(1d0-mu_k)/(1d0-mu_i**2-mu_j**2-mu_k**2)'
           write(ln,14)'CALL TRANSFORM(P,P1,',info_ijk(i,1,1),',',info_ijk(i,1,2),',',
     &                   info_ijk(i,1,3),',',info_ijk(i,1,4),',',info_ijk(i,1,5),',',
     &                   mass(info_ijk(i,1,1)),',',mass(info_ijk(i,1,2)),',',mass(info_ijk(i,1,3)),')'
            write(ln,10)'IF(PASSCUTSDIP(P1)) THEN'
            write(ln,10)'IF (alpha_ff*yp .gt.y_ijk) THEN'
           else
            write(ln,14)'CALL TRANSFORM(P,P1,',info_ijk(i,1,1),',',info_ijk(i,1,2),',',
     &                   info_ijk(i,1,3),',',info_ijk(i,1,4),',',info_ijk(i,1,5),',',
     &                   mass(info_ijk(i,1,1)),',',mass(info_ijk(i,1,2)),',',mass(info_ijk(i,1,3)),')'
            write(ln,10)'IF(PASSCUTSDIP(P1)) THEN'
            write(ln,10)'IF (alpha_ff .gt.y_ijk) THEN'
           endif
            write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &           info_p(5,inverse(iline_orig(1))),',',
     &           info_p(5,inverse(iline_orig(2))),',pdfwgt)'
            if (i.le.9) then
               write(ln,40)
     .          'CALL SDIPOL'//char(48)//char(48)//char(48+i)//'(P,ANS(',i,'))'
            elseif(i.le.99) then
               write(ln,40)'CALL SDIPOL'//char(48)//char(i/10+48)
     &              //char(mod(i,10)+48)//'(P,ANS(',i,'))'
            elseif(i.le.999) then
               write(ln,41)'CALL SDIPOL'//char(i/100+48)
     &              //char((i-(i/100)*100)/10+48)
     &              //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
            else
               write (*,*) 'ERROR: write_dipoles: too many dipoles'
            endif
            write(ln,51)'SUBTRACT = SUBTRACT - ANS(',i,')*PDFWGT'
            write(ln,51)'CALL WRITEHIST(P1,ANS(',i,')*PSWGT*PDFWGT)'
            write(ln,10)'ENDIF'
            write(ln,10)'ENDIF'
         endif
c       final-initial theta function
         if(info_ijk(i,1,1).gt.nincoming .and. info_ijk(i,1,3) .le.nincoming) then
            write(ln,70)'sij=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,2),'))'
            write(ln,70)'sjk=dot(p(0,',info_ijk(i,1,2)
     .           ,'),p(0,',info_ijk(i,1,3),'))'
            write(ln,70)'sik=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,3),'))'
            if (mass(info_ijk(i,1,1)).ne.'ZERO'.or.
     .          mass(info_ijk(i,1,2)).ne.'ZERO') then
               write(ln,71)'mass_i=',mass(info_ijk(i,1,1))
               write(ln,71)'mass_j=',mass(info_ijk(i,1,2))
c               write(ln,71)'mass_k=',mass(info_ijk(i,1,3))
               write(ln,10)'msq_i=mass_i**2'
               write(ln,10)'msq_j=mass_j**2'
c            write(ln,10)'msq_k=mass_k**2'
               write(ln,10)'if(mass_i .eq.mass_j) then'
               write(ln,10)' msq_ij=0d0'
               write(ln,10)'else'
               write(ln,10)' msq_ij=msq_i'
               write(ln,10)'endif'
               write(ln,10)'omx_ija=(sij-0.5d0*(msq_ij-msq_i-msq_j))/(sjk+sik)'
            else
               write(ln,10)'omx_ija=sij/(sjk+sik)'
            endif
            write(ln,14)'CALL TRANSFORM(P,P1,',info_ijk(i,1,1),',',info_ijk(i,1,2),',',
     &                   info_ijk(i,1,3),',',info_ijk(i,1,4),',',info_ijk(i,1,5),',',
     &                   mass(info_ijk(i,1,1)),',',mass(info_ijk(i,1,2)),',',mass(info_ijk(i,1,3)),')'
            write(ln,10)'IF(PASSCUTSDIP(P1)) THEN'
            write(ln,10)'IF (alpha_fi .gt.omx_ija) THEN'
            write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &           info_p(5,inverse(iline_orig(1))),',',
     &           info_p(5,inverse(iline_orig(2))),',pdfwgt)'
            if (i.le.9) then
               write(ln,40)'CALL SDIPOL'//char(48)//char(48)
     .              //char(48+i)//'(P,ANS(',i,'))'
            elseif(i.le.99) then
               write(ln,40)'CALL SDIPOL'//char(48)//char(i/10+48)
     &              //char(mod(i,10)+48)//'(P,ANS(',i,'))'
            elseif(i.le.999) then
               write(ln,41)'CALL SDIPOL'//char(i/100+48)
     &              //char((i-(i/100)*100)/10+48)
     &              //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
            else
               write (*,*) 'ERROR: write_dipoles: too many dipoles'
            endif
            write(ln,51)'SUBTRACT = SUBTRACT - ANS(',i,')*PDFWGT'
            write(ln,51)'CALL WRITEHIST(P1,ANS(',i,')*PSWGT*PDFWGT)'
            write(ln,10)'ENDIF'
            write(ln,10)'ENDIF'
         endif
c       initial-final theta function
         if(info_ijk(i,1,1).le.nincoming .and. info_ijk(i,1,3) .gt.nincoming) then
            write(ln,70)'sij=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,2),'))'
            write(ln,70)'sik=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,3),'))'
            write(ln,10)'u_i=sij/(sij+sik)'
            write(ln,14)'CALL TRANSFORM(P,P1,',info_ijk(i,1,1),',',info_ijk(i,1,2),',',
     &                   info_ijk(i,1,3),',',info_ijk(i,1,4),',',info_ijk(i,1,5),',',
     &                   mass(info_ijk(i,1,1)),',',mass(info_ijk(i,1,2)),',',mass(info_ijk(i,1,3)),')'
            write(ln,10)'IF(PASSCUTSDIP(P1)) THEN'
            write(ln,10)'IF (alpha_if .gt.u_i) THEN'
            write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &           info_p(5,inverse(iline_orig(1))),',',
     &           info_p(5,inverse(iline_orig(2))),',pdfwgt)'
            if (i.le.9) then
               write(ln,40)'CALL SDIPOL'//char(48)//char(48)
     .              //char(48+i)//'(P,ANS(',i,'))'
            elseif(i.le.99) then
               write(ln,40)'CALL SDIPOL'//char(48)//char(i/10+48)
     &              //char(mod(i,10)+48)//'(P,ANS(',i,'))'
            elseif(i.le.999) then
               write(ln,41)'CALL SDIPOL'//char(i/100+48)
     &              //char((i-(i/100)*100)/10+48)
     &              //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
            else
               write (*,*) 'ERROR: write_dipoles: too many dipoles'
            endif
            write(ln,51)'SUBTRACT = SUBTRACT - ANS(',i,')*PDFWGT'
            write(ln,51)'CALL WRITEHIST(P1,ANS(',i,')*PSWGT*PDFWGT)'
            write(ln,10)'ENDIF'
            write(ln,10)'ENDIF'
         endif
c     initial-initial theta function
         if(info_ijk(i,1,1).le.nincoming .and. info_ijk(i,1,3) .le.nincoming) then
            write(ln,70)'sij=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,2),'))'
            write(ln,70)'sik=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,3),'))'
            write(ln,10)'vtilde_i=sij/sik'
            write(ln,14)'CALL TRANSFORM(P,P1,',info_ijk(i,1,1),',',info_ijk(i,1,2),',',
     &                   info_ijk(i,1,3),',',info_ijk(i,1,4),',',info_ijk(i,1,5),',',
     &                   mass(info_ijk(i,1,1)),',',mass(info_ijk(i,1,2)),',',mass(info_ijk(i,1,3)),')'
            write(ln,10)'IF(PASSCUTSDIP(P1)) THEN'
            write(ln,10)'IF (alpha_ii .gt.vtilde_i) THEN'
            if (i.le.9) then
               write(ln,40)'CALL SDIPOL'//char(48)//char(48)
     .              //char(48+i)//'(P,ANS(',i,'))'
            elseif(i.le.99) then
               write(ln,40)'CALL SDIPOL'//char(48)//char(i/10+48)
     &              //char(mod(i,10)+48)//'(P,ANS(',i,'))'
            elseif(i.le.999) then
               write(ln,41)'CALL SDIPOL'//char(i/100+48)
     &              //char((i-(i/100)*100)/10+48)
     &              //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
            else
               write (*,*) 'ERROR: write_dipoles: too many dipoles'
            endif
            write(ln,51)'SUBTRACT = SUBTRACT - ANS(',i,')*PDFWGT'
            write(ln,51)'CALL WRITEHIST(P1,ANS(',i,')*PSWGT*PDFWGT)'
            write(ln,10)'ENDIF'
            write(ln,10)'ENDIF'
         endif
      endif
      enddo

      write(ln,10) 'RETURN'
      write(ln,10) 'END'
      write(ln,5) ' '
      write(ln,5) ' '
      write(ln,5) ' '
      write(ln,5) ' '
      write(ln,5) ' '

      write(ln,10)'SUBROUTINE DIPOLEPDF(x1,x2,leg1,leg2,WGT)'
      write(ln,10)'IMPLICIT NONE'
      write(ln,5) 'This is a template for a routine to calculate the pdfs of the'
      write(ln,5) 'single dipoles. The output variable is assumed to'
      write(ln,5) 'be wgt. The parton fractions are given by x1 and x2.'
      write(ln,5)
     &     'leg1 and leg2 are the PDG codes for two incoming partons.'
      write(ln,5)'The factorization scale is included via a common block in dipole.inc'
      write(ln,10) ' '
      write(ln,10)'include "nexternal.inc"'
      write(ln,10)'include "dipole.inc"'
      write(ln,5)'include "run.inc"'
      write(ln,10)'double precision X1,X2,WGT'
      write(ln,10)'INTEGER leg1,leg2'
      write(ln,10)' '
      write(ln,10)'IF (x1.gt.1d0 .or. x2.gt.1d0) THEN'
      write(ln,10)'   WGT=0d0'
      write(ln,10)'   RETURN'
      write(ln,10)'ENDIF'
      write(ln,5)'This is an empty routine, so return 1'
      write(ln,10) 'WGT=1d0'
      write(ln,10) ' '
      write(ln,10) 'END'



      write(ln,10) ' '
      write(ln,10) ' '
      write(ln,10) ' '


      write(ln,10)'SUBROUTINE WRITEHIST(P,WGT)'
      write(ln,10)'IMPLICIT NONE'
      write(ln,5) 'This is a template that can be used for writing events or histograms.'
      write(ln,5) 'The input variables are the phase space point and the according'
      write(ln,5) 'weight.'
      write(ln,10)' '
      write(ln,10)'include "nexternal.inc"'
      write(ln,10)'DOUBLE PRECISION P(0:3,NEXTERNAL),WGT'
      write(ln,10)' '
      write(ln,10)' '
      write(ln,10)'END'

      write(ln,10) ' '
      write(ln,10) ' '
      write(ln,10) ' '

      write(ln,10)'LOGICAL FUNCTION PASSCUTSDIP(P)'
      write(ln,10)'IMPLICIT NONE'
      write(ln,5) 'This is a template for a routine that puts cuts on each'
      write(ln,5) 'individual dipole, e.g. a jet algorithm. The input variable'
      write(ln,5) 'is the phase space point P. If event passes the cut, this'
      write(ln,5) 'function has to be set TRUE.'
      write(ln,10) ' '
      write(ln,10)'include "nexternal.inc"'
      write(ln,10)'DOUBLE PRECISION P(0:3,NEXTERNAL)'
      write(ln,10)' '
      write(ln,10)' '
      write(ln,5) 'Set to TRUE by default.'
      write(ln,10)'PASSCUTSDIP=.TRUE.'
      write(ln,10)' '
      write(ln,10)'END'      
      




 5    format('C ',A)
 10   format(6X,A)
 11   format(6X,A,I2,A)
 12   format(5X,'&',A,I2,A)
 13   format(6X,A,I2,A)
 14   format(6X,A,5(I2,A),A8,A,A8,A,A8,A)
 15   format('C',6X,A)
 17   format(6X,A,I3,A,I3,A)
 20   format(', P',I1)
 30   format(',P',I1,'(0:3)')
 40   format(6X,A22,I4,A2)
 41   format(6X,A23,I4,A2)
 50   format(7X,A26,I4,A1)
 51   format(7X,A,I4,A)
 60   format(6X,A11,I4,A1)
 70   format(6X,A12,I2,A6,I2,A2)
 71   format(6X,A,A5)
 72   format(6X,A15,2(I1,A7),I1,A)
 73   format(6X,A,I2,A,I2,A,I2,A)
      close(ln)
      end



      Subroutine write_intdipoles()


      implicit none

      include 'params.inc'

c     Global
      integer dipollist(maxdipoles,2,2),nrdipols
      common /list/ dipollist,nrdipols

      integer info_ijk(maxdipoles,2,5),nr
      common /ijk/ info_ijk,nr

      character*70       cversion
      common/to_version/ cversion

      character*60 proc
      integer iproc
      common/to_write/iproc,proc

      integer         nincoming,nincfirst
      common/to_proc2/nincoming,nincfirst

      character*(8) mass(max_particles)
      common /masses/ mass

      integer info_p(5,max_particles),iposx(3,max_particles)
      character*(max_string) iwave(max_particles),owave(max_particles)
      character*(8) str(3,max_particles)
      common/to_external/iwave,owave,iposx,info_p,str

      integer iline_orig(-maxlines:maxlines)
      common /ciline_orig/iline_orig

      integer iline_copy(-maxlines:maxlines,maxdipoles)
      common /to_writedipole/iline_copy

      integer        iline(-maxlines:maxlines),idir(-maxlines:maxlines)
      integer        this_coup(max_coup) ,goal_coup(max_coup)
      common/to_proc/iline,idir,this_coup,goal_coup

      character*(max_string) dipprocess(maxdipoles)
      common /process/ dipprocess

      character*(4) part_type(maxdipoles,2)
      common /part/ part_type

      character*(4*max_particles) particle(4)
      integer                               charge_c(max_particles)
      integer iparticle(0:max_particles,0:4),inverse(max_particles)
      common/to_model/iparticle,  particle,  inverse, charge_c

c     Local
      integer ln,i,ascii_emitter,ascii_spectator,ascii_nr
      integer nrdipol
      parameter (ln=82)

      nrdipol=nr
      if(nrdipol.eq.0) nrdipol=1

      open(unit=ln,file='intdipoles.f',status='unknown')

      write(ln,10) 'SUBROUTINE INTDIPOLES(P,X,Z,PSWGT,EPS,EPSSQ,FIN)'
      write(ln,5) ' '
      write(ln,5) cversion
      write(ln,5) 
     .     'RETURNS COEFFICIENTS OF 1/eps^2,1/eps AND FINITE TERMS'
      write(ln,5) 'FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL) '
      write(ln,5) ' '
      write(ln,5) 'FOR PROCESS '//proc(1:iproc)
      write(ln,5) 'WHERE Z IS THE MOMENTUM FRACTION '//
     &     'OF THE INITAL EMITTER/SPECTATOR.'
      write(ln,5) ' '

      write(ln,10) 'IMPLICIT NONE'
      write(ln,10) 'include "nexternal.inc" '
      write(ln,10) 'include "dipole.inc" '
      write(ln,10) 'include "coupl.inc"'

      write(ln,5) ' '
      write(ln,5) 'ARGUMENTS '
      write(ln,5) ' '
      write(ln,10) 'REAL*8 P(0:3,NEXTERNAL),X(2),Z,PSWGT'
      write(ln,10) 'REAL*8 EPS,EPSSQ,FIN'
      write(ln,5) ' '
      write(ln,5) 'GLOBAL '
      write(ln,5) ' '
      write(ln,10) 'REAL*8 dot,lambda_tr'
      write(ln,10) 'LOGICAL PASSCUTSDIP'
      write(ln,10) 'EXTERNAL dot,lambda_tr,PASSCUTSDIP'
      write(ln,5) ' '
      write(ln,5) 'LOCAL'
      write(ln,5) ' '
      write(ln,10) 'REAL*8 sij,sik,sjk,zp,e,esq,pdfwgt(3),FINITE(3)'
      write(ln,10) 'REAL*8 y_ijk,omx_ija,u_i,vtilde_i,zplus'
      write(ln,10) 'REAL*8 mass_i,mass_j,mass_k,mass_ij'
      write(ln,10) 'REAL*8 msq_i,msq_j,msq_k,msq_ij,sikzone'
      write(ln,10) 'REAL*8 musq_i,musq_j,musq_k,musq_ij'
      write(ln,10) 'REAL*8 vij,vtildeij,rho_i,rho_k,rho,qsq'
      write(ln,10) 'REAL*8 one,two,three,four,cf,ca,qtot(0:3),pi'
      write(ln,10) 'REAL*8 gsq,gammaE,ff,fi(9)'
      write(ln,10) 'parameter (one=1d0,two=2d0,three=3d0,four=4d0)'
      write(ln,10) 'parameter (cf=4d0/3d0,ca=3d0)'
      write(ln,60) 'REAL*8 ANS(',nrdipol,')'
      write(ln,10) 'REAL*8 ZERO'
      write(ln,10) 'INTEGER NU'
      write(ln,10) 'PARAMETER (ZERO=0D0)'
      write(ln,10) 'PARAMETER (gammaE=0.57721566d0)'
      write(ln,10) 'PARAMETER (pi=3.1415926d0)'
      write(ln,5) ' '
      write(ln,5) 'BEGIN CODE'
      write(ln,5) ' '
      write(ln,10) 'INTDIPOLE=.TRUE.'
      write(ln,10) 'gsq=GG(1)**2'
      write(ln,10) 'eps=0d0'
      write(ln,10) 'epssq=0d0'
      write(ln,10) 'fin=0d0'
      write(ln,10) 'e=0d0'
      write(ln,10) 'esq=0d0'
      write(ln,10) 'DO NU=1,3'
      write(ln,10) '   finite(NU)=0d0'
      write(ln,10) 'ENDDO'
      write(ln,60) 'DO NU=1,',nrdipol
      write(ln,10) '   ANS(NU)=0d0'
      write(ln,10) 'ENDDO'
      write(ln,10) 'IF (PASSCUTSDIP(P)) THEN'
      do i=1,nr
c     First check if unresolved is massless. If not do not print it here.
         if (mass(info_ijk(i,1,2)).eq.'ZERO')then
c**************************************************************
c       final-final theta funtion
c**************************************************************
         if(info_ijk(i,1,1).gt.nincoming
     &        .and. info_ijk(i,1,3).gt.nincoming) then
            write(ln,10)' '
            write(ln,71)'mass_i=',mass(info_ijk(i,1,1))
            write(ln,71)'mass_k=',mass(info_ijk(i,1,3))
            write(ln,10) 'DO NU=0,3'
            write(ln,17) '   qtot(NU)=p(NU,',info_ijk(i,1,4),
     &           ')+p(NU,',info_ijk(i,1,5),')'
            write(ln,10) 'ENDDO'
            write(ln,10)'qsq=dot(qtot,qtot)'

            write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &           info_p(5,inverse(iline_copy(1,i))),',',
     &           info_p(5,inverse(iline_copy(2,i))),',pdfwgt(2))'
            if (i.le.9) then
               write(ln,40)
     &          'CALL SDIPOL'//char(48)//char(48)//char(48+i)
     &              //'(P,ANS(',i,'))'
            elseif(i.le.99) then
               write(ln,40)'CALL SDIPOL'//char(48)//char(i/10+48)
     &              //char(mod(i,10)+48)//'(P,ANS(',i,'))'
            elseif(i.le.999) then
               write(ln,42)'CALL SDIPOL'//char(i/100+48)
     &              //char((i-(i/100)*100)/10+48)
     &              //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
            else
               write (*,*) 'ERROR: write_dipoles: too many dipoles'
            endif
            write(ln,14)
     &           'CALL EPSFF(mass_i,mass_k,qsq,'
     &           ,info_ijk(i,2,4),',',info_ijk(i,2,1),',esq,e)'
            write(ln,11) 'epssq=epssq-esq*ANS(',i,')*pdfwgt(2)'
            write(ln,11) 'eps=eps-e*ANS(',i,')*pdfwgt(2)'
            write(ln,14)
     &           'CALL FINITEFF(mass_i,mass_k,qsq,'
     &           ,info_ijk(i,2,4),',',info_ijk(i,2,1),',fi)'
            write(ln,11)
     &           'finite(2)=finite(2)-fi(2)*ANS(',i,')*pdfwgt(2)'

c**************************************************************
c       final-initial theta function
c**************************************************************
         elseif(info_ijk(i,1,1).gt.nincoming
     &           .and.info_ijk(i,1,3).le.nincoming) then
            write(ln,10)' '
            write(ln,71)'mass_i=',mass(info_ijk(i,1,1))
            write(ln,74)'sik=two*dot(p(0,',info_ijk(i,1,4),'),p(0,'
     &           ,info_ijk(i,1,5),'))/z'
            if(info_ijk(i,2,4).eq.0 .and. (
     &           mass(info_ijk(i,1,1)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,2)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,3)).ne.'ZERO')) then
               write (ln,11) 'zplus=1d0-4d0*mass_i**2/sik'
               write(ln,74)'sikzone=two*dot(p(0,',info_ijk(i,1,4),
     &              '),p(0,',info_ijk(i,1,5),
     &              '))/zplus'
            else
               write(ln,74)'sikzone=two*dot(p(0,',info_ijk(i,1,4),
     &              '),p(0,',info_ijk(i,1,5),'))'
            endif
            if(info_ijk(i,1,3).eq.1) then
               write(ln,17)'CALL DIPOLEPDF(x(1)/z,x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(1))'
            elseif(info_ijk(i,1,3).eq.2) then
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2)/z,',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(1))'
            else
               write (ln,*)'Error #145 in write_intdipole'
               write (*,*)'Error #145 in write_intdipole'
            endif
            if(info_ijk(i,2,4).eq.0 .and. (
     &           mass(info_ijk(i,1,1)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,2)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,3)).ne.'ZERO')) then
               if(info_ijk(i,1,3).eq.1) then
                  write(ln,17)'CALL DIPOLEPDF(x(1)/zplus,x(2),',
     &                 info_p(5,inverse(iline_copy(1,i))),',',
     &                 info_p(5,inverse(iline_copy(2,i))),',pdfwgt(3))'
               elseif(info_ijk(i,1,3).eq.2) then
                  write(ln,17)'CALL DIPOLEPDF(x(1),x(2)/zplus,',
     &                 info_p(5,inverse(iline_copy(1,i))),',',
     &                 info_p(5,inverse(iline_copy(2,i))),',pdfwgt(3))'
               endif
            else
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(2))'
            endif
            if (i.le.9) then
               write(ln,40)'CALL SDIPOL'//char(48)//char(48)
     &              //char(48+i)//'(P,ANS(',i,'))'
            elseif(i.le.99) then
               write(ln,40)'CALL SDIPOL'//char(48)//char(i/10+48)
     &              //char(mod(i,10)+48)//'(P,ANS(',i,'))'
            elseif(i.le.999) then
               write(ln,42)'CALL SDIPOL'//char(i/100+48)
     &              //char((i-(i/100)*100)/10+48)
     &              //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
            else
               write (*,*) 'ERROR: write_dipoles: too many dipoles'
            endif
            write(ln,14)
     &           'CALL EPSFI(',info_ijk(i,2,4),',',
     &           info_ijk(i,2,1),',sikzone,mass_i,esq,e)'
            write(ln,11) 'epssq=epssq-esq*ANS(',i,')*pdfwgt(2)'
            write(ln,11) 'eps=eps-e*ANS(',i,')*pdfwgt(2)'
            write(ln,14) 'CALL FINITEFI(mass_i,sik,sikzone,z,'
     &           ,info_ijk(i,2,4),',',info_ijk(i,2,1),',fi)'

            if(info_ijk(i,2,4).eq.0 .and. (
     &           mass(info_ijk(i,1,1)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,2)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,3)).ne.'ZERO')) then
               write(ln,11)
     &            'finite(2)=finite(2)-(-fi(8))*ANS(',i,')*pdfwgt(2)'
               write(ln,11) 'finite(3)=finite(3)-(fi(7)+fi(9))'//
     &              '*ANS(',i,')*pdfwgt(3)'
            else
               write(ln,11) 'finite(1)=finite(1)-(fi(1)+fi(3))'//
     &              '*ANS(',i,')*pdfwgt(1)'
               write(ln,11) 'finite(2)=finite(2)-(fi(2)-fi(4))'//
     &              '*ANS(',i,')*pdfwgt(2)'
            endif

c**************************************************************
c       initial-final theta function
c**************************************************************
         elseif(info_ijk(i,1,1).le.nincoming .and.
     &           info_ijk(i,1,3).gt.nincoming) then
            write(ln,10)' '
            write(ln,71)'mass_k=',mass(info_ijk(i,1,3))
            write(ln,74)'sik=two*dot(p(0,',info_ijk(i,1,4)
     &           ,'),p(0,',info_ijk(i,1,5),'))/z'
            write(ln,74)'sikzone=two*dot(p(0,',info_ijk(i,1,4)
     &           ,'),p(0,',info_ijk(i,1,5),'))'
            if (info_ijk(i,1,1).eq.1) then
               write(ln,17)'CALL DIPOLEPDF(x(1)/z,x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(1))'
            elseif(info_ijk(i,1,1).eq.2) then
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2)/z,',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(1))'
            else
               write (ln,*)'Error #146 in write_intdipole'
               write (*,*)'Error #146 in write_intdipole'
            endif
            if(((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.1).or.
     &           (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.0))) then
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(2))'
            endif
            if (i.le.9) then
               write(ln,41)'CALL SDIPOL'//char(48)//char(48)
     &              //char(48+i)//'(P,ANS(',i,'))'
            elseif(i.le.99) then
               write(ln,41)'CALL SDIPOL'//char(48)//char(i/10+48)
     &              //char(mod(i,10)+48)//'(P,ANS(',i,'))'
            elseif(i.le.999) then
               write(ln,41)'CALL SDIPOL'//char(i/100+48)
     &              //char((i-(i/100)*100)/10+48)
     &              //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
            else
               write (*,*) 'ERROR: write_dipoles: too many dipoles'
            endif

            if(((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.1).or.
     &           (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.0))) then
               write(ln,14) 'CALL EPSIF(',info_ijk(i,2,1),','
     &              ,info_ijk(i,2,4),',sikzone,mass_k,esq,e)'
               write(ln,11) 'epssq=epssq-esq*ANS(',i,')*pdfwgt(2)'
               write(ln,11) 'eps=eps-e*ANS(',i,')*pdfwgt(2)'
            endif
            write(ln,14)
     &           'CALL FINITEIF(mass_k,sik,sikzone,z,'
     &           ,info_ijk(i,2,1),',',info_ijk(i,2,4),',fi)'
c
c quark-quark splitting (i:quark, ij~: quark) and
c gluon-gluon splitting (i:gluon, ij~: gluon)
c
            if(((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.1).or.
     &           (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.0))) then
               if (mass(info_ijk(i,1,3)).eq.'ZERO') then
                  write(ln,11) 'finite(1)=finite(1)-(fi(1)+fi(3))'//
     &                 '*ANS(',i,')*pdfwgt(1)'
                  write(ln,11) 'finite(2)=finite(2)-(fi(2)-fi(4))'//
     &                 '*ANS(',i,')*pdfwgt(2)'
               else
                  write(ln,11) 'finite(1)=finite(1)-'//
     &                 '(fi(1)+fi(3)+fi(5))*ANS(',i,')*pdfwgt(1)'
                  write(ln,11) 'finite(2)=finite(2)-'//
     &                 '(fi(2)-fi(4)-fi(6))*ANS(',i,')*pdfwgt(2)'
               endif
c
c quark-gluon splitting (i:quark, ij~: gluon) and
c gluon-quark splitting (i:gluon, ij~: gluon)
c
            elseif((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.0) .or.
     &             (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.1)) then
             write(ln,11) 'finite(1)=finite(1)-fi(1)'//
     &              '*ANS(',i,')*pdfwgt(1)'
            else
               write (*,*) 'Error #142 in dipolcolor'
               write (ln,*) 'Error #142 in dipolcolor'
            endif

c**************************************************************
c     initial-initial theta function
c**************************************************************
         elseif(info_ijk(i,1,1).le.nincoming .and.
     &           info_ijk(i,1,3) .le.nincoming) then
            write(ln,10)' '
            write(ln,74)'sik=two*dot(p(0,',info_ijk(i,1,4),
     &           '),p(0,',info_ijk(i,1,5),'))/z'
            write(ln,74)'sikzone=two*dot(p(0,',info_ijk(i,1,4),
     &           '),p(0,',info_ijk(i,1,5),'))'
            if (info_ijk(i,1,1).eq.1) then
               write(ln,17)'CALL DIPOLEPDF(x(1)/z,x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(1))'
            elseif(info_ijk(i,1,1).eq.2) then
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2)/z,',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(1))'
            else
               write (ln,*)'Error #147 in write_intdipole'
               write (*,*)'Error #147 in write_intdipole'
            endif
            if(((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.1).or.
     &         (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.0))) then
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(2))'
            endif
            if (i.le.9) then
               write(ln,41)'CALL SDIPOL'//char(48)//char(48)
     &           //char(48+i)//'(P,ANS(',i,'))'
            elseif(i.le.99) then
               write(ln,41)'CALL SDIPOL'//char(48)//char(i/10+48)
     &           //char(mod(i,10)+48)//'(P,ANS(',i,'))'
            elseif(i.le.999) then
               write(ln,41)'CALL SDIPOL'//char(i/100+48)
     &           //char((i-(i/100)*100)/10+48)
     &           //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
            else
               write (*,*) 'ERROR: write_dipoles: too many dipoles'
            endif
            if(((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.1).or.
     &         (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.0))) then
               write(ln,14) 'CALL EPSII(',info_ijk(i,2,1),','
     &              ,info_ijk(i,2,4),',sikzone,esq,e)'
               write(ln,11) 'epssq=epssq-esq*ANS(',i,')*pdfwgt(2)'
               write(ln,11) 'eps=eps-e*ANS(',i,')*pdfwgt(2)'
            endif
            write(ln,14) 'CALL FINITEII(sik,sikzone,z,'
     &           ,info_ijk(i,2,1),',',info_ijk(i,2,4),',fi)'
c
c quark-quark splitting (i:quark, ij~: quark) and
c gluon-gluon splitting (i:gluon, ij~: gluon)
c
            if(((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.1).or.
     &         (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.0))) then
               if (mass(info_ijk(i,1,3)).eq.'ZERO') then
                  write(ln,11) 'finite(1)=finite(1)-(fi(1)+fi(3))'//
     &                 '*ANS(',i,')*pdfwgt(1)'
                  write(ln,11) 'finite(2)=finite(2)-(fi(2)-fi(4))'//
     &                 '*ANS(',i,')*pdfwgt(2)'
               else
                  write(ln,11) 'finite(1)=finite(1)-'//
     &                 '(fi(1)+fi(3)+fi(5))*ANS(',i,')*pdfwgt(1)'
                  write(ln,11) 'finite(2)=finite(2)-'//
     &                 '(fi(2)-fi(4)-fi(6))*ANS(',i,')*pdfwgt(2)'
               endif
c
c quark-gluon splitting (i:quark, ij~: gluon) and
c gluon-quark splitting (i:gluon, ij~: gluon)
c
            elseif((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.0) .or.
     &             (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.1)) then
               write(ln,11)
     &              'finite(1)=finite(1)-fi(1)*ANS(',i,')*pdfwgt(1)'
            else
               write (*,*) 'Error #144 in dipolcolor'
               write (ln,*) 'Error #144 in dipolcolor'
            endif
         else
            write (*,*) 'Error #143 in dipolcolor'
            write (ln,*) 'Error #143 in dipolcolor'
         endif
         endif ! massless unresolved
      enddo
      write(ln,10) ' '
      write(ln,10) 'ENDIF'
      write(ln,10) 'FIN=(FINITE(1)/z+FINITE(2))*PSWGT'
      write(ln,10) 'CALL WRITEHIST(P,FIN)'
      write(ln,10) 'INTDIPOLE=.FALSE.'
      write(ln,10) 'RETURN'
      write(ln,10) 'END'
      write(ln,10) ' '
      write(ln,10) ' '
      write(ln,10) ' '
      write(ln,10) ' '
      write(ln,10) ' '

      write(ln,10)
     &     'SUBROUTINE INTDIPOLESFINITE(P,X,Z,PSWGT,EPS,EPSSQ,FIN)'
      write(ln,5) ' '
      write(ln,5) cversion
      write(ln,5) 
     .     'RETURNS COEFFICIENTS OF 1/eps^2,1/eps AND FINITE TERMS'
      write(ln,5) 'FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL) '
      write(ln,5) ' '
      write(ln,5) 'FOR PROCESS '//proc(1:iproc)
      write(ln,5) 'WHERE Z IS THE MOMENTUM FRACTION '//
     &     'OF THE INITAL EMITTER/SPECTATOR.'
      write(ln,5) ' '

      write(ln,10) 'IMPLICIT NONE'
      write(ln,10) 'include "nexternal.inc" '
      write(ln,10) 'include "dipole.inc" '
      write(ln,10) 'include "coupl.inc"'

      write(ln,5) ' '
      write(ln,5) 'ARGUMENTS '
      write(ln,5) ' '
      write(ln,10) 'REAL*8 P(0:3,NEXTERNAL),X(2),Z,PSWGT'
      write(ln,10) 'REAL*8 EPS,EPSSQ,FIN'
      write(ln,5) ' '
      write(ln,5) 'GLOBAL '
      write(ln,5) ' '
      write(ln,10) 'REAL*8 dot,lambda_tr'
      write(ln,10) 'LOGICAL PASSCUTSDIP'
      write(ln,10) 'EXTERNAL dot,lambda_tr,PASSCUTSDIP'
      write(ln,5) ' '
      write(ln,5) 'LOCAL'
      write(ln,5) ' '
      write(ln,10) 'REAL*8 sij,sik,sjk,zp,e,esq,pdfwgt(3),FINITE(3)'
      write(ln,10) 'REAL*8 y_ijk,omx_ija,u_i,vtilde_i,zplus'
      write(ln,10) 'REAL*8 mass_i,mass_j,mass_k,mass_ij'
      write(ln,10) 'REAL*8 msq_i,msq_j,msq_k,msq_ij,sikzone'
      write(ln,10) 'REAL*8 musq_i,musq_j,musq_k,musq_ij'
      write(ln,10) 'REAL*8 vij,vtildeij,rho_i,rho_k,rho,qsq'
      write(ln,10) 'REAL*8 one,two,three,four,cf,ca,qtot(0:3),pi'
      write(ln,10) 'REAL*8 gsq,gammaE,ff,fi(9)'
      write(ln,10) 'parameter (one=1d0,two=2d0,three=3d0,four=4d0)'
      write(ln,10) 'parameter (cf=4d0/3d0,ca=3d0)'
      write(ln,60) 'REAL*8 ANS(',nrdipol,')'
      write(ln,10) 'REAL*8 ZERO'
      write(ln,10) 'INTEGER NU'
      write(ln,10) 'PARAMETER (ZERO=0D0)'
      write(ln,10) 'PARAMETER (gammaE=0.57721566d0)'
      write(ln,10) 'PARAMETER (pi=3.1415926d0)'
      write(ln,5) ' '
      write(ln,5) 'BEGIN CODE'
      write(ln,5) ' '
      write(ln,10) 'INTDIPOLE=.TRUE.'
      write(ln,10) 'gsq=GG(1)**2'
      write(ln,10) 'eps=0d0'
      write(ln,10) 'epssq=0d0'
      write(ln,10) 'fin=0d0'
      write(ln,10) 'e=0d0'
      write(ln,10) 'esq=0d0'
      write(ln,10) 'zplus=1d0'
      write(ln,10) 'DO NU=1,3'
      write(ln,10) '   finite(NU)=0d0'
      write(ln,10) 'ENDDO'
      write(ln,60) 'DO NU=1,',nrdipol
      write(ln,10) '   ANS(NU)=0d0'
      write(ln,10) 'ENDDO'
      write(ln,10) 'IF (PASSCUTSDIP(P)) THEN'
      do i=1,nr
c     First check if unresolved is massless. If not do print it here.
         if (mass(info_ijk(i,1,2)).ne.'ZERO')then
c**************************************************************
c       final-final theta funtion
c**************************************************************
         if(info_ijk(i,1,1).gt.nincoming
     &        .and. info_ijk(i,1,3).gt.nincoming) then
            write(ln,10)' '
            write(ln,71)'mass_i=',mass(info_ijk(i,1,1))
            write(ln,71)'mass_k=',mass(info_ijk(i,1,3))
            write(ln,10) 'DO NU=0,3'
            write(ln,17) '   qtot(NU)=p(NU,',info_ijk(i,1,4),
     &           ')+p(NU,',info_ijk(i,1,5),')'
            write(ln,10) 'ENDDO'
            write(ln,10)'qsq=dot(qtot,qtot)'

            write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &           info_p(5,inverse(iline_copy(1,i))),',',
     &           info_p(5,inverse(iline_copy(2,i))),',pdfwgt(2))'
            if (i.le.9) then
               write(ln,40)
     &          'CALL SDIPOL'//char(48)//char(48)//char(48+i)
     &              //'(P,ANS(',i,'))'
            elseif(i.le.99) then
               write(ln,40)'CALL SDIPOL'//char(48)//char(i/10+48)
     &              //char(mod(i,10)+48)//'(P,ANS(',i,'))'
            elseif(i.le.999) then
               write(ln,42)'CALL SDIPOL'//char(i/100+48)
     &              //char((i-(i/100)*100)/10+48)
     &              //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
            else
               write (*,*) 'ERROR: write_dipoles: too many dipoles'
            endif
            write(ln,14)
     &           'CALL EPSFF(mass_i,mass_k,qsq,'
     &           ,info_ijk(i,2,4),',',info_ijk(i,2,1),',esq,e)'
            write(ln,11) 'epssq=epssq-esq*ANS(',i,')*pdfwgt(2)'
            write(ln,11) 'eps=eps-e*ANS(',i,')*pdfwgt(2)'
            write(ln,14)
     &           'CALL FINITEFF(mass_i,mass_k,qsq,'
     &           ,info_ijk(i,2,4),',',info_ijk(i,2,1),',fi)'
            write(ln,11)
     &           'finite(2)=finite(2)-fi(2)*ANS(',i,')*pdfwgt(2)'

c**************************************************************
c       final-initial theta function
c**************************************************************
         elseif(info_ijk(i,1,1).gt.nincoming
     &           .and.info_ijk(i,1,3).le.nincoming) then
            write(ln,10)' '
            write(ln,71)'mass_i=',mass(info_ijk(i,1,1))
            write(ln,74)'sik=two*dot(p(0,',info_ijk(i,1,4),'),p(0,'
     &           ,info_ijk(i,1,5),'))/z'
            if(info_ijk(i,2,4).eq.0 .and. (
     &           mass(info_ijk(i,1,1)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,2)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,3)).ne.'ZERO')) then
               write (ln,11) 'zplus=1d0-4d0*mass_i**2/sik'
               write(ln,74)'sikzone=two*dot(p(0,',info_ijk(i,1,4),
     &              '),p(0,',info_ijk(i,1,5),
     &              '))/zplus'
            else
               write(ln,74)'sikzone=two*dot(p(0,',info_ijk(i,1,4),
     &              '),p(0,',info_ijk(i,1,5),'))'
            endif
            if(info_ijk(i,1,3).eq.1) then
               write(ln,17)'CALL DIPOLEPDF(x(1)/z,x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(1))'
            elseif(info_ijk(i,1,3).eq.2) then
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2)/z,',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(1))'
            else
               write (ln,*)'Error #145 in write_intdipole'
               write (*,*)'Error #145 in write_intdipole'
            endif
            if(info_ijk(i,2,4).eq.0 .and. (
     &           mass(info_ijk(i,1,1)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,2)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,3)).ne.'ZERO')) then
               if(info_ijk(i,1,3).eq.1) then
                  write(ln,17)'CALL DIPOLEPDF(x(1)/zplus,x(2),',
     &                 info_p(5,inverse(iline_copy(1,i))),',',
     &                 info_p(5,inverse(iline_copy(2,i))),',pdfwgt(3))'
               elseif(info_ijk(i,1,3).eq.2) then
                  write(ln,17)'CALL DIPOLEPDF(x(1),x(2)/zplus,',
     &                 info_p(5,inverse(iline_copy(1,i))),',',
     &                 info_p(5,inverse(iline_copy(2,i))),',pdfwgt(3))'
               endif
            else
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(2))'
            endif
            if (i.le.9) then
               write(ln,40)'CALL SDIPOL'//char(48)//char(48)
     &              //char(48+i)//'(P,ANS(',i,'))'
            elseif(i.le.99) then
               write(ln,40)'CALL SDIPOL'//char(48)//char(i/10+48)
     &              //char(mod(i,10)+48)//'(P,ANS(',i,'))'
            elseif(i.le.999) then
               write(ln,42)'CALL SDIPOL'//char(i/100+48)
     &              //char((i-(i/100)*100)/10+48)
     &              //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
            else
               write (*,*) 'ERROR: write_dipoles: too many dipoles'
            endif
            write(ln,14)
     &           'CALL EPSFI(',info_ijk(i,2,4),',',
     &           info_ijk(i,2,1),',sikzone,mass_i,esq,e)'
            write(ln,11) 'epssq=epssq-esq*ANS(',i,')*pdfwgt(2)'
            write(ln,11) 'eps=eps-e*ANS(',i,')*pdfwgt(2)'
            write(ln,14) 'CALL FINITEFI(mass_i,sik,sikzone,z,'
     &           ,info_ijk(i,2,4),',',info_ijk(i,2,1),',fi)'

            if(info_ijk(i,2,4).eq.0 .and. (
     &           mass(info_ijk(i,1,1)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,2)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,3)).ne.'ZERO')) then
               write(ln,11)
     &            'finite(2)=finite(2)-(-fi(8))*ANS(',i,')*pdfwgt(2)'
               write(ln,11) 'finite(3)=finite(3)-(fi(7)+fi(9))'//
     &              '*ANS(',i,')*pdfwgt(3)/zplus'
            else
               write(ln,11) 'finite(1)=finite(1)-(fi(1)+fi(3))'//
     &              '*ANS(',i,')*pdfwgt(1)'
               write(ln,11) 'finite(2)=finite(2)-(fi(2)-fi(4))'//
     &              '*ANS(',i,')*pdfwgt(2)'
            endif

c**************************************************************
c       initial-final theta function
c**************************************************************
         elseif(info_ijk(i,1,1).le.nincoming .and.
     &           info_ijk(i,1,3).gt.nincoming) then
            write(ln,10)' '
            write(ln,71)'mass_k=',mass(info_ijk(i,1,3))
            write(ln,74)'sik=two*dot(p(0,',info_ijk(i,1,4)
     &           ,'),p(0,',info_ijk(i,1,5),'))/z'
            write(ln,74)'sikzone=two*dot(p(0,',info_ijk(i,1,4)
     &           ,'),p(0,',info_ijk(i,1,5),'))'
            if (info_ijk(i,1,1).eq.1) then
               write(ln,17)'CALL DIPOLEPDF(x(1)/z,x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(1))'
            elseif(info_ijk(i,1,1).eq.2) then
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2)/z,',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(1))'
            else
               write (ln,*)'Error #146 in write_intdipole'
               write (*,*)'Error #146 in write_intdipole'
            endif
            if(((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.1).or.
     &           (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.0))) then
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(2))'
            endif
            if (i.le.9) then
               write(ln,41)'CALL SDIPOL'//char(48)//char(48)
     &              //char(48+i)//'(P,ANS(',i,'))'
            elseif(i.le.99) then
               write(ln,41)'CALL SDIPOL'//char(48)//char(i/10+48)
     &              //char(mod(i,10)+48)//'(P,ANS(',i,'))'
            elseif(i.le.999) then
               write(ln,41)'CALL SDIPOL'//char(i/100+48)
     &              //char((i-(i/100)*100)/10+48)
     &              //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
            else
               write (*,*) 'ERROR: write_dipoles: too many dipoles'
            endif

            write(ln,14) 'CALL EPSIF(',info_ijk(i,2,1),','
     &           ,info_ijk(i,2,4),',sikzone,mass_k,esq,e)'
            write(ln,11) 'epssq=epssq-esq*ANS(',i,')*pdfwgt(2)'
            write(ln,11) 'eps=eps-e*ANS(',i,')*pdfwgt(2)'
            write(ln,14)
     &           'CALL FINITEIF(mass_k,sik,sikzone,z,'
     &           ,info_ijk(i,2,1),',',info_ijk(i,2,4),',fi)'
c
c quark-quark splitting (i:quark, ij~: quark) and
c gluon-gluon splitting (i:gluon, ij~: gluon)
c
            if(((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.1).or.
     &           (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.0))) then
               if (mass(info_ijk(i,1,3)).eq.'ZERO') then
                  write(ln,11) 'finite(1)=finite(1)-(fi(1)+fi(3))'//
     &                 '*ANS(',i,')*pdfwgt(1)'
                  write(ln,11) 'finite(2)=finite(2)-(fi(2)-fi(4))'//
     &                 '*ANS(',i,')*pdfwgt(2)'
               else
                  write(ln,11) 'finite(1)=finite(1)-'//
     &                 '(fi(1)+fi(3)+fi(5))*ANS(',i,')*pdfwgt(1)'
                  write(ln,11) 'finite(2)=finite(2)-'//
     &                 '(fi(2)-fi(4)-fi(6))*ANS(',i,')*pdfwgt(2)'
               endif
c
c quark-gluon splitting (i:quark, ij~: gluon) and
c gluon-quark splitting (i:gluon, ij~: gluon)
c
            elseif((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.0) .or.
     &             (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.1)) then
             write(ln,11) 'finite(1)=finite(1)-fi(1)'//
     &              '*ANS(',i,')*pdfwgt(1)'
            else
               write (*,*) 'Error #142 in dipolcolor'
               write (ln,*) 'Error #142 in dipolcolor'
            endif

c**************************************************************
c     initial-initial theta function
c**************************************************************
         elseif(info_ijk(i,1,1).le.nincoming .and.
     &           info_ijk(i,1,3) .le.nincoming) then
            write(ln,10)' '
            write(ln,74)'sik=two*dot(p(0,',info_ijk(i,1,4),
     &           '),p(0,',info_ijk(i,1,5),'))/z'
            write(ln,74)'sikzone=two*dot(p(0,',info_ijk(i,1,4),
     &           '),p(0,',info_ijk(i,1,5),'))'
            if (info_ijk(i,1,1).eq.1) then
               write(ln,17)'CALL DIPOLEPDF(x(1)/z,x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(1))'
            elseif(info_ijk(i,1,1).eq.2) then
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2)/z,',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(1))'
            else
               write (ln,*)'Error #147 in write_intdipole'
               write (*,*)'Error #147 in write_intdipole'
            endif
            if(((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.1).or.
     &         (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.0))) then
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(2))'
            endif
            if (i.le.9) then
               write(ln,41)'CALL SDIPOL'//char(48)//char(48)
     &           //char(48+i)//'(P,ANS(',i,'))'
            elseif(i.le.99) then
               write(ln,41)'CALL SDIPOL'//char(48)//char(i/10+48)
     &           //char(mod(i,10)+48)//'(P,ANS(',i,'))'
            elseif(i.le.999) then
               write(ln,41)'CALL SDIPOL'//char(i/100+48)
     &           //char((i-(i/100)*100)/10+48)
     &           //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
            else
               write (*,*) 'ERROR: write_dipoles: too many dipoles'
            endif
            write(ln,14) 'CALL EPSII(',info_ijk(i,2,1),','
     &           ,info_ijk(i,2,4),',sikzone,esq,e)'
            write(ln,11) 'epssq=epssq-esq*ANS(',i,')*pdfwgt(2)'
            write(ln,11) 'eps=eps-e*ANS(',i,')*pdfwgt(2)'
            write(ln,14) 'CALL FINITEII(sik,sikzone,z,'
     &           ,info_ijk(i,2,1),',',info_ijk(i,2,4),',fi)'
c
c quark-quark splitting (i:quark, ij~: quark) and
c gluon-gluon splitting (i:gluon, ij~: gluon)
c
            if(((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.1).or.
     &         (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.0))) then
               if (mass(info_ijk(i,1,3)).eq.'ZERO') then
                  write(ln,11) 'finite(1)=finite(1)-(fi(1)+fi(3))'//
     &                 '*ANS(',i,')*pdfwgt(1)'
                  write(ln,11) 'finite(2)=finite(2)-(fi(2)-fi(4))'//
     &                 '*ANS(',i,')*pdfwgt(2)'
               else
                  write(ln,11) 'finite(1)=finite(1)-'//
     &                 '(fi(1)+fi(3)+fi(5))*ANS(',i,')*pdfwgt(1)'
                  write(ln,11) 'finite(2)=finite(2)-'//
     &                 '(fi(2)-fi(4)-fi(6))*ANS(',i,')*pdfwgt(2)'
               endif
c
c quark-gluon splitting (i:quark, ij~: gluon) and
c gluon-quark splitting (i:gluon, ij~: gluon)
c
            elseif((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.0) .or.
     &             (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.1)) then
               write(ln,11)
     &              'finite(1)=finite(1)-fi(1)*ANS(',i,')*pdfwgt(1)'
            else
               write (*,*) 'Error #144 in dipolcolor'
               write (ln,*) 'Error #144 in dipolcolor'
            endif
         else
            write (*,*) 'Error #143 in dipolcolor'
            write (ln,*) 'Error #143 in dipolcolor'
         endif
         endif ! massless unresolved
      enddo
      write(ln,10) ' '
      write(ln,10) 'ENDIF'
      write(ln,10) 'INTDIPOLE=.FALSE.'
      write(ln,10) 'FIN=(FINITE(1)/z+FINITE(2)+FINITE(3))*PSWGT'
      write(ln,10) 'CALL WRITEHIST(P,FIN)'
      write(ln,10) 'RETURN'
      write(ln,10) 'END'




 5    format('C ',A)
 10   format(6X,A)
 11   format(6X,A,I3,A)
 17   format(6X,A,I3,A,I3,A)
 12   format(5X,'&',A,I3,A)
 13   format(6X,A,I2,A)
 14   format(6X,A,I3,A,I3,A)
 15   format('C',6X,A)
 16   format(6X,A,A,A1,A4,A)
 20   format(', P',I1)
 30   format(',P',I1,'(0:3)')
 40   format(6X,A21,I4,A2)
 41   format(6X,A,I4,A2)
 42   format(6X,A23,I4,A2)
 50   format(7X,A26,I4,A1)
 60   format(6X,A,I4,A1)
 70   format(6X,A12,I2,A6,I2,A2)
 71   format(6X,A,A5)
 72   format(6X,A15,2(I1,A7),I1,A)
 73   format(6X,A15,I1,A7,I1,A)
 74   format(6X,2(A,I1),A)
 75   format(6X,4(A,I1),A)
      close(ln)







      end



      function rli2(x)
c
c real dilogarithm function  0 <= x <= 1
c
      implicit real*8(a-h,o-z)
      parameter(a1 = -0.250000000000000d0)
      parameter(a2 = -0.111111111111111d0)
      parameter(a3 = -0.010000000000000d0)
      parameter(a4 = -0.0170068027210884d0)
      parameter(a5 = -0.0194444444444444d0)
      parameter(a6 = -0.0206611570247934d0)
      parameter(a7 = -0.0214173006480699d0)
      parameter(a8 = -0.02194886637723d0)
      parameter(a9 = -0.0220893589994137d0)
      parameter(a10 = -0.0229303207720760d0)
      parameter(zeta2 =  1.644934066848226d0)
      
      if(x.gt.1.d0)then
        write(*,*)' argument greater than 1 passed to li2'
        rli2=0.d0
        return
      elseif(x.lt.0.d0)then
        write(*,*)' argument less than 0 passed to li2'
        rli2=0.d0
        return
      elseif(x.eq.1.d0)then
        rli2=zeta2
        return
      elseif(x.eq.0.d0)then
        rli2=0d0
        return
      elseif(x.gt.0.5d0)then
        y=1.d0-x
        z=-log(1.d0-y)
        z2=z*z
        rli2=-z*(1.d0+a1*z*(1.d0+a2*z*(1.d0+a3*z2*(1.d0+a4*z2*
     1 (1.d0+a5*z2*(1.d0+a6*z2*(1.d0+a7*z2*(1.d0+a8*z2*(1.d0+a9*z2*
     2 (1.d0+a10*z2))))))))))
     3 +zeta2-log(x)*log(1.d0-x)
       return
      elseif(x.le.0.5d0)then
        y=x
        z=-log(1.d0-y)
        z2=z*z
        rli2=z*(1.d0+a1*z*(1.d0+a2*z*(1.d0+a3*z2*(1.d0+a4*z2*
     1 (1.d0+a5*z2*(1.d0+a6*z2*(1.d0+a7*z2*(1.d0+a8*z2*(1.d0+a9*z2*
     2 (1.d0+a10*z2))))))))))
        return
      endif
      rli2=0d0
      return
      end

      Subroutine removelines(iline_copy,nr)
!***********************************
! Removes dipoles that are identical to dipoles
! already calculated.
!**********************************
      implicit none
      include 'params.inc'

c     Arguments
      integer iline_copy(-maxlines:maxlines,maxdipoles),nr

c    Global
      integer info_ijk(maxdipoles,2,5),nrx
      common /ijk/ info_ijk, nrx

c     Local
      integer i,j,k,l
      logical same,identlines

      identlines=.true.
      do while(identlines)
      identlines=.false.
      do i=1,nr-1
       do j=i+1,nr
        same=.true.
        do k=-maxlines,maxlines
         if(iline_copy(k,i).ne.iline_copy(k,j)) then
           same=.false.
         endif
        enddo
        if((info_ijk(i,1,1).ne.info_ijk(j,1,1)).or.
     &     (info_ijk(i,1,2).ne.info_ijk(j,1,2)).or.
     &     (info_ijk(i,1,3).ne.info_ijk(j,1,3))) then
           same=.false.
        endif
c  if  ij=ji so ji needs to be removed
c        if((info_ijk(i,2,1).eq.0.and. info_ijk(i,2,2).eq.0).and.
c     &     (info_ijk(j,2,1).eq.0.and. info_ijk(j,2,2).eq.0).and.
c     &     (info_ijk(i,1,1).eq.info_ijk(j,1,2)).and.(info_ijk(i,1,2).eq.
c     &     info_ijk(j,1,1)).and.(info_ijk(i,1,3).eq.info_ijk(j,1,3)))then
        if((info_ijk(i,1,1).eq.info_ijk(j,1,2)).and.(info_ijk(i,1,2).eq.
     &     info_ijk(j,1,1)).and.(info_ijk(i,1,3).eq.info_ijk(j,1,3)))then
           same=.true.
        endif

        if(same.eqv..true.) then
         identlines=.true.
         do l=j,nr
          do k=-maxlines,maxlines
           iline_copy(k,l)=iline_copy(k,l+1)
          enddo
          do k=1,5
           info_ijk(l,1,k)=info_ijk(l+1,1,k)
           info_ijk(l,2,k)=info_ijk(l+1,2,k)
          enddo
         enddo
         nr=nr-1
        endif
       enddo
      enddo
      enddo

      end


      Subroutine printdipolcolor(color)
! **************************************
! Prints relevant parts of color matrix
! **************************************
      implicit none

      include 'params.inc'

c    Arguments
      integer color(0:2*maxlines,0:maxfactors,0:maxterms)

c      print*, 'color(0,1,term)'
      print*, color(0,1,1),color(0,1,2),color(0,1,3),color(0,1,4)
      print*, color(1,1,1),color(2,1,1)
      print*, color(1,1,2),color(2,1,2)
      print*, color(1,1,3),color(2,1,3)
      print*, color(1,1,4),color(2,1,4)
c      print*, color(0,2,1),color(0,2,2),color(0,2,3),color(0,2,4)
c      print*,'maxifactor'
c      print*,max_ifactor(1),max_ifactor(2),max_ifactor(3),max_ifactor(4)
      print*, ' '
      print*,color(0,0,0)
      write(*,'(8i5)') color(0,0,1),color(0,0,2),color(0,0,3),color(0,0,4),color(0,0,5)
     &       ,color(0,0,6),color(0,0,7),color(0,0,8)
      write(*,'(8i5)')color(0,2,1),color(0,2,2),color(0,2,3),color(0,2,4),color(0,2,5)
     &       ,color(0,2,6),color(0,2,7),color(0,2,8)
      write(*,'(8i5)')color(1,2,1),color(1,2,2),color(1,2,3),color(1,2,4),color(1,2,5)
     &       ,color(1,2,6),color(1,2,7),color(1,2,8)
      write(*,'(8i5)')color(2,2,1),color(2,2,2),color(2,2,3),color(2,2,4),color(2,2,5)
     &       ,color(2,2,6),color(2,2,7),color(2,2,8)
      write(*,'(8i5)')color(3,2,1),color(3,2,2),color(3,2,3),color(3,2,4),color(3,2,5)
     &       ,color(3,2,6),color(3,2,7),color(3,2,8)
      write(*,'(8i5)')color(4,2,1),color(4,2,2),color(4,2,3),color(4,2,4),color(4,2,5)
     &       ,color(4,2,6),color(4,2,7),color(4,2,8)
      write(*,'(8i5)')color(0,3,1),color(0,3,2),color(0,3,3),color(0,3,4),color(0,3,5)
     &       ,color(0,3,6),color(0,3,7),color(0,3,8)
      write(*,'(8i5)')color(1,3,1),color(1,3,2),color(1,3,3),color(1,3,4),color(1,3,5)
     &       ,color(1,3,6),color(1,3,7),color(1,3,8)
      write(*,'(8i5)')color(2,3,1),color(2,3,2),color(2,3,3),color(2,3,4),color(2,3,5)
     &       ,color(2,3,6),color(2,3,7),color(2,3,8)
      write(*,'(8i5)')color(3,3,1),color(3,3,2),color(3,3,3),color(3,3,4),color(3,3,5)
     &       ,color(3,3,6),color(3,3,7),color(3,3,8)
      write(*,'(8i5)')color(4,3,1),color(4,3,2),color(4,3,3),color(4,3,4),color(4,3,5)
     &       ,color(4,3,6),color(4,3,7),color(4,3,8)
      write(*,'(8i5)')color(0,4,1),color(0,4,2),color(0,4,3),color(0,4,4),color(0,4,5)
     &       ,color(0,4,6),color(0,4,7),color(0,4,8)
      write(*,'(8i5)')color(1,4,1),color(1,4,2),color(1,4,3),color(1,4,4),color(1,4,5)
     &       ,color(1,4,6),color(1,4,7),color(1,4,8)
      write(*,'(8i5)')color(2,4,1),color(2,4,2),color(2,4,3),color(2,4,4),color(2,4,5)
     &       ,color(2,4,6),color(2,4,7),color(2,4,8)
      write(*,'(8i5)')color(3,4,1),color(3,4,2),color(3,4,3),color(3,4,4),color(3,4,5)
     &       ,color(3,4,6),color(3,4,7),color(3,4,8)
      write(*,'(8i5)')color(4,4,1),color(4,4,2),color(4,4,3),color(4,4,4),color(4,4,5)
     &       ,color(4,4,6),color(4,4,7),color(4,4,8)
      write(*,'(8i5)')color(0,5,1),color(0,5,2),color(0,5,3),color(0,5,4),color(0,5,5)
     &       ,color(0,5,6),color(0,5,7),color(0,5,8)
      write(*,'(8i5)')color(1,5,1),color(1,5,2),color(1,5,3),color(1,5,4),color(1,5,5)
     &       ,color(1,5,6),color(1,5,7),color(1,5,8)
      write(*,'(8i5)')color(2,5,1),color(2,5,2),color(2,5,3),color(2,5,4),color(2,5,5)
     &       ,color(2,5,6),color(2,5,7),color(2,5,8)
      write(*,'(8i5)')color(3,5,1),color(3,5,2),color(3,5,3),color(3,5,4),color(3,5,5)
     &       ,color(3,5,6),color(3,5,7),color(3,5,8)
      write(*,'(8i5)')color(4,5,1),color(4,5,2),color(4,5,3),color(4,5,4),color(4,5,5)
     &       ,color(4,5,6),color(4,5,7),color(4,5,8)
      write(*,'(8i5)')color(0,6,1),color(0,6,2),color(0,6,3),color(0,6,4),color(0,6,5)
     &       ,color(0,6,6),color(0,6,7),color(0,6,8)
      write(*,'(8i5)')color(1,6,1),color(1,6,2),color(1,6,3),color(1,6,4),color(1,6,5)
     &       ,color(1,6,6),color(1,6,7),color(1,6,8)
      write(*,'(8i5)')color(2,6,1),color(2,6,2),color(2,6,3),color(2,6,4),color(2,6,5)
     &       ,color(2,6,6),color(2,6,7),color(2,6,8)
      write(*,'(8i5)')color(3,6,1),color(3,6,2),color(3,6,3),color(3,6,4),color(3,6,5)
     &       ,color(3,6,6),color(3,6,7),color(3,6,8)
      write(*,'(8i5)')color(4,6,1),color(4,6,2),color(4,6,3),color(4,6,4),color(4,6,5)
     &       ,color(4,6,6),color(4,6,7),color(4,6,8)
      write(*,'(8i5)')color(0,7,1),color(0,7,2),color(0,7,3),color(0,7,4),color(0,7,5)
     &       ,color(0,7,6),color(0,7,7),color(0,7,8)
      write(*,'(8i5)')color(1,7,1),color(1,7,2),color(1,7,3),color(1,7,4),color(1,7,5)
     &       ,color(1,7,6),color(1,7,7),color(1,7,8)
      write(*,'(8i5)')color(2,7,1),color(2,7,2),color(2,7,3),color(2,7,4),color(2,7,5)
     &       ,color(2,7,6),color(2,7,7),color(2,7,8)
      write(*,'(8i5)')color(3,7,1),color(3,7,2),color(3,7,3),color(3,7,4),color(3,7,5)
     &       ,color(3,7,6),color(3,7,7),color(3,7,8)
      write(*,'(8i5)')color(4,7,1),color(4,7,2),color(4,7,3),color(4,7,4),color(4,7,5)
     &       ,color(4,7,6),color(4,7,7),color(4,7,8)
      write(*,'(8i5)')color(0,8,1),color(0,8,2),color(0,8,3),color(0,8,4),color(0,8,5)
     &       ,color(0,8,6),color(0,8,7),color(0,8,8)
      write(*,'(8i5)')color(1,8,1),color(1,8,2),color(1,8,3),color(1,8,4),color(1,8,5)
     &       ,color(1,8,6),color(1,8,7),color(1,8,8)
      write(*,'(8i5)')color(2,8,1),color(2,8,2),color(2,8,3),color(2,8,4),color(2,8,5)
     &       ,color(2,8,6),color(2,8,7),color(2,8,8)
      write(*,'(8i5)')color(3,8,1),color(3,8,2),color(3,8,3),color(3,8,4),color(3,8,5)
     &       ,color(3,8,6),color(3,8,7),color(3,8,8)
      write(*,'(8i5)')color(4,8,1),color(4,8,2),color(4,8,3),color(4,8,4),color(4,8,5)
     &       ,color(4,8,6),color(4,8,7),color(4,8,8)



      end


      Subroutine write_dipolsqed()
!***************************************
! writes all dipoles in the file dipolsum.f
! and includes the theta function of the
! alpha parameter
!***************************************

      implicit none

      include 'params.inc'

c     Global
      integer dipollist(maxdipoles,2,2),nrdipols
      common /list/ dipollist,nrdipols

      integer info_ijk(maxdipoles,2,5),nr
      common /ijk/ info_ijk,nr

      character*70       cversion
      common/to_version/ cversion

      character*60 proc
      integer iproc
      common/to_write/iproc,proc

      integer         nincoming,nincfirst
      common/to_proc2/nincoming,nincfirst

      character*(8) mass(max_particles)
      common /masses/ mass

      integer iline_copy(-maxlines:maxlines,maxdipoles)
      common /to_writedipole/iline_copy
      character*(4*max_particles) particle(4)
      integer                               charge_c(max_particles)
      integer iparticle(0:max_particles,0:4),inverse(max_particles)
      common/to_model/iparticle,  particle,  inverse, charge_c
      integer info_p(5,max_particles),iposx(3,max_particles)
      character*(max_string) iwave(max_particles),owave(max_particles)
      character*(8) str(3,max_particles)
      common/to_external/iwave,owave,iposx,info_p,str


      integer iline_orig(-maxlines:maxlines)
      common /ciline_orig/iline_orig

c     Local
      integer ln,i,j,ascii_emitter,ascii_spectator,ascii_nr
      integer nrdipol
      parameter (ln=81)
      logical first(max_particles,max_particles)

      do i=1,max_particles
       do j=1,max_particles
        first(i,j)=.true.
       enddo
      enddo

      nrdipol=nr
      if(nrdipol.eq.0) nrdipol=1

      open(unit=ln,file='dipolsumqed.f',status='unknown')

      write(ln,10) 'SUBROUTINE DIPOLSUMQED(P,X,PSWGT,SUBTRACT)'
      write(ln,5) ' '
      write(ln,5) cversion
      write(ln,5) 'RETURNS SUM OVER ALL DIPOLS'
      write(ln,5) 'FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)'
      write(ln,5) ' '
      write(ln,5) 'FOR PROCESS '//proc(1:iproc)
      write(ln,5) ' '

      write(ln,10) 'IMPLICIT NONE'
      write(ln,10) 'include "nexternal.inc" '
      write(ln,10) 'include "dipole.inc" '
      write(ln,10) 'include "coupl.inc"'

      write(ln,5) ' '
      write(ln,5) 'ARGUMENTS '
      write(ln,5) ' '
      write(ln,10) 'REAL*8 P(0:3,NEXTERNAL),X(2),PSWGT,SUBTRACT'
      write(ln,5) ' '
      write(ln,5) 'GLOBAL '
      write(ln,5) ' '
      write(ln,10) 'REAL*8 dot'
      write(ln,10) 'LOGICAL PASSCUTSDIP'
      write(ln,10) 'EXTERNAL dot, PASSCUTSDIP'
      write(ln,5) ' '
      write(ln,5) 'LOCAL'
      write(ln,5) ' '
      write(ln,10) 'REAL*8 sij,sik,sjk,q(0:3),qsq,yp,P1(0:3,NEXTERNAL),pdfwgt'
      write(ln,10) 'REAL*8 y_ijk,omx_ija,u_i,vtilde_i'
      write(ln,10) 'REAL*8 mass_i,mass_j,mass_k,mass_ij'
      write(ln,10) 'REAL*8 msq_i,msq_j,msq_k,msq_ij,mu_i,mu_j,mu_k'
      write(ln,10) 'REAL*8 one,two,three,four'
      write(ln,10) 'parameter (one=1d0,two=2d0,three=3d0,four=4d0)'
      write(ln,60) 'REAL*8 ANS(',nrdipol,')'
      write(ln,10) 'REAL*8 ZERO'
      write(ln,10) 'INTEGER NU'
      write(ln,10) 'PARAMETER (ZERO=0D0)'
      write(ln,5) ' '
      write(ln,5) 'BEGIN CODE'
      write(ln,5) ' '
      write(ln,10) 'INTDIPOLE=.FALSE.'
      write(ln,10) 'SUBTRACT=0d0'
      do i=1,nr
c     First check if unresolved is massless. If not do not print it here.
         if (mass(info_ijk(i,1,2)).eq.'ZERO')then
c       final-final theta funtion
         if(info_ijk(i,1,1).gt.nincoming .and. info_ijk(i,1,3) .gt.nincoming) then
            write(ln,10) ' '
            write(ln,70)'sij=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,2),'))'
            write(ln,70)'sjk=dot(p(0,',info_ijk(i,1,2)
     .           ,'),p(0,',info_ijk(i,1,3),'))'
            write(ln,70)'sik=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,3),'))'
            write(ln,10)'y_ijk=sij/(sij+sjk+sik)'
           if(mass(info_ijk(i,1,3)).ne.'ZERO') then
            write(ln,10)'DO NU=0,3'
            write(ln,73)' q(nu)=p(nu,',info_ijk(i,1,1),')+p(nu,',info_ijk(i,1,2),
     &                    ')+p(nu,',info_ijk(i,1,3),')'
            write(ln,10)'ENDDO'
            write(ln,10)'qsq=dot(q,q)'
            write(ln,71)'mass_i=',mass(info_ijk(i,1,1))
            write(ln,71)'mass_j=',mass(info_ijk(i,1,2))
            write(ln,71)'mass_k=',mass(info_ijk(i,1,3))
            write(ln,10)'mu_i=mass_i/SQRT(qsq)'
            write(ln,10)'mu_j=mass_j/SQRT(qsq)'
            write(ln,10)'mu_k=mass_k/SQRT(qsq)'
            write(ln,10)'yp=1d0-2d0*mu_k*(1d0-mu_k)/(1d0-mu_i**2-mu_j**2-mu_k**2)'
            write(ln,14)'CALL TRANSFORMQED(P,P1,',info_ijk(i,1,1),',',info_ijk(i,1,2),',',
     &                   info_ijk(i,1,3),',',info_ijk(i,1,4),',',info_ijk(i,1,5),',',
     &                   mass(info_ijk(i,1,1)),',',mass(info_ijk(i,1,2)),',',mass(info_ijk(i,1,3)),')'
            write(ln,10)'IF(PASSCUTSDIP(P1)) THEN'
            write(ln,10)'IF (alpha_ff*yp .gt.y_ijk) THEN'
           else
            write(ln,14)'CALL TRANSFORMQED(P,P1,',info_ijk(i,1,1),',',info_ijk(i,1,2),',',
     &                   info_ijk(i,1,3),',',info_ijk(i,1,4),',',info_ijk(i,1,5),',',
     &                   mass(info_ijk(i,1,1)),',',mass(info_ijk(i,1,2)),',',mass(info_ijk(i,1,3)),')'
            write(ln,10)'IF(PASSCUTSDIP(P1)) THEN'
            write(ln,10)'IF (alpha_ff .gt.y_ijk) THEN'
           endif
            write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &           info_p(5,inverse(iline_orig(1))),',',
     &           info_p(5,inverse(iline_orig(2))),',pdfwgt)'
            if (i.le.9) then
               write(ln,40)
     .          'CALL SDIPOLQED'//char(48)//char(48)//char(48+i)//'(P,ANS(',i,'))'
            elseif(i.le.99) then
               write(ln,40)'CALL SDIPOLQED'//char(48)//char(i/10+48)
     &              //char(mod(i,10)+48)//'(P,ANS(',i,'))'
            elseif(i.le.999) then
               write(ln,41)'CALL SDIPOLQED'//char(i/100+48)
     &              //char((i-(i/100)*100)/10+48)
     &              //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
            else
               write (*,*) 'ERROR: write_dipoles: too many dipoles'
            endif
            if(info_ijk(i,2,1).eq.info_ijk(i,2,2).and..not.
     &          first(info_ijk(i,1,1),info_ijk(i,1,2))) then
             write(ln,52)'SUBTRACT = SUBTRACT - ANS(',i,')*PDFWGT'
             write(ln,52)'CALL WRITEHIST(P1,ANS(',i,')*PSWGT*PDFWGT)'
            else
             write(ln,51)'SUBTRACT = SUBTRACT - ANS(',i,')*PDFWGT'
             write(ln,51)'CALL WRITEHIST(P1,ANS(',i,')*PSWGT*PDFWGT)'
             first(info_ijk(i,1,1),info_ijk(i,1,2))=.false.
            endif
            write(ln,10)'ENDIF'
            write(ln,10)'ENDIF'
         endif
c       final-initial theta function
         if(info_ijk(i,1,1).gt.nincoming .and. info_ijk(i,1,3) .le.nincoming) then
            write(ln,10) ' '
            write(ln,70)'sij=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,2),'))'
            write(ln,70)'sjk=dot(p(0,',info_ijk(i,1,2)
     .           ,'),p(0,',info_ijk(i,1,3),'))'
            write(ln,70)'sik=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,3),'))'
            if (mass(info_ijk(i,1,1)).ne.'ZERO'.or.
     .          mass(info_ijk(i,1,2)).ne.'ZERO') then
               write(ln,71)'mass_i=',mass(info_ijk(i,1,1))
               write(ln,71)'mass_j=',mass(info_ijk(i,1,2))
c               write(ln,71)'mass_k=',mass(info_ijk(i,1,3))
               write(ln,10)'msq_i=mass_i**2'
               write(ln,10)'msq_j=mass_j**2'
c            write(ln,10)'msq_k=mass_k**2'
               write(ln,10)'if(mass_i .eq.mass_j) then'
               write(ln,10)' msq_ij=0d0'
               write(ln,10)'else'
               write(ln,10)' msq_ij=msq_i'
               write(ln,10)'endif'
               write(ln,10)'omx_ija='//
     &              '(sij-0.5d0*(msq_ij-msq_i-msq_j))/(sjk+sik)'
            else
               write(ln,10)'omx_ija=sij/(sjk+sik)'
            endif
            write(ln,14)'CALL TRANSFORMQED(P,P1,',info_ijk(i,1,1),',',info_ijk(i,1,2),',',
     &                   info_ijk(i,1,3),',',info_ijk(i,1,4),',',info_ijk(i,1,5),',',
     &                   mass(info_ijk(i,1,1)),',',mass(info_ijk(i,1,2)),',',mass(info_ijk(i,1,3)),')'
            write(ln,10)'IF(PASSCUTSDIP(P1)) THEN'
            write(ln,10)'IF (alpha_fi .gt.omx_ija) THEN'
            write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &           info_p(5,inverse(iline_orig(1))),',',
     &           info_p(5,inverse(iline_orig(2))),',pdfwgt)'
            if (i.le.9) then
               write(ln,40)'CALL SDIPOLQED'//char(48)//char(48)
     .              //char(48+i)//'(P,ANS(',i,'))'
            elseif(i.le.99) then
               write(ln,40)'CALL SDIPOLQED'//char(48)//char(i/10+48)
     &              //char(mod(i,10)+48)//'(P,ANS(',i,'))'
            elseif(i.le.999) then
               write(ln,41)'CALL SDIPOLQED'//char(i/100+48)
     &              //char((i-(i/100)*100)/10+48)
     &              //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
            else
               write (*,*) 'ERROR: write_dipoles: too many dipoles'
            endif
            if(info_ijk(i,2,1).eq.info_ijk(i,2,2).and..not.
     &         first(info_ijk(i,1,1),info_ijk(i,1,2))) then
             write(ln,52)'SUBTRACT = SUBTRACT - ANS(',i,')*PDFWGT'
             write(ln,52)'CALL WRITEHIST(P1,ANS(',i,')*PSWGT*PDFWGT)'
            else
             write(ln,51)'SUBTRACT = SUBTRACT - ANS(',i,')*PDFWGT'
             write(ln,51)'CALL WRITEHIST(P1,ANS(',i,')*PSWGT*PDFWGT)'
             first(info_ijk(i,1,1),info_ijk(i,1,2))=.false.
            endif
            write(ln,10)'ENDIF'
            write(ln,10)'ENDIF'
         endif
c       initial-final theta function
         if(info_ijk(i,1,1).le.nincoming .and. info_ijk(i,1,3) .gt.nincoming) then
            write(ln,10) ' '
            write(ln,70)'sij=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,2),'))'
            write(ln,70)'sik=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,3),'))'
            write(ln,10)'u_i=sij/(sij+sik)'
            write(ln,14)'CALL TRANSFORMQED(P,P1,',info_ijk(i,1,1),',',info_ijk(i,1,2),',',
     &                   info_ijk(i,1,3),',',info_ijk(i,1,4),',',info_ijk(i,1,5),',',
     &                   mass(info_ijk(i,1,1)),',',mass(info_ijk(i,1,2)),',',mass(info_ijk(i,1,3)),')'
            write(ln,10)'IF(PASSCUTSDIP(P1)) THEN'
            write(ln,10)'IF (alpha_if .gt.u_i) THEN'
            write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &           info_p(5,inverse(iline_orig(1))),',',
     &           info_p(5,inverse(iline_orig(2))),',pdfwgt)'
            if (i.le.9) then
               write(ln,40)'CALL SDIPOLQED'//char(48)//char(48)
     .              //char(48+i)//'(P,ANS(',i,'))'
            elseif(i.le.99) then
               write(ln,40)'CALL SDIPOLQED'//char(48)//char(i/10+48)
     &              //char(mod(i,10)+48)//'(P,ANS(',i,'))'
            elseif(i.le.999) then
               write(ln,41)'CALL SDIPOLQED'//char(i/100+48)
     &              //char((i-(i/100)*100)/10+48)
     &              //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
            else
               write (*,*) 'ERROR: write_dipoles: too many dipoles'
            endif
             if((info_ijk(i,2,1).eq.info_ijk(i,2,2).and..not.
     &           first(info_ijk(i,1,1),info_ijk(i,1,2))).or.
     &          (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,2).eq.1.and..not.
     &          first(info_ijk(i,1,1),info_ijk(i,1,2)))) then
             write(ln,52)'SUBTRACT = SUBTRACT - ANS(',i,')*PDFWGT'
             write(ln,52)'CALL WRITEHIST(P1,ANS(',i,')*PSWGT*PDFWGT)'
            else
             write(ln,51)'SUBTRACT = SUBTRACT - ANS(',i,')*PDFWGT'
             write(ln,51)'CALL WRITEHIST(P1,ANS(',i,')*PSWGT*PDFWGT)'
             first(info_ijk(i,1,1),info_ijk(i,1,2))=.false.
            endif
            write(ln,10)'ENDIF'
            write(ln,10)'ENDIF'
         endif
c     initial-initial theta function
         if(info_ijk(i,1,1).le.nincoming .and. info_ijk(i,1,3) .le.nincoming) then
            write(ln,10) ' '
            write(ln,70)'sij=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,2),'))'
            write(ln,70)'sik=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,3),'))'
            write(ln,10)'vtilde_i=sij/sik'
            write(ln,14)'CALL TRANSFORMQED(P,P1,',info_ijk(i,1,1),',',info_ijk(i,1,2),',',
     &                   info_ijk(i,1,3),',',info_ijk(i,1,4),',',info_ijk(i,1,5),',',
     &                   mass(info_ijk(i,1,1)),',',mass(info_ijk(i,1,2)),',',mass(info_ijk(i,1,3)),')'
            write(ln,10)'IF(PASSCUTSDIP(P1)) THEN'
            write(ln,10)'IF (alpha_ii .gt.vtilde_i) THEN'
            write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &           info_p(5,inverse(iline_orig(1))),',',
     &           info_p(5,inverse(iline_orig(2))),',pdfwgt)'
            if (i.le.9) then
               write(ln,40)'CALL SDIPOLQED'//char(48)//char(48)
     .              //char(48+i)//'(P,ANS(',i,'))'
            elseif(i.le.99) then
               write(ln,40)'CALL SDIPOLQED'//char(48)//char(i/10+48)
     &              //char(mod(i,10)+48)//'(P,ANS(',i,'))'
            elseif(i.le.999) then
               write(ln,41)'CALL SDIPOLQED'//char(i/100+48)
     &              //char((i-(i/100)*100)/10+48)
     &              //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
            else
               write (*,*) 'ERROR: write_dipoles: too many dipoles'
            endif
             if((info_ijk(i,2,1).eq.info_ijk(i,2,2).and..not.
     &          first(info_ijk(i,1,1),info_ijk(i,1,2))).or.
     &          (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,2).eq.1.and..not.
     &          first(info_ijk(i,1,1),info_ijk(i,1,2)))) then
             write(ln,52)'SUBTRACT = SUBTRACT - ANS(',i,')*PDFWGT'
             write(ln,52)'CALL WRITEHIST(P1,ANS(',i,')*PSWGT*PDFWGT)'
            else
             write(ln,51)'SUBTRACT = SUBTRACT - ANS(',i,')*PDFWGT'
             write(ln,51)'CALL WRITEHIST(P1,ANS(',i,')*PSWGT*PDFWGT)'
             first(info_ijk(i,1,1),info_ijk(i,1,2))=.false.
            endif
            write(ln,10)'ENDIF'
            write(ln,10)'ENDIF'
         endif
      endif
      enddo

      write(ln,10) 'RETURN'
      write(ln,10) 'END'
      write(ln,5) ' '
      write(ln,5) ' '
      write(ln,5) ' '
      write(ln,5) ' '
      write(ln,5) ' '

      write(ln,10) 'SUBROUTINE DIPOLSUMFINITEQED(P,X,PSWGT,SUBTRACT)'
      write(ln,5) ' '
      write(ln,5) cversion
      write(ln,5) 'RETURNS SUM OVER "NON-DIVERGENT" DIPOLS,'
      write(ln,5) 'I.E. DIPOLES FOR WHICH THE UNRESOLVED IS'
      write(ln,5) 'MASSIVE (IN particles.dat MASS IS NOT'
      write(ln,5) 'EQUAL TO "ZERO")'
      write(ln,5) 'FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)'
      write(ln,5) ' '
      write(ln,5) 'FOR PROCESS '//proc(1:iproc)
      write(ln,5) ' '

      write(ln,10) 'IMPLICIT NONE'
      write(ln,10) 'include "nexternal.inc" '
      write(ln,10) 'include "dipole.inc" '
      write(ln,10) 'include "coupl.inc"'

      write(ln,5) ' '
      write(ln,5) 'ARGUMENTS '
      write(ln,5) ' '
      write(ln,10) 'REAL*8 P(0:3,NEXTERNAL),X(2),PSWGT,SUBTRACT'
      write(ln,5) ' '
      write(ln,5) 'GLOBAL '
      write(ln,5) ' '
      write(ln,10) 'REAL*8 dot'
      write(ln,10) 'LOGICAL PASSCUTSDIP'
      write(ln,10) 'EXTERNAL dot, PASSCUTSDIP'
      write(ln,5) ' '
      write(ln,5) 'LOCAL'
      write(ln,5) ' '
      write(ln,10) 'REAL*8 sij,sik,sjk,q(0:3),qsq,yp,P1(0:3,NEXTERNAL),pdfwgt'
      write(ln,10) 'REAL*8 y_ijk,omx_ija,u_i,vtilde_i'
      write(ln,10) 'REAL*8 mass_i,mass_j,mass_k,mass_ij'
      write(ln,10) 'REAL*8 msq_i,msq_j,msq_k,msq_ij,mu_i,mu_j,mu_k'
      write(ln,10) 'REAL*8 one,two,three,four'
      write(ln,10) 'parameter (one=1d0,two=2d0,three=3d0,four=4d0)'
      write(ln,60) 'REAL*8 ANS(',nrdipol,')'
      write(ln,10) 'REAL*8 ZERO'
      write(ln,10) 'INTEGER NU'
      write(ln,10) 'PARAMETER (ZERO=0D0)'
      write(ln,5) ' '
      write(ln,5) 'BEGIN CODE'
      write(ln,5) ' '
      write(ln,10) 'INTDIPOLE=.FALSE.'
      write(ln,10) 'SUBTRACT=0d0'
      do i=1,nr
         if (mass(info_ijk(i,1,2)).ne.'ZERO')then
c       final-final theta funtion
         if(info_ijk(i,1,1).gt.nincoming .and. info_ijk(i,1,3) .gt.nincoming) then
            write(ln,70)'sij=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,2),'))'
            write(ln,70)'sjk=dot(p(0,',info_ijk(i,1,2)
     .           ,'),p(0,',info_ijk(i,1,3),'))'
            write(ln,70)'sik=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,3),'))'
            write(ln,10)'y_ijk=sij/(sij+sjk+sik)'
           if(mass(info_ijk(i,1,3)).ne.'ZERO') then
            write(ln,10)'DO NU=0,3'
            write(ln,73)' q(nu)=p(nu,',info_ijk(i,1,1),')+p(nu,',info_ijk(i,1,2),
     &                    ')+p(nu,',info_ijk(i,1,3),')'
            write(ln,10)'ENDDO'
            write(ln,10)'qsq=dot(q,q)'
            write(ln,71)'mass_i=',mass(info_ijk(i,1,1))
            write(ln,71)'mass_j=',mass(info_ijk(i,1,2))
            write(ln,71)'mass_k=',mass(info_ijk(i,1,3))
            write(ln,10)'mu_i=mass_i/SQRT(qsq)'
            write(ln,10)'mu_j=mass_j/SQRT(qsq)'
            write(ln,10)'mu_k=mass_k/SQRT(qsq)'
            write(ln,10)'yp=1d0-2d0*mu_k*(1d0-mu_k)/(1d0-mu_i**2-mu_j**2-mu_k**2)'
           write(ln,14)'CALL TRANSFORMQED(P,P1,',info_ijk(i,1,1),',',info_ijk(i,1,2),',',
     &                   info_ijk(i,1,3),',',info_ijk(i,1,4),',',info_ijk(i,1,5),',',
     &                   mass(info_ijk(i,1,1)),',',mass(info_ijk(i,1,2)),',',mass(info_ijk(i,1,3)),')'
            write(ln,10)'IF(PASSCUTSDIP(P1)) THEN'
            write(ln,10)'IF (alpha_ff*yp .gt.y_ijk) THEN'
           else
            write(ln,14)'CALL TRANSFORMQED(P,P1,',info_ijk(i,1,1),',',info_ijk(i,1,2),',',
     &                   info_ijk(i,1,3),',',info_ijk(i,1,4),',',info_ijk(i,1,5),',',
     &                   mass(info_ijk(i,1,1)),',',mass(info_ijk(i,1,2)),',',mass(info_ijk(i,1,3)),')'
            write(ln,10)'IF(PASSCUTSDIP(P1)) THEN'
            write(ln,10)'IF (alpha_ff .gt.y_ijk) THEN'
           endif
            write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &           info_p(5,inverse(iline_orig(1))),',',
     &           info_p(5,inverse(iline_orig(2))),',pdfwgt)'
            if (i.le.9) then
               write(ln,40)
     .          'CALL SDIPOLQED'//char(48)//char(48)//char(48+i)//'(P,ANS(',i,'))'
            elseif(i.le.99) then
               write(ln,40)'CALL SDIPOLQED'//char(48)//char(i/10+48)
     &              //char(mod(i,10)+48)//'(P,ANS(',i,'))'
            elseif(i.le.999) then
               write(ln,41)'CALL SDIPOLQED'//char(i/100+48)
     &              //char((i-(i/100)*100)/10+48)
     &              //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
            else
               write (*,*) 'ERROR: write_dipoles: too many dipoles'
            endif
            if(info_ijk(i,2,1).eq.info_ijk(i,2,2).and..not.
     &          first(info_ijk(i,1,1),info_ijk(i,1,2))) then
             write(ln,52)'SUBTRACT = SUBTRACT - ANS(',i,')*PDFWGT'
             write(ln,52)'CALL WRITEHIST(P1,ANS(',i,')*PSWGT*PDFWGT)'
            else
             write(ln,51)'SUBTRACT = SUBTRACT - ANS(',i,')*PDFWGT'
             write(ln,51)'CALL WRITEHIST(P1,ANS(',i,')*PSWGT*PDFWGT)'
             first(info_ijk(i,1,1),info_ijk(i,1,2))=.false.
            endif
            write(ln,10)'ENDIF'
            write(ln,10)'ENDIF'
         endif
c       final-initial theta function
         if(info_ijk(i,1,1).gt.nincoming .and. info_ijk(i,1,3) .le.nincoming) then
            write(ln,70)'sij=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,2),'))'
            write(ln,70)'sjk=dot(p(0,',info_ijk(i,1,2)
     .           ,'),p(0,',info_ijk(i,1,3),'))'
            write(ln,70)'sik=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,3),'))'
            if (mass(info_ijk(i,1,1)).ne.'ZERO'.or.
     .          mass(info_ijk(i,1,2)).ne.'ZERO') then
               write(ln,71)'mass_i=',mass(info_ijk(i,1,1))
               write(ln,71)'mass_j=',mass(info_ijk(i,1,2))
c               write(ln,71)'mass_k=',mass(info_ijk(i,1,3))
               write(ln,10)'msq_i=mass_i**2'
               write(ln,10)'msq_j=mass_j**2'
c            write(ln,10)'msq_k=mass_k**2'
               write(ln,10)'if(mass_i .eq.mass_j) then'
               write(ln,10)' msq_ij=0d0'
               write(ln,10)'else'
               write(ln,10)' msq_ij=msq_i'
               write(ln,10)'endif'
               write(ln,10)'omx_ija=(sij-0.5d0*(msq_ij-msq_i-msq_j))/(sjk+sik)'
            else
               write(ln,10)'omx_ija=sij/(sjk+sik)'
            endif
            write(ln,14)'CALL TRANSFORMQED(P,P1,',info_ijk(i,1,1),',',info_ijk(i,1,2),',',
     &                   info_ijk(i,1,3),',',info_ijk(i,1,4),',',info_ijk(i,1,5),',',
     &                   mass(info_ijk(i,1,1)),',',mass(info_ijk(i,1,2)),',',mass(info_ijk(i,1,3)),')'
            write(ln,10)'IF(PASSCUTSDIP(P1)) THEN'
            write(ln,10)'IF (alpha_fi .gt.omx_ija) THEN'
            write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &           info_p(5,inverse(iline_orig(1))),',',
     &           info_p(5,inverse(iline_orig(2))),',pdfwgt)'
            if (i.le.9) then
               write(ln,40)'CALL SDIPOLQED'//char(48)//char(48)
     .              //char(48+i)//'(P,ANS(',i,'))'
            elseif(i.le.99) then
               write(ln,40)'CALL SDIPOLQED'//char(48)//char(i/10+48)
     &              //char(mod(i,10)+48)//'(P,ANS(',i,'))'
            elseif(i.le.999) then
               write(ln,41)'CALL SDIPOLQED'//char(i/100+48)
     &              //char((i-(i/100)*100)/10+48)
     &              //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
            else
               write (*,*) 'ERROR: write_dipoles: too many dipoles'
            endif
            if(info_ijk(i,2,1).eq.info_ijk(i,2,2).and..not.
     &          first(info_ijk(i,1,1),info_ijk(i,1,2))) then
             write(ln,52)'SUBTRACT = SUBTRACT - ANS(',i,')*PDFWGT'
             write(ln,52)'CALL WRITEHIST(P1,ANS(',i,')*PSWGT*PDFWGT)'
            else
             write(ln,51)'SUBTRACT = SUBTRACT - ANS(',i,')*PDFWGT'
             write(ln,51)'CALL WRITEHIST(P1,ANS(',i,')*PSWGT*PDFWGT)'
             first(info_ijk(i,1,1),info_ijk(i,1,2))=.false.
            endif
            write(ln,10)'ENDIF'
            write(ln,10)'ENDIF'
         endif
c       initial-final theta function
         if(info_ijk(i,1,1).le.nincoming .and. info_ijk(i,1,3) .gt.nincoming) then
            write(ln,70)'sij=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,2),'))'
            write(ln,70)'sik=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,3),'))'
            write(ln,10)'u_i=sij/(sij+sik)'
            write(ln,14)'CALL TRANSFORMQED(P,P1,',info_ijk(i,1,1),',',info_ijk(i,1,2),',',
     &                   info_ijk(i,1,3),',',info_ijk(i,1,4),',',info_ijk(i,1,5),',',
     &                   mass(info_ijk(i,1,1)),',',mass(info_ijk(i,1,2)),',',mass(info_ijk(i,1,3)),')'
            write(ln,10)'IF(PASSCUTSDIP(P1)) THEN'
            write(ln,10)'IF (alpha_if .gt.u_i) THEN'
            write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &           info_p(5,inverse(iline_orig(1))),',',
     &           info_p(5,inverse(iline_orig(2))),',pdfwgt)'
            if (i.le.9) then
               write(ln,40)'CALL SDIPOLQED'//char(48)//char(48)
     .              //char(48+i)//'(P,ANS(',i,'))'
            elseif(i.le.99) then
               write(ln,40)'CALL SDIPOLQED'//char(48)//char(i/10+48)
     &              //char(mod(i,10)+48)//'(P,ANS(',i,'))'
            elseif(i.le.999) then
               write(ln,41)'CALL SDIPOLQED'//char(i/100+48)
     &              //char((i-(i/100)*100)/10+48)
     &              //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
            else
               write (*,*) 'ERROR: write_dipoles: too many dipoles'
            endif
            if(info_ijk(i,2,1).eq.info_ijk(i,2,2).and..not.
     &          first(info_ijk(i,1,1),info_ijk(i,1,2))) then
             write(ln,52)'SUBTRACT = SUBTRACT - ANS(',i,')*PDFWGT'
             write(ln,52)'CALL WRITEHIST(P1,ANS(',i,')*PSWGT*PDFWGT)'
            else
             write(ln,51)'SUBTRACT = SUBTRACT - ANS(',i,')*PDFWGT'
             write(ln,51)'CALL WRITEHIST(P1,ANS(',i,')*PSWGT*PDFWGT)'
             first(info_ijk(i,1,1),info_ijk(i,1,2))=.false.
            endif
            write(ln,10)'ENDIF'
            write(ln,10)'ENDIF'
         endif
c     initial-initial theta function
         if(info_ijk(i,1,1).le.nincoming .and. info_ijk(i,1,3) .le.nincoming) then
            write(ln,70)'sij=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,2),'))'
            write(ln,70)'sik=dot(p(0,',info_ijk(i,1,1)
     .           ,'),p(0,',info_ijk(i,1,3),'))'
            write(ln,10)'vtilde_i=sij/sik'
            write(ln,14)'CALL TRANSFORMQED(P,P1,',info_ijk(i,1,1),',',info_ijk(i,1,2),',',
     &                   info_ijk(i,1,3),',',info_ijk(i,1,4),',',info_ijk(i,1,5),',',
     &                   mass(info_ijk(i,1,1)),',',mass(info_ijk(i,1,2)),',',mass(info_ijk(i,1,3)),')'
            write(ln,10)'IF(PASSCUTSDIP(P1)) THEN'
            write(ln,10)'IF (alpha_ii .gt.vtilde_i) THEN'
            if (i.le.9) then
               write(ln,40)'CALL SDIPOLQED'//char(48)//char(48)
     .              //char(48+i)//'(P,ANS(',i,'))'
            elseif(i.le.99) then
               write(ln,40)'CALL SDIPOLQED'//char(48)//char(i/10+48)
     &              //char(mod(i,10)+48)//'(P,ANS(',i,'))'
            elseif(i.le.999) then
               write(ln,41)'CALL SDIPOLQED'//char(i/100+48)
     &              //char((i-(i/100)*100)/10+48)
     &              //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
            else
               write (*,*) 'ERROR: write_dipoles: too many dipoles'
            endif
            if(info_ijk(i,2,1).eq.info_ijk(i,2,2).and..not.
     &          first(info_ijk(i,1,1),info_ijk(i,1,2))) then
             write(ln,52)'SUBTRACT = SUBTRACT - ANS(',i,')*PDFWGT'
             write(ln,52)'CALL WRITEHIST(P1,ANS(',i,')*PSWGT*PDFWGT)'
            else
             write(ln,51)'SUBTRACT = SUBTRACT - ANS(',i,')*PDFWGT'
             write(ln,51)'CALL WRITEHIST(P1,ANS(',i,')*PSWGT*PDFWGT)'
             first(info_ijk(i,1,1),info_ijk(i,1,2))=.false.
            endif
            write(ln,10)'ENDIF'
            write(ln,10)'ENDIF'
         endif
      endif
      enddo

      write(ln,10) 'RETURN'
      write(ln,10) 'END'
      write(ln,5) ' '
      write(ln,5) ' '
      write(ln,5) ' '
      write(ln,5) ' '
      write(ln,5) ' '

      write(ln,10)'SUBROUTINE DIPOLEPDF(x1,x2,leg1,leg2,WGT)'
      write(ln,10)'IMPLICIT NONE'
      write(ln,5) 'This is a template for a routine to calculate the pdfs of the'
      write(ln,5) 'single dipoles. The output variable is assumed to'
      write(ln,5) 'be wgt. The parton fractions are given by x1 and x2.'
      write(ln,5)
     &     'leg1 and leg2 are the PDG codes for two incoming partons.'
      write(ln,5)'The factorization scale is included via a common block in dipole.inc'
      write(ln,10) ' '
      write(ln,10)'include "nexternal.inc"'
      write(ln,10)'include "dipole.inc"'
      write(ln,5)'include "run.inc"'
      write(ln,10)'double precision X1,X2,WGT'
      write(ln,10)'INTEGER leg1,leg2'
      write(ln,10)' '
      write(ln,10)'IF (x1.gt.1d0 .or. x2.gt.1d0) THEN'
      write(ln,10)'   WGT=0d0'
      write(ln,10)'   RETURN'
      write(ln,10)'ENDIF'
      write(ln,5)'This is an empty routine, so return 1'
      write(ln,10) 'WGT=1d0'
      write(ln,10) ' '
      write(ln,10) 'END'



      write(ln,10) ' '
      write(ln,10) ' '
      write(ln,10) ' '


      write(ln,10)'SUBROUTINE WRITEHIST(P,WGT)'
      write(ln,10)'IMPLICIT NONE'
      write(ln,5) 'This is a template that can be used for writing events or histograms.'
      write(ln,5) 'The input variables are the phase space point and the according'
      write(ln,5) 'weight.'
      write(ln,10)' '
      write(ln,10)'include "nexternal.inc"'
      write(ln,10)'DOUBLE PRECISION P(0:3,NEXTERNAL),WGT'
      write(ln,10)' '
      write(ln,10)' '
      write(ln,10)'END'

      write(ln,10) ' '
      write(ln,10) ' '
      write(ln,10) ' '

      write(ln,10)'LOGICAL FUNCTION PASSCUTSDIP(P)'
      write(ln,10)'IMPLICIT NONE'
      write(ln,5) 'This is a template for a routine that puts cuts on each'
      write(ln,5) 'individual dipole, e.g. a jet algorithm. The input variable'
      write(ln,5) 'is the phase space point P. If event passes the cut, this'
      write(ln,5) 'function has to be set TRUE.'
      write(ln,10) ' '
      write(ln,10)'include "nexternal.inc"'
      write(ln,10)'DOUBLE PRECISION P(0:3,NEXTERNAL)'
      write(ln,10)' '
      write(ln,10)' '
      write(ln,5) 'Set to TRUE by default.'
      write(ln,10)'PASSCUTSDIP=.TRUE.'
      write(ln,10)' '
      write(ln,10)'END'      
      




 5    format('C ',A)
 10   format(6X,A)
 11   format(6X,A,I2,A)
 12   format(5X,'&',A,I2,A)
 13   format(6X,A,I2,A)
 14   format(6X,A,5(I2,A),A8,A,A8,A,A8,A)
 15   format('C',6X,A)
 17   format(6X,A,I3,A,I3,A)
 20   format(', P',I1)
 30   format(',P',I1,'(0:3)')
 40   format(6X,A,I4,A)
 41   format(6X,A23,I4,A2)
 50   format(7X,A26,I4,A1)
 51   format(7X,A,I4,A)
 52   format('C ',7X,A,I4,A)
 60   format(6X,A11,I4,A1)
 70   format(6X,A12,I2,A6,I2,A2)
 71   format(6X,A,A5)
 72   format(6X,A15,2(I1,A7),I1,A)
 73   format(6X,A,I2,A,I2,A,I2,A)
      close(ln)
      end


      Subroutine write_intdipolesqed()


      implicit none

      include 'params.inc'

c     Global
      integer dipollist(maxdipoles,2,2),nrdipols
      common /list/ dipollist,nrdipols

      integer info_ijk(maxdipoles,2,5),nr
      common /ijk/ info_ijk,nr

      character*70       cversion
      common/to_version/ cversion

      character*60 proc
      integer iproc
      common/to_write/iproc,proc

      integer         nincoming,nincfirst
      common/to_proc2/nincoming,nincfirst

      character*(8) mass(max_particles)
      common /masses/ mass

      integer info_p(5,max_particles),iposx(3,max_particles)
      character*(max_string) iwave(max_particles),owave(max_particles)
      character*(8) str(3,max_particles)
      common/to_external/iwave,owave,iposx,info_p,str

      integer iline_orig(-maxlines:maxlines)
      common /ciline_orig/iline_orig

      integer iline_copy(-maxlines:maxlines,maxdipoles)
      common /to_writedipole/iline_copy

      integer        iline(-maxlines:maxlines),idir(-maxlines:maxlines)
      integer        this_coup(max_coup) ,goal_coup(max_coup)
      common/to_proc/iline,idir,this_coup,goal_coup

      character*(max_string) dipprocess(maxdipoles)
      common /process/ dipprocess

      character*(4) part_type(maxdipoles,2)
      common /part/ part_type

      real*8 chargeinfo(max_particles)
      integer chargint(maxdipoles,2)
      common /charge_info/ chargeinfo,chargint

      character*(4*max_particles) particle(4)
      integer                               charge_c(max_particles)
      integer iparticle(0:max_particles,0:4),inverse(max_particles)
      common/to_model/iparticle,  particle,  inverse, charge_c

      logical photonpdf
      common/photon/ photonpdf


c     Local
      integer ln,i,j,ascii_emitter,ascii_spectator,ascii_nr
      integer nrdipol,leg, part(max_particles)
      parameter (ln=82)
      logical first(max_particles,max_particles)
      logical lepton,photon


C  Check if there are leptons present in the initial state
      lepton=.false.
      photon=.false.
      do i=1, nr
        if( ( abs(info_p(5,inverse(iline_copy(1,i)))).gt.10.and.
     &     info_p(5,inverse(iline_copy(1,i))).ne.21.and.
     &     info_p(5,inverse(iline_copy(1,i))).ne.22) .or.
     &    (abs(info_p(5,inverse(iline_copy(2,i)))).gt.10.and.
     &     info_p(5,inverse(iline_copy(2,i))).ne.21.and.
     &     info_p(5,inverse(iline_copy(2,i))).ne.22) ) then
         lepton=.true.
        endif
        if( ( abs(info_p(5,inverse(iline_copy(1,i)))).eq.22) .or.
     &    ( abs(info_p(5,inverse(iline_copy(2,i)))).eq.22)) then
         photon=.true.
        endif
      enddo


      do i=1,nr
       do j=1,nr
        first(i,j)=.true.
       enddo
      enddo



      nrdipol=nr
      if(nrdipol.eq.0) nrdipol=1

      open(unit=ln,file='intdipolesqed.f',status='unknown')

      write(ln,10) 'SUBROUTINE INTDIPOLESQED(P,P1,P2,X,Z,ZI,PSWGT,EPS,EPSSQ,FIN)'
      write(ln,5) ' '
      write(ln,5) cversion
      write(ln,5) 
     .     'RETURNS COEFFICIENTS OF 1/eps^2,1/eps AND FINITE TERMS'
      write(ln,5) 'FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL) '
      write(ln,5) ' '
      write(ln,5) 'FOR PROCESS '//proc(1:iproc)
      write(ln,5) 'WHERE Z IS THE MOMENTUM FRACTION '//
     &     'OF THE INITAL EMITTER/SPECTATOR.'
      write(ln,5) ' '

      write(ln,10) 'IMPLICIT NONE'
      write(ln,10) 'include "nexternal.inc" '
      write(ln,10) 'include "dipole.inc" '
      write(ln,10) 'include "coupl.inc"'

      write(ln,5) ' '
      write(ln,5) 'ARGUMENTS '
      write(ln,5) ' '
      write(ln,10) 'REAL*8 P(0:3,NEXTERNAL),P1(0:3,NEXTERNAL),P2(0:3,NEXTERNAL),X(2),Z,PSWGT(3)'
      write(ln,10) 'REAL*8 ZI,EPS,EPSSQ,FIN'
      write(ln,5) ' '
      write(ln,5) 'GLOBAL '
      write(ln,5) ' '
      write(ln,10) 'REAL*8 dot,lambda_tr,DFRAG'
      write(ln,10) 'LOGICAL PASSCUTSDIP,THETACUT'
      write(ln,10) 'EXTERNAL dot,lambda_tr,PASSCUTSDIP,THETACUT,DFRAG'
      write(ln,5) ' '
      write(ln,5) 'LOCAL'
      write(ln,5) ' '
c      write(ln,10)
c     &     'CHARACTER*(4) u,ubar,d,dbar,c,cbar,s,sbar,b,bar,t,tbar,glu,nqcd'
      write(ln,10) 'REAL*8 sij,sik,sjk,zp,e,esq,pdfwgt(3),FINITE(3)'
      write(ln,10) 'REAL*8 y_ijk,omx_ija,u_i,vtilde_i,zplus'
      write(ln,10) 'REAL*8 mass_i,mass_j,mass_k,mass_ij'
      write(ln,10) 'REAL*8 msq_i,msq_j,msq_k,msq_ij,sikzone'
      write(ln,10) 'REAL*8 musq_i,musq_j,musq_k,musq_ij'
      write(ln,10) 'REAL*8 vij,vtildeij,rho_i,rho_k,rho,qsq'
      write(ln,10) 'REAL*8 one,two,three,four,qtot(0:3),pi'
      write(ln,10) 'REAL*8 gsq,gammaE,ff,fi(9),sab,sabt'
      write(ln,10) 'parameter (one=1d0,two=2d0,three=3d0,four=4d0)'
c      write(ln,10) 'parameter (cf=4d0/3d0,ca=3d0)'
      write(ln,60) 'REAL*8 ANS(',nrdipol,')'
      if(lepton.or.(photon.and.(photonpdf.eqv..false.))) then
       write(ln,60) 'REAL*8 ANSZONE(',nrdipol,')'
      endif
      write(ln,10) 'REAL*8 ZERO'
      write(ln,10) 'INTEGER NU'
      write(ln,10) 'PARAMETER (ZERO=0D0)'
      write(ln,10) 'PARAMETER (gammaE=0.57721566d0)'
      write(ln,10) 'PARAMETER (pi=3.1415926d0)'
      write(ln,5) ' '
      write(ln,5) 'BEGIN CODE'
      write(ln,5) ' '
      write(ln,10) 'INTDIPOLE=.TRUE.'
c      write(ln,10) 'gsq=GG(1)**2'
      write(ln,10) 'eps=0d0'
      write(ln,10) 'epssq=0d0'
      write(ln,10) 'fin=0d0'
      write(ln,10) 'e=0d0'
      write(ln,10) 'esq=0d0'
      write(ln,10) 'DO NU=1,3'
      write(ln,10) '   finite(NU)=0d0'
      write(ln,10) 'ENDDO'
      write(ln,60) 'DO NU=1,',nrdipol
      write(ln,10) '   ANS(NU)=0d0'
      if(lepton.or.(photon.and.(photonpdf.eqv..false.))) then
       write(ln,10) '   ANSZONE(NU)=0d0'
      endif
      write(ln,10) 'ENDDO'
      if(.not.lepton.and..not.photon.or.(photon.and.photonpdf)) then
       write(ln,10) 'IF (PASSCUTSDIP(P)) THEN'
      endif 
      do i=1,nr



c     First check if unresolved is massless. If not do not print it here.
         if (mass(info_ijk(i,1,2)).eq.'ZERO')then
c**************************************************************
c       final-final theta funtion
c**************************************************************
         if(info_ijk(i,1,1).gt.nincoming
     &        .and. info_ijk(i,1,3).gt.nincoming) then
c  check which particle is a lepton and which hadron
         if((abs(info_p(5,inverse(iline_orig(info_ijk(i,1,1))))).gt.10.and.
     &      info_p(5,inverse(iline_orig(info_ijk(i,1,1)))).ne.21) .or.
     &      (abs(info_p(5,inverse(iline_orig(info_ijk(i,1,2))))).ne.22) ) then
           part(info_ijk(i,1,1))=0
         elseif((abs(info_p(5,inverse(iline_orig(info_ijk(i,1,1))))).le.10.or.
     &      abs(info_p(5,inverse(iline_orig(info_ijk(i,1,1))))).eq.21).and.
     &      abs(info_p(5,inverse(iline_orig(info_ijk(i,1,2))))).eq.22) then
           part(info_ijk(i,1,1))=1
        endif

            write(ln,10)' '
            write(ln,71)'mass_i=',mass(info_ijk(i,1,1))
            write(ln,71)'mass_k=',mass(info_ijk(i,1,3))
            write(ln,10) 'DO NU=0,3'
            write(ln,73) ' qtot(NU)=p(NU,',info_ijk(i,1,4),')+p(NU,',info_ijk(i,1,5),')'
c            write(ln,10) ' IF (NCS) THEN'
c            write(ln,76) ' qtot(NU)=p(NU,',info_ijk(i,1,1),')+p(NU,',info_ijk(i,1,2),
c     &                    ')+p(NU,',info_ijk(i,1,3),')'
c            write(ln,10) ' ENDIF'
            write(ln,10) 'ENDDO'
            write(ln,10)'qsq=dot(qtot,qtot)'
            write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &           info_p(5,inverse(iline_copy(1,i))),',',
     &           info_p(5,inverse(iline_copy(2,i))),',pdfwgt(2))'
           if(info_ijk(i,2,1).eq.info_ijk(i,2,2).and..not.
     &          first(info_ijk(i,1,1),info_ijk(i,1,2))) then
             write(ln,15)'IF (PASSCUTSDIP(P)) THEN'
             if (i.le.9) then
                write(ln,61)' CALL SDIPOLQED'//char(48)//char(48)
     &            //char(48+i)//'(P,ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,61)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &            //char(mod(i,10)+48)//'(P,ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,61)' CALL SDIPOLQED'//char(i/100+48)
     &            //char((i-(i/100)*100)/10+48)
     &            //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             write(ln,15)'ENDIF'
           else
             write(ln,10)'IF (PASSCUTSDIP(P)) THEN'
             if (i.le.9) then
                write(ln,62)
     &           ' CALL SDIPOLQED'//char(48)//char(48)//char(48+i)
     &               //'(P,ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,62)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &               //char(mod(i,10)+48)//'(P,ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,62)' CALL SDIPOLQED'//char(i/100+48)
     &               //char((i-(i/100)*100)/10+48)
     &               //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             first(info_ijk(i,1,1),info_ijk(i,1,2))=.false.
             write(ln,10) 'ENDIF'
           endif

            write(ln,19) ' CALL EPSFFQED(mass_i,mass_k,qsq,zi,',info_ijk(i,2,4),
     &           ',',info_ijk(i,2,1),',',chargint(i,1),',',chargint(i,2),',',
     &            part(info_ijk(i,1,1)),',esq,e)'
            write(ln,11) 'epssq=epssq-esq*ANS(',i,')*pdfwgt(2)'
            write(ln,11) 'eps=eps-e*ANS(',i,')*pdfwgt(2)'
            if(info_ijk(i,2,1).eq.0) then
              leg=info_ijk(i,1,2)
            elseif(info_ijk(i,2,2).eq.0) then
              leg=info_ijk(i,1,1)
            elseif(info_ijk(i,2,1).eq.info_ijk(i,2,2)) then
              leg=info_ijk(i,1,1)
            endif
            write(ln,77) ' CALL FINITEFFQED(mass_i,mass_k,qsq,zi,',info_ijk(i,2,4),
     &           ',',info_ijk(i,2,1),',',chargint(i,1),',',chargint(i,2),',',leg,',',
     &            part(info_ijk(i,1,1)),',fi)'
            if(lepton.or.(photon.and..not.photonpdf)) then
             write(ln,11) ' finite(2)=finite(2)-fi(2)*ANS(',i,')*pdfwgt(2)*pswgt(1)'
             write(ln,11 ) ' CALL WRITEHIST(P,-fi(2)*ANS(',i,')*pdfwgt(2)*pswgt(1))'
            if(abs(info_p(5,inverse(iline_orig(info_ijk(i,1,2))))).eq.22) then
            write(ln,10) 'IF(NCS) THEN'
            if(part(info_ijk(i,1,1)).eq.1) then
            write(ln,10) ' IF(PHOTONFRAG) THEN'
            write(ln,10) '  fi(3)=fi(3)+DFRAG(ZI)'
            write(ln,10) ' ENDIF'
            endif
            write(ln,10) ' IF(THETACUT(P,zi)) THEN'
            write(ln,11) '  finite(2)=finite(2)-fi(3)*ANS(',i,')*pdfwgt(2)*pswgt(1)'
            write(ln,11) '  CALL WRITEHIST(P,-fi(3)*ANS(',i,')*pdfwgt(2)*pswgt(1))'
            write(ln,10) ' ENDIF'
            write(ln,10) ' IF(THETACUT(P,1d0)) THEN'
            write(ln,11) '  finite(2)=finite(2)+fi(3)*ANS(',i,')*pdfwgt(2)*pswgt(1)'
            write(ln,11) '  CALL WRITEHIST(P,+fi(3)*ANS(',i,')*pdfwgt(2)*pswgt(1))'
            write(ln,10) ' ENDIF'            
            write(ln,10) 'ENDIF'
            endif
            else
             write(ln,11) ' finite(2)=finite(2)-fi(2)*ANS(',i,')*pdfwgt(2)'
            if(abs(info_p(5,inverse(iline_orig(info_ijk(i,1,2))))).eq.22) then
            write(ln,10) 'IF(NCS) THEN'
            if(part(info_ijk(i,1,1)).eq.1) then
            write(ln,10) ' IF(PHOTONFRAG) THEN'
            write(ln,10) '  fi(3)=fi(3)+DFRAG(ZI)'
            write(ln,10) ' ENDIF'
            endif
            write(ln,10) ' IF(THETACUT(P,zi)) THEN'
            write(ln,11) '  finite(2)=finite(2)-fi(3)*ANS(',i,')*pdfwgt(2)'
            write(ln,10) ' ENDIF'
            write(ln,10) ' IF(THETACUT(P,1d0)) THEN'
            write(ln,11) '  finite(2)=finite(2)+fi(3)*ANS(',i,')*pdfwgt(2)'
            write(ln,10) ' ENDIF'            
            write(ln,10) 'ENDIF'
            endif
            endif



c**************************************************************
c       final-initial theta function
c**************************************************************
         elseif(info_ijk(i,1,1).gt.nincoming
     &           .and.info_ijk(i,1,3).le.nincoming) then

c  check which particle is a lepton and which hadron
         if((abs(info_p(5,inverse(iline_orig(info_ijk(i,1,1))))).gt.10.and.
     &      info_p(5,inverse(iline_orig(info_ijk(i,1,1)))).ne.21) .or.
     &      (abs(info_p(5,inverse(iline_orig(info_ijk(i,1,2))))).ne.22) ) then
           part(info_ijk(i,1,1))=0
         elseif((abs(info_p(5,inverse(iline_orig(info_ijk(i,1,1))))).le.10.or.
     &      abs(info_p(5,inverse(iline_orig(info_ijk(i,1,1))))).eq.21).and.
     &      abs(info_p(5,inverse(iline_orig(info_ijk(i,1,2))))).eq.22) then
           part(info_ijk(i,1,1))=1
        endif
c        print*, 'i',i,'info_ijk(i,1,1)',info_ijk(i,1,1),abs(info_p(5,inverse(iline_orig(info_ijk(i,1,1)))))
c            if(info_ijk(i,2,1).eq.0) then
c              leg=info_ijk(i,1,2)
c            elseif(info_ijk(i,2,2).eq.0) then
c              leg=info_ijk(i,1,1)
c            elseif(info_ijk(i,2,1).eq.info_ijk(i,2,2)) then
c              leg=info_ijk(i,1,1)
c            endif
            leg=info_ijk(i,1,3)
            write(ln,10)' '
            write(ln,71)'mass_i=',mass(info_ijk(i,1,1))
            if(info_ijk(i,2,4).eq.0 .and. (
     &           mass(info_ijk(i,1,1)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,2)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,3)).ne.'ZERO')) then
               write (ln,11) 'zplus=1d0-4d0*mass_i**2/sik'
               write(ln,74)'sikzone=two*dot(p(0,',info_ijk(i,1,4),
     &              '),p(0,',info_ijk(i,1,5),
     &              '))/zplus'
            else
               write(ln,74)'sikzone=two*dot(p(0,',info_ijk(i,1,4),
     &              '),p(0,',info_ijk(i,1,5),'))'
            endif

           if(.not.lepton.and..not.photon.or.(photon.and.photonpdf)) then

            write(ln,74)'sik=two*dot(p(0,',info_ijk(i,1,4),'),p(0,'
     &           ,info_ijk(i,1,5),'))/z'
            if(info_ijk(i,1,3).eq.1) then
               write(ln,17)'CALL DIPOLEPDF(x(1)/z,x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(1))'
            elseif(info_ijk(i,1,3).eq.2) then
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2)/z,',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(1))'
            else
               write (ln,*)'Error #145 in write_intdipole'
               write (*,*)'Error #145 in write_intdipole'
            endif
            if(info_ijk(i,2,4).eq.0 .and. (
     &           mass(info_ijk(i,1,1)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,2)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,3)).ne.'ZERO')) then
               if(info_ijk(i,1,3).eq.1) then
                  write(ln,17)'CALL DIPOLEPDF(x(1)/zplus,x(2),',
     &                 info_p(5,inverse(iline_copy(1,i))),',',
     &                 info_p(5,inverse(iline_copy(2,i))),',pdfwgt(3))'
               elseif(info_ijk(i,1,3).eq.2) then
                  write(ln,17)'CALL DIPOLEPDF(x(1),x(2)/zplus,',
     &                 info_p(5,inverse(iline_copy(1,i))),',',
     &                 info_p(5,inverse(iline_copy(2,i))),',pdfwgt(3))'
               endif
            else
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(2))'
            endif
           if(info_ijk(i,2,1).eq.info_ijk(i,2,2).and..not.
     &          first(info_ijk(i,1,1),info_ijk(i,1,2))) then
             if (i.le.9) then
                write(ln,61)'CALL SDIPOLQED'//char(48)//char(48)
     &            //char(48+i)//'(P,ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,61)'CALL SDIPOLQED'//char(48)//char(i/10+48)
     &            //char(mod(i,10)+48)//'(P,ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,61)'CALL SDIPOLQED'//char(i/100+48)
     &            //char((i-(i/100)*100)/10+48)
     &            //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
           else
             if (i.le.9) then
                write(ln,62)'CALL SDIPOLQED'//char(48)//char(48)//char(48+i)
     &               //'(P,ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,62)'CALL SDIPOLQED'//char(48)//char(i/10+48)
     &               //char(mod(i,10)+48)//'(P,ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,62)'CALL SDIPOLQED'//char(i/100+48)
     &               //char((i-(i/100)*100)/10+48)
     &               //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             first(info_ijk(i,1,1),info_ijk(i,1,2))=.false.
           endif


             write(ln,19)' CALL EPSFIQED(',info_ijk(i,2,4),',',info_ijk(i,2,1),
     &         ',',chargint(i,1),',',chargint(i,2),',sikzone,mass_i,zi,',
     &          part(info_ijk(i,1,1)),',esq,e)'
            write(ln,11) 'epssq=epssq-esq*ANS(',i,')*pdfwgt(2)'
            write(ln,11) 'eps=eps-e*ANS(',i,')*pdfwgt(2)'
             write(ln,10) ' sab=massint(1)**2+massint(2)**2+two*dot(p(0,1),p(0,2))'
             write(ln,10) ' sabt=massint(1)**2+massint(2)**2+two*dot(p(0,1),p(0,2))*z'
c            if (info_ijk(i,1,3).eq.1) then
             write(ln,77) ' CALL FINITEFIQED(mass_i,sik,sikzone,z,',info_ijk(i,2,4),
     &         ',',info_ijk(i,2,1),',',chargint(i,1),',',chargint(i,2),
     &         ',',leg,',sab,sabt,',part(info_ijk(i,1,1)),',fi)'
c            elseif(info_ijk(i,1,3).eq.2) then
c             write(ln,19) ' CALL FINITEFIQED(mass_i,sik,sikzone,z,',info_ijk(i,2,4),
c     &         ',',info_ijk(i,2,1),',',chargint(i,1),',',chargint(i,2),
c     &         ',',info_ijk(i,1,3),',sab,sabt,fi)'
c            endif
            if(info_ijk(i,2,4).eq.0 .and. (
     &           mass(info_ijk(i,1,1)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,2)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,3)).ne.'ZERO')) then
               write(ln,11)
     &            'finite(2)=finite(2)-(-fi(8))*ANS(',i,')*pdfwgt(2)'
               write(ln,11) 'finite(3)=finite(3)-(fi(7)+fi(9))'//
     &              '*ANS(',i,')*pdfwgt(3)'
            else
               write(ln,11) 'finite(1)=finite(1)-(fi(1)+fi(3))'//
     &              '*ANS(',i,')*pdfwgt(1)'
               write(ln,11) 'finite(2)=finite(2)-(fi(2)-fi(4))'//
     &              '*ANS(',i,')*pdfwgt(2)'
             if(abs(info_p(5,inverse(iline_orig(info_ijk(i,1,2))))).eq.22) then
               write(ln,10) 'IF(NCS) THEN'
               if(part(info_ijk(i,1,1)).eq.1) then
               write(ln,10) ' IF(PHOTONFRAG) THEN'
               write(ln,10) '  fi(5)=fi(5)+DFRAG(ZI)'
               write(ln,10) '  fi(6)=fi(5)'
               write(ln,10) ' ENDIF'
               endif
               write(ln,10) ' IF(THETACUT(P,zi)) THEN'
               write(ln,11) '  finite(2)=finite(2)-fi(5)*ANS(',i,')*pdfwgt(2)'
               write(ln,11) '  finite(1)=finite(1)-fi(7)*ANS(',i,')*pdfwgt(1)'
               write(ln,10) ' ENDIF'
               write(ln,10) ' IF(THETACUT(P,1d0)) THEN'
               write(ln,11) '  finite(2)=finite(2)+fi(6)*ANS(',i,')*pdfwgt(2)'
               write(ln,11) '  finite(2)=finite(2)+fi(8)*ANS(',i,')*pdfwgt(1)'
               write(ln,10) ' ENDIF'
               write(ln,10) 'ENDIF'
              endif
            endif


           elseif(lepton.or.(photon.and..not.photonpdf)) then

            write(ln,75)'sik=two*dot(p',info_ijk(i,1,5),'(0,',info_ijk(i,1,4),
     &              '),p',info_ijk(i,1,5),'(0,',info_ijk(i,1,5),'))'
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(2))'
           if(info_ijk(i,2,1).eq.info_ijk(i,2,2).and..not.
     &          first(info_ijk(i,1,1),info_ijk(i,1,2))) then
             write(ln,15)'IF (PASSCUTSDIP(P)) THEN'
             if (i.le.9) then
                write(ln,61)' CALL SDIPOLQED'//char(48)//char(48)
     &            //char(48+i)//'(P,ANSZONE(',i,'))'
             elseif(i.le.99) then
                write(ln,61)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &            //char(mod(i,10)+48)//'(P,ANSZONE(',i,'))'
             elseif(i.le.999) then
                write(ln,61)' CALL SDIPOLQED'//char(i/100+48)
     &            //char((i-(i/100)*100)/10+48)
     &            //char(mod(i-(i/100)*100,10)+48)//'(P,ANSZONE(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             write(ln,15)'ENDIF'
c             write(ln,6)' CALL BOOST(P,Z,',info_ijk(i,1,3),',P1)'
             write(ln,9)'IF (PASSCUTSDIP(P',info_ijk(i,1,5),')) THEN'
             if (i.le.9) then
                write(ln,51)' CALL SDIPOLQED'//char(48)//char(48)
     &            //char(48+i)//'(P',info_ijk(i,1,5),',ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,51)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &            //char(mod(i,10)+48)//'(P',info_ijk(i,1,5),',ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,51)' CALL SDIPOLQED'//char(i/100+48)
     &            //char((i-(i/100)*100)/10+48)
     &            //char(mod(i-(i/100)*100,10)+48)//'(P',info_ijk(i,1,5),',ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             write(ln,15)'ENDIF'
           else
             write(ln,10)'IF (PASSCUTSDIP(P)) THEN'
             if (i.le.9) then
                write(ln,62)' CALL SDIPOLQED'//char(48)//char(48)//char(48+i)
     &               //'(P,ANSZONE(',i,'))'
             elseif(i.le.99) then
                write(ln,62)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &               //char(mod(i,10)+48)//'(P,ANSZONE(',i,'))'
             elseif(i.le.999) then
                write(ln,62)' CALL SDIPOLQED'//char(i/100+48)
     &               //char((i-(i/100)*100)/10+48)
     &               //char(mod(i-(i/100)*100,10)+48)//'(P,ANSZONE(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             write(ln,10) 'ENDIF'
c             write(ln,13)' CALL BOOST(P,Z,',info_ijk(i,1,3),',P1)'
             write(ln,7)'IF (PASSCUTSDIP(P',info_ijk(i,1,5),')) THEN'
             if (i.le.9) then
                write(ln,41)' CALL SDIPOLQED'//char(48)//char(48)//char(48+i)
     &               //'(P',info_ijk(i,1,5),',ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,41)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &               //char(mod(i,10)+48)//'(P',info_ijk(i,1,5),',ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,41)' CALL SDIPOLQED'//char(i/100+48)
     &               //char((i-(i/100)*100)/10+48)
     &               //char(mod(i-(i/100)*100,10)+48)//'(P',info_ijk(i,1,5),',ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             write(ln,10) 'ENDIF' 
             first(info_ijk(i,1,1),info_ijk(i,1,2))=.false.
           endif

             write(ln,19)' CALL EPSFIQED(',info_ijk(i,2,4),',',info_ijk(i,2,1),
     &         ',',chargint(i,1),',',chargint(i,2),',sikzone,mass_i,zi,',
     &          part(info_ijk(i,1,1)),',esq,e)'
            write(ln,11) 'epssq=epssq-esq*ANSZONE(',i,')*pdfwgt(2)'
            write(ln,11) 'eps=eps-e*ANSZONE(',i,')*pdfwgt(2)'
             write(ln,10) ' sab=massint(1)**2+massint(2)**2+two*dot(p(0,1),p(0,2))'
             write(ln,74) ' sabt=massint(1)**2+massint(2)**2+two*dot(p'
     &                      ,info_ijk(i,1,5),'(0,1),p',info_ijk(i,1,5),'(0,2))'
c            if (info_ijk(i,1,1).eq.1) then
c             write(ln,19) ' CALL FINITEFIQED(mass_i,sik,sikzone,z,',info_ijk(i,2,4),
c     &         ',',info_ijk(i,2,1),',',chargint(i,1),',',chargint(i,2),
c     &         ',',info_ijk(i,1,1),',sab,sabt,fi)'
c            elseif(info_ijk(i,1,1).eq.2) then
             write(ln,77) ' CALL FINITEFIQED(mass_i,sik,sikzone,z,',info_ijk(i,2,4),
     &         ',',info_ijk(i,2,1),',',chargint(i,1),',',chargint(i,2),
     &         ',',leg,',sab,sabt,',part(info_ijk(i,1,1)),',fi)'
c            endif
            if(info_ijk(i,2,4).eq.0 .and. (
     &           mass(info_ijk(i,1,1)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,2)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,3)).ne.'ZERO')) then
               write(ln,11)
     &            'finite(2)=finite(2)-(-fi(8))*ANSZONE(',i,')*pdfwgt(2)*pswgt(1)'
               write(ln,80) 'finite(3)=finite(3)-(fi(7)+fi(9))'//
     &              '*ANS(',i,')*pdfwgt(2)*pswgt(',info_ijk(i,1,5)+1,')'
               write(ln,11) ' CALL WRITEHIST(P,fi(8)*ANSZONE(',i,')*pdfwgt(2)*pswgt(1))'
               write(ln,8) ' CALL WRITEHIST(P',info_ijk(i,1,5),',-(fi(7)+fi(9))*ANS(',i,
     &                      ')*pdfwgt(2)*pswgt(',info_ijk(i,1,5)+1,'))'
            else
               write(ln,80) 'finite(1)=finite(1)-(fi(1)+fi(3))'//
     &              '*ANS(',i,')*pdfwgt(2)*pswgt(',info_ijk(i,1,5)+1,')'
               write(ln,11) 'finite(2)=finite(2)-(fi(2)-fi(4))'//
     &              '*ANSZONE(',i,')*pdfwgt(2)*pswgt(1)'
               write(ln,11) ' CALL WRITEHIST(P,-(fi(2)-fi(4))*ANSZONE(',i,')*pdfwgt(2)*pswgt(1))'
               write(ln,8) ' CALL WRITEHIST(P',info_ijk(i,1,5),',-(fi(1)+fi(3))/z*ANS(',i,
     &                      ')*pdfwgt(3)*pswgt(',info_ijk(i,1,5)+1,'))'
             if(abs(info_p(5,inverse(iline_orig(info_ijk(i,1,2))))).eq.22) then
               write(ln,10) 'IF(NCS) THEN'
               if(part(info_ijk(i,1,1)).eq.1) then
               write(ln,10) ' IF(PHOTONFRAG) THEN'
               write(ln,10) '  fi(5)=fi(5)+DFRAG(ZI)'
               write(ln,10) '  fi(6)=fi(5)'
               write(ln,10) ' ENDIF'
               endif
               write(ln,10) ' IF(THETACUT(P,zi)) THEN'
               write(ln,11) '  finite(2)=finite(2)-fi(5)*ANSZONE(',i,')*pdfwgt(2)*pswgt(1)'
               write(ln,80)  '  finite(1)=finite(1)-fi(7)*ANS(',i,')*pdfwgt(2)*pswgt(',info_ijk(i,1,5)+1,')'
               write(ln,11) '  CALL WRITEHIST(P,-fi(5)*ANSZONE(',i,')*pdfwgt(2)*pswgt(1))'
               write(ln,8)  '  CALL WRITEHIST(P',info_ijk(i,1,5),',-fi(7)/z*ANS(',i,
     &                      ')*pdfwgt(3)*pswgt(',info_ijk(i,1,5)+1,'))'
               write(ln,10) ' ENDIF'
               write(ln,10) ' IF(THETACUT(P,1d0)) THEN'
               write(ln,11) '  finite(2)=finite(2)+fi(6)*ANSZONE(',i,')*pdfwgt(2)*pswgt(1)'
               write(ln,11) '  finite(2)=finite(2)+fi(8)*ANSZONE(',i,')*pdfwgt(2)*pswgt(1)'
               write(ln,11) '  CALL WRITEHIST(P,(fi(6)+fi(8))*ANSZONE(',i,')*pdfwgt(2)*pswgt(1))'
               write(ln,10) ' ENDIF'
               write(ln,10) 'ENDIF'
              endif
            endif

          

           endif

c**************************************************************
c       initial-final theta function
c**************************************************************
         elseif(info_ijk(i,1,1).le.nincoming .and.
     &           info_ijk(i,1,3).gt.nincoming) then
            if(info_ijk(i,2,1).eq.0) then
              leg=info_ijk(i,1,2)
            elseif(info_ijk(i,2,2).eq.0) then
              leg=info_ijk(i,1,1)
            elseif(info_ijk(i,2,1).eq.info_ijk(i,2,2)) then
              leg=info_ijk(i,1,1)
            endif
            write(ln,10)' '
            write(ln,71)'mass_k=',mass(info_ijk(i,1,3))
            write(ln,74)'sikzone=two*dot(p(0,',info_ijk(i,1,4)
     &           ,'),p(0,',info_ijk(i,1,5),'))'


          if(.not.lepton.and..not.photon.or.(photon.and.photonpdf)) then
  
            write(ln,74)'sik=two*dot(p(0,',info_ijk(i,1,4)
     &           ,'),p(0,',info_ijk(i,1,5),'))/z'
            if (info_ijk(i,1,1).eq.1) then
               write(ln,17)'CALL DIPOLEPDF(x(1)/z,x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(1))'
            elseif(info_ijk(i,1,1).eq.2) then
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2)/z,',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(1))'
            else
               write (ln,*)'Error #146 in write_intdipole'
               write (*,*)'Error #146 in write_intdipole'
            endif
            if(((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.1).or.
     &           (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.0))) then
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(2))'
            endif
          if(info_ijk(i,2,1).eq.info_ijk(i,2,2).and..not.
     &          first(info_ijk(i,1,1),info_ijk(i,1,2))) then
             if (i.le.9) then
                write(ln,61)'CALL SDIPOLQED'//char(48)//char(48)
     &            //char(48+i)//'(P,ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,61)'CALL SDIPOLQED'//char(48)//char(i/10+48)
     &            //char(mod(i,10)+48)//'(P,ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,61)'CALL SDIPOLQED'//char(i/100+48)
     &            //char((i-(i/100)*100)/10+48)
     &            //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
           else
             if (i.le.9) then
                write(ln,62)'CALL SDIPOLQED'//char(48)//char(48)
     &               //char(48+i)//'(P,ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,62)'CALL SDIPOLQED'//char(48)//char(i/10+48)
     &               //char(mod(i,10)+48)//'(P,ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,62)'CALL SDIPOLQED'//char(i/100+48)
     &               //char((i-(i/100)*100)/10+48)
     &               //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             first(info_ijk(i,1,1),info_ijk(i,1,2))=.false.
           endif

              write(ln,19) ' CALL EPSIFQED(',info_ijk(i,2,1),',',info_ijk(i,2,4),
     &          ',',chargint(i,1),',',chargint(i,2),',',info_ijk(i,1,1),',sikzone,mass_k,esq,e)'
               write(ln,11) 'epssq=epssq-esq*ANS(',i,')*pdfwgt(2)'
               write(ln,11) 'eps=eps-e*ANS(',i,')*pdfwgt(2)'
             write(ln,10) ' sab=massint(1)**2+massint(2)**2+two*dot(p(0,1),p(0,2))/z'
             write(ln,10) ' sabt=massint(1)**2+massint(2)**2+two*dot(p(0,1),p(0,2))'
            if (info_ijk(i,1,1).eq.1) then
             write(ln,77) ' CALL FINITEIFQED(mass_k,sik,sikzone,z,',info_ijk(i,2,1),
     &          ',',info_ijk(i,2,4),',',chargint(i,1),',',chargint(i,2),','
     &          ,leg,',',info_p(5,inverse(iline_copy(1,i))),',sab,sabt,fi)'
            elseif(info_ijk(i,1,1).eq.2) then
             write(ln,77) ' CALL FINITEIFQED(mass_k,sik,sikzone,z,',info_ijk(i,2,1),
     &          ',',info_ijk(i,2,4),',',chargint(i,1),',',chargint(i,2),','
     &          ,leg,',',info_p(5,inverse(iline_copy(2,i))),',sab,sabt,fi)'
            endif
c
c quark-quark splitting (i:quark, ij~: quark) and
c gluon-gluon splitting (i:gluon, ij~: gluon)
c
            if(((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.1).or.
     &           (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.0))) then
               if (mass(info_ijk(i,1,3)).eq.'ZERO') then
                  write(ln,11) 'finite(1)=finite(1)-(fi(1)+fi(3))'//
     &                 '*ANS(',i,')*pdfwgt(1)'
                  write(ln,11) 'finite(2)=finite(2)-(fi(2)-fi(4))'//
     &                 '*ANS(',i,')*pdfwgt(2)'
               else
                  write(ln,11) 'finite(1)=finite(1)-'//
     &                 '(fi(1)+fi(3)+fi(5))*ANS(',i,')*pdfwgt(1)'
                  write(ln,11) 'finite(2)=finite(2)-'//
     &                 '(fi(2)-fi(4)-fi(6))*ANS(',i,')*pdfwgt(2)'
               endif
c
c quark-gluon splitting (i:quark, ij~: gluon) and
c gluon-quark splitting (i:gluon, ij~: gluon)
c
            elseif((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.0) .or.
     &             (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.1)) then
             write(ln,11) 'finite(1)=finite(1)-fi(1)'//
     &              '*ANS(',i,')*pdfwgt(1)'
            else
               write (*,*) 'Error #142 in dipolcolor'
               write (ln,*) 'Error #142 in dipolcolor'
            endif

          elseif(lepton.or.(photon.and..not.photonpdf)) then

            write(ln,75)'sik=two*dot(p',info_ijk(i,1,4),'(0,',info_ijk(i,1,4)
     &           ,'),p',info_ijk(i,1,4),'(0,',info_ijk(i,1,5),'))'
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(2))'
           if(info_ijk(i,2,1).eq.info_ijk(i,2,2).and..not.
     &          first(info_ijk(i,1,1),info_ijk(i,1,2))) then
             write(ln,15)'IF PASSCUTSDIP(P) THEN'
             if (i.le.9) then
                write(ln,61)' CALL SDIPOLQED'//char(48)//char(48)
     &            //char(48+i)//'(P,ANSZONE(',i,'))'
             elseif(i.le.99) then
                write(ln,61)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &            //char(mod(i,10)+48)//'(P,ANSZONE(',i,'))'
             elseif(i.le.999) then
                write(ln,61)' CALL SDIPOLQED'//char(i/100+48)
     &            //char((i-(i/100)*100)/10+48)
     &            //char(mod(i-(i/100)*100,10)+48)//'(P,ANSZONE(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             write(ln,15)'ENDIF'
c             write(ln,6)' CALL BOOST(P,Z,',info_ijk(i,1,1),',P1)'
             write(ln,9)'IF (PASSCUTSDIP(P',info_ijk(i,1,4),')) THEN'
             if (i.le.9) then
                write(ln,51)' CALL SDIPOLQED'//char(48)//char(48)
     &            //char(48+i)//'(P',info_ijk(i,1,4),',ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,51)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &            //char(mod(i,10)+48)//'(P',info_ijk(i,1,4),',ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,51)' CALL SDIPOLQED'//char(i/100+48)
     &            //char((i-(i/100)*100)/10+48)
     &            //char(mod(i-(i/100)*100,10)+48)//'(P',info_ijk(i,1,4),',ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             write(ln,15)'ENDIF'
           else
             write(ln,10)'IF (PASSCUTSDIP(P)) THEN'
             if (i.le.9) then
                write(ln,62)' CALL SDIPOLQED'//char(48)//char(48)//char(48+i)
     &               //'(P,ANSZONE(',i,'))'
             elseif(i.le.99) then
                write(ln,62)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &               //char(mod(i,10)+48)//'(P,ANSZONE(',i,'))'
             elseif(i.le.999) then
                write(ln,62)' CALL SDIPOLQED'//char(i/100+48)
     &               //char((i-(i/100)*100)/10+48)
     &               //char(mod(i-(i/100)*100,10)+48)//'(P,ANSZONE(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             write(ln,10) 'ENDIF'
c             write(ln,13)' CALL BOOST(P,Z,',info_ijk(i,1,1),',P1)'
             write(ln,7)'IF (PASSCUTSDIP(P',info_ijk(i,1,4),')) THEN'
             if (i.le.9) then
                write(ln,41)' CALL SDIPOLQED'//char(48)//char(48)//char(48+i)
     &               //'(P',info_ijk(i,1,4),',ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,41)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &               //char(mod(i,10)+48)//'(P',info_ijk(i,1,4),',ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,41)' CALL SDIPOLQED'//char(i/100+48)
     &               //char((i-(i/100)*100)/10+48)
     &               //char(mod(i-(i/100)*100,10)+48)//'(P',info_ijk(i,1,4),',ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             write(ln,10) 'ENDIF' 
             first(info_ijk(i,1,1),info_ijk(i,1,2))=.false.
           endif

              write(ln,19) ' CALL EPSIFQED(',info_ijk(i,2,1),',',info_ijk(i,2,4),
     &          ',',chargint(i,1),',',chargint(i,2),',',info_ijk(i,1,1),',sikzone,mass_k,esq,e)'
               write(ln,11) 'epssq=epssq-esq*ANSZONE(',i,')*pdfwgt(2)'
               write(ln,11) 'eps=eps-e*ANSZONE(',i,')*pdfwgt(2)'
             write(ln,10) ' sab=massint(1)**2+massint(2)**2+two*dot(p(0,1),p(0,2))'
             write(ln,74) ' sabt=massint(1)**2+massint(2)**2+two*dot(p'
     &                      ,info_ijk(i,1,4),'(0,1),p',info_ijk(i,1,4),'(0,2))'
            if (info_ijk(i,1,1).eq.1) then
             write(ln,77) ' CALL FINITEIFQED(mass_k,sik,sikzone,z,',info_ijk(i,2,1),
     &          ',',info_ijk(i,2,4),',',chargint(i,1),',',chargint(i,2),','
     &          ,leg,',',info_p(5,inverse(iline_copy(1,i))),',sab,sabt,fi)'
            elseif(info_ijk(i,1,1).eq.2) then
             write(ln,77) ' CALL FINITEIFQED(mass_k,sik,sikzone,z,',info_ijk(i,2,1),
     &          ',',info_ijk(i,2,4),',',chargint(i,1),',',chargint(i,2),','
     &          ,leg,',',info_p(5,inverse(iline_copy(2,i))),',sab,sabt,fi)'
            endif
c
c quark-quark splitting (i:quark, ij~: quark) and
c gluon-gluon splitting (i:gluon, ij~: gluon)
c
            if(((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.1).or.
     &           (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.0))) then
               if (mass(info_ijk(i,1,3)).eq.'ZERO') then
                  write(ln,80) 'finite(1)=finite(1)-(fi(1)+fi(3))'//
     &                 '*ANS(',i,')*pdfwgt(2)*pswgt(',info_ijk(i,1,4)+1,')'
                  write(ln,11) 'finite(2)=finite(2)-(fi(2)-fi(4))'//
     &                 '*ANSZONE(',i,')*pdfwgt(2)*pswgt(1)'
               write(ln,11) ' CALL WRITEHIST(P,-(fi(2)-fi(4))*ANSZONE(',i,')*pdfwgt(2)*pswgt(1))'
               write(ln,8) ' CALL WRITEHIST(P',info_ijk(i,1,4),',-(fi(1)+fi(3))/z*ANS(',i,
     &                      ')*pdfwgt(2)*pswgt(',info_ijk(i,1,4)+1,'))'
               else
                  write(ln,80) 'finite(1)=finite(1)-'//
     &                 '(fi(1)+fi(3)+fi(5))*ANS(',i,')*pdfwgt(2)*pswgt(',info_ijk(i,1,4)+1,')'
                  write(ln,11) 'finite(2)=finite(2)-'//
     &                 '(fi(2)-fi(4)-fi(6))*ANSZONE(',i,')*pdfwgt(2)*pswgt(1)'
               write(ln,11) ' CALL WRITEHIST(P,-(fi(2)-fi(4)-fi(6))*ANSZONE(',i,')*pdfwgt(2)*pswgt(1))'
               write(ln,8) ' CALL WRITEHIST(P',info_ijk(i,1,4),',-(fi(1)+fi(3)+fi(5))/z*ANS(',i,
     &                      ')*pdfwgt(2)*pswgt(',info_ijk(i,1,4)+1,'))'
               endif
c
c quark-gluon splitting (i:quark, ij~: gluon) and
c gluon-quark splitting (i:gluon, ij~: gluon)
c
            elseif((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.0) .or.
     &             (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.1)) then
             write(ln,80) 'finite(1)=finite(1)-fi(1)'//
     &              '*ANS(',i,')*pdfwgt(2)*pswgt(',info_ijk(i,1,4)+1,')'
             write(ln,8) ' CALL WRITEHIST(P',info_ijk(i,1,4),',-fi(1)/z*ANS(',i,
     &                      ')*pdfwgt(2)*pswgt(',info_ijk(i,1,4)+1,'))'
            else
               write (*,*) 'Error #142 in dipolcolor'
               write (ln,*) 'Error #142 in dipolcolor'
            endif


          endif


c**************************************************************
c     initial-initial theta function
c**************************************************************
         elseif(info_ijk(i,1,1).le.nincoming .and.
     &           info_ijk(i,1,3) .le.nincoming) then
            if(info_ijk(i,2,1).eq.0) then
              leg=info_ijk(i,1,2)
            elseif(info_ijk(i,2,2).eq.0) then
              leg=info_ijk(i,1,1)
            elseif(info_ijk(i,2,1).eq.info_ijk(i,2,2)) then
              leg=info_ijk(i,1,1)
            endif
            write(ln,10)' '
            write(ln,74)'sikzone=two*dot(p(0,',info_ijk(i,1,4),
     &           '),p(0,',info_ijk(i,1,5),'))'


          if(.not.lepton.and..not.photon.or.(photon.and.photonpdf)) then


            write(ln,74)'sik=two*dot(p(0,',info_ijk(i,1,4),
     &           '),p(0,',info_ijk(i,1,5),'))/z'
            if (info_ijk(i,1,1).eq.1) then
               write(ln,17)'CALL DIPOLEPDF(x(1)/z,x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(1))'
            elseif(info_ijk(i,1,1).eq.2) then
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2)/z,',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(1))'
            else
               write (ln,*)'Error #147 in write_intdipole'
               write (*,*)'Error #147 in write_intdipole'
            endif
            if(((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.1).or.
     &         (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.0))) then
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(2))'
            endif
            if(info_ijk(i,2,1).eq.info_ijk(i,2,2).and..not.
     &          first(info_ijk(i,1,1),info_ijk(i,1,2))) then
             if (i.le.9) then
                write(ln,51)'CALL SDIPOLQED'//char(48)//char(48)
     &            //char(48+i)//'(P,ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,51)'CALL SDIPOLQED'//char(48)//char(i/10+48)
     &            //char(mod(i,10)+48)//'(P,ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,51)'CALL SDIPOLQED'//char(i/100+48)
     &            //char((i-(i/100)*100)/10+48)
     &            //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
           else
             if (i.le.9) then
                write(ln,41)'CALL SDIPOLQED'//char(48)//char(48)
     &            //char(48+i)//'(P,ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,41)'CALL SDIPOLQED'//char(48)//char(i/10+48)
     &            //char(mod(i,10)+48)//'(P,ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,41)'CALL SDIPOLQED'//char(i/100+48)
     &            //char((i-(i/100)*100)/10+48)
     &            //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             first(info_ijk(i,1,1),info_ijk(i,1,2))=.false.
            endif

             write(ln,19) ' CALL EPSIIQED(',info_ijk(i,2,1),',',info_ijk(i,2,4),',',
     &             chargint(i,1),',',chargint(i,2),',',info_ijk(i,1,1),',sikzone,esq,e)'
               write(ln,11) 'epssq=epssq-esq*ANS(',i,')*pdfwgt(2)'
               write(ln,11) 'eps=eps-e*ANS(',i,')*pdfwgt(2)'
            if (info_ijk(i,1,1).eq.1) then
             write(ln,77) ' CALL FINITEIIQED(sik,sikzone,z,',info_ijk(i,2,1),','
     &            ,info_ijk(i,2,4),',',chargint(i,1),',',chargint(i,2),','
     &            ,leg,',',info_p(5,inverse(iline_copy(1,i))),',fi)'
            elseif(info_ijk(i,1,1).eq.2) then
             write(ln,77) ' CALL FINITEIIQED(sik,sikzone,z,',info_ijk(i,2,1),','
     &            ,info_ijk(i,2,4),',',chargint(i,1),',',chargint(i,2),','
     &            ,leg,',',info_p(5,inverse(iline_copy(2,i))),',fi)'
            endif
c
c quark-quark splitting (i:quark, ij~: quark) and
c gluon-gluon splitting (i:gluon, ij~: gluon)
c
            if(((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.1).or.
     &         (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.0))) then
               if (mass(info_ijk(i,1,3)).eq.'ZERO') then
                  write(ln,11) 'finite(1)=finite(1)-(fi(1)+fi(3))'//
     &                 '*ANS(',i,')*pdfwgt(1)'
                  write(ln,11) 'finite(2)=finite(2)-(fi(2)-fi(4))'//
     &                 '*ANS(',i,')*pdfwgt(2)'
               else
                  write(ln,11) 'finite(1)=finite(1)-'//
     &                 '(fi(1)+fi(3)+fi(5))*ANS(',i,')*pdfwgt(1)'
                  write(ln,11) 'finite(2)=finite(2)-'//
     &                 '(fi(2)-fi(4)-fi(6))*ANS(',i,')*pdfwgt(2)'
               endif
c
c quark-gluon splitting (i:quark, ij~: gluon) and
c gluon-quark splitting (i:gluon, ij~: gluon)
c
            elseif((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.0) .or.
     &             (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.1)) then
               write(ln,11)
     &              'finite(1)=finite(1)-fi(1)*ANS(',i,')*pdfwgt(1)'
            else
               write (*,*) 'Error #144 in dipolcolor'
               write (ln,*) 'Error #144 in dipolcolor'
            endif



         elseif(lepton.or.(photon.and..not.photonpdf)) then


            write(ln,75)'sik=two*dot(p',info_ijk(i,1,4),'(0,',info_ijk(i,1,4)
     &           ,'),p',info_ijk(i,1,4),'(0,',info_ijk(i,1,5),'))'
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(2))'
           if(info_ijk(i,2,1).eq.info_ijk(i,2,2).and..not.
     &          first(info_ijk(i,1,1),info_ijk(i,1,2))) then
             write(ln,15)'IF (PASSCUTSDIP(P)) THEN'
             if (i.le.9) then
                write(ln,61)' CALL SDIPOLQED'//char(48)//char(48)
     &            //char(48+i)//'(P,ANSZONE(',i,'))'
             elseif(i.le.99) then
                write(ln,61)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &            //char(mod(i,10)+48)//'(P,ANSZONE(',i,'))'
             elseif(i.le.999) then
                write(ln,61)' CALL SDIPOLQED'//char(i/100+48)
     &            //char((i-(i/100)*100)/10+48)
     &            //char(mod(i-(i/100)*100,10)+48)//'(P,ANSZONE(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             write(ln,15)'ENDIF'
c             write(ln,6)' CALL BOOST(P,Z,',info_ijk(i,1,1),',P1)'
             write(ln,9)'IF (PASSCUTSDIP(P',info_ijk(i,1,4),')) THEN'
             if (i.le.9) then
                write(ln,51)' CALL SDIPOLQED'//char(48)//char(48)
     &            //char(48+i)//'(P',info_ijk(i,1,4),',ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,51)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &            //char(mod(i,10)+48)//'(P',info_ijk(i,1,4),',ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,51)' CALL SDIPOLQED'//char(i/100+48)
     &            //char((i-(i/100)*100)/10+48)
     &            //char(mod(i-(i/100)*100,10)+48)//'(P',info_ijk(i,1,4),',ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             write(ln,15)'ENDIF'
           else
             write(ln,10)'IF (PASSCUTSDIP(P)) THEN'
             if (i.le.9) then
                write(ln,62)' CALL SDIPOLQED'//char(48)//char(48)//char(48+i)
     &               //'(P,ANSZONE(',i,'))'
             elseif(i.le.99) then
                write(ln,62)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &               //char(mod(i,10)+48)//'(P,ANSZONE(',i,'))'
             elseif(i.le.999) then
                write(ln,62)' CALL SDIPOLQED'//char(i/100+48)
     &               //char((i-(i/100)*100)/10+48)
     &               //char(mod(i-(i/100)*100,10)+48)//'(P,ANSZONE(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             write(ln,10) 'ENDIF'
c             write(ln,13)' CALL BOOST(P,Z,',info_ijk(i,1,1),',P1)'
             write(ln,7)'IF (PASSCUTSDIP(P',info_ijk(i,1,4),')) THEN'
             if (i.le.9) then
                write(ln,41)' CALL SDIPOLQED'//char(48)//char(48)//char(48+i)
     &               //'(P',info_ijk(i,1,4),',ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,41)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &               //char(mod(i,10)+48)//'(P',info_ijk(i,1,4),',ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,41)' CALL SDIPOLQED'//char(i/100+48)
     &               //char((i-(i/100)*100)/10+48)
     &               //char(mod(i-(i/100)*100,10)+48)//'(P',info_ijk(i,1,4),',ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             write(ln,10) 'ENDIF' 
             first(info_ijk(i,1,1),info_ijk(i,1,2))=.false.
           endif

             write(ln,19) ' CALL EPSIIQED(',info_ijk(i,2,1),',',info_ijk(i,2,4),',',
     &             chargint(i,1),',',chargint(i,2),',',info_ijk(i,1,1),',sikzone,esq,e)'
               write(ln,11) 'epssq=epssq-esq*ANSZONE(',i,')*pdfwgt(2)'
               write(ln,11) 'eps=eps-e*ANSZONE(',i,')*pdfwgt(2)'
            if (info_ijk(i,1,1).eq.1) then
             write(ln,77) ' CALL FINITEIIQED(sik,sikzone,z,',info_ijk(i,2,1),','
     &            ,info_ijk(i,2,4),',',chargint(i,1),',',chargint(i,2),','
     &            ,leg,',',info_p(5,inverse(iline_copy(1,i))),',fi)'
            elseif(info_ijk(i,1,1).eq.2) then
             write(ln,77) ' CALL FINITEIIQED(sik,sikzone,z,',info_ijk(i,2,1),','
     &            ,info_ijk(i,2,4),',',chargint(i,1),',',chargint(i,2),','
     &            ,leg,',',info_p(5,inverse(iline_copy(2,i))),',fi)'
            endif
c
c quark-quark splitting (i:quark, ij~: quark) and
c gluon-gluon splitting (i:gluon, ij~: gluon)
c
            if(((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.1).or.
     &         (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.0))) then
               if (mass(info_ijk(i,1,3)).eq.'ZERO') then
                  write(ln,80) 'finite(1)=finite(1)-(fi(1)+fi(3))'//
     &                 '*ANS(',i,')*pdfwgt(2)*pswgt(',info_ijk(i,1,4)+1,')'
                  write(ln,11) 'finite(2)=finite(2)-(fi(2)-fi(4))'//
     &                 '*ANSZONE(',i,')*pdfwgt(2)*pswgt(1)'
                  write(ln,11) ' CALL WRITEHIST(P,-(fi(2)-fi(4))*ANSZONE(',i,')*pdfwgt(2)*pswgt(1))'
                  write(ln,8) ' CALL WRITEHIST(P',info_ijk(i,1,4),',-(fi(1)+fi(3))/z*ANS(',i,
     &                      ')*pdfwgt(2)*pswgt(',info_ijk(i,1,4)+1,'))'
               else
                  write(ln,80) 'finite(1)=finite(1)-'//
     &                 '(fi(1)+fi(3)+fi(5))*ANS(',i,')*pdfwgt(2)*pswgt(',info_ijk(i,1,4)+1,')'
                  write(ln,11) 'finite(2)=finite(2)-'//
     &                 '(fi(2)-fi(4)-fi(6))*ANSZONE(',i,')*pdfwgt(2)*pswgt(1)'
                  write(ln,11) ' CALL WRITEHIST(P,-(fi(2)-fi(4)-fi(6))*ANSZONE(',i,')*pdfwgt(2)*pswgt(1))'
                  write(ln,8) ' CALL WRITEHIST(P',info_ijk(i,1,4),',-(fi(1)+fi(3)+fi(5))/z*ANS(',i,
     &                      ')*pdfwgt(2)*pswgt(',info_ijk(i,1,4)+1,'))'
               endif
c
c quark-gluon splitting (i:quark, ij~: gluon) and
c gluon-quark splitting (i:gluon, ij~: gluon)
c
            elseif((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.0) .or.
     &             (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.1)) then
               write(ln,80)
     &              'finite(1)=finite(1)-fi(1)*ANS(',i,')*pdfwgt(2)*pswgt(',info_ijk(i,1,4)+1,')'
               write(ln,8) ' CALL WRITEHIST(P',info_ijk(i,1,4),',-fi(1)/z*ANS(',i,
     &                      ')*pdfwgt(2)*pswgt(',info_ijk(i,1,4)+1,'))'
            else
               write (*,*) 'Error #144 in dipolcolor'
               write (ln,*) 'Error #144 in dipolcolor'
            endif








         endif


         else
            write (*,*) 'Error #143 in dipolcolor'
            write (ln,*) 'Error #143 in dipolcolor'
         endif
         endif ! massless unresolved
      enddo
      write(ln,10) ' '
      if(.not.lepton.and..not.photon.or.(photon.and.photonpdf)) then
       write(ln,10) 'ENDIF' 
      endif
      if(.not.lepton.and..not.photon.or.(photon.and.photonpdf)) then
       write(ln,10) 'FIN=(FINITE(1)/z+FINITE(2))*pswgt(1)'
       write(ln,10) 'CALL WRITEHIST(P,FIN)'
      else
       write(ln,10) 'FIN=(FINITE(1)/z+FINITE(2))'
      endif
      write(ln,10) 'INTDIPOLE=.FALSE.'
      write(ln,10) 'RETURN'
      write(ln,10) 'END'
      write(ln,10) ' '
      write(ln,10) ' '
      write(ln,10) ' '
      write(ln,10) ' '
      write(ln,10) ' '

      write(ln,10)
     &     'SUBROUTINE INTDIPOLESFINITEQED(P,P1,P2,X,Z,PSWGT,EPS,EPSSQ,FIN)'
      write(ln,5) ' '
      write(ln,5) cversion
      write(ln,5) 
     .     'RETURNS COEFFICIENTS OF 1/eps^2,1/eps AND FINITE TERMS'
      write(ln,5) 'FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL) '
      write(ln,5) ' '
      write(ln,5) 'FOR PROCESS '//proc(1:iproc)
      write(ln,5) 'WHERE Z IS THE MOMENTUM FRACTION '//
     &     'OF THE INITAL EMITTER/SPECTATOR.'
      write(ln,5) ' '

      write(ln,10) 'IMPLICIT NONE'
      write(ln,10) 'include "nexternal.inc" '
      write(ln,10) 'include "dipole.inc" '
      write(ln,10) 'include "coupl.inc"'

      write(ln,5) ' '
      write(ln,5) 'ARGUMENTS '
      write(ln,5) ' '
      write(ln,10) 'REAL*8 P(0:3,NEXTERNAL),P1(0:3,NEXTERNAL),P2(0:3,NEXTERNAL),X(2),Z,PSWGT(3)'
      write(ln,10) 'REAL*8 EPS,EPSSQ,FIN'
      write(ln,5) ' '
      write(ln,5) 'GLOBAL '
      write(ln,5) ' '
      write(ln,10) 'REAL*8 dot,lambda_tr'
      write(ln,10) 'LOGICAL PASSCUTSDIP'
      write(ln,10) 'EXTERNAL dot,lambda_tr,PASSCUTSDIP'
      write(ln,5) ' '
      write(ln,5) 'LOCAL'
      write(ln,5) ' '
      write(ln,10) 'REAL*8 sij,sik,sjk,zp,e,esq,pdfwgt(3),FINITE(3)'
      write(ln,10) 'REAL*8 y_ijk,omx_ija,u_i,vtilde_i,zplus'
      write(ln,10) 'REAL*8 mass_i,mass_j,mass_k,mass_ij'
      write(ln,10) 'REAL*8 msq_i,msq_j,msq_k,msq_ij,sikzone'
      write(ln,10) 'REAL*8 musq_i,musq_j,musq_k,musq_ij'
      write(ln,10) 'REAL*8 vij,vtildeij,rho_i,rho_k,rho,qsq'
      write(ln,10) 'REAL*8 one,two,three,four,cf,ca,qtot(0:3),pi'
      write(ln,10) 'REAL*8 gsq,gammaE,ff,fi(9),sab,sabt'
      write(ln,10) 'parameter (one=1d0,two=2d0,three=3d0,four=4d0)'
c      write(ln,10) 'parameter (cf=4d0/3d0,ca=3d0)'
      write(ln,60) 'REAL*8 ANS(',nrdipol,')'
      if(lepton.or.(photon.and..not.photonpdf)) then
       write(ln,60) 'REAL*8 ANSZONE(',nrdipol,')'
      endif
      write(ln,10) 'REAL*8 ZERO'
      write(ln,10) 'INTEGER NU'
      write(ln,10) 'PARAMETER (ZERO=0D0)'
      write(ln,10) 'PARAMETER (gammaE=0.57721566d0)'
      write(ln,10) 'PARAMETER (pi=3.1415926d0)'
      write(ln,5) ' '
      write(ln,5) 'BEGIN CODE'
      write(ln,5) ' '
      write(ln,10) 'INTDIPOLE=.TRUE.'
c      write(ln,10) 'gsq=GG(1)**2'
      write(ln,10) 'eps=0d0'
      write(ln,10) 'epssq=0d0'
      write(ln,10) 'fin=0d0'
      write(ln,10) 'e=0d0'
      write(ln,10) 'esq=0d0'
      write(ln,10) 'zplus=1d0'
      write(ln,10) 'DO NU=1,3'
      write(ln,10) '   finite(NU)=0d0'
      write(ln,10) 'ENDDO'
      write(ln,60) 'DO NU=1,',nrdipol
      write(ln,10) '   ANS(NU)=0d0'
      if(lepton.or.(photon.and.(photonpdf.eqv..false.))) then
       write(ln,10) '   ANSZONE(NU)=0d0'
      endif
      write(ln,10) 'ENDDO'
      if(.not.lepton.and..not.photon.or.(photon.and.photonpdf)) then
       write(ln,10) 'IF (PASSCUTSDIP(P)) THEN'
      endif 
      do i=1,nr
c     First check if unresolved is massive. If not do print it here.
         if (mass(info_ijk(i,1,2)).ne.'ZERO')then
c**************************************************************
c       final-final theta funtion
c**************************************************************
         if(info_ijk(i,1,1).gt.nincoming
     &        .and. info_ijk(i,1,3).gt.nincoming) then
            if(info_ijk(i,2,1).eq.0) then
              leg=info_ijk(i,1,2)
            elseif(info_ijk(i,2,2).eq.0) then
              leg=info_ijk(i,1,1)
            elseif(info_ijk(i,2,1).eq.info_ijk(i,2,2)) then
              leg=info_ijk(i,1,1)
            endif
            write(ln,10)' '
            write(ln,71)'mass_i=',mass(info_ijk(i,1,1))
            write(ln,71)'mass_k=',mass(info_ijk(i,1,3))
            write(ln,10) 'DO NU=0,3'
            write(ln,73) ' qtot(NU)=p(NU,',info_ijk(i,1,4),')+p(NU,',info_ijk(i,1,5),')'
c            write(ln,10) ' IF (NCS) THEN'
c            write(ln,76) ' qtot(NU)=p(NU,',info_ijk(i,1,1),')+p(NU,',info_ijk(i,1,2),
c     &                    ')+p(NU,',info_ijk(i,1,3),')'
c            write(ln,10) ' ENDIF'
            write(ln,10) 'ENDDO'
            write(ln,10)'qsq=dot(qtot,qtot)'
            write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &           info_p(5,inverse(iline_copy(1,i))),',',
     &           info_p(5,inverse(iline_copy(2,i))),',pdfwgt(2))'
           if(info_ijk(i,2,1).eq.info_ijk(i,2,2).and..not.
     &          first(info_ijk(i,1,1),info_ijk(i,1,2))) then
             write(ln,15)'IF (PASSCUTSDIP(P)) THEN'
             if (i.le.9) then
                write(ln,51)' CALL SDIPOLQED'//char(48)//char(48)
     &            //char(48+i)//'(P,ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,51)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &            //char(mod(i,10)+48)//'(P,ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,51)' CALL SDIPOLQED'//char(i/100+48)
     &            //char((i-(i/100)*100)/10+48)
     &            //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             write(ln,15)'ENDIF'
           else
             write(ln,10)'IF (PASSCUTSDIP(P)) THEN'
             if (i.le.9) then
                write(ln,41)
     &           ' CALL SDIPOLQED'//char(48)//char(48)//char(48+i)
     &               //'(P,ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,41)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &               //char(mod(i,10)+48)//'(P,ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,41)' CALL SDIPOLQED'//char(i/100+48)
     &               //char((i-(i/100)*100)/10+48)
     &               //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             first(info_ijk(i,1,1),info_ijk(i,1,2))=.false.
             write(ln,10) 'ENDIF'
           endif
            write(ln,77) ' zi=dot(p(0,',info_ijk(i,1,1),'),p(0,',info_ijk(i,1,3),
     &           '))/(dot(p(0,',info_ijk(i,1,2),'),p(0,',info_ijk(i,1,3),'))+dot(p(0,'
     &            ,info_ijk(i,1,1),'),p(0,',info_ijk(i,1,3),')))'
            write(ln,18) ' CALL EPSFFQED(mass_i,mass_k,qsq,zi,',info_ijk(i,2,4),
     &           ',',info_ijk(i,2,1),',',chargint(i,1),',',chargint(i,2),',esq,e)'
            write(ln,11) 'epssq=epssq-esq*ANS(',i,')*pdfwgt(2)'
            write(ln,11) 'eps=eps-e*ANS(',i,')*pdfwgt(2)'
            write(ln,19) ' CALL FINITEFFQED(mass_i,mass_k,qsq,zi,',info_ijk(i,2,4),
     &           ',',info_ijk(i,2,1),',',chargint(i,1),',',chargint(i,2),',',leg,',fi)'
            if(lepton.or.(photon.and..not.photonpdf)) then   
             write(ln,11) ' finite(2)=finite(2)-fi(2)*ANS(',i,')*pdfwgt(2)*pswgt(1)'
             write(ln,11 ) ' CALL WRITEHIST(P,-fi(2)*ANS(',i,')*pdfwgt(2)*pswgt(1))'
            else
             write(ln,11) ' finite(2)=finite(2)-fi(2)*ANS(',i,')*pdfwgt(2)'
            endif

c**************************************************************
c       final-initial theta function
c**************************************************************
         elseif(info_ijk(i,1,1).gt.nincoming
     &           .and.info_ijk(i,1,3).le.nincoming) then
            if(info_ijk(i,2,1).eq.0) then
              leg=info_ijk(i,1,2)
            elseif(info_ijk(i,2,2).eq.0) then
              leg=info_ijk(i,1,1)
            elseif(info_ijk(i,2,1).eq.info_ijk(i,2,2)) then
              leg=info_ijk(i,1,1)
            endif
            write(ln,10)' '
            write(ln,71)'mass_i=',mass(info_ijk(i,1,1))
            if(info_ijk(i,2,4).eq.0 .and. (
     &           mass(info_ijk(i,1,1)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,2)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,3)).ne.'ZERO')) then
               write (ln,11) 'zplus=1d0-4d0*mass_i**2/sik'
               write(ln,74)'sikzone=two*dot(p(0,',info_ijk(i,1,4),
     &              '),p(0,',info_ijk(i,1,5),
     &              '))/zplus'
            else
               write(ln,74)'sikzone=two*dot(p(0,',info_ijk(i,1,4),
     &              '),p(0,',info_ijk(i,1,5),'))'
            endif

           if(.not.lepton.and..not.photon.or.(photon.and.photonpdf)) then


            write(ln,74)'sik=two*dot(p(0,',info_ijk(i,1,4),'),p(0,'
     &           ,info_ijk(i,1,5),'))/z'
            if(info_ijk(i,1,3).eq.1) then
               write(ln,17)'CALL DIPOLEPDF(x(1)/z,x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(1))'
            elseif(info_ijk(i,1,3).eq.2) then
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2)/z,',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(1))'
            else
               write (ln,*)'Error #145 in write_intdipole'
               write (*,*)'Error #145 in write_intdipole'
            endif
            if(info_ijk(i,2,4).eq.0 .and. (
     &           mass(info_ijk(i,1,1)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,2)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,3)).ne.'ZERO')) then
               if(info_ijk(i,1,3).eq.1) then
                  write(ln,17)'CALL DIPOLEPDF(x(1)/zplus,x(2),',
     &                 info_p(5,inverse(iline_copy(1,i))),',',
     &                 info_p(5,inverse(iline_copy(2,i))),',pdfwgt(3))'
               elseif(info_ijk(i,1,3).eq.2) then
                  write(ln,17)'CALL DIPOLEPDF(x(1),x(2)/zplus,',
     &                 info_p(5,inverse(iline_copy(1,i))),',',
     &                 info_p(5,inverse(iline_copy(2,i))),',pdfwgt(3))'
               endif
            else
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(2))'
            endif
           if(info_ijk(i,2,1).eq.info_ijk(i,2,2).and..not.
     &          first(info_ijk(i,1,1),info_ijk(i,1,2))) then
             if (i.le.9) then
                write(ln,51)'CALL SDIPOLQED'//char(48)//char(48)
     &            //char(48+i)//'(P,ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,51)'CALL SDIPOLQED'//char(48)//char(i/10+48)
     &            //char(mod(i,10)+48)//'(P,ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,51)'CALL SDIPOLQED'//char(i/100+48)
     &            //char((i-(i/100)*100)/10+48)
     &            //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
           else
             if (i.le.9) then
                write(ln,41)'CALL SDIPOLQED'//char(48)//char(48)//char(48+i)
     &               //'(P,ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,41)'CALL SDIPOLQED'//char(48)//char(i/10+48)
     &               //char(mod(i,10)+48)//'(P,ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,41)'CALL SDIPOLQED'//char(i/100+48)
     &               //char((i-(i/100)*100)/10+48)
     &               //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             first(info_ijk(i,1,1),info_ijk(i,1,2))=.false.
           endif


             write(ln,18)' CALL EPSFIQED(',info_ijk(i,2,4),',',info_ijk(i,2,1),
     &         ',',chargint(i,1),',',chargint(i,2),',sikzone,mass_i,zi,esq,e)'
            write(ln,11) 'epssq=epssq-esq*ANS(',i,')*pdfwgt(2)'
            write(ln,11) 'eps=eps-e*ANS(',i,')*pdfwgt(2)'
             write(ln,10) ' sab=massint(1)**2+massint(2)**2+two*dot(p(0,1),p(0,2))'
             write(ln,10) ' sabt=massint(1)**2+massint(2)**2+two*dot(p(0,1),p(0,2))*z'
c            if (info_ijk(i,1,1).eq.1) then
             write(ln,19) ' CALL FINITEFIQED(mass_i,sik,sikzone,z,',info_ijk(i,2,4),
     &         ',',info_ijk(i,2,1),',',chargint(i,1),',',chargint(i,2),
     &         ',',leg,',sab,sabt,fi)'
c            elseif(info_ijk(i,1,1).eq.2) then
c             write(ln,19) ' CALL FINITEFIQED(mass_i,sik,sikzone,z,',info_ijk(i,2,4),
c     &         ',',info_ijk(i,2,1),',',chargint(i,1),',',chargint(i,2),
c     &         ',',info_ijk(i,1,3),',sab,sabt,fi)'
c            endif
            if(info_ijk(i,2,4).eq.0 .and. (
     &           mass(info_ijk(i,1,1)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,2)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,3)).ne.'ZERO')) then
               write(ln,11)
     &            'finite(2)=finite(2)-(-fi(8))*ANS(',i,')*pdfwgt(2)'
               write(ln,11) 'finite(3)=finite(3)-(fi(7)+fi(9))'//
     &              '*ANS(',i,')*pdfwgt(3)'
            else
               write(ln,11) 'finite(1)=finite(1)-(fi(1)+fi(3))'//
     &              '*ANS(',i,')*pdfwgt(1)'
               write(ln,11) 'finite(2)=finite(2)-(fi(2)-fi(4))'//
     &              '*ANS(',i,')*pdfwgt(2)'
            endif


           elseif(lepton.or.(photon.and..not.photonpdf)) then


            write(ln,75)'sik=two*dot(p',info_ijk(i,1,5),'(0,',info_ijk(i,1,4),
     &              '),p',info_ijk(i,1,5),'(0,',info_ijk(i,1,5),'))'
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(2))'
           if(info_ijk(i,2,1).eq.info_ijk(i,2,2).and..not.
     &          first(info_ijk(i,1,1),info_ijk(i,1,2))) then
             write(ln,15)'IF (PASSCUTSDIP(P)) THEN'
             if (i.le.9) then
                write(ln,61)' CALL SDIPOLQED'//char(48)//char(48)
     &            //char(48+i)//'(P,ANSZONE(',i,'))'
             elseif(i.le.99) then
                write(ln,61)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &            //char(mod(i,10)+48)//'(P,ANSZONE(',i,'))'
             elseif(i.le.999) then
                write(ln,61)' CALL SDIPOLQED'//char(i/100+48)
     &            //char((i-(i/100)*100)/10+48)
     &            //char(mod(i-(i/100)*100,10)+48)//'(P,ANSZONE(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             write(ln,15)'ENDIF'
c             write(ln,6)' CALL BOOST(P,Z,',info_ijk(i,1,3),',P1)'
             write(ln,9)'IF (PASSCUTSDIP(P',info_ijk(i,1,5),')) THEN'
             if (i.le.9) then
                write(ln,51)' CALL SDIPOLQED'//char(48)//char(48)
     &            //char(48+i)//'(P',info_ijk(i,1,5),',ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,51)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &            //char(mod(i,10)+48)//'(P',info_ijk(i,1,5),',ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,51)' CALL SDIPOLQED'//char(i/100+48)
     &            //char((i-(i/100)*100)/10+48)
     &            //char(mod(i-(i/100)*100,10)+48)//'(P',info_ijk(i,1,5),',ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             write(ln,15)'ENDIF'
           else
             write(ln,10)'IF (PASSCUTSDIP(P)) THEN'
             if (i.le.9) then
                write(ln,62)' CALL SDIPOLQED'//char(48)//char(48)//char(48+i)
     &               //'(P,ANSZONE(',i,'))'
             elseif(i.le.99) then
                write(ln,62)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &               //char(mod(i,10)+48)//'(P,ANSZONE(',i,'))'
             elseif(i.le.999) then
                write(ln,62)' CALL SDIPOLQED'//char(i/100+48)
     &               //char((i-(i/100)*100)/10+48)
     &               //char(mod(i-(i/100)*100,10)+48)//'(P,ANSZONE(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             write(ln,10) 'ENDIF'
c             write(ln,13)' CALL BOOST(P,Z,',info_ijk(i,1,3),',P1)'
             write(ln,7)'IF (PASSCUTSDIP(P',info_ijk(i,1,5),')) THEN'
             if (i.le.9) then
                write(ln,41)' CALL SDIPOLQED'//char(48)//char(48)//char(48+i)
     &               //'(P',info_ijk(i,1,5),',ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,41)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &               //char(mod(i,10)+48)//'(P',info_ijk(i,1,5),',ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,41)' CALL SDIPOLQED'//char(i/100+48)
     &               //char((i-(i/100)*100)/10+48)
     &               //char(mod(i-(i/100)*100,10)+48)//'(P',info_ijk(i,1,5),',ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             write(ln,10) 'ENDIF' 
             first(info_ijk(i,1,1),info_ijk(i,1,2))=.false.
           endif

             write(ln,18)' CALL EPSFIQED(',info_ijk(i,2,4),',',info_ijk(i,2,1),
     &         ',',chargint(i,1),',',chargint(i,2),',sikzone,mass_i,zi,esq,e)'
            write(ln,11) 'epssq=epssq-esq*ANSZONE(',i,')*pdfwgt(2)'
            write(ln,11) 'eps=eps-e*ANSZONE(',i,')*pdfwgt(2)'
             write(ln,10) ' sab=massint(1)**2+massint(2)**2+two*dot(p(0,1),p(0,2))'
             write(ln,74) ' sabt=massint(1)**2+massint(2)**2+two*dot(p'
     &                      ,info_ijk(i,1,5),'(0,1),p',info_ijk(i,1,5),'(0,2))'
c            if (info_ijk(i,1,1).eq.1) then
c             write(ln,19) ' CALL FINITEFIQED(mass_i,sik,sikzone,z,',info_ijk(i,2,4),
c     &         ',',info_ijk(i,2,1),',',chargint(i,1),',',chargint(i,2),
c     &         ',',info_ijk(i,1,3),',sab,sabt,fi)'
c            elseif(info_ijk(i,1,1).eq.2) then
             write(ln,19) ' CALL FINITEFIQED(mass_i,sik,sikzone,z,',info_ijk(i,2,4),
     &         ',',info_ijk(i,2,1),',',chargint(i,1),',',chargint(i,2),
     &         ',',leg,',sab,sabt,fi)'
c            endif
            if(info_ijk(i,2,4).eq.0 .and. (
     &           mass(info_ijk(i,1,1)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,2)).ne.'ZERO'.or.
     &           mass(info_ijk(i,1,3)).ne.'ZERO')) then
               write(ln,11)
     &            'finite(2)=finite(2)-(-fi(8))*ANSZONE(',i,')*pdfwgt(2)*pswgt(1)'
               write(ln,80) 'finite(3)=finite(3)-(fi(7)+fi(9))'//
     &              '*ANS(',i,')*pdfwgt(3)*pswgt(',info_ijk(i,1,5)+1,')'
               write(ln,11) ' CALL WRITEHIST(P,fi(8)*ANSZONE(',i,')*pdfwgt(2)*pswgt(1))'
               write(ln,8) ' CALL WRITEHIST(P',info_ijk(i,1,5),',-(fi(7)+fi(9))*ANS(',i,
     &                      ')*pdfwgt(3)*pswgt(',info_ijk(i,1,5)+1,'))'
            else
               write(ln,80) 'finite(1)=finite(1)-(fi(1)+fi(3))'//
     &              '*ANS(',i,')*pdfwgt(2)*pswgt(',info_ijk(i,1,5)+1,')'
               write(ln,11) 'finite(2)=finite(2)-(fi(2)-fi(4))'//
     &              '*ANSZONE(',i,')*pdfwgt(2)*pswgt(1)'
               write(ln,11) ' CALL WRITEHIST(P,-(fi(2)-fi(4))*ANSZONE(',i,')*pdfwgt(2)*pswgt(1))'
               write(ln,8) ' CALL WRITEHIST(P',info_ijk(i,1,5),',-(fi(1)+fi(3))/z*ANS(',i,
     &                      ')*pdfwgt(3)*pswgt(',info_ijk(i,1,5)+1,'))'
            endif

          

           endif

c**************************************************************
c       initial-final theta function
c**************************************************************
         elseif(info_ijk(i,1,1).le.nincoming .and.
     &           info_ijk(i,1,3).gt.nincoming) then
            if(info_ijk(i,2,1).eq.0) then
              leg=info_ijk(i,1,2)
            elseif(info_ijk(i,2,2).eq.0) then
              leg=info_ijk(i,1,1)
            elseif(info_ijk(i,2,1).eq.info_ijk(i,2,2)) then
              leg=info_ijk(i,1,1)
            endif
            write(ln,10)' '
            write(ln,71)'mass_k=',mass(info_ijk(i,1,3))
            write(ln,74)'sikzone=two*dot(p(0,',info_ijk(i,1,4)
     &           ,'),p(0,',info_ijk(i,1,5),'))'


          if(.not.lepton.and..not.photon.or.(photon.and.photonpdf)) then
  

            write(ln,74)'sik=two*dot(p(0,',info_ijk(i,1,4)
     &           ,'),p(0,',info_ijk(i,1,5),'))/z'
            if (info_ijk(i,1,1).eq.1) then
               write(ln,17)'CALL DIPOLEPDF(x(1)/z,x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(1))'
            elseif(info_ijk(i,1,1).eq.2) then
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2)/z,',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(1))'
            else
               write (ln,*)'Error #146 in write_intdipole'
               write (*,*)'Error #146 in write_intdipole'
            endif
            if(((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.1).or.
     &           (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.0))) then
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(2))'
            endif
          if(info_ijk(i,2,1).eq.info_ijk(i,2,2).and..not.
     &          first(info_ijk(i,1,1),info_ijk(i,1,2))) then
             if (i.le.9) then
                write(ln,51)'CALL SDIPOLQED'//char(48)//char(48)
     &            //char(48+i)//'(P,ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,51)'CALL SDIPOLQED'//char(48)//char(i/10+48)
     &            //char(mod(i,10)+48)//'(P,ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,51)'CALL SDIPOLQED'//char(i/100+48)
     &            //char((i-(i/100)*100)/10+48)
     &            //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
           else
             if (i.le.9) then
                write(ln,41)'CALL SDIPOLQED'//char(48)//char(48)
     &               //char(48+i)//'(P,ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,41)'CALL SDIPOLQED'//char(48)//char(i/10+48)
     &               //char(mod(i,10)+48)//'(P,ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,41)'CALL SDIPOLQED'//char(i/100+48)
     &               //char((i-(i/100)*100)/10+48)
     &               //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             first(info_ijk(i,1,1),info_ijk(i,1,2))=.false.
           endif

              write(ln,19) ' CALL EPSIFQED(',info_ijk(i,2,1),',',info_ijk(i,2,4),
     &          ',',chargint(i,1),',',chargint(i,2),',',info_ijk(i,1,1),',sikzone,mass_k,esq,e)'
               write(ln,11) 'epssq=epssq-esq*ANS(',i,')*pdfwgt(2)'
               write(ln,11) 'eps=eps-e*ANS(',i,')*pdfwgt(2)'
             write(ln,10) ' sab=massint(1)**2+massint(2)**2+two*dot(p(0,1),p(0,2))/z'
             write(ln,10) ' sabt=massint(1)**2+massint(2)**2+two*dot(p(0,1),p(0,2))'
            if (info_ijk(i,1,1).eq.1) then
             write(ln,77) ' CALL FINITEIFQED(mass_k,sik,sikzone,z,',info_ijk(i,2,1),
     &          ',',info_ijk(i,2,4),',',chargint(i,1),',',chargint(i,2),','
     &          ,leg,',',info_p(5,inverse(iline_copy(1,i))),',sab,sabt,fi)'
            elseif(info_ijk(i,1,1).eq.2) then
             write(ln,77) ' CALL FINITEIFQED(mass_k,sik,sikzone,z,',info_ijk(i,2,1),
     &          ',',info_ijk(i,2,4),',',chargint(i,1),',',chargint(i,2),','
     &          ,leg,',',info_p(5,inverse(iline_copy(2,i))),',sab,sabt,fi)'
            endif
c
c quark-quark splitting (i:quark, ij~: quark) and
c gluon-gluon splitting (i:gluon, ij~: gluon)
c
            if(((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.1).or.
     &           (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.0))) then
               if (mass(info_ijk(i,1,3)).eq.'ZERO') then
                  write(ln,11) 'finite(1)=finite(1)-(fi(1)+fi(3))'//
     &                 '*ANS(',i,')*pdfwgt(1)'
                  write(ln,11) 'finite(2)=finite(2)-(fi(2)-fi(4))'//
     &                 '*ANS(',i,')*pdfwgt(2)'
               else
                  write(ln,11) 'finite(1)=finite(1)-'//
     &                 '(fi(1)+fi(3)+fi(5))*ANS(',i,')*pdfwgt(1)'
                  write(ln,11) 'finite(2)=finite(2)-'//
     &                 '(fi(2)-fi(4)-fi(6))*ANS(',i,')*pdfwgt(2)'
               endif
c
c quark-gluon splitting (i:quark, ij~: gluon) and
c gluon-quark splitting (i:gluon, ij~: gluon)
c
            elseif((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.0) .or.
     &             (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.1)) then
             write(ln,11) 'finite(1)=finite(1)-fi(1)'//
     &              '*ANS(',i,')*pdfwgt(1)'
            else
               write (*,*) 'Error #142 in dipolcolor'
               write (ln,*) 'Error #142 in dipolcolor'
            endif

          elseif(lepton.or.(photon.and..not.photonpdf)) then


            write(ln,75)'sik=two*dot(p',info_ijk(i,1,4),'(0,',info_ijk(i,1,4)
     &           ,'),p',info_ijk(i,1,4),'(0,',info_ijk(i,1,5),'))'
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(2))'
           if(info_ijk(i,2,1).eq.info_ijk(i,2,2).and..not.
     &          first(info_ijk(i,1,1),info_ijk(i,1,2))) then
             write(ln,15)'IF PASSCUTSDIP(P) THEN'
             if (i.le.9) then
                write(ln,61)' CALL SDIPOLQED'//char(48)//char(48)
     &            //char(48+i)//'(P,ANSZONE(',i,'))'
             elseif(i.le.99) then
                write(ln,61)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &            //char(mod(i,10)+48)//'(P,ANSZONE(',i,'))'
             elseif(i.le.999) then
                write(ln,61)' CALL SDIPOLQED'//char(i/100+48)
     &            //char((i-(i/100)*100)/10+48)
     &            //char(mod(i-(i/100)*100,10)+48)//'(P,ANSZONE(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             write(ln,15)'ENDIF'
c             write(ln,6)' CALL BOOST(P,Z,',info_ijk(i,1,1),',P1)'
             write(ln,9)'IF (PASSCUTSDIP(P',info_ijk(i,1,4),')) THEN'
             if (i.le.9) then
                write(ln,51)' CALL SDIPOLQED'//char(48)//char(48)
     &            //char(48+i)//'(P',info_ijk(i,1,4),',ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,51)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &            //char(mod(i,10)+48)//'(P',info_ijk(i,1,4),',ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,51)' CALL SDIPOLQED'//char(i/100+48)
     &            //char((i-(i/100)*100)/10+48)
     &            //char(mod(i-(i/100)*100,10)+48)//'(P',info_ijk(i,1,4),',ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             write(ln,15)'ENDIF'
           else
             write(ln,10)'IF (PASSCUTSDIP(P)) THEN'
             if (i.le.9) then
                write(ln,62)' CALL SDIPOLQED'//char(48)//char(48)//char(48+i)
     &               //'(P,ANSZONE(',i,'))'
             elseif(i.le.99) then
                write(ln,62)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &               //char(mod(i,10)+48)//'(P,ANSZONE(',i,'))'
             elseif(i.le.999) then
                write(ln,62)' CALL SDIPOLQED'//char(i/100+48)
     &               //char((i-(i/100)*100)/10+48)
     &               //char(mod(i-(i/100)*100,10)+48)//'(P,ANSZONE(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             write(ln,10) 'ENDIF'
c             write(ln,13)' CALL BOOST(P,Z,',info_ijk(i,1,1),',P1)'
             write(ln,7)'IF (PASSCUTSDIP(P',info_ijk(i,1,4),')) THEN'
             if (i.le.9) then
                write(ln,41)' CALL SDIPOLQED'//char(48)//char(48)//char(48+i)
     &               //'(P',info_ijk(i,1,4),',ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,41)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &               //char(mod(i,10)+48)//'(P',info_ijk(i,1,4),',ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,41)' CALL SDIPOLQED'//char(i/100+48)
     &               //char((i-(i/100)*100)/10+48)
     &               //char(mod(i-(i/100)*100,10)+48)//'(P',info_ijk(i,1,4),',ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             write(ln,10) 'ENDIF' 
             first(info_ijk(i,1,1),info_ijk(i,1,2))=.false.
           endif

              write(ln,19) ' CALL EPSIFQED(',info_ijk(i,2,1),',',info_ijk(i,2,4),
     &          ',',chargint(i,1),',',chargint(i,2),',',info_ijk(i,1,1),',sikzone,mass_k,esq,e)'
               write(ln,11) 'epssq=epssq-esq*ANSZONE(',i,')*pdfwgt(2)'
               write(ln,11) 'eps=eps-e*ANSZONE(',i,')*pdfwgt(2)'
             write(ln,10) ' sab=massint(1)**2+massint(2)**2+two*dot(p(0,1),p(0,2))'
             write(ln,74) ' sabt=massint(1)**2+massint(2)**2+two*dot(p'
     &                      ,info_ijk(i,1,4),'(0,1),p',info_ijk(i,1,4),'(0,2))'
            if (info_ijk(i,1,1).eq.1) then
             write(ln,77) ' CALL FINITEIFQED(mass_k,sik,sikzone,z,',info_ijk(i,2,1),
     &          ',',info_ijk(i,2,4),',',chargint(i,1),',',chargint(i,2),','
     &          ,leg,',',info_p(5,inverse(iline_copy(1,i))),',sab,sabt,fi)'
            elseif(info_ijk(i,1,1).eq.2) then
             write(ln,77) ' CALL FINITEIFQED(mass_k,sik,sikzone,z,',info_ijk(i,2,1),
     &          ',',info_ijk(i,2,4),',',chargint(i,1),',',chargint(i,2),','
     &          ,leg,',',info_p(5,inverse(iline_copy(2,i))),',sab,sabt,fi)'
            endif
c
c quark-quark splitting (i:quark, ij~: quark) and
c gluon-gluon splitting (i:gluon, ij~: gluon)
c
            if(((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.1).or.
     &           (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.0))) then
               if (mass(info_ijk(i,1,3)).eq.'ZERO') then
                  write(ln,80) 'finite(1)=finite(1)-(fi(1)+fi(3))'//
     &                 '*ANS(',i,')*pdfwgt(2)*pswgt',info_ijk(i,1,4)+1,')'
                  write(ln,11) 'finite(2)=finite(2)-(fi(2)-fi(4))'//
     &                 '*ANSZONE(',i,')*pdfwgt(2)*pswgt(1)'
               write(ln,11) ' CALL WRITEHIST(P,-(fi(2)-fi(4))*ANSZONE(',i,')*pdfwgt(2)*pswgt(1))'
               write(ln,8) ' CALL WRITEHIST(P',info_ijk(i,1,4),',-(fi(1)+fi(3))/z*ANS(',i,
     &                      ')*pdfwgt(2)*pswgt(',info_ijk(i,1,4)+1,'))'
               else
                  write(ln,80) 'finite(1)=finite(1)-'//
     &                 '(fi(1)+fi(3)+fi(5))*ANS(',i,')*pdfwgt(2)*pswgt(',info_ijk(i,1,4)+1,')'
                  write(ln,11) 'finite(2)=finite(2)-'//
     &                 '(fi(2)-fi(4)-fi(6))*ANSZONE(',i,')*pdfwgt(2)*pswgt(1)'
               write(ln,11) ' CALL WRITEHIST(P,-(fi(2)-fi(4)-fi(6))*ANSZONE(',i,')*pdfwgt(2)*pswgt(1))'
               write(ln,8) ' CALL WRITEHIST(P',info_ijk(i,1,4),',-(fi(1)+fi(3)+fi(5))/z*ANS(',i,
     &                      ')*pdfwgt(2)*pswgt(',info_ijk(i,1,4)+1,'))'
               endif
c
c quark-gluon splitting (i:quark, ij~: gluon) and
c gluon-quark splitting (i:gluon, ij~: gluon)
c
            elseif((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.0) .or.
     &             (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.1)) then
             write(ln,80) 'finite(1)=finite(1)-fi(1)'//
     &              '*ANS(',i,')*pdfwgt(2)*pswgt',info_ijk(i,1,4)+1,')'
             write(ln,8) ' CALL WRITEHIST(P',info_ijk(i,1,4),',-fi(1)/z*ANS(',i,
     &                      ')*pdfwgt(2)*pswgt(',info_ijk(i,1,4)+1,'))'
            else
               write (*,*) 'Error #142 in dipolcolor'
               write (ln,*) 'Error #142 in dipolcolor'
            endif


          endif


c**************************************************************
c     initial-initial theta function
c**************************************************************
         elseif(info_ijk(i,1,1).le.nincoming .and.
     &           info_ijk(i,1,3) .le.nincoming) then
            if(info_ijk(i,2,1).eq.0) then
              leg=info_ijk(i,1,2)
            elseif(info_ijk(i,2,2).eq.0) then
              leg=info_ijk(i,1,1)
            elseif(info_ijk(i,2,1).eq.info_ijk(i,2,2)) then
              leg=info_ijk(i,1,1)
            endif
            write(ln,10)' '
            write(ln,74)'sikzone=two*dot(p(0,',info_ijk(i,1,4),
     &           '),p(0,',info_ijk(i,1,5),'))'


          if(.not.lepton.and..not.photon.or.(photon.and.photonpdf)) then


            write(ln,74)'sik=two*dot(p(0,',info_ijk(i,1,4),
     &           '),p(0,',info_ijk(i,1,5),'))/z'
            if (info_ijk(i,1,1).eq.1) then
               write(ln,17)'CALL DIPOLEPDF(x(1)/z,x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(1))'
            elseif(info_ijk(i,1,1).eq.2) then
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2)/z,',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(1))'
            else
               write (ln,*)'Error #147 in write_intdipole'
               write (*,*)'Error #147 in write_intdipole'
            endif
            if(((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.1).or.
     &         (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.0))) then
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(2))'
            endif
            if(info_ijk(i,2,1).eq.info_ijk(i,2,2).and..not.
     &          first(info_ijk(i,1,1),info_ijk(i,1,2))) then
             if (i.le.9) then
                write(ln,51)'CALL SDIPOLQED'//char(48)//char(48)
     &            //char(48+i)//'(P,ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,51)'CALL SDIPOLQED'//char(48)//char(i/10+48)
     &            //char(mod(i,10)+48)//'(P,ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,51)'CALL SDIPOLQED'//char(i/100+48)
     &            //char((i-(i/100)*100)/10+48)
     &            //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
           else
             if (i.le.9) then
                write(ln,41)'CALL SDIPOLQED'//char(48)//char(48)
     &            //char(48+i)//'(P,ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,41)'CALL SDIPOLQED'//char(48)//char(i/10+48)
     &            //char(mod(i,10)+48)//'(P,ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,41)'CALL SDIPOLQED'//char(i/100+48)
     &            //char((i-(i/100)*100)/10+48)
     &            //char(mod(i-(i/100)*100,10)+48)//'(P,ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             first(info_ijk(i,1,1),info_ijk(i,1,2))=.false.
            endif

             write(ln,19) ' CALL EPSIIQED(',info_ijk(i,2,1),',',info_ijk(i,2,4),',',
     &             chargint(i,1),',',chargint(i,2),',',info_ijk(i,1,1),',sikzone,esq,e)'
               write(ln,11) 'epssq=epssq-esq*ANS(',i,')*pdfwgt(2)'
               write(ln,11) 'eps=eps-e*ANS(',i,')*pdfwgt(2)'
            if (info_ijk(i,1,1).eq.1) then
             write(ln,77) ' CALL FINITEIIQED(sik,sikzone,z,',info_ijk(i,2,1),','
     &            ,info_ijk(i,2,4),',',chargint(i,1),',',chargint(i,2),','
     &            ,leg,',',info_p(5,inverse(iline_copy(1,i))),',fi)'
            elseif(info_ijk(i,1,1).eq.2) then
             write(ln,77) ' CALL FINITEIIQED(sik,sikzone,z,',info_ijk(i,2,1),','
     &            ,info_ijk(i,2,4),',',chargint(i,1),',',chargint(i,2),','
     &            ,leg,',',info_p(5,inverse(iline_copy(2,i))),',fi)'
            endif
c
c quark-quark splitting (i:quark, ij~: quark) and
c gluon-gluon splitting (i:gluon, ij~: gluon)
c
            if(((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.1).or.
     &         (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.0))) then
               if (mass(info_ijk(i,1,3)).eq.'ZERO') then
                  write(ln,11) 'finite(1)=finite(1)-(fi(1)+fi(3))'//
     &                 '*ANS(',i,')*pdfwgt(1)'
                  write(ln,11) 'finite(2)=finite(2)-(fi(2)-fi(4))'//
     &                 '*ANS(',i,')*pdfwgt(2)'
               else
                  write(ln,11) 'finite(1)=finite(1)-'//
     &                 '(fi(1)+fi(3)+fi(5))*ANS(',i,')*pdfwgt(1)'
                  write(ln,11) 'finite(2)=finite(2)-'//
     &                 '(fi(2)-fi(4)-fi(6))*ANS(',i,')*pdfwgt(2)'
               endif
c
c quark-gluon splitting (i:quark, ij~: gluon) and
c gluon-quark splitting (i:gluon, ij~: gluon)
c
            elseif((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.0) .or.
     &             (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.1)) then
               write(ln,11)
     &              'finite(1)=finite(1)-fi(1)*ANS(',i,')*pdfwgt(1)'
            else
               write (*,*) 'Error #144 in dipolcolor'
               write (ln,*) 'Error #144 in dipolcolor'
            endif



         elseif(lepton.or.(photon.and..not.photonpdf)) then


            write(ln,75)'sik=two*dot(p',info_ijk(i,1,4),'(0,',info_ijk(i,1,4)
     &           ,'),p',info_ijk(i,1,4),'(0,',info_ijk(i,1,5),'))'
               write(ln,17)'CALL DIPOLEPDF(x(1),x(2),',
     &              info_p(5,inverse(iline_copy(1,i))),',',
     &              info_p(5,inverse(iline_copy(2,i))),',pdfwgt(2))'
           if(info_ijk(i,2,1).eq.info_ijk(i,2,2).and..not.
     &          first(info_ijk(i,1,1),info_ijk(i,1,2))) then
             write(ln,15)'IF (PASSCUTSDIP(P)) THEN'
             if (i.le.9) then
                write(ln,61)' CALL SDIPOLQED'//char(48)//char(48)
     &            //char(48+i)//'(P,ANSZONE(',i,'))'
             elseif(i.le.99) then
                write(ln,61)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &            //char(mod(i,10)+48)//'(P,ANSZONE(',i,'))'
             elseif(i.le.999) then
                write(ln,61)' CALL SDIPOLQED'//char(i/100+48)
     &            //char((i-(i/100)*100)/10+48)
     &            //char(mod(i-(i/100)*100,10)+48)//'(P,ANSZONE(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             write(ln,15)'ENDIF'
c             write(ln,6)' CALL BOOST(P,Z,',info_ijk(i,1,1),',P1)'
             write(ln,9)'IF (PASSCUTSDIP(P',info_ijk(i,1,4),')) THEN'
             if (i.le.9) then
                write(ln,51)' CALL SDIPOLQED'//char(48)//char(48)
     &            //char(48+i)//'(P',info_ijk(i,1,4),',ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,51)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &            //char(mod(i,10)+48)//'(P',info_ijk(i,1,4),',ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,51)' CALL SDIPOLQED'//char(i/100+48)
     &            //char((i-(i/100)*100)/10+48)
     &            //char(mod(i-(i/100)*100,10)+48)//'(P',info_ijk(i,1,4),',ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             write(ln,15)'ENDIF'
           else
             write(ln,10)'IF (PASSCUTSDIP(P)) THEN'
             if (i.le.9) then
                write(ln,41)' CALL SDIPOLQED'//char(48)//char(48)//char(48+i)
     &               //'(P,ANSZONE(',i,'))'
             elseif(i.le.99) then
                write(ln,41)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &               //char(mod(i,10)+48)//'(P,ANSZONE(',i,'))'
             elseif(i.le.999) then
                write(ln,41)' CALL SDIPOLQED'//char(i/100+48)
     &               //char((i-(i/100)*100)/10+48)
     &               //char(mod(i-(i/100)*100,10)+48)//'(P,ANSZONE(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             write(ln,10) 'ENDIF'
c             write(ln,13)' CALL BOOST(P,Z,',info_ijk(i,1,1),',P1)'
             write(ln,7)'IF (PASSCUTSDIP(P',info_ijk(i,1,4),')) THEN'
             if (i.le.9) then
                write(ln,41)' CALL SDIPOLQED'//char(48)//char(48)//char(48+i)
     &               //'(P',info_ijk(i,1,4),',ANS(',i,'))'
             elseif(i.le.99) then
                write(ln,41)' CALL SDIPOLQED'//char(48)//char(i/10+48)
     &               //char(mod(i,10)+48)//'(P',info_ijk(i,1,4),',ANS(',i,'))'
             elseif(i.le.999) then
                write(ln,41)' CALL SDIPOLQED'//char(i/100+48)
     &               //char((i-(i/100)*100)/10+48)
     &               //char(mod(i-(i/100)*100,10)+48)//'(P',info_ijk(i,1,4),',ANS(',i,'))'
             else
                write (*,*) 'ERROR: write_dipoles: too many dipoles'
             endif
             write(ln,10) 'ENDIF' 
             first(info_ijk(i,1,1),info_ijk(i,1,2))=.false.
           endif

             write(ln,19) ' CALL EPSIIQED(',info_ijk(i,2,1),',',info_ijk(i,2,4),',',
     &             chargint(i,1),',',chargint(i,2),',',info_ijk(i,1,1),',sikzone,esq,e)'
               write(ln,11) 'epssq=epssq-esq*ANSZONE(',i,')*pdfwgt(2)'
               write(ln,11) 'eps=eps-e*ANSZONE(',i,')*pdfwgt(2)'
            if (info_ijk(i,1,1).eq.1) then
             write(ln,77) ' CALL FINITEIIQED(sik,sikzone,z,',info_ijk(i,2,1),','
     &            ,info_ijk(i,2,4),',',chargint(i,1),',',chargint(i,2),','
     &            ,leg,',',info_p(5,inverse(iline_copy(1,i))),',fi)'
            elseif(info_ijk(i,1,1).eq.2) then
             write(ln,77) ' CALL FINITEIIQED(sik,sikzone,z,',info_ijk(i,2,1),','
     &            ,info_ijk(i,2,4),',',chargint(i,1),',',chargint(i,2),','
     &            ,leg,',',info_p(5,inverse(iline_copy(2,i))),',fi)'
            endif
c
c quark-quark splitting (i:quark, ij~: quark) and
c gluon-gluon splitting (i:gluon, ij~: gluon)
c
            if(((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.1).or.
     &         (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.0))) then
               if (mass(info_ijk(i,1,3)).eq.'ZERO') then
                  write(ln,80) 'finite(1)=finite(1)-(fi(1)+fi(3))'//
     &                 '*ANS(',i,')*pdfwgt(2)*pswgt(',info_ijk(i,1,4)+1,')'
                  write(ln,11) 'finite(2)=finite(2)-(fi(2)-fi(4))'//
     &                 '*ANSZONE(',i,')*pdfwgt(2)*pswgt(1)'
                  write(ln,11) ' CALL WRITEHIST(P,-(fi(2)-fi(4))*ANSZONE(',i,')*pdfwgt(2)*pswgt(1))'
                  write(ln,8) ' CALL WRITEHIST(P',info_ijk(i,1,4),',-(fi(1)+fi(3))/z*ANS(',i,
     &                      ')*pdfwgt(2)*pswgt(',info_ijk(i,1,4)+1,'))'
               else
                  write(ln,80) 'finite(1)=finite(1)-'//
     &                 '(fi(1)+fi(3)+fi(5))*ANS(',i,')*pdfwgt(2)*pswgt(',info_ijk(i,1,4)+1,')'
                  write(ln,11) 'finite(2)=finite(2)-'//
     &                 '(fi(2)-fi(4)-fi(6))*ANSZONE(',i,')*pdfwgt(2)*pswgt(1)'
                  write(ln,11) ' CALL WRITEHIST(P,-(fi(2)-fi(4)-fi(6))*ANSZONE(',i,')*pdfwgt(2)*pswgt(1))'
                  write(ln,8) ' CALL WRITEHIST(P',info_ijk(i,1,4),',-(fi(1)+fi(3)+fi(5))/z*ANS(',i,
     &                      ')*pdfwgt(2)*pswgt(',info_ijk(i,1,4)+1,'))'
               endif
c
c quark-gluon splitting (i:quark, ij~: gluon) and
c gluon-quark splitting (i:gluon, ij~: gluon)
c
            elseif((info_ijk(i,2,1).eq.1.and.info_ijk(i,2,4).eq.0) .or.
     &             (info_ijk(i,2,1).eq.0.and.info_ijk(i,2,4).eq.1)) then
               write(ln,80)
     &              'finite(1)=finite(1)-fi(1)*ANS(',i,')*pdfwgt(2)*pswgt(',info_ijk(i,1,4)+1,')'
               write(ln,8) ' CALL WRITEHIST(P',info_ijk(i,1,4),',-fi(1)/z*ANS(',i,
     &                      ')*pdfwgt(2)*pswgt(',info_ijk(i,1,4)+1,'))'
            else
               write (*,*) 'Error #144 in dipolcolor'
               write (ln,*) 'Error #144 in dipolcolor'
            endif





         endif


         else
            write (*,*) 'Error #143 in dipolcolor'
            write (ln,*) 'Error #143 in dipolcolor'
         endif
         endif ! massless unresolved
      enddo
      write(ln,10) ' '
      if(.not.lepton.and..not.photon.or.(photon.and.photonpdf)) then
       write(ln,10) 'ENDIF' 
      endif
      write(ln,10) 'INTDIPOLE=.FALSE.'
      if(lepton.or.(photon.and..not.photonpdf)) then
       write(ln,10) 'FIN=(FINITE(1)/z+FINITE(2)+FINITE(3))'
      else
       write(ln,10) 'FIN=(FINITE(1)/z+FINITE(2)+FINITE(3))*pswgt(1)'
       write(ln,10) 'CALL WRITEHIST(P,FIN)'
      endif
      write(ln,10) 'RETURN'
      write(ln,10) 'END'
      write(ln,10) ' '
      write(ln,10) ' '



c   Cut Routine for non-collinear-safe observables

      write(ln,10) 'LOGICAL FUNCTION THETACUT(P,ZI)'
      write(ln,10) 'IMPLICIT NONE'
      write(ln,5) 'This is a template for a Theta function to include'
      write(ln,5) 'correction terms for non-collinear-safe observables.'
      write(ln,10) ' '
      write(ln,10) 'include "nexternal.inc"'
      write(ln,10) 'REAL*8 P(0:3,nexternal),zi'
      write(ln,10) ' '
      write(ln,5) 'This is an empty routine so return TRUE'
      write(ln,10) 'THETACUT=.TRUE.'
      write(ln,10) ' '
      write(ln,10) 'END'

c   Routine for calling nonperturbative parts of fragmentation function

      write(ln,10) ' '
      write(ln,10) ' '
      write(ln,10) 'REAL*8 FUNCTION DFRAG(ZI)'
      write(ln,10) 'IMPLICIT NONE'
      write(ln,5) 'This is a template for calculating the non-perturbative'
      write(ln,5) 'parts of the quark-to-photon fragmentation function'
      write(ln,10) ' '
      write(ln,10) 'include "nexternal.inc"'
      write(ln,10) 'REAL*8 zi'
      write(ln,10) ' '
      write(ln,5)  'This is an empty routine so return zero'
      write(ln,10) 'DFRAG=0d0'
      write(ln,10) ' '
      write(ln,10) 'END'


 5    format('C ',A)
 6    format('C ',7X,A,I2,A)
 7    format(6X,A,I1,A)
 8    format(6X,A,I1,A,I4,A,I1,A)
 9    format('C',6X,A,I1,A)
 10   format(6X,A)
 11   format(6X,A,I3,A)
 17   format(6X,A,I3,A,I3,A)
 12   format(5X,'&',A,I3,A)
 13   format(6X,A,I2,A)
 14   format(6X,A,I3,A,I3,A)
 15   format('C',6X,A)
 16   format(6X,A,A,A1,A4,A)
 18   format(6X,A,4(I3,A))
 19   format(6X,A,5(I3,A))
 20   format(', P',I1)
 30   format(',P',I1,'(0:3)')
 40   format(6X,A21,I4,A2)
 41   format(6X,A,I1,A,I4,A2)
 42   format(6X,A23,I4,A2)
 50   format(7X,A26,I4,A1)
 51   format('C ',6X,A,I1,A,I4,A2)
 60   format(6X,A,I4,A1)
 61   format('C ',6X,A,I4,A2)
 62   format(6X,A,I4,A2)
 70   format(6X,A12,I2,A6,I2,A2)
 71   format(6X,A,A5)
 72   format(6X,A15,2(I1,A7),I1,A)
 73   format(6X,A15,I1,A7,I1,A)
 74   format(6X,2(A,I1),A)
 75   format(6X,4(A,I1),A)
 76   format(6X,3(A,I3),A)
 77   format(6X,6(A,I3),A)
 78   format(6X,3(A,I1),A)
 80   format(6X,A,I4,A,I1,A)
      close(ln)







      end


      Subroutine write_dipoleinc()
      implicit none

      include 'nexternal.inc'
      logical photonpdf
      common/photon/ photonpdf
      integer ln,ln1,i,j,k
      logical done, pdffound, massfound
      character*120 buff

      done=.false.
      pdffound=.false.
      massfound=.false.
      ln=31
      ln1=32
      open(unit=ln,file='../SubProcesses/dipole.inc',status='unknown')
      open(unit=ln1,file='../SubProcesses/dipole.inc.tmp',status='unknown')
      rewind(ln)
      rewind(ln1)

      do while (.not.done)
 10     read(ln,'(a)',end=11),buff

       do i=1,120

        if(.not.pdffound) then
        if((index(buff(i:i+10),"photonpdf=").gt.0).and.photonpdf) then
         buff="      parameter (photonpdf=.TRUE.)"
         pdffound=.true.
         write(ln1,'(a)') buff
         goto 10
        elseif((index(buff(i:i+10),"photonpdf=").gt.0).and..not.photonpdf) then
         buff="      parameter (photonpdf=.FALSE.)"
         pdffound=.true.
         write(ln1,'(a)') buff
         goto 10
        endif
        endif


        if(.not.massfound) then
        if((index(buff(i:i+8),"massint(").gt.0)) then
          buff(1:22)="      real*8 massint("
          buff(23:24)=char(48+nexternal)
          buff(24:25)=")"
          write(ln1,'(a)') buff
          goto 10
        endif

         if(index(buff(i:i+20),"      data massint/").gt.0) then
          buff(1:20)="      data massint/"
          j=21
           do k=1, nexternal-1
             buff(j:j+4)="1d0,"
             j=j+4
           enddo
           buff(j:j+4)="1d0/"
          buff=buff(1:j+4)
          write(ln1,'(a)') buff
          massfound=.true.
          goto 10
        endif
        endif 

       enddo
c        print*, buff
        write(ln1,'(a)') buff
      enddo

 11   continue
      close(ln)
      close(ln1)
      end

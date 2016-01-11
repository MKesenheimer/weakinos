      program plot_events
c***********************************************************************
c     1. reads banner and finds the particles classes                  *
c     2. reads events in a applies further cuts if needed              *
c     3. write out plots.top                                           *    
c----------------------------------------------------------------------*
c     FIRST  VERSION 16-May-2003  BY FM                                *
c     SECOND VERSION 25-Aug-2005  BY FM                                *
c----------------------------------------------------------------------*
c     LAST UPDATE  23-FEB-2006 BY RIKKERT FREDERIX                     *
c     many improvements:                                               *
c     - events can be different                                        *
c     - inclusive plots added                                          *
c     - pgs events can be added                                        *
c     - sa_card.dat for input classes, cuts and plot options           *
c***********************************************************************
      implicit none

      include 'info.inc'
c     
c     Local
c
      integer nexternal,nfinal,ic(7,MaxParticles),ictemp(MaxParticles)
      double precision P(0:4,MaxParticles),wgt,pmiss(0:3)
      real*8 sum,mxwgt,mET
      real*8 sump,mxwgtp
      integer i,j,l,k,m,lmax,nevent,neventp
      integer ievent,checkwght
      double precision scale,aqcd,aqed,sumw,sumpw
      logical done,found,firsttime
      data firsttime/.true./
      character*130 infile,outfile,logfile
      data outfile/'plots.top'/
      data logfile/'plots.log'/
      integer npar,idn(maxparticles),linfile
      double precision ptv(maxparticles)
      character*140 buff2
c
c     Global
c
      integer maxclass,maxevparticles,partoclass(-300:300,0:Maxparticles)
      integer maxevpart(0:maxclasses)
      common/classes/maxclass,maxevparticles,partoclass,maxevpart
      logical lhe,unweighted
      common/fileformat/lhe,unweighted
      double precision awgt
      common/ne/sum,awgt,nevent

c
c     External
c
      logical passcuts,decayed
      double precision pt,et,ordering
      common/plot_decayed/decayed
c
c
c-----
c   Begin Code
c-----

      write(*,*) '*****************************************************'
      write(*,*) '*                 SimpleAnalysis                    *'
      write(*,*) '*               a MadEvent program                  *'
      write(*,*) '*    for making plots and do simple analysis for    *'
      write(*,*) '*            LHE and LHCO event files               *'  
      write(*,*) '*        ---------------------------------          *'
      write(*,*) '*        Version compliant with MG_ME_V4.0          *'
      write(*,*) '*                                                   *'
      write(*,*) '*       Authors: R. Frederix and F. Maltoni         *'
      write(*,*) '*                                                   *'
      write(*,*) '*         Last Update by R.F., 23-Feb-2006          *'
      write(*,*) '*****************************************************'
      write(*,*) ' '


      write(*,*) 'input name of event file:'
      write(*,*) '-------------------------'

      read(*,'(a)')  infile
      open(unit=98,file=infile,status='old')


      open(unit=52,file=logfile,status='unknown')

      call no_spaces_long(infile,linfile)
      if (infile(linfile-3:linfile).eq.'.lhe')then
         LHE=.true.
         write (*,'(a)') 'Found Les Houches Event File '
      elseif (infile(linfile-3:linfile).eq.'lhco')then
         LHE=.false.
         write (*,'(a)') 'Found LHC Olympics (PGS4) Event File '
      else
         write (*,*) 'File extension unknown'
         write (*,*) 'Can only plot .lhe and .lhco files.'
         return
      endif

      outfile='plots.top'
      write(*,*) '   output file: ',outfile
      open(unit=99,file=outfile,status='unknown')
c
c--   Create all classes :events will be read a first time
c
      call get_classes

c--   Set-up histograms

      call histo_init
c
c--   Initialize cuts
c
      call setcuts
c
c now start reading the events to fill histograms
c

      write (*,*) "Now parsing events to fill histograms"
      done=.false.
      nevent=0      
      sumw=0d0
      sum=0d0
      mxwgt=-1d0
      neventp=0
      sump=0d0
      mxwgtp=-1d0
      checkwght=1
      rewind(98)
      do while (.not. done)
         if (lhe)then
            call read_event(98,P,wgt,nexternal,ic,ievent,scale,aqcd,aqed,buff2,done)
         else
            call read_lhco_event(98,P,wgt,nexternal,ic,done,sum)
         endif
	 if(.not.done) then
            nevent=nevent+1
            sumw=sumw+wgt
            if (lhe) then
               sum=sumw
            endif
            mxwgt = max(wgt,mxwgt)
c     
c     remove decayed resonances from particle listing if .not.decayed and parton level events
c     
            nfinal=nexternal
            if (.not.decayed.and.lhe)then
               do i=nexternal-2,3,-1
                  if(ic(6,i).ne.1) then
                     do j=i+1,nfinal
                        do k=1,7
                           ic(k,j-1)=ic(k,j)
                        enddo
                        do k=0,4
                           P(k,j-1)=P(k,j)
                        enddo
                     enddo
                     nfinal=nfinal-1
                  endif
               enddo
            endif
c For SUSY particles
            do i=3,nfinal
               ictemp(i)=MOD(ic(1,i),1000000)+ic(1,i)/1000000*100
            enddo
c
c     put particle class information in eventinfo vector
c
            do i=0,maxclasses
               do j=1, maxparticles
                  if (eventinfo(i,j).ne.0) eventinfo(i,j)=0
               enddo
            enddo
            do i=3,nfinal
              if (partoclass(ictemp(i),0).ge.1)then
                  do j=1,partoclass(ictemp(i),0)
                     k=1
                     do while(eventinfo(partoclass(ictemp(i),j),k).ne.0)
                        k=k+1
                     enddo
                     eventinfo(partoclass(ictemp(i),j),k)=i  
                  enddo
               endif
            enddo
c
c     Order particles in vector per class -> Uses the ordering function to order
c
            do i=0,maxclass         
               npar=1
               do while (npar.le.maxparticles.and.eventinfo(i,npar).ne.0)
                  ptv(npar)=ordering(p(0,eventinfo(i,npar)))
                  idn(npar)=eventinfo(i,npar)
                  npar=npar+1
               enddo
               if (npar-1.gt.1)then
                  call sort2(ptv,idn,npar-1)
                  do k=1,npar-1
                     eventinfo(i,k)=idn(k)
                  enddo
               endif
            enddo
c
c For the missing ET class (mET) we have to calculated the missing ET
c
            if (maxevpart(0).ge.1) then
               do j=0,3
                  pmiss(j)=0d0
               enddo
               do k=1,maxevpart(0)
                  do j=0,3
                     pmiss(j)=pmiss(j)+p(j,eventinfo(0,k))
                  enddo
               enddo
               mET=et(pmiss) !Missing Et
            else
               mET=0d0
            endif

c            do i=0,maxclass
c            write (*,*)'class:',i,', with particles:',(eventinfo(i,k),k=1,maxparticles)
c            enddo

            if(passcuts(p,nfinal)) then ! check cuts 
               if (unweighted)then
                  call graph_event(p,1d0,mET)
               else
                  call graph_event(p,wgt,mET)
               endif
               mxwgtp = max(wgt,mxwgtp)
               neventp=neventp+1
               if (lhe)then
                  sump=sump+wgt
                  sumpw=sump
                  awgt=1
                  if (.not.unweighted) awgt=sump/neventp
               else
                  sump=neventp/nevent*sum
                  sumpw=neventp
                  awgt=1d0
               endif
            endif
         endif
      enddo
c      firsttime=.false.
c      enddo

      write(*,*) ' '
      write(*,*) '----------------------------------------------'
      write(*,*) '              events information              '
      write(*,*) '----------------------------------------------'
      write(*,*) ' '
      write(*,*) '  Events in file               :  ',nevent
      write(*,*) '  Integrated weight (pb)       :  ',sum
      write(*,*) '  Max wgt                      :  ',mxwgt
      write(*,*) '  Average wgt                  :  ',sumw/nevent
      write(*,*) ' '
      write(*,*) '  Passing the cuts (plotted)'
      write(*,*) '  --------------------------'
      write(*,*) ' '
      write(*,*) '  Events                       :  ',neventp
      write(*,*) '  Integrated weight (pb)       :  ',sump
      write(*,*) '  Max wgt                      :  ',mxwgtp
      write(*,*) '  Average wgt                  :  ',sumpw/neventp
c
c     finalize histograms
c
      call histo_final
c
c     close datafiles
c
      close(52)
      close(97)
      close(98)
      close(99)
 99   continue
      end




      subroutine get_classes
c###################################################################
c Get all classes and max number of particles
c and matrix with all particle info
c###################################################################
      implicit none
      include 'info.inc'

      integer nexternal,nfinal,ic(7,MaxParticles),ictemp(MaxParticles)
      double precision P(0:4,MaxParticles),wgt,sum,wgt1
      integer ievent
      double precision scale,aqcd,aqed
      logical done
      data done/.false./

      character*25 classesfile
      character*4 classname(0:maxclasses),buffclass
      character*15 plotname
      character*7 test,output
      character*140 buff2,buffer
      integer iname(maxparticles),l1
      logical foundallclasses, foundallparticles,firsttime
      integer nclass,nparticles,class(0:maxclasses,maxparticles)
      integer maxparticlesinclass(0:maxclasses),evpartinclass(0:maxclasses)
      integer partoclass(-300:300,0:maxparticles),maxevpart(0:maxclasses)
      integer maxevparticles,maxclass,avpart
      integer beginclasses,endclasses,beginoptions,endoptions
      integer plot,plotclass
      integer i,j,k,l,m
      logical lhe,unweighted
      common/classes/maxclass,maxevparticles,partoclass,maxevpart
      common/fileformat/lhe,unweighted
      logical decayed
      common/plot_decayed/decayed
      character*4 id(-300:300)
      common/to_id/id
      data classesfile/'sa_card.dat'/


      write (*,*) "Classes will be read from file: ./sa_card.dat"
      open(unit=97,file=classesfile,status='old')

      call set_id




      rewind(97)
 61   read(97,'(a)',end=90,err=90)buffer
      if(buffer(1:18).eq. "# Begin PlotOutput ")then
 60      read(97,'(a)')buffer
            if (buffer(1:16).eq. "# End PlotOutput ") goto 63 ! End of PlotOptions
            l1=40
            if(index(buffer,"#").ne.0) l1=index(buffer,"#")-1 ! ignore comments
            read(buffer(1:l1),*,end=62,err=62) plotname,output
 62         continue
            if (plotname(1:6).eq.'output') then
               if (output(1:7).eq.'gnuplot')then
                  gnu=.true.
               else
                  gnu=.false.
               endif
            elseif (plotname(1:12).eq.'plot_decayed') then
               if (output(1:3).eq.'yes')then
                  decayed=.true.
               else
                  decayed=.false.
               endif
            endif
         goto 60  !Not yet at end of PlotOptions -> Read next line
 63      continue !Found all PlotOptions ->  continue here
      endif     
      goto 61   !Go back to beginning to read next line
 90   continue  !End of file -> continue






      do i=0,maxclasses
         do j=1,maxparticles
            class(i,j)=0
         enddo
      enddo
      nclass=0
      m=0
      rewind(97)
 10   read (97,'(a)',end=76,err=76) buffer
      if(buffer(1:16).eq. "# Begin Classes ")then
 43      read (97,'(a)') buffer
            if (buffer(1:14).eq. "# End Classes ") goto 44 ! End of classes list
            l1=40
            if(index(buffer,"#").ne.0) l1=index(buffer,"#")-1 ! ignore comments
            nclass=nclass+1
            read(buffer(1:l1),*,end=11,err=11) classname(nclass),(class(nclass,i),i=1,maxparticles)
 11         continue
ccc   Missing ET class (mET) is special. Also remove empty classes from the list
            if (classname(nclass).eq.'mET')then !missing ET class is special
               classname(0)='mET'
               do i=1,maxparticles
                  class(0,i)=class(nclass,i)
               enddo
               if (class(0,1).eq.0) then !empty missing ET class
                  maxparticlesinclass(0)=0
               else
                  l=1
                  do while (class(0,l).ne.0) !find all particles in missing ET class
                     maxparticlesinclass(0)=l
                     l=l+1
                  enddo
               endif
               nclass=nclass-1
            elseif (class(nclass,1).eq.0 !empty class, remove from list
     &              .and.classname(nclass).ne.'mET')then !if it is not the mET class
               nclass=nclass-1
            endif
            l=1
            do while (class(nclass,l).ne.0) !find all particles in class
               maxparticlesinclass(nclass)=l
               l=l+1
            enddo
         goto 43   ! Go back to read class from list
 44      continue  ! Found all classes in list -> continue here
      endif
      goto 10      ! Go back to beginning to read next line of the file 
 76   continue     ! Continue here when at end of file
c Make sure that the missing ET class exists
      if (classname(0).ne.'mET')then
         classname(0)='mET'
         maxparticlesinclass(0)=0
         if (.not.lhe) then     ! for the .lhco file add the 'neutrino' to this class
            maxparticlesinclass(0)=1
            class(0,1)=12
         endif
      endif
c      write (*,*) nclass

c      do i=1,nclass
c         write (*,*) classname(i),(class(i,j),j=1,maxparticlesinclass(i))
c      enddo
      
 45   foundallclasses=.true.

      write (*,*) 'Found the following classes in input classes file:'
      if (nclass.ne.0)then
         do i=0,nclass
            write (*,'(2x,a7,i2,a3,3x,a4,30(i4))')
     &           'class #',i,' is',classname(i),(class(i,j),j=1,maxparticlesinclass(i))
         enddo
      else
         write (*,*)'None, no classes specified in ',classesfile
      endif

c
c read all events to fill the classes
c
      write (*,*) "Now start parsing events a first time"
      do i=0,maxclasses
         maxevpart(i)=0
      enddo
      firsttime=.true.
      rewind(98)
      do while (.not.done)
         if (lhe)then
            call read_event(98,P,wgt,nexternal,ic,ievent,scale,aqcd,aqed,buff2,done)
         else
            call read_lhco_event(98,P,wgt,nexternal,ic,done,sum)
         endif
         if (done) goto 65
c
c     Check if this are weighted or unweighted events
c
         if (firsttime)then
            wgt1=wgt
            unweighted=.true.
            firsttime=.false.
         else
            if (wgt.eq.wgt1.and.unweighted)then
               unweighted=.true.
            elseif(wgt.ne.wgt1)then
               unweighted=.false.
            endif
         endif
c
c     remove decayed resonances from particle listing if .not.decayed and parton level events
c
         nfinal=nexternal
         if (.not.decayed.and.lhe)then
            do i=nexternal-2,3,-1
               if(ic(6,i).ne.1) then
                  do j=i+1,nfinal
                     do k=1,7
                        ic(k,j-1)=ic(k,j)
                     enddo
                     do k=0,4
                        P(k,j-1)=P(k,j)
                     enddo
                  enddo
                  nfinal=nfinal-1
               endif
            enddo
         endif
c For Susy particles:
         do i=3,nfinal
            ictemp(i)=MOD(ic(1,i),1000000)+ic(1,i)/1000000*100
         enddo
c
c Compare the particles in the events with the classes
c
         do i=0,nclass
            evpartinclass(i)=0
         enddo
         do i=3,nfinal
            k=0
            do j=0,nclass
               do l=1,maxparticlesinclass(j)
c                  write (*,*) ictemp(i),class(j,l)
                  if (ictemp(i).eq.class(j,l)) then
                     k=k+1
                     partoclass(ictemp(i),0) = k
                     partoclass(ictemp(i),k) = j
                     evpartinclass(j) = evpartinclass(j)+1
c                     write (*,*)ictemp(i),partoclass(ictemp(i),1)
                     maxevpart(j)=max(maxevpart(j),evpartinclass(j))
c                     if (maxevpart(j).ge.)write (*,*) j,maxevpart(j)
                  endif
               enddo
            enddo
            if (k.eq.0) then !particle does not fit in existing class -> create new class
               nclass=nclass+1
               if (nclass.gt.maxclasses)then
                  write (*,*) 'Error! Too many classes, max is',MaxClasses
               endif
               class(nclass,1)=ictemp(i)
               classname(nclass)=id(ictemp(i))
               partoclass(ictemp(i),0) = 1
               partoclass(ictemp(i),1) = nclass
               maxparticlesinclass(nclass)=1
               maxevpart(nclass)=1
        write (*,*)'created new class ',classname(nclass),'for particle',ictemp(i)
            endif
         enddo
      enddo
 65   continue
      
      maxclass=nclass

      maxevparticles=0
      avpart=0
      do i=0,nclass
         maxevparticles=maxevparticles+maxevpart(i)
         avpart=max(maxevpart(i),avpart)
c      write (*,*) i,maxevpart(i)
      enddo
      if (maxevparticles.gt.maxparticles)then
         write (*,*) 'Error! Too many particles, max is',MaxParticles
         write (*,*) maxevparticles
      else
      write(52,*)'Final states contain',maxevparticles,
     &                                ' different particles.'
      endif

c
c     Give all classes and particles appropriate names
c
      k=1
      do i=0,maxclasses
         call no_spaces(classname(i),iname(i))
         do j=1,maxevpart(i)
            lname(i,j)=iname(i)
            name(i,j)(1:lname(i,j))=classname(i)
            if(maxevpart(i).ge.1.and.j.le.9)then
               lname(i,j)=lname(i,j)+1
               name(i,j)(lname(i,j):lname(i,j))=char(48+j)
            elseif (j.gt.9)then
               lname(i,j)=lname(i,j)+1
               name(i,j)(lname(i,j):lname(i,j))=char(49)
               lname(i,j)=lname(i,j)+1
               name(i,j)(lname(i,j):lname(i,j))=char(38+j)
            endif
      write (52,'(a10,x,i2,a2,a6,3x,a15,i2)')'particle #',k,': ',name(i,j)(1:lname(i,j)),
     &           'fits in class #',i
            k=k+1
         enddo
      enddo


c
c     Read plotting options from sa_card.dat
c
      plot=0
      plotclass=0
      plotname=''
      do i=1,maxclasses
         ptplot(i)=0
         etplot(i)=0
         eplot(i)=0
         momplot(i)=0
         etaplot(i)=0
         yplot(i)=0
         costhplot(i)=0
         phiplot(i)=0
         X1plot(i)=0
         X2plot(i)=0
         X3plot(i)=0
         dRijplot(i)=0
         mijplot(i)=0
         delta_phiplot(i)=0
         delta_etaplot(i)=0
         kTijplot(i)=0
         XY1plot(i)=0
         XY2plot(i)=0
         XY3plot(i)=0
         XYZ1plot(i)=0
         XYZ2plot(i)=0
         XYZ3plot(i)=0
         XYZA1plot(i)=0
         XYZA2plot(i)=0
         XYZA3plot(i)=0
      enddo
      rewind(97)
 51   read(97,'(a)',end=80,err=80)buffer
      if(buffer(1:17).eq. "# Begin PlotDefs ")then
 50      read(97,'(a)')buffer
            if (buffer(1:15).eq. "# End PlotDefs ") goto 53 ! End of PlotOptions
            l1=40
            if(index(buffer,"#").ne.0) l1=index(buffer,"#")-1 ! ignore comments
            read(buffer(1:l1),*,end=52,err=52) plotname,plotclass,plot
 52         continue
            if (plotname(1:3).eq.'pt ')then
               ptplot(plotclass)=plot
            elseif (plotname(1:3).eq.'et ')then
               etplot(plotclass)=plot
            elseif (plotname(1:2).eq.'e ')then
               eplot(plotclass)=plot
            elseif (plotname(1:4).eq.'mom ')then
               momplot(plotclass)=plot
            elseif (plotname(1:4).eq.'eta ')then
               etaplot(plotclass)=plot
            elseif (plotname(1:2).eq.'y ')then
               yplot(plotclass)=plot
            elseif (plotname(1:6).eq.'costh ')then
               costhplot(plotclass)=plot
            elseif (plotname(1:4).eq.'phi ')then
               phiplot(plotclass)=plot
            elseif (plotname(1:3).eq.'X1 ')then
               X1plot(plotclass)=plot
            elseif (plotname(1:3).eq.'X2 ')then
               X2plot(plotclass)=plot
            elseif (plotname(1:3).eq.'X3 ')then
               X3plot(plotclass)=plot
            elseif (plotname(1:5).eq.'dRij ')then
               dRijplot(plotclass)=plot
            elseif (plotname(1:4).eq.'mij ')then
               mijplot(plotclass)=plot
            elseif (plotname(1:10).eq.'delta_phi ')then
               delta_phiplot(plotclass)=plot
            elseif (plotname(1:10).eq.'delta_eta ')then
               delta_etaplot(plotclass)=plot
            elseif (plotname(1:5).eq.'kTij ')then
               kTijplot(plotclass)=plot
            elseif (plotname(1:4).eq.'XY1 ')then
               XY1plot(plotclass)=plot
            elseif (plotname(1:4).eq.'XY2 ')then
               XY2plot(plotclass)=plot
            elseif (plotname(1:4).eq.'XY3 ')then
               XY3plot(plotclass)=plot
            elseif (plotname(1:5).eq.'XYZ1 ')then
               XYZ1plot(plotclass)=plot
            elseif (plotname(1:5).eq.'XYZ2 ')then
               XYZ2plot(plotclass)=plot
            elseif (plotname(1:5).eq.'XYZ3 ')then
               XYZ3plot(plotclass)=plot
            elseif (plotname(1:6).eq.'XYZA1 ')then
               XYZA1plot(plotclass)=plot
            elseif (plotname(1:6).eq.'XYZA2 ')then
               XYZA2plot(plotclass)=plot
            elseif (plotname(1:6).eq.'XYZA3 ')then
               XYZA3plot(plotclass)=plot
            endif
         goto 50  !Not yet at end of PlotOptions -> Read next line
 53      continue !Found all PlotOptions ->  continue here
      endif     
      goto 51   !Go back to beginning to read next line
 80   continue  !End of file -> continue


      end







      subroutine histo_init
c*************************************************************************
c     Set up graphing
c*************************************************************************
      implicit none
      include 'info.inc'
c
c     LOCAL
c
      integer i,j,l,npar,ntot,k,n,m,a,b,c,d
      character*1 cparti,cpartj
      character*20 title
c
c     Global
c
      integer nexternal
      character*4 id(-300:300)
      common/to_id/id
      integer maxclass,maxevparticles,partoclass(-300:300,0:MaxParticles)
      integer maxevpart(0:maxclasses),maxpart(0:maxclasses)
      common/classes/maxclass,nexternal,partoclass,maxevpart

      call inihist

      call setplotrange
C
      print*,'Setting up graphs'
c

      write (52,*)''
      write (52,*)'plots set up:'

C
C     PLOTS 
C
      l=0      !counter of the number of plots
c-- weight
      l=l+1
      title='Weights'
      write (52,*) title
      call mbook(l,title,0.02d0,0d0,2d0)
c-- counts
      do i=1,maxclass
         if (maxevpart(i).ge.1)then
            l=l+1
            title='# '//name(i,1)(1:lname(i,1)-1)
      write (52,*) title
            call mbook(l,title,1d0,-0.5d0,10.5d0)
         endif
      enddo
c-- Ht
      l=l+1
      title='Ht'
      write (52,*) title
      call mbook(l,title,10d0,0d0,500d0)

c-- Missing ET
      if (maxevpart(0).ge.1) then
         l=l+1
         title='Missing ET'
      write (52,*) title
         call mbook(l,title,etmissrange(1),etmissrange(2),etmissrange(3))
      endif
c-- pt
      do i=1, maxclass 
         if (ptplot(i).ne.0) then
            j=1
            do while (j.le.ptplot(i).and.j.le.maxevpart(i))
               l=l+1
               title='pt('//name(i,j)(1:lname(i,j))//')'
      write (52,*) title
               call mbook(l,title,ptrange(1),ptrange(2),ptrange(3))
               j=j+1
            enddo
         endif
      enddo
c-- et
      do i=1, maxclass 
         if (etplot(i).ne.0) then
            j=1
            do while (j.le.etplot(i).and.j.le.maxevpart(i))
               l=l+1
               title='et('//name(i,j)(1:lname(i,j))//')'
      write (52,*) title
               call mbook(l,title,etrange(1),etrange(2),etrange(3))
               j=j+1
            enddo
         endif
      enddo
c-- e
      do i=1, maxclass 
         if (eplot(i).ne.0) then
            j=1
            do while (j.le.eplot(i).and.j.le.maxevpart(i))
               l=l+1
               title='e('//name(i,j)(1:lname(i,j))//')'
      write (52,*) title
               call mbook(l,title,erange(1),erange(2),erange(3))
               j=j+1
            enddo
         endif
      enddo
c-- momentum
      do i=1, maxclass 
         if (momplot(i).ne.0) then
            j=1
            do while (j.le.momplot(i).and.j.le.maxevpart(i))
               l=l+1
               title='mom|'//name(i,j)(1:lname(i,j))//'|'
      write (52,*) title
               call mbook(l,title,momrange(1),momrange(2),momrange(3))
               j=j+1
            enddo
         endif
      enddo
c-- pseudo-rapidity eta
      do i=1, maxclass 
         if (etaplot(i).ne.0) then
            j=1
            do while (j.le.etaplot(i).and.j.le.maxevpart(i))
               l=l+1
               title='eta('//name(i,j)(1:lname(i,j))//')'
      write (52,*) title
               call mbook(l,title,etarange(1),etarange(2),etarange(3))
               j=j+1
            enddo
         endif
      enddo
c-- rapidity y
      do i=1, maxclass 
         if (yplot(i).ne.0) then
            j=1
            do while (j.le.yplot(i).and.j.le.maxevpart(i))
               l=l+1
               title='y('//name(i,j)(1:lname(i,j))//')'
      write (52,*) title
               call mbook(l,title,yrange(1),yrange(2),yrange(3))
               j=j+1
            enddo
         endif
      enddo
c-- cos theta
      do i=1, maxclass 
         if (costhplot(i).ne.0) then
            j=1
            do while (j.le.costhplot(i).and.j.le.maxevpart(i))
               l=l+1
               title='costh('//name(i,j)(1:lname(i,j))//')'
      write (52,*) title
               call mbook(l,title,costhrange(1),costhrange(2),costhrange(3))
               j=j+1
            enddo
         endif
      enddo
c-- azymuthal phi
      do i=1, maxclass 
         if (phiplot(i).ne.0) then
            j=1
            do while (j.le.phiplot(i).and.j.le.maxevpart(i))
               l=l+1
               title='phi('//name(i,j)(1:lname(i,j))//')'
      write (52,*) title
               call mbook(l,title,phirange(1),phirange(2),phirange(3))
               j=j+1
            enddo
         endif
      enddo
c-- User defined
      do i=1, maxclass 
         if (X1plot(i).ne.0) then
            j=1
            do while (j.le.X1plot(i).and.j.le.maxevpart(i))
               l=l+1
               title='X1('//name(i,j)(1:lname(i,j))//')'
      write (52,*) title
               call mbook(l,title,X1range(1),X1range(2),X1range(3))
               j=j+1
            enddo
         endif
      enddo
      do i=1, maxclass 
         if (X2plot(i).ne.0) then
            j=1
            do while (j.le.X2plot(i).and.j.le.maxevpart(i))
               l=l+1
               title='X2('//name(i,j)(1:lname(i,j))//')'
      write (52,*) title
               call mbook(l,title,X2range(1),X2range(2),X2range(3))
               j=j+1
            enddo
         endif
      enddo
      do i=1, maxclass 
         if (X3plot(i).ne.0) then
            j=1
            do while (j.le.X3plot(i).and.j.le.maxevpart(i))
               l=l+1
               title='X3('//name(i,j)(1:lname(i,j))//')'
      write (52,*) title
               call mbook(l,title,X3range(1),X3range(2),X3range(3))
               j=j+1
            enddo
         endif
      enddo
c-- DeltaR
      do i=1, maxclass
         if (dRijplot(i).ne.0) then
            do m=i, maxclass
               if (dRijplot(m).ne.0) then
                  j=1
                  do while (j.le.dRijplot(i).and.j.le.maxevpart(i))
                     n=1
                     do while (n.le.dRijplot(m).and.n.le.maxevpart(m))
                        if (i.ne.m.or.j.lt.n)then ! Removes double counting
                           l=l+1
                           title='R('//name(i,j)(1:lname(i,j))//','//name(m,n)(1:lname(m,n))//')'
      write (52,*) title
                           call mbook(l,title,dRijrange(1),dRijrange(2),dRijrange(3))
                        endif
                        n=n+1
                     enddo
                     j=j+1
                  enddo
               endif
            enddo
         endif
      enddo
c-- invariant mass
      do i=1, maxclass
         if (mijplot(i).ne.0) then
            do m=i, maxclass
               if (mijplot(m).ne.0) then
                  j=1
                  do while (j.le.mijplot(i).and.j.le.maxevpart(i))
                     n=1
                     do while (n.le.mijplot(m).and.n.le.maxevpart(m))
                        if (i.ne.m.or.j.lt.n)then ! Removes double counting
                           l=l+1
                           title='m('//name(i,j)(1:lname(i,j))//','//name(m,n)(1:lname(m,n))//')'
      write (52,*) title
                           call mbook(l,title,mijrange(1),mijrange(2),mijrange(3))
                        endif
                        n=n+1
                     enddo
                     j=j+1
                  enddo
               endif
            enddo
         endif
      enddo
c-- delta_phi
      do i=1, maxclass
         if (delta_phiplot(i).ne.0) then
            do m=i, maxclass
               if (delta_phiplot(m).ne.0) then
                  j=1
                  do while (j.le.delta_phiplot(i).and.j.le.maxevpart(i))
                     n=1
                     do while (n.le.delta_phiplot(m).and.n.le.maxevpart(m))
                        if (i.ne.m.or.j.lt.n)then ! Removes double counting
                           l=l+1
                           title='Dphi('//name(i,j)(1:lname(i,j))//','//name(m,n)(1:lname(m,n))//')'
      write (52,*) title
                           call mbook(l,title,delta_phirange(1),delta_phirange(2),delta_phirange(3))
                        endif
                        n=n+1
                     enddo
                     j=j+1
                  enddo
               endif
            enddo
         endif
      enddo
c-- delta_eta
      do i=1, maxclass
         if (delta_etaplot(i).ne.0) then
            do m=i, maxclass
               if (delta_etaplot(m).ne.0) then
                  j=1
                  do while (j.le.delta_etaplot(i).and.j.le.maxevpart(i))
                     n=1
                     do while (n.le.delta_etaplot(m).and.n.le.maxevpart(m))
                        if (i.ne.m.or.j.lt.n)then ! Removes double counting
                           l=l+1
                           title='Deta('//name(i,j)(1:lname(i,j))//','//name(m,n)(1:lname(m,n))//')'
      write (52,*) title
                           call mbook(l,title,delta_etarange(1),delta_etarange(2),delta_etarange(3))
                        endif
                        n=n+1
                     enddo
                     j=j+1
                  enddo
               endif
            enddo
         endif
      enddo
c-- kTij
      do i=1, maxclass
         if (kTijplot(i).ne.0) then
            do m=i, maxclass
               if (kTijplot(m).ne.0) then
                  j=1
                  do while (j.le.kTijplot(i).and.j.le.maxevpart(i))
                     n=1
                     do while (n.le.kTijplot(m).and.n.le.maxevpart(m))
                        if (i.ne.m.or.j.lt.n)then ! Removes double counting
                           l=l+1
                           title='kTij('//name(i,j)(1:lname(i,j))//','//name(m,n)(1:lname(m,n))//')'
      write (52,*) title
                           call mbook(l,title,kTijrange(1),kTijrange(2),kTijrange(3))
                        endif
                        n=n+1
                     enddo
                     j=j+1
                  enddo
               endif
            enddo
         endif
      enddo
c-- User defined
      do i=1, maxclass
         if (XY1plot(i).ne.0) then
            do m=i, maxclass
               if (XY1plot(m).ne.0) then
                  j=1
                  do while (j.le.XY1plot(i).and.j.le.maxevpart(i))
                     n=1
                     do while (n.le.XY1plot(m).and.n.le.maxevpart(m))
                        if (i.ne.m.or.j.lt.n)then ! Removes double counting
                           l=l+1
                           title='XY1('//name(i,j)(1:lname(i,j))//','//name(m,n)(1:lname(m,n))//')'
      write (52,*) title
                           call mbook(l,title,XY1range(1),XY1range(2),XY1range(3))
                        endif
                        n=n+1
                     enddo
                     j=j+1
                  enddo
               endif
            enddo
         endif
      enddo
      do i=1, maxclass
         if (XY2plot(i).ne.0) then
            do m=i, maxclass
               if (XY2plot(m).ne.0) then
                  j=1
                  do while (j.le.XY2plot(i).and.j.le.maxevpart(i))
                     n=1
                     do while (n.le.XY2plot(m).and.n.le.maxevpart(m))
                        if (i.ne.m.or.j.lt.n)then ! Removes double counting
                           l=l+1
                           title='XY2('//name(i,j)(1:lname(i,j))//','//name(m,n)(1:lname(m,n))//')'
      write (52,*) title
                           call mbook(l,title,XY2range(1),XY2range(2),XY2range(3))
                        endif
                        n=n+1
                     enddo
                     j=j+1
                  enddo
               endif
            enddo
         endif
      enddo
      do i=1, maxclass
         if (XY3plot(i).ne.0) then
            do m=i, maxclass
               if (XY3plot(m).ne.0) then
                  j=1
                  do while (j.le.XY3plot(i).and.j.le.maxevpart(i))
                     n=1
                     do while (n.le.XY3plot(m).and.n.le.maxevpart(m))
                        if (i.ne.m.or.j.lt.n)then ! Removes double counting
                           l=l+1
                           title='XY3('//name(i,j)(1:lname(i,j))//','//name(m,n)(1:lname(m,n))//')'
      write (52,*) title
                           call mbook(l,title,XY3range(1),XY3range(2),XY3range(3))
                        endif
                        n=n+1
                     enddo
                     j=j+1
                  enddo
               endif
            enddo
         endif
      enddo
c-- User defined
      do i=1, maxclass
         if (XYZ1plot(i).ne.0) then
            do m=i, maxclass
               if (XYZ1plot(m).ne.0) then
                  do a=m, maxclass
                     if (XYZ1plot(a).ne.0) then
                        j=1
                        do while (j.le.XYZ1plot(i).and.j.le.maxevpart(i))
                           n=1
                           do while (n.le.XYZ1plot(m).and.n.le.maxevpart(m))
                              b=1
                              do while (b.le.XYZ1plot(a).and.b.le.maxevpart(a))
                                 if ((i.ne.m.or.j.lt.n).and.(m.ne.a.or.n.lt.b))then ! Removes double counting
                                    l=l+1
                                    title='XYZ1('//name(i,j)(1:lname(i,j))//
     &                                    ','//name(m,n)(1:lname(m,n))//','//name(a,b)(1:lname(a,b))//')'
      write (52,*) title
                                    call mbook(l,title,XYZ1range(1),XYZ1range(2),XYZ1range(3))
                                 endif
                                 b=b+1
                              enddo
                              n=n+1
                           enddo
                           j=j+1
                        enddo
                     endif
                  enddo
               endif
            enddo
         endif
      enddo
      do i=1, maxclass
         if (XYZ2plot(i).ne.0) then
            do m=i, maxclass
               if (XYZ2plot(m).ne.0) then
                  do a=m, maxclass
                     if (XYZ2plot(a).ne.0) then
                        j=1
                        do while (j.le.XYZ2plot(i).and.j.le.maxevpart(i))
                           n=1
                           do while (n.le.XYZ2plot(m).and.n.le.maxevpart(m))
                              b=1
                              do while (b.le.XYZ2plot(a).and.b.le.maxevpart(a))
                                 if ((i.ne.m.or.j.lt.n).and.(m.ne.a.or.n.lt.b))then ! Removes double counting
                                    l=l+1
                                    title='XYZ2('//name(i,j)(1:lname(i,j))//
     &                                    ','//name(m,n)(1:lname(m,n))//','//name(a,b)(1:lname(a,b))//')'
      write (52,*) title
                                    call mbook(l,title,XYZ2range(1),XYZ2range(2),XYZ2range(3))
                                 endif
                                 b=b+1
                              enddo
                              n=n+1
                           enddo
                           j=j+1
                        enddo
                     endif
                  enddo
               endif
            enddo
         endif
      enddo
      do i=1, maxclass
         if (XYZ3plot(i).ne.0) then
            do m=i, maxclass
               if (XYZ3plot(m).ne.0) then
                  do a=m, maxclass
                     if (XYZ3plot(a).ne.0) then
                        j=1
                        do while (j.le.XYZ3plot(i).and.j.le.maxevpart(i))
                           n=1
                           do while (n.le.XYZ3plot(m).and.n.le.maxevpart(m))
                              b=1
                              do while (b.le.XYZ3plot(a).and.b.le.maxevpart(a))
                                 if ((i.ne.m.or.j.lt.n).and.(m.ne.a.or.n.lt.b))then ! Removes double counting
                                    l=l+1
                                    title='XYZ3('//name(i,j)(1:lname(i,j))//
     &                                    ','//name(m,n)(1:lname(m,n))//','//name(a,b)(1:lname(a,b))//')'
      write (52,*) title
                                    call mbook(l,title,XYZ3range(1),XYZ3range(2),XYZ3range(3))
                                 endif
                                 b=b+1
                              enddo
                              n=n+1
                           enddo
                           j=j+1
                        enddo
                     endif
                  enddo
               endif
            enddo
         endif
      enddo
      do i=1, maxclass
         if (XYZA1plot(i).ne.0) then
            do m=i, maxclass
               if (XYZA1plot(m).ne.0) then
                  do a=m, maxclass
                     if (XYZA1plot(a).ne.0) then
                        do c=a, maxclass
                           if (XYZA1plot(c).ne.0) then
                              j=1
                              do while (j.le.XYZA1plot(i).and.j.le.maxevpart(i))
                                 n=1
                                 do while (n.le.XYZA1plot(m).and.n.le.maxevpart(m))
                                    b=1
                                    do while (b.le.XYZA1plot(a).and.b.le.maxevpart(a))
                                       d=1
                                       do while (d.le.XYZA1plot(c).and.d.le.maxevpart(c))
                                          if ((i.ne.m.or.j.lt.n).and.(m.ne.a.or.n.lt.b).and.(a.ne.c.or.b.lt.d))then
                                             l=l+1
                                             title='XYZA1('//name(i,j)(1:lname(i,j))//','//name(m,n)(1:lname(m,n))//
     &                                            ','//name(a,b)(1:lname(a,b))//','//name(c,d)(1:lname(c,d))//')'
      write (52,*) title
                                             call mbook(l,title,XYZA1range(1),XYZA1range(2),XYZA1range(3))
                                          endif
                                          d=d+1
                                       enddo
                                       b=b+1
                                    enddo
                                    n=n+1
                                 enddo
                                 j=j+1
                              enddo
                           endif
                        enddo
                     endif
                  enddo
               endif
            enddo
         endif
      enddo
      do i=1, maxclass
         if (XYZA2plot(i).ne.0) then
            do m=i, maxclass
               if (XYZA2plot(m).ne.0) then
                  do a=m, maxclass
                     if (XYZA2plot(a).ne.0) then
                        do c=a, maxclass
                           if (XYZA2plot(c).ne.0) then
                              j=1
                              do while (j.le.XYZA2plot(i).and.j.le.maxevpart(i))
                                 n=1
                                 do while (n.le.XYZA2plot(m).and.n.le.maxevpart(m))
                                    b=1
                                    do while (b.le.XYZA2plot(a).and.b.le.maxevpart(a))
                                       d=1
                                       do while (d.le.XYZA2plot(c).and.d.le.maxevpart(c))
                                          if ((i.ne.m.or.j.lt.n).and.(m.ne.a.or.n.lt.b).and.(a.ne.c.or.b.lt.d))then
                                             l=l+1
                                             title='XYZA2('//name(i,j)(1:lname(i,j))//','//name(m,n)(1:lname(m,n))//
     &                                            ','//name(a,b)(1:lname(a,b))//','//name(c,d)(1:lname(c,d))//')'
      write (52,*) title
                                             call mbook(l,title,XYZA2range(1),XYZA2range(2),XYZA2range(3))
                                          endif
                                          d=d+1
                                       enddo
                                       b=b+1
                                    enddo
                                    n=n+1
                                 enddo
                                 j=j+1
                              enddo
                           endif
                        enddo
                     endif
                  enddo
               endif
            enddo
         endif
      enddo
      do i=1, maxclass
         if (XYZA3plot(i).ne.0) then
            do m=i, maxclass
               if (XYZA3plot(m).ne.0) then
                  do a=m, maxclass
                     if (XYZA3plot(a).ne.0) then
                        do c=a, maxclass
                           if (XYZA3plot(c).ne.0) then
                              j=1
                              do while (j.le.XYZA3plot(i).and.j.le.maxevpart(i))
                                 n=1
                                 do while (n.le.XYZA3plot(m).and.n.le.maxevpart(m))
                                    b=1
                                    do while (b.le.XYZA3plot(a).and.b.le.maxevpart(a))
                                       d=1
                                       do while (d.le.XYZA3plot(c).and.d.le.maxevpart(c))
                                          if ((i.ne.m.or.j.lt.n).and.(m.ne.a.or.n.lt.b).and.(a.ne.c.or.b.lt.d))then
                                             l=l+1
                                             title='XYZA3('//name(i,j)(1:lname(i,j))//','//name(m,n)(1:lname(m,n))//
     &                                            ','//name(a,b)(1:lname(a,b))//','//name(c,d)(1:lname(c,d))//')'
      write (52,*) title
                                             call mbook(l,title,XYZA3range(1),XYZA3range(2),XYZA3range(3))
                                          endif
                                          d=d+1
                                       enddo
                                       b=b+1
                                    enddo
                                    n=n+1
                                 enddo
                                 j=j+1
                              enddo
                           endif
                        enddo
                     endif
                  enddo
               endif
            enddo
         endif
      enddo

     
      write(*,*) l," plots set up"

      end


      subroutine graph_event(p,dwgt,mET)
c***************************************************************************
c     fill historgrams
c
c     Use eventinfo(i,j) to specify particles.
c     i              class  (i=0 for missing Et class)
c     j              particle in class (ordered from large to small pt)
c     eventinfo(i,j) is equal to 0 if particle does not exists
c***************************************************************************
      implicit none
      include 'info.inc'
      double precision  pi              , to_deg
      parameter        (pi = 3.1415927d0, to_deg=180d0/pi)
c
c     Arguments
c
      double precision dwgt,mET
      REAL*8 P(0:4,MaxParticles)
c
c     Local
c
      real*8 www,mij,R,ptv(MaxParticles)
      integer i,j,k,l,lmax,m,n,a,b,c,d
      integer idn(MaxParticles),npart(MaxParticles),ntot,npar
      real*8 dj,Ht
c
c     External
c
      double precision sumdot,dot,et,pt,eta,dRij,delta_phi,y,e,mom,costh,phi,delta_eta,kTij
      double precision X1,X2,X3,XY1,XY2,XY3,XYZ1,XYZ2,XYZ3,XYZA1,XYZA2,XYZA3
c
c     Global
c
      integer maxclass,maxevparticles,partoclass(-300:300,0:MaxParticles)
      integer maxevpart(0:maxclasses),maxpart(0:maxclasses)
      common/classes/maxclass,maxevparticles,partoclass,maxevpart

c-----
c  Begin Code
c-----

      www=dwgt


      l=0 !counter of the number of plots


c-- weights
      l=l+1
      call mfill(l,www,www)

c-- counts
      do i=1,maxclass
         if (maxevpart(i).ge.1)then
         l=l+1
         j=1
         do while (eventinfo(i,j).ne.0)
            j=j+1
         enddo
         dj=real(j-1)
         call mfill(l,dj,www)
         endif
      enddo
c-- Ht
      l=l+1
      Ht=0d0
      do i=1,maxclass  ! Loop over all classes
         do j=1,maxevpart(i)   ! Loop from 1 to max # particles in any event (belonging to class i)
            if (eventinfo(i,j).ne.0)then       ! Do not sum if particle does not exist in an event
               Ht=Ht+et(p(0,eventinfo(i,j)))
            endif
         enddo
      enddo
      Ht=Ht+mET
      call mfill(l,Ht,www)
c-- Missing ET
      if (maxevpart(0).ge.1) then
         l=l+1
         call mfill(l,mET,www)
      endif
c-- pt
      do i=1, maxclass 
         if (ptplot(i).ne.0) then
            j=1
            do while (j.le.ptplot(i).and.j.le.maxevpart(i))
               l=l+1
               if (eventinfo(i,j).ne.0)then ! Do not plot if particle does not exist in an event
                  call mfill(l,pt(p(0,eventinfo(i,j))),www)
               endif
               j=j+1
            enddo
         endif
      enddo
c-- et
      do i=1, maxclass 
         if (etplot(i).ne.0) then
            j=1
            do while (j.le.etplot(i).and.j.le.maxevpart(i))
               l=l+1
               if (eventinfo(i,j).ne.0)then ! Do not plot if particle does not exist in an event
                  call mfill(l,et(p(0,eventinfo(i,j))),www)
               endif
               j=j+1
            enddo
         endif
      enddo
c-- e
      do i=1, maxclass 
         if (eplot(i).ne.0) then
            j=1
            do while (j.le.eplot(i).and.j.le.maxevpart(i))
               l=l+1
               if (eventinfo(i,j).ne.0)then ! Do not plot if particle does not exist in an event
                  call mfill(l,e(p(0,eventinfo(i,j))),www)
               endif
               j=j+1
            enddo
         endif
      enddo
c-- mom
      do i=1, maxclass 
         if (momplot(i).ne.0) then
            j=1
            do while (j.le.momplot(i).and.j.le.maxevpart(i))
               l=l+1
               if (eventinfo(i,j).ne.0)then ! Do not plot if particle does not exist in an event
                  call mfill(l,mom(p(0,eventinfo(i,j))),www)
               endif
               j=j+1
            enddo
         endif
      enddo
c-- pseudo-rapidity eta
      do i=1, maxclass 
         if (etaplot(i).ne.0) then
            j=1
            do while (j.le.etaplot(i).and.j.le.maxevpart(i))
               l=l+1
               if (eventinfo(i,j).ne.0)then ! Do not plot if particle does not exist in an event
                  call mfill(l,eta(p(0,eventinfo(i,j))),www)
               endif
               j=j+1
            enddo
         endif
      enddo
c-- rapidity y
      do i=1, maxclass 
         if (yplot(i).ne.0) then
            j=1
            do while (j.le.yplot(i).and.j.le.maxevpart(i))
               l=l+1
               if (eventinfo(i,j).ne.0)then ! Do not plot if particle does not exist in an event
                  call mfill(l,y(p(0,eventinfo(i,j))),www)
               endif
               j=j+1
            enddo
         endif
      enddo
c-- cos theta
      do i=1, maxclass 
         if (costhplot(i).ne.0) then
            j=1
            do while (j.le.costhplot(i).and.j.le.maxevpart(i))
               l=l+1
               if (eventinfo(i,j).ne.0)then ! Do not plot if particle does not exist in an event
                  call mfill(l,costh(p(0,eventinfo(i,j))),www)
               endif
               j=j+1
            enddo
         endif
      enddo
c-- azymuthal phi
      do i=1, maxclass 
         if (phiplot(i).ne.0) then
            j=1
            do while (j.le.phiplot(i).and.j.le.maxevpart(i))
               l=l+1
               if (eventinfo(i,j).ne.0)then ! Do not plot if particle does not exist in an event
                  call mfill(l,phi(p(0,eventinfo(i,j))),www)
               endif
               j=j+1
            enddo
         endif
      enddo
c-- User defined
      do i=1, maxclass 
         if (X1plot(i).ne.0) then
            j=1
            do while (j.le.X1plot(i).and.j.le.maxevpart(i))
               l=l+1
               if (eventinfo(i,j).ne.0)then ! Do not plot if particle does not exist in an event
                  call mfill(l,X1(p(0,eventinfo(i,j))),www)
               endif
               j=j+1
            enddo
         endif
      enddo
      do i=1, maxclass 
         if (X2plot(i).ne.0) then
            j=1
            do while (j.le.X2plot(i).and.j.le.maxevpart(i))
               l=l+1
               if (eventinfo(i,j).ne.0)then ! Do not plot if particle does not exist in an event
                  call mfill(l,X2(p(0,eventinfo(i,j))),www)
               endif
               j=j+1
            enddo
         endif
      enddo
      do i=1, maxclass 
         if (X3plot(i).ne.0) then
            j=1
            do while (j.le.X3plot(i).and.j.le.maxevpart(i))
               l=l+1
               if (eventinfo(i,j).ne.0)then ! Do not plot if particle does not exist in an event
                  call mfill(l,X3(p(0,eventinfo(i,j))),www)
               endif
               j=j+1
            enddo
         endif
      enddo
c-- DeltaR
      do i=1, maxclass
         if (dRijplot(i).ne.0) then
            do m=i, maxclass
               if (dRijplot(m).ne.0) then
                  j=1
                  do while (j.le.dRijplot(i).and.j.le.maxevpart(i))
                     n=1
                     do while (n.le.dRijplot(m).and.n.le.maxevpart(m))
                        if (i.ne.m.or.j.lt.n)then ! Removes double counting
                           l=l+1
                           if (eventinfo(i,j).ne.0.and.eventinfo(m,n).ne.0)then
                              call mfill(l,dRij(p(0,eventinfo(i,j)),p(0,eventinfo(m,n))),www)
                           endif
                        endif
                        n=n+1
                     enddo
                     j=j+1
                  enddo
               endif
            enddo
         endif
      enddo
c-- invariant mass
      do i=1, maxclass
         if (mijplot(i).ne.0) then
            do m=i, maxclass
               if (mijplot(m).ne.0) then
                  j=1
                  do while (j.le.mijplot(i).and.j.le.maxevpart(i))
                     n=1
                     do while (n.le.mijplot(m).and.n.le.maxevpart(m))
                        if (i.ne.m.or.j.lt.n)then ! Removes double counting
                           l=l+1
                           if (eventinfo(i,j).ne.0.and.eventinfo(m,n).ne.0)then
                              call mfill(l,mij(p(0,eventinfo(i,j)),p(0,eventinfo(m,n))),www)
                           endif
                        endif
                        n=n+1
                     enddo
                     j=j+1
                  enddo
               endif
            enddo
         endif
      enddo
c-- delta_phi
      do i=1, maxclass
         if (delta_phiplot(i).ne.0) then
            do m=i, maxclass
               if (delta_phiplot(m).ne.0) then
                  j=1
                  do while (j.le.delta_phiplot(i).and.j.le.maxevpart(i))
                     n=1
                     do while (n.le.delta_phiplot(m).and.n.le.maxevpart(m))
                        if (i.ne.m.or.j.lt.n)then ! Removes double counting
                           l=l+1
                           if (eventinfo(i,j).ne.0.and.eventinfo(m,n).ne.0)then
                              call mfill(l,delta_phi(p(0,eventinfo(i,j)),p(0,eventinfo(m,n))),www)
                           endif
                        endif
                        n=n+1
                     enddo
                     j=j+1
                  enddo
               endif
            enddo
         endif
      enddo
c-- delta_eta
      do i=1, maxclass
         if (delta_etaplot(i).ne.0) then
            do m=i, maxclass
               if (delta_etaplot(m).ne.0) then
                  j=1
                  do while (j.le.delta_etaplot(i).and.j.le.maxevpart(i))
                     n=1
                     do while (n.le.delta_etaplot(m).and.n.le.maxevpart(m))
                        if (i.ne.m.or.j.lt.n)then ! Removes double counting
                           l=l+1
                           if (eventinfo(i,j).ne.0.and.eventinfo(m,n).ne.0)then
                              call mfill(l,delta_eta(p(0,eventinfo(i,j)),p(0,eventinfo(m,n))),www)
                           endif
                        endif
                        n=n+1
                     enddo
                     j=j+1
                  enddo
               endif
            enddo
         endif
      enddo
c-- kTij
      do i=1, maxclass
         if (kTijplot(i).ne.0) then
            do m=i, maxclass
               if (kTijplot(m).ne.0) then
                  j=1
                  do while (j.le.kTijplot(i).and.j.le.maxevpart(i))
                     n=1
                     do while (n.le.kTijplot(m).and.n.le.maxevpart(m))
                        if (i.ne.m.or.j.lt.n)then ! Removes double counting
                           l=l+1
                           if (eventinfo(i,j).ne.0.and.eventinfo(m,n).ne.0)then
                              call mfill(l,kTij(p(0,eventinfo(i,j)),p(0,eventinfo(m,n))),www)
                           endif
                        endif
                        n=n+1
                     enddo
                     j=j+1
                  enddo
               endif
            enddo
         endif
      enddo
c-- user defined
      do i=1, maxclass
         if (XY1plot(i).ne.0) then
            do m=i, maxclass
               if (XY1plot(m).ne.0) then
                  j=1
                  do while (j.le.XY1plot(i).and.j.le.maxevpart(i))
                     n=1
                     do while (n.le.XY1plot(m).and.n.le.maxevpart(m))
                        if (i.ne.m.or.j.lt.n)then ! Removes double counting
                           l=l+1
                           if (eventinfo(i,j).ne.0.and.eventinfo(m,n).ne.0)then
                              call mfill(l,XY1(p(0,eventinfo(i,j)),p(0,eventinfo(m,n))),www)
                           endif
                        endif
                        n=n+1
                     enddo
                     j=j+1
                  enddo
               endif
            enddo
         endif
      enddo
      do i=1, maxclass
         if (XY2plot(i).ne.0) then
            do m=i, maxclass
               if (XY2plot(m).ne.0) then
                  j=1
                  do while (j.le.XY2plot(i).and.j.le.maxevpart(i))
                     n=1
                     do while (n.le.XY2plot(m).and.n.le.maxevpart(m))
                        if (i.ne.m.or.j.lt.n)then ! Removes double counting
                           l=l+1
                           if (eventinfo(i,j).ne.0.and.eventinfo(m,n).ne.0)then
                              call mfill(l,XY2(p(0,eventinfo(i,j)),p(0,eventinfo(m,n))),www)
                           endif
                        endif
                        n=n+1
                     enddo
                     j=j+1
                  enddo
               endif
            enddo
         endif
      enddo
      do i=1, maxclass
         if (XY3plot(i).ne.0) then
            do m=i, maxclass
               if (XY3plot(m).ne.0) then
                  j=1
                  do while (j.le.XY3plot(i).and.j.le.maxevpart(i))
                     n=1
                     do while (n.le.XY3plot(m).and.n.le.maxevpart(m))
                        if (i.ne.m.or.j.lt.n)then ! Removes double counting
                           l=l+1
                           if (eventinfo(i,j).ne.0.and.eventinfo(m,n).ne.0)then
                              call mfill(l,XY3(p(0,eventinfo(i,j)),p(0,eventinfo(m,n))),www)
                           endif
                        endif
                        n=n+1
                     enddo
                     j=j+1
                  enddo
               endif
            enddo
         endif
      enddo
c-- User defined
      do i=1, maxclass
         if (XYZ1plot(i).ne.0) then
            do m=i, maxclass
               if (XYZ1plot(m).ne.0) then
                  do a=m, maxclass
                     if (XYZ1plot(a).ne.0) then
                        j=1
                        do while (j.le.XYZ1plot(i).and.j.le.maxevpart(i))
                           n=1
                           do while (n.le.XYZ1plot(m).and.n.le.maxevpart(m))
                              b=1
                              do while (b.le.XYZ1plot(a).and.b.le.maxevpart(a))
                                 if ((i.ne.m.or.j.lt.n).and.(m.ne.a.or.n.lt.b))then ! Removes double counting
                                    l=l+1
                                    if (eventinfo(i,j).ne.0.and.eventinfo(m,n).ne.0.and.eventinfo(a,b).ne.0)then
                                       call mfill(l,XYZ1(p(0,eventinfo(i,j)),p(0,eventinfo(m,n)),p(0,eventinfo(a,b))),www)
                                    endif
                                 endif
                                 b=b+1
                              enddo
                              n=n+1
                           enddo
                           j=j+1
                        enddo
                     endif
                  enddo
               endif
            enddo
         endif
      enddo
      do i=1, maxclass
         if (XYZ2plot(i).ne.0) then
            do m=i, maxclass
               if (XYZ2plot(m).ne.0) then
                  do a=m, maxclass
                     if (XYZ2plot(a).ne.0) then
                        j=1
                        do while (j.le.XYZ2plot(i).and.j.le.maxevpart(i))
                           n=1
                           do while (n.le.XYZ2plot(m).and.n.le.maxevpart(m))
                              b=1
                              do while (b.le.XYZ2plot(a).and.b.le.maxevpart(a))
                                 if ((i.ne.m.or.j.lt.n).and.(m.ne.a.or.n.lt.b))then ! Removes double counting
                                    l=l+1
                                    if (eventinfo(i,j).ne.0.and.eventinfo(m,n).ne.0.and.eventinfo(a,b).ne.0)then
                                       call mfill(l,XYZ2(p(0,eventinfo(i,j)),p(0,eventinfo(m,n)),p(0,eventinfo(a,b))),www)
                                    endif
                                 endif
                                 b=b+1
                              enddo
                              n=n+1
                           enddo
                           j=j+1
                        enddo
                     endif
                  enddo
               endif
            enddo
         endif
      enddo
      do i=1, maxclass
         if (XYZ3plot(i).ne.0) then
            do m=i, maxclass
               if (XYZ3plot(m).ne.0) then
                  do a=m, maxclass
                     if (XYZ3plot(a).ne.0) then
                        j=1
                        do while (j.le.XYZ3plot(i).and.j.le.maxevpart(i))
                           n=1
                           do while (n.le.XYZ3plot(m).and.n.le.maxevpart(m))
                              b=1
                              do while (b.le.XYZ3plot(a).and.b.le.maxevpart(a))
                                 if ((i.ne.m.or.j.lt.n).and.(m.ne.a.or.n.lt.b))then ! Removes double counting
                                    l=l+1
                                    if (eventinfo(i,j).ne.0.and.eventinfo(m,n).ne.0.and.eventinfo(a,b).ne.0)then
                                       call mfill(l,XYZ3(p(0,eventinfo(i,j)),p(0,eventinfo(m,n)),p(0,eventinfo(a,b))),www)
                                    endif
                                 endif
                                 b=b+1
                              enddo
                              n=n+1
                           enddo
                           j=j+1
                        enddo
                     endif
                  enddo
               endif
            enddo
         endif
      enddo
      do i=1, maxclass
         if (XYZA1plot(i).ne.0) then
            do m=i, maxclass
               if (XYZA1plot(m).ne.0) then
                  do a=m, maxclass
                     if (XYZA1plot(a).ne.0) then
                        do c=a, maxclass
                           if (XYZA1plot(c).ne.0) then
                              j=1
                              do while (j.le.XYZA1plot(i).and.j.le.maxevpart(i))
                                 n=1
                                 do while (n.le.XYZA1plot(m).and.n.le.maxevpart(m))
                                    b=1
                                    do while (b.le.XYZA1plot(a).and.b.le.maxevpart(a))
                                       d=1
                                       do while (d.le.XYZA1plot(c).and.d.le.maxevpart(c))
                                          if ((i.ne.m.or.j.lt.n).and.(m.ne.a.or.n.lt.b).and.(a.ne.c.or.b.lt.d))then
                                             l=l+1
                                             if (eventinfo(i,j).ne.0.and.eventinfo(m,n).ne.0.and.eventinfo(a,b).ne.0)then
                                                call mfill(l,XYZA1(p(0,eventinfo(i,j)),p(0,eventinfo(m,n)),
     &                                               p(0,eventinfo(a,b)),p(0,eventinfo(c,d))),www)
                                             endif
                                          endif
                                          d=d+1
                                       enddo
                                       b=b+1
                                    enddo
                                    n=n+1
                                 enddo
                                 j=j+1
                              enddo
                           endif
                        enddo
                     endif
                  enddo
               endif
            enddo
         endif
      enddo
      do i=1, maxclass
         if (XYZA2plot(i).ne.0) then
            do m=i, maxclass
               if (XYZA2plot(m).ne.0) then
                  do a=m, maxclass
                     if (XYZA2plot(a).ne.0) then
                        do c=a, maxclass
                           if (XYZA2plot(c).ne.0) then
                              j=1
                              do while (j.le.XYZA2plot(i).and.j.le.maxevpart(i))
                                 n=1
                                 do while (n.le.XYZA2plot(m).and.n.le.maxevpart(m))
                                    b=1
                                    do while (b.le.XYZA2plot(a).and.b.le.maxevpart(a))
                                       d=1
                                       do while (d.le.XYZA2plot(c).and.d.le.maxevpart(c))
                                          if ((i.ne.m.or.j.lt.n).and.(m.ne.a.or.n.lt.b).and.(a.ne.c.or.b.lt.d))then
                                             l=l+1
                                             if (eventinfo(i,j).ne.0.and.eventinfo(m,n).ne.0.and.eventinfo(a,b).ne.0)then
                                                call mfill(l,XYZA2(p(0,eventinfo(i,j)),p(0,eventinfo(m,n)),
     &                                               p(0,eventinfo(a,b)),p(0,eventinfo(c,d))),www)
                                             endif
                                          endif
                                          d=d+1
                                       enddo
                                       b=b+1
                                    enddo
                                    n=n+1
                                 enddo
                                 j=j+1
                              enddo
                           endif
                        enddo
                     endif
                  enddo
               endif
            enddo
         endif
      enddo
      do i=1, maxclass
         if (XYZA3plot(i).ne.0) then
            do m=i, maxclass
               if (XYZA3plot(m).ne.0) then
                  do a=m, maxclass
                     if (XYZA3plot(a).ne.0) then
                        do c=a, maxclass
                           if (XYZA3plot(c).ne.0) then
                              j=1
                              do while (j.le.XYZA3plot(i).and.j.le.maxevpart(i))
                                 n=1
                                 do while (n.le.XYZA3plot(m).and.n.le.maxevpart(m))
                                    b=1
                                    do while (b.le.XYZA3plot(a).and.b.le.maxevpart(a))
                                       d=1
                                       do while (d.le.XYZA3plot(c).and.d.le.maxevpart(c))
                                          if ((i.ne.m.or.j.lt.n).and.(m.ne.a.or.n.lt.b).and.(a.ne.c.or.b.lt.d))then
                                             l=l+1
                                             if (eventinfo(i,j).ne.0.and.eventinfo(m,n).ne.0.and.eventinfo(a,b).ne.0)then
                                                call mfill(l,XYZA3(p(0,eventinfo(i,j)),p(0,eventinfo(m,n)),
     &                                               p(0,eventinfo(a,b)),p(0,eventinfo(c,d))),www)
                                             endif
                                          endif
                                          d=d+1
                                       enddo
                                       b=b+1
                                    enddo
                                    n=n+1
                                 enddo
                                 j=j+1
                              enddo
                           endif
                        enddo
                     endif
                  enddo
               endif
            enddo
         endif
      enddo


      return
      end

      subroutine histo_final
c*************************************************************************
c     Stores graphs
c*************************************************************************
      implicit none
      include 'info.inc'
c
c     Local
c
      integer k,i,j,l,n,m,a,b,c,d
      character*20 xtit,ytit
      character*3 scale
      double precision norm
c
c     Global
c
      double precision sum,awgt
      integer nevent
      common/ne/sum,awgt,nevent
      integer maxclass,nexternal,partoclass(-300:300,0:MaxParticles)
      integer maxevpart(0:maxclasses),maxpart(0:maxclasses)
      common/classes/maxclass,nexternal,partoclass,maxevpart
      
      do k=1,200
         call mfinal(k)
      enddo 
c
c Calculate normalization for the plots:
c
      norm=1d0/awgt


      scale='LOG'
      ytit ='# events/bin'

      l=0  !counter of the number of plots
c-- weights
      l=l+1
      xtit='wgt'
      call mtop(l,99,xtit,ytit,scale,nevent,sum,norm,gnu)
c-- counts
      do i=1,maxclass
         if (maxevpart(i).ge.1) then
         l=l+1
         xtit='# particles'
         call mtop(l,99,xtit,ytit,scale,nevent,sum,norm,gnu)
         endif
      enddo
c-- Ht
      l=l+1
      xtit='ET'
      call mtop(l,99,xtit,ytit,scale,nevent,sum,norm,gnu)
c-- Missing ET
      if (maxevpart(0).ge.1) then
         l=l+1
         xtit='ET'
         call mtop(l,99,xtit,ytit,scale,nevent,sum,norm,gnu)
      endif
c-- pt
      do i=1, maxclass 
         if (ptplot(i).ne.0) then
            j=1
            do while (j.le.ptplot(i).and.j.le.maxevpart(i))
               l=l+1
               xtit='pt'
               call mtop(l,99,xtit,ytit,scale,nevent,sum,norm,gnu)
               j=j+1
            enddo
         endif
      enddo
c-- et
      do i=1, maxclass 
         if (etplot(i).ne.0) then
            j=1
            do while (j.le.etplot(i).and.j.le.maxevpart(i))
               l=l+1
               xtit='et'
               call mtop(l,99,xtit,ytit,scale,nevent,sum,norm,gnu)
               j=j+1
            enddo
         endif
      enddo
c-- e
      do i=1, maxclass 
         if (eplot(i).ne.0) then
            j=1
            do while (j.le.eplot(i).and.j.le.maxevpart(i))
               l=l+1
               xtit='E'
               call mtop(l,99,xtit,ytit,scale,nevent,sum,norm,gnu)
               j=j+1
            enddo
         endif
      enddo
c-- mom
      do i=1, maxclass 
         if (momplot(i).ne.0) then
            j=1
            do while (j.le.momplot(i).and.j.le.maxevpart(i))
               l=l+1
               xtit='|momentum|'
               call mtop(l,99,xtit,ytit,scale,nevent,sum,norm,gnu)
               j=j+1
            enddo
         endif
      enddo
c-- pseudo-rapidity eta
      do i=1, maxclass 
         if (etaplot(i).ne.0) then
            j=1
            do while (j.le.etaplot(i).and.j.le.maxevpart(i))
               l=l+1
               xtit='eta'
               call mtop(l,99,xtit,ytit,'lin',nevent,sum,norm,gnu)
               j=j+1
            enddo
         endif
      enddo
c-- rapidity y
      do i=1, maxclass 
         if (yplot(i).ne.0) then
            j=1
            do while (j.le.yplot(i).and.j.le.maxevpart(i))
               l=l+1
               xtit='y'
               call mtop(l,99,xtit,ytit,'lin',nevent,sum,norm,gnu)
               j=j+1
            enddo
         endif
      enddo
c-- cos theta
      do i=1, maxclass 
         if (costhplot(i).ne.0) then
            j=1
            do while (j.le.costhplot(i).and.j.le.maxevpart(i))
               l=l+1
               xtit='cos theta'
               call mtop(l,99,xtit,ytit,'lin',nevent,sum,norm,gnu)
               j=j+1
            enddo
         endif
      enddo
c-- azymuthal phi
      do i=1, maxclass 
         if (phiplot(i).ne.0) then
            j=1
            do while (j.le.phiplot(i).and.j.le.maxevpart(i))
               l=l+1
               xtit='phi'
               call mtop(l,99,xtit,ytit,'lin',nevent,sum,norm,gnu)
               j=j+1
            enddo
         endif
      enddo
c-- user defined
      do i=1, maxclass 
         if (X1plot(i).ne.0) then
            j=1
            do while (j.le.X1plot(i).and.j.le.maxevpart(i))
               l=l+1
               xtit='X1'
               call mtop(l,99,xtit,ytit,'lin',nevent,sum,norm,gnu)
               j=j+1
            enddo
         endif
      enddo
      do i=1, maxclass 
         if (X2plot(i).ne.0) then
            j=1
            do while (j.le.X2plot(i).and.j.le.maxevpart(i))
               l=l+1
               xtit='X2'
               call mtop(l,99,xtit,ytit,'lin',nevent,sum,norm,gnu)
               j=j+1
            enddo
         endif
      enddo
      do i=1, maxclass 
         if (X3plot(i).ne.0) then
            j=1
            do while (j.le.X3plot(i).and.j.le.maxevpart(i))
               l=l+1
               xtit='X3'
               call mtop(l,99,xtit,ytit,'lin',nevent,sum,norm,gnu)
               j=j+1
            enddo
         endif
      enddo
c-- DeltaR
      do i=1, maxclass
         if (dRijplot(i).ne.0) then
            do m=i, maxclass
               if (dRijplot(m).ne.0) then
                  j=1
                  do while (j.le.dRijplot(i).and.j.le.maxevpart(i))
                     n=1
                     do while (n.le.dRijplot(m).and.n.le.maxevpart(m))
                        if (i.ne.m.or.j.lt.n)then ! Removes double counting
                           l=l+1
                           xtit='DeltaR'
                           call mtop(l,99,xtit,ytit,scale,nevent,sum,norm,gnu)
                        endif
                        n=n+1
                     enddo
                     j=j+1
                  enddo
               endif
            enddo
         endif
      enddo
c-- invariant mass
      do i=1, maxclass
         if (mijplot(i).ne.0) then
            do m=i, maxclass
               if (mijplot(m).ne.0) then
                  j=1
                  do while (j.le.mijplot(i).and.j.le.maxevpart(i))
                     n=1
                     do while (n.le.mijplot(m).and.n.le.maxevpart(m))
                        if (i.ne.m.or.j.lt.n)then ! Removes double counting
                           l=l+1
                           xtit='mij'
                           call mtop(l,99,xtit,ytit,scale,nevent,sum,norm,gnu)
                        endif
                        n=n+1
                     enddo
                     j=j+1
                  enddo
               endif
            enddo
         endif
      enddo
c-- delta_phi
      do i=1, maxclass
         if (delta_phiplot(i).ne.0) then
            do m=i, maxclass
               if (delta_phiplot(m).ne.0) then
                  j=1
                  do while (j.le.delta_phiplot(i).and.j.le.maxevpart(i))
                     n=1
                     do while (n.le.delta_phiplot(m).and.n.le.maxevpart(m))
                        if (i.ne.m.or.j.lt.n)then ! Removes double counting
                           l=l+1
                           xtit='DeltaPhi'
                           call mtop(l,99,xtit,ytit,scale,nevent,sum,norm,gnu)
                        endif
                        n=n+1
                     enddo
                     j=j+1
                  enddo
               endif
            enddo
         endif
      enddo
c-- delta_eta
      do i=1, maxclass
         if (delta_etaplot(i).ne.0) then
            do m=i, maxclass
               if (delta_etaplot(m).ne.0) then
                  j=1
                  do while (j.le.delta_etaplot(i).and.j.le.maxevpart(i))
                     n=1
                     do while (n.le.delta_etaplot(m).and.n.le.maxevpart(m))
                        if (i.ne.m.or.j.lt.n)then ! Removes double counting
                           l=l+1
                           xtit='Delta_eta'
                           call mtop(l,99,xtit,ytit,scale,nevent,sum,norm,gnu)
                        endif
                        n=n+1
                     enddo
                     j=j+1
                  enddo
               endif
            enddo
         endif
      enddo
c-- kTij
      do i=1, maxclass
         if (kTijplot(i).ne.0) then
            do m=i, maxclass
               if (kTijplot(m).ne.0) then
                  j=1
                  do while (j.le.kTijplot(i).and.j.le.maxevpart(i))
                     n=1
                     do while (n.le.kTijplot(m).and.n.le.maxevpart(m))
                        if (i.ne.m.or.j.lt.n)then ! Removes double counting
                           l=l+1
                           xtit='kTij'
                           call mtop(l,99,xtit,ytit,scale,nevent,sum,norm,gnu)
                        endif
                        n=n+1
                     enddo
                     j=j+1
                  enddo
               endif
            enddo
         endif
      enddo
c-- User defined
      do i=1, maxclass
         if (XY1plot(i).ne.0) then
            do m=i, maxclass
               if (XY1plot(m).ne.0) then
                  j=1
                  do while (j.le.XY1plot(i).and.j.le.maxevpart(i))
                     n=1
                     do while (n.le.XY1plot(m).and.n.le.maxevpart(m))
                        if (i.ne.m.or.j.lt.n)then ! Removes double counting
                           l=l+1
                           xtit='XY1'
                           call mtop(l,99,xtit,ytit,scale,nevent,sum,norm,gnu)
                        endif
                        n=n+1
                     enddo
                     j=j+1
                  enddo
               endif
            enddo
         endif
      enddo
      do i=1, maxclass
         if (XY2plot(i).ne.0) then
            do m=i, maxclass
               if (XY2plot(m).ne.0) then
                  j=1
                  do while (j.le.XY2plot(i).and.j.le.maxevpart(i))
                     n=1
                     do while (n.le.XY2plot(m).and.n.le.maxevpart(m))
                        if (i.ne.m.or.j.lt.n)then ! Removes double counting
                           l=l+1
                           xtit='XY2'
                           call mtop(l,99,xtit,ytit,scale,nevent,sum,norm,gnu)
                        endif
                        n=n+1
                     enddo
                     j=j+1
                  enddo
               endif
            enddo
         endif
      enddo
      do i=1, maxclass
         if (XY3plot(i).ne.0) then
            do m=i, maxclass
               if (XY3plot(m).ne.0) then
                  j=1
                  do while (j.le.XY3plot(i).and.j.le.maxevpart(i))
                     n=1
                     do while (n.le.XY3plot(m).and.n.le.maxevpart(m))
                        if (i.ne.m.or.j.lt.n)then ! Removes double counting
                           l=l+1
                           xtit='XY3'
                           call mtop(l,99,xtit,ytit,scale,nevent,sum,norm,gnu)
                        endif
                        n=n+1
                     enddo
                     j=j+1
                  enddo
               endif
            enddo
         endif
      enddo
c-- User defined
      do i=1, maxclass
         if (XYZ1plot(i).ne.0) then
            do m=i, maxclass
               if (XYZ1plot(m).ne.0) then
                  do a=m, maxclass
                     if (XYZ1plot(a).ne.0) then
                        j=1
                        do while (j.le.XYZ1plot(i).and.j.le.maxevpart(i))
                           n=1
                           do while (n.le.XYZ1plot(m).and.n.le.maxevpart(m))
                              b=1
                              do while (b.le.XYZ1plot(a).and.b.le.maxevpart(a))
                                 if ((i.ne.m.or.j.lt.n).and.(m.ne.a.or.n.lt.b))then ! Removes double counting
                                    l=l+1
                                    xtit='XYZ1'
                                    call mtop(l,99,xtit,ytit,scale,nevent,sum,norm,gnu)
                                 endif
                                 b=b+1
                              enddo
                              n=n+1
                           enddo
                           j=j+1
                        enddo
                     endif
                  enddo
               endif
            enddo
         endif
      enddo
      do i=1, maxclass
         if (XYZ2plot(i).ne.0) then
            do m=i, maxclass
               if (XYZ2plot(m).ne.0) then
                  do a=m, maxclass
                     if (XYZ2plot(a).ne.0) then
                        j=1
                        do while (j.le.XYZ2plot(i).and.j.le.maxevpart(i))
                           n=1
                           do while (n.le.XYZ2plot(m).and.n.le.maxevpart(m))
                              b=1
                              do while (b.le.XYZ2plot(a).and.b.le.maxevpart(a))
                                 if ((i.ne.m.or.j.lt.n).and.(m.ne.a.or.n.lt.b))then ! Removes double counting
                                    l=l+1
                                    xtit='XYZ2'
                                    call mtop(l,99,xtit,ytit,scale,nevent,sum,norm,gnu)
                                 endif
                                 b=b+1
                              enddo
                              n=n+1
                           enddo
                           j=j+1
                        enddo
                     endif
                  enddo
               endif
            enddo
         endif
      enddo
      do i=1, maxclass
         if (XYZ3plot(i).ne.0) then
            do m=i, maxclass
               if (XYZ3plot(m).ne.0) then
                  do a=m, maxclass
                     if (XYZ3plot(a).ne.0) then
                        j=1
                        do while (j.le.XYZ3plot(i).and.j.le.maxevpart(i))
                           n=1
                           do while (n.le.XYZ3plot(m).and.n.le.maxevpart(m))
                              b=1
                              do while (b.le.XYZ1plot(a).and.b.le.maxevpart(a))
                                 if ((i.ne.m.or.j.lt.n).and.(m.ne.a.or.n.lt.b))then ! Removes double counting
                                    l=l+1
                                    xtit='XYZ3'
                                    call mtop(l,99,xtit,ytit,scale,nevent,sum,norm,gnu)
                                 endif
                                 b=b+1
                              enddo
                              n=n+1
                           enddo
                           j=j+1
                        enddo
                     endif
                  enddo
               endif
            enddo
         endif
      enddo
      do i=1, maxclass
         if (XYZA1plot(i).ne.0) then
            do m=i, maxclass
               if (XYZA1plot(m).ne.0) then
                  do a=m, maxclass
                     if (XYZA1plot(a).ne.0) then
                        do c=a, maxclass
                           if (XYZA1plot(c).ne.0) then
                              j=1
                              do while (j.le.XYZA1plot(i).and.j.le.maxevpart(i))
                                 n=1
                                 do while (n.le.XYZA1plot(m).and.n.le.maxevpart(m))
                                    b=1
                                    do while (b.le.XYZA1plot(a).and.b.le.maxevpart(a))
                                       d=1
                                       do while (d.le.XYZA1plot(c).and.d.le.maxevpart(c))
                                          if ((i.ne.m.or.j.lt.n).and.(m.ne.a.or.n.lt.b).and.(a.ne.c.or.b.lt.d))then
                                             l=l+1
                                             xtit='XYZA1'
                                             call mtop(l,99,xtit,ytit,scale,nevent,sum,norm,gnu)
                                          endif
                                          d=d+1
                                       enddo
                                       b=b+1
                                    enddo
                                    n=n+1
                                 enddo
                                 j=j+1
                              enddo
                           endif
                        enddo
                     endif
                  enddo
               endif
            enddo
         endif
      enddo
      do i=1, maxclass
         if (XYZA2plot(i).ne.0) then
            do m=i, maxclass
               if (XYZA2plot(m).ne.0) then
                  do a=m, maxclass
                     if (XYZA2plot(a).ne.0) then
                        do c=a, maxclass
                           if (XYZA2plot(c).ne.0) then
                              j=1
                              do while (j.le.XYZA2plot(i).and.j.le.maxevpart(i))
                                 n=1
                                 do while (n.le.XYZA2plot(m).and.n.le.maxevpart(m))
                                    b=1
                                    do while (b.le.XYZA2plot(a).and.b.le.maxevpart(a))
                                       d=1
                                       do while (d.le.XYZA2plot(c).and.d.le.maxevpart(c))
                                          if ((i.ne.m.or.j.lt.n).and.(m.ne.a.or.n.lt.b).and.(a.ne.c.or.b.lt.d))then
                                             l=l+1
                                             xtit='XYZA2'
                                             call mtop(l,99,xtit,ytit,scale,nevent,sum,norm,gnu)
                                          endif
                                          d=d+1
                                       enddo
                                       b=b+1
                                    enddo
                                    n=n+1
                                 enddo
                                 j=j+1
                              enddo
                           endif
                        enddo
                     endif
                  enddo
               endif
            enddo
         endif
      enddo
      do i=1, maxclass
         if (XYZA3plot(i).ne.0) then
            do m=i, maxclass
               if (XYZA3plot(m).ne.0) then
                  do a=m, maxclass
                     if (XYZA3plot(a).ne.0) then
                        do c=a, maxclass
                           if (XYZA3plot(c).ne.0) then
                              j=1
                              do while (j.le.XYZA3plot(i).and.j.le.maxevpart(i))
                                 n=1
                                 do while (n.le.XYZA3plot(m).and.n.le.maxevpart(m))
                                    b=1
                                    do while (b.le.XYZA3plot(a).and.b.le.maxevpart(a))
                                       d=1
                                       do while (d.le.XYZA3plot(c).and.d.le.maxevpart(c))
                                          if ((i.ne.m.or.j.lt.n).and.(m.ne.a.or.n.lt.b).and.(a.ne.c.or.b.lt.d))then
                                             l=l+1
                                             xtit='XYZA3'
                                             call mtop(l,99,xtit,ytit,scale,nevent,sum,norm,gnu)
                                          endif
                                          d=d+1
                                       enddo
                                       b=b+1
                                    enddo
                                    n=n+1
                                 enddo
                                 j=j+1
                              enddo
                           endif
                        enddo
                     endif
                  enddo
               endif
            enddo
         endif
      enddo

      end




      subroutine set_id
c****************************************
c     to set the particles ids
*****************************************
      implicit none
c
c    Particle's ID
c     	
      character*4 id(-300:300)
      common/to_id/id
c
c
c-Quarks: d u s c b t d~ u~ s~ c~ b~ t~ 
      id(1) ='d'
      id(2) ='u'
      id(3) ='s'
      id(4) ='c'
      id(5) ='b'
      id(6) ='t'
      id(-1)='d~'
      id(-2)='u~'
      id(-3)='s~'
      id(-4)='c~'
      id(-5)='b~'
      id(-6)='t~'  
c-Leptons: e- mu- ta- ve vm vt e+ mu+ ta+ ve~ vm~ vt~ 
      id(11) ='e-'
      id(12) ='ve'
      id(13) ='mu-'
      id(14) ='vm'
      id(15) ='ta-'
      id(16) ='vt'
      id(-11)='e+'
      id(-12)='ve~'
      id(-13)='mu+'
      id(-14)='vm~'
      id(-15)='ta+'
      id(-16)='vt~'
c-Bosons: g a z w+ w- h        
      id(21)  ='j'
      id(22)  ='a'
      id(23)  ='z'
      id(24)  ='w+'
      id(-24) ='w-'
      id(25)  ='h'

c-SUSY particles
c     squarks
      id(101) = 'dl'
      id(-101) = 'dl~'
      id(201) = 'dr'
      id(-201) = 'dr~'
      id(102) = 'ul'
      id(-102) = 'ul~'
      id(202) = 'ur'
      id(-202) = 'ur~'
      id(103) = 'sl'
      id(-103) = 'sl~'
      id(203) = 'sr'
      id(-203) = 'sr~'
      id(104) = 'cl'
      id(-104) = 'cl~'
      id(204) = 'cr'
      id(-204) = 'cr~'
      id(105) = 'b1'
      id(-105) = 'b1~'
      id(205) = 'b2'
      id(-205) = 'b2~'
      id(106) = 't1'
      id(-106) = 't1~'
      id(206) = 't2'
      id(-206) = 't2~'

c     sleptons
      id(111) = 'el-'
      id(-111) = 'el+'
      id(211) = 'er-'
      id(-211) = 'er+'
      id(113) = 'mul-'
      id(-113) = 'mul+'
      id(213) = 'mur-'
      id(-213) = 'mur+'
      id(115) = 'ta1-'
      id(-115) = 'ta1+'
      id(215) = 'ta2-'
      id(-215) = 'ta2+'
      id(112) = 'sve'
      id(-112) = 'sve~'
      id(114) = 'svm'
      id(-114) = 'svm~'
      id(116) = 'svt'
      id(-116) = 'svt~'

c     higgs bosons
      id(25) = 'h1'
      id(35) = 'h2'
      id(36) = 'h3'
      id(37) = 'h-'
      id(-37) = 'h+'

c     inos
      id(121) = 'go'
      id(122) = 'n1'
      id(123) = 'n2'
      id(125) = 'n3'
      id(135) = 'n4'
      id(124) = 'x1-'
      id(-124) = 'x1+'
      id(137) = 'x2-'
      id(-137) = 'x2+'

      return
      end

      subroutine sort2(array,aux1,n)
      implicit none
! Arguments
      integer n
      integer aux1(n)
      double precision array(n)
!  Local Variables
      integer i,k
      double precision temp
      logical done

!-----------
! Begin Code
!-----------
      do i=n-1,1,-1
         done = .true.
         do k=1,i
            if (array(k) .lt. array(k+1)) then
               temp = array(k)
               array(k) = array(k+1)
               array(k+1) = temp
               temp = aux1(k)
               aux1(k) = aux1(k+1)
               aux1(k+1) = temp
               done = .false.
            end if
         end do
         if (done) return
      end do
      end 


      subroutine no_spaces(buff,nchars)
c**********************************************************************
c     Given buff a buffer of words separated by spaces
c     returns it where all space are moved to the right
c     returns also the length of the single word.
c     maxlength is the length of the buffer
c**********************************************************************
      implicit none
c
c     Constants
c
      integer    maxline
      parameter (maxline=4)
      character*1 null
      parameter  (null=' ')
c
c     Arguments
c
      character*(maxline) buff
      integer nchars,maxlength
c
c     Local
c
      integer i,j
      character*(maxline) temp
      
c-----
c  Begin Code
c-----
      nchars=0
c      write (*,*) "buff=",buff(1:maxlength)
      do i=1,maxline
         if(buff(i:i).ne.null) then
            nchars=nchars+1
            temp(nchars:nchars)=buff(i:i)
         endif
c         write(*,*) i,":",buff(1:maxlength),":",temp(1:nchars),":"
      enddo
      buff=temp      
      end

      subroutine no_spaces_long(buff,nchars)
c**********************************************************************
c     Given buff a buffer of words separated by spaces
c     returns it where all space are moved to the right
c     returns also the length of the single word.
c     maxlength is the length of the buffer
c**********************************************************************
      implicit none
c
c     Constants
c
      integer    maxline
      parameter (maxline=130)
      character*1 null
      parameter  (null=' ')
c
c     Arguments
c
      character*(maxline) buff
      integer nchars,maxlength
c
c     Local
c
      integer i,j
      character*(maxline) temp
      
c-----
c  Begin Code
c-----
      nchars=0
c      write (*,*) "buff=",buff(1:maxlength)
      do i=1,maxline
         if(buff(i:i).ne.null) then
            nchars=nchars+1
            temp(nchars:nchars)=buff(i:i)
         endif
c         write(*,*) i,":",buff(1:maxlength),":",temp(1:nchars),":"
      enddo
      buff=temp      
      end




      subroutine read_lhco_event(lun,p,wgt,nexternal,ic,done,sum)
      implicit none
      include 'info.inc'
      integer lun,nexternal,ic(7,maxparticles),final(maxparticles)
      double precision p(0:4,maxparticles),wgt,sum
      logical done

      character*256 buff

      logical lhco_banner_open
      integer lhco_lun_ban
      common/to_lhco_banner/lhco_banner_open, lhco_lun_ban
      data lhco_lun_ban/37/
      data lhco_banner_open/.false./

      integer number,type(maxparticles),i,j
      double precision eta(maxparticles),phi(maxparticles),e_y(maxparticles)
      double precision pt(maxparticles),jmas(maxparticles),ntrk(maxparticles)
      double precision btag(maxparticles),had_em(maxparticles),dum(maxparticles)
      

      done=.false.
      if (.not.lhco_banner_open)then
         read(lun,'(a256)',end=33,err=33) buff
         do while(index(buff,"  #  typ      eta ") .eq. 0)
            if (index(buff,"##  Integrated weight (pb)  :").eq. 1)then
c               buff(31:31)='0'
                write (lhco_lun_ban,'(a256)')buff
               rewind (lhco_lun_ban)
               read(lhco_lun_ban,'(a30,e11.5)')buff,sum
c               write (*,*) buff,sum
            endif
            read(lun,'(a256)') buff
         enddo
         read (lun,'(a256)')buff
         lhco_banner_open=.true.
      endif
      number=1
      i=0
      do while(number.ne.0)
         i=i+1
         read(lun,*,end=33) number,type(i),eta(i),phi(i),pt(i),
     &        jmas(i),ntrk(i),btag(i),had_em(i)
         if (number.eq.0)then
            do j=1,1
               backspace(lun)
            enddo
         endif
c     write (*,*) number, i
      enddo
 99   continue
      number=i-1
      if (number.gt.maxparticles-2)then
         write (*,*) 'Too many particles in final state'
         write (*,*) 'Increase MaxParticles to at least',number+2
      endif
      do i=1,number
         e_y(i)=(dsqrt(jmas(i)**2+pt(i)**2*dcosh(eta(i))**2)+pt(i)*dsinh(eta(i))) / dsqrt(jmas(i)**2+pt(i)**2)
         p(0,i+2)=sqrt(pt(i)**2*0.25d0*(e_y(i)+1d0/e_y(i))**2+jmas(i)**2)
         p(1,i+2)=pt(i)*cos( phi(i))
         p(2,i+2)=pt(i)*sin( phi(i))
         p(3,i+2)=pt(i)*0.5d0*(e_y(i)-1d0/e_y(i))
         p(4,i+2)=jmas(i)
         if     (type(i).eq.0)then
            ic(1,i+2)=22        !photon
         elseif (type(i).eq.1.and.ntrk(i).eq.-1.0)then
            ic(1,i+2)=11        !electron
         elseif (type(i).eq.1.and.ntrk(i).eq.1.0)then 
            ic(1,i+2)=-11       !positron
         elseif (type(i).eq.2.and.ntrk(i).eq.-1.0)then
            ic(1,i+2)=13        !mu-
         elseif (type(i).eq.2.and.ntrk(i).eq.1.0)then 
            ic(1,i+2)=-13       !mu+
         elseif (type(i).eq.3.and.ntrk(i).le.-1.0)then
            ic(1,i+2)=15        !tau-
         elseif (type(i).eq.3.and.ntrk(i).ge.1.0)then 
            ic(1,i+2)=-15       !tau+
         elseif (type(i).eq.4)then 
            if (btag(i).ge.1) then
               ic(1,i+2)=5      !bjet
            else
               ic(1,i+2)=21     !non-bjet (or gluon)
            endif
         elseif (type(i).eq.6)then
            ic(1,i+2)=12        !Missing E_t (or neutrino)
         endif
      enddo
      nexternal=number+2
      do i=1,nexternal
         do j=2,7
            ic(j,i)=1
         enddo
      enddo
      wgt=1d0
      return
 33   done=.true.
      lhco_banner_open=.false.
      return
      end




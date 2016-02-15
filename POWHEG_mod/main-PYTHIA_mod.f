      program main_pythia
      implicit none
      include 'LesHouches.h'
      include 'hepevt.h'
      integer iev,temp,i
      external pydata
      character * 6 WHCPRG
      common/cWHCPRG/WHCPRG
      integer maxpr
      parameter (maxpr=6)
      integer maxev
      common/mcmaxev/maxev
      character*100 slhafilename

      ! WHCPRG tells the analysis subroutine which program is calling the
      ! analysis
      WHCPRG='PYTHIA'

      ! how many events do we have?  
      call opencount(maxev)

      ! read in the slha file
      call powheginputstring("SLHA",slhafilename)
      call pythia_read_slha(trim(slhafilename))
      
      ! initialize pythia
      call pythia_init
      
      call init_hist
      nevhep=0
      do iev=1,maxev
         call pythia_event
         if(nup.eq.0) then
            write(*,*) ' no event generated; skipping'
            goto 111
         endif
         call pyanalysis
         if(nevhep.gt.0.and.mod(nevhep,20000).EQ.0) then
            WRITE(*,*)'# of events processed=',iev
            WRITE(*,*)'# of events generated=',nevhep
            call pyaend
         endif 
      enddo
 111  continue
      write(*,*) 'At the end nevhep is ',nevhep
      call pyaend
      end
      
c     jumps to next event and calls analysis
      subroutine pyanalysis
      implicit none
      include 'hepevt.h'
      include 'LesHouches.h'
      nevhep=nevhep+1
      if(abs(idwtup).eq.3) xwgtup=xwgtup*xsecup(1)
      call analysis(xwgtup)
      call pwhgaccumup
      end
      
c     writes PYTHIA analysis output into .top file
      subroutine pyaend
      implicit none
      character *20 pwgprefix
      character *100 filename
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      filename=pwgprefix(1:lprefix)//'POWHEG+PYTHIA-output'
      call pwhgsetout
      call pwhgtopout(filename)
      close(99)
      end
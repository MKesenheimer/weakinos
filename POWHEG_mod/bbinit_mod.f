c MK: copied and modified version of bbinit.f, revision 3154
c modified on the basis of disquark/bbinit.f
c changes marked with "! CH, MK:" or "! MK:"
c bigger changes over a whole section are marked with ! ===...

      subroutine bbinit
      implicit none
      integer iret
      real * 8 powheginput
      integer parallelstages
      external powheginput
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'cgengrids.h'
      include 'pwhg_rad.h'
      include 'pwhg_rnd.h'
c MK: added
#include "osres.h"
#include "pwhg_flst_add.h"
#include "Flags.h"
#include "cgengrids_add.h"
#include "pwhg_rad_add.h"
      real * 8 xx(ndiminteg)
      integer mcalls,icalls,j
      integer btilde,sigremnant
      external btilde,sigremnant
      integer ichan
      double precision rad_totosres_sum ! MK: added
      integer sigosres ! CH, MK: added
      external sigosres ! CH, MK: added
c parallelstages:
c 1   prepare the importance sampling grids
c 2   prepare the upper bounding envelopes for the
c     generation of the b_tilde function
c 3   prepare the upper bound for the generation of radiation
c 4   generate events
      parallelstages =  powheginput('#parallelstage')
      if(flg_newweight .and. parallelstages .gt. 0) then
         write(*,*) ' Since we are running in reweighting mode '
         write(*,*) ' we set parallelstages to 4 '
         parallelstages = 4
      endif
      if(parallelstages.gt.0.and.rnd_cwhichseed.eq.'none') then
         write(*,*) ' with parallelstage also manyseeds '
         write(*,*) ' must be set'
         call pwhg_exit(-1)
      endif
      if(parallelstages.le.0) then
c Do all stages in one go if needed;
c first look for <prefix>xgrid.dat file      
         call loadxgrids(iret)
         if(iret.ne.0) then
c if not there look for pwggridinfo-* files from parallel runs
            call loadlatestxgridinfo(iret)
            if(iret.ne.0) then
c if not there generate the xgrid; the argument 0 means:
c generate *xgrid.dat file
               call bbinitxgrids(0)
            endif
         endif
c make sure the ifold arrays are read in
         do j=1,ndiminteg
            ifold(j)=1
            ifoldrm(j)=1
            ifoldosres(j)=1 ! CH, MK: added
         enddo
c Override real integration parameters with powheg.input values         
         ifold(ndiminteg-2) = powheginput("foldcsi")
         ifold(ndiminteg-1) = powheginput("foldy")
         ifold(ndiminteg)   = powheginput("foldphi")
         if(flg_storemintupb) then
            ! CH, MK: changed the following call
            call loadgrids(iret,xgrid,ymax,ymaxrat,xgridrm,ymaxrm,
     #           ymaxratrm,xgridosres,ymaxosres,ymaxratosres,
     #           ifold,ifoldrm,ifoldosres,'fullgrid')
         else
            ! CH, MK: changed the following call
            call loadgrids(iret,xgrid,ymax,ymaxrat,xgridrm,ymaxrm,
     #           ymaxratrm,xgridosres,ymaxosres,ymaxratosres,
     #           ifold,ifoldrm,ifoldosres,'grid')
         endif
         if(iret.eq.0) then
            call apply_totarr ! MK: added
            write(*,*) 'upper bound grids successfully loaded'
            write(*,*) 'btilde pos.   weights:', rad_totposbtl,' +-',
     1           rad_etotposbtl
            write(*,*) 'btilde |neg.| weights:', rad_totnegbtl,' +-',
     1           rad_etotnegbtl
            write(*,*)
     1           'btilde total (pos.-|neg.|):', rad_totbtl,' +-',
     2           rad_etotbtl
         else
            call bbinitgrids
            !print*,"rad_totosresgen(1)",rad_totosresgen(1)
            !print*,"rad_osres_totarr(1,5,1)",rad_osres_totarr(1,5,1)
            !stop
            ! MK: gets called only if storemintupb in powheg.input is 
            ! set to 1
            if(flg_storemintupb) then
               call loadmintupb(ndiminteg,'btildeupb',ymax,ymaxrat)
               write(*,*) ' Upper bounding envelope for btilde computed'
               write(*,*)
     1              ' Efficiency for btilde generation is printed above'
               if((flg_withreg.or.flg_withdamp)
     1              .and..not.flg_bornonly) then
                  call loadmintupb(ndiminteg,'remnupb',ymaxrm,ymaxratrm)
                  write(*,*)
     1                 ' Upper bounding envelope for remnants computed'
                  write(*,*)
     1             ' Efficiency for remnant generation is printed above'
               endif
               ! CH, MK added or changed the following lines:
               !========================================================
#ifdef DSUB_II
               if(.not.flg_bornonly) then
                  call loadmintupb(ndiminteg,'osresupb',ymaxosres,
     &                             ymaxratosres)
                  write(*,*)
     1               ' Upper bounding envelope for osres-terms computed'
                  write(*,*)
     1             ' Efficiency for osres generation is printed above'
               endif
#endif
!                call storegrids(xgrid,ymax,ymaxrat,xgridrm,ymaxrm,
!      1              ymaxratrm,ifold,ifoldrm,-1,-1,'fullgrid')
               call storegrids(xgrid,ymax,ymaxrat,xgridrm,ymaxrm,
     1              ymaxratrm,xgridosres,ymaxosres,ymaxratosres,
     2              ifold,ifoldrm,ifoldosres,-1,-1,-1,-1,'fullgrid')
               !========================================================
            endif
         endif

c initialize gen; the array xmmm is set up at this stage.
         call gen(btilde,ndiminteg,xgrid,ymax,ymaxrat,xmmm,ifold,0,
     1    mcalls,icalls,xx)
         
         if (.not.flg_LOevents) then
c load or compute normalization of upper bounding function for radiation
c (iret ignored)
            call do_maxrat(mcalls,icalls,-1,iret)
         endif
c print statistics
         call gen(btilde,ndiminteg,xgrid,ymax,ymaxrat,xmmm,ifold,3,
     1        mcalls,icalls,xx)
         if(xx(1).gt.0)
     1        write(*,*) 'POWHEG: efficiency in the generation'
     2        //' of the Born variables =',xx(1)
         if((flg_withreg.or.flg_withdamp).and..not.flg_bornonly) then
c initialize gen for remnants
            call gen(sigremnant,ndiminteg,xgridrm,ymaxrm,ymaxratrm,
     1           xmmmrm,ifoldrm,0,mcalls,icalls,xx)
c     save random number seeds
            ! CH, MK: changed this to include regular terms, too
            !if(rad_totrm/rad_totgen.gt.1d-4.and.
            if((rad_totrem+rad_totreg)/rad_totgen.gt.1d-4.and.
     1           powheginput('#skipextratests').lt.0) then
               call randomsave
c     generate few events from remnants, just to determine the generation efficiency
               do j=1,int(min(powheginput('nubound'),10d0))
                  call gen(sigremnant,ndiminteg,xgridrm,ymaxrm,
     1                ymaxratrm,xmmmrm,ifoldrm,1,mcalls,icalls,xx)
               enddo
c     restore  random number seeds
               call randomrestore
c     print statistics
               call gen(sigremnant,ndiminteg,xgridrm,ymaxrm,ymaxratrm,
     1              xmmmrm,ifoldrm,3,mcalls,icalls,xx)
               if(xx(1).gt.0) then
                  write(*,*) 'POWHEG: efficiency in the generation'
     1                 //' of the remnant variables =',xx(1)
               endif
            endif
         endif
         ! CH, MK: added the following lines
         !==============================================================
         ! CH, MK: the osres-parts
#ifdef DSUB_II
         if(.not.flg_bornonly) then
            ! initialize gen for osress
            ! CH, MK: set the flg_btilde to false for osres-stuff
            flg_btilde=.false.
            call gen(sigosres,ndiminteg,xgridosres,ymaxosres,
     1           ymaxratosres,xmmmosres,ifoldosres,0,
     2           mcalls,icalls,xx)
            rad_totosres_sum = 0D0
            do ichan=1,nosres
              rad_totosres_sum = rad_totosres_sum + rad_totosres(ichan)
            enddo
            ! save random number seeds    
            if(dabs(rad_totosres_sum)/rad_totgen.gt.1d-4.and.
     1           powheginput('#skipextratests').lt.0) then
               call randomsave
               ! generate few events from osres-terms, just to determine
               ! the generation efficiency
               do j=1,int(min(powheginput('nubound'),10d0))
                  call gen(sigosres,ndiminteg,xgridosres,ymaxosres,
     1               ymaxratosres,xmmmosres,ifoldosres,1,
     2               mcalls,icalls,xx)
               enddo
               ! restore  random number seeds
               call randomrestore
               ! print statistics
               call gen(sigosres,ndiminteg,xgridosres,ymaxosres,
     1              ymaxratosres,xmmmosres,ifoldosres,3,
     2              mcalls,icalls,xx)
               if(xx(1).gt.0) then
                  write(*,*) 'POWHEG: efficiency in the generation'
     1                 //' of the osres variables =',xx(1)
               endif
            endif
            flg_btilde=.true.
         endif
#endif

         ! CH, MK: if the manyseeds-flag is set and we want to generate events using
         ! an old grid: simply set the random-number-generator to 
         ! the appropriate seed at this point
         if(rnd_cwhichseed.ne.'none')
     1     call setrandom(rnd_initialseed,rnd_i1,rnd_i2)
         !==============================================================
         return
      endif
      if(parallelstages.eq.1) then
         call bbinitxgrids(1)
         call pwhg_exit(0)
      endif
      call loadxgrids(iret)
      if(iret.ne.0) then
c if not there look for pwggridinfo-* files from parallel runs
         call loadlatestxgridinfo(iret)
         if(iret.ne.0) then
            write(*,*) ' cannot load xgrid or gridinfo files;'
            write(*,*) ' cannot perform stage 2'
            call pwhg_exit(-1)
         endif
      endif
c make sure the ifold arrays are read in
      do j=1,ndiminteg
         ifold(j)=1
         ifoldrm(j)=1
         ifoldosres(j)=1 ! CH, MK: added
      enddo
      ifold(ndiminteg-2) = powheginput("foldcsi")
      ifold(ndiminteg-1) = powheginput("foldy")
      ifold(ndiminteg)   = powheginput("foldphi")
      if(parallelstages.eq.2) then
         call bbinitgrids
         call writestat('st2')
         call pwhg_exit(0)
      endif
c in all cases pwggrids files must be loaded, to include information
c on total cross sections components (rad_tot* and rad_etot* variables)
      ! CH, MK: changed the following call
      call loadgrids(iret,xgrid,ymax,ymaxrat,xgridrm,ymaxrm,
     #           ymaxratrm,xgridosres,ymaxosres,ymaxratosres,
     #           ifold,ifoldrm,ifoldosres,'grid')
      if(iret.ne.0) then
         write(*,*) ' cannot load grid files'
         call pwhg_exit(-1)
      endif
      if(flg_storemintupb) then
c     Set up better upper bounding envelopes for btilde
         write(*,*) ' Loading the upper bounding envelope for btilde:'
         ! CH, MK: changed the following call
         call loadgrids(iret,xgrid,ymax,ymaxrat,xgridrm,ymaxrm,
     #           ymaxratrm,xgridosres,ymaxosres,ymaxratosres,
     #           ifold,ifoldrm,ifoldosres,'fullgrid')
         if(iret.ne.0) then
c this aborts if no files are found
            call loadmintupb(ndiminteg,'btildeupb',ymax,ymaxrat)
            write(*,*) ' Upper bounding envelope for btilde computed'
            write(*,*)
     1           ' Efficiency for btilde generation is printed above'
            if((flg_withreg.or.flg_withdamp).and..not.flg_bornonly) then
               call loadmintupb(ndiminteg,'remnupb',ymaxrm,ymaxratrm)
               write(*,*)
     1             ' Upper bounding envelope for remnants computed'
               write(*,*)
     1             ' Efficiency for remnant generation is printed above'
            endif
            ! CH, MK: added the following lines
            !===========================================================
#ifdef DSUB_II
            if(.not.flg_bornonly) then
               call loadmintupb(ndiminteg,'osresupb',ymaxosres,ymaxratosres)
               write(*,*)
     1             ' Upper bounding envelope for osres-temrs computed'
               write(*,*)
     1             ' Efficiency for osres generation is printed above'
            endif
#endif
!             call storegrids(xgrid,ymax,ymaxrat,xgridrm,ymaxrm,
!      1        ymaxratrm,ifold,ifoldrm,-1,-1,'fullgrid')
               call storegrids(xgrid,ymax,ymaxrat,xgridrm,ymaxrm,
     1              ymaxratrm,xgridosres,ymaxosres,ymaxratosres,
     2              ifold,ifoldrm,ifoldosres,-1,-1,-1,-1,'fullgrid')

            !===========================================================
         endif
      endif

c initialize gen; the array xmmm is set up at this stage.
      call gen(btilde,ndiminteg,xgrid,ymax,ymaxrat,xmmm,ifold,0,
     1     mcalls,icalls,xx)
      
      if((flg_withreg.or.flg_withdamp).and..not.flg_bornonly) then
c initialize gen for remnants
         call gen(sigremnant,ndiminteg,xgridrm,ymaxrm,ymaxratrm,
     1        xmmmrm,ifoldrm,0,mcalls,icalls,xx)
      endif

      ! CH, MK: added the following lines:
      !=================================================================
#ifdef DSUB_II
      if(.not.flg_bornonly) then
c initialize gen for osress
         call gen(sigosres,ndiminteg,xgridosres,ymaxosres,ymaxratosres,
     1        xmmmosres,ifoldosres,0,mcalls,icalls,xx)
      endif
#endif
      !=================================================================

      if(parallelstages.eq.3) then
         call writestat('st3')
c force compute upper bound for radiation (iret ignored)
         if(rnd_cwhichseed.ne.'none')
     1     call setrandom(rnd_initialseed,rnd_i1,rnd_i2)
         if (.not.flg_LOevents) then
            call do_maxrat(mcalls,icalls,0,iret)
         endif
         call pwhg_exit(0)
      endif

c force loading bound for radiation (iret ignored)
      iret=0
      if (.not.flg_LOevents) then
         call do_maxrat(mcalls,icalls,1,iret)
      endif
      if(iret.ne.0) then
         call pwhg_exit(0)
      endif
      if(rnd_cwhichseed.ne.'none')
     1     call setrandom(rnd_initialseed,rnd_i1,rnd_i2)
      end

      subroutine bbinitgrids
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_rnd.h'
      include 'pwhg_rad.h'
      include 'multigrid.h'
      include 'cgengrids.h' ! MK: changed position
c MK: added
#include "osres.h"
#include "pwhg_flst_add.h"
#include "Flags.h"
#include "cgengrids_add.h"
#include "pwhg_rad_add.h"
      integer iret1,iret2,iun
      ! CH, MK: changed
      !=================================================================
      real * 8 sigbtl,errbtl,sigrm,errrm,sigos,erros
      real * 8 btilde,sigremnant,sigosres,totneg
      integer ncall1,ncall1rm,ncall2,ncall2rm,ncall1osres,ncall2osres
      integer itmx1,itmx1rm,itmx1osres,itmx2,itmx2rm,itmx2osres
      !=================================================================
      real * 8 xx(ndiminteg)
      character * 40 mergelabels
      character * 40 filename
      character * 20 pwgprefix
      integer lprefix
      real*8 temparrosres(2,ntot_osres) ! CH, MK: added
      common/cpwgprefix/pwgprefix,lprefix
      integer j,k,mcalls,icalls,imode,iunstat
      logical savewithnegweights,multigrid
      real * 8 powheginput
      external btilde,sigremnant,powheginput,sigosres ! CH, MK: added
      double precision rad_totosresgen_sum, rad_etotosresgen_sum2 ! MK: added
      double precision rad_totosres_sum, rad_etotosres_sum2 ! MK: added
      double precision rad_totnegosres_sum ! MK: added
      integer ichan ! MK: added 
      call newunit(iunstat)
      open(unit=iunstat,file=mergelabels(pwgprefix,rnd_cwhichseed,
     1     'stat.dat',' '),status='unknown')
      if(rnd_cwhichseed.ne.'none')
     1     call setrandom(rnd_initialseed,rnd_i1,rnd_i2)
      if (powheginput('#testplots').eq.1d0) call init_hist
      if(flg_withnegweights) then
         write(*,*)' POWHEG: Computing pos.+|neg.| '
     1        //' weight contribution to inclusive cross section' 
      else
         write(*,*)' POWHEG: Computing positive weight'
     1        //' contribution to inclusive cross section' 
      endif
      ncall2=powheginput('ncall2')
      itmx2=powheginput('itmx2')
      ncall2rm=powheginput('#ncall2rm')
      if(ncall2rm.lt.0) ncall2rm=ncall2
      itmx2rm=powheginput('#itmx2rm')
      if(itmx2rm.lt.0) itmx2rm=itmx2
      ! CH, MK: added the following lines
      !=================================================================
      ncall2osres=powheginput('ncall2osres')
      if(ncall2osres.lt.0) ncall2osres=ncall2
      itmx2osres=powheginput('itmx2osres')
      if(itmx2osres.lt.0) itmx2osres=itmx2
      !=================================================================
      if(ncall2*itmx2.le.0) then
         write(*,*) 'ncall2=',ncall2,'itmx2=',itmx2
         write(*,*) 'Cannot continue'
         call pwhg_exit(0)
      endif
      flg_nlotest=.true.
      imode=1
c Totals will also be made available in the rad_tot???btl variables.
c The output in sigbtl is: positive weight only (flg_withnegweights=.false.)
c                          pos-|neg|            (flg_withnegweights=.true.)
c Results in rad_tot???btl do not depend upon flg_withnegweights
      call resettotals
      if(flg_storemintupb) call startstoremintupb('btildeupb')
      call setstage2init
      call mint(btilde,ndiminteg,ncall2,itmx2,ifold,imode,iun,
     1        xgrid,xint,xacc,nhits,ymax,ymaxrat,sigbtl,errbtl)
      if(flg_storemintupb) call stopstoremintupb
      call finaltotals
c finalize btilde output in histograms
      call pwhgaddout
      flg_nlotest=.false.
      write(*,*) 'btilde pos.   weights:', rad_totposbtl,' +-',
     1     rad_etotposbtl
      write(*,*) 'btilde |neg.| weights:', rad_totnegbtl,' +-',
     1     rad_etotnegbtl
      write(*,*) 'btilde total (pos.-|neg.|):', rad_totbtl,' +-',
     1     rad_etotbtl
      write(iunstat,*) 'btilde pos.   weights:', rad_totposbtl,' +-',
     1     rad_etotposbtl
      write(iunstat,*) 'btilde |neg.| weights:', rad_totnegbtl,' +-',
     1     rad_etotnegbtl
      write(iunstat,*) 'btilde Total (pos.-|neg.|):', rad_totbtl,
     1     ' +-',rad_etotbtl
c Now compute the remnant contributions
      call resettotalsrem ! CH: added
      if((flg_withreg.or.flg_withdamp).and..not.flg_bornonly) then
         write(*,*)' POWHEG: Computing remnant'//
     1        ' and/or regular remnants'
         flg_nlotest=.true.
         imode=1
         if(flg_storemintupb) call startstoremintupb('remnupb')
         call samegridasbtilde(ndiminteg,
     1        ncall2,itmx2,ifold,xgrid,
     1        ncall2rm,itmx2rm,ifoldrm,xgridrm)
         call mint(sigremnant,ndiminteg,ncall2rm,itmx2rm,ifoldrm,imode,
     1        iun,xgridrm,xintrm,xaccrm,nhitsrm,
     2        ymaxrm,ymaxratrm,sigrm,errrm)
         if(flg_storemintupb) call stopstoremintupb
c     add finalized remnant contributions in histograms
         call pwhgaddout
         flg_nlotest=.false.
         ! CH: commented
         !else
         !  sigrm=0
         !  errrm=0
      endif
      ! CH, MK: added or changed the following lines
      !=================================================================
      call finaltotalsrem ! CH: added, rad_totrm is deprecated from now on
      !rad_totrm=sigrm
      !rad_etotrm=errrm
      ! for the osres-parts: similar contribution as remnants
      call resettotalsosres ! CH: added
#ifdef DSUB_II
      if(.not.flg_bornonly) then
            write(*,*)' POWHEG: Computing ' ! CH, MK: modified the textoutput
            write(*,*) 'remaining terms from on-shell R,'//
     &                 ' incl. on-shell-subtraction'
            flg_nlotest=.true.
            imode=1
            ! set the flg_btilde to false for "osres remanants"
            flg_btilde=.false.
            if(flg_storemintupb) call startstoremintupb('osresupb')
            ! CH: changes: see integration of sigremnant
            call mint(sigosres,ndiminteg,ncall2osres,itmx2osres,
     1        ifoldosres,imode,iun,xgridosres,dabs(xintosres),xaccosres, ! MK: added dabs (without there occured many NaNs)
     2        nhitsosres,ymaxosres,ymaxratosres,sigos,erros)
            if(flg_storemintupb) call stopstoremintupb
c add finalized remnant contributions in histograms
            call pwhgaddout
            flg_nlotest=.false.
      endif
#endif
      call finaltotalsosres ! CH: added
      flg_btilde=.true.
      ! added: write-out the obtained results:
      if(flg_withdamp) then
        write(6,*) 'remnant total (pos.-|neg.|):', rad_totrem,
     &                                      ' +-', rad_etotrem
        write(iunstat,*) 'remnant total (pos.-|neg.|):', rad_totrem,
     &                                            ' +-', rad_etotrem
      endif
      if(flg_withreg) then
        write(6,*) 'regular total (pos.-|neg.|):', rad_totreg,
     &                                      ' +-', rad_etotreg
        write(iunstat,*) 'regular total (pos.-|neg.|):', rad_totreg,
     &                                            ' +-', rad_etotreg
      endif

#ifdef DSUB_II

      ! MK: combine the results of all on shell resonances:
      ! this is called only once and is mandatory
      call update_totarr ! MK: OK!
      do k=1,2
        do j=1,ntot_osres
          temparrosres(k,j)=0D0
          do ichan=1,nosres
            temparrosres(k,j)=temparrosres(k,j)+
     &                        rad_osres_totarr(k,j,ichan)
          enddo
        enddo
      enddo     
      write(6,*) 'real_osres pos.   weights:', 
     1        temparrosres(1,3),' +-',temparrosres(2,3)
      write(6,*) 'real_osres |neg.| weights:', 
     1        temparrosres(1,4),' +-',temparrosres(2,4)
      write(6,*) 'real_osres total (pos.-|neg.|):',
     1        temparrosres(1,1),' +-',temparrosres(2,1)
      write(iunstat,*) 'real_osres pos.   weights:', 
     1        temparrosres(1,3),' +-',temparrosres(2,3)
      write(iunstat,*) 'real_osres |neg.| weights:', 
     1        temparrosres(1,4),' +-',temparrosres(2,4)
      write(iunstat,*) 'real_osres total (pos.-|neg.|):', 
     1        temparrosres(1,1),' +-',temparrosres(2,1)
#ifdef DEBUGQ
      print*,"[DEBUG]: in bbinit_mod:596"
      print*,"temparrosres(k,j)",temparrosres
      print*,"rad_osres_totarr",rad_osres_totarr
      !stop
#endif 
#endif
      !=================================================================

c rad_totgen is used for the generation of the events.
c btilde and remnant event are chosen in proportion to
c rad_totbtlgen and rad_totrm.
c           if(flg_withnegweights) then
c              rad_totbtlgen=rad_totabsbtl
c              rad_etotbtlgen=rad_etotabsbtl
c           else
c     c notice: this is correct only if the negative fraction is
c     c negligible
c              rad_totbtlgen=rad_totbtl
c              rad_etotbtlgen=rad_etotbtl
c           endif
c           rad_totgen=rad_totrm+rad_totbtlgen
c           rad_etotgen=sqrt(rad_etotbtlgen**2+rad_etotrm**2)
c           rad_tot=rad_totrm+rad_totbtl
c           rad_etot=sqrt(rad_etotbtl**2+rad_etotrm**2)
c CH, MK: add analogues to totbtlgen (totremgen,totreggen,totosresgen(res)),
c which contain the remnants, the reg. and the osres-contributions 
c respectively...
      if(flg_withnegweights) then
         rad_totbtlgen=rad_totabsbtl
         rad_etotbtlgen=rad_etotabsbtl
         ! CH, MK: added the following lines
         !==============================================================
         rad_totreggen  = rad_totabsreg
         rad_etotreggen = rad_etotabsreg
         rad_totremgen  = rad_totabsrem
         rad_etotremgen = rad_etotabsrem
         do ichan=1,nosres
           rad_totosresgen(ichan)  = rad_totabsosres(ichan)
           rad_etotosresgen(ichan) = rad_etotabsosres(ichan)
         enddo
         !==============================================================
      else
c notice: this is correct only if the negative fraction is
c negligible
         rad_totbtlgen=rad_totbtl
         rad_etotbtlgen=rad_etotbtl
         ! CH, MK: added the following lines
         !==============================================================
         rad_totreggen  = rad_totreg
         rad_etotreggen = rad_etotreg
         rad_totremgen  = rad_totrem
         rad_etotremgen = rad_etotrem
         do ichan=1,nosres 
           rad_totosresgen(ichan)  = rad_totosres(ichan)
           rad_etotosresgen(ichan) = rad_etotosres(ichan)
         enddo
         !==============================================================
      endif

      ! CH, MK: added or changed the following lines
      !=================================================================
      rad_totosresgen_sum   = 0D0
      rad_etotosresgen_sum2 = 0D0
      rad_totosres_sum      = 0D0
      rad_etotosres_sum2    = 0D0
      rad_totnegosres_sum   = 0D0
      
      do ichan=1,nosres
        rad_totosresgen_sum   = rad_totosresgen_sum + 
     &                          rad_totosresgen(ichan)
        rad_etotosresgen_sum2 = rad_etotosresgen_sum2 + 
     &                          rad_etotosresgen(ichan)**2
        rad_totosres_sum      = rad_totosres_sum + 
     &                          rad_totosres(ichan)
        rad_etotosres_sum2    = rad_etotosres_sum2 +
     &                          rad_etotosres(ichan)**2
        rad_totnegosres_sum   = rad_totnegosres_sum + 
     &                          rad_totnegosres(ichan)
      enddo

#ifdef DEBUGQ
      print*,"[DEBUG]: in bbinit_mod:678"
      print*,"rad_totosresgen_sum",rad_totosresgen_sum
      print*,"rad_etotosresgen_sum2",rad_etotosresgen_sum2
      print*,"rad_totosres_sum",rad_totosres_sum
      print*,"rad_etotosres_sum2",rad_etotosres_sum2
      print*,"rad_totnegosres_sum",rad_totnegosres_sum
      !stop
#endif
      
      rad_totgen  = rad_totreggen + rad_totremgen + 
     &              rad_totosresgen_sum + rad_totbtlgen
      rad_etotgen = dsqrt(rad_etotbtlgen**2 + rad_etotreggen**2 +
     &                    rad_etotremgen**2 + rad_etotosresgen_sum2)

      ! physical xs
      rad_tot  = rad_totreg + rad_totrem + rad_totosres_sum + rad_totbtl
      rad_etot = dsqrt(rad_etotreg**2 + rad_etotrem**2 +
     &                 rad_etotosres_sum2 + rad_etotbtl**2)

      !=================================================================
      
c Grids are stored in all cases; they contain informations in
c rad_tot* and rad_etot* variables. The upper bound informations
c in pwggrid files is used only if mintupb files have not been saved.
c CH, MK: changed the subroutine call
c      call storegrids(xgrid,ymax,ymaxrat,
c    1 xgridrm,ymaxrm,ymaxratrm,ifold,ifoldrm,ncall2,itmx2,'grid')
      call storegrids(xgrid,ymax,ymaxrat,xgridrm,ymaxrm,
     1              ymaxratrm,xgridosres,ymaxosres,ymaxratosres,
     2              ifold,ifoldrm,ifoldosres,ncall2,itmx2,ncall2osres,
     3              itmx2osres,'grid')
c Output NLO histograms
      if (powheginput('#testplots').eq.1d0) then
         filename=mergelabels(pwgprefix,rnd_cwhichseed,'NLO',' ')
         call pwhgtopout(filename)
      endif
c CH, MK: already written out earlier...         
c      if(flg_withreg.or.flg_withdamp) then
c      write(iunstat,*) ' Remnant cross section in pb',
c     1        rad_totrm,'+-',rad_etotrm
c      endif
      ! CH, MK: modified the following statements
      !=================================================================
      if (flg_weightedev) then
         write(iunstat,*)
     1 ' total (btilde+remnants+regulars+osresR) cross section times '
         write(iunstat,*) ' suppression factor in pb',
     1        rad_tot,'+-',rad_etot         
      else
         write(iunstat,*)
     1 ' total (btilde+remnants+regulars+osresR) cross section in pb',
     2     rad_tot,'+-',rad_etot
      endif
      totneg=rad_totnegreg+rad_totnegrem+rad_totnegosres_sum+
     &          rad_totnegbtl ! CH, MK: added
      write(iunstat,*) ' negative weight fraction:',
     1              totneg/(2*totneg+rad_tot)
c     1     rad_totnegbtl/(2*rad_totnegbtl+rad_tot)
      
c      if(flg_withreg.or.flg_withdamp) then
c      write(*,*) ' Remnant cross section in pb',
c     1        rad_totrm,'+-',rad_etotrm
c      endif
      if (flg_weightedev) then
         write(*,*)
     1 ' total (btilde+remnants+regulars+osresR) cross section times '
         write(*,*) ' suppression factor in pb',
     1        rad_tot,'+-',rad_etot               
      else
         write(*,*)
     1 ' total (btilde+remnants+regulars+osresR) cross section in pb',
     2     rad_tot,'+-',rad_etot
      endif
      write(*,*) ' negative weight fraction:',
     1              totneg/(2*totneg+rad_tot)
c     1     rad_totnegbtl/(2*rad_totnegbtl+rad_tot)
      close(iunstat)
      !=================================================================
      end ! end bbinitgrids
      
      subroutine samegridasbtilde(ndiminteg,
     1        ncall2,itmx2,ifold,xgrid,
     1        ncall2rm,itmx2rm,ifoldrm,xgridrm)
c if the flag stage2init is set to 1 in powheg.input, the remnant
c cross section is computed with the same grid as the btilde cross
c section. In this way, the total cross section (times the suppression
c factor) computed with powheg (when using the same importance sampling
c grids) is identical whether withdamp is set or
c not, irrespective of the statistics. Useful for debugging.
      implicit none
      integer nintervals
      parameter (nintervals=50)
      integer ndiminteg,
     1        ncall2,itmx2,ifold(ndiminteg),
     1        ncall2rm,itmx2rm,ifoldrm(ndiminteg)
      real * 8 xgrid(0:nintervals,ndiminteg),
     1       xgridrm(0:nintervals,ndiminteg)
      real * 8 powheginput
      if(powheginput("#stage2init").eq.1d0) then
         ncall2rm = ncall2
         itmx2rm = itmx2
         ifoldrm = ifold
         xgridrm = xgrid
         call resetrandom
      endif
      end

      subroutine setstage2init
      implicit none
      real * 8 powheginput
      if(powheginput("#stage2init").eq.1d0) then
         call resetrandom
      endif
      end

c CH, MK: added stuff concerning osres parts
      subroutine bbinitxgrids(iparallel)
c iparallel = 0: generate xgrid file
c iparallel = 1: generate gridinfo files for parallel grid      
      implicit none
      character * 6 tag,prevtag
      integer iparallel,itmp
      ! CH, MK: added   
      integer iteration,imode,j,itmx1,ncall1,ncall1rm,itmx1rm
      integer ncall1osres,itmx1osres,iun,iret
      integer btilde,sigremnant,sigosres
      include 'nlegborn.h'
      include 'pwhg_flg.h'
      include 'pwhg_rnd.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'pwhg_par.h'
      include 'cgengrids.h'
c MK: added
#include "osres.h"
#include "pwhg_flst_add.h"
#include "Flags.h"
#include "cgengrids_add.h"
#include "pwhg_rad_add.h"
      character * 20 pwgprefix
      character * 40 mergelabels
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      real * 8 sigbtl,errbtl,sigrm,errrm,sigos,erros ! CH, MK: added
      real * 8 random,powheginput
      logical savewithnegweights
      external btilde,sigremnant,random,powheginput,sigosres ! CH, MK: added
c No folding while generating importance sampling grids
      do j=1,ndiminteg
         ifold(j)=1
         ifoldrm(j)=1
         ifoldosres(j)=1 ! CH, MK: added
      enddo
c         
      if(iparallel.eq.1) then
         iteration = powheginput('#xgriditeration')
         if(iteration.gt.par_maxxgriditerations) then
            write(*,*) ' POWHEG is compiled with a maximum'
            write(*,*) ' number of x-grid iterations=',
     1           par_maxxgriditerations
            write(*,*) ' increase the par_maxxgriditerations parameter'
            write(*,*) ' in the pwhg_par.h file to use more.'
            write(*,*) ' (BUT YOU SHOULD NOT NEED MORE THAN 3 or 4!'
            write(*,*) ' increase ncall1 instead)'
            call pwhg_exit(-1)
         endif
         if(iteration.le.0) iteration = 1
         write(tag,'(i3)') iteration
         tag=adjustl(tag)
         tag='xg'//tag
         if(iteration.gt.1) then
            write(prevtag,'(i3)') iteration-1
            prevtag=adjustl(prevtag)
            prevtag='xg'//prevtag
         endif
      else
         iteration = 0
         tag=' '
         prevtag=' '
      endif
c In this block we compute the importance sampling grid
      write(*,*)
      write(*,*)' POWHEG: initialization'
      write(*,*)' Computing the integral of the absolute value'
      if(flg_weightedev) then
         write(*,*)' of the cross section times the suppression' 
         write(*,*)' factor to set up the adaptive grid'
      else
         write(*,*)' of the cross section to set up the adaptive grid'
      endif
c If parallel operations, set a different random number also for
c different grid iterations
      if(iparallel.eq.1) then
         itmp=0
         call resetrandom
         do j=1,iteration
            itmp = random()*1d9
         enddo
         call setrandom(rnd_initialseed+itmp,rnd_i1,rnd_i2)
      endif
      if(iteration.le.1) then
         call initxgrid(xgrid,ndiminteg)
      else
         call loadgridinfo(mergelabels('btl',prevtag,' ',' '),.false.,
     1        xgrid,xint,iret)
         if(iret.ne.0) then
            write(*,*) 'cannot find level '//prevtag//'gridinfo files'
            call pwhg_exit()
         endif
      endif
c The actual value of the grid is the one to be saved in the
c gridinfo file, since loadgridinfo updates the grid by itself
      xgrid0 = xgrid
      ncall1 = powheginput("ncall1")
      ncall1rm = powheginput("#ncall1rm")
      ! ncall1rm is additional, if it is not present use the standard one:
      if (ncall1rm.lt.0d0) ncall1rm = ncall1
      ! CH, MK: added
      !=================================================================
      ncall1osres = powheginput("ncall1osres")
      ! ncall1osres is additional, if it is not present use the standard one:
      if (ncall1osres.lt.0d0) ncall1osres = ncall1
      !=================================================================
      itmx1 = powheginput("itmx1")
      itmx1rm = powheginput("#itmx1rm")
      if(itmx1rm.lt.0) itmx1rm=itmx1
      ! CH, MK: added
      !=================================================================
      itmx1osres = powheginput("itmx1osres")
      if(itmx1osres.lt.0) itmx1osres=itmx1
      !=================================================================
c with parallel grids only one iteration is allowed
      if(iparallel.eq.1) then
         itmx1 = 1
         itmx1rm = 1
         itmx1osres = 1
      endif
      call newunit(iun)
      call regridplotopen(mergelabels(pwgprefix(1:lprefix),tag,
     1     rnd_cwhichseed,'btlgrid.top'))
      write(*,*)' result +- errtot (picobarn) for each iteration'
      flg_nlotest=.false.
      imode=0
      savewithnegweights=flg_withnegweights
      flg_withnegweights=.true.
c ********** CALL to mint for btilde
      call mint(btilde,ndiminteg,ncall1,itmx1,ifold,imode,iun,
     1     xgrid,xint,xacc,nhits,ymax,ymaxrat,sigbtl,errbtl)

c **********
      call regridplotclose
      if(iteration.ge.1) then
         call storegridinfo('btl-'//trim(tag),xgrid0,xint,
     1        xacc,nhits,ndiminteg)
      endif
      flg_withnegweights=savewithnegweights
      close(iun)
      if((flg_withreg.or.flg_withdamp).and..not.flg_bornonly) then
         write(*,*) ' Computing the integral of the'//
     1        ' remnant cross section'
         if(flg_weightedev) then
            write(*,*) 'times the suppression factor'
         endif
         write(*,*) ' to set up the adaptive grid'
         flg_nlotest=.false.
         imode=0
         if(iteration.le.1) then
            call initxgrid(xgridrm,ndiminteg)
         else
            call loadgridinfo(mergelabels('rmn',prevtag,' ',' '),
     1           .false.,xgridrm,xintrm,iret)
         endif
         call newunit(iun)
         call regridplotopen(mergelabels(pwgprefix(1:lprefix),tag,
     1     rnd_cwhichseed,'rmngrid.top'))
         xgrid0rm=xgridrm
c ********** CALL to mint for remnants
         call mint(sigremnant,ndiminteg,ncall1rm,itmx1rm,ifoldrm,imode,
     1        iun,xgridrm,xintrm,xaccrm,nhitsrm,
     1        ymaxrm,ymaxratrm,sigrm,errrm)

c **********
         call regridplotclose
         if(iteration.ge.1) then
            call storegridinfo('rmn-'//trim(tag),xgrid0rm,xintrm,
     1           xaccrm,nhitsrm,ndiminteg)
         endif
      endif
      ! CH, MK: added the osres-part
      !=================================================================
#ifdef DSUB_II
      if(.not.flg_bornonly) then
          flg_nlotest=.false.
          imode=0
          if(iteration.le.1) then
            call initxgrid(xgridosres,ndiminteg)
          else
            call loadgridinfo(mergelabels('osres',prevtag,' ',' '),
     1           .false.,xgridosres,xintosres,iret)
         endif
         call newunit(iun)
         call regridplotopen(mergelabels(pwgprefix(1:lprefix),tag,
     1     rnd_cwhichseed,'osresgrid.top'))
         xgrid0osres=xgridosres
c ********** CALL to mint for osres remnants
         ! CH, MK: set the flg_btilde to false for "osres remnants"
         flg_btilde=.false.
         call mint(sigosres,ndiminteg,ncall1osres,itmx1osres,ifoldosres,imode,
     1        iun,xgridosres,xintosres,xaccosres,nhitsosres,
     1        ymaxosres,ymaxratosres,sigos,erros)
         flg_btilde=.true.
c **********
         call regridplotclose
         if(iteration.ge.1) then
            call storegridinfo('osres-'//trim(tag),xgrid0osres,xintosres,
     1           xaccosres,nhitsosres,ndiminteg)
         endif
      endif
#endif
      !=================================================================
      if(iparallel.eq.0) then
         !call storexgrid(xgrid,xint,xgridrm,xintrm)
         call storexgrid(xgrid,xint,xgridrm,xintrm,xgridosres,xintosres) ! CH: changed
c      else
c         call newunit(iun)
c         open(unit=iun,
c     1        file=mergelabels('stage-xiterations',tag,'.dat',' '),
c     1        status='unknown')
c         write(iun,*)'done'
c         close(iun)
      endif
      write(*,*)' Importance sampling x grids generated and stored'
      end

      subroutine pwhg_openoutput(iun,string,suffix)
      implicit none
      integer iun
      character * (*) string
      character * (*) suffix
      include 'pwhg_rnd.h'
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      character * 40 mergelabels
      open(unit=iun,file=mergelabels(pwgprefix(1:lprefix)
     1     //trim(string),rnd_cwhichseed,suffix,' '),status='unknown')
      end


      subroutine gen_btilde(mcalls,icalls)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer mcalls,icalls
      include 'cgengrids.h'
c MK: added
#include "cgengrids_add.h"
#include "pwhg_flst_add.h"
      real * 8 xx(ndiminteg)      
      real * 8 btilde
      external btilde
      call gen(btilde,ndiminteg,xgrid,ymax,ymaxrat,xmmm,ifold,1,
     1     mcalls,icalls,xx)
      end

      subroutine gen_sigremnant
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'cgengrids.h'
c MK: added
#include "pwhg_flst_add.h"
#include "Flags.h"
#include "cgengrids_add.h"
      real * 8 xx(ndiminteg)
      integer mcalls,icalls
      logical savelogical
      real * 8 sigremnant
      external sigremnant
c communicate file to load upper bound data
      savelogical=flg_fastbtlbound
      flg_fastbtlbound=.false.
      call gen(sigremnant,ndiminteg,xgridrm,ymaxrm,ymaxratrm,
     1    xmmmrm,ifoldrm,1,mcalls,icalls,xx)
      flg_fastbtlbound=savelogical
      end


c CH, MK: added this routine, in principle a copy of gen_sigremnant
      subroutine gen_sigosres
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'cgengrids.h'
c MK: added
#include "pwhg_flst_add.h"
#include "Flags.h"
#include "cgengrids_add.h"
      real * 8 xx(ndiminteg)
      integer mcalls,icalls
      logical savelogical
      real * 8 sigosres
      external sigosres
c communicate file to load upper bound data
      savelogical=flg_fastbtlbound
      flg_fastbtlbound=.false.
      call gen(sigosres,ndiminteg,xgridosres,ymaxosres,ymaxratosres,
     &    xmmmosres,ifoldosres,1,mcalls,icalls,xx)
      flg_fastbtlbound=savelogical
      end


      ! CH, MK: changed the routine arguments
      subroutine storegrids(xgrid,ymax,ymaxrat,xgridrm,ymaxrm,
     1                ymaxratrm,xgridosres,ymaxosres,ymaxratosres,
     2                ifold,ifoldrm,ifoldosres,ncall2,itmx2,ncall2osres,
     3                itmx2osres,gridtag)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_pdf.h'
      include 'pwhg_rad.h'
      include 'pwhg_rnd.h'
c MK: added
#include "osres.h"
#include "pwhg_rad_add.h"
#include "pwhg_flst_add.h"
      real * 8 xgrid(0:50,ndiminteg),ymax(50,ndiminteg),
     1     ymaxrat(50,ndiminteg),xgridrm(0:50,ndiminteg),
     2     ymaxrm(50,ndiminteg),ymaxratrm(50,ndiminteg),
     3     xgridosres(0:50,ndiminteg),ymaxosres(50,ndiminteg), ! CH, MK: added
     4     ymaxratosres(50,ndiminteg) ! CH, MK: added
      integer nbins
      parameter (nbins=50)
      integer ifold(ndiminteg),ifoldrm(ndiminteg),ifoldosres(ndiminteg), ! CH, MK: added osres
     #        ncall2,itmx2,ncall2osres,itmx2osres ! CH, MK: added _reg
      integer ichan ! MK: added
      character *(*) gridtag
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      integer j,k,iun
      call newunit(iun)
      if(rnd_cwhichseed.eq.'none') then
         open(unit=iun,file=pwgprefix(1:lprefix)//gridtag//'.dat',
     #     form='unformatted',status='unknown')
      else
         open(unit=iun,
     1        file=pwgprefix(1:lprefix)//gridtag//'-'//rnd_cwhichseed//
     2        '.dat',form='unformatted',status='unknown')
      endif
      write(iun) ((xgrid(j,k),k=1,ndiminteg),j=0,nbins)
      write(iun) ((ymax(j,k),k=1,ndiminteg),j=1,nbins)
      write(iun) ((ymaxrat(j,k),k=1,ndiminteg),j=1,nbins)
      write(iun) ((xgridrm(j,k),k=1,ndiminteg),j=0,nbins)
      write(iun) ((ymaxrm(j,k),k=1,ndiminteg),j=1,nbins)
      write(iun) ((ymaxratrm(j,k),k=1,ndiminteg),j=1,nbins)
      ! CH, MK: added
      !=================================================================
      write(iun) ((xgridosres(j,k),k=1,ndiminteg),j=0,nbins)
      write(iun) ((ymaxosres(j,k),k=1,ndiminteg),j=1,nbins)
      write(iun) ((ymaxratosres(j,k),k=1,ndiminteg),j=1,nbins)
      !=================================================================
      write(iun) (ifold(k),k=1,ndiminteg)
      write(iun) (ifoldrm(k),k=1,ndiminteg)
      write(iun) (ifoldosres(k),k=1,ndiminteg) ! CH, MK: added
      write(iun) ncall2*itmx2
      write(iun) ncall2osres*itmx2osres ! CH, MK: added
      write(iun) kn_sbeams, pdf_ih1, pdf_ih2, pdf_ndns1, pdf_ndns2
      ! MK: rad_totarr contains all the different parts,
      ! via the subroutines update_totarr and apply_totarr they are
      ! related tot rad_totbtl, ...
      !=================================================================
      ! this is called only once and is mandatory
      call update_totarr ! MK: OK!
      write(iun) ((rad_totarr(j,k),j=1,2),k=1,ntot)
      write(iun) (((rad_osres_totarr(j,k,ichan),j=1,2),k=1,ntot_osres)
     &                                              ,ichan=1,nosres)
      !=================================================================
      close(iun)
      end

      ! CH, MK: changed routine arguments
      subroutine loadgrids(iret,xgrid,ymax,ymaxrat,xgridrm,ymaxrm,
     #           ymaxratrm,xgridosres,ymaxosres,ymaxratosres,
     #           ifold,ifoldrm,ifoldosres,gridtag)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_pdf.h'
      include 'pwhg_rad.h'
      include 'pwhg_rnd.h'
      include 'pwhg_par.h'
c MK: added
#include "osres.h"
#include "pwhg_rad_add.h"
#include "pwhg_flst_add.h"
      real * 8 xgrid(0:50,ndiminteg),ymax(50,ndiminteg),
     1     ymaxrat(50,ndiminteg),xgridrm(0:50,ndiminteg),
     2     ymaxrm(50,ndiminteg),ymaxratrm(50,ndiminteg)
      real * 8 xxgrid(0:50,ndiminteg),xymax(50,ndiminteg),
     1     xymaxrat(50,ndiminteg),xxgridrm(0:50,ndiminteg),
     2     xymaxrm(50,ndiminteg),xymaxratrm(50,ndiminteg)
      ! CH, MK: added or modified the following lines
      !=================================================================
      real * 8 xgridosres(0:50,ndiminteg),ymaxosres(50,ndiminteg),
     1     ymaxratosres(50,ndiminteg),xxgridosres(0:50,ndiminteg),
     2     xymaxosres(50,ndiminteg),xymaxratosres(50,ndiminteg)
!          real * 8 tot(2,8),rtot(2,8) ! CH: replaced with following lines
      real*8 tot(2,ntot) ! CH, MK: contains the rad_totarr-entries
      real*8 tot_osres(2,ntot_osres,cnosres) ! MK: contains the rad_osres_totarr-entries
      integer ichan ! MK: added to distinguish the on-shell resonances
      integer ifold(ndiminteg),ifoldrm(ndiminteg),ifoldosres(ndiminteg)
      integer iifold(ndiminteg),iifoldrm(ndiminteg),iifoldosres(ndiminteg)
      integer iret,iretcheck,jfound
      !=================================================================
      character *(*) gridtag
c
      integer ios
      integer nbins
      parameter (nbins=50)
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      character * 4 chseed, firstfound
      real * 8 shx
      integer ih1x, ih2x, ndns1x, ndns2x
      integer j,k,iun,jfile,nfiles,ncall2,itmx2,ncall2osres,itmx2osres ! CH, MK: added _reg
      logical lpresent,manyfiles,filefound
      double precision, allocatable :: totarr(:,:)
      logical, allocatable :: goodentries(:)
      logical firsttime
      logical, save :: check_bad_st2
      real * 8 rjfound, rncall2
      real * 8 powheginput
      external powheginput
      if(powheginput('use-old-grid').eq.0) then
         iret=1
         return
      endif
      iret=0
      call newunit(iun)
      open(unit=iun,file=pwgprefix(1:lprefix)//gridtag//'.dat',
     #     form='unformatted',status='old',iostat=ios)
      if(ios.eq.0) then
         nfiles=1
         check_bad_st2 = .false.
      else
         nfiles=par_maxseeds
         manyfiles=.true.
         check_bad_st2 =
     1        powheginput("#check_bad_st2") == 1
         if(check_bad_st2) then
            allocate(totarr(2,nfiles),goodentries(nfiles))
            goodentries = .true.
            totarr = 0
         endif
      endif
c Try to open and merge a set of grid files, generated with different
c random seeds
      firsttime = .true.
 11   continue
      filefound=.false.
      jfound=0
      do jfile=1,nfiles
         if(manyfiles) then
            if(check_bad_st2) then
               if(.not.goodentries(jfile)) cycle
            endif
            write(chseed,'(i4)') jfile
            do k=1,4
               if(chseed(k:k).eq.' ') chseed(k:k)='0'
            enddo
            inquire(file=pwgprefix(1:lprefix)//gridtag//'-'//
     1           chseed//'.dat',exist=lpresent)
            if(.not.lpresent) goto 111
            open(unit=iun,file=pwgprefix(1:lprefix)//gridtag//'-'//
     1           chseed//'.dat',
     2           form='unformatted',status='old',iostat=ios)
            if(ios.ne.0) then
               iret=-1
               return
            else
               write(*,*)
     1              ' Opened ',pwgprefix(1:lprefix)//gridtag//'-'//
     2              chseed//'.dat'
            endif
         endif
         filefound=.true.
         read(iun,iostat=ios) ((xxgrid(j,k),k=1,ndiminteg),j=0,nbins)
         if(ios.ne.0) goto 998
         read(iun,iostat=ios) ((xymax(j,k),k=1,ndiminteg),j=1,nbins)
         if(ios.ne.0) goto 998
         read(iun,iostat=ios) ((xymaxrat(j,k),k=1,ndiminteg),j=1,nbins)
         if(ios.ne.0) goto 998
         read(iun,iostat=ios) ((xxgridrm(j,k),k=1,ndiminteg),j=0,nbins)
         if(ios.ne.0) goto 998
         read(iun,iostat=ios) ((xymaxrm(j,k),k=1,ndiminteg),j=1,nbins)
         if(ios.ne.0) goto 998
         read(iun,iostat=ios)((xymaxratrm(j,k),k=1,ndiminteg),j=1,nbins)
         if(ios.ne.0) goto 998
         ! CH, MK: added the following lines
         !==============================================================
         read(iun,iostat=ios) 
     &                   ((xxgridosres(j,k),k=1,ndiminteg),j=0,nbins)
         if(ios.ne.0) goto 998
         read(iun,iostat=ios) 
     &                   ((xymaxosres(j,k),k=1,ndiminteg),j=1,nbins)
         if(ios.ne.0) goto 998    
         read(iun,iostat=ios)  
     &                   ((xymaxratosres(j,k),k=1,ndiminteg),j=1,nbins)
         if(ios.ne.0) goto 998  
         !==============================================================
         read(iun,iostat=ios) (iifold(k),k=1,ndiminteg)
         if(ios.ne.0) goto 998
         read(iun,iostat=ios) (iifoldrm(k),k=1,ndiminteg)
         if(ios.ne.0) goto 998
         ! CH, MK: added or modified the following lines
         !==============================================================
         read(iun,iostat=ios) (iifoldosres(k),k=1,ndiminteg)
         if(ios.ne.0) goto 998
         read(iun,iostat=ios) ncall2
         read(iun,iostat=ios) ncall2osres ! CH, MK: added
         if(powheginput("#ncallfrominput").eq.1) then
            ncall2=powheginput("ncall2")
            itmx2=powheginput("itmx2")
            ncall2=ncall2*itmx2
            ! CH, MK: read-in the different statistic-parameters 
            ! for the regulars (if present!) ->
            ncall2osres=powheginput("#ncall2osres")
            itmx2osres=powheginput("#itmx2osres")
            if (ncall2osres.lt.0) ncall2osres=ncall2
            if (itmx2osres.lt.0) itmx2osres=itmx2
            ncall2osres=ncall2osres*itmx2osres
         !==============================================================
         endif
         if(ios.ne.0) goto 998
         read(iun,iostat=ios) shx, ih1x, ih2x, ndns1x, ndns2x
         if(ios.ne.0) goto 998
         if(shx.ne.kn_sbeams.or.ih1x.ne.pdf_ih1.or.ih2x.ne.pdf_ih2
     1      .or.ios.ne.0)
     2        goto 998
!          read(iun,iostat=ios) ((tot(k,j),k=1,2),j=1,8)
         read(iun,iostat=ios) ((tot(k,j),k=1,2),j=1,ntot) ! CH, MK: changed
         read(iun,iostat=ios) (((tot_osres(k,j,ichan),k=1,2),
     &                           j=1,ntot_osres),ichan=1,nosres) ! MK: added
         if(ios.ne.0) goto 998
         jfound=jfound+1
         if(jfound.lt.2) then
            firstfound = chseed
            do k=1,ndiminteg
               do j=0,nbins
                  xgrid(j,k)=xxgrid(j,k)
                  xgridrm(j,k)=xxgridrm(j,k)
                  xgridosres(j,k)=xxgridosres(j,k) ! CH, MK: added
               enddo
               ifold(k)=iifold(k)
               ifoldrm(k)=iifoldrm(k)
               ifoldosres(k)=iifoldosres(k) ! CH, MK: added
            enddo
            do k=1,ndiminteg
               do j=1,nbins
                  ymax(j,k)=xymax(j,k)
                  ymaxrm(j,k)=xymaxrm(j,k)
                  ymaxosres(j,k)=xymaxosres(j,k) ! CH, MK: added
                  ymaxrat(j,k)=xymaxrat(j,k)
                  ymaxratrm(j,k)=xymaxratrm(j,k)
                  ymaxratosres(j,k)=xymaxratosres(j,k) ! CH, MK: added
               enddo
            enddo
            do k=1,2
               ! CH, MK: if we have only one grid: simply save the results
               ! contained in it
               do j=1,ntot
                   rad_totarr(k,j)=tot(k,j)
               enddo
               ! MK: added the following statements
               do j=1,ntot_osres
                 do ichan=1,nosres
                   rad_osres_totarr(k,j,ichan)=tot_osres(k,j,ichan)
                 enddo
               enddo
            enddo
            ! this is called only once and is mandatory
            call apply_totarr ! MK: save the array totarr to the individual variables -> OK!
         else
            do k=1,ndiminteg
               do j=0,nbins
                  if(xgrid(j,k).ne.xxgrid(j,k).or.
     1                 xgridrm(j,k).ne.xxgridrm(j,k).or.
     2                 xgridosres(j,k).ne.xxgridosres(j,k)) then ! CH, MK: added
                     write(*,*) ' error loading grids: '
                     write(*,*)  pwgprefix(1:lprefix)//gridtag//'-'//
     1          rnd_cwhichseed//'.dat does not have the same importance'
                    write(*,*) 'sampling grid as ',pwgprefix(1:lprefix)
     1              //gridtag//'-'//firstfound//'.dat'
                     call pwhg_exit(-1)
                  endif
               enddo
               if(ifold(k).ne.iifold(k)
     1              .or.ifoldrm(k).ne.iifoldrm(k).or.
     2                  ifoldosres(k).ne.ifoldosres(k)) then ! CH, MK: added
                  write(*,*) ' error loading grids: '
                  write(*,*)  pwgprefix(1:lprefix)//gridtag//'-'//
     1                 rnd_cwhichseed//
     2                 '.dat does not have the same folding as '
                  write(*,*) pwgprefix(1:lprefix)//gridtag//'-'
     1                 //firstfound//'.dat'
                  call pwhg_exit(-1)
               endif
            enddo
            do k=1,ndiminteg
               do j=1,nbins
                  ymax(j,k)=max(ymax(j,k),xymax(j,k))
                  ymaxrm(j,k)=max(ymaxrm(j,k),xymaxrm(j,k))
                  ymaxosres(j,k)=max(ymaxosres(j,k),xymaxosres(j,k)) ! CH, MK: added
                  ymaxrat(j,k)=max(ymaxrat(j,k),xymaxrat(j,k))
                  ymaxratrm(j,k)=max(ymaxratrm(j,k),xymaxratrm(j,k))
                  ymaxratosres(j,k)=max(ymaxratosres(j,k),
     &                                  xymaxratosres(j,k)) ! CH, MK: added
               enddo
            enddo
c turn these to reals; very large integer can overflow in fortran
c               rjfound = jfound
c               rncall2 = ncall2
c               rtot(2,j)=sqrt((rtot(2,j)**2*(rjfound-1)**2+tot(2,j)**2)
c     1              /rjfound**2+(rjfound-1)*(rtot(1,j)-tot(1,j))**2/
c     2              (rjfound**3*rncall2))
c               rtot(1,j)=(rtot(1,j)*(rjfound-1)+tot(1,j))/rjfound
c CH: here we have to take the average of all entries found in the grid file...
c again assign these values directly to the rad_totarr
c if we use different VEGAS-parameters for the osres-terms: make sure to take the
c appropriate ncall2 here (this is not possible for rad_tot(gen), as here btilde and
c osres. parts are naturally "mixed"-> simply recalculate these 2 entries at the end
            do j=1,ntot
              if(j.le.15) then! CH, MK: these are the btilde/remnant-entries
                rad_totarr(2,j)=dsqrt((rad_totarr(2,j)**2*(jfound-1)**2+
     &                          tot(2,j)**2)/jfound**2+(jfound-1)*
     &                          (rad_totarr(1,j)-tot(1,j))**2/
     &                          (jfound**3*ncall2))
              else
                rad_totarr(2,j)=dsqrt((rad_totarr(2,j)**2*(jfound-1)**2+
     &                          tot(2,j)**2)/jfound**2+(jfound-1)*
     &                          (rad_totarr(1,j)-tot(1,j))**2/
     &                          (jfound**3*ncall2osres))
              endif
              rad_totarr(1,j)=(rad_totarr(1,j)*(jfound-1)+tot(1,j))/
     &                         jfound
            enddo

            ! MK: added the on-shell contributions
            do j=1,ntot_osres
              do ichan=1,nosres
                rad_osres_totarr(2,j,ichan)=
     &            dsqrt((rad_osres_totarr(2,j,ichan)**2*(jfound-1)**2+
     &            tot_osres(2,j,ichan)**2)/jfound**2+(jfound-1)*
     &            (rad_osres_totarr(1,j,ichan)-tot_osres(1,j,ichan))**2/
     &            (jfound**3*ncall2osres))
                rad_osres_totarr(1,j,ichan)=
     &            (rad_osres_totarr(1,j,ichan)*(jfound-1)+
     &            tot_osres(1,j,ichan))/jfound
              enddo
            enddo
            call apply_totarr ! MK: added -> OK!
         endif
         close(iun)
 111     continue
      enddo
      if(filefound) then
         if(manyfiles .and. check_bad_st2) then
c check that the different runs are more or less consistent
            if(firsttime) then
               call check_stat_consistency(nfiles,totarr,
     1              goodentries,iretcheck)
               if(iretcheck.eq.-1) then
                  firsttime = .false.
                  goto 11
               endif
            endif
            deallocate(goodentries,totarr)
         endif
         return
      endif
 998  continue
      iret=-1
      end

      subroutine  check_stat_consistency(nentries,res,goodentries,iret)
      implicit none
      integer nentries,iret
      integer indices(nentries)
      double precision res(2,nentries)
      logical goodentries(nentries)
      double precision average,weight,tmpav,tmpweight
      double precision tmp2(2),ow,oav
      integer nonz,j,k,itmp
      logical :: ex=.true.
c DEBUG
c      res(1,7) = 3.7
c      res(2,7) = 1
c      res(1,5) = 3.7
c      res(2,5) = 1
c      res(1,3) = 3.7
c      res(2,3) = 1
c end DEBUG
      do j=1,nentries
         indices(j)=j
      enddo
c bubblesort
      do while(ex)
         ex = .false.
         do j=1,nentries-1
c swap in growing order, but put zeros at the end
            if((res(2,j+1).ne.0.and.res(2,j+1).lt.res(2,j)).or.
     1        (res(2,j+1).ne.0.and.res(2,j).eq.0)) then
               tmp2 = res(:,j)
               res(:,j) = res(:,j+1)
               res(:,j+1) = tmp2
               itmp = indices(j)
               indices(j) = indices(j+1)
               indices(j+1) = itmp
               ex = .true.
            endif
         enddo
      enddo
      do j=1,nentries
         if(res(2,j).eq.0) then
            nonz = j - 1
            exit
         endif
      enddo
c Compute the average. Neglect zero entries.
      average = 0
      weight = 0
      do j=1,nonz
         if(res(1,j).ne.0) then
            oav = average
            ow = weight
            average = (average*(j-1) + res(1,j))/j
            weight = sqrt((weight*(j-1))**2 + res(2,j)**2)/j
            if(j.gt.1) then
c               write(*,*) ' old,new av.',oav,average
c               write(*,*) ' old,new err.',ow,weight
c               write(*,*) ' deviation:',abs(oav-average)/ow
c after half of the runs
               if(abs(oav-average)/ow.gt.10) then
                  write(*,*) ' check_stat_consistency:'
                  write(*,*)
     1                 ' The program has detected inconsistent results'
                  write(*,*) ' among different runs. The runs:'
                  write(*,*) (indices(k),k=j,nonz)
                  write(*,*) ' look suspicious. '
                  write(*,*) ' Inspect your runs at stage 2 '
                  if(nonz-j+1.gt. nonz/10) then
                     write(*,*)
     1               ' The fraction of inconsistent runs is too large'
                     write(*,*) ' exiting ...'
                     call exit(-1)
                  else
                     write(*,*)
     1               ' The fraction of inconsistent file is < 10%'
                     write(*,*) ' We discard the following files:'
                     write(*,*) (indices(k),k=j,nonz)
                     write(*,*) ' and reload the others'                     
                     do k=j,nonz
                        goodentries(indices(k)) = .false.
                     enddo
                     iret = -1
                     return
                  endif
               endif
            endif
         endif
      enddo
      iret = 0
      end

c CH, MK: changed routine arguments
      subroutine storexgrid(xgrid,xint,xgridrm,xintrm,xgridosres,xintosres)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_pdf.h'
      include 'pwhg_rad.h'
c MK: added
#include "osres.h"
#include "pwhg_rad_add.h"
#include "pwhg_flst_add.h"
      real * 8 xgrid(0:50,ndiminteg),xgridrm(0:50,ndiminteg),
     1     xint,xintrm,xgridosres(0:50,ndiminteg),xintosres ! CH, MK: added
      integer nbins
      parameter (nbins=50)
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      integer j,k,iun
      call newunit(iun)
      open(unit=iun,file=pwgprefix(1:lprefix)//'xgrid.dat',
     #     form='unformatted',status='unknown')
      write(iun) ((xgrid(j,k),k=1,ndiminteg),j=0,nbins),xint
      write(iun) ((xgridrm(j,k),k=1,ndiminteg),j=0,nbins),xintrm
      write(iun) ((xgridosres(j,k),k=1,ndiminteg),j=0,nbins),xintosres ! CH, MK: added
      write(iun) kn_sbeams, pdf_ih1, pdf_ih2, pdf_ndns1, pdf_ndns2
      close(iun)
      end

      subroutine loadxgrids(iret)
c loads an integration grid from file <prefix>xgrid.dat, if found
c (returns iret=0). If not found iret=-1. It checks that some key
c parameters are the same as for the current run
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_pdf.h'
      include 'pwhg_rad.h'
      include 'cgengrids.h'
c MK: added
#include "osres.h"
#include "cgengrids_add.h"
#include "pwhg_rad_add.h"
#include "pwhg_flst_add.h"
      integer iret
c
      integer ios
      integer nbins
      parameter (nbins=50)
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      real * 8 shx
      integer ih1x, ih2x, ndns1x, ndns2x
      integer j,k,iun
      real * 8 powheginput
      external powheginput
      if(powheginput('use-old-grid').eq.0) then
         iret=1
         return
      endif
      call newunit(iun)
      open(unit=iun,file=pwgprefix(1:lprefix)//'xgrid.dat',
     #     form='unformatted',status='old',iostat=ios)
      if(ios.ne.0) then
         iret=-1
         return
      endif
      read(iun,iostat=ios) ((xgrid(j,k),k=1,ndiminteg),j=0,nbins),xint
      read(iun,iostat=ios) ((xgridrm(j,k),k=1,ndiminteg),j=0,nbins),
     1     xintrm
      read(iun,iostat=ios) ((xgridosres(j,k),k=1,ndiminteg),j=0,nbins),
     1     xintosres ! CH, MK: added
      read(iun,iostat=ios) shx, ih1x, ih2x, ndns1x, ndns2x
      if(shx.ne.kn_sbeams.or.ih1x.ne.pdf_ih1.or.ih2x.ne.pdf_ih2
     #  .or.ios.ne.0) then
         iret=-1
         close(iun)
         return
      endif
      close(iun)
      iret=0
      end


      subroutine loadlatestxgridinfo(iret)
c loads x information from gridinfo files, and computes the xgrid.
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flg.h'
      include 'pwhg_par.h'
      include 'cgengrids.h' ! MK: changed position
      integer iret
      character * 40 mergelabels
      character * 6 chnum
      integer j,iteration
      logical probe
      real * 8 ans
c MK: added
#include "Flags.h"
#include "cgengrids_add.h"
      probe=.true.
      do j=par_maxxgriditerations,1,-1
         write(chnum,'(i4)') j
         chnum='xg'//adjustl(chnum)
c remember, tags look like '-1', '-2', etc.
         call loadgridinfo
     1        (mergelabels('btl',chnum,' ',' '),probe,xgrid,ans,iret)
         if(iret.eq.0) then
            iteration = j
            exit
         endif
      enddo
      if(j.eq.0) then
         iret = -1
         return
      endif
      probe = .false.
      write(*,*) ' loading '//chnum//' iteration file'
      call loadgridinfo
     1     (mergelabels('btl',chnum,' ',' '),probe,xgrid,xint,iret)
      if(iret.ne.0) return
      if((flg_withreg.or.flg_withdamp).and..not.flg_bornonly) then
         call loadgridinfo
     1     (mergelabels('rmn',chnum,' ',' '),probe,xgridrm,xintrm,iret)
      endif
c CH, MK: added this for on-shell resonant case
#ifdef DSUB_II
      if(.not.flg_bornonly) then
         call loadgridinfo
     1     (mergelabels('osres',chnum,' ',' '),probe,xgridosres,xintosres,iret)
      endif
#endif
      end
      
      subroutine loadgridinfo(storelabel,probe,xgrid,ans,iret)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_rnd.h'
      include 'pwhg_par.h'
      character *(*) storelabel
      logical probe
      integer nintervals,ndimmax
      parameter (nintervals=50,ndimmax=ndiminteg)
      real * 8 xgrid(0:nintervals,*),ans
      integer iret
      real * 8 xacc(0:nintervals,ndimmax)
      integer nhits(1:nintervals,ndimmax)
      real * 8 xacc0(0:nintervals,ndimmax),tmp
      integer nhits0(1:nintervals,ndimmax),ndim0
      character * 4 chseed
      character * 20 pwgprefix
      character * 100 file
      character * 40 mergelabels
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      integer j,k,iun,ios,kdim
      logical filefound,lpresent
      integer nfiles,jfile,jfound

c first probe if there are files to load
      nfiles=par_maxseeds
      do jfile=1,nfiles
         write(chseed,'(i4)') jfile
         do k=1,4
            if(chseed(k:k).eq.' ') chseed(k:k)='0'
         enddo
         file = mergelabels(pwgprefix(1:lprefix)//'gridinfo',
     1        storelabel,chseed,'.dat')
         inquire(file= file,exist=lpresent)
         if(lpresent) then
            if(probe) then
               iret = 0
               return
            endif
            exit
         endif
      enddo
      if(probe) then
         iret = -1
         return
      endif

      nfiles=par_maxseeds
      filefound=.false.
      jfound=0
      xacc(:,1:ndiminteg)=0
      nhits(:,1:ndiminteg)=0
      ans=0
      do jfile=1,nfiles
         write(chseed,'(i4)') jfile
         do k=1,4
            if(chseed(k:k).eq.' ') chseed(k:k)='0'
         enddo
         file = mergelabels(pwgprefix(1:lprefix)//'gridinfo',
     1        storelabel,chseed,'.dat')
         inquire(file= file,exist=lpresent)
         if(.not.lpresent) cycle
         call newunit(iun)
         open(unit=iun,file= file,
     2        form='unformatted',status='old',iostat=ios)
         if(ios.ne.0) then
            iret=-1
            return
         else
            write(*,*) ' Opened ', file
         endif
         filefound=.true.
         jfound=jfound+1
         read(iun,iostat=ios) ndim0
         if(ios.ne.0.or.ndim0.ne.ndiminteg) goto 111
         read(iun,iostat=ios) tmp
         ans = ans + tmp
         if(ios.ne.0) goto 111
         read(iun,iostat=ios) xgrid(0:nintervals,1:ndiminteg)
         if(ios.ne.0) goto 111
         read(iun,iostat=ios) xacc0(0:nintervals,1:ndiminteg)
         if(ios.ne.0) goto 111
         xacc(:,1:ndiminteg) =
     1        xacc(:,1:ndiminteg) + xacc0(:,1:ndiminteg)
         read(iun,iostat=ios) nhits0(1:nintervals,1:ndiminteg)
         if(ios.ne.0) goto 111
         nhits(:,1:ndiminteg) =
     1        nhits(:,1:ndiminteg) + nhits0(:,1:ndiminteg)
         close(iun)
      enddo
c      write(*,*) ' loadgridinfo: found ',jfound,
c     1     ' files with grid information'
      if(filefound) then
         if(rnd_cwhichseed.eq.'none') then
            call regridplotopen(pwgprefix(1:lprefix)//'-'//
     1          storelabel(1:3) //'grid.top')
         else
            call regridplotopen(pwgprefix(1:lprefix)//'-'//
     1           rnd_cwhichseed//'-'//storelabel(1:3)//'grid.top')
         endif
         ans = ans/jfound
         do kdim=1,ndiminteg
            call regrid(xacc(0,kdim),xgrid(0,kdim),
     1           nhits(1,kdim),kdim,nintervals)
         enddo
         iret = 0
         call regridplotclose
         return
      endif

      iret = -1
      return
 111  continue
      write(*,*) ' loadgridinfo: problems loading files'
      call pwhg_exit(-1)
      end


      subroutine storegridinfo(storelabel,xgrid,xint,xacc,nhits,ndim)
c stores accumulated results and number of hits in gridinfo files
      implicit none
      include 'nlegborn.h'
      include 'pwhg_rnd.h'
      character *(*) storelabel
      integer nintervals,ndimmax
      parameter (nintervals=50,ndimmax=ndiminteg)
      real * 8 xgrid(0:nintervals,*),xint
      real * 8 xacc(0:nintervals,*)
      integer nhits(1:nintervals,*),ndim
      character * 20 pwgprefix
      character * 40 mergelabels
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      integer j,k,iun
      call newunit(iun)
      open(unit=iun,
     1     file=mergelabels(pwgprefix(1:lprefix)//'gridinfo',
     2     storelabel,rnd_cwhichseed,'.dat'),
     3     form='unformatted',status='unknown')
      write(iun) ndim
      write(iun) xint
      write(iun) xgrid(0:nintervals,1:ndim)
      write(iun) xacc(0:nintervals,1:ndim)
      write(iun) nhits(1:nintervals,1:ndim)
      close(iun)
      end


      function mergelabels(lab1,lab2,lab3,lab4)
c puts together up to 4 labels, separating them with '-'.
c Empty labels are ignored.
      character * 40 mergelabels
      character *(*) lab1,lab2,lab3,lab4
      integer where
      mergelabels=' '
      if(lab1.ne.' '.and.lab1.ne.'none') then
         mergelabels=adjustl(lab1)
         where=index(mergelabels,' ')
         if(where.eq.0) goto 999
      endif
      if(lab2.ne.' '.and.lab2.ne.'none') then
         mergelabels(where:)='-'//adjustl(lab2)
         where=index(mergelabels,' ')
         if(where.eq.0) goto 999
      endif
      if(lab3.ne.' '.and.lab3.ne.'none') then
         mergelabels(where:)='-'//adjustl(lab3)
         where=index(mergelabels,' ')
         if(where.eq.0) goto 999
      endif
      if(lab4.ne.' '.and.lab4.ne.'none') then
         mergelabels(where:)='-'//adjustl(lab4)
         where=index(mergelabels,' ')
         if(where.eq.0) goto 999
      endif
c get rid of hiphen before extension
      where=index(mergelabels,'-.',.true.)
      if(where.ne.0) then
         mergelabels(where:)=mergelabels(where+1:)
      endif
      return
 999  write(*,*) ' mergelabels: strings too long'
      call pwhg_exit(-1)
      end

      subroutine writestat(stage)
      implicit none
      character *(*) stage
      character * 40 mergelabels
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'pwhg_rnd.h'
      include 'pwhg_flg.h'
c MK: added
#include "osres.h"
#include "pwhg_rad_add.h"
#include "Flags.h"
      real*8 temparrosres(2,ntot_osres), totneg ! MK: added
      real*8 rad_totnegosres_sum ! MK: added 
      integer ichan ! MK: added to distinguish the on-shell resonances
      integer j,k
      integer iunstat
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      call newunit(iunstat)
      open(unit=iunstat,file=mergelabels(pwgprefix,stage,rnd_cwhichseed,
     1     'stat.dat'),status='unknown')

      write(iunstat,*) 'btilde pos.   weights:', rad_totposbtl,' +-',
     1     rad_etotposbtl
      write(iunstat,*) 'btilde |neg.| weights:', rad_totnegbtl,' +-',
     1     rad_etotnegbtl
      write(iunstat,*) 'btilde Total (pos.-|neg.|):', rad_totbtl,
     1     ' +-',rad_etotbtl
      ! MK: addded or modified
      ! ================================================================
#ifdef DSUB_II
      ! MK: combine the results of all on shell resonances:
      ! this is called only once and is mandatory
      call update_totarr ! MK: OK!
      do k=1,2
        do j=1,ntot_osres
          temparrosres(k,j)=0D0
          do ichan=1,nosres
            temparrosres(k,j)=temparrosres(k,j)+
     &                        rad_osres_totarr(k,j,ichan)
          enddo
        enddo
      enddo     
      write(iunstat,*) 'real_osres pos.   weights:', 
     1        temparrosres(1,3),' +-',temparrosres(2,3)
      write(iunstat,*) 'real_osres |neg.| weights:',
     1         temparrosres(1,4),' +-',temparrosres(2,4)
      write(iunstat,*) 'real_osres total (pos.-|neg.|):', 
     1        temparrosres(1,1),' +-',temparrosres(2,1)
#endif
      if(flg_withreg.or.flg_withdamp) then
         write(iunstat,*) ' Remnant cross section in pb',
!      1        rad_totrm,'+-',rad_etotrm
     1        rad_totrem,'+-',rad_etotrem
      endif
      if (flg_weightedev) then
         write(iunstat,*)
     1 ' total (btilde+remnants+regulars+osresR) cross section times '
         write(iunstat,*) ' suppression factor in pb',
     1        rad_tot,'+-',rad_etot         
      else
         write(iunstat,*)
     1 ' total (btilde+remnants+regulars+osresR) cross section in pb',
     2     rad_tot,'+-',rad_etot
      endif
      rad_totnegosres_sum = 0D0
      do ichan=1,nosres
        rad_totnegosres_sum   = rad_totnegosres_sum + 
     &                          rad_totnegosres(ichan)
      enddo
      totneg=rad_totnegreg+rad_totnegrem+rad_totnegosres_sum+
     &          rad_totnegbtl ! CH, MK: added
      write(iunstat,*) ' negative weight fraction:',
     1              totneg/(2*totneg+rad_tot)
c     1     rad_totnegbtl/(2*rad_totnegbtl+rad_tot)
      ! ================================================================
      close(iunstat)
      end

c MK: added this routine
c save all variables in array rad_totarr and rad_osres_totarr
      subroutine update_totarr
        implicit none

#include "nlegborn.h"
#include "pwhg_flst.h"
#include "pwhg_rad.h"
#include "osres.h"
#include "pwhg_rad_add.h"

        integer ichan

        ! btilde contributions 1:4
        rad_totarr(1,1) = rad_totbtl
        rad_totarr(2,1) = rad_etotbtl
        rad_totarr(1,2) = rad_totabsbtl
        rad_totarr(2,2) = rad_etotabsbtl
        rad_totarr(1,3) = rad_totposbtl
        rad_totarr(2,3) = rad_etotposbtl
        rad_totarr(1,4) = rad_totnegbtl
        rad_totarr(2,4) = rad_etotnegbtl

        ! regular contributions 5:8
        rad_totarr(1,5) = rad_totreg
        rad_totarr(2,5) = rad_etotreg
        rad_totarr(1,6) = rad_totabsreg
        rad_totarr(2,6) = rad_etotabsreg
        rad_totarr(1,7) = rad_totposreg
        rad_totarr(2,7) = rad_etotposreg
        rad_totarr(1,8) = rad_totnegreg
        rad_totarr(2,8) = rad_etotnegreg

        ! remnant contributions 9:15
        rad_totarr(1,9)  = rad_totrem
        rad_totarr(2,9)  = rad_etotrem
        rad_totarr(1,10) = rad_totabsrem
        rad_totarr(2,10) = rad_etotabsrem
        rad_totarr(1,11) = rad_totposrem
        rad_totarr(2,11) = rad_etotposrem
        rad_totarr(1,12) = rad_totnegrem
        rad_totarr(2,12) = rad_etotnegrem

        rad_totarr(1,13) = rad_totbtlgen
        rad_totarr(2,13) = rad_etotbtlgen
        rad_totarr(1,14) = rad_totreggen
        rad_totarr(2,14) = rad_etotreggen
        rad_totarr(1,15) = rad_totremgen
        rad_totarr(2,15) = rad_etotremgen

        ! totals 16:17
        rad_totarr(1,16) = rad_tot
        rad_totarr(2,16) = rad_etot
        rad_totarr(1,17) = rad_totgen
        rad_totarr(2,17) = rad_etotgen

        ! on-shell contributions
        do ichan=1,nosres
          rad_osres_totarr(1,1,ichan) = rad_totosres(ichan)
          rad_osres_totarr(2,1,ichan) = rad_etotosres(ichan)
          rad_osres_totarr(1,2,ichan) = rad_totabsosres(ichan)
          rad_osres_totarr(2,2,ichan) = rad_etotabsosres(ichan)
          rad_osres_totarr(1,3,ichan) = rad_totpososres(ichan)
          rad_osres_totarr(2,3,ichan) = rad_etotpososres(ichan)
          rad_osres_totarr(1,4,ichan) = rad_totnegosres(ichan)
          rad_osres_totarr(2,4,ichan) = rad_etotnegosres(ichan)
          rad_osres_totarr(1,5,ichan) = rad_totosresgen(ichan)
          rad_osres_totarr(2,5,ichan) = rad_etotosresgen(ichan)
        enddo
      end

c MK: added this routine
c write back the array rad_totarr and rad_osres_totarr to the original values
      subroutine apply_totarr
        implicit none

#include "nlegborn.h"
#include "pwhg_flst.h"
#include "pwhg_rad.h"
#include "osres.h"
#include "pwhg_rad_add.h"

        integer ichan

        ! btilde contributions 1:4
        rad_totbtl     = rad_totarr(1,1) 
        rad_etotbtl    = rad_totarr(2,1) 
        rad_totabsbtl  = rad_totarr(1,2) 
        rad_etotabsbtl = rad_totarr(2,2) 
        rad_totposbtl  = rad_totarr(1,3) 
        rad_etotposbtl = rad_totarr(2,3) 
        rad_totnegbtl  = rad_totarr(1,4) 
        rad_etotnegbtl = rad_totarr(2,4) 

        ! regular contributions 5:8
        rad_totreg     = rad_totarr(1,5)
        rad_etotreg    = rad_totarr(2,5)
        rad_totabsreg  = rad_totarr(1,6)
        rad_etotabsreg = rad_totarr(2,6)
        rad_totposreg  = rad_totarr(1,7)
        rad_etotposreg = rad_totarr(2,7)
        rad_totnegreg  = rad_totarr(1,8)
        rad_etotnegreg = rad_totarr(2,8)

        ! remnant contributions 9:15
        rad_totrem     = rad_totarr(1,9) 
        rad_etotrem    = rad_totarr(2,9) 
        rad_totabsrem  = rad_totarr(1,10)
        rad_etotabsrem = rad_totarr(2,10)
        rad_totposrem  = rad_totarr(1,11)
        rad_etotposrem = rad_totarr(2,11)
        rad_totnegrem  = rad_totarr(1,12)
        rad_etotnegrem = rad_totarr(2,12)

        rad_totbtlgen  = rad_totarr(1,13)
        rad_etotbtlgen = rad_totarr(2,13)
        rad_totreggen  = rad_totarr(1,14)
        rad_etotreggen = rad_totarr(2,14)
        rad_totremgen  = rad_totarr(1,15)
        rad_etotremgen = rad_totarr(2,15)

        ! totals 16:17
        rad_tot     = rad_totarr(1,16)
        rad_etot    = rad_totarr(2,16)
        rad_totgen  = rad_totarr(1,17)
        rad_etotgen = rad_totarr(2,17)

        ! on-shell contributions
        do ichan=1,nosres
          rad_totosres(ichan)     = rad_osres_totarr(1,1,ichan)
          rad_etotosres(ichan)    = rad_osres_totarr(2,1,ichan)
          rad_totabsosres(ichan)  = rad_osres_totarr(1,2,ichan)
          rad_etotabsosres(ichan) = rad_osres_totarr(2,2,ichan)
          rad_totpososres(ichan)  = rad_osres_totarr(1,3,ichan)
          rad_etotpososres(ichan) = rad_osres_totarr(2,3,ichan)
          rad_totnegosres(ichan)  = rad_osres_totarr(1,4,ichan)
          rad_etotnegosres(ichan) = rad_osres_totarr(2,4,ichan)
          rad_totosresgen(ichan)  = rad_osres_totarr(1,5,ichan)
          rad_etotosresgen(ichan) = rad_osres_totarr(2,5,ichan)
        enddo
      end
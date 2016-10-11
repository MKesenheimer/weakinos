c MK: copied and modified version of gen_radiation.f, revision 3154
c modified on the basis of disquark/gen_radiation.f
c changes marked with "! CH, MK:" or "! MK:"
c bigger changes over a whole section are marked with !===...

      subroutine pwhgevent
      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      include 'pwhg_rwl.h'
      include 'LesHouches.h'
c MK: added
#include "osres.h"
#include "pwhg_flst_add.h"
#include "Flags.h"
#include "pwhg_rad_add.h"
      integer iret,iun
      real * 8 suppfact
      real * 8 random,powheginput
      external random,powheginput
      integer mcalls,icalls
      data mcalls,icalls/0,0/
      save mcalls,icalls
      real * 8 pwhg_pt2,pt2max_regular,pt2max_osres ! CH: added pt2max_osres
      external pwhg_pt2,pt2max_regular,pt2max_osres
      integer ichan ! MK: added
      double precision rad_totosresgen_sum ! MK: added
      real * 8 weight
      real*8 ran1 ! CH: added
      integer i
      real * 8 seconds
      integer lh_seed,lh_n1,lh_n2
      common/lhseeds/lh_seed,lh_n1,lh_n2
      logical notfinite_kin
      flg_monitorubound = .true.
c If at the end the event is not generated for some reason (nup=0)
c restart from here
 1    continue
      if(idwtup.eq.3) then
         weight=1
      elseif(idwtup.eq.-4) then
         weight=rad_totgen * rad_branching
      else
         write(*,*) ' only 3 and -4 are allowed for idwtup'
         call exit(-1)
      endif
c     store current random seeds. To be used to restart at problematic events
      call readcurrentrandom(lh_seed,lh_n1,lh_n2)
c      if(random().gt.rad_totrm/rad_totgen) then
      ! CH: changed this: regular part now osres in several parts
      ran1=random() ! CH: added
      ! MK: added the following lines
      !=================================================================
      rad_totosresgen_sum = 0D0
      do ichan=1,nosres
        rad_totosresgen_sum   = rad_totosresgen_sum + 
     &                          rad_totosresgen(ichan)
      enddo
      !=================================================================
      if(ran1.gt.(rad_totremgen+rad_totreggen+rad_totosresgen_sum)/ ! MK: changed
     &                rad_totgen) then
           ! CH: set the flg_btilde to true for btilde part
           flg_btilde=.true.

c     generate underlying Born kinematics
         call reset_timer
         call gen_btilde(mcalls,icalls)
         if(notfinite_kin('Born')) goto 1
         call get_timer(seconds)
         call addtocnt('btilde time (sec)',seconds)
c     generate underlying Born flavour
         call gen_uborn_idx
c
         if(powheginput('#testsuda').eq.1) then
            call testsuda
         endif
         if(.not.flg_LOevents) then
c generate radiation
            call reset_timer
            call gen_radiation
            if(notfinite_kin('Real')) goto 1
            call get_timer(seconds)
            call addtocnt('radiation time (sec)',seconds)
            rad_pt2max=pwhg_pt2()
         else
            kn_csi = 0
         endif
c add a random azimuthal rotation around beam axis
         call add_azimuth
c --- set up les houches interface
         call gen_leshouches
c if negative weight, flip the sign of weight
         if(rad_btilde_sign(rad_ubornidx).eq.-1) then
            weight=-weight
         endif
c rad_type=1 for btilde events (used only for debugging purposes)
         rad_type=1
         call increasecnt("btilde event")
         rwl_type = rad_type
         rwl_index = rad_ubornidx
         rwl_weight = rad_btilde_arr(rad_ubornidx)
     $        *rad_btilde_sign(rad_ubornidx)
      ! regular or remnant contribution
      elseif(ran1.gt.(rad_totosresgen_sum/rad_totgen)) then ! MK: changed
         ! CH: set the flg_btilde to true for normal remnant part
         flg_btilde=.true.

c generate remnant n+1 body cross section
         call reset_timer
         call gen_sigremnant
         if(notfinite_kin('Real')) goto 1
         call get_timer(seconds)
         call addtocnt("remnant time (sec)",seconds)
c pick a configuration according to its cross section
c iret=1: rem contribution (leftover from damping factor on R)
c iret=2: reg contribution (real graphs without singular regions)
c and regenerate real phase space accordingly
         call gen_remnant(iret)
         include 'post_gen_remnant_hook.h'
         if(notfinite_kin('Real')) goto 1
c         if (pwhg_pt2().lt.rad_ptsqmin) then
c            write(*,*) '****************************************'
c            write(*,*) 'WARNING in gen_remnant'
c            write(*,*) 'pwhg_pt2 < rad_ptsqmin ',
c     #           pwhg_pt2(),' < ',rad_ptsqmin
c            write(*,*) (flst_alr(i,rad_realalr),i=1,nlegreal)
c            write(*,*) 'To generate this event, use the following seeds'
c            call printcurrentrandom
c            write(*,*) '****************************************'
c         endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         call add_azimuth
         if(iret.eq.1) then
c     set st_muren2 equal to pt2 for scalup value
            rad_pt2max=max(rad_ptsqmin,pwhg_pt2())
            call set_rad_scales(rad_pt2max)
            call gen_leshouches
c     rad_type=2 for remnants
            rad_type=2
            if(rad_damp_rem_sign(rad_realalr).eq.-1) then ! CH: remnant might be negative, if we split the reals
              weight=-weight
            endif
            call increasecnt("remnant event")
            rwl_type = rad_type
            rwl_index = rad_realalr
            rwl_weight = rad_damp_rem_arr(rad_realalr)
         else
c     set st_muren2 for scalup value for regular contributions
            rad_pt2max=max(rad_ptsqmin,pt2max_regular())
            call set_rad_scales(rad_pt2max)
            call gen_leshouches_reg
c rad_type=3 for regular contributions
            rad_type=3
            if(rad_reg_sign(rad_realreg).eq.-1) then ! CH: should never occur, but anyway...
              weight=-weight
            endif
            call increasecnt("regular event")
            rwl_type = rad_type
            rwl_index = rad_realreg
            rwl_weight = rad_reg_arr(rad_realreg)
         endif
c=======================================================================        
c CH, MK: new part here:
      else
         ! CH: set the flg_btilde to false for osres part
         flg_btilde=.false.
         ! generate osres n+1 body cross section
         call gen_sigosres
         ! MK:
         ! pick a configuration according to its cross section
         ! iret=3: first on-shell resonant real contribution
         ! iret=4: second on-shell resonant real contribution
         ! remember:     iret = ichan + 2 (see gen_index_mod.f) and
         !           rad_type = ichan + 3 (see pwhgreweight_mod.f)
         ! therefor  rad_type = iret + 1
         ! this is in full analogy with disquark
         call gen_osres(iret)
         ! MK: changed the following lines
         ! if negative weight, flip the sign of xwgtup
         if( (iret.ge.3) .and. (iret.le.(nosres+2)) ) then
           if(rad_osres_sign(rad_realosres,iret-2).eq.-1) then
             weight = -weight
           endif
         else
            print*, 'Wrong return-statement of gen_osres: iret=',iret
            stop        
         endif

         call add_azimuth

         if( (iret.ge.3) .and. (iret.le.(nosres+2)) ) then
            ! set st_muren2 for scalup value for osres contributions 
            ! (alphas should be set correctly)
            ! pt2max_osres-fct in Real_osres_phsp.f
            rad_pt2max=max(rad_ptsqmin,pt2max_osres())
            call set_rad_scales(rad_pt2max)
            call gen_leshouches_osres
            rad_type = iret + 1
         else
            print*, "error in gen_radiation_mod.f: invalid iret-type."
            print*, "iret: ", iret
            print*, "nosres: ", nosres
            stop
         endif   
         call increasecnt("osres event")
         rwl_type = rad_type
         rwl_index = rad_realosres
         rwl_weight = rad_osres_arr(rad_realosres,iret-2)
     &               *rad_osres_sign(rad_realosres,iret-2)
c======================================================================= 
      endif
      if(flg_weightedev) then
         call born_suppression(suppfact)
         if(suppfact.eq.0) then
            write(*,*) ' 0 suppression factor in event generation'
            write(*,*) ' aborting'
            call exit(-1)
         endif
         weight=weight/suppfact
      endif
c     correct for bound violations
      if(flg_ubexcess_correct) then
         weight = weight * rad_genubexceeded
      endif
c If at the end the event is not generated for some reason (nup=0)
c restart from here
      if(nup.eq.0) goto 1
      xwgtup = weight
      end

      logical function notfinite_kin(BornOrReal)
      implicit none
      character * 4 BornOrReal
      logical pwhg_isfinite
      integer j,mu
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      if(BornOrReal.eq.'Born') then
         do j=1,nlegborn
            do mu=0,3
               if(.not.pwhg_isfinite(kn_pborn(mu,j))) then
                  call increasecnt("not_finite_kin in Born")
                  notfinite_kin = .true.
                  return
               endif
            enddo
         enddo
      elseif(BornOrReal.eq.'Real') then
         do j=1,nlegreal
            do mu=0,3
               if(.not.pwhg_isfinite(kn_preal(mu,j))) then
                  call increasecnt("not_finite_kin in Real")
                  notfinite_kin = .true.
                  return
               endif
            enddo
         enddo
      else
         write(*,*) 'not_finite_kin: error, argument should be'
     1        //'either Born or Real, got ',BornOrReal
         call exit(-1)
      endif
      notfinite_kin = .false.
      end




      subroutine gen_radiation
      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flg.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      real * 8 t,csi,y,azi,sig,born
      real * 8 tmax
      common/ctmax/tmax
      integer kinreg,firstreg,lastreg,fl1,fl2,flemitter
      logical ini
      data ini/.true./
      real * 8 pwhg_pt2,powheginput
      external pwhg_pt2,powheginput
      save ini,firstreg,lastreg
      logical is_charged,is_coloured
      external is_charged,is_coloured
      if(ini) then
         firstreg=powheginput("#radregion")
         if(firstreg.le.0) then
            firstreg=1
            lastreg=rad_nkinreg
         else
            lastreg=firstreg
         endif
         ini=.false.
      endif
c Use highest bid procedure (see appendix B of FNO2006)
      tmax=0
      kinreg=0
      do rad_kinreg=firstreg,lastreg
         if(rad_kinreg_on(rad_kinreg)) then
            if(rad_kinreg.eq.1) then
c     initial state radiation
c kn_emitter may be 0,1,2 depending upon the flavour
c of the process, which is undefined here.
c Set it to a value less than 2, to avoid problems later.
               kn_emitter = 0
               fl1=flst_born(1,rad_ubornidx)
               fl2=flst_born(2,rad_ubornidx)
               if((.not.is_coloured(fl1).and..not.is_coloured(fl2))
     1          .and.(is_charged(fl1).or.is_charged(fl2))) then
                  flg_em_rad = .true.
               else
                  flg_em_rad = .false.
               endif
               call gen_rad_isr(t)
            else
c     final state radiation
               kn_emitter=flst_lightpart+rad_kinreg-2
               flemitter=flst_born(kn_emitter,rad_ubornidx)
               if(.not.is_coloured(flemitter).and.is_charged(flemitter))
     1              then
                  flg_em_rad = .true.
               else
                  flg_em_rad = .false.
               endif
               call gen_rad_fsr(t)
            endif
            include 'pwhg_gen_radiation_hook.h'
            if(t.gt.tmax) then
               tmax=t
               kinreg=rad_kinreg
               csi=kn_csi
               y=kn_y
               azi=kn_azi
            endif
         endif
      enddo
c Set up radiation kinematics
      if(kinreg.eq.0) then
c Generate a Born like event
         kn_csi=0
         rad_kinreg=0
         return
      else         
         rad_kinreg=kinreg
         kn_csi=csi
         kn_y=y
         kn_azi=azi
         if(rad_kinreg.eq.1) then
            call gen_real_phsp_isr_rad
         else
            call gen_real_phsp_fsr_rad
         endif
         t=pwhg_pt2()
         call set_rad_scales(t)
c We call sigborn_rad now, becayse the real may depend
c upon the Born through the soft and collinear terms,
c that are used in the real if bornzerodamp is used.
c Failing to do so may cause problems in picking the
c flavour
         call sigborn_rad(born)
         call sigreal_rad(sig)
         call gen_real_idx
      endif
      end



      function pwhg_pt2()
      implicit none
      real * 8 pwhg_pt2
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      real * 8 pres(0:3),q2
      integer em,ires
      if(rad_kinreg.eq.1) then
         pwhg_pt2=(kn_sreal/4)*(1-kn_y**2)*kn_csi**2
      else
         em=flst_lightpart+rad_kinreg-2
         if(kn_masses(em).eq.0) then
            if(flst_bornres(em,1).ne.0) then
               ires=flst_bornres(em,1)
               pres=kn_cmpborn(:,ires)
               q2=pres(0)**2-pres(1)**2-pres(2)**2-pres(3)**2
            else
               q2=kn_sreal
            endif
            pwhg_pt2=(q2/2)*(1-kn_y)*kn_csi**2
         else
            call comppt2fsrmv(kn_y,kn_csi,pwhg_pt2)
         endif
      endif
      end

      function pwhg_upperb_rad()
      implicit none
      real * 8 pwhg_upperb_rad
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      include 'pwhg_st.h'
      include 'pwhg_em.h'
      real * 8 x,y,csi
      integer em
      csi=kn_csi
      x=1-csi
      y=kn_y
      if(rad_kinreg.eq.1) then
         if(rad_iupperisr.eq.1) then
            pwhg_upperb_rad = 1/((1-x)*(1-y**2))
c Possible alternatives:
c rad_iupper=2   pwhg_upperb_rad = 1/(x*(1-x)*(1-y**2))
c 
c rad_iupper=3:  pwhg_upperb_rad = 1/(x**2*(1-x)*(1-y**2))
         else
            write(*,*) ' rad_iupper=',rad_iupperisr,
     1        'alternative not implemented'
            call exit(1)
         endif
      else
c Final state radiation
         em=flst_lightpart+rad_kinreg-2
         if(kn_masses(em).eq.0) then
c for now use the same
            if(rad_iupperfsr.eq.1) then
               pwhg_upperb_rad = 1/(csi*(1-y))
            elseif(rad_iupperfsr.eq.2) then
               pwhg_upperb_rad = 1/(csi**2*(1-y)*(1-csi/2*(1-y))**2)
     2              *csi
            elseif(rad_iupperfsr.eq.3) then
               pwhg_upperb_rad = 1/(csi*(1-y)*
     2              (1-csi/2*(1-y)))
            else
               write(*,*) ' rad_iupper=',rad_iupperfsr,
     1              'alternative not implemented'
               call exit(1)
            endif
         else
c     massive emitter
            call compubradmv(y,csi,pwhg_upperb_rad)
         endif
      endif
      if(flg_em_rad) then
         pwhg_upperb_rad = pwhg_upperb_rad * em_alpha
      else
         pwhg_upperb_rad = pwhg_upperb_rad * st_alpha
      endif

      end



      function pt2solve(pt2,i)
c Returns  xlr - log(Delta^{(tilde{V})}) , see eq. D14, D15 in ZZ paper
c We use it to find its zero in pt2.
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_kn.h'
      include 'pwhg_em.h'
      include 'pwhg_rad.h'
      include 'pwhg_math.h'
      real * 8 pt2solve,pt2
c i set by dzero: 1 for first call, 2 for subsequent calls, 3 for last call
c before a normal exit; not used here
      integer i,em
      real * 8 xlr,q2,xlam2c,kt2max,unorm,cunorm,sborn
      integer nlc
      common/cpt2solve/xlr,q2,kt2max,xlam2c,unorm,sborn,nlc
      real * 8 b0,xm,p,tmp
      b0=(11*CA-4*TF*nlc)/(12*pi)
      if(flg_em_rad) then
         cunorm=unorm*em_alpha
      else
         cunorm=unorm
      endif
      if(rad_kinreg.eq.1) then
         if(rad_iupperisr.eq.1) then
c see Notes/upperbounding-isr.pdf
            if(pt2.lt.sborn) then
               if(sborn.lt.kt2max) then
                  pt2solve=cunorm*pi/b0*(
     #      (log(2*sborn/xlam2c)*log(log(sborn/xlam2c)/log(pt2/xlam2c))
     #        - log(sborn/pt2)) +
     #           log(2d0)*log(log(kt2max/xlam2c)/log(sborn/xlam2c)))
     #           + xlr
               else
                  pt2solve=cunorm*pi/b0*(
     #     (log(2*sborn/xlam2c)*log(log(kt2max/xlam2c)/log(pt2/xlam2c))
     #        - log(kt2max/pt2)) )
     #           + xlr
               endif
            else
               pt2solve=cunorm*pi/b0*(log(2d0)
     #           *log(log(kt2max/xlam2c)/log(pt2/xlam2c)))
     #           + xlr
            endif
         else
            write(*,*) ' rad_iupper=',rad_iupperisr,' not implemented'
            call exit(1)
c Alternatives: rad_iupper=2
c         pt2solve=cunorm*pi/b0/2
c     #        *(log(q2/xlam2c)*log(log(kt2max/xlam2c)/log(pt2/xlam2c))
c     #        - log(kt2max/pt2)) + xlr
         endif
      else
         em = flst_lightpart+rad_kinreg-2
         if(kn_masses(em).ne.0) then
            call compintub(pt2,pt2solve)
c The following lines are used to test the analytic integration
c versus a vegas one; uncomment to test
c            call compintubveg(pt2,tmp)
c            write(*,'(a,3(1x,d10.4))') ' testintub:',pt2,pt2solve,tmp
            pt2solve=cunorm*pt2solve+xlr
         else
            if(rad_iupperfsr.eq.1) then
c final state radiation
               pt2solve=cunorm*pi/b0*(
     1       (log(kt2max/xlam2c)*log(log(kt2max/xlam2c)/log(pt2/xlam2c))
     2              - log(kt2max/pt2)) )
     3              + xlr
            elseif(rad_iupperfsr.eq.2) then
               xm=kn_csimax
               p=sqrt(pt2/sborn)
               pt2solve=cunorm*2*pi*2*(
     3        (log(xm-xm**2)+(2*xm-2)*log(xm)-2*log(1-xm)*xm-2)/xm/2.d+0
     1     -(p*log(xm-p**2)+(2*p*log(p)-2*log(1-p)*p-2)*xm-2*p*log(p))
     2    /(p*xm)/2.d+0) + xlr
            elseif(rad_iupperfsr.eq.3) then
               xm=kn_csimax
               p=sqrt(pt2/sborn)
               pt2solve=cunorm*2*pi*2*(
     3     (log(xm-xm**2)+(2*xm-2)*log(xm)-2*log(1-xm)*xm-2)/xm/2.d+0
     1   -(p*log(xm-p**2)+(2*p*log(p)-2*log(1-p)*p-2)*xm-2*p*log(p))
     2   /(p*xm)/2.d+0) + xlr
            else
               write(*,*)
     1 ' rad_iupper=',rad_iupperfsr,' not implemented'
               call exit(1)
            endif            
         endif
      endif
      end

      subroutine gen_rad_isr(t)
c Generates hard radiation kinematics according to
c appendix D in ZZ paper.
c
c  common/cptmin/ptminsq: minimum pt^2 accepted
c 
      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      include 'pwhg_st.h'
      real * 8 t
      real * 8 x,y,x1b,x2b
      real * 8 xlr,q2,xlam2c,kt2max,unorm,sborn
      integer nlc
      common/cpt2solve/xlr,q2,kt2max,xlam2c,unorm,sborn,nlc
      real * 8 xmin,rv,xp,xm,chi,tk,uk,ubound,ufct,
     #   value,err,tmp1,tmp2,tmp,rvalue,born,sig
      common/cdfxmin/xmin
      real * 8 tmax
      common/ctmax/tmax
      real * 8 random,pt2solve,dfxmin,pwhg_alphas0,pwhg_upperb_rad
      external random,pt2solve,dfxmin,pwhg_alphas0,pwhg_upperb_rad
      unorm=rad_norms(rad_kinreg,rad_ubornidx)
      sborn=kn_sborn
      x1b=kn_xb1
      x2b=kn_xb2
c See Notes/kt2max.pdf
      kt2max = sborn*(1-x2b**2)*(1-x1b**2)/(x1b+x2b)**2
      if(kt2max.lt.rad_ptsqmin.or.kt2max.lt.tmax) then
         t=-1
         goto 3
      endif
c upper bound is log(q2/t)
      if(rad_iupperisr.eq.1) then
         q2=2*sborn
      elseif(rad_iupperisr.eq.2) then
         write(*,*) ' rad_iupper=',rad_iupperisr,' not implemented'
         call exit(1)
c Alternative rad_iupper=2
c         q2=4*sborn/min(x1b,x2b)**2
      endif
c see section 4 in ZZ paper, last paragraph
      xlam2c=rad_lamll**2
      nlc=5
      xlr=0
 1    continue
      xlr=xlr+log(random())
c CERNLIB voodoo:
      call KERSET('C205. ',0,0,101)
c solve for zero of pt2solve
c dzero(xmin,xmax,x,err,eps,maxcalls,function)
c err: on exit if no error occours: |y-y0|<err 
c      error C205.1 function(xmin)*function(xmax)>0,
c                   x=0 and r=-2(ymax-ymin)
c      error C205.2 Number of calls to F exceeds maxcalls,
c                   x=0 and r=-(xmax-xmin)/2
c eps: required accuracy
      call dzero(rad_ptsqmin,kt2max,t,err,1d-8,1000000,pt2solve)
c error conditions
      if(t.eq.0.and.err.lt.0d0
     # .and.err.gt.rad_ptsqmin-kt2max) then
         write(*,*) 'DZERO fails'
         write(*,*) ' number of calls exceeded'
         call exit(1)
      endif
 3    if(t.lt.rad_ptsqmin.or.t.lt.tmax) then
c below cut (either below absolute minimum, or below previously generated
c radiation in highest bid loop): generate a born event
         t=-1
         kn_csi=0
         return
      endif
c vetoes:
      rv=random()
      xp=(sqrt(1+t/sborn)+sqrt(t/sborn))**2
      xm=(sqrt(1+t/sborn)-sqrt(t/sborn))**2
c tmp1: V(t)/tilde{V}(t) in appendix D of ZZ paper;
c (typo: in D.13, log log -> log
      xmin=min(x1b,x2b)/(2*sqrt(1+t/sborn))
      if(rad_iupperisr.eq.1) then
         tmp1=log((sqrt(xp-xmin)+sqrt(xm-xmin))
     #       /(sqrt(xp-xmin)-sqrt(xm-xmin)))
         if(t.lt.sborn) then
            tmp1=tmp1/(log(2*sborn/t)/2)
         else
            tmp1=tmp1/(log(2d0)/2)
         endif
      elseif(rad_iupperisr.eq.2) then
         tmp1=log(2/xmin*(sqrt((xp-xmin)*(xm-xmin))
     # +1-xmin/2*(xp+xm))/(xp-xm)) /(log(q2/t)/2)
      endif
c compare with D.11-D.12
c to set xmuren2:
      call set_rad_scales(t)
      tmp2=st_alpha / pwhg_alphas0(t,rad_lamll,nlc)
      tmp=tmp1*tmp2
      if(tmp.gt.1) then
         write(*,*) ' Error: upper bound lower than actual value',
     #        tmp,tmp1,tmp2,t,st_alpha
         call exit(1)
      endif
      if(rv.gt.tmp) then
         goto 1
      endif
c At this stage: pt generated according to D.2
c generate x proportional to 1/(x sqrt((xp-x)*(xm-x)))
c in the range xmin < x < xm (cf. D.5)
c Generate in d sqrt(xm-x) /sqrt(xp-x)  (rad_iupper=1) or d sqrt(xm-x) /(x sqrt(xp-x)) (rad_iupper=2)
c using       d sqrt(xm-x) /sqrt(xp-xm) (rad_iupper=1) or d sqrt(xm-x) /(xmin sqrt(xp-xm)) (rad_iupper=2) as upper bound using hit and miss
 2    chi=sqrt(xm-xmin)*random()
      x=xm-chi**2
      if(rad_iupperisr.eq.1) then
         if(random().gt.sqrt(xp-xm)/sqrt(xp-x)) goto 2
      elseif(rad_iupperisr.eq.2) then
         if(random().gt.(xmin*sqrt(xp-xm))/(x*sqrt(xp-x))) goto 2
      endif
c get y (abs to avoid tiny negative values)
      y=sqrt(abs(1-4*x/(1-x)**2*t/sborn))
      if(random().gt.0.5d0) y=-y
c At this point an x-y pair is generated according to the
c distribution upper().
c
c Veto if out of range (x1>1 or x2>1)
      tk=-1d0/2*(1-x)*(1-y)
      uk=-1d0/2*(1-x)*(1+y)
      if(   x1b*sqrt((1+tk)/(1+uk)/x) .gt. 1
     # .or. x2b*sqrt((1+uk)/(1+tk)/x) .gt. 1) then
         goto 1
      endif
c extra suppression factor of upper bounding function (may depend upon radiation variables)
      call uboundfct(ufct,1-x,y)
      if(random().gt.ufct) goto 1
c Veto from upper bound to real value. Count how many vetoes,
c since these may be expensive.
      call sigborn_rad(born)
      if(born.lt.0) then
         born=0
      endif
      if(born.eq.0) then
c bizarre situation that may arise when the scale gets so low
c that some pdf vanish (typically heavy flavour pdf's)
         t=-1
         goto 3
      endif
      kn_y=y
      kn_csi=1-x
      kn_azi=2*pi*random()
      ubound=born*pwhg_upperb_rad()*unorm*ufct
      call gen_real_phsp_isr_rad
      call sigreal_rad(sig)
      value=sig*kn_jacreal
      if(value.gt.ubound) then
         call increasecnt(
     # 'upper bound failures in generation of radiation')
      endif
      rvalue=random()*ubound
      if(rvalue.gt.value) then
         call increasecnt('vetoed radiation')
         goto 1
      endif
      end


      subroutine getkt2maxands(kt2,s)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_mvem.h'
      real * 8 kt2,s
c setupmvemitter fixed to works also in the massless case
      call setupmvemitter
      kt2=kt2max
      s=q**2
      end

      subroutine gen_rad_fsr(t)
c Generates final state hard radiation kinematics according to
c Notes/upperbounding-fsr.pdf
      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      include 'pwhg_st.h'
      real * 8 t
      real * 8 csi,y
      real * 8 xlr,q2,xlam2c,kt2max,unorm,sborn
      integer nlc
      common/cpt2solve/xlr,q2,kt2max,xlam2c,unorm,sborn,nlc
      real * 8 xmin,rv,ubound,ufct,
     #   s,value,err,tmp,rvalue,born,sig
      common/cdfxmin/xmin
      real * 8 tmax
      common/ctmax/tmax
      real *8 kt2minqed
      common/showerqed/kt2minqed
      real * 8 random,pt2solve,pwhg_alphas0,pwhg_upperb_rad,pwhg_pt2
      external random,pt2solve,pwhg_alphas0,pwhg_upperb_rad,pwhg_pt2
      unorm=rad_norms(rad_kinreg,rad_ubornidx)
c kn_sborn=kn_sreal:
      call getkt2maxands(kt2max,s)
      sborn=s
c below is for the QED case; it will never hit that limit anyhow ...
      if(flg_em_rad) then
         if(rad_iupperfsr.eq.1) then
            write(*,*)
     1  'error gen_rad_fsr: '
            write(*,*) 
     1  'cannot use iupperfsr = 1 for electromagnetic radiation'
            call exit(-1)
         endif
         if(kt2max.lt.kt2minqed.or.kt2max.lt.tmax) then
            t=-1
            goto 3
         endif
      else
         if(kt2max.lt.rad_ptsqmin.or.kt2max.lt.tmax) then
            t=-1
            goto 3
         endif
      endif
c see section 4 in ZZ paper, last paragraph
      xlam2c=rad_lamll**2
      nlc=5
      xlr=0
 1    continue
      xlr=xlr+log(random())
c CERNLIB voodoo:
      call KERSET('C205. ',0,0,101)
c solve for zero of pt2solve
c dzero(xmin,xmax,x,err,eps,maxcalls,function)
c err: on exit if no error occours: |y-y0|<err 
c      error C205.1 function(xmin)*function(xmax)>0,
c                   x=0 and r=-2(ymax-ymin)
c      error C205.2 Number of calls to F exceeds maxcalls,
c                   x=0 and r=-(xmax-xmin)/2
c eps: required accuracy
      call dzero(rad_ptsqmin,kt2max,t,err,1d-8,1000000,pt2solve)
c error conditions
      if(t.eq.0.and.err.lt.0d0
     # .and.err.gt.rad_ptsqmin-kt2max) then
         write(*,*) 'DZERO fails'
         write(*,*) ' number of calls exceeded'
         call exit(1)
      endif
 3    if(t.lt.rad_ptsqmin.or.t.lt.tmax) then
c below cut (either below absolute minimum, or below previously generated
c radiation in highest bid loop): generate a born event
         t=-1
         kn_csi=0
         return
      endif
c vetoes:
      rv=random()
      call set_rad_scales(t)
      if(kn_masses(kn_emitter).eq.0) then
         if(rad_iupperfsr.eq.1) then
            tmp=st_alpha / pwhg_alphas0(t,rad_lamll,nlc)
         elseif(rad_iupperfsr.eq.2) then
            tmp=st_alpha
         elseif(rad_iupperfsr.eq.3) then
            tmp=st_alpha
         endif
      else
         tmp=st_alpha
      endif
c Only for pp ->W, to account for em radiation from the electron
      if(flg_em_rad) then
c This should be equivalent at setting tmp=1
         tmp=tmp/st_alpha
      endif
      if(tmp.gt.1.000000001d0) then
         write(*,*) ' Error: upper bound lower than actual value',
     1        tmp,t
         call exit(1)
      endif
      if(rv.gt.tmp) then
         goto 1
      endif

      if(kn_masses(kn_emitter).eq.0) then
         if(rad_iupperfsr.eq.1) then         
c At this stage: pt generated according to (1) of upperbounding-fsr.pdf;
c generate csi uniformly in 1/csi
c in the range t/s < csi^2 < csimax^2
            rv=random()
            csi=exp(rv*log(t/s)/2+(1-rv)*log(kn_csimax))
c get y
            y=1-2*t/(s*csi**2)
c At this point a csi-y pair is generated according to the
c distribution upper(). It is automatically within range.
         elseif(rad_iupperfsr.eq.2) then
c     csi distributed uniformly in 1/(csi-t/s)
            rv=random()
            csi=1/(rv/(sqrt(t/s)-t/s)+(1-rv)/(kn_csimax-t/s))+t/s
c extra csi dependent factor
            if(random().gt.csi) goto 1
c get y
            y=1-2*t/(s*csi**2)
c At this point a csi-y pair is generated according to the
c distribution upper(). It is automatically within range,
c unless we have a massive emitter
         elseif(rad_iupperfsr.eq.3) then
c     csi distributed uniformly in 1/(csi-t/s)
            rv=random()
            csi=1/(rv/(sqrt(t/s)-t/s)+(1-rv)/(kn_csimax-t/s))+t/s
c get y
            y=1-2*t/(s*csi**2)
            if(random().gt.(csi-t/s)) goto 1
         else
            write(*,*) ' gen_rad_fsr:  rad_iupper=',rad_iupperfsr,
     1           ' invalid'
         endif
      else
c massive emitter case
         rv=random()
         call gencsiymv(t,rv,csi,y)
c Now veto if we are out of range
         if(csi.gt.1) goto 1
      endif
c
c extra suppression factor of upper bounding function (may depend upon radiation variables)
      call uboundfct(ufct,csi,y)
      if(random().gt.ufct) goto 1
c Veto from upper bound to real value. Count how many vetoes,
c since these may be expensive.
c      write(*,*) ' genrad_fsr: y and csi ',y,csi
      call sigborn_rad(born)
      if(born.lt.0) then
         born=0
      endif
      if(born.eq.0) then
c bizarre situation that may arise when the scale gets so low
c that some pdf vanish (typically heavy flavour pdf's)
         t=-1
         goto 3
      endif
      kn_y=y
      kn_csi=csi
      kn_azi=2*pi*random()
      ubound=born*pwhg_upperb_rad()*unorm*ufct
      call gen_real_phsp_fsr_rad
      call sigreal_rad(sig)
      value=sig*kn_jacreal
      if(value.gt.ubound) then
         call increasecnt(
     # 'upper bound failures in generation of radiation')
      endif
      rvalue=random()*ubound
      if(rvalue.gt.value) then
         call increasecnt('vetoed radiation')
         goto 1
      endif
      end


      subroutine add_azimuth
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_math.h'
      include 'pwhg_kn.h'
      integer ileg
      real * 8 azi,sazi,cazi
      real * 8 dir(3)
      data dir/0d0,0d0,1d0/
      real * 8 random
      external random
      azi=2d0*pi*random()
      sazi = sin(azi)
      cazi = cos(azi)
      if (kn_csi.eq.0d0) then
         do ileg=1, nlegborn
            call mrotate(dir,sazi,cazi,kn_pborn(1,ileg))
         enddo
      else
         do ileg=1, nlegreal
            call mrotate(dir,sazi,cazi,kn_preal(1,ileg))
         enddo
      endif
      end
      
      

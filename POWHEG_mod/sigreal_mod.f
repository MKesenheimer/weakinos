c MK: copied and modified version of sigreal.f, revision 3154
c modified on the basis of disquark/sigreal.f
c changes marked with "! CH, MK:" or "! MK:"
c bigger changes over a whole section are marked with !===...

      subroutine btildereal(xrad,resreal,www)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_flg.h'
      include 'pwhg_par.h'
      real * 8 xrad(3),resreal(maxprocborn),www
      real * 8 rr(maxalr),rc(maxalr),rp(maxalr),rm(maxalr),
     # rs(maxalr),rcs(maxalr),rps(maxalr),rms(maxalr),xl,xlp,xlm,
     # jac_over_csi,jac_over_csi_coll,jac_over_csi_soft,
     # jac_over_csi_p,jac_over_csi_m,rrr,rrrc,rrrs,rrrcs,
     # rrrp,rrrps,rrrm,rrrms,remnant,out0,out1
      real * 8, allocatable :: out1arr(:),out0arr(:)
      integer j,iuborn
      logical valid_emitter
      external valid_emitter
      logical pwhg_isfinite
      external pwhg_isfinite
      logical, save :: ini=.true.

      if(ini) then
         if(flg_analysisextrainfo) then
            allocate(out1arr(maxalr),out0arr(maxalr))
         endif
         ini = .false.
      endif

      do j=1,flst_nborn
         resreal(j)=0
      enddo
      do kn_emitter=0,nlegborn
c output values for analysis_driver
         out0=0
         out1=0
         if(flg_analysisextrainfo) then
            out0arr = 0
            out1arr = 0
         endif
c check that emitter is valid
         if(valid_emitter(kn_emitter)) then
            if(kn_emitter.gt.2) then
c     final state radiation
               call gen_real_phsp_fsr(xrad,jac_over_csi,
     #              jac_over_csi_coll,jac_over_csi_soft)
c This subroutine may set the scales with values depending
c upon the real emission kinematics
               call setscalesbtlreal
c sigreal fills the array rr with the value of the R_alpha contribution
c that have emitter equal to kn_emitter. All other contributions are set
c to zero. 
               call sigreal_btl(rr)
               if(flg_withsubtr) then
c We may prefer to set the counterterms scales different from the real scales
                  call setscalesbtlct
                  call collfsr(rc)
c     soft subtraction
                  call soft(rs)
                  call softcollfsr(rcs)
c     in final state radiation csimax is independent of y
                  xl=log(kn_csimax/par_csicut)
               endif
               do j=1,flst_nalr
                  iuborn=flst_alr2born(j)
                  rrr=rr(j)*kn_jacborn
     #                 *jac_over_csi/(1-kn_y)/kn_csitilde
                  if(flg_withsubtr) then
                     rrrc=rc(j)*kn_jacborn
     #                 *jac_over_csi_coll/(1-kn_y)/kn_csitilde
                     rrrs=rs(j)*kn_jacborn
     #                 *jac_over_csi_soft/(1-kn_y)/kn_csitilde
                     rrrcs=rcs(j)*kn_jacborn
     #                 *jac_over_csi_soft/(1-kn_y)/kn_csitilde
                     remnant=(rrrs-rrrcs)*xl*kn_csitilde
                  endif
                  if(flg_withsubtr) then
                     resreal(iuborn)= resreal(iuborn)+rrr-rrrc
     #-rrrs+rrrcs+remnant
                  else
c     provide a damping factor for the singular region,
c     to avoid divergent integral (25 is an ad hoc value
                     resreal(iuborn)= resreal(iuborn)
     #               +rrr*(1-kn_y**2)*kn_csi/
     #                  (25/kn_sbeams+(1-kn_y**2)*kn_csi)
                  endif
                  if(flg_nlotest) then
                     out1=out1+rrr
                     if(flg_analysisextrainfo) then
                        out1arr(j) = rrr
                     endif
                     if(flg_withsubtr) then
                        out0=out0-rrrc-rrrs+rrrcs+remnant
                        if(flg_analysisextrainfo) then
                           out0arr(j) = -rrrc-rrrs+rrrcs+remnant
                        endif
                     endif
                  endif
               enddo
            else
c     initial state singularities.
c     Regions that have only + (-) collinear singularity should return
c     zero rm (rp).
               call gen_real_phsp_isr
     #(xrad,jac_over_csi,jac_over_csi_p,jac_over_csi_m,
     #jac_over_csi_soft)
               call setscalesbtlreal
               call sigreal_btl(rr)
               if(flg_withsubtr) then
                  call setscalesbtlct
                  call soft(rs)
                  if(kn_emitter.ne.2) then
                     call collisrp(rp)
                     call softcollisrp(rps)
                  endif
                  if(kn_emitter.ne.1) then
                     call collisrm(rm)
                     call softcollisrm(rms)
                  endif
c     remnants (see xscaled.pdf in docs directory)
                  xl =log(kn_csimax/par_csicut)
                  xlp=log(kn_csimaxp/par_csicut)
                  xlm=log(kn_csimaxm/par_csicut)
               endif
               do j=1,flst_nalr
                  rrr=rr(j)*kn_jacborn
     #                 *jac_over_csi/(1-kn_y**2)/kn_csitilde
                  if(flg_withsubtr) then
                     rrrs=rs(j)*kn_jacborn
     #                 *jac_over_csi_soft/(1-kn_y**2)/kn_csitilde
                     remnant=rrrs*xl*kn_csitilde
                     if(kn_emitter.ne.2) then
                        rrrp=rp(j)*kn_jacborn
     #                 *jac_over_csi_p/(1-kn_y)/kn_csitilde/2
                        rrrps=rps(j)*kn_jacborn
     #                 *jac_over_csi_soft/(1-kn_y)/kn_csitilde/2
                        remnant=remnant-rrrps*xlp*kn_csitilde
                     else
                        rrrp=0
                        rrrps=0
                     endif
                     if(kn_emitter.ne.1) then
                        rrrm=rm(j)*kn_jacborn
     #                 *jac_over_csi_m/(1+kn_y)/kn_csitilde/2
                        rrrms=rms(j)*kn_jacborn
     #                 *jac_over_csi_soft/(1+kn_y)/kn_csitilde/2
                        remnant=remnant-rrrms*xlm*kn_csitilde
                     else
                        rrrm=0
                        rrrms=0
                     endif
                  endif
                  iuborn=flst_alr2born(j)
                  if(flg_withsubtr) then
                     resreal(iuborn)= resreal(iuborn)+rrr
     #              -rrrs-rrrp-rrrm+rrrps+rrrms+remnant
                  else
c     provide a damping factor for the singular region,
c     to avoid divergent integral (25 is an ad hoc value)
                     resreal(iuborn)= resreal(iuborn)
     #               +rrr*(1-kn_y**2)*kn_csi/
     #                  (25/kn_sbeams+(1-kn_y**2)*kn_csi)
                  endif
                  if(flg_nlotest) then
                     out1=out1+rrr
                     if(flg_analysisextrainfo) then
                        out1arr(j) = rrr
                     endif
                     if(flg_withsubtr) then
                        out0=out0-rrrs-rrrp-rrrm+rrrps+rrrms+remnant
                        if(flg_analysisextrainfo) then
                           out0arr(j) = -rrrs-rrrp-rrrm+rrrps+rrrms
     1                          +remnant
                        endif
                     endif
                  endif
               enddo
            endif
         endif
         if (.not.pwhg_isfinite(out0).or..not.pwhg_isfinite(out1)) then
            out0 = 0d0 
            out1 = 0d0 
            resreal = 0d0 
            if(flg_analysisextrainfo) then
               out1arr = 0
               out0arr = 0
            endif
         endif
         if(flg_nlotest) then
            out0=out0*www
            out1=out1*www
            call analysis_extrainfo('realct',flst_nalr,out0arr,www)
            if(out0.ne.0d0) call analysis_driver(out0,0)
            call analysis_extrainfo('real',flst_nalr,out1arr,www)
            if(out1.ne.0d0) call analysis_driver(out1,1)
         endif
      enddo
      end

      subroutine checklims(iun)
      implicit none
      integer iun
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_dbg.h'
      include 'Flags.h' ! MK: added
      external sigreal_btl,soft,collfsr,softcollfsr, collisrp,
     $     softcollisrp,collisrm,softcollisrm
      call randomsave
      flg_btilde=.true. ! MK: added
      if(dbg_softtest) then
         write(iun,*) '******************************************'   
         write(iun,*) '           CHECK  SOFT LIMITS             '     
         write(iun,*)
         do kn_emitter=0,nlegborn
            call checksoft(sigreal_btl,soft,' soft',iun)
         enddo
         write(iun,*) '******************************************'
         write(iun,*)
      endif
      
      if(dbg_colltest) then
         write(iun,*) '******************************************'   
         write(iun,*) '      CHECK  COLL. LIMITS FOR FSR       '         
         write(iun,*)
         do kn_emitter=3,nlegborn
            call checkcoll(sigreal_btl,collfsr,1,' coll',iun)
         enddo
         write(iun,*) '******************************************'   
         write(iun,*)
         write(iun,*) '******************************************'   
         write(iun,*) '      CHECK  COLL. LIMITS FOR ISR       '         
         write(iun,*)
         do kn_emitter=0,2
c            call randomsave
            if(kn_emitter.ne.2) call checkcoll(sigreal_btl,collisrp,1
     $           ,' coll-plus',iun)
c            call randomrestore
            if(kn_emitter.ne.1) call checkcoll(sigreal_btl,collisrm,-1
     $           ,' coll-minus',iun)
         enddo
         write(iun,*) '******************************************'   
         write(iun,*)
      endif
      
      if(dbg_softtest.and.dbg_colltest) then  
         write(iun,*) '******************************************'   
         write(iun,*) '   CHECK  SOFT-COLL. LIMITS FOR FSR     '         
         write(iun,*)
         do kn_emitter=3,nlegborn
            call checksoft(collfsr,softcollfsr,' soft-coll',iun)
         enddo
         write(iun,*) '******************************************'
         write(iun,*)
         write(iun,*) '******************************************'   
         write(iun,*) '   CHECK  SOFT-COLL. LIMITS FOR ISR +   '         
         write(iun,*)
         do kn_emitter=0,2
            if(kn_emitter.ne.2)call checksoft(collisrp,softcollisrp,
     $           ' soft-coll-plus',iun)
         enddo
         write(iun,*) '******************************************'
         write(iun,*)
         write(iun,*) '******************************************'   
         write(iun,*) '   CHECK  SOFT-COLL. LIMITS FOR ISR -   '         
         write(iun,*)
         do kn_emitter=0,2
            if(kn_emitter.ne.1)call checksoft(collisrm,softcollisrm,
     $           ' soft-coll-minus',iun)
         enddo
         write(iun,*) '******************************************'
         write(iun,*)

         
         write(iun,*) '******************************************'   
         write(iun,*) '   CHECK  COLL.-SOFT LIMITS FOR FSR     '         
         write(iun,*)
         do kn_emitter=3,nlegborn
            call checkcoll(soft,softcollfsr,1,' coll-soft',iun)
         enddo
         write(iun,*) '******************************************'
         write(iun,*)
         write(iun,*) '******************************************'   
         write(iun,*) '   CHECK  COLL.-SOFT LIMITS FOR ISR     '         
         write(iun,*)
         do kn_emitter=0,2
            if(kn_emitter.ne.2) call checkcoll(soft,softcollisrp,1
     $           ,' coll-plus-soft',iun)
            if(kn_emitter.ne.1) call checkcoll(soft,softcollisrm,-1
     $           ,' coll-minus-soft',iun)
         enddo
      endif
      call randomrestore
      end

      subroutine checkborn(iun)
c Check if Born, colour correlated born and spin correlated Born
c are consistent with total Born
      implicit none
      integer iun
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_br.h'
      include 'pwhg_st.h'
      integer  iborn,j,k,mu,kres,ires
      real * 8 tot
      real * 8 gtens(0:3,0:3),ap
      data gtens/1d0, 0d0, 0d0, 0d0,
     #           0d0,-1d0, 0d0, 0d0,
     #           0d0, 0d0,-1d0, 0d0,
     #           0d0, 0d0, 0d0,-1d0/
      save gtens
      logical colcorr
      external colcorr
      do ires=1,flst_nreson         
         do iborn=1,flst_nborn
            kres=flst_reslist(ires)
            do j=1,nlegborn
               if(colcorr(j,iborn,kres)) then
                  tot=0
                  do k=1,nlegborn
                     if(colcorr(k,iborn,kres)) then
                        if(k.ne.j) then
                           tot=tot+br_bornjk(j,k,iborn)
                        endif
                     endif
                  enddo
                  if(flst_born(j,iborn).eq.0) then
                     tot=tot/(ca*br_born(iborn))
                  else
                     tot=tot/(cf*br_born(iborn))
                  endif
                  if(abs((tot-1)/tot).gt.1d-8) then
                     write(iun,'(f6.3,a,20(i8,1x))') tot, ! MK: changed i3 -> i8
     1                    ' colour check fails for flav. struct:',kres,
     2                    (flst_born(k,iborn),k=1,nlegborn)
                  endif
               endif
            enddo
         enddo
      enddo
      do iborn=1,flst_nborn
         do j=1,nlegborn
            if(flst_born(j,iborn).eq.0) then
               tot=0
               do mu=0,3
                  tot=tot-gtens(mu,mu)*br_bmunu(mu,mu,j,iborn)
               enddo
               tot=tot/br_born(iborn)
               if(abs((tot-1)/tot).gt.1d-8) then
                  write(iun,'(f6.3,a,i2,a,20(i8,1x))') ! MK: changed i3 -> i8
     1    tot, ' spin correlated amplitude'//
     2 ' wrong for leg', j, ' flavour struct:',
     3                 (flst_born(k,iborn),k=1,nlegborn)
               endif
            endif
         enddo
      enddo
      end



      subroutine checksoft(sig,sigs,label,iun)
      implicit none
      include 'pwhg_dbg.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_st.h'
      character *(*) label
      integer iun
      real * 8 xborn(ndiminteg-3),xrad(3)
      integer nexp
      parameter (nexp=5)
      real * 8 jac_over_csi,
     #jac_over_csi_coll,jac_over_csi_soft,rr(maxalr,nexp),
     #rs(maxalr,nexp),jac_over_csi_p,jac_over_csi_m
      integer j,jexp,alr,alrp
      character * 15 flag
      character * 32 fff
      logical ident(maxalr)
      real * 8 random,dotp
      external random,dotp
      logical valid_emitter,iszero,isnonzero,isequal
      external valid_emitter
      do j=1,ndiminteg-3
         xborn(j)=random()
      enddo
      call gen_born_phsp(xborn)
      call setscalesbtilde
      call allborn
      call checkborn(iun)
c      write(iun,*)' mass',sqrt(2*dotp(kn_pborn(0,3),kn_pborn(0,4)))
      do j=1,3
         xrad(j)=random()
      enddo
c Check soft limits
      if(valid_emitter(kn_emitter)) then
         write(iun,*) ' Random Born variables ====> ',xborn
          write(iun,*) ' Random radiation variables ====> ',xrad
          do jexp=1,nexp
            xrad(1)=10d0**(-jexp)
            if(kn_emitter.gt.2) then
               call gen_real_phsp_fsr(xrad,jac_over_csi,
     $              jac_over_csi_coll,jac_over_csi_soft)
            else
               call gen_real_phsp_isr (xrad,jac_over_csi,jac_over_csi_p,
     $              jac_over_csi_m,jac_over_csi_soft)
            endif
            write(iun,*) '### Check soft',xrad(1)
            call sig(rr(1,jexp))
            call sigs(rs(1,jexp))
         enddo
         do alr=1,flst_nalr
            ident(alr)=.false.
         enddo
         do alr=1,flst_nalr
c only radiated gluons or photons
            if(flst_alr(nlegreal,alr).ne.0 .and.
     1         flst_alr(nlegreal,alr).ne.22  ) cycle
            if(ident(alr)) cycle
c     if one rr is zero, all others must be zero
            iszero=.false.
            isnonzero=.false.
            do jexp=1,nexp
c               if(rs(alr,jexp).ne.0) isnonzero=.true.
c               if(rs(alr,jexp).eq.0) iszero=.true.
               if(rr(alr,jexp).ne.0) isnonzero=.true.
               if(rr(alr,jexp).eq.0) iszero=.true.
            enddo
            if(iszero.and.isnonzero) then
               write(iun,*) ' some vanish and some do not'
            endif
            if(isnonzero.and..not.iszero) then
               ! MK: Changes for weakino made here (neutralino
               ! number 1000022 not correctly written out in realequiv
               ! file, i3 -> i8):
               fff = '(a,1x,i8,1x,a, 20(1x,i8),a,a,a)'
               write(fff(15:17),'(i3)') nlegreal
               write(iun,fff)
     $              ' emitter ',kn_emitter, ', process ',
     $              (flst_alr(j,alr),j=1,nlegreal)!,', ',label,':'
               do alrp=alr+1,flst_nalr
                  isequal=.true.
                  do jexp=1,nexp
                     if(rr(alr,jexp).ne.rr(alrp,jexp).or. rs(alr,jexp)
     $                    .ne.rs(alrp,jexp)) isequal=.false.
                  enddo
                  if(isequal) then
c                     write(iun,'(a,1x,i3,1x,a,20(1x,i3))')
c     $                    ' emitter ',kn_emitter, ', process ',
c     $                    (flst_alr(j,alrp),j=1,nlegreal) !,', ',label,':'
                     ident(alrp)=.true.
                  endif
               enddo
               do jexp=2,nexp
                  call setwarnflag
     1           (abs(rs(alr,jexp)/rr(alr,jexp)-1),jexp,2,flag)
                  write(iun,*) rs(alr,jexp)/rr(alr,jexp),flag
               enddo
            endif
         enddo
      endif
      end


      subroutine checkcoll(sig,sigc,idir,label,iun)
      implicit none
      integer iun
      include 'pwhg_dbg.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      character *(*) label
      integer idir
      real * 8 xborn(ndiminteg-3),xrad(3)
      integer nexp
      parameter (nexp=8)
      real * 8 jac_over_csi,
     #jac_over_csi_coll,jac_over_csi_soft,rr(maxalr,nexp),
     #rc(maxalr,nexp),jac_over_csi_p,jac_over_csi_m
      integer j,jexp,jexpfirst,alr,alrp
      real * 8 random
      external random
      logical ident(maxalr)
      character * 15 flag
      logical valid_emitter,iszero,isnonzero,isequal
      external valid_emitter
      do j=1,ndiminteg-3
         xborn(j)=random()
      enddo
      call gen_born_phsp(xborn)
      call setscalesbtilde
      call allborn
      do j=1,3
         xrad(j)=random()
      enddo
      if(valid_emitter(kn_emitter)) then
         write(iun,*) ' Random Born variables ====> ',xborn
         write(iun,*) ' Random radiation variables ====> ',xrad
         do jexp=1,nexp
            if(idir.ne.-1) then
               xrad(2)=10d0**(-jexp)
            else
               xrad(2)=1-10d0**(-jexp)
            endif
            if(kn_emitter.gt.2) then
               call gen_real_phsp_fsr(xrad,jac_over_csi,
     #jac_over_csi_coll,jac_over_csi_soft)
            else
               call gen_real_phsp_isr
     #(xrad,jac_over_csi,jac_over_csi_p,jac_over_csi_m,
     #jac_over_csi_soft)
            endif
            write(iun,*) '######### Check coll',xrad(2)
            call sig(rr(1,jexp))
            call sigc(rc(1,jexp))
         enddo
         do alr=1,flst_nalr
            ident(alr)=.false.
         enddo
         do alr=1,flst_nalr
            do jexp=1,nexp
               if(rc(alr,jexp).ne.rc(alr,1)) then
                  write(iun,*)
     #' checklims error : coll lim depends upon coll variable'
               endif
            enddo
         enddo
         do alr=1,flst_nalr
            if(ident(alr)) cycle
c     if one rr is zero, all others must be zero
            iszero=.false.
            isnonzero=.false.
            jexpfirst=2
            do jexp=nexp,1,-1
c               if(rc(alr,jexp).ne.0) isnonzero=.true.
c               if(rc(alr,jexp).eq.0) iszero=.true.
c               if(rr(alr,jexp).ne.0) isnonzero=.true.
               if(rr(alr,jexp).eq.0) then
                  if(rc(alr,jexp).ne.0) then
                     write(iun,*) ' some vanish and some do not'
                  endif
                  jexpfirst=jexp+1
                  goto 111
               endif
            enddo
 111        continue
            if(jexpfirst.le.nexp) then
               ! MK: Changes made in this line, i3 -> i8
               write(iun,'(2a,i2,a,20(1x,i8))') label,'  emitter ',
     #kn_emitter,', process ',(flst_alr(j,alr),j=1,nlegreal)
               do alrp=alr+1,flst_nalr
                  isequal=.true.
                  do jexp=1,nexp
                     if(rr(alr,jexp).ne.rr(alrp,jexp).or.
     #rc(alr,jexp).ne.rc(alrp,jexp)) isequal=.false.
                  enddo
                  if(isequal) then
c                     write(iun,'(2a,i2,a,20(1x,i3))') label,' emitter ',
c     # kn_emitter,', process ',(flst_alr(j,alrp),j=1,nlegreal)
                     ident(alrp)=.true.
                  endif
               enddo
               do jexp=jexpfirst,nexp
                  call setwarnflag(abs(rc(alr,jexp)/rr(alr,jexp)-1),
     1                 jexp,jexpfirst,flag)
c     Added this 'if' to be sure that no division by zero occurs
c     if((rr(alr,jexp-1)-rc(alr,jexp-1)).ne.0d0) then
c     write(iun,*) (rr(alr,jexp)-rc(alr,jexp))/
c     #(rr(alr,jexp-1)-rc(alr,jexp-1)),rc(alr,jexp)/rr(alr,jexp),flag
c     endif
                  write(iun,*) rc(alr,jexp)/rr(alr,jexp),flag
 1             enddo
            endif
         enddo
      endif
      end

      subroutine setwarnflag(dist,jexp,jexpfirst,flag)
      implicit none
      character * 15 flag
      real * 8 dist
      integer jexp,jexpfirst
      flag=' '
      if(dist.gt.0.01) then
         if(jexp.eq.jexpfirst) then
            if(dist.lt.0.1) then
               flag='*-WARN-*'
            else
               flag='*-WWARN-*'
            endif
         elseif(jexp.eq.3) then
            if(dist.lt.0.1) then
               flag='*-WWARN-*'
            else
               flag='*-WWWARN-*'
            endif
         elseif(jexp.ge.4) then
            if(dist.lt.0.1) then
               flag='*-WWWARN-*'
            elseif(dist.lt.0.3) then
               flag='*-WWWWARN-*'
            else
               flag='*-WWWWWARN-*'
            endif
         endif
      endif
      end



      subroutine sigreal_rad(sig)
      implicit none
      real * 8 sig
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      include 'pwhg_flg.h'
      include 'pwhg_par.h'
      include 'pwhg_pdf.h'
      include 'Flags.h' ! MK: added
      real * 8 rr(maxalr),rc(maxalr),rs(maxalr),rcs(maxalr)
      integer alr,alrpr,iret,em
      integer nmomset,emitter
      parameter (nmomset=10)
      real * 8 res(nmomset,maxalr),preal(0:3,nlegreal,nmomset),cprop
      integer equivto(maxalr)
      real * 8 equivcoef(maxalr)
      integer j,k
      real * 8 sumdijinv,dampfac,r
      real * 8 pdf1(-pdf_nparton:pdf_nparton),
     1         pdf2(-pdf_nparton:pdf_nparton)
      real * 8 ptsq,pwhg_pt2,dijterm
      logical computed(maxalr)
      logical condition
      logical ini
      data ini/.true./
      save ini,equivto,equivcoef
      external pwhg_pt2,dijterm
      flg_rad=.true. ! MK: added
      if(ini) then
         do alr=1,flst_nalr
            equivto(alr)=-1
         enddo
         if(flg_smartsig) then
            flg_in_smartsig = .true.
            call randomsave
c     generate "nmomset" random real-phase space configurations
            call fillmomenta(nlegreal,nmomset,kn_masses,preal)
            do alr=1,flst_nalr
               if(kn_emitter.eq.0) then
                  kn_resemitter=0
               else
                  kn_resemitter=flst_alrres(nlegreal,alr)
               endif
               flst_cur_alr=alr ! MK: added
               do j=1,nmomset
                  call realgr(
     1                 flst_alr(1,alr),preal(0,1,j),res(j,alr))
               enddo
               call compare_vecs(nmomset,alr,res,1,alrpr,cprop,iret)
               if(iret.eq.0) then
c     they are equal
                  equivto(alr)=alrpr
                  equivcoef(alr)=1
               elseif(iret.eq.1) then
c     they are proportional
                  equivto(alr)=alrpr
                  equivcoef(alr)=cprop
               else
c     < 0 for unequal:
                  equivto(alr)=-1
               endif
            enddo
            call randomrestore
         endif
         flg_in_smartsig = .false.
         ini=.false.
      endif
c End initialization phase; compute graphs
      do alr=1,flst_nalr
         rr(alr)=0
      enddo
      call pdfcall(1,kn_x1,pdf1)
      call pdfcall(2,kn_x2,pdf2)
      if(flg_withdamp) then
         call collrad(rc)
         call collsoftrad(rcs)
         call softrad(rs)
      endif
      do alr=1,flst_nalr
         computed(alr)=.false.
      enddo
      do j=1,rad_alr_nlist
         alr=rad_alr_list(j)
         em=flst_emitter(alr)
c check if emitter corresponds to current radiation region (i.e. rad_kinreg):
         if((rad_kinreg.eq.1.and.em.le.2).or.(em.gt.2.and.
     #       flst_lightpart+rad_kinreg-2.eq.em))then
c check if we have a g -> Q Qbar splitting below threshold:
            if(em.gt.0) then
               if(flst_alr(em,alr)+flst_alr(nlegreal,alr).eq.0.and.
     #abs(flst_alr(em,alr)).ge.4) then
                  ptsq=pwhg_pt2()
                  if(abs(flst_alr(em,alr)).eq.4
     #  .and.ptsq.lt.rad_charmthr2.or.
     # abs(flst_alr(em,alr)).eq.5.and.ptsq.lt.rad_bottomthr2) then
                     rr(alr)=0
                     goto 995
                  endif
               endif
            endif
c ----------------
c Gimnastic to avoid problem with non-lazy evaluation of logical
c expressions in gfortran; replaces the line
c            if(equivto(alr).lt.0.or..not.computed(equivto(alr))) then
            if(equivto(alr).lt.0) then
               condition=.true.
            elseif(.not.computed(equivto(alr))) then
               condition=.true.
            else
               condition=.false.
            endif
            if(condition) then
               if(kn_emitter.eq.0) then
                  kn_resemitter=0
               else
                  kn_resemitter=flst_alrres(nlegreal,alr)
               endif
               flst_cur_alr = alr
               call realgr(flst_alr(1,alr),kn_cmpreal,rr(alr))
               sumdijinv=0
               do k=1,flst_allreg(1,0,alr)
                  sumdijinv=sumdijinv
     #+1/dijterm(flst_allreg(1,k,alr),flst_allreg(2,k,alr),alr)
               enddo
               rr(alr)=rr(alr)/dijterm(em,nlegreal,alr)/sumdijinv
               if(em.gt.2) then
                  if(flg_doublefsr) then
c    supply a factor E_em/(E_em+E_rad), times 2 if both gluons
                     rr(alr)=rr(alr)
     1                    *kn_cmpreal(0,kn_emitter)**par_2gsupp/
     2                    (kn_cmpreal(0,kn_emitter)**par_2gsupp
     3                    +kn_cmpreal(0,nlegreal)**par_2gsupp)
                     if(flst_alr(kn_emitter,alr).eq.0.and.
     1                    flst_alr(nlegreal,alr).eq.0) then
                        rr(alr)=rr(alr)*2
                     endif
                  else
c     If the emitter is in the final state, and if the emitted and emitter
c     are both gluons, supply a factor E_em/(E_em+E_rad) * 2
                     if(flst_alr(kn_emitter,alr).eq.0.and.
     1                    flst_alr(nlegreal,alr).eq.0) then
                        rr(alr)=rr(alr)*2
     1                       *kn_cmpreal(0,kn_emitter)**par_2gsupp/
     2                       (kn_cmpreal(0,kn_emitter)**par_2gsupp
     3                       +kn_cmpreal(0,nlegreal)**par_2gsupp)
                     endif
                  endif
               endif
               rr(alr)=rr(alr)*flst_mult(alr)
c supply Born zero damping factor, if required
               if(flg_withdamp) then
                  r=rr(alr)
                  call bornzerodamp(alr,r,rc(alr),rs(alr),rcs(alr),
     1                 dampfac)
                  rr(alr)=rr(alr) * dampfac
               endif
               computed(alr)=.true.
               if(equivto(alr).gt.0) then
                  rr(equivto(alr))=rr(alr)/equivcoef(alr)
                  computed(equivto(alr))=.true.
               endif
            else
               rr(alr)=rr(equivto(alr))*equivcoef(alr)
            endif
         else
            rr(alr)=0
         endif
 995     continue
      enddo
      sig=0
      do j=1,rad_alr_nlist
         alr=rad_alr_list(j)
         if(rr(alr).ne.0) then
            rr(alr)=rr(alr)*pdf1(flst_alr(1,alr))*pdf2(flst_alr(2,alr))
            sig=sig+rr(alr)
            rad_real_arr(j)=rr(alr)
         else
            rad_real_arr(j)=0
         endif
      enddo
      end


      subroutine sigreal_btl(rr)
      implicit none
      real * 8 rr(*)
      call sigreal_btl0(rr,0)
      end

c Real cross section, required by btilde;
c fills the array rr(alr) with the invariant cross section, multiplied
c by csi^2 (1-y^2) for ISR regions
c    csi^2 (1-y)   for FSR regions
      subroutine sigreal_btl0(rr,imode)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_flg.h'
      include 'pwhg_par.h'
      include 'pwhg_pdf.h'
      include 'Flags.h' ! MK: added
      integer imode
      real * 8 rr(maxalr)
      real * 8 rc(maxalr),rs(maxalr),rcs(maxalr),r
      integer alr,alrpr,iret
      integer nmomset
      parameter (nmomset=10)
      real * 8 res(nmomset,maxalr),preal(0:3,nlegreal,nmomset),cprop
      integer equivto(maxalr),markused(maxalr)
      real * 8 equivcoef(maxalr)
      common/cequivtoreal/equivto
      common/cequivcoefreal/equivcoef
      integer j,k
      real * 8 sumdijinv,dampfac
      real * 8 pdf1(-pdf_nparton:pdf_nparton),
     1         pdf2(-pdf_nparton:pdf_nparton)
      real * 8 rescfac
      logical ini
      data ini/.true./
      save ini,/cequivtoreal/,/cequivcoefreal/
      real * 8 dijterm
      external dijterm
      flg_rad=.false. ! MK: added flags
      if(ini) then
         do alr=1,flst_nalr
            equivto(alr)=-1
         enddo
         if(flg_smartsig) then
            flg_in_smartsig = .true.
            call  printrealequiv         
            call randomsave
c     generate "nmomset" random real-phase space configurations
            call fillmomenta(nlegreal,nmomset,kn_masses,preal)
            do alr=1,flst_nalr
               do j=1,nmomset
                  if(kn_emitter.eq.0) then
                     kn_resemitter=0
                  else
                     kn_resemitter=flst_alrres(nlegreal,alr)
                     flst_cur_alr=alr ! MK: added
                  endif
                  call realgr(
     1                 flst_alr(1,alr),preal(0,1,j),res(j,alr))
               enddo
               call compare_vecs(nmomset,alr,res,0,alrpr,cprop,iret)
               if(iret.eq.0) then
c     they are equal:
                  equivto(alr)=alrpr
                  equivcoef(alr)=1
               elseif(iret.eq.1) then
c     they are proportional:
                  equivto(alr)=alrpr
                  equivcoef(alr)=cprop
               else
c     < 0 for unequal:
                  equivto(alr)=-1
               endif
            enddo
            call randomrestore
         endif
         call printrealequivregions
         flg_in_smartsig = .false.
         ini=.false.
      endif
c End initialization phase; compute graphs
      do alr=1,flst_nalr
         rr(alr)=0
         markused(alr)=0
      enddo
      if(flg_withdamp) then
         call collbtl(rc)
         call collsoftbtl(rcs)
         call softbtl(rs)
      endif
      do alr=1,flst_nalr
c Only R_alpha (namely alr) with the current emitter: 
         if(flst_emitter(alr).eq.kn_emitter) then
            if(equivto(alr).lt.0) then
c Not equal to any previous one, compute explicitly.
c First mark as being computed
               markused(alr)=1
               if(kn_emitter.eq.0) then
                  kn_resemitter=0
               else
                  kn_resemitter=flst_alrres(nlegreal,alr)
               endif
               flst_cur_alr = alr
               call realgr(flst_alr(1,alr),kn_preal,rr(alr))
c Supply FKS factor to separate singular region:
               sumdijinv=0
c Loop over all singular regions of the given contribution
               do k=1,flst_allreg(1,0,alr)
c flst_allreg({1,2},...) are the two legs that identify the k'th region
                  sumdijinv=sumdijinv
     1      +1/dijterm(flst_allreg(1,k,alr),flst_allreg(2,k,alr),alr)
               enddo
               rr(alr)=rr(alr)/dijterm(kn_emitter,nlegreal,alr)
     1              /sumdijinv
c If the emitter is in the final state, and if the emitted and emitter
c are both gluons, supply a factor E_em/(E_em+E_rad) * 2
               if(kn_emitter.gt.2) then
                  if(flg_doublefsr) then
c    supply a factor E_em/(E_em+E_rad), times 2 if both gluons
                     rr(alr)=rr(alr)
     1                    *kn_cmpreal(0,kn_emitter)**par_2gsupp/
     2                    (kn_cmpreal(0,kn_emitter)**par_2gsupp
     3                    +kn_cmpreal(0,nlegreal)**par_2gsupp)
                     if(flst_alr(kn_emitter,alr).eq.0.and.
     1                    flst_alr(nlegreal,alr).eq.0) then
                        rr(alr)=rr(alr)*2
                     endif
                  else
c     If the emitter is in the final state, and if the emitted and emitter
c     are both gluons, supply a factor E_em/(E_em+E_rad) * 2
                     if(flst_alr(kn_emitter,alr).eq.0.and.
     1                    flst_alr(nlegreal,alr).eq.0) then
                        rr(alr)=rr(alr)*2
     1                       *kn_cmpreal(0,kn_emitter)**par_2gsupp/
     2                       (kn_cmpreal(0,kn_emitter)**par_2gsupp
     3                       +kn_cmpreal(0,nlegreal)**par_2gsupp)
                     endif
                  endif
               endif
c supply Born zero damping factor, if required
               if(flg_withdamp) then
                  if(kn_emitter.gt.2) then
                     r=rr(alr)*(1-kn_y)*kn_csi**2
                  else
                     r=rr(alr)*(1-kn_y**2)*kn_csi**2
                  endif
                  r=r*flst_mult(alr)
                  call bornzerodamp(alr,r,rc(alr),rs(alr),rcs(alr),
     1                 dampfac)
                  if(imode.eq.0) then
                     rr(alr) =rr(alr) * dampfac
                  elseif(imode.eq.1) then
                     rr(alr) =rr(alr) * (1-dampfac)
                  else
                     write(*,*) ' sigreal_btl0: improper call'
                  endif
               endif
            else
               if(markused(equivto(alr)).ne.1) then
                  write(*,*) ' error: sigreal_btl flg_smartsig bug'
                  call exit(1)
               endif
               rr(alr)=rr(equivto(alr))*equivcoef(alr)
               markused(alr)=1
            endif
         endif
      enddo
      if(.not.flg_minlo) then
         rescfac = 1
         call pdfcall(1,kn_x1,pdf1)
         call pdfcall(2,kn_x2,pdf2)
      endif
      do alr=1,flst_nalr
         if(flg_minlo) then
            flg_minlo_real=.true.
            call setlocalscales(flst_alr2born(alr),2,rescfac)
            flg_minlo_real=.false.
            call pdfcall(1,kn_x1,pdf1)
            call pdfcall(2,kn_x2,pdf2)
         endif
         rr(alr)=rr(alr)*flst_mult(alr)
         if(kn_emitter.gt.2) then
            rr(alr)=rr(alr)*(1-kn_y)*kn_csi**2
         else
            rr(alr)=rr(alr)*(1-kn_y**2)*kn_csi**2
         endif
c include pdf's
         rr(alr)=rr(alr)*pdf1(flst_alr(1,alr))*pdf2(flst_alr(2,alr))
     1        *rescfac
      enddo
      end





      subroutine fillmomenta(nparticles,nmomset,masses,p)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      integer nparticles,nmomset
      real * 8 p(0:3,nparticles,nmomset),masses(nparticles)
      integer mu,j,k,i
      real * 8 pl,en,scale,ptki
      real * 8 random
      external random
      real * 8 ptmin,pt,ptcut,beta, vec(3),plab(0:3,nparticles),
     $     pcm(0:3,nparticles)!,costh,modk,modi
      integer last,sign
      logical debug
      parameter (debug=.false.)
      real * 8 dotp
      external dotp

c     produce ONLY events with minimum pt larger than ptcut
      ptcut=sqrt(kn_sbeams)/((nparticles+2)*10)
      if(nparticles.eq.3) ptcut = 0d0
      if (debug)  write(*,*) 'ptcut min in fillmomenta = ',ptcut

c     average value of the generated momentum.
      scale = sqrt(kn_sbeams)/(nparticles+2)
      do j=1,nmomset
 1       ptmin=1d30
         last = (nparticles-2)*random()+1
         last=last+2
         do k=3,nparticles
            if (k.ne.last) then
               do mu=1,3
                  plab(mu,k)=scale*(random()-0.5d0)
               enddo
            endif
         enddo
         plab(3,last)=scale*(random()-0.5d0)
         plab(1,last)=0d0
         plab(2,last)=0d0
         do k=3,nparticles
            if (k.ne.last) then
               plab(1,last)=plab(1,last)-plab(1,k)
               plab(2,last)=plab(2,last)-plab(2,k)
            endif
         enddo
         pl=0d0
         en=0d0
         do k=3,nparticles
            plab(0,k)=sqrt(masses(k)**2+
     #         plab(1,k)**2+plab(2,k)**2+plab(3,k)**2)
            pl=pl+plab(3,k)
            en=en+plab(0,k)
         enddo
         plab(0,1)=(en+pl)/2
         plab(0,2)=(en-pl)/2
         plab(3,1)=plab(0,1)
         plab(3,2)=-plab(0,2)
         plab(1,1)=0d0
         plab(1,2)=0d0
         plab(2,1)=0d0
         plab(2,2)=0d0

c         do k=1,nparticles
c            do mu=0,3
c               plab(mu,k)=plab(mu,k)*scale
c            enddo
c         enddo
c     boost momenta in the center-of-mass frame

         beta=-(plab(0,1)-plab(0,2))/(plab(0,1)+plab(0,2))
         vec(1)=0d0
         vec(2)=0d0
         vec(3)=1d0
         call mboost(nparticles,vec,beta,plab(0,1),pcm(0,1))
         
c     compute the minimum pt wrt the beam axis
         do k=3,nparticles
            pt=sqrt(pcm(1,k)**2+pcm(2,k)**2)
            ptmin=min(ptmin,pt)
         enddo
         if (ptmin.lt.ptcut) then
c            write(*,*) 'ptmin = ',ptmin            
            goto 1
         endif
c     compute the minimum pt among each pair of final-state momenta
         do k=3,nparticles-1
            do i=k+1,nparticles
c               modk=sqrt(pcm(1,k)**2+pcm(2,k)**2+pcm(3,k)**2)
c               modi=sqrt(pcm(1,i)**2+pcm(2,i)**2+pcm(3,i)**2)
c               costh=
c     $          (pcm(1,k)*pcm(1,i)+pcm(2,k)*pcm(2,i)+pcm(3,k)*pcm(3,i))/
c     $              modk/modi
c               ptmin=min(ptmin,modk*sqrt(abs(1-costh**2)),
c     $                         modi*sqrt(abs(1-costh**2)))
               ptki = sqrt(abs(2*dotp(pcm(0,k),pcm(0,i))*
     $              pcm(0,k)*pcm(0,i)/(pcm(0,k)+pcm(0,i))**2))               
               ptmin=min(ptmin,ptki)
            enddo
         enddo
         if (ptmin.lt.ptcut) then
c            write(*,*) 'ptmin = ',ptmin            
            goto 1
         endif

         if (debug) then
            write(*,*) 'set of momenta # ',j
            do k=1,nparticles
               write(*,'(10d12.4)') pcm(:,k), 
     $              sqrt(abs(pcm(0,k)**2-pcm(1,k)**2-pcm(2,k)**2-
     $              pcm(3,k)**2))
            enddo
            write(*,*) 'last particle ',last
         endif         
         p(:,:,j)=pcm(:,:)
      enddo
      end
 

      subroutine compare_vecs(nmomset,alr,res,imode,alrpr,cprop,iret)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      real * 8 ep
      parameter (ep=1d-8)
      integer nmomset,alr,alrpr,imode,iret,j,k
      real * 8 res(nmomset,*),cprop,rat
c imode=0 when called from btilde,
c imode=1 when called for radiation. In the latter
c case, graphs that do not have the same underlying Born
c are not considered.
      do j=1,alr-1
         if(flst_emitter(j).ne.flst_emitter(alr)) goto 10
         if(imode.eq.1.and.flst_alr2born(j).ne.flst_alr2born(alr))
     1        goto 10
         rat=res(1,alr)/res(1,j)
         do k=1,nmomset
            if(abs(1-res(k,alr)/res(k,j)/rat).gt.ep) goto 10
         enddo
         if(abs(1-rat).lt.ep) then
            iret=0
            cprop=1
         else
            iret=1
            cprop=rat
         endif
         alrpr=j
         return
 10      continue
      enddo
      iret=-1
      end

      subroutine realgr(rflav,p,res)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer rflav(nlegreal)
      real * 8 p(0:3,nlegreal),res
      call real_ampsq(p,rflav,res)
c flux factor
      res=res/(8*p(0,1)*p(0,2))
      if(res.eq.0) then
c         write(*,*) 'realgr:', rflav
         continue
      endif
      end



      subroutine real_ampsq(p,rflav,amp2)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      real * 8 p(0:3,nlegreal)
      integer rflav(nlegreal)
      real * 8 amp2 
      logical pwhg_isfinite
      external pwhg_isfinite
      call setreal(p,rflav,amp2)
c     check if amp2 is finite
      if (.not.pwhg_isfinite(amp2)) amp2=0d0
      amp2 = amp2*st_alpha/(2*pi)
      end






      subroutine printrealequivregions
c it prints the set of equivalent alr regions
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer equivto(maxalr)
      common/cequivtoreal/equivto
      real * 8 equivcoef(maxalr)
      common/cequivcoefreal/equivcoef
      integer j,k,iun,count
      save count
      data count/0/
      call newunit(iun)
      write(*,*) 'Writing realequivregions file...'
      open(unit=iun,file='realequivregions',status='unknown')
      do j=1,flst_nalr
         if(equivto(j).eq.-1) then
            write(iun,'(a)')
     1           'Beginning sequence of equivalent amplitudes'
            write(iun,100) 1d0,j, flst_alr(:,j)
            do k=1,flst_nalr
               if(equivto(k).eq.j) then
                  write(iun,100) equivcoef(k),k,flst_alr(:,k)
               endif
            enddo
            count=count+1
         endif
      enddo
      write(iun,*) ''
      write(iun,'(a,i4,a)') 'Found ',count, ' equivalent groups'
      close(iun)
      write(*,*) 'Done'
c MK: Changes for dineutral_minimal made here (neutralino number 1000022
c not correctly written out in bornequiv file, i4 -> i7):
 100  format(d11.4,5x,i8,5x,100(i8,1x))
      end

      subroutine printrealequiv
c it prints the set of equivalent real configurations
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_kn.h'
      integer nmomset
      parameter (nmomset=10)
      real * 8 res(nmomset,maxprocreal),preal(0:3,nlegreal,nmomset)
      real * 8 cprop
      integer equivto(maxprocreal)
      real * 8 equivcoef(maxprocreal)
      Integer j,k,iun,count
      integer ireal,irealpr,iret
      save count
      data count/0/      
      do ireal=1,flst_nreal
         equivto(ireal)=-1
      enddo
      call randomsave
c     generate "nmomset" random real-phase space configurations
      call fillmomenta(nlegreal,nmomset,kn_masses,preal)
      do ireal=1,flst_nreal
         do j=1,nmomset            
            call setreal(preal(0,1,j),flst_real(1,ireal),res(j,ireal))
         enddo
         call compare_vecs(nmomset,ireal,res,0,irealpr,cprop,iret)
         if(iret.eq.0) then
c     they are equal:
            equivto(ireal)=irealpr
            equivcoef(ireal)=1
         elseif(iret.eq.1) then
c     they are proportional:
            equivto(ireal)=irealpr
            equivcoef(ireal)=cprop
         else
c     < 0 for unequal:
            equivto(ireal)=-1
         endif
      enddo
      call randomrestore

      call newunit(iun)
      open(unit=iun,file='realequiv',status='unknown')
      write(*,*) 'Writing realequiv file...'
      do j=1,flst_nreal
         if(equivto(j).eq.-1) then
            write(iun,'(a)')
     1           'Beginning sequence of equivalent amplitudes'
            write(iun,100) 1d0,j, flst_real(:,j)
            do k=1,flst_nreal
               if(equivto(k).eq.j) then
                  write(iun,100) equivcoef(k),k,flst_real(:,k)
               endif
            enddo
            count=count+1
         endif
      enddo
      write(iun,*) ''
      write(iun,'(a,i4,a)') 'Found ',count, ' equivalent groups'
      close(iun)
      write(*,*) 'Done'
      ! MK: Changes for dineutral_minimal made here (neutralino number 
      ! 1000022 not correctly written out in bornequiv file, i4 -> i7):
 100  format(d11.4,5x,i8,5x,100(i8,1x))
      end




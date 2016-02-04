c MK: copied and modified version of sigremnants.f, revision 3154
c modified on the basis of disquark/sigremnants.f
c changes marked with "! CH, MK:" or "! MK:"
c bigger changes over a whole section are marked with ! =====...


c Function to be integrated with mint, used to generate
c non-singular contributions using standard mint-gen method
c These contributions arise from real graphs without singular regions,
c or, when damping of R/B value is adopted, from the remnant of the
c damping
      function sigremnant(xx,ww,ifirst,imode,retval,retval0)
c retval is the function return value
c retvavl0 is an 'avatar' function the has similar value, but is much
c easier to compute (i.e. the Born term in this case)
c imode = 0 compute retval0 only.
c imode = 1 compute retval, retval0
c return value: output, 0: success; 1: retval0 was not computed
c                 (this function does not support an avatar function)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      include 'pwhg_flg.h'
      include 'pwhg_math.h'
c MK: added
#include "osres.h"
#include "pwhg_flst_add.h"
#include "pwhg_rad_add.h"
#include "Flags.h"
      integer sigremnant,imode
      real * 8 retval,retval0,xx(ndiminteg),ww
      integer ifirst
      real * 8 xrad(3)
      real * 8 xborn(ndiminteg-3)
      integer j,alr
      real * 8 ttt
      real * 8 jac_over_csi,jac_over_csi_p,jac_over_csi_m,
     #         jac_over_csi_s,jac_over_csi_coll
      real * 8 xjac,suppfact,totreg,totrem
      logical valid_emitter
      external valid_emitter
      logical pwhg_isfinite
      external pwhg_isfinite
      sigremnant = 1
      if(ifirst.eq.2) then
         !retval=rad_reg_tot+rad_damp_rem_tot
         ! CH: changed: calculate the rad_... values
         ! in seperate routine, sim. to btilde, get sigremnant 
         ! directly from there
         call addupweightsrem(retval)
         if(flg_nlotest) call pwhgaccumup
         return
      endif
      do j=1,ndiminteg-3
         xborn(j)=xx(j)
      enddo
      do j=1,3
         xrad(j)=xx(ndiminteg-3 + j)
         rad_xradremn(j)=xrad(j)
      enddo
c regular contributions; any phase space parametrization should be OK
      kn_emitter=0
      call gen_born_phsp(xborn)
c set scales
      call setscalesbtilde
c the following is needed to compute soft and collinear limits
      call allborn
      call born_suppression(suppfact)
      if(flg_withreg.or.(flg_withdamp.and.(
     # valid_emitter(0).or.valid_emitter(1).or.valid_emitter(2)))) then
         call gen_real_phsp_isr(xrad,
     #    jac_over_csi,jac_over_csi_p,jac_over_csi_m,jac_over_csi_s)
         xjac=jac_over_csi*kn_csi*kn_csimax*kn_jacborn*ww*hc2
      endif
      totreg = 0d0 ! CH: added
      if(flg_withreg) then
c This subroutine may set the scales with values depending
c upon the real emission kinematics
         call setscalesbtlreal
         call sigreal_reg(xjac,totreg,rad_reg_arr) ! CH: replaced rad_reg_tot->totreg
         ! CH: added the following lines
         ! =============================================================
         if (.not.pwhg_isfinite(totreg)) then 
            totreg = 0d0 
            retval = 0d0
         endif
         call transfersign(rad_reg_arr,rad_reg_sign,flst_nregular)
         ! =============================================================
         if(flg_nlotest) then
            call analysis_extrainfo('reg',flst_nregular,rad_reg_arr,1d0)
            call analysis_driver(totreg/suppfact,1) ! CH: changed
         endif
c           else
c              rad_reg_tot=0
      endif

      totrem=0d0 ! CH: added
      if(flg_withdamp) then
c              rad_damp_rem_tot=0 !CH: commented out (removed this entry from pwhg_rad.h)
         do alr=1,flst_nalr
            rad_damp_rem_arr(alr)=0
         enddo
         do kn_emitter=0,nlegborn
            if(valid_emitter(kn_emitter)) then
               if(kn_emitter.le.2) then
c     No need to generate phase space; it is already available
                  call setscalesbtlreal
                  call sigreal_damp_rem(xjac,ttt,rad_damp_rem_arr)
                  ! CH: added the following lines
                  ! ====================================================
                  if (.not.pwhg_isfinite(ttt)) then 
                     ttt = 0d0 
                     retval = 0d0
                  endif
                  ! ====================================================
                  if(flg_nlotest) then
                     call analysis_extrainfo('remn',
     1                    flst_nalr,rad_damp_rem_arr,1d0)
                     call analysis_driver(ttt,1)
                  endif
c                       rad_damp_rem_tot=rad_damp_rem_tot+ttt
               else
                  call gen_real_phsp_fsr(xrad,
     #jac_over_csi,jac_over_csi_coll,jac_over_csi_s)
                  xjac=jac_over_csi*kn_csi*kn_csimax
     #                *kn_jacborn*ww*hc2
                  call setscalesbtlreal
                  call sigreal_damp_rem(xjac,ttt,rad_damp_rem_arr)
                  ! CH: added the following lines
                  ! ====================================================
                  if (.not.pwhg_isfinite(ttt)) then 
                     ttt = 0d0 
                     retval = 0d0
                  endif
                  ! ====================================================
                  if(flg_nlotest) then
                     call analysis_extrainfo('remn',
     1                    flst_nalr,rad_damp_rem_arr,1d0)
                     call analysis_driver(ttt,1)
                  endif
c                       rad_damp_rem_tot=rad_damp_rem_tot+ttt
               endif
               totrem=totrem+ttt ! CH: added
            endif
         enddo
         ! CH: added
         call transfersign(rad_damp_rem_arr,rad_damp_rem_sign,flst_nalr)
         ! CH: commented
c           else
c              rad_damp_rem_tot=0
      endif

c        if (.not.pwhg_isfinite(rad_reg_tot)) then
c           rad_reg_tot=0d0
c           rad_reg_arr=0d0
c        endif
c        if (.not.pwhg_isfinite(rad_damp_rem_tot)) then
c           rad_damp_rem_tot=0d0
c           rad_damp_rem_arr=0d0
c        endif
c        retval=rad_reg_tot+rad_damp_rem_tot
         retval=totrem+totreg

#ifdef DEBUGQ
         print*,"sigremnants.f:167: retval",retval
         stop
#endif
      end

      subroutine sigreal_reg(xjac,sig,r0)
c contributions from real graphs that do not have a singular region
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      include 'pwhg_flg.h'
      include 'pwhg_pdf.h'
c MK: added
#include "osres.h"
#include "pwhg_flst_add.h"
#include "pwhg_rad_add.h"
#include "Flags.h"
      real * 8 xjac,sig,r0(maxprocreal)
      integer lreg,lregpr,iret
      integer nmomset
      parameter (nmomset=10)
      real * 8 res(nmomset,maxprocreal),preal(0:3,nlegreal,nmomset),
     #   cprop
      integer equivto(maxprocreal)
      real * 8 equivcoef(maxprocreal)
      integer j
      real * 8 pdf1(-pdf_nparton:pdf_nparton),
     1         pdf2(-pdf_nparton:pdf_nparton)
      logical ini
      data ini/.true./
      save ini,equivto,equivcoef
      if(ini) then
         do lreg=1,flst_nregular
            equivto(lreg)=-1
         enddo
         if(flg_smartsig) then
            flg_in_smartsig = .true.
            call randomsave
c     generate "nmomset" random real-phase space configurations
            call fillmomenta(nlegreal,nmomset,kn_masses,preal)
            do lreg=1,flst_nregular
               do j=1,nmomset
                  call realgr(
     1                 flst_regular(1,lreg),preal(0,1,j),res(j,lreg))
               enddo
               call compare_vecs_reg(nmomset,lreg,res,lregpr,cprop,iret)
               if(iret.eq.0) then
c     they are equal
                  equivto(lreg)=lregpr
                  equivcoef(lreg)=1
               elseif(iret.eq.1) then
c     they are proportional
                  equivto(lreg)=lregpr
                  equivcoef(lreg)=cprop
               endif
            enddo
            call randomrestore
         endif
         flg_in_smartsig = .false.
         ini=.false.
      endif
c End initialization phase; compute graphs
      call pdfcall(1,kn_x1,pdf1)
      call pdfcall(2,kn_x2,pdf2)
      do lreg=1,flst_nregular
c ----------------
         if(equivto(lreg).lt.0) then
            flst_cur_alr = -1
            call realgr(flst_regular(1,lreg),kn_cmpreal,r0(lreg))
         else
            r0(lreg)=r0(equivto(lreg))*equivcoef(lreg)
         endif
      enddo
      sig=0
      do lreg=1,flst_nregular
         r0(lreg)=xjac*r0(lreg)*
     #      pdf1(flst_regular(1,lreg))*pdf2(flst_regular(2,lreg))
         sig=sig+r0(lreg)
      enddo
      end

      subroutine compare_vecs_reg(nmomset,lreg,res,lregpr,cprop,iret)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
c MK: added
#include "osres.h"
#include "pwhg_flst_add.h"
      real * 8 ep
      parameter (ep=1d-12)
      integer nmomset,lreg,lregpr,iret,j,k
      real * 8 res(nmomset,*),cprop,rat
      do j=1,lreg-1
         rat=res(1,lreg)/res(1,j)
         do k=1,nmomset
            if(abs(1-res(k,lreg)/res(k,j)/rat).gt.ep) goto 10
         enddo
         if(abs(1-rat).lt.ep) then
            iret=0
            cprop=1
         else
            iret=1
            cprop=rat
         endif
         lregpr=j
         return
 10      continue
      enddo
      iret=-1
      end

      subroutine sigreal_damp_rem(xjac,sig,r1)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
c MK: added
#include "osres.h"
#include "pwhg_flst_add.h"
      real * 8 xjac,sig
      real * 8 r0(maxalr),r1(maxalr)
      integer alr
      sig=0
      call sigreal_btl0(r0,1)
      do alr=1,flst_nalr
         if(kn_emitter.eq.flst_emitter(alr)) then
            if(kn_emitter.le.2) then
               r0(alr)=r0(alr)/((1-kn_y**2)*kn_csi**2)
            else
               r0(alr)=r0(alr)/((1-kn_y)*kn_csi**2)
            endif
            r0(alr)=r0(alr)*xjac
            r1(alr)=r1(alr)+r0(alr)
            sig=sig+r0(alr)
         endif
      enddo
      end

c the following routines are similar to the btilde-case
       subroutine addupweightsrem(sigremnant)
       implicit none
       include 'nlegborn.h'
       include 'pwhg_flst.h'
       include 'pwhg_rad.h'
c MK: added
#include "osres.h"
#include "pwhg_flst_add.h"
#include "pwhg_rad_add.h"
       real*8 sigremnant

       real*8 dtotreg,dtotabsreg,dtotposreg,dtotnegreg,
     &       dtotrem,dtotabsrem,dtotposrem,dtotnegrem
      
       real*8 totreg,etotreg,totabsreg,etotabsreg,totposreg,etotposreg,
     &       totnegreg,etotnegreg
       real*8 totrem,etotrem,totabsrem,etotabsrem,totposrem,etotposrem,
     &       totnegrem,etotnegrem
       integer ncall,i

       common/cadduptotalsrem/totreg,etotreg,totabsreg,etotabsreg,
     &             totposreg,etotposreg,totnegreg,etotnegreg,
     &             totrem,etotrem,totabsrem,etotabsrem,
     &             totposrem,etotposrem,totnegrem,etotnegrem,
     &             ncall
      
       ncall      = ncall + 1
       dtotreg    = 0
       dtotabsreg = 0
       dtotposreg = 0
       dtotnegreg = 0
       dtotrem    = 0
       dtotabsrem = 0
       dtotposrem = 0
       dtotnegrem = 0

       do i=1,flst_nregular
         dtotreg    = dtotreg + rad_reg_arr(i)*rad_reg_sign(i)
         dtotabsreg = dtotabsreg + rad_reg_arr(i)
         if(rad_reg_sign(i).eq.1) then
           dtotposreg = dtotposreg+rad_reg_arr(i)
         else
           dtotnegreg = dtotnegreg+rad_reg_arr(i)
         endif
       enddo

       do i=1,flst_nalr
         dtotrem    = dtotrem + rad_damp_rem_arr(i)*rad_damp_rem_sign(i)
         dtotabsrem = dtotabsrem+rad_damp_rem_arr(i)
         if(rad_damp_rem_sign(i).eq.1) then
           dtotposrem = dtotposrem+rad_damp_rem_arr(i)
         else
           dtotnegrem = dtotnegrem+rad_damp_rem_arr(i)
         endif
       enddo

       totreg     = totreg     + dtotreg
       etotreg    = etotreg    + dtotreg**2
       totabsreg  = totabsreg  + dtotabsreg
       etotabsreg = etotabsreg + dtotabsreg**2
       totposreg  = totposreg  + dtotposreg
       etotposreg = etotposreg + dtotposreg**2
       totnegreg  = totnegreg  + dtotnegreg
       etotnegreg = etotnegreg + dtotnegreg**2
       totrem     = totrem     + dtotrem
       etotrem    = etotrem    + dtotrem**2
       totabsrem  = totabsrem  + dtotabsrem
       etotabsrem = etotabsrem + dtotabsrem**2
       totposrem  = totposrem  + dtotposrem
       etotposrem = etotposrem + dtotposrem**2
       totnegrem  = totnegrem  + dtotnegrem
       etotnegrem = etotnegrem + dtotnegrem**2
       
       sigremnant = dtotabsreg+dtotabsrem
       end

c   set all new totals concerning regulars/remnants/split to 0
       subroutine resettotalsrem
       implicit none
       real*8 totreg,etotreg,totabsreg,etotabsreg,totposreg,etotposreg,
     &       totnegreg,etotnegreg
       real*8 totrem,etotrem,totabsrem,etotabsrem,totposrem,etotposrem,
     &       totnegrem,etotnegrem
       integer ncall

       common/cadduptotalsrem/totreg,etotreg,totabsreg,etotabsreg,
     &             totposreg,etotposreg,totnegreg,etotnegreg,
     &             totrem,etotrem,totabsrem,etotabsrem,
     &             totposrem,etotposrem,totnegrem,etotnegrem,
     &             ncall

       ncall      = 0
       totreg     = 0d0
       etotreg    = 0d0
       totabsreg  = 0d0
       etotabsreg = 0d0
       totposreg  = 0d0
       etotposreg = 0d0
       totnegreg  = 0d0
       etotnegreg = 0d0
       totrem     = 0d0
       etotrem    = 0d0
       totabsrem  = 0d0
       etotabsrem = 0d0
       totposrem  = 0d0
       etotposrem = 0d0
       totnegrem  = 0d0
       etotnegrem = 0d0
       end

c   similar to corresponding routine in btilde
       subroutine finaltotalsrem
       implicit none
       include 'nlegborn.h'
       include 'pwhg_flst.h'
       include 'pwhg_rad.h'
c MK: added
#include "osres.h"
#include "pwhg_flst_add.h"
#include "pwhg_rad_add.h"
       real*8 totreg,etotreg,totabsreg,etotabsreg,totposreg,etotposreg,
     &       totnegreg,etotnegreg
       real*8 totrem,etotrem,totabsrem,etotabsrem,totposrem,etotposrem,
     &       totnegrem,etotnegrem
       integer ncall
       real*8 calc_error
       external calc_error
       common/cadduptotalsrem/totreg,etotreg,totabsreg,etotabsreg,
     &             totposreg,etotposreg,totnegreg,etotnegreg,
     &             totrem,etotrem,totabsrem,etotabsrem,
     &             totposrem,etotposrem,totnegrem,etotnegrem,
     &             ncall
       if(ncall.eq.0) ncall = 1 ! if we never call the remnant-routines: avoid NaNs
       rad_totreg     = totreg/ncall
       rad_etotreg    = calc_error(totreg,etotreg,ncall)
       rad_totabsreg  = totabsreg/ncall
       rad_etotabsreg = calc_error(totabsreg,etotabsreg,ncall)
       rad_totposreg  = totposreg/ncall
       rad_etotposreg = calc_error(totposreg,etotposreg,ncall)
       rad_totnegreg  = totnegreg/ncall
       rad_etotnegreg = calc_error(totnegreg,etotnegreg,ncall)
       rad_totrem     = totrem/ncall
       rad_etotrem    = calc_error(totrem,etotrem,ncall)
       rad_totabsrem  = totabsrem/ncall
       rad_etotabsrem = calc_error(totabsrem,etotabsrem,ncall)
       rad_totposrem  = totposrem/ncall
       rad_etotposrem = calc_error(totposrem,etotposrem,ncall)
       rad_totnegrem  = totnegrem/ncall
       rad_etotnegrem = calc_error(totnegrem,etotnegrem,ncall)
       end

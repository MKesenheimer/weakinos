c MK: copied and modified version of gen_index.f, revision 3154
c modified on the basis of disquark/gen_index.f
c changes marked with "! CH, MK:" or "! MK:"
c bigger changes over a whole section are marked with !===...

      subroutine gen_uborn_idx
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      integer j,alr,em
      do j=1,rad_nkinreg
         rad_kinreg_on(j)=.false.
      enddo
      call pick_random(flst_nborn,rad_btilde_arr,rad_ubornidx)
      rad_alr_nlist=flst_born2alr(0,rad_ubornidx)
      do j=1,rad_alr_nlist
         alr=flst_born2alr(j,rad_ubornidx)
         em=flst_emitter(alr)
         rad_alr_list(j)=alr
         if(em.le.2) then
            rad_kinreg_on(1)=.true.
         else
            rad_kinreg_on(em-flst_lightpart+2)=.true.
         endif
      enddo
      end

      subroutine gen_real_idx
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      ! MK: picks a weighted random index according to real-xsec (rad_real_arr)
      ! of the resonant region and sets the appropriate emitter of this
      ! region
      call pick_random(rad_alr_nlist,rad_real_arr,rad_realidx)
      rad_realalr=rad_alr_list(rad_realidx)
      kn_emitter = flst_emitter(rad_realalr)
#ifdef DEBUGQ
      print*,"rad_real_arr",rad_real_arr(1:rad_alr_nlist)
      print*,"rad_alr_list",rad_alr_list(1:rad_alr_nlist)
      !print*,"rad_realidx",rad_realidx
      print*,"rad_realalr",rad_realalr
      print*,"flst_emitter",flst_emitter(1:rad_alr_nlist)
      print*,"kn_emitter",kn_emitter
      !print*,"maxalr",maxalr
      print*
#endif
      end

      subroutine pick_random(n,r,jret)
      implicit none
      integer n,jret
      real * 8 r(n)
      real * 8 tot(0:n)
      integer j
      real * 8 ran,random
      tot(0)=0
      do j=1,n
         tot(j)=tot(j-1)+r(j)
      enddo
      ran=tot(n)*random()
      do j=1,n
         if(tot(j).gt.ran) then
            jret=j
            return
         endif
      enddo
      write(*,*) ' ***********************************'
      write(*,*) ' POWHEGBOX:pick_random: could not pick a value'
      write(*,*) ' set choice to 1'
      write(*,*) ' ***********************************'
      jret=1
c      call exit(-1)
      end
      

      subroutine gen_remnant(iret)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
#include "osres.h"
#include "pwhg_flst_add.h"
#include "pwhg_rad_add.h"
      integer iret
      real * 8 dum1,dum2,dum3,dum4
      real * 8 random
      external random
c Remnant or regular?
c      if(random().gt.rad_reg_tot/(rad_reg_tot+rad_damp_rem_tot)) then
      ! CH: just modified the naming...
      if(random().gt.rad_totreggen/(rad_totreggen+rad_totremgen)) then
c remnant
         iret=1
         call pick_random(flst_nalr,rad_damp_rem_arr,rad_realalr)
         rad_ubornidx=flst_alr2born(rad_realalr)
         kn_emitter=flst_emitter(rad_realalr)
         if(kn_emitter.le.2) then
            call gen_real_phsp_isr(rad_xradremn,dum1,dum2,dum3,dum4)
         else
            call gen_real_phsp_fsr(rad_xradremn,dum1,dum2,dum3)
         endif
      else
c regular
         iret=2
         call pick_random(flst_nregular,rad_reg_arr,rad_realreg)
         call gen_real_phsp_isr(rad_xradremn,dum1,dum2,dum3,dum4)
      endif
      end
      
! MK: added the following routine for the osres contributions
      subroutine gen_osres(iret)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      include 'pwhg_flg.h'
c MK: added
#include "PhysPars.h"
#include "osres.h"
#include "pwhg_flst_add.h"
#include "pwhg_rad_add.h"
      integer iret
      real * 8 random, ran1
      external random
      integer ichan ! MK: added
      double precision rad_totosres_sum ! MK: added
      double precision rad_totosres_frac_sum ! MK: added

#ifdef DEBUGQ
      print*,"[DEBUG] in gen_index_mod:124"
      print*,"uncomment to continue"
      stop
#endif
      
      ! sum
      rad_totosres_sum = 0D0
      do ichan=1,nosres
        rad_totosres_sum = rad_totosres_sum + rad_totosres(ichan)
      enddo
      
      rad_totosres_frac_sum = 0D0
      
      ! We want to switch randomly between all on-shell resonant
      ! contributions.
      ! So we go through all entries of the array rad_totosres
      ! and sum them up step by step. Then we check:
      ! random()<=rad_totosres(1)/rad_totosres_sum? -> iret = 3
      ! random()<=(rad_totosres(1)+rad_totosres(2))/rad_totosres_sum? -> iret = 4
      ! ...
      ! random()<=rad_totosres_sum/rad_totosres_sum? -> iret = 2+nosres
      !
      ! figuratively, if ran1 lies e.g. in intervall 3:
      !
      ! contr.   1    2     3   4   5    ... nosres
      ! sum    |---|-----|----|--|-----| ...
      ! step1: |---|                                 -> if(.false.) -> next step
      ! step2: |---|-----|                           -> if(.false.) -> next step
      ! step3: |---|-----|----|                      -> if(.true.)  -> iret = 2+3 = 5
      ! ran1               ^
      !  
      !
      ! Note: since genremnant uses indices iret=1,2 we have to shift
      ! the indices for the on-shell resonances by 2.
      ! => iret = ichan + 2 
      
c which version is better?
#define GENOSRES1
#ifdef GENOSRES1
      ran1 = random()
      do ichan=1,nosres  
        rad_totosres_frac_sum = rad_totosres_frac_sum
     &                           + rad_totosres(ichan)
        if(ran1.le.rad_totosres_frac_sum/rad_totosres_sum) then
          iret=2+ichan
          call pick_random(flst_nosres,rad_osres_arr(:,ichan),
     &                                                  rad_realosres)
          call real_osres_phsp(rad_xradosres,ichan)
          rad_ubornidx=flst_alr2born(rad_realosres)
          kn_emitter=flst_emitter(rad_realosres)
          goto 10
        endif
      enddo
 10   continue
#endif

#ifdef GENOSRES2
      ran1 = random()
      rad_totosres_frac_sum = rad_totosres_sum
      do ichan=1,nosres  
        rad_totosres_frac_sum = rad_totosres_frac_sum
     &                           - rad_totosres(ichan)
        if(ran1.gt.rad_totosres_frac_sum/rad_totosres_sum) then
          iret=2+ichan
          call pick_random(flst_nosres,rad_osres_arr(:,ichan),
     &                                                  rad_realosres)
          call real_osres_phsp(rad_xradosres,ichan)
          rad_ubornidx=flst_alr2born(rad_realosres)
          kn_emitter=flst_emitter(rad_realosres)
          goto 10
        endif
      enddo
 10   continue
#endif
      end

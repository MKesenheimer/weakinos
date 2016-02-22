c############### sigosres.f ############################################
c last modified by MK, 05.12.2015
c weakino pair productions
c new file on the basis of sigremnants.f, revision 3154
c changes not marked

c############### function sigosres #####################################
c Function to be integrated with mint, used to generate
c non-singular contributions using standard mint-gen method
c These contributions arise from real graphs without singular regions,
c or, when damping of R/B value is adopted, from the remnant of the
c damping
c retval is the function return value
c retvavl0 is an 'avatar' function the has similar value, but is much
c easier to compute (i.e. the Born term in this case)
c TODO:
c imode = 0 compute retval0 only.
c imode = 1 compute retval, retval0
c return value: output, 0: success; 1: retval0 was not computed
c                   (this function does not support an avatar function)
      function sigosres(xx,ww1,ifirst,imode,retval,retval0)
        implicit none

#include "nlegborn.h"
#include "pwhg_flst.h"
#include "pwhg_kn.h"
#include "pwhg_rad.h"
#include "pwhg_flg.h"
#include "pwhg_math.h"
#include "PhysPars.h"

c keep this order
#include "osres.h"
#include "pwhg_flst_add.h"
#include "pwhg_rad_add.h"

        integer sigosres,imode
        double precision retval,retval0,xx(ndiminteg),ww1
        double precision yy(ndiminteg)
        integer ifirst,ichan,idi
        double precision xjac
        ! temporary results for calc. of resonant contributions
        double precision sigosres_contr
        
        ! -> return succes
        sigosres = 1
        
        ! set this to zero: radiation is always generated from production 
        ! process for osres contributions
        kn_resemitter = 0
     
        if(ifirst.eq.2) then
          call addupweightsosres(retval)
          !retval = rad_reg_tot+rad_damp_rem_tot
          if(flg_nlotest) call pwhgaccumup
          return
        endif

        call setscalesbtilde
        ! we need the last random-numbers for the radiation-generation
        do idi=1,ndiminteg
          rad_xradosres(idi) = xx(idi)
        enddo
        
        ! reset result
        retval = 0D0

        ! use the generic phase-space here,
        ! provide tan-mapping for the resonant particles
        ! sum over the resonances
        do ichan=1,nosres
          ! first create a momentum-config where the particles i and j are
          ! tan-mapped on a resonant squark
          ! phase space that builds the 2->3 PS by using only 1->2 sub PS
          call real_osres_phsp(xx,ichan)
          xjac = kn_jacreal*ww1*hc2
          call sigreal_osres(xjac,sigosres_contr,
     &                      rad_osres_arr(:,ichan),ichan)         
          call transfersign(rad_osres_arr(:,ichan),
     &                      rad_osres_sign(:,ichan),flst_nosres)
          if(flg_nlotest) then
            call analysis_driver(sigosres_contr,1)
          endif
          ! TODO: ausprobieren:
          retval = retval + dabs(sigosres_contr)
          ! old:
          !retval = retval - sigosres_contr
        enddo

        ! retval0 = 0D0
        
#ifdef DEBUGQ
         print*,"retval",retval
         print*,"sigosres_contr",sigosres_contr
         print*,"xjac",xjac
         print*
         !stop
#endif
      end
c############### end function sigosres #################################

c############### subroutine sigreal_osres ##############################
c contributions from real graphs that do not have a singular region
      subroutine sigreal_osres(xjac,sig,r0,ichan)
        implicit none
        
#include "nlegborn.h"
#include "pwhg_flst.h"
#include "pwhg_kn.h"
#include "pwhg_rad.h"
#include "pwhg_flg.h"
#include "pwhg_pdf.h"

c keep this order
#include "osres.h"
#include "pwhg_flst_add.h"
#include "pwhg_rad_add.h"

        double precision xjac,sig,r0(maxprocreal)
        integer lset,lsetpr,iret
        integer nmomset
        parameter (nmomset=10)
        double precision res(nmomset,maxprocreal)
        double precision preal(0:3,nlegreal,nmomset)
        double precision cprop
        integer j
        double precision pdf1(-pdf_nparton:pdf_nparton)
        double precision pdf2(-pdf_nparton:pdf_nparton)
        ! we need to initialize this routine for every on-shell resonance
        integer equivto(maxprocreal,cnosres)
        double precision equivcoef(maxprocreal,cnosres)
        logical init(cnosres)
        data init/cnosres*.true./
        save equivto,equivcoef
        save init
        integer ichan

        ! first check if the jacobian is 0
        ! if yes return immediately and save time
        if(xjac.eq.0D0) then 
          sig = 0
          do lset=1,flst_nosres
            r0(lset) = 0D0
          enddo
          return
        endif

        ! initialization phase
        if(init(ichan)) then
          do lset=1,flst_nosres
            equivto(lset,ichan) = -1
          enddo
          flg_in_smartsig = .true. ! causes trouble with subtraction
          if(flg_smartsig) then
            call randomsave
            ! generate "nmomset" random real-phase space configurations
            call fillmomenta(nlegreal,nmomset,kn_masses,preal)
            do lset=1,flst_nosres
              do j=1,nmomset
                ! call realgr(flst_osres(1,lset),preal(0,1,j),res(j,lset))
                call setosresreal(preal(0,1,j),flst_osres(1,lset),
     &                            osresID(ichan),res(j,lset))
              enddo
              
              call compare_vecs_osres(nmomset,lset,res,lsetpr,
     &                                cprop,iret)
              if(iret.eq.0) then
                ! they are equal
                equivto(lset,ichan)   = lsetpr
                equivcoef(lset,ichan) = 1D0
              elseif(iret.eq.1) then
                ! they are proportional
                equivto(lset,ichan)   = lsetpr
                equivcoef(lset,ichan) = cprop
              endif
            enddo
            call randomrestore
          endif
          flg_in_smartsig = .false.
          init(ichan)     = .false.
        endif
      
        ! End initialization phase; compute graphs
        call pdfcall(1,kn_x1,pdf1)
        call pdfcall(2,kn_x2,pdf2)
        do lset=1,flst_nosres
          if(equivto(lset,ichan).lt.0) then
            flst_cur_alr = -1
            !call realgr(flst_osres(1,lset),kn_cmpreal,r0(lset))
            call setosresreal(kn_cmpreal,flst_osres(1,lset),
     &                        osresID(ichan),r0(lset))
          else
            r0(lset)=r0(equivto(lset,ichan))*equivcoef(lset,ichan)
          endif
        enddo
        sig = 0D0
        do lset=1,flst_nosres
          r0(lset) = xjac*r0(lset)*
     &               pdf1(flst_osres(1,lset))*pdf2(flst_osres(2,lset))
          sig = sig+r0(lset)
        enddo
      end
c############### end subroutine sigreal_osres ##########################


c############### subroutine compare_vecs_osres #########################
c slightliy modified copy of compare_vecs_reg
c make sure to avoid 0-amplitudes correctly
      subroutine compare_vecs_osres(nmomset,lset,res,lsetpr,cprop,iret)
        implicit none

#include "nlegborn.h"
#include "pwhg_flst.h"

        double precision ep
        parameter (ep=1d-12)
        integer nmomset,lset,lsetpr,iret,j,k
        double precision res(nmomset,*),cprop,rat

        ! added this section
        !===============================================================
        if(lset.eq.1) then
          do j=1,nmomset
            if(res(j,lset).ne.0) then
              iret=-1
              return
            endif
          enddo
          iret   = 1 ! this is important
          lsetpr = 1 ! by default, take it prop to first flst with prop-fact 0
          cprop  = 0D0
          return
        endif
        !===============================================================

        do j=1,lset-1
          ! added this section
          !=============================================================
          ! if the res of the amp-routine is 0:
          ! make sure not to divide by 0 (is independent of momentum-config):
          if(res(1,lset).eq.0) then
            iret   = 1 ! this is important
            lsetpr = 1 ! by default, take it prop to first flst with prop-fact 0
            cprop  = 0D0
            return
          endif
          if(res(1,j).eq.0D0) goto 10 !no need to compare to a 0-result
          !=============================================================

          rat = res(1,lset)/res(1,j)
          do k=1,nmomset
            if(abs(1D0-res(k,lset)/res(k,j)/rat).gt.ep) goto 10
          enddo
          if(abs(1D0-rat).lt.ep) then
            iret  = 0
            cprop = 1D0
          else
            iret  = 1
            cprop = rat
          endif
          lsetpr = j
          return
 10       continue
        enddo
        iret = -1
      end
c############### end subroutine compare_vecs_osres #####################

c############### subroutine addupweightsosres ##########################
c the following routines are similar to the btilde-case...
      subroutine addupweightsosres(sigosres)
        implicit none

#include "nlegborn.h"
#include "pwhg_flst.h"
# include "pwhg_rad.h"

c keep this order
#include "osres.h"
#include "pwhg_flst_add.h"
#include "pwhg_rad_add.h"

        double precision sigosres

        double precision dtotosres(cnosres),dtotabsosres(cnosres)
        double precision dtotpososres(cnosres),dtotnegosres(cnosres)
        double precision totosres(cnosres),etotosres(cnosres)
        double precision totabsosres(cnosres),etotabsosres(cnosres)
        double precision totpososres(cnosres),etotpososres(cnosres)
        double precision totnegosres(cnosres),etotnegosres(cnosres)
        integer ncalls,i,j

        common/cadduptotalsosres/totosres,etotosres,totabsosres
        common/cadduptotalsosres/etotabsosres,totpososres
        common/cadduptotalsosres/etotpososres,totnegosres,etotnegosres
        common/cadduptotalsosres/ncalls

        ! keep track of the number of calls
        ncalls = ncalls + 1
        !print*,ncalls

        ! reset all contributions
        do j=1,nosres
          dtotosres(j)    = 0D0
          dtotabsosres(j) = 0D0
          dtotpososres(j) = 0D0
          dtotnegosres(j) = 0D0
        enddo

        ! sum over resonances
        do j=1,nosres
          ! sum over processes
          do i=1,flst_nosres
            dtotosres(j)    = dtotosres(j) + rad_osres_arr(i,j)
     &                                     * rad_osres_sign(i,j)
            dtotabsosres(j) = dtotabsosres(j) + rad_osres_arr(i,j)
            if(rad_osres_sign(i,j).eq.1) then
              dtotpososres(j) = dtotpososres(j) + rad_osres_arr(i,j)
            else
              dtotnegosres(j) = dtotnegosres(j) + rad_osres_arr(i,j)
            endif
          enddo

          totosres(j)     = totosres(j)     + dtotosres(j)
          etotosres(j)    = etotosres(j)    + dtotosres(j)**2
          totabsosres(j)  = totabsosres(j)  + dtotabsosres(j)
          etotabsosres(j) = etotabsosres(j) + dtotabsosres(j)**2
          totpososres(j)  = totpososres(j)  + dtotpososres(j)
          etotpososres(j) = etotpososres(j) + dtotpososres(j)**2
          totnegosres(j)  = totnegosres(j)  + dtotnegosres(j)
          etotnegosres(j) = etotnegosres(j) + dtotnegosres(j)**2

          sigosres = sigosres + dtotabsosres(i)
        enddo

#ifdef DEBUGQ
        print*,"sigosres",sigosres
#endif
      end

c############### end subroutine addupweightsosres ######################


c############### subroutine addupweightsosres ##########################
c set all new totals concerning regulars/remnants/osres to zero
      subroutine resettotalsosres
        implicit none

#include "osres.h"

        double precision totosres(cnosres),etotosres(cnosres)
        double precision totabsosres(cnosres),etotabsosres(cnosres)
        double precision totpososres(cnosres),etotpososres(cnosres)
        double precision totnegosres(cnosres),etotnegosres(cnosres)
        integer ncalls,j

        common/cadduptotalsosres/totosres,etotosres,totabsosres
        common/cadduptotalsosres/etotabsosres,totpososres
        common/cadduptotalsosres/etotpososres,totnegosres,etotnegosres
        common/cadduptotalsosres/ncalls

        ncalls = 0

        do j=1,nosres
          totosres(j)     = 0D0
          etotosres(j)    = 0D0
          totabsosres(j)  = 0D0
          etotabsosres(j) = 0D0
          totpososres(j)  = 0D0
          etotpososres(j) = 0D0
          totnegosres(j)  = 0D0
          etotnegosres(j) = 0D0
        enddo
      end

c############### end subroutine addupweightsosres ######################

c############### subroutine finaltotalsosres ###########################
c similar to corresponding routine in btilde
      subroutine finaltotalsosres
        implicit none

#include "nlegborn.h"
#include "pwhg_flst.h"
#include "pwhg_rad.h"

c keep this order
#include "osres.h"
#include "pwhg_flst_add.h"
#include "pwhg_rad_add.h"

        double precision calc_error
        external calc_error

        double precision totosres(cnosres),etotosres(cnosres)
        double precision totabsosres(cnosres),etotabsosres(cnosres)
        double precision totpososres(cnosres),etotpososres(cnosres)
        double precision totnegosres(cnosres),etotnegosres(cnosres)
        integer ncalls,j

        common/cadduptotalsosres/totosres,etotosres,totabsosres
        common/cadduptotalsosres/etotabsosres,totpososres
        common/cadduptotalsosres/etotpososres,totnegosres,etotnegosres
        common/cadduptotalsosres/ncalls

        ! if we never call the remnant-routines: avoid NaNs
        if(ncalls.eq.0) ncalls = 1

        do j=1,nosres
          rad_totosres(j)     = totosres(j)/ncalls
          rad_etotosres(j)    = calc_error(totosres(j),
     &                                     etotosres(j),ncalls)
          rad_totabsosres(j)  = totabsosres(j)/ncalls
          rad_etotabsosres(j) = calc_error(totabsosres(j),
     &                                     etotabsosres(j),ncalls)
          rad_totpososres(j)  = totpososres(j)/ncalls
          rad_etotpososres(j) = calc_error(totpososres(j),
     &                                     etotpososres(j),ncalls)
          rad_totnegosres(j)  = totnegosres(j)/ncalls
          rad_etotnegosres(j) = calc_error(totnegosres(j),
     &                                     etotnegosres(j),ncalls)
        enddo
      end
c############### end subroutine finaltotalsosres #######################

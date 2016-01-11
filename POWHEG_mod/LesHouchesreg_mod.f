c MK: copied and modified version of LesHouchesreg.f, revision 3154
c modified on the basis of disquark/LesHouchesreg.f
c changes marked with "! CH, MK:" or "! MK:"
c bigger changes over a whole section are marked with !===...

      subroutine gen_leshouches_reg
      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      include 'LesHouches.h'
      integer ireg
      nup=nlegreal
      scalup=sqrt(rad_pt2max)
      call momenta_lh(kn_preal,nlegreal) !CH: need to set momenta explicitly!
      do ireg=1,nup
c Remember: gluons are marked 0 here!
         idup(ireg)=flst_regular(ireg,rad_realreg)
         if (idup(ireg).eq.0) idup(ireg)=21
         if(ireg.le.2) then
            istup(ireg)=-1
            mothup(1,ireg)=0
            mothup(2,ireg)=0            
         else
            istup(ireg)=1
            mothup(1,ireg)=1
            mothup(2,ireg)=2            
         endif
         spinup(ireg)=9
         vtimup(ireg)=0
      enddo
c add resonances, perform decays, put particles on shell, etc.(or nothing!)
      call finalize_lh 
c no remnants for now!
c     Don't forget to set scale for scalup equal to the pt of the 
c     radiation (whatever it is now!)
c     set color connections for all particles
      write(*,*) 'gen_leshouches_reg: dummy interface to regular'//
     1           ' contributions'
      write(*,*) ' Replace with your own process-dependent one, to be'
      write(*,*) ' put in the process-specific directory (e.g. /W,'//'
     #     /Z, /VBF_H...)'
      write(*,*) ' The Makefile will automatically compile the version'
      write(*,*) ' in the process-specific directory'
      call exit(1)
      end

c=======================================================================
c CH, MK: added the following routine specifically for the osres-contributions
c relevant for the treatment of the real CFs
      subroutine gen_leshouches_osres
      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      include 'LesHouches.h'
c MK: added
#include "osres.h"
#include "pwhg_flst_add.h"
#include "pwhg_rad_add.h"
      integer ireg,moth,j,iub,i
      logical flag_res45
      nup=nlegreal
      call momenta_lh(kn_preal,nlegreal) !always have radiation here
      ! id of the event
      idprup=lprup(1)

      do ireg=1,nup
         ! Remember: gluons are marked 0 here!
         idup(ireg)=flst_osres(ireg,rad_realosres)
         if (idup(ireg).eq.0) idup(ireg)=21
         if(ireg.le.2) then
            istup(ireg)=-1
            mothup(1,ireg)=0
            mothup(2,ireg)=0
         else
            istup(ireg)=1
            mothup(1,ireg)=1
            mothup(2,ireg)=2
         endif
         spinup(ireg)=9
         vtimup(ireg)=0
      enddo
      call lh_resonances
      call finalize_lh ! flg_btilde is in this case false: we create a regular...

      ! radiated parton always from born
      mothup(1,nup)=1
      mothup(2,nup)=2
        
      ! Don't forget to set scale for scalup equal to the pt of the
      ! radiation (whatever it is now!)
      scalup=sqrt(rad_pt2max) !rad_pt2max was set using the pt2max_osres-routine!

      end
c=======================================================================

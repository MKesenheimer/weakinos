#!/bin/bash

cat <<EOF
 Insert in Makefile:
 MADLIB=-L MadGraphStuff/MadGraph_POWHEG/my_proc/lib/ -lmadgraph -lmodel -ldhelas3
 USER=testsigmas.o ....
 add \$(MADLIB) to the make target
 edit Born.f, insert:
 call checkbornwithmg((p,flav,born,bornjk,bmunu)
 before return
 edit real.f; insert
 call checkrealwithmg(p,flav,real)
 before return.
 copy the coupl.inc file to the compilation directory
 copy the Cards/param_card.dat into a run directory,
 under a subdirectory named Cards
EOF


cat Born.f | sed 's/subroutine setborn/subroutine msetborn/ ; s/subroutine borncolour_lh/subroutine mborncolour_lh/ ; s/subroutine finalize_lh/subroutine mfinalize_lh/' > testsigmas.f

cat real.f | sed 's/subroutine setreal/subroutine msetreal/' >> testsigmas.f

cat init_couplings.f | sed 's/subroutine init_couplings/subroutine minit_couplings/' >> testsigmas.f


cat <<EOF >> testsigmas.f

      subroutine checkbornwithmg(p,flav,born,bornjk,bmunu)
      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      integer nlegs
      parameter (nlegs=nlegborn)
      integer flav(nlegs)
      real * 8 p(0:3,nlegs)
      real * 8 born,bornjk(nlegs,nlegs),bmunu(0:3,0:3,nlegs)
      real * 8 mborn,mbornjk(nlegs,nlegs),mbmunu(0:3,0:3,nlegs)
      integer j,k,mu,nu

      call mminit_couplings

      call msetborn(p,flav,mborn,mbornjk,mbmunu)

      if(mborn.ne.0) then
         if(abs(born/mborn-1).gt.1d-7) then
            write(*,*) flav
            write(*,*) ' born/mborn=',born/mborn
            return
         endif
      else
         if(born.ne.0) then
            write(*,*) flav
            write(*,*) ' born/mborn=Inf'
            return
         endif
      endif
      do j=1,nlegborn
      do mu=0,3
      do nu=0,3
         if(mbmunu(mu,nu,j).ne.0) then
            if(abs(bmunu(mu,nu,j)/mbmunu(mu,nu,j)-1).gt.1d-7) then
               write(*,*) flav
               write(*,*) ' bmunu(',mu,',',nu,',',j,')/mbmunu=',
     1              bmunu(mu,nu,j)/mbmunu(mu,nu,j)
            endif
         else
            if(bmunu(mu,nu,j).ne.0) then
               write(*,*) flav
               write(*,*) ' bmunu:(',mu,',',nu,',',j,')/mbmunu=Inf'
            endif
         endif
      enddo
      enddo
      enddo
      do j=1,nlegborn
         do k=1,nlegborn
            if(j.ne.k.or.kn_masses(j).ne.0) then
               if(mbornjk(j,k).ne.0) then
                  if(abs(bornjk(j,k)/mbornjk(j,k)-1).gt.1d-7) then
                     write(*,*) flav
                     write(*,*) ' bjk(',j,',',k,')/mbjk=',
     1                    bornjk(j,k)/mbornjk(j,k)
                  endif
               else
                  if(bornjk(j,k).ne.0) then
                     write(*,*) flav
                     write(*,*) ' bjk(',j,',',k,')/mbjk=Inf'
                  endif
               endif
            endif
         enddo
      enddo
      end

      subroutine checkrealwithmg(p,flav,real)
      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer nlegs
      parameter (nlegs=nlegreal)
      integer flav(nlegs)
      real * 8 p(0:3,nlegs)
      real * 8 real,mreal

      call mminit_couplings

      call msetreal(p,flav,mreal)

      if(mreal.ne.0) then
         if(abs(real/mreal-1).gt.1d-7) then
            write(*,*) flav
            write(*,*) ' real/mreal=',real/mreal
         endif
      else
         if(real.ne.0) then
            write(*,*) flav
            write(*,*) ' real/mreal=Inf'
         endif
      endif
      end

      subroutine mminit_couplings
      logical ini
      data ini/.true./
      save ini

      if(ini) then
         call minit_couplings
         ini=.false.
      endif

      end



EOF

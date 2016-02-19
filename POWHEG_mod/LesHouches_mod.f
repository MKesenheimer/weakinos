c MK: copied and modified version of LesHouches_mod.f, revision 3154
c modified on the basis of disquark/LesHouches.f
c changes marked with "! CH, MK:" or "! MK:"
c bigger changes over a whole section are marked with !===...

      subroutine gen_leshouches
      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      include 'LesHouches.h'
      integer alr,em,rad,flem,flrad,coluborn(2)
      if(kn_csi.eq.0d0) then
c it is a Born event or a LO event
c first assign flavours and colour to the underlying Born process
         call momenta_lh(kn_pborn,nlegborn)
         call born_lh
         if (flg_LOevents) then
c     set default scalup to the sqrt of the center-of-mass energy of the event.
c     For processes that need different values, the scalup value can be
c     overridden in the finalize_lh user subroutine.
            scalup=2*kn_cmpborn(0,1)
         else
            scalup=sqrt(rad_ptsqmin)
         endif
      else
c it is an event with radiation
         call momenta_lh(kn_preal,nlegreal)
         call born_lh
         scalup=sqrt(rad_pt2max)
         nup=nlegreal
         istup(nup)=1
         spinup(nup)=9
         vtimup(nup)=0
         mothup(1,nup)=1
         mothup(2,nup)=2
         alr=rad_realalr
         em=flst_emitter(alr)
         rad=nlegreal
         if(em.eq.0) then
            if(kn_y.gt.0) then
               em=1
            else
               em=2
            endif
         endif
         flem=flst_alr(em,alr)
         flrad=flst_alr(rad,alr)
         coluborn(1)=icolup(1,em)
         coluborn(2)=icolup(2,em)
         if(em.le.2) then
c setcolour_rad works with all incoming (or all outgoing) flavours and colours;
c for ISR make everything outgoing:
c thus in input change the sign of the emitter flavour, to make it outgoing
c and conjugate its colour in output, to make it incoming.
            call setcolour_rad(coluborn,-flem,flrad,
     #                     icolup(1,em),icolup(1,rad))
            call colour_conj(icolup(1,em))
         else
c For ISR make everything incoming;
c Change the sign of the emitter and radiated flavours in input,
c to make it incoming;
c conjugate their colours in the output, to make them outgoing.
            call setcolour_rad(coluborn,-flem,-flrad,
     #                     icolup(1,em),icolup(1,rad))
            call colour_conj(icolup(1,rad))
            call colour_conj(icolup(1,em))
         endif
#ifdef DEBUGQ
         print*,"[DEBUG] in LesHouches_mod.f:73"
         print*,"flst_alr(:,alr)",flst_alr(:,alr)
         print*,"icolup(:,",1,")",icolup(:,1)
         print*,"icolup(:,",2,")",icolup(:,2)
         print*,"icolup(:,",5,")",icolup(:,5)
         print*
#endif
         ! MK: added gluon color checks
         if((flem.eq.0) .and. 
     &      ((icolup(1,em).eq.0) .or. (icolup(2,em).eq.0))) then
           print*,"error in LesHouches_mod.f:83"
           print*,"emitter is gluon, but has too less icolup entries:"
           print*,"icolup(:,",em,")",icolup(:,em)
           stop
         endif
         if((flrad.eq.0) .and. 
     &      ((icolup(1,rad).eq.0) .or. (icolup(2,rad).eq.0))) then
           print*,"error in LesHouches_mod.f:90"
           print*,"radiator is gluon, but has too less icolup entries:"
           print*,"icolup(:,",rad,")",icolup(:,rad)
           stop
         endif  
         if(flem.eq.0) flem=21
         if(flrad.eq.0) flrad=21
         idup(em)=flem
         idup(rad)=flrad
      endif
c if there are resonances, setup up appropriate LH pointers
      call lh_resonances
c add resonances, perform decays, put particles on shell, etc.(or nothing!)
      call finalize_lh
      end


c     i1<i2
      subroutine add_resonance(idpdg,i1,i2)
      implicit none
      include 'LesHouches.h'
      real * 8 iduptmp(maxnup),istuptmp(maxnup),mothuptmp(2,maxnup),
     #     icoluptmp(2,maxnup),puptmp(5,maxnup),vtimuptmp(maxnup),
     #     spinuptmp(maxnup)
      integer idpdg,i1,i2,i,j
      if (i1.ge.i2) then
         write(*,*) 'wrong sequence in add_resonance'
         stop
      endif
      do i=i1,nup
         iduptmp(i) = idup(i)
         istuptmp(i) = istup(i)
         mothuptmp(1,i) = mothup(1,i)
         mothuptmp(2,i) = mothup(2,i)
         icoluptmp(1,i) = icolup(1,i)
         icoluptmp(2,i) = icolup(2,i)
         do j=1,5
            puptmp(j,i)=pup(j,i)
         enddo
         vtimuptmp(i)=vtimup(i)
         spinuptmp(i)=spinup(i)
      enddo
      idup(i1)=idpdg
      istup(i1)=+2
      mothup(1,i1)=1
      mothup(2,i1)=2
      icolup(1,i1)=0
      icolup(2,i1)=0
      do j=1,4
         pup(j,i1)=puptmp(j,i1)+puptmp(j,i2)
      enddo
      pup(5,i1)=sqrt(pup(4,i1)**2-
     #     (pup(1,i1)**2+pup(2,i1)**2+pup(3,i1)**2))
      vtimup(i1)=0
      spinup(i1)=9
c     change mothers of decaying particles
      mothuptmp(1,i1)=i1
      mothuptmp(2,i1)=i1
      mothuptmp(1,i2)=i1
      mothuptmp(2,i2)=i1
      nup=nup+1      
      do i=i1+1, nup
         idup(i) = iduptmp(i-1)
         istup(i) = istuptmp(i-1)
         mothup(1,i) = mothuptmp(1,i-1)
         mothup(2,i) = mothuptmp(2,i-1)
         icolup(1,i) = icoluptmp(1,i-1)
         icolup(2,i) = icoluptmp(2,i-1)
         do j=1,5
            pup(j,i)=puptmp(j,i-1)
         enddo
         vtimup(i)=vtimuptmp(i-1)
         spinup(i)=spinuptmp(i-1)         
      enddo
      end
      


      subroutine colour_conj(icol)
      implicit none
      integer icol(2)
      integer itmp
      itmp=icol(2)
      icol(2)=icol(1)
      icol(1)=itmp
      end

      subroutine momenta_lh(p,n)
      implicit none
      integer n
      real * 8 p(0:3,n)
      integer k,mu
      include 'LesHouches.h'
      do k=1,n
         do mu=1,3
            pup(mu,k)=p(mu,k)
         enddo
         pup(4,k)=p(0,k)
         pup(5,k)=sqrt(abs(p(0,k)**2-p(1,k)**2
     #       -p(2,k)**2-p(3,k)**2))
      enddo
      end

      subroutine born_lh
      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      include 'LesHouches.h'
c      include 'PhysPars.h'
      include 'pwhg_st.h'
      integer k,fl
      nup=nlegborn
c id of the event
      idprup=lprup(1)
ccccc CAVEAT: aqedup must not be set here
ccccc since otherwise a dependence on the
ccccc process-dependent PhysPars.h must be
ccccc introduced
      aqedup=-1                 ! ph_alphaem
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      aqcdup=st_alpha

      do k=1,nlegborn
         fl=flst_born(k,rad_ubornidx)
c gluons are numbered 21 in pdg
         if(fl.eq.0) fl=21
         if(fl.eq.1000000) fl=1000021 ! CH: added for gluinos
         idup(k)=fl
         if(k.gt.2) then
            istup(k)=1
            mothup(1,k)=1
            mothup(2,k)=2
         else
            istup(k)=-1
            mothup(1,k)=0
            mothup(2,k)=0
         endif
         spinup(k)=9
         vtimup(k)=0
      enddo
      call borncolour_lh
      end

      subroutine lh_resonances
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'LesHouches.h'
      integer moth,j,iub,em
      iub = rad_ubornidx
      do j=1,nlegborn
         moth=flst_bornres(j,iub)
         if(moth.ne.0) then
            mothup(1,j)=moth
            mothup(2,j)=moth
            istup(moth)=2
         endif
      enddo
      if(nup.eq.nlegreal) then
         em = flst_emitter(rad_realalr)
         if(em.ne.0) then
            moth=flst_bornres(flst_emitter(rad_realalr),iub)
            if(moth.ne.0) then
               mothup(1,nlegreal)=moth
               mothup(2,nlegreal)=moth
            endif
         endif
      endif
      end

      subroutine getnewcolor(newcolor)
      implicit none
      integer newcolor
      include 'LesHouches.h'
      integer j,k
      newcolor=511
 1    continue
      do j=1,nup
         do k=1,2
            if(icolup(k,j).eq.newcolor) then
               newcolor=newcolor+1
               goto 1
            endif
         enddo
      enddo
      end

c Given a vertex with 3 partons, the first parton with incoming colour
c colin, the second and third partons with incoming flavours fl2, fl3,
c it assigns the incoming colours of the second and third partons
c colout2, colout3. The colour assignment is unambiguous in all cases, but in
c the one where both fl2 and fl3 are gluons, where a 50% random choice is made
c over the two possible colour connections.
      subroutine setcolour_rad(colin,fl2,fl3,colout2,colout3)
      implicit none
      integer colin(2),fl2,fl3,colout2(2),colout3(2)
      integer newcolor
      real * 8 random
      logical is_coloured
      external random,is_coloured
      call getnewcolor(newcolor)
c First handle photon case
      if(fl2.eq.22) then
         colout2=0
         colout3(1)=colin(2)
         colout3(2)=colin(1)
      elseif(fl3.eq.22) then
         colout3=0
         colout2(1)=colin(2)
         colout2(2)=colin(1)
      elseif(colin(1).eq.0.and.colin(2).eq.0) then
         if(is_coloured(fl2)) then
            if(fl2.gt.0) then
               colout2(1)=newcolor
               colout2(2)=0
               colout3(2)=newcolor
               colout3(1)=0
            elseif(fl3.gt.0) then
               colout2(1)=0
               colout2(2)=newcolor
               colout3(2)=0
               colout3(1)=newcolor
            else
               goto 998
            endif
         else
            colout2=0
            colout3=0
         endif
      else
c now handle all coloured particles
         if(colin(1).ne.0.and.colin(2).ne.0) then
            if(fl2.eq.0.and.fl3.eq.0) then
               if(random().gt.0.5d0) then
                  colout2(1)=colin(2)
                  colout2(2)=newcolor
                  colout3(1)=newcolor
                  colout3(2)=colin(1)
               else
                  colout3(1)=colin(2)
                  colout3(2)=newcolor
                  colout2(1)=newcolor
                  colout2(2)=colin(1)
               endif
            elseif(fl2.gt.0.and.fl3.lt.0) then
               colout2(2)=0
               colout2(1)=colin(2)
               colout3(1)=0
               colout3(2)=colin(1)
            elseif(fl2.lt.0.and.fl3.gt.0) then
               colout2(1)=0
               colout2(2)=colin(1)
               colout3(2)=0
               colout3(1)=colin(2)
            else
               goto 998
            endif
         elseif(colin(2).eq.0) then
            if(fl2.eq.0) then
               colout3(1)=0
               colout3(2)=newcolor
               colout2(1)=newcolor
               colout2(2)=colin(1)
            elseif(fl3.eq.0) then
               colout2(1)=0
               colout2(2)=newcolor
               colout3(1)=newcolor
               colout3(2)=colin(1)
            else
               goto 998
            endif
         elseif(colin(1).eq.0) then
            if(fl2.eq.0) then
               colout3(1)=newcolor
               colout3(2)=0
               colout2(2)=newcolor
               colout2(1)=colin(2)
            elseif(fl3.eq.0) then
               colout2(1)=newcolor
               colout2(2)=0
               colout3(1)=colin(2)
               colout3(2)=newcolor
            else
               goto 998
            endif
         else
            goto 998
         endif
      endif
      return
 998  continue
      write(*,*)
     #     ' Error in setcolour_rad: inconsistent colour flavour input'
      write(*,*) 'colin[1:2] ',colin(1),' ',colin(2)
      write(*,*) 'fl2, fl3 ',fl2,' ',fl3
      write(*,*) 'newcolor ',newcolor
      call pwhg_exit(1)
      end


      subroutine displeshouches
      implicit none
      include 'LesHouches.h'
      integer pdgfac
      parameter(pdgfac=1000000)
      integer mark(200),icol(200),j,k,ilist,icur,imarked
      write(*,*) 'incoming beams:', idbmup(1), idbmup(2)
      write(*,*) 'number of partons in subprocess',nup
      write(*,'(a,10(i3,1x))') 'ist:',(istup(j),j=1,nup)
      write(*,'(a,10(i8,1x))') 'ids:',(idup(j),j=1,nup) ! CH: i3->i8, to include SUSY-PDGs as well
      do j=1,nup
         write(*,'(5(d10.4,1x))') (pup(k,j),k=1,5)
      enddo
c check that quarks have zero anticolour, and antiquarks
c have zero colour
      do j=1,nup
         if(idup(j).ge.1.and.idup(j).le.6.and.icolup(2,j).ne.0) then
            write(*,*) ' error: anticolor in quark'
         elseif(idup(j).ge.-6.and.idup(j).le.-1.and.icolup(1,j).ne.0)
     #           then
            write(*,*) ' error: color in antiquark'
         endif
           if(abs(idup(j)).ge.pdgfac) then
           ! CH: change this check to include squarks, antisquarks as well
           if(mod(idup(j),pdgfac).ge.1.and.mod(idup(j),pdgfac).le.6.and.icolup(2,j).ne.0) then
             write(*,*) ' error: anticolor in squark'
           elseif(mod(idup(j),pdgfac).ge.-6.and.mod(idup(j),pdgfac).le.-1.and.icolup(1,j).ne.0)
     #           then
             write(*,*) ' error: color in antisquark' ! CH: changed to anti 
           endif
             if(idup(j).eq.pdgfac+21) then
               if((icolup(1,j)).eq.0.or.icolup(2,j).eq.0) then
                write(*,*) ' error: no color or anticolor in gluino: ',icolup(1:2,j)
               endif
             endif
           endif
           ! CH: check if g/go have color and anticolor:
           if(idup(j).eq.21) then
             if((icolup(1,j)).eq.0.or.icolup(2,j).eq.0) then
             write(*,*) ' error: no color or anticolor in gluon: ',icolup(1:2,j)
             endif
           endif
      enddo
c Print color linked lists
c first conjugate incoming colors
      call colour_conj(icolup(1,1))
      call colour_conj(icolup(1,2))
      do j=1,nup
         mark(j)=1
      enddo
      imarked=0
      do j=1,nup
         if(istup(j).eq.2.or.
     #        (icolup(1,j).eq.0.and.icolup(2,j).eq.0)) then
            mark(j)=0
            imarked=imarked+1
         endif
         if(istup(j).ne.2.and.(((idup(j).gt.0.and.idup(j).lt.6)).or.
     &         (idup(j).gt.pdgfac.and.mod(idup(j),pdgfac).lt.6))) then !CH: added this line
            icur=j
         endif
      enddo
      mark(icur)=0
      ilist=1
      imarked=imarked+1
      icol(ilist)=icur
 12   continue
      do j=1,nup
         if(istup(j).ne.2.and.
     #        (icolup(1,j).ne.0.or.icolup(2,j).ne.0)) then
         if(mark(j).ne.0.and.icolup(1,icur).eq.icolup(2,j)) then
            ilist=ilist+1
            imarked=imarked+1
            icol(ilist)=j
            mark(j)=0
            icur=j
            if(imarked.eq.nup) goto 22
            if(icolup(1,icur).eq.0) then
               do k=1,nup
                  if(mark(k).ne.0.and.icolup(2,k).eq.0) then
                     ilist=ilist+2
                     icol(ilist)=k
                     icol(ilist-1)=0
                     imarked=imarked+1
                     icur=k
                     mark(k)=0
                     goto 12
                  endif
               enddo
               write(*,*) ' inconsistent colors!'
            endif
            goto 12
         endif
         endif
      enddo
      write(*,*) ' inconsistent colors!'
 22   continue
      call colour_conj(icolup(1,1))
      call colour_conj(icolup(1,2))
      write(*,*) (icol(j),j=1,ilist)
      end



      subroutine pdfreweightinfo(id1,id2,x1,x2,xmufact,xf1,xf2)
      implicit none
      integer id1,id2
      real * 8 x1,x2,xf1,xf2,xmufact
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_rad.h'
      include 'pwhg_pdf.h'
      include 'pwhg_st.h'
      include 'LesHouches.h'
c MK: added
#include "osres.h"
#include "pwhg_flst_add.h"
#include "Flags.h"
#include "pwhg_rad_add.h"
      real * 8 pdf(-pdf_nparton:pdf_nparton)
      if(rad_type.eq.1) then
c Btilde event: pass x1 and x2, id1, id2 etc. of the current underlying Born
         id1=flst_born(1,rad_ubornidx)
         id2=flst_born(2,rad_ubornidx)
         x1=kn_xb1
         x2=kn_xb2
         call setscalesbtilde
      elseif(rad_type.eq.2.or.rad_type.eq.3.or.
       ! MK: added 4<=rad_type<=3+nosres
     & ((rad_type.ge.4).and.(rad_type.le.(nosres+3)))) then
         id1=idup(1)
         id2=idup(2)
         if(id1.eq.21) id1=0
         if(id2.eq.21) id2=0         
         x1=kn_x1
         x2=kn_x2
         call setscalesbtilde
      endif
      xmufact=sqrt(st_mufact2)
      call pdfcall(1,x1,pdf)
      xf1=x1*pdf(id1)
      call pdfcall(2,x2,pdf)
      xf2=x2*pdf(id2)
      if(id1.eq.0) id1=21
      if(id2.eq.0) id2=21
      end

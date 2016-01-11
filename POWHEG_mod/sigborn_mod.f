c MK: copied and modified version of gen_real_phsp.f, revision 3154
c modified on the basis of disquark/gen_real_phsp.f
c changes marked with "! CH, MK:" or "! MK:"
c bigger changes over a whole section are marked with !===...

      subroutine btildeborn(res)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_br.h'
      include 'pwhg_flg.h'
      include 'pwhg_pdf.h'
      real * 8 pdf1(-pdf_nparton:pdf_nparton),
     1         pdf2(-pdf_nparton:pdf_nparton)
      real * 8 res(flst_nborn),tot,rescfac
      integer j
      if(.not.flg_minlo) then
         rescfac = 1
         call pdfcall(1,kn_xb1,pdf1)
         call pdfcall(2,kn_xb2,pdf2)
      endif
      tot=0
      do j=1,flst_nborn
         if(flg_minlo) then
            call setlocalscales(j,1,rescfac)
            call pdfcall(1,kn_xb1,pdf1)
            call pdfcall(2,kn_xb2,pdf2)
         endif
         res(j)=br_born(j) *
     #  pdf1(flst_born(1,j))*pdf2(flst_born(2,j))*kn_jacborn*rescfac
         tot=tot+res(j)
c     if one wants to do only the integration over the whole phase space, then
c     uncomment the following (with or without the flux factor 1/(2*kn_sborn)
c        
c         res(j)=kn_jacborn/(2*kn_sborn)/flst_nborn
      enddo
c If born rescaling of some sort is needed, the user should provide
c its own rescaling procedure in a file named pwhg_born_rescaling_hook.h.
c The default one is in the include directory, and is empty
      include 'pwhg_born_rescaling_hook.h'
      end

      subroutine sigborn_rad(born)
      implicit none
      real * 8 born
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_br.h'
      include 'pwhg_rad.h'
      include 'pwhg_pdf.h'
      real * 8 pdf1(-pdf_nparton:pdf_nparton),
     1         pdf2(-pdf_nparton:pdf_nparton),
     2         bornjk(nlegborn,nlegborn),bmunu(0:3,0:3,nlegborn)
      call pdfcall(1,kn_xb1,pdf1)
      call pdfcall(2,kn_xb2,pdf2)
      flst_cur_iborn = rad_ubornidx
      call setborn0(kn_cmpborn,flst_born(1,rad_ubornidx),born,
     #        bornjk,bmunu)
c Store in the br_born arrays; they are used to compute
c collinear and soft approximations, and affect the subtraction
c of the remnant component.
      br_born(rad_ubornidx)=born
      br_bornjk(:,:,rad_ubornidx)=bornjk
      br_bmunu(:,:,:,rad_ubornidx)=bmunu

      born=born *
     #  pdf1(flst_born(1,rad_ubornidx))*pdf2(flst_born(2,rad_ubornidx))
      end

     

      subroutine allborn
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_br.h'
      include 'pwhg_flg.h'
      integer equivto(maxprocborn)
      common/cequivtoborn/equivto
      real * 8 equivcoef(maxprocborn)
      common/cequivcoefborn/equivcoef
      integer nmomset
      parameter (nmomset=10)
      real * 8 pborn(0:3,nlegborn,nmomset),cprop
      real * 8 born(nmomset,maxprocborn)
      real * 8 bornjk(nlegborn,nlegborn,nmomset,maxprocborn)
      real * 8 bmunu(0:3,0:3,nlegborn,nmomset,maxprocborn)
      integer iborn,ibornpr,mu,nu,k,j,iret
      logical ini
      data ini/.true./
      save ini,/cequivtoborn/,/cequivcoefborn/
      if(ini) then
         do iborn=1,flst_nborn
            equivto(iborn)=-1
         enddo
         if(flg_smartsig) then
            flg_in_smartsig = .true.
            call randomsave
            call fillmomenta(nlegborn,nmomset,kn_masses,pborn)
            do iborn=1,flst_nborn
               do j=1,nmomset
                  flst_cur_iborn = iborn
                  call setborn0(pborn(0,1,j),flst_born(1,iborn),
     1                 born(j,iborn),bornjk(1,1,j,iborn),
     2                 bmunu(0,0,1,j,iborn))
               enddo
               call compare_vecsb(nmomset,iborn,born,bornjk,bmunu,
     1              ibornpr,cprop,iret)
               if(iret.eq.0) then
                  equivto(iborn)=ibornpr
                  equivcoef(iborn)=1
               elseif(iret.eq.1) then
                  equivto(iborn)=ibornpr
                  equivcoef(iborn)=cprop
               endif
            enddo
            call randomrestore
         endif
         flg_in_smartsig = .false.
         ini=.false.
      endif
      do iborn=1,flst_nborn
         if(equivto(iborn).lt.0) then
            flst_cur_iborn = iborn
            call setborn0(kn_cmpborn,flst_born(1,iborn),br_born(iborn),
     #        br_bornjk(1,1,iborn),br_bmunu(0,0,1,iborn))
         else
            br_born(iborn)=br_born(equivto(iborn))*equivcoef(iborn)
            do j=1,nlegborn
               do k=1,nlegborn
                  br_bornjk(j,k,iborn)=br_bornjk(j,k,equivto(iborn))
     #*equivcoef(iborn)
               enddo
            enddo
            do mu=0,3
               do nu=0,3
                  do j=1,nlegborn
                     br_bmunu(mu,nu,j,iborn)=
     #                   br_bmunu(mu,nu,j,equivto(iborn))
     #                   *equivcoef(iborn)
                  enddo
               enddo
            enddo
         endif
      enddo
      end

      subroutine compare_vecsb(nmomset,iborn,born,bornjk,bmunu,
     1     ibornpr,cprop,iret)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      real * 8 ep
      parameter (ep=1d-12)
      integer nmomset,iborn,ibornpr,iret,jborn,k,jleg,kleg,mu,nu
      real * 8 born(nmomset,iborn),cprop,rat,resi,resj
      real * 8 bornjk(nlegborn,nlegborn,nmomset,maxprocborn)
      real * 8 bmunu(0:3,0:3,nlegborn,nmomset,maxprocborn)
      do jborn=1,iborn-1
         rat=born(1,iborn)/born(1,jborn)
         do k=1,nmomset
            resi=born(k,iborn)
            resj=born(k,jborn)
            if(abs(1-resi/resj/rat).gt.ep) goto 10
         enddo
c totals are identical; see if also colour correlated are
         do jleg=1,nlegborn
            do kleg=1,nlegborn
               do k=1,nmomset
                  resi=bornjk(jleg,kleg,k,iborn)
                  resj=bornjk(jleg,kleg,k,jborn)
c some of these are zero, be careful
                  if(resi.ne.resj*rat) then
                     if(abs(1-resi/resj/rat).gt.ep) goto 10
                  endif
               enddo
            enddo
         enddo
c totals are identical; see if also spin correlated are
         do jleg=1,nlegborn
            do mu=0,3
               do nu=0,3
                  do k=1,nmomset
                     resi=bmunu(mu,nu,jleg,k,iborn)
                     resj=bmunu(mu,nu,jleg,k,jborn)
c some of these are zero, be careful
                     if(resi.ne.resj*rat) then
                        if(abs(1-resi/resj/rat).gt.ep) goto 10
                     endif
                  enddo
               enddo
            enddo
         enddo
         if(abs(1-rat).lt.ep) then
            iret=0
            cprop=1
         else
            iret=1
            cprop=rat
         endif
         ibornpr=jborn
         return
 10      continue
      enddo
      iret=-1
      end

      subroutine setborn0(p,bflav,born,bornjk,bmunu)
c provide the flux factor to the user Born routine
      implicit none
      include 'pwhg_math.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      integer nlegs
      parameter (nlegs=nlegborn)
      real * 8 p(0:3,nlegs)
      integer bflav(nlegs)
      real * 8 born,bornjk(nlegs,nlegs),bmunu(0:3,0:3,nlegs)
      integer j,k,mu,nu
      logical pwhg_isfinite
      external pwhg_isfinite
      call setborn(p,bflav,born,bornjk,bmunu)
c     check if born, bornjk and bmunu are finite
      if (.not.pwhg_isfinite(born)) born=0d0
      do j=1,nlegs
         do mu=0,3
            do nu=0,3
               if (.not.pwhg_isfinite(bmunu(mu,nu,j))) 
     $               bmunu(mu,nu,j)=0d0
            enddo
         enddo
      enddo
      do j=1,nlegs
         do k=1,nlegs
            if (.not.pwhg_isfinite(bornjk(j,k))) bornjk(j,k)=0d0
         enddo
      enddo
c
      born=born/(2*kn_sborn)
      do j=1,nlegs
         do mu=0,3
            do nu=0,3
               bmunu(mu,nu,j)=bmunu(mu,nu,j)/(2*kn_sborn)
            enddo
         enddo
      enddo
      do j=1,nlegs
         do k=1,nlegs
            bornjk(j,k)=bornjk(j,k)/(2*kn_sborn)
         enddo
      enddo
      end


      subroutine printbornequiv
c When invoked after the first call to allborn,
c it prints the set of equivalent Born configurations
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer equivto(maxprocborn)
      common/cequivtoborn/equivto
      real * 8 equivcoef(maxprocborn)
      common/cequivcoefborn/equivcoef
      integer j,k,iun,count
      save count
      data count/0/
      call newunit(iun)
      open(unit=iun,file='bornequiv',status='unknown')
      write(*,*) 'Writing bornequiv file...'
      do j=1,flst_nborn
         if(equivto(j).eq.-1) then
            write(iun,'(a)')
     1           'Beginning sequence of equivalent amplitudes'
            write(iun,100) 1d0,j, flst_born(:,j)
            do k=1,flst_nborn
               if(equivto(k).eq.j) then
                  write(iun,100) equivcoef(k),k,flst_born(:,k)
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
      ! 1000022 not correctly written out in bornequiv file, i3 -> i7):
 100  format(d11.4,5x,i8,5x,100(i8,1x))
      end

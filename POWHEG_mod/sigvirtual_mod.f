c MK: copied and modified version of sigvirtual.f, revision 3154
c modified on the basis of disquark/sigvirtual.f
c changes marked with "! CH, MK:" or "! MK:"
c bigger changes over a whole section are marked with !===...

      subroutine sigvirtual(virt_arr)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_kn.h'
      include 'pwhg_br.h'
      include 'pwhg_flg.h'
      real * 8 virt_arr(maxprocborn)
      integer equivto(maxprocborn)
      common/cequivtovirt/equivto
      real * 8 equivcoef(maxprocborn)
      common/cequivcoefvirt/equivcoef
      integer nmomset
      parameter (nmomset=10)
      real * 8 pborn(0:3,nlegborn,nmomset),cprop
      real * 8 virtual(nmomset,maxprocborn)
      integer iborn,ibornpr,mu,nu,k,j,iret
      logical ini
      data ini/.true./
      save ini,/cequivtovirt/,/cequivcoefvirt/
      logical pwhg_isfinite
      external pwhg_isfinite
      if(ini) then
         do iborn=1,flst_nborn
            equivto(iborn)=-1
         enddo
         if(flg_smartsig.and..not.flg_novirtual) then
            flg_in_smartsig = .true.
            call randomsave
            call fillmomenta(nlegborn,nmomset,kn_masses,pborn)
            do iborn=1,flst_nborn
               do j=1,nmomset
                  flst_cur_iborn = iborn 
                  call setvirtual(pborn(0,1,j),flst_born(1,iborn),
     #                 virtual(j,iborn))
c     check if virtual(j,iborn) is finite
                  if (.not.pwhg_isfinite(virtual(j,iborn))) 
     #                 virtual(j,iborn)=0d0
               enddo
               call compare_vecsv(nmomset,iborn,virtual,ibornpr,
     #              cprop,iret)
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
         call printvirtequiv
      endif
      do iborn=1,flst_nborn
         flst_cur_iborn = iborn ! MK: added
         if(equivto(iborn).lt.0) then
            flst_cur_iborn = iborn
            call setvirtual(kn_cmpborn,flst_born(1,iborn),
     #           virt_arr(iborn))
c     check if virt_arr(iborn) is finite
                  if (.not.pwhg_isfinite(virt_arr(iborn))) 
     #                 virt_arr(iborn)=0d0
            virt_arr(iborn)=virt_arr(iborn)/(2*kn_sborn)
         else
            virt_arr(iborn)=virt_arr(equivto(iborn))*equivcoef(iborn)
         endif
      enddo
      end

      subroutine compare_vecsv(nmomset,iborn,res,ibornpr,cprop,iret)
      implicit none
      real * 8 ep
      parameter (ep=1d-8)
      integer nmomset,iborn,ibornpr,iret,j,k
      real * 8 res(nmomset,iborn),cprop,rat
      do j=1,iborn-1
         rat=res(1,iborn)/res(1,j)
         do k=1,nmomset
            if(abs(1-res(k,iborn)/res(k,j)/rat).gt.ep) goto 10
         enddo
         if(abs(1-rat).lt.ep) then
            iret=0
            cprop=1
         else
            iret=1
            cprop=rat
         endif
         ibornpr=j
         return
 10      continue
      enddo
      iret=-1
      end


      subroutine printvirtequiv
c When invoked after the first call to virtuals,
c it prints the set of equivalent virtual configurations
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer equivto(maxprocborn)
      common/cequivtovirt/equivto
      real * 8 equivcoef(maxprocborn)
      common/cequivcoefvirt/equivcoef
      integer j,k,iun,count
      save count
      data count/0/
      call newunit(iun)
      open(unit=iun,file='virtequiv',status='unknown')
      write(*,*) 'Writing virtequiv file...'
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
      ! 1000022 not correctly written out in bornequiv file, i4 -> i7):
 100  format(d11.4,5x,i8,5x,100(i8,1x))
      end

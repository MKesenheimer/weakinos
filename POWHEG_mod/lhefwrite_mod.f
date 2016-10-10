c MK: copied and modified version of lhefwrite.f, revision 3154
c modified on the basis of disquark/lhefwrite.f
c changes marked with "! CH, MK:" or "! MK:"
c bigger changes over a whole section are marked with !===...

c...lhefheader(nlf)
c...writes initialization information to a les houches events file on unit nlf. 
      subroutine lhefwritehdr(nlf)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'pwhg_flg.h'
      include 'pwhg_lhrwgt.h'
      integer nlf
      real * 8 version
      real * 8 powheginput ! CH: added
      external powheginput
      common/cversion/version
      data version/1.0/
      integer ipr,iran,n1ran,n2ran,iun,iret
      character * 3 whichpdfpk
      external whichpdfpk
      character(100) buffer
      include 'LesHouches.h'
      call  pwhg_io_write(nlf,'<LesHouchesEvents version="3.0">')
c      write(nlf,'(a)') '<LesHouchesEvents version="3.0">'
      call  pwhg_io_write(nlf,'<!--')
c     write(nlf,'(a)') '<!--'
      call  pwhg_io_write(nlf,'file generated with POWHEG-BOX-V2')
c      write(nlf,'(a,f3.1)') 'file generated with POWHEG-BOX-V2'
      iun = nlf
      include 'svn.version'
      call  pwhg_io_write(nlf,'Input file powheg.input contained:')
c      write(nlf,'(a)') 'Input file powheg.input contained:'
      call wrtpowheginput(nlf)
      call  pwhg_io_write(nlf,'End of powheg.input content')
c      write(nlf,'(a)') 'End of powheg.input content'
      write(buffer,'(a)') 'PDF package: '//whichpdfpk()
      call  pwhg_io_write(nlf,trim(buffer))
      call rm48ut(iran,n1ran,n2ran)
      write(buffer,*) 'Random number generator initialized with: ',
     1     iran,' ',n1ran,' ',n2ran
      call  pwhg_io_write(nlf,trim(buffer))
      call  pwhg_io_write(nlf,'-->')
c     write(nlf,'(a)') '-->'
      ! CH, MK: add output of event-number
      ! ================================================================
      call  pwhg_io_write(nlf,'<header>')
      call  pwhg_io_write(nlf,'<MGGenerationInfo>')
      write(buffer,*) '#  Number of Events        :',
     1                 powheginput('numevts')
      call  pwhg_io_write(nlf,trim(buffer))
      call  pwhg_io_write(nlf,'</MGGenerationInfo>')
      call  pwhg_io_write(nlf,'</header>')
      ! ================================================================
      if(flg_rwl .or. (flg_reweight.and.lhrwgt_id.ne.' ')) then
         call  pwhg_io_write(nlf,'<header>')
c     write(nlf,'(a)') '<header>'
         if(flg_rwl) then
            call rwl_write_rwgt_info(nlf)
         else
            call printrwghthdr(nlf)
         endif
         call  pwhg_io_write(nlf,'</header>')
c     write(nlf,'(a)') '</header>'
      endif
      call  pwhg_io_write(nlf,'<init>')
c      write(nlf,'(a)') '<init>'
      write(buffer,110) idbmup(1),idbmup(2),ebmup(1),ebmup(2),
     1     pdfgup(1),pdfgup(2),pdfsup(1),pdfsup(2),idwtup,nprup
      call  pwhg_io_write(nlf,trim(buffer))
      do 100 ipr=1,nprup
         write(buffer,120) xsecup(ipr),xerrup(ipr),xmaxup(ipr),
     &        lprup(ipr)
         call  pwhg_io_write(nlf,trim(buffer))
 100  continue
      call  pwhg_io_write(nlf,'</init>')
c      write(nlf,'(a)') '</init>'
 110  format(1p,2(1x,i8),2(1x,e12.5),6(1x,i6))
 120  format(1p,3(1x,e12.5),1x,i7) ! MK: changed here i6 -> i7
      end


      subroutine printleshouches
c useful for debugging
      call lhefwritev(6)
      end

c...lhefeader(nlf)
c...writes event information to a les houches events file on unit nlf. 
      subroutine lhefwritev(nlf)
      implicit none
      integer nlf
      include 'LesHouches.h'
      include 'pwhg_flg.h'
      include 'pwhg_lhrwgt.h'
      integer i,j,iret
      integer, save :: counter=0
      character(200) buffer
      if(flg_noevents) then
c do not write events, write only the event count
         counter = counter + 1
         if((counter/1000)*1000.eq.counter) then
            write(nlf,*) counter
         endif
         return
      endif
      call  pwhg_io_write(nlf,'<event>')
c      write(nlf,'(a)')'<event>'
      write(buffer,210) nup,idprup,xwgtup,scalup,aqedup,aqcdup
      call  pwhg_io_write(nlf,trim(buffer))
      do 200 i=1,nup
         ! MK: added gluon color check
         if((idup(i).eq.21) .and. 
     &      ((icolup(1,i).eq.0) .or. (icolup(2,i).eq.0))) then
           print*,"error in lhefwrite_mod.f:94"
           print*,"idup is gluon, but has too less icolup entries:"
           print*,"icolup(:,",i,")",icolup(:,i)
           stop
         endif
         write(buffer,220) idup(i),istup(i),mothup(1,i),
     & mothup(2,i),icolup(1,i),icolup(2,i),(pup(j,i),j=1,5),
     & vtimup(i),spinup(i)
         call  pwhg_io_write(nlf,trim(buffer))
 200  continue
      if(flg_reweight) then
         call lhefwriteevrw(nlf)
         if(flg_rwl) then
            call rwl_write_weights(nlf)
         elseif(lhrwgt_id.ne.' ') then
            call  pwhg_io_write(nlf,'<rwgt>')
c            write(nlf,'(a)')'<rwgt>'
            call printrwgtev(nlf,xwgtup)
            call  pwhg_io_write(nlf,'</rwgt>')
c            write(nlf,'(a)')'</rwgt>'
         endif
      endif
      if(flg_pdfreweight) call lhefwritepdfrw(nlf)
      if(flg_debug) call lhefwritextra(nlf)
      call  pwhg_io_write(nlf,'</event>')
c      write(nlf,'(a)')'</event>'
 210  format(1p,2(1x,i7),4(1x,e12.5)) ! MK: changed here i6 -> i7
 220  format(1p,i8,5(1x,i5),5(1x,e16.9),1x,e12.5,1x,e10.3)
      end

c...lheftrailer(nlf)
c...writes trailer to a les houches events file on unit nlf. 
      subroutine lhefwritetrailer(nlf)
      implicit none
      integer nlf,iran,n1ran,n2ran,iret
      character(100) buffer
      call pwhg_io_write(nlf,'</LesHouchesEvents>')
c     save last random number
      call rm48ut(iran,n1ran,n2ran)
      write(buffer,*) '#Random number generator exit values: ',
     # iran,' ',n1ran,' ',n2ran
      call  pwhg_io_write(nlf,trim(buffer))
      end

      subroutine lhefwritextra(nlf)
      implicit none
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      include 'pwhg_flg.h'
      integer nlf
      integer lh_seed,lh_n1,lh_n2,iret
      common/lhseeds/lh_seed,lh_n1,lh_n2
      character(100) buffer
c      write(nlf,'(a)') '# Start extra-info-previous-event'
      call pwhg_io_write(nlf,'# Start extra-info-previous-event')
      write(buffer,*) '# ',rad_kinreg,'       rad_kinreg'
      call pwhg_io_write(nlf,trim(buffer))
      write(buffer,*) '# ',rad_type,'         rad_type'
      call pwhg_io_write(nlf,trim(buffer))
      write(buffer,*) '# ', lh_seed,' ',lh_n1,' ',lh_n2,
     #     "    previous event's random seeds "
      call pwhg_io_write(nlf,trim(buffer))
c      write(nlf,'(a)') '# End extra-info-previous-event'
      call pwhg_io_write(nlf,'# End extra-info-previous-event')
      end

      subroutine lhefwritepdfrw(nlf)
      implicit none
      integer nlf
      integer id1,id2,iret
      real * 8 x1,x2,xf1,xf2,xmufact
      character(100) buffer
      call pdfreweightinfo(id1,id2,x1,x2,xmufact,xf1,xf2)
      write(buffer,111)'#pdf ',id1,id2,x1,x2,xmufact,xf1,xf2
      call pwhg_io_write(nlf,trim(buffer))
 111  format(a,2(1x,i2),5(1x,e14.8))
      end


      subroutine lhefwriteevrw(nlf)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'pwhg_rwl.h'
      integer nlf
      character * 132 string
c     rad_type=1,2,3 for btilde,remnants,regulars, respectively
      write(string,*)'#rwgt ',rwl_type,
     $        rwl_index,rwl_weight,rwl_seed,rwl_n1,rwl_n2
c This gymnastics to avoid some fortran compiler going automatically to a new line
c when writing too long records with fmt=*
      call pwhg_io_write(nlf,trim(adjustl(string)))
      end

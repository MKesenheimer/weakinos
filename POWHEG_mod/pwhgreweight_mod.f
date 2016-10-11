c MK: copied and modified version of pwhgreweight.f, revision 3154
c modified on the basis of disquark/pwhgreweight.f
c changes marked with "! CH, MK:" or "! MK:"
c bigger changes over a whole section are marked with !===...

      subroutine pwhgnewweight(iunin,iunrwgt)
      implicit none
      integer iunin,iunrwgt
      character * 500 stringin
      integer iret
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'pwhg_flg.h'
      include 'pwhg_st.h'
      include 'pwhg_pdf.h'
      include 'LesHouches.h'
c ! MK: added
#include "osres.h"
#include "Flags.h"
#include "pwhg_flst_add.h"
#include "pwhg_rad_add.h"
      integer maxev
      integer k
      integer gen_seed,gen_n1,gen_n2
      common/cgenrand/gen_seed,gen_n1,gen_n2
      real * 8 newweight
      logical pwhg_isfinite
      external pwhg_isfinite
      character * 3 whichpdfpk
      character * 200 string
      real * 8 powheginput
      external powheginput
      call lhefreadevlhrwgt(iunin,iunrwgt,iret,stringin)
      if(iret.lt.0) then
         write(*,*) ' End of event file! Aborting ...'
         call exit(-1)
      endif
      call setrandom(gen_seed,gen_n1,gen_n2)
      if(rad_type.eq.1) then
         call gen_btilderw
         newweight=rad_btilde_arr(rad_ubornidx)*
     1        rad_btilde_sign(rad_ubornidx)
         if(flg_fullrwgt) then
            call fullrwgt(newweight)
         endif
      elseif(rad_type.eq.2) then
         call gen_sigremnantrw
         newweight=rad_damp_rem_arr(rad_realalr)
      elseif(rad_type.eq.3) then
         call gen_sigremnantrw
         newweight=rad_reg_arr(rad_realreg)
      ! MK: added the following statements
      !=================================================================
      ! Note: since pwhgnewweight uses indices rad_type=1,2,3 we have to shift
      ! the indices for the on-shell resonances by 3.
      ! => rad_type = ichan + 3 
      elseif((rad_type.ge.4) .and. (rad_type.le.(nosres+3))) then
         call gen_sigosresrw
         newweight=rad_osres_arr(rad_realosres,rad_type-3)
      !=================================================================
      else
         write(*,*) 'Error in pwhgnewweight, invalid rad_type: ',
     $        rad_type
         call exit(-1)
      endif

      if(.not.pwhg_isfinite(newweight)) newweight=0d0
      write(string,*) '#new weight,renfact,facfact,pdf1,pdf2',
     1        xwgtup*newweight/rad_currentweight,st_renfact,
     2        st_facfact,pdf_ndns1,pdf_ndns2,' ',whichpdfpk()
      call pwhg_io_write(iunrwgt,trim(adjustl(string)))
c      write(iunrwgt,'(a)') trim(adjustl(string))

      if(adjustl(stringin).eq.'<rwgt>') then
         call pwhg_io_write(iunrwgt,trim(adjustl(stringin)))
c         write(iunrwgt,'(a)') trim(stringin)
         do k=1,1000000
            call pwhg_io_read(iunin,stringin,iret)
            if(iret /= 0) goto 999
c            read(unit=iunin,fmt='(a)',end=999,err=999) stringin
            if(stringin.ne.'</rwgt>') then
               call pwhg_io_write(iunrwgt,trim(stringin))
c               write(iunrwgt,'(a)') trim(stringin)
            else
               call lhrwgt_writeweight
     1              (iunrwgt,xwgtup*newweight/rad_currentweight)
               call pwhg_io_write(iunrwgt,trim(stringin))
c               write(iunrwgt,'(a)') trim(stringin)
               return
            endif
         enddo
      else
         call pwhg_io_write(iunrwgt,trim(stringin))
c         write(iunrwgt,'(a)') trim(stringin)
      endif
      return
 999  continue
      write(*,*) 'Did not find </rwgt> after <rwgt>'
      write(*,*) 'in Les Houches file. Exiting ...'
      call exit(-1)
      end


      subroutine gen_btilderw
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      integer mcalls,icalls
      include 'cgengrids.h'
      real * 8 xx(ndiminteg)      
      real * 8 btilde
      external btilde
      mcalls=0
      icalls=0
      call gen(btilde,ndiminteg,xgrid,ymax,ymaxrat,xmmm,ifold,2,
     #    mcalls,icalls,xx)
      end

      subroutine gen_sigremnantrw
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'cgengrids.h'
      real * 8 xx(ndiminteg)
      integer mcalls,icalls
      real * 8 sigremnant
      external sigremnant
      mcalls=0
      icalls=0
      call gen(sigremnant,ndiminteg,xgridrm,ymaxrm,ymaxratrm,
     1 xmmmrm,ifoldrm,2,mcalls,icalls,xx)
      end

      ! CH, MK: added this routine
      subroutine gen_sigosresrw
      implicit none
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'cgengrids.h'
c MK: added
#include "pwhg_flst_add.h"
#include "Flags.h"
#include "cgengrids_add.h"
      real * 8 xx(ndiminteg)
      integer mcalls,icalls
      real * 8 sigosres
      external sigosres
      mcalls=0
      icalls=0
      call gen(sigosres,ndiminteg,xgridosres,ymaxosres,ymaxratosres,
     &         xmmmosres,ifoldosres,2,mcalls,icalls,xx)
      end      
      
      subroutine openoutputrw(iunrwgt)
      implicit none
      include 'pwhg_flg.h'
      include 'pwhg_rnd.h'
      integer iunrwgt
      character * 20 pwgprefix
      character * 200 file
      integer lprefix,iret
      common/cpwgprefix/pwgprefix,lprefix
      real * 8 powheginput
      logical exist
      if(rnd_cwhichseed.ne.'none') then
         file=pwgprefix(1:lprefix)//'events-rwgt-'
     1        //rnd_cwhichseed//'.lhe'
      else
         file=pwgprefix(1:lprefix)//'events-rwgt.lhe'
      endif
      if(powheginput("#clobberlhe").ne.1) then
         inquire(file=file,exist=exist)
         if(exist) then
            write(*,*) 'pwhg_main: error, file ',trim(file),
     1           ' exists! will not overwrite, exiting ...'
            call exit(-1)
         endif
      endif
      call pwhg_io_open_write(trim(file),iunrwgt,
     1     flg_compress_lhe,iret)
      if(iret /= 0) then
         write(*,*) ' could not open ',trim(file),' for writing'
         write(*,*) ' exiting ...'
         call exit(-1)
      endif
      end

c...reads event information from a les houches events file on unit nlf. 
      subroutine lhefreadevlhrwgt(nlf,nuo,iret,string)
      implicit none
      integer nlf,nuo,iret
      character * 500 string
      include 'LesHouches.h'
      integer i,j
      iret=0
 1    continue
      string=' '
      call pwhg_io_read(nlf,string,iret)
      if(iret /= 0) goto 666
c     read(nlf,fmt='(a)',err=777,end=666) string
      call pwhg_io_write(nuo,trim(string))
c      write(nuo,'(a)') trim(string)
      if(string.eq.'</LesHouchesEvents>') then
         goto 998
      endif
      if(string(1:6).eq.'<event') then
c on error try next event. The error may be cause by merging
c truncated event files. On EOF return with no event found
         call pwhg_io_read(nlf,string,iret)
         if(iret /= 0) goto 666
c         read(nlf,fmt='(a)',err=777,end=666) string
         call pwhg_io_write(nuo,trim(string))
c         write(nuo,'(a)') trim(string)
         read(string,fmt=*,end=998,err=1)
     1        nup,idprup,xwgtup,scalup,aqedup,aqcdup
         do i=1,nup
            call pwhg_io_read(nlf,string,iret)
            if(iret /= 0) goto 666
c            read(nlf,'(a)',err=777,end=666) string
            call pwhg_io_write(nuo,trim(string))
c            write(nuo,'(a)') trim(string)
            read(string,fmt=*,end=998,err=1)
     1           idup(i),istup(i),mothup(1,i),
     &           mothup(2,i),icolup(1,i),icolup(2,i),(pup(j,i),j=1,5),
     &           vtimup(i),spinup(i)
         enddo
         call lhefreadextrarw(nlf,nuo,iret,string)
c Go on reading, find the end of the weights section
         goto 999
      else
         goto 1
      endif
c no event found:
 777  continue
      write(*,*) "Error in reading"
      write(*,*) string
      !call exit(-1) ! MK: changed according to lhefread.f
      nup=0
      return
 666  continue
      !iret=-1 ! MK: changed
      !return
      print *,"reached EOF"
      print *,string
      nup=0
      return
 998  continue
      print *,"read </LesHouchesEvents>"
      iret=-1
      nup=0      
 999  end


      subroutine lhefwritetrailernw(iunin,iunout)
      implicit none
      integer iunin,iunout
      character * 500 string
      integer iret
 1    continue
      string=' '
      call pwhg_io_read(iunin,string,iret)
      if(iret /= 0) goto 998
c      read(iunin,fmt='(a)',err=998,end=998) string
      call pwhg_io_write(iunout,trim(string))
c      write(iunout,'(a)') trim(string)
      if(string.eq.'</LesHouchesEvents>') then
         goto 998
      endif
      goto 1
 998  continue
      end



      subroutine lhefreadextrarw(nlf,nou,iret,string)
      implicit none
      include 'LesHouches.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'pwhg_st.h'
      include 'pwhg_kn.h'
      include 'pwhg_flg.h'
c !  MK: added
#include "osres.h"
#include "pwhg_flst_add.h"
#include "pwhg_rad_add.h"
      character * 500 string
      integer nlf,nou,iret
      logical readrw
      integer gen_seed,gen_n1,gen_n2
      common/cgenrand/gen_seed,gen_n1,gen_n2
      readrw = .false.
 1    continue
      call pwhg_io_read(nlf,string,iret)
      if(iret /= 0) goto 998
c      read(unit=nlf,fmt='(a)',end=998) string
      if(string.eq.'<rwgt>'.or.
     1     string.eq.'</event>') then
c Don't write the end event record; first we must output the new weight
         return
      endif
      if(string.eq.'<event>') then
         if(.not.readrw) then
            write(*,*) 
     $ 'Error in lhefreadextra, while reading rwg informations'
            write(*,*)'Abort run'
            call exit(-1)
         endif
         backspace nlf
         return
      endif
      if(string.eq.'# Start extra-info-previous-event') then
         call pwhg_io_write(nou,trim(string))
c         write(nou,'(a)') trim(string)
         call pwhg_io_read(nlf,string,iret)
c         read(nlf,'(a)') string
         call pwhg_io_write(nou,trim(string))
c         write(nou,'(a)') trim(string)
         read(string(3:),*) rad_kinreg
         call pwhg_io_read(nlf,string,iret)
c         read(nlf,'(a)') string
         read(string(3:),*) rad_type
      endif
      call pwhg_io_write(nou,trim(string))
c      write(nou,'(a)') trim(string)
      if(flg_newweight) then
c read a string; if it starts with #rwgt, read first rad_type from the
c string, then all other information, depending upon rad_type.
c set readrw to true
         string=adjustl(string)
         if(string(1:5).eq.'#rwgt') then
            string(1:5)=' '
c     do things
c            print*, 'FOUND'
            read(string,*) rad_type
            if(rad_type.eq.1) then
c     btilde
               read(string,*)rad_type,
     $              rad_ubornidx,rad_currentweight,
     $              gen_seed,gen_n1,gen_n2
            elseif(rad_type.eq.2) then
c     remnant
               read(string,*)rad_type,
     $              rad_realalr,rad_currentweight,
     $              gen_seed,gen_n1,gen_n2
            elseif(rad_type.eq.3) then
c     regular
               read(string,*)rad_type,
     $              rad_realreg,rad_currentweight,
     $              gen_seed,gen_n1,gen_n2
            ! MK: added the following statements
            !===========================================================
            elseif((rad_type.ge.4) .and. (rad_type.le.(3+nosres))) then
               read(string,*)rad_type,
     $              rad_realosres,rad_currentweight,
     $              gen_seed,gen_n1,gen_n2
            !===========================================================
            else
               write(*,*) 'Invalid rad_type in lhefwriteevrw: ',rad_type
               call exit(-1)
            endif
c     if all went ok, set readrw to true
            readrw=.true.
         endif
      endif
      goto 1
 998  continue
      iret=-1
      end


c The following are routines for the implementation of the
c Les Houches standard for reweighting


      subroutine lhrwgt_writeweight(iun,weight)
      implicit none
      integer iun
      include 'pwhg_lhrwgt.h'
      double precision weight
      character * 100 string
      integer iret
      write(string,'(a,e11.5,a)')
     1     "<wgt id='"//trim(adjustl(lhrwgt_id))//"'> ",
     2     weight,' </wgt>'
      call pwhg_io_write(iun,trim(string))
      end

      subroutine lhrwgt_writeheader(iun)
      implicit none
      integer iun
      include 'pwhg_lhrwgt.h'
      character * 200 tmpstring
      character * 200 string
      integer j,iret
      logical written,in_group
c This flag is set to true when the current weight has
c been written out
      written = .false.
      call pwhg_io_write(iun,'<initrwgt>')
c      write(iun,'(a)') '<initrwgt>'
      if(lhrwgt_group_name.eq.' ') then
c Write all lines firs
         do j=1,lhrwgt_nheader
            call pwhg_io_write(iun,trim(lhrwgt_header(j)))
c           write(iun,'(a)') trim(lhrwgt_header(j))
         enddo
      else
c this flag is set to true when we are writing out weights in
c the current group
         in_group = .false.
         do j=1,lhrwgt_nheader
            if(.not.written.and.in_group
     1           .and.adjustl(lhrwgt_header(j)).eq.'</weightgroup>')
     2           then
               write(string,'(a)') "<weight id='"//
     1              trim(adjustl(lhrwgt_id))//
     2              "'>"//' '//trim(adjustl(lhrwgt_descr))
     3              //' </weight>'
               call pwhg_io_write(iun,trim(string))
               in_group = .false.
               written = .true.
            endif
            call pwhg_io_write(iun,trim(lhrwgt_header(j)))
c            write(iun,'(a)') trim(lhrwgt_header(j))
            tmpstring = adjustl(lhrwgt_header(j))
            if(tmpstring(1:18).eq.'<weightgroup name=') then
               tmpstring = tmpstring(19:)
               call getquotedstring(tmpstring,tmpstring,iret)
               if(tmpstring.eq.lhrwgt_group_name) in_group = .true.
            endif
         enddo
      endif
      if(.not.written) then
         if(lhrwgt_group_name.ne.' ') then
            if(lhrwgt_group_combine.ne.' ') then
               write(string,'(a)') "<weightgroup name='"
     1              //trim(lhrwgt_group_name)//"' combine='"
     1              //trim(lhrwgt_group_combine)//"'>"
               call pwhg_io_write(iun,trim(string))
            else
               write(string,'(a)') "<weightgroup name='"
     1              //trim(lhrwgt_group_name)//"'>"
               call pwhg_io_write(iun,trim(string))
            endif
         endif
         write(string,'(a)') "<weight id='"//
     1        trim(adjustl(lhrwgt_id))//
     1        "'>"//' '//trim(adjustl(lhrwgt_descr))//' </weight>'
         call pwhg_io_write(iun,trim(string))
         if(lhrwgt_group_name.ne.' ') then
            call pwhg_io_write(iun,'</weightgroup>')
         endif
      endif
      call pwhg_io_write(iun,'</initrwgt>')
c      write(iun,'(a)') '</initrwgt>'
      end

      subroutine lhrwgt_checkheader
      implicit none
      include 'pwhg_lhrwgt.h'
      character * (lhrwgt_max_header_columns) string
      integer j,iret
      do j=1,lhrwgt_nheader
         string = adjustl(lhrwgt_header(j))
         if(string(1:11).eq.'<weight id=') then
            call getquotedstring(string(12:),string,iret)
            if(string.eq.lhrwgt_id) then
               write(*,*) ' lhrwgt_checkheader: ERROR'
               write(*,*) ' lhrwgt_id ',trim(string),' already in use'
               call exit(-1)
            endif
         endif
      enddo
      end

      subroutine lhrwgt_readheader_rw(iun)
      implicit none
      include 'pwhg_lhrwgt.h'
      integer iun
      character * (lhrwgt_max_header_columns+200) string,string0
      integer j,k,iret
      do j=1,1000000
         call pwhg_io_read(iun,string0,iret)
c     read(iun,'(a)') string0
         string=adjustl(string0)
         if(string.eq.'<initrwgt>') then
            lhrwgt_nheader = 0
            do k=1,1000000
               call pwhg_io_read(iun,string0,iret)
c               read(iun,'(a)') string0
               string=adjustl(string0)
               if(string.ne.'</initrwgt>') then
                  lhrwgt_nheader = lhrwgt_nheader + 1
                  if(lhrwgt_nheader.gt.lhrwgt_maxnheader) then
                     write(*,*) '******** ERROR ******'
                     write(*,*) ' lhrwgt_readheader_rw: '
                     write(*,*) 
     1                    ' found too many strings in header'
                     write(*,*) ' increase lhrwgt_maxnheader'
     1                    //' in include/pwg_lhrwgt.h'
                     call exit(-1)
                  endif                    
c check that the string fits
                  if(len(trim(string0)).gt.
     1                 lhrwgt_max_header_columns) then
                     write(*,*) '******** ERROR ******'
                     write(*,*) ' lhrwgt_readheader_rw: '
                     write(*,*) 
     1                    ' found too long a string in the header'
                     write(*,*) 
     1                    ' increase lhrwgt_max_header_columns'
     1                    //' in include/pwg_lhrwgt.h'
                     call exit(-1)
                  endif
                  lhrwgt_header(lhrwgt_nheader) = trim(string0)
               else
                  return
               endif
            enddo
            write(*,*) '******** ERROR ******'
            write(*,*) ' lhrwgt_readheader_rw: '
            write(*,*) ' didn;t find the end of the header'
            call exit(-1)
         endif
      enddo
      write(*,*) '******** ERROR ******'
      write(*,*) ' lhrwgt_readheader_rw: '
      write(*,*) ' didn;t find the header'
      call exit(-1)
      end

      subroutine lhrwgt_copyheader(iunin,iunout)
      implicit none
      include 'pwhg_lhrwgt.h'
      integer iunin,iunout
      character * (400) string,string0
      integer j,k,iret
      logical written,in_group
      do j=1,1000000
         call pwhg_io_read(iunin,string0,iret)
         if(iret /= 0) goto 999
c         read(unit=iunin,fmt='(a)',end=999,err=999) string0
         string=adjustl(string0)
         if(string.eq.'<initrwgt>') then
            call pwhg_io_backspace(iunin)
            call lhrwgt_readheader_rw(iunin)
            call lhrwgt_checkheader
            call lhrwgt_writeheader(iunout)
            return
         else
            call pwhg_io_write(iunout,trim(string0))
c     write(iunout,'(a)') trim(string0)
         endif
      enddo
 999  continue
      write(*,*) ' ERROR: did not find <initrwgt> in lhe file'
      write(*,*) ' Most likely, the initial run did not specify'
      write(*,*) ' a rwgt_id line'
      call exit(-1)
      end


      subroutine getreweightinput
      implicit none
      include 'pwhg_lhrwgt.h'
      call powheginputstring("#lhrwgt_group_name",lhrwgt_group_name)
      call powheginputstring("#lhrwgt_group_combine",
     1     lhrwgt_group_combine)
      call powheginputstring("#lhrwgt_id",lhrwgt_id)
      call powheginputstring("#lhrwgt_descr",lhrwgt_descr)
      end

      subroutine printrwghthdr(iun)
      implicit none
      include 'pwhg_lhrwgt.h'
      INTEGER IUN
      if(lhrwgt_id.ne.' ') then
c     write(iun,'(a)') '<initrwgt>'
         call pwhg_io_write(iun,'<initrwgt>')
         if(lhrwgt_group_name.ne.' ') then
            if(lhrwgt_group_combine.ne.' ') then
               call pwhg_io_write(iun,"<weightgroup name='"
     1              //trim(lhrwgt_group_name)//"' combine='"
     1              //trim(lhrwgt_group_combine)//"'>")
c               write(iun,'(a)') "<weightgroup name='"
c     1              //trim(lahrwgt_group_name)//"' combine='"
c     1              //trim(lhrwgt_group_combine)//"'>"
            else
               call pwhg_io_write(iun,"<weightgroup name='"
     1              //trim(lhrwgt_group_name)//"'>")
c               write(iun,'(a)') "<weightgroup name='"
c     1              //trim(lhrwgt_group_name)//"'>"
            endif
         endif
         call pwhg_io_write(iun,"<weight id='"//trim(lhrwgt_id)
     1       //"'> "//trim(lhrwgt_descr)//" </weight>")
c         write(iun,'(a)') "<weight id='"//trim(lhrwgt_id)
c     1       //"'> "//trim(lhrwgt_descr)//" </weight>"
         if(lhrwgt_group_name.ne.' ') then
             call pwhg_io_write(iun,"</weightgroup>")
c     write(iun,'(a)') "</weightgroup>"
          endif
          call pwhg_io_write(iun,"</initrwgt>")
c         write(iun,'(a)') "</initrwgt>"
      endif
      end

      subroutine printrwgtev(nlf,weight)
      implicit none
      integer nlf
      double precision weight
      character(len=300) string
      include 'pwhg_lhrwgt.h'
      write(string,'(a,e16.9,a)')"<wgt id='"//trim(lhrwgt_id)//"'> ",
     1     weight,' </wgt>'
      call pwhg_io_write(nlf,trim(string))
      end
c From here to the end of file are the subroutines needed for full reweighting,
c that also includes the effect of the change in pdf's in the hardest radiation

      subroutine readpowheginputinfo(nlf)
c If any information from the current lhe file is desired, get it here.
c In this case, we set the pdf's used for generating the file, that may
c be needed for reweighting.
      implicit none
      integer nlf
      character * 300 string
      character * 3 whichpdf
      include 'pwhg_pdf.h'
      integer iret,j,ndns1,ndns2,lhans1,lhans2
      call pwhg_io_rewind(nlf)
      do j=1,1000000
         call pwhg_io_read(nlf,string,iret)
         if(iret /= 0) goto 999
         string = adjustl(string)
         if(string(1:5) .eq. 'ndns1') then
            read(string(6:),*) ndns1
         elseif(string(1:5) .eq. 'ndns2') then
            read(string(6:),*) ndns2
         elseif(string(1:6) .eq. 'lhans1') then
            read(string(7:),*) lhans1
         elseif(string(1:6) .eq. 'lhans2') then
            read(string(7:),*) lhans2
         elseif(string(1:12) .eq. 'PDF package:') then
            if(adjustl(string(13:)) .eq. 'lha') then
               pdf_ndns1lhe = lhans1
               pdf_ndns2lhe = lhans1
               whichpdf = 'lha'
            elseif(adjustl(string(13:)) .eq. 'mlm') then
               pdf_ndns1lhe = ndns1
               pdf_ndns2lhe = ndns2
               whichpdf = 'mlm'
            else
               write(*,*) ' readpowheginputinfo: cannot identify'
               write(*,*) ' psf package ',trim(adjustl(string))
               call exit(-1)
            endif
            call pwhg_io_rewind(nlf)
            if(pdf_ndns1lhe .ne. pdf_ndns2lhe) then
               write(*,*)
     1      ' readpowheginputinfo: got two different pdf sets',
     2     ' for the two imcoming hadrons: ',
     3              pdf_ndns1lhe, pdf_ndns2lhe
               write(*,*) ' cannot handle that! exiting ...'
               call exit(-1)
            endif
            write(*,*) 
     1         ' readpowheginputinfo: got ',whichpdf,' set ',
     1         pdf_ndns1lhe,' from input lhe file: '
            return
         endif
      enddo
 999  continue
      write(*,*) ' readpowheginputinfo: something wrong reading'
      call exit(-1)
      end

      subroutine fullrwgt(weight)
c This subroutines reweights also the R/B exp[-int R/B] term in the
c Sudakov form factor.
c It needs to now the pdf's used in the original .lhe
c file. These are provided in a powheg.input "pdforig" entry.
c it computes all reweighting factors needed to update the R/B
c and Sudakov form factor (the last one in an approximate way).
      implicit none
      real * 8 weight
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      include 'pwhg_flg.h'
      include 'pwhg_st.h'
      include 'pwhg_pdf.h'
      include 'LesHouches.h'
c These are needed to call generipdfpar.
      integer iord,ih,iret
      real * 8 lam5l
      character * 2 scheme
c
      real * 8 mu2
      integer fl1b,fl2b,fl1r,fl2r
      real * 8 pdf(-pdf_nparton:pdf_nparton)
      real * 8 lumorig,lum,as,asorig,exponent,exponentorig,
     1         lumbornorig,lumborn,lumreal,lumrealorig,q2l,
     2         save_st_lambda5MSB,save_st_mufact2,x1r,x2r,fact
      real * 8 powheginput,dotp
      integer npdforig,npdf
      logical ini
      data ini/.true./
      integer mode
      save ini,npdforig,npdf,mode
      logical pwhg_isfinite
      external pwhg_isfinite

c These variables are the only ones accessed by the contained subroutines
      double precision mm,mm2,ym,ptm,q,ptm2,q2,mtm,mtm2,
     1     logtaumin,logtaumax,taumin,taumax,taubmin
      integer ngauss
      parameter (ngauss=16)
      double precision xgauss(ngauss),wgauss(ngauss)
      double precision shat,z,em,kmom,aslim,ycmmin,ycmmax,
     1  alphads

      double precision lam5,q2low,q2high
c end variables accessed by the contained subroutines

c two possibilities:
c A) HERWIG style Sudakov
c Radiation factor = Real(kt2)/Born(kt2) * Lumborn(kt2)/Lumborn(q2l) * Delta_Herwig
c The correction factor is
c alpha(kt2)/alpha_old(kt2) * (Lumreal(kt2)/Lumrealold(kt2))/(Lumborn(kt2)/Lumborn_old(kt2))
c Times (Lumborn(kt2)/Lumborn(q2l))/(Lumborn_old(kt2)/Lumborn_old(q2l)) Delta_Herwig/Delta_Herwig_old
c Equals to
c alpha(kt2)/alpha_old(kt2) * (Lumreal(kt2)/Lumreal_old(kt2)) Lumborn_old(q2l)/Lumborn(q2l)
c times Delta_Herwig/Delta_Herwig_old;
c B) Sjostrand style Sudakov
c Radiation factor = Real(kt2)/Born(kt2) * Delta_Sjostrand
c Correction factor:
c alpha(kt2)/alpha_old(kt2) * (Lumreal(kt2)/Lumreal_old(kt2))/(Lumborn(kt2)/Lumborn_old(kt2))
c Times Delta_Sjostrand/Delta_Sjostrand_old
c Equals to
c alpha(kt2)/alpha_old(kt2) * (Lumreal(kt2)/Lumreal_old(kt2)) * (Lumborn_old(kt2)/Lumborn(kt2))
c 

      if(ini) then
         mode = powheginput("#fullrwgtmode")
         if(mode.lt.0) mode = 4
         if(pdf_ndns1.ne.pdf_ndns2) then
            write(*,*) " fullpdfrwgt now works only for identical pdf's"
            write(*,*) " for the two beams ! exiting ..."
            call exit(-1)
         endif
         ini = .false.
      endif
      npdf = pdf_ndns1
      npdforig = pdf_ndns1lhe

      x1r = pup(4,1)/ebmup(1)
      x2r = pup(4,2)/ebmup(2)

      fl1b = flst_born(1,rad_ubornidx)
      fl2b = flst_born(2,rad_ubornidx)

      fl1r = idup(1)
      if(fl1r.eq.21) fl1r = 0
      fl2r = idup(2)
      if(fl2r.eq.21) fl2r = 0

c If it has Born kinematics, it is already reweighted by the standard POWHEG reweighter.
c However, it must be reweighted with the Sudakov form factor for not emitting below the cutoff scale
c Check that it has Born kinematics by requiring the same flavours and same x_1/2 for real and born.
c Reals are written to the LH file with a precision of 1d-9
      if(fl1r .eq. fl1b .and. fl2r .eq. fl2b .and.
     1     abs((kn_xb1-x1r)/x1r) .lt. 2d-9   .and. 
     2     abs((kn_xb2-x2r)/x2r) .lt. 2d-9)   then
         if(mode.eq.4) then
c we must compute the Sudakov form factor at the current scale with current pdf's
            mu2=rad_ptsqmin
            save_st_lambda5MSB = st_lambda5MSB
            save_st_mufact2 = st_mufact2
c     Sudakov, Sjostrand kind, exact phase space
            call sudakovxxx(mu2,exponent)
c     Compute the same for original pdf's
            pdf_ndns1 = npdforig
            pdf_ndns2 = npdforig
            call genericpdfpar(npdforig,ih,lam5l,scheme,iord,iret)
            st_lambda5MSB = lam5l
c     Sudakov, Sjostrand kind, exact phase space
            call sudakovxxx(mu2,exponentorig)
            fact = exp(exponent-exponentorig)
            if(pwhg_isfinite(fact)) then
               weight = weight * fact
            endif
c Set parameters affecting pdf's as in the current settings
            st_lambda5MSB = save_st_lambda5MSB
            st_mufact2 = save_st_mufact2
            pdf_ndns1 = npdf
            pdf_ndns2 = npdf
            call genericpdfpar(npdf,ih,lam5l,scheme,iord,iret)
         endif
         return
      endif
c Set scale equal to boson pt; take real kinematics from LH event;
c the kn_*real variables are undefined during reweighting.
      mu2 = pup(1,nup)**2 + pup(2,nup)**2
      q2l = dotp(kn_cmpborn(:,3)+kn_cmpborn(:,4),
     1          kn_cmpborn(:,3)+kn_cmpborn(:,4))

      save_st_mufact2 = st_mufact2
c Compute Born lum for current pdf's

      call genericpdfpar(npdf,ih,lam5l,scheme,iord,iret)

      if(mode.ne.3.and.mode.ne.4) then
c in these mode, it is assumed that HERWIG style Sudakovs is a good approximation
         st_mufact2 = max(pdf_q2min,q2l)
      else
c Sjostrand style Sudakov; must supply the ratio of born luminosities computed at kt2
         st_mufact2 = max(pdf_q2min,mu2)
         call pdfcall(1,kn_xb1,pdf)
         lumborn = pdf(fl1b)
         call pdfcall(2,kn_xb2,pdf)
         lumborn = lumborn * pdf(fl2b)
      endif

c Compute Real lum for current pdf's

      st_mufact2 =  max(pdf_q2min,mu2)
      call pdfcall(1,x1r,pdf)
      lumreal = pdf(fl1r)
      call pdfcall(2,x2r,pdf)
      lumreal = lumreal * pdf(fl2r)

c Get as to compute radiation for current pdf's
      as = mcalphas(mu2,lam5l)
      save_st_lambda5MSB = st_lambda5MSB
      st_lambda5MSB = lam5l
      if(q2l.gt.mu2) then
         if(mode.eq.1) then
c Sudakov MiNLO style (HERWIG kind)
            call sudakov_exponent(mu2,q2l,q2l,exponent,.true.,
     1                            2,.false.)
         elseif(mode.eq.2) then
c Sudakov sikmplified (HERWIG kind)
            exponent = simplesudakov(lam5l,mu2,q2l)
         elseif(mode.eq.3) then
c Sudakov sikmplified (Sjostrand kind)
            exponent =
     1           simplesudakovx(lam5l,mu2,q2l,kn_xb1,kn_xb2,fl1b,fl2b)
         endif
      else
         exponent = 0
      endif
      if(mode.eq.4) then
c Sudakov, Sjostrand kind, exact phase space
         call sudakovxxx(mu2,exponent)
         exponent = exponent/2
      endif
c Compute the same for original pdf's
      pdf_ndns1 = npdforig
      pdf_ndns2 = npdforig
      call genericpdfpar(npdforig,ih,lam5l,scheme,iord,iret)

      if(mode.ne.3.and.mode.ne.4) then
         st_mufact2 = max(pdf_q2min,q2l)
         call pdfcall(1,kn_xb1,pdf)
         lumbornorig = pdf(fl1b)
         call pdfcall(2,kn_xb2,pdf)
         lumbornorig = lumbornorig * pdf(fl2b)
      else
c Sjostrand style Sudakov; must supply the ratio of luminosities computed at kt2
         st_mufact2 = max(pdf_q2min,mu2)
         call pdfcall(1,kn_xb1,pdf)
         lumbornorig = pdf(fl1b)
         call pdfcall(2,kn_xb2,pdf)
         lumbornorig = lumbornorig * pdf(fl2b)
      endif

      st_mufact2 = max(pdf_q2min,mu2)
      call pdfcall(1,x1r,pdf)
      lumrealorig = pdf(fl1r)
      call pdfcall(2,x2r,pdf)
      lumrealorig = lumrealorig * pdf(fl2r)

c Get as to compute radiation for original pdf's
      asorig = mcalphas(mu2,lam5l)
      st_lambda5MSB = lam5l

      if(q2l.gt.mu2) then
         if(mode.eq.1) then
c Sudakov MiNLO style (HERWIG kind)
            call sudakov_exponent(mu2,q2l,q2l,exponentorig,.true.,
     1           2,.false.)
         elseif(mode.eq.2) then
c Sudakov sikmplified (HERWIG kind)
            exponentorig = simplesudakov(lam5l,mu2,q2l)
         elseif(mode.eq.3) then
c Sudakov sikmplified (Sjostrand kind)
            exponentorig =
     1           simplesudakovx(lam5l,mu2,q2l,kn_xb1,kn_xb2,fl1b,fl2b)
         endif
      else
         exponentorig = 0
      endif
      if(mode.eq.4) then
c Sudakov, Sjostrand kind, exact phase space
         call sudakovxxx(mu2,exponentorig)
         exponentorig = exponentorig/2
      endif
c Set parameters affecting pdf's as in the current settings
      st_lambda5MSB = save_st_lambda5MSB
      st_mufact2 = save_st_mufact2
      pdf_ndns1 = npdf
      pdf_ndns2 = npdf
      call genericpdfpar(npdf,ih,lam5l,scheme,iord,iret)
c Compute reweighting factor:
      
      fact =   as/asorig                               ! as for radiation updated
      if(mode.ne.3.and.mode.ne.4) then
         fact = fact * lumbornorig/lumborn             ! pdf ratio arising
                                                       ! from Sudakov form factors
      else
         fact = fact * lumbornorig/lumborn             ! Factor due to R/B in radiation
                                                       ! this time evaluated at kt2
      endif

      fact = fact * lumreal/lumrealorig                ! new pdf's in radiation
      fact = fact * exp(2*(exponent-exponentorig))     ! Sudakov form factors

      if(pwhg_isfinite(fact)) then
         weight = weight * fact
      endif

      contains


      double precision function
     1     simplesudakovx(lam5,q2low,q2high,x1,x2,fl1,fl2)
      implicit none
      double precision lam5,q2low,q2high,x1,x2
      integer fl1,fl2
      integer ngauss
      parameter (ngauss=16)
      double precision xgauss(ngauss),wgauss(ngauss),xxx(2)
      integer j,k
c      double precision simplesudakovx0,
      double precision res
      logical ini
      data ini/.true./
      save ini,xgauss,wgauss
c Setup gaussian weights and points
      if(ini) then
         call dgset(0d0,1d0,ngauss,xgauss,wgauss)
         ini = .false.
      endif
      res = 0
      do j=1,ngauss
         xxx(1) = xgauss(j)
         do k=1,ngauss
            xxx(2) = xgauss(k)
            res = res + wgauss(j)*wgauss(k)*
     1           simplesudakovx0(lam5,q2low,q2high,x1,x2,fl1,fl2,xxx)
         enddo
      enddo
      simplesudakovx = res
      end function simplesudakovx

      double precision function
     1     simplesudakovx0(lam5,q2low,q2high,x1,x2,fl1,fl2,xxx)
      implicit none
      double precision lam5,q2low,q2high,x1,x2,xxx(2)
      integer fl1,fl2,fl,j
      include 'pwhg_math.h'
      include 'pwhg_pdf.h'
      include 'pwhg_st.h'
      real * 8 kt2,x,xjac,lq2,as,eta,xi,ximax,ximin,y,ymax,z,xx
      real * 8 fxr(-pdf_nparton:pdf_nparton),
     1         fxb(-pdf_nparton:pdf_nparton)
c      real * 8 mcalphas
c      external mcalphas
      xjac = 1
      lq2 = log(q2high/q2low)
      kt2 = exp(xxx(1)*lq2)*q2low
      xjac = xjac * kt2 * lq2
      eta = 4*kt2/q2high
      simplesudakovx0 = 0
      do j=1,2
         if(j.eq.1) then
            xx = x1
            fl = fl1
         else
            xx = x2
            fl = fl2
         endif
         if((1-xx)**2.gt.eta) then
c define y=sqrt((1-x)^2-eta);
c y d y = (1-x) dx, dx = y dy /(1-x)
c y goes from 0 to sqrt((1-x1)**2-eta). Needs importance sampling
c near the upper bound, and also near the lower bound if x1 is small.
c Take xi = log((y+eta)/(1-y)) as sampling variable xi
            ymax = sqrt((1-xx)**2-eta)
            ximax = log((ymax+eta)/(1-ymax))
            ximin = log(eta)
            xi = (ximax-ximin)*xxx(2) + ximin
            xjac = xjac*(ximax-ximin)
            z=exp(xi)
            xjac = xjac * z
            y = (z-eta)/(1+z)
            xjac = xjac * (1+eta) / (1+z)**2
            x = 1 - sqrt( y**2 + eta )
            xjac = xjac * y/(1-x)
            st_mufact2 = max(pdf_q2min,kt2)
            as =  mcalphas(kt2,lam5)
            call pdfcall(1,xx,fxb)
            call pdfcall(1,xx/x,fxr)
            simplesudakovx0 = simplesudakovx0 + (
     1           ( cf * (1+x**2)/(1-x) * fxr(fl) +
     2           tf * (x**2+(1-x)**2) * fxr(0)  )/fxb(fl)   )
     3           *  as/(2*pi) /kt2/x *xjac
     4           / sqrt(1-eta/(1-x)**2)
         endif
      enddo
c we divide by two simply because we multiply by 2 in the main program
      simplesudakovx0 = - simplesudakovx0 / 2
      end function simplesudakovx0

      double precision function simplesudakov(lam50,q2low0,q2high0)
      implicit none
      double precision lam50,q2low0,q2high0
      double precision dgauss
      ! MK: added
      double precision simplesudakov0
      !external simplesudakov0
      lam5 = lam50
      q2low = q2low0
      q2high = q2high0
      simplesudakov = dgauss(simplesudakov0,0d0,1d0,1d-4)
      end function simplesudakov

      double precision function simplesudakov0(x)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'pwhg_rad.h'
      double precision x
      double precision dgauss      
      real * 8 as,lqlow,l,q2
      integer nf
      lqlow = log(q2high/q2low)
      l = x*lqlow
      q2 = q2low*exp(l)
      as = mcalphas(q2,lam5)
c Sudakov exponent is (for each leg)
c              d q2   C_f       [      Q2       ]
c      \int   -----  ---- as(q) [ log --- + B1  ],  B1 = -3/2 = -1.5
c               q2   2 pi       [      q2       ]

c lqlow is the jacobian
      simplesudakov0 = - cf/(2*pi)*as*(l-1.5d0) * lqlow

      end function simplesudakov0

      

      double precision function mcalphas(q2,lam5)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_math.h'
      include 'pwhg_rad.h'
      include 'pwhg_st.h'
      real * 8 q2
      real * 8 lam5
      integer nf
      real * 8 as,pwhg_alphas
      external pwhg_alphas
      if(q2.lt.rad_charmthr2) then
         nf=3
      elseif(st_muren2.lt.rad_bottomthr2) then
         nf=4
      else
         nf=5
      endif
      as = pwhg_alphas(q2,lam5,-1)
      as = as * (1+as/(2*pi)*((67d0/18-pi**2/6)*ca-5d0/9*nf))
      mcalphas = as
      end function mcalphas

      
      subroutine sudakovxxx(kt2max,result)
      implicit none
      real * 8 result,kt2max
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      include 'pwhg_st.h'
      include 'pwhg_pdf.h'
      real * 8 pdf1(-pdf_nparton:pdf_nparton),
     1         pdf2(-pdf_nparton:pdf_nparton)
      double precision deltasud,tau0,tau1,eta,pm(0:3)
      real * 8 dotp
      logical ini
      data ini/.true./
      save ini
c Setup gaussian weights and points
      if(ini) then
         call dgset(0d0,1d0,ngauss,xgauss,wgauss)
         ini = .false.
      endif
c the vector meson rapidity and virtuality is the same in the Born and Real kinematics
      pm = kn_pborn(:,3)+kn_pborn(:,4)
      mm2 = dotp(pm,pm)
      mm=sqrt(mm2)
      ym = 0.5d0*log((pm(0)+pm(3))/(pm(0)-pm(3)))
      ptm2 = kt2max
      ptm = sqrt(kt2max)
      q2 = kn_sbeams
      q  = sqrt(q2)
      taubmin=mm2/q2
      mtm2 = mm2+ptm2
      mtm = sqrt(mtm2)
      eta=exp(-2*abs(ym))
      tau0 = (ptm+mtm)**2/q2
      tau1 = ((eta*q2-mm2)*sqrt(eta*q2*mtm2)+eta*q2*ptm2)
     1     /(q2*(eta*q2-mtm2))
      if(tau0.lt.eta) then
         if(tau1.lt.eta) then
            taumin=tau0
            taumax=1
            logtaumin = log(tau0)
            logtaumax = 0
         else
            write(*,*) 'bmass_in_minlo: tau0<eta and tau1>eta !!!'
            write(*,*) "can't be, exiting ..."
            call exit(-1)
         endif
      else
         if(tau1.lt.tau0) then
            write(*,*) 'bmass_in_minlo: tau1<tau0 !!!'
            write(*,*) "can't be, exiting ..."
            call exit(-1)
         endif
         if(tau1.gt.1) then
            write(*,*) 'bmass_in_minlo: tau1>1 !!!'
            write(*,*) "can't be, exiting ..."
            call exit(-1)
         endif
         taumin=tau1
         taumax=1
         logtaumin = log(tau1)
         logtaumax = 0
      endif
C      alphads = st_alpha ! fixed coupling 
      call dointdeltasud(result)

c      call plotdeltasud

      end subroutine sudakovxxx

      function deltasud0(ycm,tau) result(res)
      implicit none
      include 'LesHouches.h'
      include 'pwhg_math.h'
      include 'pwhg_st.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_rad.h'
      double precision, intent(in) :: ycm,tau
      double precision res
      double precision kt2,xjac,x1,x2,taub,gg,gq,qg,att,auu,
     1     kl,alphas,x,x1b,x2b,k0
      include 'pwhg_pdf.h'
c sing1 and sing2 will stand for the singlet on the 1 and 2 incoming particles
      real * 8 pdf1(-pdf_nparton:pdf_nparton),
     1         pdf2(-pdf_nparton:pdf_nparton),sing1,sing2
      integer j,fl1r,fl2r

      real * 8 born
      integer  bflav(4)
      real * 8 ppp(0:3,4) ! Partonic COM frame momenta
      real * 8 R_t,R_tb,R_tb_int,R_tb_int_fact,R_tb_int_non_fact
      real * 8 B_t,B_tb,B_t_inf,B_tb_inf
      real * 8 powheginput
      external powheginput
      real * 8 savemfer
      real * 8 approx_result
      integer  saveafer
      logical ini
      data    ini/.true./
      save    ini
      integer approx
      save    approx
c      double precision mcalphas,as
      double precision as

      kl = em * tanh(ym-ycm)
      kt2 = kmom**2-kl**2

      as = mcalphas(kt2,st_lambda5MSB)
c
      xjac = z/(4*pi)*(mm2+kt2)/(mm2*(1+z))
c
      x1 = sqrt(tau)*exp(ycm)
      x2 = sqrt(tau)*exp(-ycm)

c      shat = q2*tau

      x = mm2/shat

      k0 = sqrt(shat)*(1-x)/2
c
C      call checkkin(x1,x2,kt2)

      st_mufact2 = max(kt2,pdf_q2min)
      call pdfcall(1,x1,pdf1)
      call pdfcall(2,x2,pdf2)

c Get underlyinng Born flavours for current event
      fl1r = flst_born(1,rad_ubornidx)
      fl2r = flst_born(2,rad_ubornidx)

c the factors (k0-kl)/(2k0) and (k0+kl)/(2k0) select the regions with backward and forward
c emission respectively

      if(fl1r .ne. 0 .and. fl2r .ne. 0) then
         res = - ( cf * (1+x**2)/(1-x) * (pdf1(fl1r)*pdf2(fl2r))
     2     +   tf * (x**2+(1-x)**2) *
     3 (pdf1(0)*pdf2(fl2r)*(k0+kl)+pdf1(fl1r)*pdf2(0)*(k0-kl))/(2*k0))
     4     *  as*(8*pi) /(kt2/(1-x)) 
      elseif(fl1r .eq. 0 .and. fl2r .eq. 0) then
         sing1 = sum(pdf1(-6:-1)) + sum(pdf1(1:6))
         sing2 = sum(pdf2(-6:-1)) + sum(pdf2(1:6))
         res = - ( 2*ca * (x/(1-x)+(1-x)/x+x*(1-x))
     2      * (pdf1(fl1r)*pdf2(fl2r)) +   cf * (1+(1-x)**2)/x *
     3 (sing1*pdf2(fl2r)*(k0+kl)+pdf1(fl1r)*sing2*(k0-kl))/(2*k0) )
     4     *  as*(8*pi) /(kt2/(1-x))
      elseif(fl1r .eq. 0 .and. fl2r .ne. 0) then
         sing1 = sum(pdf1(-6:-1)) + sum(pdf1(1:6))
         res = - (
     1        (2*ca * (x/(1-x)+(1-x)/x+x*(1-x)) * (k0+kl)
     2         + cf * (1+x**2)/(1-x) * (k0-kl))
     3        * (pdf1(fl1r)*pdf2(fl2r))
     4     +  cf * (1+(1-x)**2)/x * sing1*pdf2(fl2r)*(k0+kl)
     5     +  tf * (x**2+(1-x)**2) * pdf1(fl1r)*pdf2(0)*(k0-kl) )/(2*k0)
     6     *  as*(8*pi) /(kt2/(1-x))
      elseif(fl1r .ne. 0 .and. fl2r .eq. 0) then
         sing2 = sum(pdf2(-6:-1)) + sum(pdf2(1:6))
         res = - (
     1        (cf * (1+x**2)/(1-x) * (k0+kl)) +
     2         2*ca * (x/(1-x)+(1-x)/x+x*(1-x)) * (k0-kl) 
     3        * (pdf1(fl1r)*pdf2(fl2r))
     4     +  tf * (x**2+(1-x)**2) * pdf1(0)*pdf2(fl2r)*(k0+kl)
     5     +  cf * (1+(1-x)**2)/x * pdf1(fl1r)*sing2*(k0-kl)    )/(2*k0)
     6     *  as*(8*pi) /(kt2/(1-x))
      endif

c 8 pi as (1-x) shat/(tk uk) according to MNR, = 8 pi as (1-x)/kt2

c
c
c Born pdf
      x1b = mm/q*exp(ym)
      x2b = mm2/q2/x1b
      call pdfcall(1,x1b,pdf1)
      call pdfcall(2,x2b,pdf2)

c divide by 2*pi/q2, that arises from the Born phase space in the denominator
      res = res/(pdf1(fl1r)*pdf2(fl2r)) * xjac /(2*pi/q2)

      end function deltasud0

      function deltasud1(xx) result(res)
      implicit none
c this constructs the kinematics starting from x(1) and (x2),
c and ym, ptm, ptm2 (in phspsudint.h), computes ycm,tau 
c and a jacobian to go from xx(1:2) to ycm,tau
      double precision, intent(in) :: xx(2)
      double precision res
      double precision xjac,tau,llmin,llmax,ll,ycm,zzz
c integration limit
c      tau = exp(xx(1)**2*(logtaumax-logtaumin)+logtaumin)
c      xjac = tau*(logtaumax-logtaumin)*2*xx(1)

      llmin = log(taumin-taubmin)
      llmax = log(taumax-taubmin)
      ll = xx(1)*(llmax-llmin)+llmin
      xjac = (llmax-llmin)
      tau = exp(ll)+taubmin
      xjac = xjac * exp(ll)

c range in ycm
      shat = tau*q2
      z = mm2/shat
      em = mm*(1+z)/(2*sqrt(z))
      kmom = mm*(1-z)/(2*sqrt(z))
      aslim = atanh(sqrt(kmom**2-ptm2)/em)
      ycmmin = max(log(tau)/2,
     1     ym-aslim)
      ycmmax = min(-log(tau)/2,
     1     ym+aslim)
      if(ycmmin.gt.ycmmax) then
         write(*,*) 'deltasud1:  Worry!'
         write(*,*) 'ycmmmin>ycmmax!!'
         call exit(-1)
      endif
      zzz = 6*(xx(2)**2/2-xx(2)**3/3)
      xjac = xjac * 6*xx(2)*(1-xx(2))
      ycm = zzz*(ycmmax-ycmmin)+ycmmin
      xjac = xjac*(ycmmax-ycmmin)
      res = deltasud0(ycm,tau)*xjac
      end function deltasud1

      subroutine dointdeltasud(res)
      implicit none
      real * 8 res
      real * 8 xx(2),random
      integer j,k,ncalls
      parameter (ncalls=100000)
      res = 0
c      do j=1,ncalls
c         xx(1)=random()
c         xx(2)=random()
c         res = res + deltasud1(xx)/ncalls
c      enddo
c      return

      do j = 1,ngauss
         xx(1) = xgauss(j)
         do k = 1,ngauss
            xx(2) = xgauss(k)
            res = res + deltasud1(xx)*wgauss(j)*wgauss(k)
         enddo
      enddo
      end subroutine dointdeltasud
c$$$
c$$$      subroutine plotdeltasud
c$$$      implicit none
c$$$      real * 8 res
c$$$      integer ncalls
c$$$      parameter (ncalls=100000)
c$$$      character * 20 str
c$$$      real * 8 xx(2),random
c$$$      integer j,k
c$$$      logical ini
c$$$      data ini/.true./
c$$$      save ini
c$$$      write(str,'(f4.1,"-",i3)') ym,int(ptm)
c$$$      do j=1,15
c$$$         do k=1,2
c$$$            if(str(j:j).eq.' ') str(j:)=str(j+1:)
c$$$         enddo
c$$$      enddo
c$$$      write(*,*) '!'//str//'!'
c$$$
c$$$      call inihists
c$$$      call bookupeqbins('x1',0.01d0,0d0,1d0)
c$$$      call bookupeqbins('x2',0.01d0,0d0,1d0)
c$$$
c$$$      res = 0
c$$$      do k = 1,ncalls
c$$$         xx(1) = random()
c$$$         xx(2) = random()
c$$$         res = deltasud1(xx)
c$$$         call filld('x1',xx(1),res)
c$$$         call filld('x2',xx(2),res)
c$$$         call pwhgaccumup
c$$$      enddo
c$$$      call pwhgsetout
c$$$      call pwhgtopout("deltasudhists-"//trim(str))
c$$$      end subroutine plotdeltasud
c$$$

      subroutine checkkin(x1,x2,kt2)
      implicit none
      real * 8 x1,x2,kt2,ycm,p,p0,ppar,ycmm,shat_lcl
      if(x1.lt.0.or.x1.gt.1) then
         write(*,*) ' checkkin: x1=',x1
         call exit(-1)
      endif
      if(x2.lt.0.or.x2.gt.1) then
         write(*,*) ' checkkin: x2=',x2
         call exit(-1)
      endif
      if(kt2.lt.0) then
         write(*,*) ' checkkin: kt2=',kt2
         call exit(-1)
      endif
      ycm = 0.5d0*log(x1/x2)
      shat_lcl = q2*x1*x2
      ycmm = ym-ycm
      p0 = (shat_lcl + mm2)/(2*sqrt(shat_lcl))
      p = (shat_lcl - mm2)/(2*sqrt(shat_lcl))
      ppar = p0*tanh(ycmm)
      if(abs(p**2-kt2-ppar**2).gt.1d-6) then
         write(*,*) ' checkkin: kinematics do not check'
         call exit(-1)
      endif
      end subroutine checkkin

      end
      


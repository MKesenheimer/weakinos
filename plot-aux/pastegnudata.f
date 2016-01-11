      implicit none
c substitute for paste command in gnuplot; can paste datasets with the same comment line
      character * 300 commandline
      character * 500 outl
      character * 150 file1, file2, linef(10),line
      integer length,status,tag1,tag2,fil1(10),fil2(10),i,k,ios,iout,
     1     firstnb,firstb,find,iblank,icur,nfil
      logical eqstrings
      call get_command(commandline,length,status)
      if(status.eq.-1) then
         write(*,*) 'commandline too long;'
         call exit(-1)
      elseif(status.gt.0) then
         write(*,*) 'error reading commandline'
         call exit(-1)
      elseif(length.eq.0) then
         call usage
         call exit(-1)
      endif
c commandline(tag1:tag2) is the tag of the data set
      tag1=find(commandline,1,'[')
      if(tag1.lt.0) call usage
      tag1=tag1+1
      tag2=find(commandline,tag1,']')
      if(tag2.lt.0) call usage
      icur=tag2+1
      tag2=tag2-1
      do i=1,10
         fil1(i)=firstnb(commandline,icur)
         icur=fil1(i)
         if(icur.gt.length) then
            nfil=i-1
            exit
         endif
         fil2(i)=firstb(commandline,icur)-1
         icur=fil2(i)+1
      enddo
      do i=1,nfil
         open(unit=10+i,status='old',file=commandline(fil1(i):fil2(i)),
     1        iostat=ios)
         if(ios.lt.0) then
            write(*,*) ' file ',commandline(fil1(i):fil2(i)),
     1           ' cannot be opened'
            call exit(-1)
         endif
         do k=1,10000000
            read(unit=10+i,fmt='(a)',iostat=ios) line
            if(ios.lt.0) then
               write(*,*) 'file ',commandline(fil1(i):fil2(i)),
     1           ' ended while looking for ',commandline(tag1:tag2)
               call exit(-1)
            endif
            if(index(line,commandline(tag1:tag2)).gt.0) exit
         enddo
      enddo
c Now all files are positioned at the tag
      write(*,'(a)') commandline(tag1:tag2)
      do k=1,10000000
         do i=1,nfil
            read(unit=10+i,fmt='(a)',iostat=ios) linef(i)
            if(ios.lt.0) then
               write(*,*) 'file ',commandline(fil1(i):fil2(i)),
     1           ' ended while looking for data '
               call exit(-1)
            endif
         enddo
         iblank=0
         do i=1,nfil
            if(linef(i).eq.' ') then
               iblank=iblank+1
            endif
         enddo
         if(iblank.eq.nfil) goto 999
         if(iblank.ne.0) then
            write(*,*) ' data ended in some files, but not all'
            call exit(-1)
         endif
         call outline(linef,nfil,outl,iout)
         write(*,'(10(a))') outl(1:iout)
      enddo
 999  write(*,*)
      end

      subroutine outline(linef,nfil,outl,iout)
      implicit none
      integer nfil
      character *(*) linef(nfil)
      character *(*) outl
      character * 1 c
      integer lout,iout,ib,i,k,l
      lout=len(outl)
      iout=0
      ib=1
      do i=1,nfil
         l=len(linef(i))
         do k=1,l
            c=linef(i)(k:k)
            if(c.ne.' '.or.ib.eq.0) then
               iout=iout+1
               if(iout.gt.lout) goto 998
               outl(iout:iout)=c
               if(c.ne.' ') then
                  ib=0
               else
                  ib=1
               endif
            endif
         enddo
      enddo
      goto 999
 998  continue
      write(*,*) ' outline: output line too short'
      call exit(-1)
 999  end
               
               
      
      function firstnb(str,icur)
      implicit none
      character *(*) str
      integer firstnb,icur,l,k
      l=len(str)
      do k=icur,l
         if(str(k:k).ne.' ') then
            firstnb=k
            return
         endif
      enddo
      firstnb=l+1
      end

      function find(str,ini,substr)
      implicit none
      character *(*) str,substr
      integer find,l,k,ini
      l=len(str)
      k=index(str(ini:),substr)
      if(k.eq.0) then
         find=-10
         return
      endif
      find=ini+k-1
      end

      function firstb(str,icur)
      implicit none
      character *(*) str
      integer firstb,icur,l,k
      l=len(str)
      do k=icur,l
         if(str(k:k).eq.' ') then
            firstb=k
            return
         endif
      enddo
      firstb=l+1
      end

      function eqstrings(str1,str2)
      implicit none
      character *(*) str1,str2
      logical eqstrings
      integer l1,l2,i1,i2,i
      character * 1 c1,c2
      l1=len(str1)
      l2=len(str2)
      i1=0
      i2=0
c previous character
      c1=' '
      do i=1,1000000
         i1=i1+1
         i2=i2+1
         if(c1.eq.' ') then
            do while(str1(i1:i1).eq.' '.and.i1.le.l1)
               i1=i1+1
            enddo
            do while(str2(i2:i2).eq.' '.and.i2.le.l2)
               i2=i2+1
            enddo
         endif
         if(i1.gt.l1.and.i2.gt.l2) exit
         if(i1.gt.l1) then
            c1=' '
         else
            c1=str1(i1:i1)
         endif
         if(i2.gt.l2) then
            c2=' '
         else
            c2=str2(i2:i2)
         endif
         if(c1.ne.c2) then
            eqstrings=.false.
            return
         endif
      enddo
      eqstrings=.true.
      end

      subroutine usage
      implicit none
      write(*,*)
     1     'pastegnudata "tag name" file1 file2 [file3 ... file10]'
      end

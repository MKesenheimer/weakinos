      program write_proc_labels
      implicit none
      character*3 s
      character*60 buff
      integer maxproc,maxpart
      parameter (maxproc=999,maxpart=20)
      integer total_proc,i,j,proc_number(maxproc),npart(maxproc),
     &     partid(maxpart,maxproc),lstr(maxproc),l
      character*140 str(maxproc)
      character*2 ls
      logical need_switching
      
      
      open(unit=10, file='proc_number_real',status='old',err=32)
      read (10,*,err=32,end=32) total_proc
      close(10)
      write (*,*) 'Number of subprocesses: ',total_proc

      if (total_proc.ge.999) then
         write (*,*) 'ERROR: Too many subprocesses, cannot continue.'
         stop
      endif

      open (unit=11,file='proc_label_real',status='old',err=33)
      do i=1,total_proc
         read(11,*,err=33,end=33) proc_number(i),npart(i),
     &        (partid(j,i),j=1,npart(i))
         call convert_to_sting(npart(i),partid(1,i),str(i),lstr(i))
         write (*,*) 'Subprocess found ',proc_number(i),npart(i),
     &        (partid(j,i),j=1,npart(i)),' ',str(i)(1:lstr(i)),lstr(i)
      enddo
      close(11)
      
      l=2
      do i=npart(1),3,-1
         if (abs(partid(i,1)).gt.9) then
            l=i
            exit
         endif
      enddo
      write(ls,'(i2)') l

      if (npart(1).gt.l+1) then
         need_switching=.true.
      else
         need_switching=.false.
      endif
      
      open (unit=12,file='sreal_proc.f',status='unknown')

      write (12,*) '     subroutine sreal_proc(p,legs,wgt)'
      write (12,*) '     implicit none'
      write (12,*) '     include "nexternal.inc"'
      write (12,*) '     include "coupl.inc"'
      write (12,*)
     &     '     double precision p(0:3,nexternal),wgt'
      write (12,*) '     integer legs(nexternal),lstr'
      write (12,*) '     character*140 str'
      if (need_switching) then
         write (12,*) '     double precision P1(0:3,nexternal)'
         write (12,*) '     integer i,ic(nexternal),'
     &        //'legs1(nexternal)'
         write (12,*) '     logical mtc,even'
         write (12,*) '     '
         write (12,*) '     do i=1,nexternal'
         write (12,*) '        ic(i)=i'
         write (12,*) '     enddo'
         write (12,*) '     mtc=.false.'
         write (12,*) '10   call nexper(nexternal-'//ls//
     &        ',ic('//ls//'+1),mtc,even)'
         write (12,*) '     do i='//ls//'+1,nexternal'
         write (12,*) '        ic(i)=ic(i)+'//ls
         write (12,*) '     enddo'
         write (12,*) '     CALL SWITCHMOM (P,P1,IC,NEXTERNAL)'
         write (12,*) '     CALL SWITCHLEGS(legs,legs1,IC,NEXTERNAL)'
         write (12,*) '     '
         write (12,*)
     &        '     call convert_to_string(nexternal,legs1,str,lstr)'
      else
         write (12,*) '     '
         write (12,*)
     &        '     call convert_to_string(nexternal,legs,str,lstr)'
      endif
      write (12,*) '     '
      do i=1,total_proc
         if (proc_number(i).le.9) then
            s='00'//char(proc_number(i)+48)
         elseif(proc_number(i).le.99)then
            s='0'//char(proc_number(i)/10+48)
     &           //char(mod(proc_number(i),10)+48)
         elseif(proc_number(i).le.999) then
            s=char(proc_number(i)/100+48)
     &           //char((proc_number(i)-(proc_number(i)/100)*100)/10+48)
     &           //char(mod(proc_number(i)-
     &           (proc_number(i)/100)*100,10)+48)
         endif
         if (i.eq.1) then
            write (12,*)
     &           '     if (str.eq."'//str(i)(1:lstr(i))//'") then'
         else
            write (12,*)
     &           '     elseif (str.eq."'//str(i)(1:lstr(i))//'") then'
         endif
         if (need_switching) then
            write (12,*) '        call srealmtrx_'//s//'(p1,wgt)'
         else
            write (12,*) '        call srealmtrx_'//s//'(p,wgt)'
         endif
         write (12,*) '        goto 20'
      enddo
      write (12,*) '     endif'
      write (12,*) '     '
      if (need_switching) then
         write (12,*) '     do while(mtc)'
         write (12,*) '        do i='//ls//'+1,nexternal'
         write (12,*) '           ic(i)=ic(i)-'//ls
         write (12,*) '        enddo'
         write (12,*) '        goto 10'
         write (12,*) '     enddo'
         write (12,*) '     if (.not.mtc) then'
         write (12,*) '        write (*,*) "Error #1,'//
     &        ' in sreal_proc.f"'
         write (12,*) '        stop'
         write (12,*) '     endif'
         write (12,*) '     '
      endif
      write (12,*) '20   continue'      
      write (12,*) '     return'
      write (12,*) '     end'
      write (12,*) '     '
      write (12,*) '     '



      write (12,*) '     subroutine real_color(legs,color)'
      write (12,*) '     implicit none'
      write (12,*) '     include "nexternal.inc"'
      write (12,*) '     integer maxamps'
      write (12,*) '     parameter (maxamps=6000)'
      do i=1,total_proc
         if (proc_number(i).le.9) then
            s='00'//char(proc_number(i)+48)
         elseif(proc_number(i).le.99)then
            s='0'//char(proc_number(i)/10+48)
     &           //char(mod(proc_number(i),10)+48)
         elseif(proc_number(i).le.999) then
            s=char(proc_number(i)/100+48)
     &           //char((proc_number(i)-(proc_number(i)/100)*100)/10+48)
     &           //char(mod(proc_number(i)-
     &           (proc_number(i)/100)*100,10)+48)
         endif
      write (12,*) '     Double Precision amp2'//s//'(maxamps), jamp2'
     &        //s//'(0:maxamps)'
      write (12,*) '     common/to_Ramps_'//s//'/amp2'//s//',jamp2'//s
      enddo
      write (12,*) '     double precision jamp2cum(0:maxamps)'
      write (12,*) '     integer ICOLUP(2,nexternal,maxamps)'
      write (12,*) '     integer color(2,nexternal),'
     &     //'color1(2,nexternal)'
      write (12,*) '     double precision random,xtarget'
      write (12,*) '     external random'
      write (12,*) '     integer legs(nexternal),lstr,i,j'
      write (12,*) '     character*140 str'
      if (need_switching) write (12,*) '     integer ic(nexternal),'
     &     //'legs1(nexternal)'
      write (12,*) '     integer iflow,ifl'

      if (need_switching) write (12,*) '     logical mtc,even'
      write (12,*) '     '
      if (need_switching) then
         write (12,*) '     do i=1,nexternal'
         write (12,*) '        ic(i)=i'
         write (12,*) '     enddo'
         write (12,*) '     mtc=.false.'
         write (12,*) '10   call nexper(nexternal-'//ls//
     &        ',ic('//ls//'+1),mtc,even)'
         write (12,*) '     do i='//ls//'+1,nexternal'
         write (12,*) '        ic(i)=ic(i)+'//ls
         write (12,*) '     enddo'
         write (12,*) '     CALL SWITCHLEGS(legs,legs1,IC,NEXTERNAL)'
         write (12,*) '     '
         write (12,*)
     &        '     call convert_to_string(nexternal,legs1,str,lstr)'
      else
         write (12,*)
     &        '     call convert_to_string(nexternal,legs,str,lstr)'
      endif
      write (12,*) '     '
      do i=1,total_proc
         if (proc_number(i).le.9) then
            s='00'//char(proc_number(i)+48)
         elseif(proc_number(i).le.99)then
            s='0'//char(proc_number(i)/10+48)
     &           //char(mod(proc_number(i),10)+48)
         elseif(proc_number(i).le.999) then
            s=char(proc_number(i)/100+48)
     &           //char((proc_number(i)-(proc_number(i)/100)*100)/10+48)
     &           //char(mod(proc_number(i)-
     &           (proc_number(i)/100)*100,10)+48)
         endif
         if (i.eq.1) then
            write (12,*)
     &           '     if (str.eq."'//str(i)(1:lstr(i))//'") then'
         else
            write (12,*)
     &           '     elseif (str.eq."'//str(i)(1:lstr(i))//'") then'
         endif
            write (12,*) '        include "leshouches_R_'//s//'.inc"'
            write (12,*) '        iflow=nint(jamp2'//s//'(0))'
            write (12,*) '        jamp2cum(0)=0d0'
            write (12,*) '        do i=1,iflow'
            write (12,*) '           jamp2cum(i)=jamp2cum(i-1)+jamp2'//
     &           s//'(i)'
            write (12,*) '        enddo'
            write (12,*) '        goto 20'
      enddo
      write (12,*) '     endif'
      write (12,*) '     '
      if (need_switching) then
         write (12,*) '     do while(mtc)'
         write (12,*) '        do i='//ls//'+1,nexternal'
         write (12,*) '           ic(i)=ic(i)-'//ls
         write (12,*) '        enddo'
         write (12,*) '        goto 10'
         write (12,*) '     enddo'
         write (12,*) '     if (.not.mtc) then'
         write (12,*) '        write (*,*) "Error #1,'//
     &        ' in sborn_proc.f"'
         write (12,*) '        stop'
         write (12,*) '     endif'
         write (12,*) '     '
      endif
      write (12,*) '20   continue'
      write (12,*) '     xtarget=jamp2cum(iflow)*random()'
      write (12,*) '     ifl=1'
      write (12,*) '     do while (jamp2cum(ifl).lt.xtarget)'
      write (12,*) '        ifl=ifl+1'
      write (12,*) '     enddo'
      write (12,*) '     do i=1,2'
      write (12,*) '        do j=1,nexternal'
      if (need_switching) then
         write (12,*) '           color1(i,j)=ICOLUP(i,j,ifl)'
      else
         write (12,*) '           color(i,j)=ICOLUP(i,j,ifl)'
      endif
      write (12,*) '        enddo'
      write (12,*) '     enddo'
      if (need_switching)
     &     write (12,*) '     call switchcolor(color1,color,'
      if (need_switching) write (12,*) '    &     ic,nexternal)'
      write (12,*) '     '
      write (12,*) '     return'
      write (12,*) '     end'
      write (12,*) '     '
      write (12,*) '     '
      write (12,*) '     '
      write (12,*) '     '



      close (12)

      open (unit=15,file='init_processes.f',status='old',
     &     access='append')
      write (15,*) ''
      write (15,*) ''
      write (15,*) '     subroutine init_processes_real'
      write (15,*) '     implicit none'
      write (15,*) '     include "nlegborn.h"'
      write (15,*) '     include "pwhg_flst.h"'
      write (15,*) ''
      do i=1,total_proc
         do j=1,npart(1)
            write (15,'(6X,a10,i4,a1,i4,a2,i12)')
     &           'flst_real(',j,',',i,')=',partid(j,i)
         enddo
         write (15,*) ''
      enddo
      write (15,*) '     flst_nreal=',total_proc
      write (15,*) ''
      write (15,*) '     return'
      write (15,*) '     end'
      write (15,*) ''
      close(15)


      return
 32   write (*,*)
     &     "ERROR: file 'proc_number' not found or not correct format."
      return
 33   write (*,*)
     &     "ERROR: file 'proc_label' not found or not correct format."
      end


      subroutine convert_to_sting(npart,id,string,lstring)
      implicit none
      integer npart,lstring,i,maxpart
      parameter (maxpart=20)
      integer id(maxpart)
      character*140 string
      character*8 s
      
      ! DEBUG
      !print*, "id = ", id

      do i=1,140
         string(i:i)=' '
      enddo
      
      lstring=0
      
      
      do i=1,npart
         if (id(i).eq.21) id(i)=0
         if (abs(id(i)).le.9) then
            s=char(abs(id(i))+48)
         elseif(abs(id(i)).le.99)then
            s=char(abs(id(i))/10+48)
     &        //char(mod(abs(id(i)),10)+48)
         elseif(abs(id(i)).le.999) then
            s=char(abs(id(i))/100+48)
     &        //char((abs(id(i))-(abs(id(i))/100)*100)/10+48)
     &        //char(mod(abs(id(i))-(abs(id(i))/100)*100,10)+48)
         elseif(abs(id(i)).le.9999) then
            s=char(abs(id(i))/1000+48)
     &        //char((abs(id(i))-(abs(id(i))/1000)*1000)/100+48)
     &        //char((abs(id(i))-(abs(id(i))/100)*100)/10+48)
     &        //char(mod(abs(id(i))-(abs(id(i))/100)*100,10)+48)
         elseif(abs(id(i)).le.99999) then
            s=char(abs(id(i))/10000+48)
     &        //char((abs(id(i))-(abs(id(i))/10000)*10000)/1000+48)
     &        //char((abs(id(i))-(abs(id(i))/1000)*1000)/100+48)
     &        //char((abs(id(i))-(abs(id(i))/100)*100)/10+48)
     &        //char(mod(abs(id(i))-(abs(id(i))/100)*100,10)+48)
         elseif(abs(id(i)).le.999999) then
            s=char(abs(id(i))/100000+48)
     &        //char((abs(id(i))-(abs(id(i))/100000)*100000)/10000+48)
     &        //char((abs(id(i))-(abs(id(i))/10000)*10000)/1000+48)
     &        //char((abs(id(i))-(abs(id(i))/1000)*1000)/100+48)
     &        //char((abs(id(i))-(abs(id(i))/100)*100)/10+48)
     &        //char(mod(abs(id(i))-(abs(id(i))/100)*100,10)+48)
         elseif(abs(id(i)).le.9999999) then
            s=char(abs(id(i))/1000000+48)
     &       //char((abs(id(i))-(abs(id(i))/1000000)*1000000)/100000+48)
     &       //char((abs(id(i))-(abs(id(i))/100000)*100000)/10000+48)
     &       //char((abs(id(i))-(abs(id(i))/10000)*10000)/1000+48)
     &       //char((abs(id(i))-(abs(id(i))/1000)*1000)/100+48)
     &       //char((abs(id(i))-(abs(id(i))/100)*100)/10+48)
     &       //char(mod(abs(id(i))-(abs(id(i))/100)*100,10)+48)
         else
            write (*,*) 'error, particle ID is too large',abs(id(i))
            stop
         endif
         
         if (id(i).ge.0) then
            if (id(i).le.9) then
               string(lstring+1:lstring+1)=s
               lstring=lstring+1
            elseif (id(i).le.99) then
               string(lstring+1:lstring+2)=s
               lstring=lstring+2
            elseif (id(i).le.999) then
               string(lstring+1:lstring+3)=s
               lstring=lstring+3
            elseif (id(i).le.9999) then
               string(lstring+1:lstring+4)=s
               lstring=lstring+4
            elseif (id(i).le.99999) then
               string(lstring+1:lstring+5)=s
               lstring=lstring+5
            elseif (id(i).le.999999) then
               string(lstring+1:lstring+6)=s
               lstring=lstring+6
            elseif (id(i).le.9999999) then
               string(lstring+1:lstring+7)=s
               lstring=lstring+7
            endif
         else
            if (abs(id(i)).le.9) then
               string(lstring+1:lstring+2)='-'//s
               lstring=lstring+2
            elseif (abs(id(i)).le.99) then
               string(lstring+1:lstring+3)='-'//s
               lstring=lstring+3
            elseif (abs(id(i)).le.999) then
               string(lstring+1:lstring+4)='-'//s
               lstring=lstring+4
            elseif (abs(id(i)).le.9999) then
               string(lstring+1:lstring+5)='-'//s
               lstring=lstring+5
            elseif (abs(id(i)).le.99999) then
               string(lstring+1:lstring+6)='-'//s
               lstring=lstring+6
            elseif (abs(id(i)).le.999999) then
               string(lstring+1:lstring+7)='-'//s
               lstring=lstring+7
            elseif (abs(id(i)).le.9999999) then
               string(lstring+1:lstring+8)='-'//s
               lstring=lstring+8
            endif
         endif
      enddo
      return
      end

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
      integer date(3), time(3)
      integer flst_gosam(maxproc),k,alpha,alphas
      
      open(unit=10, file='proc_number',status='old',err=32)
      read (10,*,err=32,end=32) total_proc
      close(10)
      write (*,*) 'Number of subprocesses: ',total_proc

      if (total_proc.ge.999) then
         write (*,*) 'ERROR: Too many subprocesses, cannot continue.'
         stop
      endif

      open (unit=11,file='proc_label',status='old',err=33)
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


      open (unit=12,file='sborn_proc.f',status='unknown')

      write (12,*) '     subroutine sborn_proc(p_born,'//
     &     'legs,wgt,wgtjk,wgtmunu)'
      write (12,*) '     implicit none'
      write (12,*) '     include "nexternal.inc"'
      write (12,*) '     include "coupl.inc"'
      write (12,*) '     double precision wgt'
      write (12,*) '     double precision p_born(0:3,nexternal-1),'//
     &     'wgt2(nexternal-1),'
      write (12,*) '    &   wgtmunu(0:3,0:3,nexternal-1),'//
     &     'wgtjk(nexternal-1,nexternal-1)'
      if (need_switching) write (12,*) '     double precision '//
     &     'p_born1(0:3,nexternal-1),wgt1(0:nexternal-1),'
      if (need_switching) write (12,*)
     &     '    &   wgtmunu1(0:3,0:3,nexternal-1),'//
     &     'wgtjk1(nexternal-1,nexternal-1)'
      write (12,*) '     integer legs(nexternal-1),lstr,i'
      write (12,*) '     character*140 str'
      if (need_switching) write (12,*) '     integer ic(nexternal-1),'
     &     //'legs1(nexternal-1)'
      if (need_switching) write (12,*) '     logical mtc,even'
      write (12,*) '     logical calculatedBorn'
      write (12,*) '     integer skip'
      write (12,*) '     common/cBorn/calculatedBorn,skip'
      write (12,*) ''
      write (12,*) '     calculatedBorn=.false.'
      write (12,*) '     '
      if (need_switching) then
         write (12,*) '     do i=1,nexternal-1'
         write (12,*) '        ic(i)=i'
         write (12,*) '     enddo'
         write (12,*) '     mtc=.false.'
         write (12,*) '10   call nexper(nexternal-1-'//ls//
     &        ',ic('//ls//'+1),mtc,even)'
         write (12,*) '     do i='//ls//'+1,nexternal-1'
         write (12,*) '        ic(i)=ic(i)+'//ls
         write (12,*) '     enddo'
         write (12,*)
     &        '     CALL SWITCHMOM (P_born,P_born1,IC,NEXTERNAL-1)'
         write (12,*) '     CALL SWITCHLEGS(legs,legs1,IC,NEXTERNAL-1)'
         write (12,*) '     '
         write (12,*)
     &        '     call convert_to_string(nexternal-1,legs1,str,lstr)'
         write (12,*) '     '
      else 
         write (12,*)
     &        '     call convert_to_string(nexternal-1,legs,str,lstr)'
         write (12,*) '     '
      endif
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
            write (12,*) '        call sborn_cl_'//s//
     &           '(p_born1,wgtmunu1,wgt1)'
            write (12,*) '        call sborn_sf_'//s//
     &           '(p_born1,wgtjk1)'
         else
            write (12,*) '        call sborn_cl_'//s//
     &           '(p_born,wgtmunu,wgt2)'
            write (12,*) '        call sborn_sf_'//s//
     &           '(p_born,wgtjk)'
         endif
         write (12,*) '        goto 20'
      enddo
      write (12,*) '     endif'
      write (12,*) '     '
      if (need_switching) then
         write (12,*) '     do while(mtc)'
         write (12,*) '        do i='//ls//'+1,nexternal-1'
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
      write (12,*) '20   wgt=0d0'

      if (need_switching) write (12,*) '     call switchborns('//
     &     'wgt1(1),wgt2,wgtjk1,wgtjk,wgtmunu1,wgtmunu,'
      if (need_switching) write (12,*) '    &     ic,nexternal-1)'
      write (12,*) '     do i=1,nexternal-1'
      write (12,*) '        if(wgt.eq.0d0 .and. wgt2(i).ne.0d0) then'
      write (12,*) '           wgt=wgt2(i)'
      write (12,*)
     &     '        elseif (wgt.ne.0d0 .and. wgt2(i).ne.0d0 .and.'
      write (12,*)
     &     '    &           abs((wgt-wgt2(i))/wgt).gt.1d-7 ) then'
      write (12,*) '           write (*,*) "Error #2 in sborn_proc'//
     &     ' ",i,wgt2'
      write (12,*) '           stop'
      write (12,*) '        endif'
      write (12,*) '     enddo'
      write (12,*) '     '
      write (12,*) '     end'
      write (12,*) '     '
      write (12,*) '     '
      write (12,*) '     '
      write (12,*) '     '
      write (12,*) '     '



      write (12,*) '     subroutine born_color(legs,color)'
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
      write (12,*) '     common/to_amps_'//s//'/amp2'//s//',jamp2'//s
      enddo
      write (12,*) '     double precision jamp2cum(0:maxamps)'
      write (12,*) '     integer ICOLUP(2,nexternal-1,maxamps)'
      write (12,*) '     integer color(2,nexternal-1),'
     &     //'color1(2,nexternal-1)'
      write (12,*) '     double precision random,xtarget'
      write (12,*) '     external random'
      write (12,*) '     integer legs(nexternal-1),lstr,i,j'
      write (12,*) '     character*140 str'
      if (need_switching) write (12,*) '     integer ic(nexternal-1),'
     &     //'legs1(nexternal-1)'
      write (12,*) '     integer iflow,ifl'
      if (need_switching) write (12,*) '     logical mtc,even'
      write (12,*) '     '
      if (need_switching) then
         write (12,*) '     do i=1,nexternal-1'
         write (12,*) '        ic(i)=i'
         write (12,*) '     enddo'
         write (12,*) '     mtc=.false.'
         write (12,*) '10   call nexper(nexternal-1-'//ls//
     &        ',ic('//ls//'+1),mtc,even)'
         write (12,*) '     do i='//ls//'+1,nexternal-1'
         write (12,*) '        ic(i)=ic(i)+'//ls
         write (12,*) '     enddo'
         write (12,*) '     CALL SWITCHLEGS(legs,legs1,IC,NEXTERNAL-1)'
         write (12,*) '     '
         write (12,*)
     &        '     call convert_to_string(nexternal-1,legs1,str,lstr)'
      else
         write (12,*)
     &        '     call convert_to_string(nexternal-1,legs,str,lstr)'
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
            write (12,*) '        include "leshouches_'//s//'.inc"'
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
         write (12,*) '        do i='//ls//'+1,nexternal-1'
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
      write (12,*) '        do j=1,nexternal-1'
      if (need_switching) then
         write (12,*) '           color1(i,j)=ICOLUP(i,j,ifl)'
      else
         write (12,*) '           color(i,j)=ICOLUP(i,j,ifl)'
      endif
      write (12,*) '        enddo'
      write (12,*) '     enddo'
      if (need_switching)
     &     write (12,*) '     call switchcolor(color1,color,'
      if (need_switching) write (12,*) '    &     ic,nexternal-1)'
      write (12,*) '     '
      write (12,*) '     end'
      write (12,*) '     '
      write (12,*) '     '
      write (12,*) '     '
      write (12,*) '     '



      write (12,*) "     subroutine"
     &     //" convert_to_string(npart,id,string,lstring)"
      write (12,*) "     implicit none"
      write (12,*) "     integer npart,lstring,i"
      write (12,*) "     integer id(npart)"
      write (12,*) "     character*140 string"
      write (12,*) "     character*8 s"
      write (12,*) "     "
      write (12,*) "     do i=1,140"
      write (12,*) "        string(i:i)=' '"
      write (12,*) "     enddo"
      write (12,*) "     lstring=0"
      write (12,*) "     do i=1,npart"
      write (12,*) "        if (id(i).eq.21) id(i)=0"
      write (12,*) "        if (abs(id(i)).le.9) then"
      write (12,*) "           s=char(abs(id(i))+48)"
      write (12,*) "        elseif(abs(id(i)).le.99)then"
      write (12,*) "           s=char(abs(id(i))/10+48)"
      write (12,*) "    &           //char(mod(abs(id(i)),10)+48)"
      write (12,*) "              elseif(abs(id(i)).le.999) then"
      write (12,*) "                 s=char(abs(id(i))/100+48)"
      write (12,*) "    &           //char((abs(id(i))"
     &     //"-(abs(id(i))/100)*100)/10+48)"
      write (12,*) "    &           //char(mod(abs(id(i))"
     &     //"-(abs(id(i))/100)*100,10)+48)"
      write (12,*) "        elseif(abs(id(i)).le.9999) then"
      write (12,*) "           s=char(abs(id(i))/1000+48)"
      write (12,*) "    &        //char((abs(id(i))"
     &     //"-(abs(id(i))/1000)*1000)/100+48)"
      write (12,*) "    &        //char((abs(id(i))"
     &     //"-(abs(id(i))/100)*100)/10+48)"
      write (12,*) "    &        //char(mod(abs(id(i))"
     &     //"-(abs(id(i))/100)*100,10)+48)"
      write (12,*) "        elseif(abs(id(i)).le.99999) then"
      write (12,*) "          s=char(abs(id(i))/10000+48)"
      write (12,*) "    &        //char((abs(id(i))"
     &     //"-(abs(id(i))/10000)*10000)/1000+48)"
      write (12,*) "    &        //char((abs(id(i))"
     &     //"-(abs(id(i))/1000)*1000)/100+48)"
      write (12,*) "    &        //char((abs(id(i))"
     &     //"-(abs(id(i))/100)*100)/10+48)"
      write (12,*) "    &        //char(mod(abs(id(i))"
     &     //"-(abs(id(i))/100)*100,10)+48)"
      write (12,*) "        elseif(abs(id(i)).le.999999) then"
      write (12,*) "           s=char(abs(id(i))/100000+48)"
      write (12,*) "    &        //char((abs(id(i))"
     &     //"-(abs(id(i))/100000)*100000)/10000+48)"
      write (12,*) "    &        //char((abs(id(i))"
     &     //"-(abs(id(i))/10000)*10000)/1000+48)"
      write (12,*) "    &        //char((abs(id(i))"
     &     //"-(abs(id(i))/1000)*1000)/100+48)"
      write (12,*) "    &        //char((abs(id(i))"
     &     //"-(abs(id(i))/100)*100)/10+48)"
      write (12,*) "    &        //char(mod(abs(id(i))"
     &     //"-(abs(id(i))/100)*100,10)+48)"
      write (12,*) "        elseif(abs(id(i)).le.9999999) then"
      write (12,*) "           s=char(abs(id(i))/1000000+48)"
      write (12,*) "    &       //char((abs(id(i))"
     &     //"-(abs(id(i))/1000000)*1000000)/100000+48)"
      write (12,*) "    &       //char((abs(id(i))"
     &     //"-(abs(id(i))/100000)*100000)/10000+48)"
      write (12,*) "    &       //char((abs(id(i))"
     &     //"-(abs(id(i))/10000)*10000)/1000+48)"
      write (12,*) "    &       //char((abs(id(i))"
     &     //"-(abs(id(i))/1000)*1000)/100+48)"
      write (12,*) "    &       //char((abs(id(i))"
     &     //"-(abs(id(i))/100)*100)/10+48)"
      write (12,*) "    &       //char(mod(abs(id(i))"
     &     //"-(abs(id(i))/100)*100,10)+48)"
      write (12,*) "        else"
      write (12,*) "           write (*,*) "
     &     //"'error, particle ID is too large',abs(id(i))"
      write (12,*) "        endif"
      write (12,*) "        if (id(i).ge.0) then"
      write (12,*) "           if (id(i).le.9) then"
      write (12,*) "              string(lstring+1:lstring+1)=s"
      write (12,*) "              lstring=lstring+1"
      write (12,*) "           elseif (id(i).le.99) then"
      write (12,*) "              string(lstring+1:lstring+2)=s"
      write (12,*) "              lstring=lstring+2"
      write (12,*) "           elseif (id(i).le.999) then"
      write (12,*) "              string(lstring+1:lstring+3)=s"
      write (12,*) "              lstring=lstring+3"
      write (12,*) "           elseif (id(i).le.9999) then"
      write (12,*) "             string(lstring+1:lstring+4)=s"
      write (12,*) "             lstring=lstring+4"
      write (12,*) "           elseif (id(i).le.99999) then"
      write (12,*) "             string(lstring+1:lstring+5)=s"
      write (12,*) "             lstring=lstring+5"
      write (12,*) "           elseif (id(i).le.999999) then"
      write (12,*) "             string(lstring+1:lstring+6)=s"
      write (12,*) "             lstring=lstring+6"
      write (12,*) "           elseif (id(i).le.9999999) then"
      write (12,*) "             string(lstring+1:lstring+7)=s"
      write (12,*) "             lstring=lstring+7"
      write (12,*) "           endif"
      write (12,*) "        else"
      write (12,*) "           if (abs(id(i)).le.9) then"
      write (12,*) "              string(lstring+1:lstring+2)='-'//s"
      write (12,*) "              lstring=lstring+2"
      write (12,*) "           elseif (abs(id(i)).le.99) then"
      write (12,*) "              string(lstring+1:lstring+3)='-'//s"
      write (12,*) "              lstring=lstring+3"
      write (12,*) "           elseif (abs(id(i)).le.999) then"
      write (12,*) "              string(lstring+1:lstring+4)='-'//s"
      write (12,*) "              lstring=lstring+4"
      write (12,*) "           elseif (abs(id(i)).le.9999) then"
      write (12,*) "              string(lstring+1:lstring+5)='-'//s"
      write (12,*) "              lstring=lstring+5"
      write (12,*) "           elseif (abs(id(i)).le.99999) then"
      write (12,*) "              string(lstring+1:lstring+6)='-'//s"
      write (12,*) "              lstring=lstring+6"
      write (12,*) "           elseif (abs(id(i)).le.999999) then"
      write (12,*) "              string(lstring+1:lstring+7)='-'//s"
      write (12,*) "              lstring=lstring+7"
      write (12,*) "           elseif (abs(id(i)).le.9999999) then"
      write (12,*) "              string(lstring+1:lstring+8)='-'//s"
      write (12,*) "              lstring=lstring+8"
      write (12,*) "           endif"
      write (12,*) "        endif"
      write (12,*) "     enddo"
      write (12,*) "     end"

      close(12)

      open (unit=13,file='nexternal.inc',status='unknown')
      write (13,*) ''
      write (13,*) '     integer nexternal'
      write (13,*) '     parameter (nexternal=',npart(1)+1,')'
      write (13,*) ''
      close(13)

      open (unit=14,file='nlegborn.h',status='unknown')
      write (14,*) '     '
      write (14,*) '     integer nlegborn,nlegreal'
      write (14,*) '     parameter (nlegborn=',npart(1),')'
      write (14,*) '     parameter (nlegreal=nlegborn+1)'
      write (14,*) '     '
      write (14,*) '     integer ndiminteg'
      write (14,*) '     parameter (ndiminteg=(nlegreal-2)*3-4+2-1)'
      write (14,*) '     '
      write (14,*) '     integer maxprocborn,maxprocreal'
      write (14,*) '     parameter (maxprocborn=999,maxprocreal=999)'
      write (14,*) '     '
      write (14,*) '     parameter '//
     $     '(maxalr=maxprocreal*nlegreal*(nlegreal-1)/2)'
      close(14)

      open (unit=15,file='init_processes.f',status='unknown')
      write (15,*) '     subroutine init_processes'
      write (15,*) '     implicit none'
      write (15,*) '     include "nlegborn.h"'
      write (15,*) '     include "pwhg_flst.h"'
      write (15,*) '     include "pwhg_st.h"'
      write (15,*) '     include "coupl.inc"'
      write (15,*) '     integer i'
      write (15,*) ''
      write (15,*) '     call init_processes_born'
      write (15,*) '     call init_processes_real'
      write (15,1) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
      write (15,1) 'C    Set here the number of light flavours'
      write (15,1) 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'
      write (15,*) '     st_nlight=5'
      write (15,*) '     call init_couplings'
      write (15,1) 'c     if (tmass.eq.0d0) then'
      write (15,1) 'c        st_nlight=6'
      write (15,1) 'c     elseif(bmass.eq.0d0) then'
      write (15,1) 'c        st_nlight=5'
      write (15,1) 'c     elseif(cmass.eq.0d0) then'
      write (15,1) 'c        st_nlight=4'
      write (15,1) 'c     else'
      write (15,1) 'c        st_nlight=3'
      write (15,1) 'c     endif'
      write (15,*) '     do i=3,nlegreal'
      write (15,*) '        if (abs(flst_real(i,1)).le.st_nlight) then'
      write (15,*) '           flst_lightpart=i'
      write (15,*) '           exit'
      write (15,*) '        endif'
      write (15,*) '     enddo'
      write (15,*) ''
      write (15,*) '     end'
      write (15,*) ''
      write (15,*) ''
      write (15,*) ''
      write (15,*) '     subroutine init_processes_born'
      write (15,*) '     implicit none'
      write (15,*) '     include "nlegborn.h"'
      write (15,*) '     include "pwhg_flst.h"'
      write (15,*) ''
      do i=1,total_proc
         do j=1,npart(1)
            write (15,'(6X,a10,i4,a1,i4,a2,i12)')
     &           'flst_born(',j,',',i,')=',partid(j,i)
         enddo
         write (15,*) ''
      enddo
      write (15,*) '     flst_nborn=',total_proc
      write (15,*) ''
      write (15,*) '     end'
      write (15,*) ''
 1    format (a)
      close(15)


      call idate(date)
      call itime(time)
      open(unit=16, file='proc_couplings',status='old',err=34)
      read (16,*,err=32,end=32) alpha
      read (16,*,err=32,end=32) alphas
      close(16)

      open(unit=18,file='orderfile.lh',err=35)
      write(*,*) '************************************************'//
     $     '******'
      write(*,*) 'Writing the orderfile...'
      
      write(18,*) '# orderfile.lh'
      write(18,*) '# Created by POWHEG-BOX'
      write(18,1000) date(1), date(2), date(3), time 
      write(18,*)
      write(18,*) 'MatrixElementSquareType CHaveraged'
      write(18,*) 'CorrectionType          QCD'
      write(18,*) 'IRregularisation        CDR'
      write(18,'(1x,a,i2)') 'AlphasPower            ',alphas-1
      write(18,'(1x,a,i2)') 'AlphaPower             ',alpha
      write(18,*) 'SubdivideSubprocess     no'
      write(18,*)
      write(18,*) '# processes  list'
      do j=1,total_proc
         do k=1,npart(1)
            flst_gosam(k)=partid(k,j)
            if (flst_gosam(k).eq.0) flst_gosam(k)=21
         enddo
         write(18,'(2i4,a3,10i4)') (flst_gosam(k),k=1,2),
     $  ' ->', (flst_gosam(k),k=3,npart(1))
      enddo
      close(18)

 1000 format ( ' # Date ', i2.2, '/', i2.2, '/', i4.4, '; time ',
     &     i2.2, ':', i2.2, ':', i2.2 )

      return
 32   write (*,*)
     &     "ERROR: file 'proc_number' not found or not correct format."
      return
 33   write (*,*)
     &     "ERROR: file 'proc_label' not found or not correct format."
 34   write (*,*)
     &  "ERROR: file 'proc_couplings' not found or not correct format."
      return
 35   write(*,*) '******************************************'
      write(*,*) 'Problem in writing the order file'
      write(*,*) 'CANNOT PROCEED. POWHEG execution abort'
      write(*,*) '******************************************'
      return
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
      
      end


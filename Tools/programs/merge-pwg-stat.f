c program to merge gnuplot data tables.
c a line starting with '#' followed by a line with 4 numbers
c is considered the beginning of a data set. All files to be
c merged must be identical in format.
c When the program starts, it expects as input a list of filenames,
c one per line, to be merged. An empty lines terminates the list.
      program merge_pwg_stat
      implicit none
      integer maxfiles,maxlines
      parameter (maxfiles=1000,maxlines=15000)
      character *(100) files(maxfiles)
      character *(300) line(maxlines,maxfiles)
      character *(300) templine
      character *(300) searchstr
      integer nlines(maxfiles)
      integer ifile,nfiles,ios,k,imethod,N_NaN
      character cmethod
      real * 8 v1,v2,v3,v4,y,err
      integer ilength, ipos
      external ilength
      integer numberlength
      parameter (numberlength=22)
      double precision temp
      double precision btildepos(maxfiles), btildeposerr(maxfiles)
      double precision btildeneg(maxfiles), btildenegerr(maxfiles)
      double precision btildetot(maxfiles), btildetoterr(maxfiles)
      double precision realospos(maxfiles), realosposerr(maxfiles)
      double precision realosneg(maxfiles), realosnegerr(maxfiles)
      double precision realostot(maxfiles), realostoterr(maxfiles)
      double precision tot(maxfiles), toterr(maxfiles)

      double precision cbtildepos, cbtildeposerr
      double precision cbtildeneg, cbtildenegerr
      double precision cbtildetot, cbtildetoterr
      double precision crealospos, crealosposerr
      double precision crealosneg, crealosnegerr
      double precision crealostot, crealostoterr
      double precision ctot, ctoterr


      ! init
      cbtildepos = 0D0; cbtildeposerr = 0D0
      cbtildeneg = 0D0; cbtildenegerr = 0D0
      cbtildetot = 0D0; cbtildetoterr = 0D0
      crealospos = 0D0; crealosposerr = 0D0
      crealosneg = 0D0; crealosnegerr = 0D0
      crealostot = 0D0; crealostoterr = 0D0
      ctot = 0D0; ctoterr = 0D0

      do ifile=1,maxfiles
         CALL getarg(ifile,files(ifile))
         if(trim(files(ifile)).eq.'') then
            nfiles=ifile-1
            write(6,*) 'mergedata found',nfiles,
     $                 'files on the command line ...'
            goto 9
         endif
      enddo
 9    continue


      ! load data
      do ifile=1,nfiles
         open(unit=11,file=files(ifile),status='old')
         do k=1,maxlines+1
            read(unit=11,fmt='(a)',end=111) line(k,ifile)
            if(k.eq.maxlines+1) then
               write(*,*) ' too many lines in file, increase maxlines'
               call exit(-1)
            endif
            goto 12
 111        nlines(ifile)=k-1
            goto 11
 12         continue
         enddo
 11      continue
      enddo

      do ifile=1,nfiles
         if(nlines(ifile).ne.nlines(1)) then
            write(*,*) ' error: file', files(ifile),
     &           ' does not match in length'
            call exit(-1)
         endif
      enddo

      do ifile=1,nfiles
        searchstr = "btilde pos.   weights:"
        call getnumber(nlines(1),line(:,ifile),searchstr,
     &                 btildepos(ifile),btildeposerr(ifile))
        searchstr = "btilde |neg.| weights:"
        call getnumber(nlines(1),line(:,ifile),searchstr,
     &                 btildeneg(ifile),btildenegerr(ifile))
        searchstr = "btilde Total (pos.-|neg.|):"
        call getnumber(nlines(1),line(:,ifile),searchstr,
     &                 btildetot(ifile),btildetoterr(ifile))
        searchstr = "real_osres pos.   weights:"
        call getnumber(nlines(1),line(:,ifile),searchstr,
     &                 realospos(ifile),realosposerr(ifile))
        searchstr = "real_osres |neg.| weights:"
        call getnumber(nlines(1),line(:,ifile),searchstr,
     &                 realosneg(ifile),realosnegerr(ifile))
        searchstr = "real_osres total (pos.-|neg.|):"
        call getnumber(nlines(1),line(:,ifile),searchstr,
     &                 realostot(ifile),realostoterr(ifile))
        searchstr = "total (btilde+remnants+regulars+osresR) cross "//
     &              "section in pb"
        call getnumber(nlines(1),line(:,ifile),searchstr,
     &                 tot(ifile),toterr(ifile))
      enddo

      !do ifile=1,2
      !print*,btildepos(ifile),btildeposerr(ifile)
      !print*,btildeneg(ifile),btildenegerr(ifile)
      !print*,btildetot(ifile),btildetoterr(ifile)
      !print*,realospos(ifile),realosposerr(ifile)
      !print*,realosneg(ifile),realosnegerr(ifile)
      !print*,realostot(ifile),realostoterr(ifile)
      !print*,tot(ifile),toterr(ifile)
      !enddo

      do ifile=1,nfiles
        cbtildepos = cbtildepos + btildepos(ifile)
        cbtildeposerr = cbtildeposerr + btildeposerr(ifile)**2

        cbtildeneg = cbtildeneg + btildeneg(ifile)
        cbtildenegerr = cbtildenegerr + btildenegerr(ifile)**2

        cbtildetot = cbtildetot + btildetot(ifile)
        cbtildetoterr = cbtildetoterr + btildetoterr(ifile)**2

        crealospos = crealospos + realospos(ifile)
        crealosposerr = crealosposerr + realosposerr(ifile)**2

        crealosneg = crealosneg + realosneg(ifile)
        crealosnegerr = crealosnegerr + realosnegerr(ifile)**2

        crealostot = crealostot + realostot(ifile)
        crealostoterr = crealostoterr + realostoterr(ifile)**2

        ctot = ctot + tot(ifile)
        ctoterr = ctoterr + toterr(ifile)**2
      enddo

      cbtildepos = cbtildepos/nfiles
      cbtildeneg = cbtildeneg/nfiles
      cbtildetot = cbtildetot/nfiles

      crealospos = crealospos/nfiles
      crealosneg = crealosneg/nfiles
      crealostot = crealostot/nfiles

      cbtildeposerr = dsqrt(cbtildeposerr/nfiles**2)
      cbtildenegerr = dsqrt(cbtildenegerr/nfiles**2)
      cbtildetoterr = dsqrt(cbtildetoterr/nfiles**2)

      crealosposerr = dsqrt(crealosposerr/nfiles**2)
      crealosnegerr = dsqrt(crealosnegerr/nfiles**2)
      crealostoterr = dsqrt(crealostoterr/nfiles**2)

      ctot = ctot/nfiles
      ctoterr = dsqrt(ctoterr/nfiles**2)

      print*,"btilde pos.   weights: ",cbtildepos," +- ",cbtildeposerr
      print*,"btilde |neg.| weights: ",cbtildeneg," +- ",cbtildenegerr
      print*,"btilde Total (pos.-|neg.|): ",cbtildetot," +- ",
     &       cbtildetoterr
      print*,"real_osres pos.   weights: ",crealospos," +- ",
     &       crealosposerr
      print*,"real_osres |neg.| weights: ",crealosneg," +- ",
     &       crealosnegerr
      print*,"real_osres total (pos.-|neg.|): ",crealostot," +- ",
     &       crealostoterr
      print*,"total (btilde+remnants+regulars+osresR) cross section "//
     &       "in pb ",ctot," +- ", ctoterr

      end


      function ilength(line)
      integer ilength
      character *(*) line
      ilength=len(line)
      do j=ilength,1,-1
         if(line(j:j).ne.' ') then
            ilength=j
            return
         endif
      enddo
      ilength=0
      end

      subroutine getnumber(nlines,line,searchstr,dnumber,derror)
      implicit none
      integer maxlines
      parameter (maxlines=15000)
      character *(300) line(maxlines)
      character *(300) templine
      character *(300) searchstr
      integer strl
      integer nlines
      integer ilength, ipos
      external ilength
      integer numberlength
      parameter (numberlength=22)
      double precision temp, dnumber, derror
      integer k
      dnumber = 0D0
      derror = 0D0
      strl = len(trim(searchstr))+3 ! accounts for the 3 spaces
      do k=1,nlines
        ipos = index(line(k),trim(searchstr))
        if(ipos .ne. 0) then
          templine = line(k)
          read(templine(ipos+strl:ipos+strl+numberlength),*) dnumber
          ipos = index(templine,"+-   ")
          if(ipos .ne. 0) then
            read(templine(ipos+5:ipos+5+numberlength),*) derror
          endif
        endif
      enddo
      end

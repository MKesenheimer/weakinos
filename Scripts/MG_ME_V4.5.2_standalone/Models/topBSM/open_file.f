
           subroutine open_file_mdl(lun,filename,fopened)
c***********************************************************************
c     opens file input-card.dat in current directory or above
c***********************************************************************
      implicit none
c
c     Arguments
c
      integer lun
      logical fopened
      character*(*) filename
      character*90  tempname
      integer fine
      integer dirup,i

c-----
c  Begin Code
c-----
c
c     first check that we will end in the main directory
c
      tempname=filename
      fine=index(tempname,' ')
      if(fine.eq.0) fine=len(tempname)
      tempname=tempname(1:fine)
c
c         if I have to read a card
c
      if(index(filename,"_card").gt.0) then
         tempname='./Cards/'//tempname
      endif


      fopened=.false.
      do i=0,5
         open(unit=lun,file=tempname,status='old',ERR=30)
         fopened=.true.
         write(*,*) 'read model file',tempname
         exit
 30      tempname='../'//tempname
         if (i.eq.5)then
            write(*,*) 'Warning: file ',tempname,' is not correct'
            stop
         endif
      enddo


      return
      end
 

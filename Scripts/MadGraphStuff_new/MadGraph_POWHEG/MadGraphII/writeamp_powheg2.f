!-----------------------
! madgraph - a Feynman Diagram package by Tim Stelzer and Bill Long 
! (c) 1993
!
! Filename: writesub.f
!-----------------------
      subroutine write_symmetry
!***************************************************************************
c
c     Here we'll write out the identical particle information
c
!***************************************************************************
      implicit none

! Constants

      include 'params.inc'

! Arguments


! Local Variables

      integer i,j,ident
      integer jmatch(0:maxlines,0:maxlines)
      logical done

! External
      
! Global Variables
      integer        iline(-maxlines:maxlines),idir(-maxlines:maxlines)
      integer        this_coup(max_coup) ,goal_coup(max_coup)
      common/to_proc/iline,idir,this_coup,goal_coup

!-----------
! Begin Code
!-----------

      jmatch(0,0) = 0
      ident = 0
      done = .false.
      do i=3,iline(0)
         j = 1
         done = .false.
         do while (j .le. jmatch(0,0) .and. .not. done)
            if (jmatch(0,j) .eq. iline(i)) then
               done=.true.
            else
               j=j+1
            endif
         enddo
         if (j .gt. jmatch(0,0)) then
            jmatch(0,0)=j
            jmatch(0,j)=iline(i)
            jmatch(j,0)=1
            jmatch(1,j)=i
         else
            jmatch(j,0)=jmatch(j,0)+1
            if (jmatch(j,0) .eq. 2) ident=ident+1  !Number of identical groups
            jmatch(jmatch(j,0),j)=i
         endif
      enddo
      end
 
      subroutine writeamp(ngraphs,nexternal,nqed,nqcd,nproc,ncolor)
!***************************************************************************
!     writes out the subroutines
!***************************************************************************
      implicit none

! Constants

      include 'params.inc'
      integer    ln,ln2   
      parameter (ln=92,ln2=58 )

! Arguments

      integer ngraphs,nqed,nqcd,nexternal,nproc,ncolor


! Local Variables

      character*75 buff
      character*20 name
c      character*4 temp
      integer i,jchar,ishift,nlength,j,k,ic,nflows
      integer idenom,inum,emit
      logical injamp

! External
      
      logical equal
      integer igcf

! Global Variables
      double precision               den(maxgraphs)
      integer           neigen,nrows,jplace(maxgraphs)
      common /to_colormat/ den,neigen,nrows,jplace

      integer flows(0:2*maxlines,0:maxfactors,0:24,0:maxflows)
      integer         graphcolor(0:2,0:maxflows,0:maxgraphs)
      common/to_color/graphcolor,flows

      integer alen(maxgraphs)
      integer wnum(-maxlines:maxlines),wcount,wlen(maxgraphs*maxlines)
      common/toopt/wcount,wnum,wlen,alen

      character*25 gname
      integer iname
      common/to_name/iname,gname
      character*60 proc
      integer iproc
      common/to_write/iproc,proc
      integer            symfact(maxgraphs)
      common /tosymmetry/symfact

      character*70       cversion
      common/to_version/ cversion

      logical              StandAlone
      common/to_StandAlone/StandAlone

      integer info_ijk(maxdipoles,2,5),nr,nrs,nrc,dipol
      logical offdiag
      common /ijk/ info_ijk,nr,nrs,nrc,offdiag
      common /dipolenumber/ dipol

      integer ncomb
      common /ihel/ ncomb

      integer         nincoming,nincfirst
      common/to_proc2/nincoming,nincfirst

      integer        iline(-maxlines:maxlines),idir(-maxlines:maxlines)
      integer        this_coup(max_coup) ,goal_coup(max_coup)
      common/to_proc/iline,idir,this_coup,goal_coup

      character*(max_string) iwave(max_particles),owave(max_particles)
      integer iposx(3,max_particles)
      integer info_p(5,max_particles)
      character*(8) str(3,max_particles)
      common/to_external/iwave,owave,iposx,info_p,str

      character*(8) mass(max_particles)
      common /masses/ mass

      character*3 cpn
      common/ccpn/cpn

!-----------
! Begin Code
!-----------
      if (nproc .gt. 0) then
         open(unit=ln,file=gname(1:iname)//'.f',status='unknown')
      else
         open(unit=ln,status='unknown')
      endif
      name = gname(1:iname)
      nlength=iname

      call upper_case(name)

      print*,'Writing function ',name(1:nlength) ,
     &     ' in file ',GNAME(1:iname)//'.f.'

      call writesum(nproc,ngraphs,ncolor)
      buff = 'REAL*8 FUNCTION '
      jchar = 17
      ishift = nlength
      buff(jchar:jchar+ishift) = name(1:nlength)
      jchar = jchar+ishift
      ishift= 15
      buff(jchar:jchar+ishift) = '(P,NHEL,HELL,IC)'
      jchar =jchar+ishift
      write(ln,'(6x,a)') buff(1:jchar)

! write comments

      write(ln,5) ' '
      write(ln,5) cversion
      write(ln,5) 'RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS'
      write(ln,5) 'FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)'
c      write(ln,5) 'FOR THE POINT IN PHASE SPACE P(0:3,N_PARTICLES)'
      write(ln,5) ' '
      write(ln,5) 'FOR PROCESS '//proc(1:iproc)
      write(ln,5) ' '

      write(ln,10) 'IMPLICIT NONE'

! write parameters

      write(ln,5) ' '
      write(ln,5) 'CONSTANTS'
      write(ln,5) ' '
      buff = 'INTEGER    NGRAPHS,    NEIGEN'
      write(ln,10) buff(1:30)
      buff = 'PARAMETER (NGRAPHS=????,NEIGEN=???)'
      write(buff(20:23),'(i4)') ngraphs
      write(buff(32:34),'(i3)') neigen
      write(ln,10) buff(1:36)
      write(ln,10) 'include "nexternal.inc"'
      buff = 'INTEGER    NWAVEFUNCS     , NCOLOR'
      write(ln,10) buff(1:34)
      buff = 'PARAMETER (NWAVEFUNCS=????, NCOLOR=????) '
      write(buff(23:26),'(i4)') WCOUNT
      write(buff(36:39),'(i4)') NCOLOR
      write(ln,10) buff(1:41)
      write(ln,10) 'REAL*8     ZERO'
      write(ln,10) 'PARAMETER (ZERO=0D0)'

! write argument declarations

      write(ln,5) ' '
      write(ln,5) 'ARGUMENTS '
      write(ln,5) ' '

      write(ln,10) 'REAL*8 P(0:3,NEXTERNAL)'
      write(ln,10) 'INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL), HELL'

! write local declarations

      write(ln,5) ' '
      write(ln,5) 'LOCAL VARIABLES '
      write(ln,5) ' '

      write(ln,10) 'INTEGER I,J'
      write(ln,10) 'COMPLEX*16 ZTEMP'
      write(ln,10) 'REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)'
      write(ln,10) 'COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)'

c by RF Jan. 2006
      write(ln,10) 'COMPLEX*16 W(18,NWAVEFUNCS)'   !,WX1(6),WX2(6)'
cend RF

! write global declarations

      write(ln,5) ' '
      write(ln,5) 'GLOBAL VARIABLES'
      write(ln,5) ' '
      write(92,10) 'integer maxamps'
      write(92,10) 'parameter (maxamps=6000)'
      write(ln,10) 'Double Precision amp2(maxamps), jamp2(0:maxamps)'
      write(ln,10) 'common/to_Ramps_'//cpn//'/  amp2,       jamp2'

      write(92,10) 'integer max_bhel'
      write(92,*) '     parameter ( max_bhel =', ncomb,')'
	  
      write(ln,10) 'include "coupl.inc"'
      write(ln,5) ' '
      write(ln,5) 'COLOR DATA'
      write(ln,5) ' '
      do while (.true.)
         read(99,'(a)',end=99) buff
         write(ln,'(a)') buff
      enddo
 99   continue

      write(ln,5) '----------'
      write(ln,5) 'BEGIN CODE'
      write(ln,5) '----------'

      do while (.true.)
         read(91,'(a)',end=999) buff
         write(ln,'(a)') buff
      enddo
 999  continue

      nflows=flows(0,0,0,0)

c      write(*,*) 'Writing flows',nflows
      do i=1,nflows
c
c     First determine a common denominator for these amplitudes
c
         idenom=1
         inum = 999999
         do j=1,ngraphs
            do k=1,graphcolor(0,0,j)         !Number of terms/flows
               if (graphcolor(0,k,j) .eq. i) then
                  inum = min(inum,abs(graphcolor(1,k,j)))
c                  write(*,*) i,j,k,idenom,graphcolor(2,k,j)
                  idenom=abs(idenom*graphcolor(2,k,j)/
     &                 igcf(idenom,graphcolor(2,k,j)))
c                  write(*,*) i,j,k,idenom,graphcolor(2,k,j)
               endif
            enddo
         enddo
         if (idenom .ne. 1) inum=1
         if (inum .eq. 1 .and. idenom .eq. 1) then
            write(ln,'(a,i4,a$)')    '      JAMP(',i,') = '
            ic = 19
         elseif (inum .ne. 1) then
            write(ln,'(a,i4,a,i6,a$)')    '      JAMP(',i,') = ',
     &           inum,'*( '
            ic = 28 
         elseif (idenom .ne. 1) then
            write(ln,'(a,i4,a$)')    '      JAMP(',i,') = ('
            ic = 20 
         endif
         do j=1,ngraphs
            do k=1,graphcolor(0,0,j)         !Number of terms/flows
            if (graphcolor(0,k,j) .eq. i) then                 
               if (ic .gt. 60) then                    
                  write(ln,*)                          
                  write(ln,'(a$)')'     &             '
                  ic = 19
               endif
               den(j)=1d0
               den(j)=dble(graphcolor(1,k,j))/graphcolor(2,k,j)*
     &              idenom*symfact(j)/inum
c               write(*,*) i,j,k,idenom,graphcolor(1,k,j),
c     &              graphcolor(2,k,j),den(j)
c               if (equal(den(j)/den(i),1d0)) then
               if (equal(den(j),1d0)) then
                  write(ln,'(a,i4,a$)') '+AMP(',j,')'
                  ic = ic+10
c               elseif (equal(den(j)/den(i),-1d0)) then
               elseif (equal(den(j),-1d0)) then
                  ic = ic+10
                  write(ln,'(a,i4,a$)') '-AMP(',j,')'
c               elseif (den(j)/den(i) .gt. 0d0) then
               elseif (den(j) .gt. 0d0) then
                  if (ic .gt. 50) then
                     write(ln,*)
                     write(ln,'(a$)')'     & '
                     ic = 8
                  endif
                  ic = ic+22
                  write(ln,'(a,e10.3,a,i4,a$)') '+'
     $                ,den(j),'*AMP(',j,')'
c               elseif (den(j)/den(i) .lt. 0d0) then
               elseif (den(j) .lt. 0d0) then
                  if (ic .gt. 50) then
                     write(ln,*)
                     write(ln,'(a$)')'     & '
                     ic = 8
                  endif
                  ic = ic+22
                  write(ln,'(a,e10.3,a,i4,a$)') '-'
     $                ,abs(den(j)),'*AMP(',j,')'
c     $                 ,den(j)/den(i),'*AMP(',j,')'
               endif
            endif
         enddo
         enddo
         if (idenom .ne. 1) then
            write(ln,'(a,i6$)') ')/',idenom 
         elseif (inum .ne. 1) then
            write(ln,'(a$)') ')'
         endif
         write(ln,*)
      enddo



      buff(1:nlength) = name(1:nlength)
      buff(nlength+1:nlength+8) = ' = 0.D0'
      write(ln,10) buff(1:nlength+8)
      write(ln,10) 'DO I = 1, NCOLOR'
      write(ln,10) '    ZTEMP = (0.D0,0.D0)'
      write(ln,10) '    DO J = 1, NCOLOR'
      write(ln,10) '        ZTEMP = ZTEMP + CF(J,I)*JAMP(J)'
      write(ln,10) '    ENDDO'
      buff(1:4) = '    '
      buff(5:nlength+5) = name(1:nlength)
      buff(nlength+6:nlength+6) = '='
      buff(nlength+7:2*nlength+7) = name(1:nlength)
      buff(2*nlength+7:2*nlength+40)=
     .     '+ZTEMP*DCONJG(JAMP(I))/DENOM(I)'
      write(ln,10) buff(1:2*nlength+40)
      write(ln,10) 'ENDDO'

      write(ln,10) 'Do I = 1, NGRAPHS'
      write(ln,10) '    amp2(i)=amp2(i)+amp(i)*dconjg(amp(i))'
      write(ln,10) 'Enddo'
      write(ln,10) 'Do I = 1, NCOLOR'
      write(ln,10) '    Jamp2(i)=Jamp2(i)+Jamp(i)*dconjg(Jamp(i))'
      write(ln,10) 'Enddo'

      write(ln,10) 'END'
      write(ln,10) ' '
      write(ln,10) ' '
            

      close(ln)



c     write the coloramps_???.inc (unit=58 is opened in mainfks.f)
      write(ln2,'(a,i4,a,i4,a)')'      LOGICAL ICOLAMP(',NGRAPHS,
     $     ',',NCOLOR,')'
      write(ln2,'(a$)')      '      DATA ICOLAMP/'
      ic = 18
      do i=1,nflows
         do j=1,ngraphs
            if (ic .gt. 60) then                    
               write(ln2,*)                          
               write(ln2,'(a$)')'     &             '
               ic = 18
            endif
            injamp=.false.
            do k=1,graphcolor(0,0,j) !Number of terms/flows
               if (graphcolor(0,k,j) .eq. i) then                 
                  injamp=.true.
               endif
            enddo
            if(injamp)then
               write(ln2,'(a$)') '.true. '
               ic = ic+8
            else
               write(ln2,'(a$)') '.false.'
               ic = ic+8
            endif
            if(i.lt.nflows.or.j.lt.ngraphs)then
               write(ln2,'(a$)') ','
            else
               write(ln2,'(a$)') '/'
            endif
         enddo
      enddo
      close(ln2)

 5    format('C ',A)
 10   format(6X,A)
 15   format('C',6X,A)
 20   format(', P',I1)
 30   format(',P',I1,'(0:3)')
 41   format(5X,'&',4X,A8,A,A8,A,A8,A)
 42   format(6X,A18,12(I2,A))
 43   format(5X,'&',4X,A8,A,A8,A,A8,A,I2,A)
       end

 
      subroutine writesum(nproc,ngraphs,ncolor)
!***************************************************************************
!     writes out the subroutine which sums over helicities and
!     and does all crossings
!***************************************************************************
      implicit none

! Constants

      include 'params.inc'

! Arguments

      integer nproc, ngraphs, ncolor

! Local Variables

      character*75 buff,fname
      character*130 buff1
      character*20 name
      character*4 str1
      integer i,jchar,ishift,j,flen,nexternal,ichar,nlen,k,l,nlength
      integer isym,jhel,initc,iden(maxcross),emit
      integer nhel(maxlines,3**maxlines),step(maxlines),ihel
      integer ncross,ipart(maxlines,maxcross),ic(maxlines,maxcross)
      integer jline(0:4),imatch(0:maxcross)
      integer icross,iboson
      logical found
      integer nwp,nwm,nz,ngot,totalhelicities
      character*3 number

! Global Variables

      integer        iline(-maxlines:maxlines),idir(-maxlines:maxlines)
      integer        this_coup(max_coup) ,goal_coup(max_coup)
      common/to_proc/iline,idir,this_coup,goal_coup

      character*25 gname
      integer iname
      common/to_name/iname,gname
      integer         nincoming,nincfirst
      common/to_proc2/nincoming,nincfirst
      character*60 proc
      integer iproc
      common/to_write/iproc,proc

      integer info_p(5,max_particles),iposx(3,max_particles)
      character*(max_string) iwave(max_particles),owave(max_particles)
      character*(8) str(3,max_particles)
      common/to_external/iwave,owave,iposx,info_p,str

      character*(4*max_particles) particle(4)
      integer                               charge_c(max_particles)
      integer iparticle(0:max_particles,0:4),inverse(max_particles)
      common/to_model/iparticle,  particle,  inverse, charge_c

      logical         lwp,lwm,lz,decay,cross
      common/to_decay/lwp,lwm,lz,decay,cross

      character*70       cversion
      common/to_version/ cversion

      logical              StandAlone
      common/to_StandAlone/StandAlone

      data StandAlone/.false./

      integer info_ijk(maxdipoles,2,5),nr,nrs,nrc,dipol
      logical offdiag
      common /ijk/ info_ijk,nr,nrs,nrc,offdiag
      common /dipolenumber/ dipol

      integer ncomb
      common /ihel/ ncomb

      integer isymn1
      common /isym/ isymn1

      character*(8) mass(max_particles)
      common /masses/ mass

      character*(max_string) dipprocess(maxdipoles)
      common /process/ dipprocess

      character*(4) part_type(maxdipoles,2)
      common /part/ part_type
      character*3 cpn
      common/ccpn/cpn

!-----------
! Begin Code
!-----------

      offdiag=.false.
      standalone=.true.
      if (cross .and. .false.) then
         call crossing(ncross,ipart,imatch)
      else
         ncross=1
         imatch(0) =1
         imatch(1) =1
         do i=1,iline(0)
            ipart(i,1)=i
         enddo
      endif
c
c set up all the different crossings which are not identical
c
      icross = 0
      do i=1,ncross
         if (imatch(i) .gt. icross) then
            call getdenom(isym,initc,jhel,step,ipart(1,i))
            iden(imatch(i)) = isym*initc*jhel
            if (imatch(i) .ne. icross+1) then
               print*,'Error in imatch writesub.f',i,imatch(i),icross+1
            endif
            icross=imatch(i)
         endif
      enddo
      if (icross .ne. imatch(0)) then
         print*,'Warning icross .ne. imatch(0)',icross,imatch(0)
      endif
      name = gname(1:iname)
      nlength=iname

      call upper_case(name)
      nexternal = iline(0)

      call sethel(nhel,step,nexternal,ihel)
      ncomb=ihel
      buff = 'SUBROUTINE S'
      jchar = 13
      ishift = iname
      buff(jchar:jchar+ishift) = name(1:nlength)
      jchar = jchar+ishift
      ishift= 7
      buff(jchar:jchar+ishift) = '(P1,ANS)'
      jchar =jchar+ishift
      write(92,10) buff(1:jchar)

      call genparton(ncross,ipart,buff(12:jchar),nproc,imatch)

      jchar=jchar-4
      FNAME= buff(13:jchar)
      FNAME(JCHAR-12:JCHAR-12+26)=',NHEL(1,IHEL),IHEL,JC(1))'
      FLEN = jchar+26

! write comments

      write(92,5) ' '
      write(92,5) cversion
	  if(standalone) write (92,5) "MadGraph StandAlone Version"
      write(92,5) 'RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS'
      write(92,5) 'AND HELICITIES'
      write(92,5) 'FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)'
      write(92,5) ' '
      write(92,5) 'FOR PROCESS '//proc(1:iproc)
      write(92,5) ' '
      write(92,10) 'IMPLICIT NONE'

! write parameters

      write(92,5) ' '
      write(92,5) 'CONSTANTS'
      write(92,5) ' '
      write(92,10) 'Include "nexternal.inc"'
      buff = 'INTEGER                 NCOMB,     NCROSS'
      write(92,10) buff(1:50)
      buff = 'PARAMETER (             NCOMB=????, NCROSS=???)'
      write(buff(31:34),'(i4)') ihel
      ncomb=ihel
      write(buff(44:46),'(i3)') icross
      write(92,10) buff(1:47)
      write(92,10) 'INTEGER    THEL'
      write(92,10) 'PARAMETER (THEL=NCOMB*NCROSS)'
      write(92,10) 'INTEGER NGRAPHS'
      buff = 'PARAMETER (NGRAPHS=????)'
      write(buff(20:23),'(i4)') ngraphs
      write(92,10) buff(1:24)

      write(92,5) ' '
      write(92,5) 'ARGUMENTS '
      write(92,5) ' '
      write(92,10) 'REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)'

      write(92,5) ' '
      write(92,5) 'LOCAL VARIABLES '
      write(92,5) ' '
      write(92,10) 'REAL*8 P(0:3,NEXTERNAL)'
      write(92,10) 'INTEGER NHEL(NEXTERNAL,NCOMB),NTRY'
      write(92,10) 'REAL*8 T'
      buff = 'REAL*8 '
      buff(8:7+nlength)=name(1:nlength)
      write(92,10) buff(1:7+nlength)

      write(92,10) 'REAL*8 ZERO'
      write(92,10) 'PARAMETER(ZERO=0d0)'
      write(92,10) 'INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)'
      write(92,10) 'INTEGER IPROC,JC(NEXTERNAL), I,L,K'
      write(92,10) 'LOGICAL GOODHEL(NCOMB,NCROSS)'

      if(standalone)  write(92,10) 'DATA NTRY/0/'
	  	  	     
      write (92,10) 'INTEGER NGOOD,igood(ncomb),jhel'
      write (92,10) 'data ngood /0/'
      write (92,10) 'save igood,jhel'
      write (92,10) 'REAL*8 hwgt'
      
      write(92,10) 'integer maxamps'
      write(92,10) 'parameter (maxamps=6000)'
      write(92,10) 'Double Precision amp2(maxamps), jamp2(0:maxamps)'
      write(92,10) 'common/to_Ramps_'//cpn//'/  amp2,       jamp2'
      write(92,10)
      if (standalone) write(92,10)'integer j,jj'
      write(92,10) 'integer max_bhel'
      write(92,*) '     parameter ( max_bhel =', ihel,')'
      write(92,10)

      write(92,10) 'INTEGER NCOLOR'
      buff = 'DATA NCOLOR   /'
      write(buff(16:24),'(i4,a)') ncolor,'/'
      write(92,10) buff(1:30)

      write(92,10) 'DATA GOODHEL/THEL*.FALSE./'
      BUFF = 'DATA (NHEL(IHEL,????),IHEL=1,??) /'
      write(buff(30:31),'(i2)') nexternal
      do i=1,ihel
         write(buff(17:20),'(i4)') i
         jchar=35
         ishift = 3
         do j=1,nexternal
            write(buff(jchar:jchar+ishift),'(I2,a)') nhel(j,i),','
            jchar=jchar+ishift
         enddo
         jchar=jchar-1
         write(buff(jchar:jchar),'(a1)') '/'
         write(92,10) buff(1:jchar)
      enddo
      BUFF = 'DATA (  IC(IHEL,???),IHEL=1,??) /'
      write(buff(29:30),'(i2)') nexternal
      icross=0
      do i=1,ncross
         if (imatch(i) .gt. icross) then
            write(buff(17:19),'(i3)') imatch(i)
            jchar=34
            ishift = 3
            do j=1,nexternal
            write(buff(jchar:jchar+ishift),'(I2,a)') ipart(j,i),','
            jchar=jchar+ishift
            enddo
            jchar=jchar-1
            write(buff(jchar:jchar),'(a1)') '/'
            write(92,10) buff(1:jchar)
            icross=icross+1
         endif
      enddo
      BUFF = 'DATA (IDEN(IHEL),IHEL=???,???) /'
      do i=1,icross-4,4
         write(buff(23:25),'(i3)') i
         write(buff(27:29),'(i3)') i+3
         jchar=33
         ishift = 5
         do j=0,3
            write(buff(jchar:jchar+ishift),'(I4,a)') iden(j+i),','
            jchar=jchar+ishift
         enddo
         jchar=jchar-1
         write(buff(jchar:jchar),'(a1)') '/'
         write(92,10) buff(1:jchar)
      enddo
      write(buff(23:25),'(i3)') i
      write(buff(27:29),'(i3)') icross
      jchar=33
      ishift = 5
      do j=0,icross-i
         write(buff(jchar:jchar+ishift),'(I4,a)') iden(j+i),','
         jchar=jchar+ishift
      enddo
      jchar=jchar-1
      write(buff(jchar:jchar),'(a1)') '/'
c      if (i .le. ncross) then
      if (i .le. icross) then
         write(92,10) buff(1:jchar)
      endif
      
      write(92,5) '----------'
      write(92,5) 'BEGIN CODE'
      write(92,5) '----------'
      write(92,10) 'NTRY=NTRY+1'
      write(92,10) 'DO IPROC=1,NCROSS'
      write(92,10) 'DO IHEL=1,NEXTERNAL'
      write(92,10) '   JC(IHEL) = +1'
      write(92,10) 'ENDDO'

      write(92,10) 'DO IHEL=1,NGRAPHS'
      write(92,10) '    amp2(ihel)=0d0'
      write(92,10) 'ENDDO'
      write(92,10) 'jamp2(0)=dble(NCOLOR)'
      write(92,10) 'DO IHEL=1,int(jamp2(0))'
      write(92,10) '    jamp2(ihel)=0d0'
      write(92,10) 'ENDDO'
      WRITE(92,10) 'ANS(IPROC) = 0D0'
      write(92,11) 'DO IHEL=1,NCOMB'
      write(92,11) '   IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN'
      buff(1:8) = '        '
      buff(8:9)= 'T='
      buff(10:10+flen-1) = fname(1:flen)  
      write(92,11),buff(1:10+flen-1)
      write(92,11) '     ANS(IPROC)=ANS(IPROC)+T'
      write(92,13)' IF (T .GT. 0D0 .AND. .NOT. '//
     &                        'GOODHEL(IHEL,IPROC)) THEN'
      write(92,11) '         GOODHEL(IHEL,IPROC)=.TRUE.'
      WRITE(92,11) '         NGOOD = NGOOD +1'
      WRITE(92,11) '         IGOOD(NGOOD) = IHEL'
      write(92,11) '     ENDIF'
      write(92,11) '   ENDIF'
      write(92,11) 'ENDDO'
      write(92,10) 'ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))'
      write(92,10) 'ENDDO'
      write(92,10) 'END'
      write(92,10) ' '
      write(92,10) ' '

 5    format('C ',A)
 10   format(6X,A)
 11   FORMAT(10X,A)
 12   FORMAT(5X,'&',4X,A)
 13   FORMAT(14X,A)
 14   FORMAT(5X,'&',8X,A)
 15   format('C',6X,A)
 16   FORMAT(5X,'&',5X,A)
 20   format(', P',I1)
 30   format(',P',I1,'(0:3)')
 40   format(6X,A,5(I2,A),A8,A,A8,A,A8,A)
 41   format(6X,A18,10(I2,A))
      end


      subroutine writesub_f90(ngraphs,nexternal,nqed,nqcd)
!***************************************************************************
!     writes out the subroutine for Fortran 90 
!***************************************************************************
      implicit none

! Constants


! Arguments

      integer ngraphs,nqed,nqcd,nexternal

! Local Variables


!-----------
! Begin Code
!-----------
      end


      subroutine writesum_f90
!***************************************************************************
!     writes out the subroutine - f90 version
!***************************************************************************
      implicit none

      end

      
      subroutine sethel(nhel,step,npart,ihel)
c************************************************************************
c     Sets up all possible helicity combinations in nhel(ipart,ihel)
c     given the number of particles, and the step size to be used
c     for each particle
c     The first state is always -1, the next is 0 or 1 depending onstep
c************************************************************************
      implicit none

!     Constants
      
      include 'params.inc'

!     Arguments
      
      integer nhel(maxlines,3**maxlines),step(maxlines),npart,ihel

!     Local

      integer ipart,iflip

c**************KEK**********************   
      integer ql,ql2,ql3
c**************KEK**********************

      logical gotone

      ihel=1
      do ipart=1,npart
c*****************KEK1**************************
         if(step(ipart).gt.0) then
         ql=0
         else 
         ql=-1
         endif 
c****************KEK1***************************
         nhel(ipart,ihel)=-1+ql
      enddo
      gotone=.true.
      do while (gotone)
c         write(*,'(8I4)'),ihel,(nhel(ipart,ihel),ipart=1,npart)
         ihel=ihel+1
         do ipart=1,npart
            nhel(ipart,ihel)=nhel(ipart,ihel-1)
         enddo
         gotone = .false.
         iflip = npart
         do while(.not. gotone .and. iflip .gt. 0)
c*****************KEK2**************************
         if(step(iflip).gt.0) then
         ql=0
         ql2=1
         else
         ql=-1  
         ql2=-1
         endif 
c****************KEK2***************************
            nhel(iflip,ihel) = nhel(iflip,ihel)+step(iflip)*ql2
            ql3=1-ql
            if (nhel(iflip,ihel) .gt. ql3) then
               nhel(iflip,ihel) = -1+ql
               iflip =iflip-1
            else
               gotone=.true.
            endif
         enddo
      enddo
      ihel=ihel-1
      end


!-------------------------------------
      subroutine lower_case(string)
      implicit none

! Constants

      integer    char_a,              char_z
c      parameter (char_a = ichar('A'), char_z = ichar('Z'))
      integer    lcshift
c      parameter (lcshift = ichar('a') - ichar('A'))


! Argument

      character*(*) string

! Local Variable

      integer i

!-----------
! Begin Code
!-----------
      char_a = ichar('A')
      char_z = ichar('Z')
      lcshift = ichar('a')-ichar('A')
      do i=1,len(string)
         if (ichar(string(i:i)) .ge. char_a .and.
     &       ichar(string(i:i)) .le. char_z) then
            string(i:i)=char(ichar(string(i:i))+lcshift)
         end if
      end do
      end

!-------------------------------------
      subroutine upper_case(string)
      implicit none

! Constants

      integer    char_a,              char_z
c      parameter (char_a = ichar('a'), char_z = ichar('z'))
      integer    ucshift
c      parameter (ucshift = ichar('A') - ichar('a'))

! Argument

      character*(*) string

! Local Variable

      integer i

!-----------
! Begin Code
!-----------
      char_A = ichar('a')
      char_Z = ichar('z')
      ucshift = ichar('A')-ichar('a')
c      write(*,*) string
      do i=1,len(string)
         if (ichar(string(i:i)) .ge. char_a .and.
     &       ichar(string(i:i)) .le. char_z) then
            string(i:i)=char(ichar(string(i:i))+ucshift)
         end if
      end do
c      write(*,*) string
      end


      subroutine writebuff2(buff,lun)

      implicit none

! Arguments

      character*(*) buff
      integer lun

! Local Variables

      integer length,i

!-----------
! Begin Code
!-----------
      length = len(buff)
      if (length .le. 60) then
         write (lun,100) buff
      else
         i = 60
         do while (buff(i:i) .ne. '+' .and. buff(i:i).ne.'-')
            i = i - 1
         end do
         if (length-i .gt. 60) then
            print*,'Sorry buff too long in writebuff2',length,i
            print*,buff(:i)
            print*,buff(i:i+60)
            print*,buff(i+60:)
         endif
         write (lun,100) buff(:i)
         write (lun,200) buff(i+1:)
      end if
 100  format ('      ',a)
 200  format ('     &     ',a)
      end

      

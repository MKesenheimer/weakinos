c----------------------------------------------------------------------
C  couplings.f 
c----------------------------------------------------------------------
c  This files takes the inputs of the standard model from a Les Houches 
c  file (param_card.dat) and calculates all the couplings that end up
c  in the Feynman rules, i.e., in the HELAS routines.
c   
c  With the convention of the New Web Generation MadEvent in this
c  part no non-trivial computation is made. The SM model secondary
c  parameters have been all calculated by the SM calculator, SMCalc
c  which produced param_card.dat.
c
c  The only exception to the above rule is for alpha_S. In the case
c  where a pdf is used in the initial state, then the values as(MZ)
c  set in the les houches card is superseeded.
c  Scales and pdf's are set in the run_card.dat.
c
c This file contains the following routines:
c 
c- madgraph original call to set the parameters
c- lh_readin is called from here.
c  It talks to the lh_readin through the common block values.
c      subroutine setpara
c
c-routine to read the LH parameters in
c      subroutine lh_readin
c
c-to set
c      subroutine set_it(npara,ivalue,value,name,id,
c      subroutine case_trap(string,length)
c      subroutine no_spaces(buff,nchars)
c---------------------------------------------------------------------- 


      subroutine setpara(param_name,readlha)
c***********************************************************************
c This subroutine sets up the HELAS couplings of the STANDARD MODEL.
c***********************************************************************
      implicit none
c
c local
c
      character*(*) param_name
      logical readlha
      integer i
      real*8 dum
c
c     common file with the couplings
c
      include 'coupl.inc'
c
c     local
c
      double precision  ee, ee2, ez, ey, sw, cw, sc2, sin2w, wm
      double precision  gwne, gwud, lambda, lam4, xt, rew, rqcd
      double precision  alphas, alfa, alfaw, mfrun
      external          alphas, alfa, alfaw, mfrun
c
c     Common to lh_readin and printout
c
      double precision  alpha, gfermi, alfas
      double precision  mtMS,mbMS,mcMS,mtaMS!MSbar masses
      double precision  Vud,Vus             !CKM matrix elements
      common/values/    alpha,gfermi,alfas,   
     &                  mtMS,mbMS,mcMS,mtaMS,
     &                  Vud
c
c constants
c
      double complex  ci
      parameter( ci = ( 0.0d0, 1.0d0 ) )
      double precision  Zero, One, Two, Three, Four, Half, Rt2
      parameter( Zero = 0.0d0, One = 1.0d0, Two = 2.0d0 )
      parameter( Three = 3.0d0, Four = 4.0d0, Half = 0.5d0 )
      parameter( Rt2   = 1.414213562d0 )
      double precision  Pi, Fourpi
      parameter( Pi = 3.14159265358979323846d0 )
      parameter( Fourpi = Four * Pi )
c
c------------------------------------------
c Start calculating the couplings for HELAS
c------------------------------------------
c
      if(readlha) then 
         call lh_readin(param_name)
         G = DSQRT(4d0*PI*ALFAS) ! use setting of the param_card.dat @ NLO
c
c auxiliary local values
c
      wm = sqrt(zmass**2/Two+
     $     sqrt(zmass**4/Four-Pi/Rt2*alpha/gfermi*zmass**2))
      sin2w  = One-(wm/zmass)**2
      cw  = sqrt( One - sin2w )
      ee2 = alpha * Fourpi
      sw  = sqrt( sin2w )
      ee  = sqrt( ee2 )
      ez  = ee/(sw*cw)
      ey  = ee*(sw/cw)
      sc2 = sin2w*( One - sin2w )
      v   = Two*wm*sw/ee   ! the wmass is used to calculate v
      lambda = hmass**2 / (Two * v**2)
c
c vector boson couplings
c
      gw   = ee/sw
      gwwa = ee
      gwwz = ee*cw/sw
c
c gauge & higgs boson coupling constants
c
      gwwh  = dcmplx( ee2/sin2w*Half*v, Zero )
      gzzh  = dcmplx( ee2/sc2*Half*v, Zero )
      ghhh  = dcmplx( -hmass**2/v*Three, Zero )
      gwwhh = dcmplx( ee2/sin2w*Half, Zero )
      gzzhh = dcmplx( ee2/sc2*Half, Zero)
      ghhhh = ghhh/v
c
c fermion-fermion-vector couplings
c
      gal(1) = dcmplx(  ee          , Zero )
      gal(2) = dcmplx(  ee          , Zero )
      gau(1) = dcmplx( -ee*Two/Three, Zero )
      gau(2) = dcmplx( -ee*Two/Three, Zero )
      gad(1) = dcmplx(  ee/Three    , Zero )
      gad(2) = dcmplx(  ee/Three    , Zero )

      gwf(1) = dcmplx( -ee/sqrt(Two*sin2w), Zero )
      gwf(2) = dcmplx(  Zero              , Zero )

      gzn(1) = dcmplx( -ez*Half                     , Zero )
      gzn(2) = dcmplx(  Zero                        , Zero )
      gzl(1) = dcmplx( -ez*(-Half + sin2w)          , Zero )
      gzl(2) = dcmplx( -ey                          , Zero )
      gzu(1) = dcmplx( -ez*( Half - sin2w*Two/Three), Zero )
      gzu(2) = dcmplx(  ey*Two/Three                , Zero )
      gzd(1) = dcmplx( -ez*(-Half + sin2w/Three)    , Zero )
      gzd(2) = dcmplx( -ey/Three                    , Zero )
c
c fermion-fermion-Higgs couplings (complex) hff(2)
c
c NOTE: the running mass is evaluated @ the same order 
c nloop of alpha_s set by the PDF choice
c 

      if(mtMS.gt.1d0) then
         ghtop(1) = dcmplx( -mtMS/v, Zero )
      else
         ghtop(1) = dcmplx( Zero,Zero)
      endif
      ghtop(2) = ghtop(1)

      if(mbMS.gt.1d0) then
         ghbot(1) = dcmplx( -mbMS/v, Zero )
      else
         ghbot(1) = dcmplx( Zero, Zero )
      endif
      ghbot(2) = ghbot(1)
      
      if(mcMS.gt.1d0) then
         ghcha(1) = dcmplx( -mcMS/v, Zero )
      else
         ghcha(1) = dcmplx( Zero, Zero )
      endif
      ghcha(2) = ghcha(1)

      ghtau(1) = dcmplx( -mtaMS/v, Zero )
      ghtau(2) = ghtau(1)

c---------------------------------------------------------
c Set Photon Width to Zero, used by symmetry optimization
c---------------------------------------------------------

      awidth = 0d0
            
      endif
c     
c     Strong coupling
c
c     As a rule we first check if a pdf has been chosen in the    
c     run_card.dat (which has been already read at this stage).
c     If there pdfs in the initial state, then the alpha_s(MZ) used
c     is set to the corresponding value.  

      GG(1) = -G
      GG(2) = -G     

c----------------------------
c end subroutine coupsm
c----------------------------


      return
      end
 
c$$$  	   subroutine lh_readin(param_name)
c$$$c----------------------------------------------------------------------
c$$$c Read the parameters from the lh file
c$$$c
c$$$c 1. Input values for the EW sector 
c$$$c 2. Higgs mass and width
c$$$c 3. Fermion masses (pole and MSbar) and widths
c$$$c---------------------------------------------------------------------- 
c$$$      implicit none
c$$$
c$$$c
c$$$c     parameters
c$$$c
c$$$      integer maxpara
c$$$      parameter (maxpara=100)
c$$$      double precision  Two, Four, Rt2, Pi
c$$$      parameter( Two = 2.0d0, Four = 4.0d0 )
c$$$      parameter( Rt2   = 1.414213562d0 )
c$$$      parameter( Pi = 3.14159265358979323846d0 )
c$$$c
c$$$c     local
c$$$c
c$$$      character*(*) param_name
c$$$      integer npara,l1,l2,id
c$$$	  character*132 buff
c$$$	  real*8 real_value 
c$$$	  real*8 value(maxpara)
c$$$	  integer ivalue(maxpara),n
c$$$	  character*20 name(maxpara),bn
c$$$	  logical block_found,done,fopened	  
c$$$      integer iunit,i,name_length,idum
c$$$c
c$$$c	block info
c$$$c      
c$$$      character*20 block_name
c$$$c
c$$$c   Common
c$$$c
c$$$	  include 'coupl.inc'
c$$$c
c$$$c     Common to lh_readin and printout
c$$$c
c$$$      double precision  alpha, gfermi, alfas
c$$$      double precision  mtMS,mbMS,mcMS,mtaMS!MSbar masses
c$$$      double precision  Vud,Vus             !CKM matrix elements
c$$$      common/values/    alpha,gfermi,alfas,   
c$$$     &                  mtMS,mbMS,mcMS,mtaMS,
c$$$     &                  Vud
c$$$c
c$$$c----------
c$$$c     start
c$$$c----------
c$$$c
c$$$c     open file
c$$$c
c$$$      iunit=14
c$$$      call open_file_mdl(iunit,param_name,fopened)
c$$$	   done=.false.
c$$$ 
c$$$       n=0
c$$$	   do while(.not.done)
c$$$	   block_found=.false.
c$$$c
c$$$c looks for the blocks or for decay
c$$$c
c$$$      do while(.not.block_found)
c$$$       read(iunit,'(a132)',end=99,err=99) buff
c$$$c--	change to lower case
c$$$	   call case_trap(buff,20)
c$$$       if(buff(1:5).eq.'block') then
c$$$         if(index(buff,"#").ne.0) l1=index(buff,"#")-1 
c$$$         block_name=buff(6:min(l1,26)) 
c$$$         call no_spaces(block_name,name_length)
c$$$c        write(*,*) block_name(1:name_length)
c$$$         block_found=.true.        
c$$$	   elseif(buff(1:5).eq.'decay') then
c$$$               n=n+1
c$$$               l1=60
c$$$               if(index(buff,"#").ne.0) l1=index(buff,"#")-1 ! ignore comments       
c$$$               read(buff(6:l1),*) ivalue(n),value(n)
c$$$               name(n)="decay" 
c$$$       endif		   
c$$$      end do
c$$$c
c$$$c      
c$$$      if(block_found) then
c$$$ 	  do while(.true.)
c$$$          read(iunit,'(a132)',end=99,err=99) buff
c$$$     	  call case_trap(buff,20)
c$$$          if(buff(1:1).eq.'b'.or.buff(1:1).eq.'d') then
c$$$          	backspace iunit
c$$$          	exit
c$$$          endif	         	
c$$$            if(buff(1:1).ne.'#') then  !if it not a comment
c$$$              n=n+1	       
c$$$              l1=60
c$$$              if(index(buff,"#").ne.0) l1=index(buff,"#")-1 ! ignore comments       
c$$$c
c$$$c  WARNING:... not all blocks have the same sintax!! You need to change it
c$$$c  depending on the block you are reading
c$$$c
c$$$             if(block_name(1:5).eq."mgckm") then
c$$$                read(buff(1:l1),*) ivalue(n),idum,value(n)
c$$$              else
c$$$                read(buff(1:l1),*) ivalue(n),value(n)
c$$$              endif  
c$$$              name(n)=block_name(1:name_length)
c$$$c              write(*,"(1x,i2,2x,e16.8,1x,a)") 
c$$$c     &        ivalue(n),value(n),name(n)
c$$$           	  endif
c$$$      end do ! do while in the block
c$$$      else
c$$$	  done=.true.
c$$$	  endif
c$$$	  end do ! do while the entire file
c$$$	  
c$$$
c$$$99	continue      
c$$$	
c$$$	   bn="sminputs"
c$$$       call set_it(n,ivalue,value,name,1,bn,alpha,128.9d0)
c$$$       alpha=1d0/alpha
c$$$       call set_it(n,ivalue,value,name,2,bn,gfermi,0.1166d-4)
c$$$       call set_it(n,ivalue,value,name,3,bn,alfas,0.119d0)
c$$$       call set_it(n,ivalue,value,name,4,bn,zmass,91.188d0)
c$$$       call set_it(n,ivalue,value,name,6,bn,tmass,174.3d0)
c$$$       call set_it(n,ivalue,value,name,7,bn,lmass,0d0)
c$$$       bn="mgyukawa"
c$$$       call set_it(n,ivalue,value,name,4,bn,mcMS,0d0)
c$$$       call set_it(n,ivalue,value,name,5,bn,mbMS,4.2d0)
c$$$       call set_it(n,ivalue,value,name,6,bn,mtMS,174d0)
c$$$       call set_it(n,ivalue,value,name,15,bn,mtaMS,0d0)
c$$$       bn="mgckm"
c$$$       call set_it(n,ivalue,value,name,1,bn,vud,1d0)
c$$$       bn="mass"
c$$$       call set_it(n,ivalue,value,name,4,bn,cmass,0d0)
c$$$       call set_it(n,ivalue,value,name,5,bn,bmass,4.7d0)
c$$$       call set_it(n,ivalue,value,name,6,bn,tmass,tmass*1d0)
c$$$       call set_it(n,ivalue,value,name,15,bn,lmass,lmass*1d0)
c$$$       call set_it(n,ivalue,value,name,25,bn,hmass,120d0)
c$$$       call set_it(n,ivalue,value,name,23,bn,zmass,zmass*1d0)
c$$$       call set_it(n,ivalue,value,name,24,bn,wmass,sqrt(zmass**2/Two+
c$$$     $     sqrt(zmass**4/Four-Pi/Rt2*alpha/gfermi*zmass**2)))
c$$$       bn="decay"
c$$$       call set_it(n,ivalue,value,name,6,bn,twidth,1.5083d0)
c$$$       call set_it(n,ivalue,value,name,25,bn,hwidth,0.0037d0)
c$$$       call set_it(n,ivalue,value,name,23,bn,zwidth,2.441d0)
c$$$       call set_it(n,ivalue,value,name,24,bn,wwidth,2.0476d0)
c$$$     
c$$$      return
c$$$      end

      
      subroutine set_it(npara,ivalue,value,name,id,
     &                  block_name,var,def_value)
c----------------------------------------------------------------------------------
c     finds the parameter value  in block_name and associate var to it.
c     If it is not found a default is given.
c----------------------------------------------------------------------------------
      implicit none

c
c     parameters
c
      integer maxpara
      parameter (maxpara=100)
c
c     arguments
c
      integer npara,ivalue(maxpara),id
      character*20  block_name,name(maxpara)
      real*8 var,def_value,value(maxpara)
c
c     local
c
      logical found
      integer i
c
c     start
c
	  found=.false.
      do i=1,npara
         found = (id.eq.ivalue(i)).and.(name(i).eq.block_name)
 	               if(found) then
         	var=value(i)
            exit
          endif	
      enddo
      
      if (.not.found) then
c         write (*,*) "Warning: parameter not found"
c         write (*,*) "         setting it to default value ",def_value
         var=def_value
      endif
      return

      end
      
      
      subroutine case_trap(string,length)
c**********************************************************    
c change string to lowercase if the input is not
c**********************************************************
      implicit none
c
c     ARGUMENT
c      
      character*(*) string
      integer length
c
c     LOCAL
c
      integer i,k

      do i=1,length
         k=ichar(string(i:i))
         if(k.ge.65.and.k.le.90) then  !upper case A-Z
            k=ichar(string(i:i))+32   
            string(i:i)=char(k)        
         endif
      enddo

      return
      end


      subroutine no_spaces(buff,nchars)
c**********************************************************************
c     Given buff a buffer of words separated by spaces
c     returns it where all space are moved to the right
c     returns also the length of the single word.
c     maxlength is the length of the buffer
c**********************************************************************
      implicit none
c
c     Constants
c
      integer    maxline
      parameter (maxline=20)
      character*1 null
      parameter  (null=' ')
c
c     Arguments
c
      character*(maxline) buff
      integer nchars,maxlength
c
c     Local
c
      integer i,j
      character*(maxline) temp
c-----
c  Begin Code
c-----
      nchars=0
c      write (*,*) "buff=",buff(1:maxlength)
      do i=1,maxline
         if(buff(i:i).ne.null) then
            nchars=nchars+1
            temp(nchars:nchars)=buff(i:i)
         endif
c         write(*,*) i,":",buff(1:maxlength),":",temp(1:nchars),":"
      enddo
      buff=temp      
      end

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
 30      continue
         if (i.eq.5)then
            write(*,*) 'Error: file ',tempname,' cannot be opened'
            call exit(-1)
         else
            write(*,*) 'Warning: file ',tempname,' cannot be opened'
         endif
         tempname='../'//tempname
      enddo


      return
      end
 

 

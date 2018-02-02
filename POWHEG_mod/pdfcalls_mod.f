c MK: copied and modified version of pdfcalls.f, revision 3154
c changes marked with "! MK:"

      subroutine pdfcall(ih,x,pdf)
      implicit none
      include 'pwhg_pdf.h'
      integer ih
      real * 8 x,pdf(-pdf_nparton:pdf_nparton)
      include 'pwhg_st.h'
      if(x .ge. 1) then
         if(x-1 .gt. 1d-4) then
            write(*,*) 'pdfcall: warning, x=',x
            write(*,*) 'returning pdf=0'
         endif
         pdf = 0
         return
      endif

#ifdef DEBUGQ
      !x = 2.0393373472149643D-3
      x = 1.6480036108839822D-3
#define CHECK
#endif

#ifdef DEBUGQ
      if(ih.eq.1) then
         call genericpdf0(pdf_ndns1,pdf_ih1,77497.7D0,x,pdf) ! mu2 = 278.384^2
      elseif(ih.eq.2) then
         call genericpdf0(pdf_ndns2,pdf_ih2,77497.7D0,x,pdf) ! mu2 = 278.384^2
      else
         write(*,*) ' pdfcall: invalid call, ih=',ih
         stop
      endif
#else
      if(ih.eq.1) then
         call genericpdf0(pdf_ndns1,pdf_ih1,st_mufact2,x,pdf)
      elseif(ih.eq.2) then
         call genericpdf0(pdf_ndns2,pdf_ih2,st_mufact2,x,pdf)
      else
         write(*,*) ' pdfcall: invalid call, ih=',ih
         stop
      endif
#endif

#ifdef DEBUGQ
      ! gluon pdf equal zero (test for no gluons in the initial state)
      !pdf(0) = 0D0
      ! u, d, s, c, b, t pdfs equal to zero
      !pdf(1) = 0D0
      !pdf(2) = 0D0      
      !pdf(3) = 0D0
      !pdf(4) = 0D0
      !pdf(5) = 0D0
      !pdf(6) = 0D0
      !pdf(-1) = 0D0  
      !pdf(-2) = 0D0  
      !pdf(-3) = 0D0
      !pdf(-4) = 0D0
      !pdf(-5) = 0D0
      !pdf(-6) = 0D0
      ! test u ubar -> u ubar
      pdf(0) = 0D0
      pdf(1) = 0D0    
      pdf(3) = 0D0
      pdf(4) = 0D0
      pdf(5) = 0D0
      pdf(6) = 0D0
      pdf(-1) = 0D0   
      pdf(-3) = 0D0
      pdf(-4) = 0D0
      pdf(-5) = 0D0
      pdf(-6) = 0D0
#endif
      
#ifdef CHECK
      !print*,"ih",ih
      print*,"mu",dsqrt(st_mufact2)
      print*,"x",x
      print*,"pdf(-6)",pdf(-6)
      print*,"pdf(-5)",pdf(-5) 
      print*,"pdf(-4)",pdf(-4) 
      print*,"pdf(-3)",pdf(-3) 
      print*,"pdf(-2)",pdf(-2) 
      print*,"pdf(-1)",pdf(-1) 
      print*,"pdf(0)",pdf(0)
      print*,"pdf(1)",pdf(1)
      print*,"pdf(2)",pdf(2)
      print*,"pdf(3)",pdf(3)
      print*,"pdf(4)",pdf(4)
      print*,"pdf(5)",pdf(5)
      print*,"pdf(6)",pdf(6)
      stop
#endif
      end


c Front end to genericpdf; it stores the arguments and return values of
c the nrec most recent calls to genericpdf. When invoked it looks in the
c stored calls; if a match its found, its return value is used.
c In this framework it is found that nrec=8 would be enough.
c This provides a remarkable increase in speed (better than a factor of 3)
c when cteq6 pdf are used.
      subroutine genericpdf0(ns,ih,xmu2,x,fx)
      implicit none
      include 'pwhg_pdf.h'
      integer maxparton
      parameter (maxparton=22)
      integer ns,ih
      real * 8 xmu2,x,fx(-pdf_nparton:pdf_nparton)
      integer nrec
      parameter (nrec=10)
      real * 8 oxmu2(nrec),ox(nrec),ofx(-maxparton:maxparton,nrec)
      integer ons(nrec),oih(nrec)
      integer irec
      save oxmu2,ox,ofx,ons,oih,irec
c set to impossible values to begin with
      data ox/nrec*-1d0/
      data irec/0/
      integer j,k
      real * 8 charmthr2,bottomthr2
      logical ini
      data ini/.true./
      save ini,charmthr2,bottomthr2
      real * 8 powheginput
      external powheginput
      if(ini) then
         charmthr2=powheginput('#charmthrpdf')
         bottomthr2=powheginput('#bottomthrpdf')
         if(charmthr2.lt.0) charmthr2=1.5
         if(bottomthr2.lt.0) bottomthr2=5
         charmthr2=charmthr2**2
         bottomthr2=bottomthr2**2
         ini=.false.
      endif
      do j=irec,1,-1
         if(x.eq.ox(j)) then
            if(xmu2.eq.oxmu2(j)) then
               if(ns.eq.ons(j).and.ih.eq.oih(j)) then
                  fx=ofx(-pdf_nparton:pdf_nparton,j)
                  return
               endif
            endif
         endif
      enddo
      do j=nrec,irec+1,-1
         if(x.eq.ox(j)) then
            if(xmu2.eq.oxmu2(j)) then
               if(ns.eq.ons(j).and.ih.eq.oih(j)) then
                  fx=ofx(-pdf_nparton:pdf_nparton,j)
                  return
               endif
            endif
         endif
      enddo
      irec=irec+1
      if(irec.gt.nrec) irec=1
      ons(irec)=ns
      oih(irec)=ih
      oxmu2(irec)=xmu2
      ox(irec)=x
      call genericpdf(ns,ih,xmu2,x,ofx(-pdf_nparton:pdf_nparton,irec))
c Flavour thresholds:
      if(xmu2.lt.bottomthr2) then
         ofx(5,irec)=0
         ofx(-5,irec)=0
      endif
      if(xmu2.lt.charmthr2) then
         ofx(4,irec)=0
         ofx(-4,irec)=0
      endif
      do k=-pdf_nparton,pdf_nparton
         if (ofx(k,irec).lt.0) then
            call increasecnt("negative pdf values");
            ofx(k,irec)=0
         endif
      enddo
      fx=ofx(-pdf_nparton:pdf_nparton,irec)
      end


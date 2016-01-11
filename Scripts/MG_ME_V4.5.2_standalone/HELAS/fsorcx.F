      subroutine fsorcx(roc,sc,gc,fmass,fwidth , fsorc)
c
c This subroutine computes an off-shell antifermion wavefunction from a
c flowing-OUT antiRarita-Schwinger fermion and a scalar boson.
c
c input:
c       complex roc(18)        : flow-out antiRSfermion           v<roc|
c       complex sc(6)          : input    scalar                   s(q)
c       complex gc(2)          : coupling constants                GFRS 
c       real    fmass          : mass  of OUTPUT fermion f'
c       real    fwidth         : width of OUTPUT fermion f'
c
c output:
c       complex fsorc(6)        : off-shell fermion          <roc,sqv,f'|
c
c- by Yoshitaro Takaesu - 2011/03/25
c
      implicit none
      double complex roc(18),sc(3),fsorc(6),gc(2),sq1(2,2),sq2(2,2),ds
      double complex rsq(0:3,4),rsqg(4),sp1(2,2),sp2(2,2),cImag
      double complex rc(0:3,4),gcc(2)
      double precision pf(0:3),fmass,fwidth,pf2,q(0:3)
      integer i

      cImag = ( 0.0d0, 1.0d0 )

#ifdef HELAS_CHECK
      double precision rZero, cZero
      parameter( rZero = 0.0d0, cZero = ( 0.0d0, 0.0d0 ) )
      integer stdo
      parameter( stdo = 6, )
#endif
c
#ifdef HELAS_CHECK
      if (abs(roc(1))+abs(roc(2))+abs(roc(3))+abs(roc(4))+abs(roc(5))
     & +abs(roc(6))+abs(roc(7))+abs(roc(8))abs(roc(9))+abs(roc(10))
     & +abs(roc(11))+abs(roc(12))+abs(roc(13))+abs(roc(14))+abs(roc(15))
     & +abs(roc(16)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : roc in fsorcx is zero spinor'
      endif
      if ( abs(roc(17))+abs(roc(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : roc in fsorcx has zero momentum'
      endif
      if ( sc(1).eq.cZero ) then
         write(stdo,*) ' helas-warn  : sc in fsorcx is zero scalar'
      endif
      if ( abs(sc(2))+abs(sc(3)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : sc in fsorcx has zero momentum'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*) ' helas-error : gc in fsorcx is zero coupling'
      endif
      if ( fmass.lt.rZero ) then
         write(stdo,*) ' helas-error : fmass in fsorcx is negative'
         write(stdo,*) '             : fmass = ',fmass
      endif
      if ( fwidth.lt.rZero ) then
         write(stdo,*) ' helas-error : fwidth in fsorcx is negative'
         write(stdo,*) '               fwidth = ',fwidth
      endif
#endif

      fsorc(5) = roc(17)+sc(2)
      fsorc(6) = roc(18)+sc(3)

      pf(0) = dble( fsorc(5))
      pf(1) = dble( fsorc(6))
      pf(2) = dimag(fsorc(6))
      pf(3) = dimag(fsorc(5))
      pf2 = pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)

#ifdef HELAS_CHECK
      if ( abs(fsorc(5))+abs(fsorc(6)).eq.rZero ) then
          write(stdo,*)
     &        ' helas-error : fsorc in fsorcx has zero momentum'
       endif
       if ( fwidth.eq.rZero .and. pf2.eq.fmass**2 ) then
          write(stdo,*)
     &         ' helas-error : fsorc in fsorcx is on fmass pole'
          write(stdo,*)
     &         '             : p     = ',pf(0),pf(1),pf(2),pf(3)
          write(stdo,*)
     &         '             : abs(p)= ',sqrt(abs(pf2))
         fsorc(1) = cZero
         fsorc(2) = cZero
         fsorc(3) = cZero
         fsorc(4) = cZero
         return
      endif
#endif

      gcc(1) = dconjg(gc(2))
      gcc(2) = dconjg(gc(1)) 

      q(0) = dble( sc(2))
      q(1) = dble( sc(3))
      q(2) = dimag(sc(3))
      q(3) = dimag(sc(2))

      rc(0,1) = roc(1)
      rc(0,2) = roc(2)
      rc(0,3) = roc(3)
      rc(0,4) = roc(4)
      rc(1,1) = roc(5)
      rc(1,2) = roc(6)
      rc(1,3) = roc(7)
      rc(1,4) = roc(8)
      rc(2,1) = roc(9)
      rc(2,2) = roc(10)
      rc(2,3) = roc(11)
      rc(2,4) = roc(12)
      rc(3,1) = roc(13)
      rc(3,2) = roc(14)
      rc(3,3) = roc(15)
      rc(3,4) = roc(16)

      sq1(1,1) = q(0)-q(3)
      sq1(1,2) = -(q(1)-cImag*q(2))
      sq1(2,1) = -(q(1)+cImag*q(2))
      sq1(2,2) = q(0)+q(3)

      sq2(1,1) = q(0)+q(3)
      sq2(1,2) = q(1)-cImag*q(2)
      sq2(2,1) = q(1)+cImag*q(2)
      sq2(2,2) = q(0)-q(3)

      do i=0,3
        rsq(i,1) = rc(i,3)*sq2(1,1)+rc(i,4)*sq2(2,1)
	rsq(i,2) = rc(i,3)*sq2(1,2)+rc(i,4)*sq2(2,2)
	rsq(i,3) = rc(i,1)*sq1(1,1)+rc(i,2)*sq1(2,1)
	rsq(i,4) = rc(i,1)*sq1(1,2)+rc(i,2)*sq1(2,2)
      enddo

      rsqg(1)=rsq(0,3)+rsq(1,4)+cImag*rsq(2,4)+rsq(3,3)
      rsqg(2)=rsq(0,4)+rsq(1,3)-cImag*rsq(2,3)-rsq(3,4)
      rsqg(3)=rsq(0,1)-rsq(1,2)-cImag*rsq(2,2)-rsq(3,1)
      rsqg(4)=rsq(0,2)-rsq(1,1)+cImag*rsq(2,1)+rsq(3,2)

      sp1(1,1) = pf(0)-pf(3)
      sp1(1,2) = -(pf(1)-cImag*pf(2))
      sp1(2,1) = -(pf(1)+cImag*pf(2))
      sp1(2,2) = pf(0)+pf(3)

      sp2(1,1) = pf(0)+pf(3)
      sp2(1,2) = pf(1)-cImag*pf(2)
      sp2(2,1) = pf(1)+cImag*pf(2)
      sp2(2,2) = pf(0)-pf(3)

      ds = -sc(1)/dcmplx( pf2-fmass**2, fmass*fwidth )
      ds = -ds

      fsorc(1)=ds*( gcc(1)*rsqg(1)*fmass
     &            +gcc(2)*(rsqg(3)*sp2(1,1)+rsqg(4)*sp2(2,1)))
      fsorc(2)=ds*( gcc(1)*rsqg(2)*fmass
     &            +gcc(2)*(rsqg(3)*sp2(1,2)+rsqg(4)*sp2(2,2)))
      fsorc(3)=ds*( gcc(1)*(rsqg(1)*sp1(1,1)+rsqg(2)*sp1(2,1))
     &            +gcc(2)*rsqg(3)*fmass)
      fsorc(4)=ds*( gcc(1)*(rsqg(1)*sp1(1,2)+rsqg(2)*sp1(2,2))
     &            +gcc(2)*rsqg(4)*fmass)
     
      return
      end

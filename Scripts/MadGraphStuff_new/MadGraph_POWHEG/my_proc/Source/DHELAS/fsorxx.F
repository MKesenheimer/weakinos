      subroutine fsorxx(ro,sc,gc,fmass,fwidth , fsor)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-OUT Rarita-Schwinger fermion and a scalar boson.
c
c input:
c       complex ro(18)         : flow-out RS fermion              v<ro|
c       complex sc(6)          : input    scalar                   s(q)
c       complex gc(2)          : coupling constants                GFRS 
c       real    fmass          : mass  of OUTPUT fermion f'
c       real    fwidth         : width of OUTPUT fermion f'
c
c output:
c       complex fsor(6)        : off-shell fermion          <ro,sqv,f'|
c
c- by Y.Takaesu - 2009/10/23
c
      implicit none
      double complex ro(18),sc(3),fsor(6),gc(2),sq1(2,2),sq2(2,2),ds
      double complex rsq(0:3,4),rsqg(4),sp1(2,2),sp2(2,2),cImag
      double complex rc(0:3,4)
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
      if (abs(ro(1))+abs(ro(2))+abs(ro(3))+abs(ro(4))+abs(ro(5))
     & +abs(ro(6))+abs(ro(7))+abs(ro(8))abs(ro(9))+abs(ro(10))
     & +abs(ro(11))+abs(ro(12))+abs(ro(13))+abs(ro(14))+abs(ro(15))
     & +abs(ro(16)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : ro in fsorxx is zero spinor'
      endif
      if ( abs(ro(17))+abs(ro(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : ro in fsorxx has zero momentum'
      endif
      if ( sc(1).eq.cZero ) then
         write(stdo,*) ' helas-warn  : sc in fsorxx is zero scalar'
      endif
      if ( abs(sc(2))+abs(sc(3)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : sc in fsorxx has zero momentum'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*) ' helas-error : gc in fsorxx is zero coupling'
      endif
      if ( fmass.lt.rZero ) then
         write(stdo,*) ' helas-error : fmass in fsorxx is negative'
         write(stdo,*) '             : fmass = ',fmass
      endif
      if ( fwidth.lt.rZero ) then
         write(stdo,*) ' helas-error : fwidth in fsorxx is negative'
         write(stdo,*) '               fwidth = ',fwidth
      endif
#endif

      fsor(5) = ro(17)+sc(2)
      fsor(6) = ro(18)+sc(3)

      pf(0) = dble( fsor(5))
      pf(1) = dble( fsor(6))
      pf(2) = dimag(fsor(6))
      pf(3) = dimag(fsor(5))
      pf2 = pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)

#ifdef HELAS_CHECK
      if ( abs(fsor(5))+abs(fsor(6)).eq.rZero ) then
          write(stdo,*)
     &        ' helas-error : fsor in fsorxx has zero momentum'
       endif
       if ( fwidth.eq.rZero .and. pf2.eq.fmass**2 ) then
          write(stdo,*)
     &         ' helas-error : fsor in fsorxx is on fmass pole'
          write(stdo,*)
     &         '             : p     = ',pf(0),pf(1),pf(2),pf(3)
          write(stdo,*)
     &         '             : abs(p)= ',sqrt(abs(pf2))
         fsor(1) = cZero
         fsor(2) = cZero
         fsor(3) = cZero
         fsor(4) = cZero
         return
      endif
#endif

      q(0) = dble( sc(2))
      q(1) = dble( sc(3))
      q(2) = dimag(sc(3))
      q(3) = dimag(sc(2))

      rc(0,1) = ro(1)
      rc(0,2) = ro(2)
      rc(0,3) = ro(3)
      rc(0,4) = ro(4)
      rc(1,1) = ro(5)
      rc(1,2) = ro(6)
      rc(1,3) = ro(7)
      rc(1,4) = ro(8)
      rc(2,1) = ro(9)
      rc(2,2) = ro(10)
      rc(2,3) = ro(11)
      rc(2,4) = ro(12)
      rc(3,1) = ro(13)
      rc(3,2) = ro(14)
      rc(3,3) = ro(15)
      rc(3,4) = ro(16)

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

      fsor(1)=ds*( gc(1)*rsqg(1)*fmass
     &            +gc(2)*(rsqg(3)*sp2(1,1)+rsqg(4)*sp2(2,1)))
      fsor(2)=ds*( gc(1)*rsqg(2)*fmass
     &            +gc(2)*(rsqg(3)*sp2(1,2)+rsqg(4)*sp2(2,2)))
      fsor(3)=ds*( gc(1)*(rsqg(1)*sp1(1,1)+rsqg(2)*sp1(2,1))
     &            +gc(2)*rsqg(3)*fmass)
      fsor(4)=ds*( gc(1)*(rsqg(1)*sp1(1,2)+rsqg(2)*sp1(2,2))
     &            +gc(2)*rsqg(4)*fmass)
     
      return
      end

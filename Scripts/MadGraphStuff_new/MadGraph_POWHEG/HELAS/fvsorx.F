      subroutine fvsorx(ro,vc,sc,gc,fmass,fwidth , fvsor)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-OUT external Rarita-Schwinger fermion, a vector boson, and 
c a scalar.
c
c input:
c       complex ro(18)         : flow-out Rarita-Schwinger fermion   <ro|
c	complex vc(6)          : input    vector boson                v
c       complex sc(6)          : input    scalar                      s
c       complex gc(2)          : coupling constants                 gvsrf
c       real    fmass          : mass  of OUTPUT fermion f'
c       real    fwidth         : width of OUTPUT fermion f'
c
c output:
c       complex fvsor(6)        : off-shell fermion             <ro,v,s,f'|
c
c- by Y.Takaesu - 2009/11/01
c
      implicit none
      double complex ro(18),vc(6),sc(3),fvsor(6),gc(2),sq1(2,2),sq2(2,2)
      double complex rsq(0:3,4),rsqg(4),sp1(2,2),sp2(2,2),cImag
      double complex rc(0:3,4),q(0:3),ds
      double precision pf(0:3),fmass,fwidth,pf2
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
         write(stdo,*) ' helas-warn  : ro in fvsorx is zero spinor'
      endif
      if ( abs(ro(17))+abs(ro(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : ro in fvsorx has zero momentum'
      endif
	if ( abs(vc(1))+abs(vc(2))+abs(vc(3))+abs(vc(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : vc in fvsirx is zero vector'
      endif
      if ( abs(vc(5))+abs(vc(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : vc in fvsirx has zero momentum'
      endif
      if ( sc(1).eq.cZero ) then
         write(stdo,*) ' helas-warn  : sc in fvsorx is zero scalar'
      endif
      if ( abs(sc(2))+abs(sc(3)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : sc in fvsorx has zero momentum'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*) ' helas-error : gc in fvsorx is zero coupling'
      endif
      if ( fmass.lt.rZero ) then
         write(stdo,*) ' helas-error : fmass in fvsorx is negative'
         write(stdo,*) '             : fmass = ',fmass
      endif
      if ( fwidth.lt.rZero ) then
         write(stdo,*) ' helas-error : fwidth in fvsorx is negative'
         write(stdo,*) '               fwidth = ',fwidth
      endif
#endif

      fvsor(5) = ro(17)+vc(5)+sc(2)
      fvsor(6) = ro(18)+vc(6)+sc(3)

      pf(0) = dble( fvsor(5))
      pf(1) = dble( fvsor(6))
      pf(2) = dimag(fvsor(6))
      pf(3) = dimag(fvsor(5))
      pf2 = pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)

#ifdef HELAS_CHECK
      if ( abs(fvsor(5))+abs(fvsor(6)).eq.rZero ) then
          write(stdo,*)
     &        ' helas-error : fvsor in fvsorx has zero momentum'
       endif
       if ( fwidth.eq.rZero .and. pf2.eq.fmass**2 ) then
          write(stdo,*)
     &         ' helas-error : fvsor in fvsorx is on fmass pole'
          write(stdo,*)
     &         '             : p     = ',pf(0),pf(1),pf(2),pf(3)
          write(stdo,*)
     &         '             : abs(p)= ',sqrt(abs(pf2))
         fvsor(1) = cZero
         fvsor(2) = cZero
         fvsor(3) = cZero
         fvsor(4) = cZero
         return
      endif
#endif
      q(0) = vc(1)
      q(1) = vc(2)
      q(2) = vc(3)
      q(3) = vc(4)

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

      ds = -1d0/dcmplx( pf2-fmass**2, fmass*fwidth )

	fvsor(1)=ds*sc(1)*(gc(1)*rsqg(1)*fmass+gc(2)
     & *(rsqg(3)*sp2(1,1)+rsqg(4)*sp2(2,1)))
	fvsor(2)=ds*sc(1)*(gc(1)*rsqg(2)*fmass+gc(2)
     & *(rsqg(3)*sp2(1,2)+rsqg(4)*sp2(2,2)))
      fvsor(3)=ds*sc(1)*(gc(1)*(rsqg(1)*sp1(1,1)+rsqg(2)
     & *sp1(2,1))+gc(2)*rsqg(3)*fmass)
	fvsor(4)=ds*sc(1)*(gc(1)*(rsqg(1)*sp1(1,2)+rsqg(2)
     & *sp1(2,2))+gc(2)*rsqg(4)*fmass)
     
      return
      end

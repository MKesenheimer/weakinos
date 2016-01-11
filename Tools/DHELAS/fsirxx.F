      subroutine fsirxx(ri,sc,gc,fmass,fwidth , fsir)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-IN Rarita-Schwinger fermion and a scalar boson.
c
c input:
c       complex ri(18)         : flow-in  RS fermion               |ri>v
c       complex sc(3)          : input    scalar                    s(q)
c       complex gc(2)          : coupling constants                 GFRS 
c       real    fmass          : mass  of OUTPUT fermion f'
c       real    fwidth         : width of OUTPUT fermion f'
c
c output:
c       complex fsir(6)        : off-shell fermion           |f',svq,ri>
c
c- by Y.Takaesu - 2010/01/18
c  
      implicit none
      double complex ri(18),sc(3),fsir(6),gc(2),ds,rc(0:3,4)
      double precision pf(0:3),fmass,fwidth,pf2,q(0:3)
      double complex sqr(0:3,4),gsqr(4),sq1(2,2),sq2(2,2)
      double complex sp1(2,2),sp2(2,2),cImag
      parameter(cImag = (0.0d0,1.0d0))
      integer i

#ifdef HELAS_CHECK
      double precision rZero, cZero
      parameter( rZero = 0.0d0, cZero = ( 0.0d0, 0.0d0 ) )
      integer stdo
      parameter( stdo = 6, cImag = ( 0.0d0, 1.0d0 ) )
#endif
c
#ifdef HELAS_CHECK
      if (abs(ri(1))+abs(ri(2))+abs(ri(3))+abs(ri(4))+abs(ri(5))
     & +abs(ri(6))+abs(ri(7))+abs(ri(8))abs(ri(9))+abs(ri(10))
     & +abs(ri(11))+abs(ri(12))+abs(ri(13))+abs(ri(14))+abs(ri(15))
     & +abs(ri(16)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : ri in fsirxx is zero spinor'
      endif
      if ( abs(ri(17))+abs(ri(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : ri in fsirxx has zero momentum'
      endif
      if ( sc(1).eq.cZero ) then
         write(stdo,*) ' helas-warn  : sc in fsirxx is zero scalar'
      endif
      if ( abs(sc(2))+abs(sc(3)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : sc in fsirxx has zero momentum'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*) ' helas-error : gc in fsirxx is zero coupling'
      endif
      if ( fmass.lt.rZero ) then
         write(stdo,*) ' helas-error : fmass in fsirxx is negative'
         write(stdo,*) '               fmass = ',fmass
      endif
      if ( fwidth.lt.rZero ) then
         write(stdo,*) ' helas-error : fwidth in fsirxx is negative'
         write(stdo,*) '               fwidth = ',fwidth
      endif
#endif

      fsir(5) = ri(17)-sc(2)
      fsir(6) = ri(18)-sc(3)

      pf(0) = dble( fsir(5))
      pf(1) = dble( fsir(6))
      pf(2) = dimag(fsir(6))
      pf(3) = dimag(fsir(5))
      pf2 = pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)

#ifdef HELAS_CHECK
      if ( abs(fsir(5))+abs(fsir(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fsir in fsirxx has zero momentum'
      endif
      if ( fwidth.eq.rZero .and. pf2.eq.fmass**2 ) then
         write(stdo,*)
     &        ' helas-error : fsir in fsirxx is on fmass pole'
         write(stdo,*)
     &        '             : p     = ',pf(0),pf(1),pf(2),pf(3)
         write(stdo,*)
     &        '             : abs(p)= ',sqrt(abs(pf2))
         fsir(1) = cZero
         fsir(2) = cZero
         fsir(3) = cZero
         fsir(4) = cZero
         return
      endif
#endif

      q(0) = dble( sc(2))
      q(1) = dble( sc(3))
      q(2) = dimag(sc(3))
      q(3) = dimag(sc(2))

      rc(0,1) = ri(1)
      rc(0,2) = ri(2)
      rc(0,3) = ri(3)
      rc(0,4) = ri(4)
      rc(1,1) = ri(5)
      rc(1,2) = ri(6)
      rc(1,3) = ri(7)
      rc(1,4) = ri(8)
      rc(2,1) = ri(9)
      rc(2,2) = ri(10)
      rc(2,3) = ri(11)
      rc(2,4) = ri(12)
      rc(3,1) = ri(13)
      rc(3,2) = ri(14)
      rc(3,3) = ri(15)
      rc(3,4) = ri(16)

      sq1(1,1) = q(0)-q(3)
      sq1(1,2) = -(q(1)-cImag*q(2))
      sq1(2,1) = -(q(1)+cImag*q(2))
      sq1(2,2) = q(0)+q(3)

      sq2(1,1) = q(0)+q(3)
      sq2(1,2) = q(1)-cImag*q(2)
      sq2(2,1) = q(1)+cImag*q(2)
      sq2(2,2) = q(0)-q(3)

      sp1(1,1) = pf(0)-pf(3)
      sp1(1,2) = -(pf(1)-cImag*pf(2))
      sp1(2,1) = -(pf(1)+cImag*pf(2))
      sp1(2,2) = pf(0)+pf(3)

      sp2(1,1) = pf(0)+pf(3)
      sp2(1,2) = pf(1)-cImag*pf(2)
      sp2(2,1) = pf(1)+cImag*pf(2)
      sp2(2,2) = pf(0)-pf(3)

      do i=0,3
        sqr(i,1) = rc(i,3)*sq1(1,1)+rc(i,4)*sq1(1,2)
	sqr(i,2) = rc(i,3)*sq1(2,1)+rc(i,4)*sq1(2,2)
	sqr(i,3) = rc(i,1)*sq2(1,1)+rc(i,2)*sq2(1,2)
	sqr(i,4) = rc(i,1)*sq2(2,1)+rc(i,2)*sq2(2,2)
      enddo

      gsqr(1)=sqr(0,3)-sqr(1,4)+cImag*sqr(2,4)-sqr(3,3)
      gsqr(2)=sqr(0,4)-sqr(1,3)-cImag*sqr(2,3)+sqr(3,4)
      gsqr(3)=sqr(0,1)+sqr(1,2)-cImag*sqr(2,2)+sqr(3,1)
      gsqr(4)=sqr(0,2)+sqr(1,1)+cImag*sqr(2,1)-sqr(3,2)

      ds = sc(1)/dcmplx( pf2-fmass**2, fmass*fwidth )

      fsir(1)=ds*dconjg(gc(1))*(sp1(1,1)*gsqr(3)+sp1(1,2)*gsqr(4))
      fsir(2)=ds*dconjg(gc(1))*(sp1(2,1)*gsqr(3)+sp1(2,2)*gsqr(4))
      fsir(3)=ds*dconjg(gc(1))*fmass*gsqr(3)
      fsir(4)=ds*dconjg(gc(1))*fmass*gsqr(4)

      fsir(1)=fsir(1)+ds*dconjg(gc(2))*fmass*gsqr(1)
      fsir(2)=fsir(2)+ds*dconjg(gc(2))*fmass*gsqr(2)
      fsir(3)=fsir(3)
     &        +ds*dconjg(gc(2))*(sp2(1,1)*gsqr(1)+sp2(1,2)*gsqr(2))
      fsir(4)=fsir(4)
     &        +ds*dconjg(gc(2))*(sp2(2,1)*gsqr(1)+sp2(2,2)*gsqr(2))

      return
      end

      subroutine fvsirc(ric,vc,sc,gc,fmass,fwidth , fvsoc)
c
c This subroutine computes an off-shell antifermion wavefunction from a
c flowing-IN antiRarita-Schwinger fermion, a vector boson and a scalar boson.
c
c input:
c       complex ric(18)        : flow-in antiRSfermion             |ric>
c       complex vc(6)          : input    vector                      v
c       complex sc(3)          : input    scalar                      s
c       complex gc(2)          : coupling constants                 gvsrf
c       real    fmass          : mass  of OUTPUT fermion f'
c       real    fwidth         : width of OUTPUT fermion f'
c
c output:
c       complex fvsoc(6)        : off-shell fermion             |f',v,s,ric>
c
c   
c- by Yoshitaro Takaesu - 2011/03/25
c  
      implicit none
      double complex ric(18),sc(3),fvsoc(6),gc(2),gcc(2),ds,rc(0:3,4),vc(6)
      double precision pf(0:3),fmass,fwidth,pf2
      double complex sqr(0:3,4),gsqr(4),sq1(2,2),sq2(2,2)
      double complex sp1(2,2),sp2(2,2),cImag,q(0:3)
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
      if (abs(ric(1))+abs(ric(2))+abs(ric(3))+abs(ric(4))+abs(ric(5))
     & +abs(ric(6))+abs(ric(7))+abs(ric(8))abs(ric(9))+abs(ric(10))
     & +abs(ric(11))+abs(ric(12))+abs(ric(13))+abs(ric(14))+abs(ric(15))
     & +abs(ric(16)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : ric in fvsirc is zero spinor'
      endif
      if ( abs(ric(17))+abs(ric(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : ric in fvsirc has zero momentum'
      endif
	if ( abs(vc(1))+abs(vc(2))+abs(vc(3))+abs(vc(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : vc in fvsirc is zero vector'
      endif
      if ( abs(vc(5))+abs(vc(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : vc in fvsirc has zero momentum'
      endif
      if ( sc(1).eq.cZero ) then
         write(stdo,*) ' helas-warn  : sc in fvsirc is zero scalar'
      endif
      if ( abs(sc(2))+abs(sc(3)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : sc in fvsirc has zero momentum'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*) ' helas-error : gc in fvsirc is zero coupling'
      endif
      if ( fmass.lt.rZero ) then
         write(stdo,*) ' helas-error : fmass in fvsirc is negative'
         write(stdo,*) '               fmass = ',fmass
      endif
      if ( fwidth.lt.rZero ) then
         write(stdo,*) ' helas-error : fwidth in fvsirc is negative'
         write(stdo,*) '               fwidth = ',fwidth
      endif
#endif

      fvsoc(5) = ric(17)-vc(5)-sc(2)
      fvsoc(6) = ric(18)-vc(6)-sc(3)

      pf(0) = dble( fvsoc(5))
      pf(1) = dble( fvsoc(6))
      pf(2) = dimag(fvsoc(6))
      pf(3) = dimag(fvsoc(5))
      pf2 = pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)

#ifdef HELAS_CHECK
      if ( abs(fvsoc(5))+abs(fvsoc(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fvsoc in fvsirc has zero momentum'
      endif
      if ( fwidth.eq.rZero .and. pf2.eq.fmass**2 ) then
         write(stdo,*)
     &        ' helas-error : fvsoc in fvsirc is on fmass pole'
         write(stdo,*)
     &        '             : p     = ',pf(0),pf(1),pf(2),pf(3)
         write(stdo,*)
     &        '             : abs(p)= ',sqrt(abs(pf2))
         fvsoc(1) = cZero
         fvsoc(2) = cZero
         fvsoc(3) = cZero
         fvsoc(4) = cZero
         return
      endif
#endif

      gcc(1) = dconjg(gc(2))
      gcc(2) = dconjg(gc(1)) 

      ds = -1d0/dcmplx( pf2-fmass**2, fmass*fwidth )

      q(0) = vc(1)
      q(1) = vc(2)
      q(2) = vc(3)
      q(3) = vc(4)

      rc(0,1) = ric(1)
      rc(0,2) = ric(2)
      rc(0,3) = ric(3)
      rc(0,4) = ric(4)
      rc(1,1) = ric(5)
      rc(1,2) = ric(6)
      rc(1,3) = ric(7)
      rc(1,4) = ric(8)
      rc(2,1) = ric(9)
      rc(2,2) = ric(10)
      rc(2,3) = ric(11)
      rc(2,4) = ric(12)
      rc(3,1) = ric(13)
      rc(3,2) = ric(14)
      rc(3,3) = ric(15)
      rc(3,4) = ric(16)

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

      fvsoc(1)=ds*sc(1)*dconjg(gcc(1))*(sp1(1,1)*gsqr(3)
     &        +sp1(1,2)*gsqr(4))
      fvsoc(2)=ds*sc(1)*dconjg(gcc(1))*(sp1(2,1)*gsqr(3)
     &        +sp1(2,2)*gsqr(4))
      fvsoc(3)=ds*sc(1)*dconjg(gcc(1))*fmass*gsqr(3)
      fvsoc(4)=ds*sc(1)*dconjg(gcc(1))*fmass*gsqr(4)

      fvsoc(1)=fvsoc(1)+ds*sc(1)*dconjg(gcc(2))*fmass*gsqr(1)
      fvsoc(2)=fvsoc(2)+ds*sc(1)*dconjg(gcc(2))*fmass*gsqr(2)
      fvsoc(3)=fvsoc(3)+ds*sc(1)*dconjg(gcc(2))*(sp2(1,1)*gsqr(1)
     &        +sp2(1,2)*gsqr(2))
      fvsoc(4)=fvsoc(4)+ds*sc(1)*dconjg(gcc(2))*(sp2(2,1)*gsqr(1)
     &        +sp2(2,2)*gsqr(2))

      return
      end

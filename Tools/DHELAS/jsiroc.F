      subroutine jsiroc(ric,foc,sc,gc,vmass,vwidth , jsirc)
c
c This subroutine computes an off-shell vector current from an anti-
c fermion, antiRarita-Schwinger fermion, and a scalar. The vector
c boson propagator is given in Feynman gauge for a massless vector and
c in unitary gauge for a massive vector.
c
c input:
c       complex ric(18)        : flow-in  antiRSfermion              |ric>
c       complex foc(6)         : flow-out antifermion                <foc|
c       complex gc(2)          : coupling constants                  gfrvs
c       real    vmass          : mass  of OUTPUT vector v
c       real    vwidth         : width of OUTPUT vector v
c
c output:
c       complex jsirc(6)        : vector current          j^mu(<foc|s|ric>)
c   
c- by Yoshitaro Takaesu - 2011/03/27
c  
      implicit none
      double complex foc(6),ric(18),gc(2),sc(3),jsirc(6)
      double precision q(0:3),vmass,vwidth,q2,vm2
      double complex gcc(2),c0,c1,c2,c3,cs,d
      double complex cm2        ! mass**2- I Gamma mass (Fabio)

      double complex sqr(0:3,4),gsqr(4),sq1(2,2),sq2(2,2),fgsqr,fr(4)
      double complex rc(0:3,4)
      integer i,j

      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )

#ifdef HELAS_CHECK
      integer stdo
      parameter( stdo = 6 )
#endif
c
#ifdef HELAS_CHECK
      if ( abs(foc(1))+abs(foc(2))+abs(foc(3))+abs(foc(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : foc in jsiroc is zero spinor'
      endif
      if ( abs(foc(5))+abs(foc(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : foc in jsiroc has zero momentum'
      endif
      if (abs(ric(1))+abs(ric(2))+abs(ric(3))+abs(ric(4))+abs(ric(5))
     & +abs(ric(6))+abs(ric(7))+abs(ric(8))abs(ric(9))+abs(ric(10))
     & +abs(ric(11))+abs(ric(12))+abs(ric(13))+abs(ric(14))+abs(ric(15))
     & +abs(ric(16)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : ric in jsiroc is zero spinor'
      endif
      if ( abs(ric(17))+abs(ric(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : ric in jsiroc has zero momentum'
      endif
	if ( sc(1).eq.cZero ) then
         write(stdo,*) ' helas-warn  : sc in jsiroc is zero scalar'
      endif
      if ( abs(sc(2))+abs(sc(3)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : sc in jsiroc has zero momentum'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*) ' helas-error : gc in jsiroc is zero coupling'
      endif
      if ( vmass.lt.rZero ) then
         write(stdo,*) ' helas-error : vmass in jsiroc is negative'
         write(stdo,*) '             : vmass = ',vmass
      endif
      if ( vwidth.lt.rZero ) then
         write(stdo,*) ' helas-error : vwidth in jsiroc is negative'
         write(stdo,*) '             : vwidth = ',vwidth
      endif
#endif

      jsirc(5) = foc(5)-ric(17)+sc(2)
      jsirc(6) = foc(6)-ric(18)+sc(3)

      q(0) = dble( jsirc(5))
      q(1) = dble( jsirc(6))
      q(2) = dimag(jsirc(6))
      q(3) = dimag(jsirc(5))
      
      
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)
      vm2 = vmass**2

#ifdef HELAS_CHECK
      if ( abs(jsirc(5))+abs(jsirc(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : jsirc in jsiroc has zero momentum'
      endif
      if ( vwidth.eq.rZero .and. q2.eq.vm2 ) then
         write(stdo,*)
     &        ' helas-error : jsirc in jsiroc is on vmass pole'
         write(stdo,*)
     &        '             : q     = ',q(0),q(1),q(2),q(3)
         write(stdo,*)
     &        '             : abs(q)= ',sqrt(abs(q2))
         jsirc(1) = cZero
         jsirc(2) = cZero
         jsirc(3) = cZero
         jsirc(4) = cZero
         return
      endif
#endif

      gcc(1) = dconjg(gc(2))
      gcc(2) = dconjg(gc(2))

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

      fr(1) = dconjg(gcc(1))*(foc(3)*rc(0,3)+foc(4)*rc(0,4))
     &       +dconjg(gcc(2))*(foc(1)*rc(0,1)+foc(2)*rc(0,2))
      fr(2) = dconjg(gcc(1))*(foc(3)*rc(1,3)+foc(4)*rc(1,4))
     &       +dconjg(gcc(2))*(foc(1)*rc(1,1)+foc(2)*rc(1,2))
      fr(3) = dconjg(gcc(1))*(foc(3)*rc(2,3)+foc(4)*rc(2,4))
     &       +dconjg(gcc(2))*(foc(1)*rc(2,1)+foc(2)*rc(2,2))
      fr(4) = dconjg(gcc(1))*(foc(3)*rc(3,3)+foc(4)*rc(3,4))
     &       +dconjg(gcc(2))*(foc(1)*rc(3,1)+foc(2)*rc(3,2))

      sq1(1,1) = q(0)-q(3)
      sq1(1,2) = -(q(1)-cImag*q(2))
      sq1(2,1) = -(q(1)+cImag*q(2))
      sq1(2,2) = q(0)+q(3)

      sq2(1,1) = q(0)+q(3)
      sq2(1,2) = q(1)-cImag*q(2)
      sq2(2,1) = q(1)+cImag*q(2)
      sq2(2,2) = q(0)-q(3)

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

      fgsqr = sc(1)*( dconjg(gcc(1))*(gsqr(3)*foc(3)+gsqr(4)*foc(4))
     &               +dconjg(gcc(2))*(gsqr(1)*foc(1)+gsqr(2)*foc(2)))
      
      
      jsirc(1) = cZero
      jsirc(2) = cZero
      jsirc(3) = cZero
      jsirc(4) = cZero

      if ( vmass.ne.rZero ) then

       d = rOne/dcmplx( q2-vm2, vmass*vwidth )

       do i=1,4
	  jsirc(i) = d*(2*sc(1)*fr(i)-q(i-1)/vmass**2*fgsqr)
       enddo

      else

       d = dcmplx( rOne/q2, rZero )
	 do i=1,4
	  jsirc(i) = d*2*sc(1)*fr(i)
	 enddo         

      end if

      return
      end

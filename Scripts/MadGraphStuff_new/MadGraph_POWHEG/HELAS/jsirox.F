      subroutine jsirox(ri,fo,sc,gc,vmass,vwidth , jsiro)
c
c This subroutine computes an off-shell vector current from an external
c fermion, internal Rarita-Schwinger fermion, and a scalar. The vector
c boson propagator is given in Feynman gauge for a massless vector and
c in unitary gauge for a massive vector.
c
c input:
c       complex ri(18)         : flow-in  Rarita-Schwinger fermion   |ri>
c       complex fo(6)         : flow-out fermion                     <fo|
c       complex gc(2)          : coupling constants                  gfrvs
c       real    vmass          : mass  of OUTPUT vector v
c       real    vwidth         : width of OUTPUT vector v
c
c output:
c       complex jsiro(6)        : vector current          j^mu(<fo|s|ri>)
c   
c- by Y.Takaesu - 2009/11/01
c  
      implicit none
      double complex fo(6),ri(18),gc(2),sc(3),jsiro(6),c0,c1,c2,c3,cs,d
      double precision q(0:3),vmass,vwidth,q2,vm2
      double complex cm2 ! mass**2- I Gamma mass (Fabio)

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
      if ( abs(fo(1))+abs(fo(2))+abs(fo(3))+abs(fo(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : fo in jsirox is zero spinor'
      endif
      if ( abs(fo(5))+abs(fo(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fo in jsirox has zero momentum'
      endif
      if (abs(ri(1))+abs(ri(2))+abs(ri(3))+abs(ri(4))+abs(ri(5))
     & +abs(ri(6))+abs(ri(7))+abs(ri(8))abs(ri(9))+abs(ri(10))
     & +abs(ri(11))+abs(ri(12))+abs(ri(13))+abs(ri(14))+abs(ri(15))
     & +abs(ri(16)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : ri in jsirox is zero spinor'
      endif
      if ( abs(ri(17))+abs(ri(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : ri in jsirox has zero momentum'
      endif
	if ( sc(1).eq.cZero ) then
         write(stdo,*) ' helas-warn  : sc in jsirox is zero scalar'
      endif
      if ( abs(sc(2))+abs(sc(3)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : sc in jsirox has zero momentum'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*) ' helas-error : gc in jsirox is zero coupling'
      endif
      if ( vmass.lt.rZero ) then
         write(stdo,*) ' helas-error : vmass in jsirox is negative'
         write(stdo,*) '             : vmass = ',vmass
      endif
      if ( vwidth.lt.rZero ) then
         write(stdo,*) ' helas-error : vwidth in jsirox is negative'
         write(stdo,*) '             : vwidth = ',vwidth
      endif
#endif

      jsiro(5) = fo(5)-ri(17)+sc(2)
      jsiro(6) = fo(6)-ri(18)+sc(3)

      q(0) = dble( jsiro(5))
      q(1) = dble( jsiro(6))
      q(2) = dimag(jsiro(6))
      q(3) = dimag(jsiro(5))
      
      
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)
      vm2 = vmass**2

#ifdef HELAS_CHECK
      if ( abs(jsiro(5))+abs(jsiro(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : jsiro in jsirox has zero momentum'
      endif
      if ( vwidth.eq.rZero .and. q2.eq.vm2 ) then
         write(stdo,*)
     &        ' helas-error : jsiro in jsirox is on vmass pole'
         write(stdo,*)
     &        '             : q     = ',q(0),q(1),q(2),q(3)
         write(stdo,*)
     &        '             : abs(q)= ',sqrt(abs(q2))
         jsiro(1) = cZero
         jsiro(2) = cZero
         jsiro(3) = cZero
         jsiro(4) = cZero
         return
      endif
#endif

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

      fr(1) = dconjg(gc(1))*(fo(3)*rc(0,3)+fo(4)*rc(0,4))
     &       +dconjg(gc(2))*(fo(1)*rc(0,1)+fo(2)*rc(0,2))
      fr(2) = dconjg(gc(1))*(fo(3)*rc(1,3)+fo(4)*rc(1,4))
     &       +dconjg(gc(2))*(fo(1)*rc(1,1)+fo(2)*rc(1,2))
      fr(3) = dconjg(gc(1))*(fo(3)*rc(2,3)+fo(4)*rc(2,4))
     &       +dconjg(gc(2))*(fo(1)*rc(2,1)+fo(2)*rc(2,2))
      fr(4) = dconjg(gc(1))*(fo(3)*rc(3,3)+fo(4)*rc(3,4))
     &       +dconjg(gc(2))*(fo(1)*rc(3,1)+fo(2)*rc(3,2))

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

      fgsqr = sc(1)*( dconjg(gc(1))*(gsqr(3)*fo(3)+gsqr(4)*fo(4))
     &               +dconjg(gc(2))*(gsqr(1)*fo(1)+gsqr(2)*fo(2)))
      
      
      jsiro(1) = cZero
      jsiro(2) = cZero
      jsiro(3) = cZero
      jsiro(4) = cZero

      if ( vmass.ne.rZero ) then

       d = rOne/dcmplx( q2-vm2, vmass*vwidth )

       do i=1,4
	  jsiro(i) = d*(2*sc(1)*fr(i)-q(i-1)/vmass**2*fgsqr)
       enddo

      else

       d = dcmplx( rOne/q2, rZero )
	 do i=1,4
	  jsiro(i) = d*2*sc(1)*fr(i)
	 enddo         

      end if

      return
      end

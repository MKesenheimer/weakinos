      subroutine jsiorx(fi,ro,sc,gc,vmass,vwidth , jsior)
c
c This subroutine computes an off-shell vector current from an external
c Rarita-Schwinger fermion, internal fermion, and a scalar. The vector
c boson propagator is given in Feynman gauge for a massless vector and
c in unitary gauge for a massive vector.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex ro(18)         : flow-out Rarita-Schwinger fermion  <ro|
c       complex gc(2)          : coupling constants                 gfrvs
c       real    vmass          : mass  of OUTPUT vector v
c       real    vwidth         : width of OUTPUT vector v
c
c output:
c       complex jsior(6)        : vector current          j^mu(<ro|s|fi>)
c   
c- by Y.Takaesu - 2009/11/01
c  
      implicit none
      double complex fi(6),ro(18),gc(2),sc(3),jsior(6),c0,c1,c2,c3,cs,d
      double precision q(0:3),vmass,vwidth,q2,vm2
      double complex cm2 ! mass**2- I Gamma mass (Fabio)

      double complex rsq(0:3,4),rsqg(4),sq1(2,2),sq2(2,2),rsqgf,rf(4)
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
      if ( abs(fi(1))+abs(fi(2))+abs(fi(3))+abs(fi(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : fi in jsiorx is zero spinor'
      endif
      if ( abs(fi(5))+abs(fi(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fi in jsiorx has zero momentum'
      endif
      if (abs(ro(1))+abs(ro(2))+abs(ro(3))+abs(ro(4))+abs(ro(5))
     & +abs(ro(6))+abs(ro(7))+abs(ro(8))abs(ro(9))+abs(ro(10))
     & +abs(ro(11))+abs(ro(12))+abs(ro(13))+abs(ro(14))+abs(ro(15))
     & +abs(ro(16)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : ro in jsiorx is zero spinor'
      endif
      if ( abs(ro(17))+abs(ro(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fo in jsiorx has zero momentum'
      endif
	if ( sc(1).eq.cZero ) then
         write(stdo,*) ' helas-warn  : sc in jsiorx is zero scalar'
      endif
      if ( abs(sc(2))+abs(sc(3)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : sc in jsiorx has zero momentum'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*) ' helas-error : gc in jsiorx is zero coupling'
      endif
      if ( vmass.lt.rZero ) then
         write(stdo,*) ' helas-error : vmass in jsiorx is negative'
         write(stdo,*) '             : vmass = ',vmass
      endif
      if ( vwidth.lt.rZero ) then
         write(stdo,*) ' helas-error : vwidth in jsiorxx is negative'
         write(stdo,*) '             : vwidth = ',vwidth
      endif
#endif

      jsior(5) = ro(17)-fi(5)+sc(2)
      jsior(6) = ro(18)-fi(6)+sc(3)

      q(0) = dble( jsior(5))
      q(1) = dble( jsior(6))
      q(2) = dimag(jsior(6))
      q(3) = dimag(jsior(5))
      
      
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)
      vm2 = vmass**2

#ifdef HELAS_CHECK
      if ( abs(jsior(5))+abs(jsior(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : jsior in jsiorx has zero momentum'
      endif
      if ( vwidth.eq.rZero .and. q2.eq.vm2 ) then
         write(stdo,*)
     &        ' helas-error : jsior in jsiorx is on vmass pole'
         write(stdo,*)
     &        '             : q     = ',q(0),q(1),q(2),q(3)
         write(stdo,*)
     &        '             : abs(q)= ',sqrt(abs(q2))
         jsior(1) = cZero
         jsior(2) = cZero
         jsior(3) = cZero
         jsior(4) = cZero
         return
      endif
#endif

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

      rf(1) = gc(1)*(rc(0,1)*fi(1)+rc(0,2)*fi(2))
     &       +gc(2)*(rc(0,3)*fi(3)+rc(0,4)*fi(4))
      rf(2) = gc(1)*(rc(1,1)*fi(1)+rc(1,2)*fi(2))
     &       +gc(2)*(rc(1,3)*fi(3)+rc(1,4)*fi(4))
      rf(3) = gc(1)*(rc(2,1)*fi(1)+rc(2,2)*fi(2))
     &       +gc(2)*(rc(2,3)*fi(3)+rc(2,4)*fi(4))
      rf(4) = gc(1)*(rc(3,1)*fi(1)+rc(3,2)*fi(2))
     &       +gc(2)*(rc(3,3)*fi(3)+rc(3,4)*fi(4))

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

      rsqgf = sc(1)*(gc(1)*(rsqg(1)*fi(1)+rsqg(2)*fi(2))
     &              +gc(2)*(rsqg(3)*fi(3)+rsqg(4)*fi(4)))
      
      jsior(1) = cZero
      jsior(2) = cZero
      jsior(3) = cZero
      jsior(4) = cZero

      if ( vmass.ne.rZero ) then

       d = rOne/dcmplx( q2-vm2, vmass*vwidth )

       do i=1,4
	  jsior(i) = d*(2*rf(i)*sc(1)-q(i-1)/vmass**2*rsqgf)
	 enddo

      else

       d = dcmplx( rOne/q2, rZero )
	 do i=1,4
	  jsior(i) = d*2*rf(i)*sc(1)
	 enddo         

      end if

      return
      end

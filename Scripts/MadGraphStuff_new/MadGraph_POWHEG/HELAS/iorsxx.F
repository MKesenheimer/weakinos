      subroutine iorsxx(fi,ro,sc,gc , vertex)
c
c This subroutine computes an amplitude of the fermion-(Rarita-Schwinger
c fermion)-scalar coupling.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex ro(18)         : flow-out RS fermion               v<ro|
c       complex sc(3)          : input    scalar                    s(q)
c       complex gc(2)          : coupling constants                 GFRS
c
c output:
c       complex vertex         : amplitude                   <ro|sqv|fi>
c  
c- by Y.Takaesu - 2009/10/30
c     
      implicit none
      double complex fi(6),ro(18),sc(3),gc(2),vertex,rc(0:3,4)
      double complex rsq(0:3,4),rsqg(4),sq1(2,2),sq2(2,2),cImag
      parameter (cImag = ( 0.0d0, 1.0d0 ))
      double precision q(0:3)
      integer i

#ifdef HELAS_CHECK
      double precision p0,p1,p2,p3,q0,q1,q2,q3,pm
      double precision epsi
      parameter( epsi = 4.0d-5 )
      double precision rZero
      parameter( rZero = 0.0d0 )
      double complex cZero
      parameter(  cZero = ( 0.0d0, 0.0d0 ) )
      integer stdo
      parameter( stdo = 6 )
#endif
c
#ifdef HELAS_CHECK
      p0 = -dble( fi(5))
      p1 = -dble( fi(6))
      p2 = -dimag(fi(6))
      p3 = -dimag(fi(5))
      q0 = dble( ro(17))
      q1 = dble( ro(18))
      q2 = dimag(ro(18))
      q3 = dimag(ro(17))
      r(0) = dble( sc(2))
      r(1) = dble( sc(3))
      r(2) = dimag(sc(3))
      r(3) = dimag(sc(2))
      if ( abs(fi(1))+abs(fi(2))+abs(fi(3))+abs(fi(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : fi in iorsxx is zero spinor'
      endif
      if ( abs(fi(5))+abs(fi(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fi in iorsxx has zero momentum'
      endif
      if (abs(ro(1))+abs(ro(2))+abs(ro(3))+abs(ro(4))+abs(ro(5))
     & +abs(ro(6))+abs(ro(7))+abs(ro(8))abs(ro(9))+abs(ro(10))
     & +abs(ro(11))+abs(ro(12))+abs(ro(13))+abs(ro(14))+abs(ro(15))
     & +abs(ro(16)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : ro in iorsxx is zero spinor'
      endif
      if ( abs(ro(17))+abs(ro(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : ro in iorsxx has zero momentum'
      endif
      if ( sc(1).eq.cZero ) then
         write(stdo,*) ' helas-warn  : sc in iorsxx is zero scalar'
      endif
      if ( abs(sc(2))+abs(sc(3)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : sc in iorsxx has zero momentum'
      endif
      pm = max( abs(p0),abs(q0),abs(r0),abs(p1),abs(q1),abs(r1),
     &          abs(p2),abs(q2),abs(r2),abs(p3),abs(q3),abs(r3) )
      if ( abs(-fi(5)+ro(17)+sc(2))+abs(-fi(6)+ro(18)+sc(3))
     &                                               .ge.pm*epsi) then
         write(stdo,*)
     &        ' helas-error : fi,ro,sc in iorsxx'
         write(stdo,*)
     &        '             :          have not balanced momenta'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*)
     &        ' helas-error : gc in iorsxx is zero coupling'
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

      vertex = sc(1)*( gc(1)*(rsqg(1)*fi(1)+rsqg(2)*fi(2))
     &                +gc(2)*(rsqg(3)*fi(3)+rsqg(4)*fi(4)))
      
      return
      end

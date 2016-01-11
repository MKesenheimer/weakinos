      subroutine iorvsc(fic,roc,vc,sc,gc , vertex)
c
c This subroutine computes an amplitude of the antifermion-(antiRarita-Schwinger fermion)
c -vector-scalar coupling.
c
c input:
c       complex fic(6)         : flow-in  antifermion               |fic>
c       complex roc(18)        : flow-out antiRSfermion             <roc|
c	complex vc(6)          : input    vector                      v
c       complex sc(3)          : input    scalar                      s
c       complex gc(2)          : coupling constants                 gvhrf
c
c output:
c       complex vertex         : amplitude                     <roc|vs|fic>
c  
c- by Yoshitaro Takaesu - 2011/03/26
c     
      implicit none
      double complex fic(6),roc(18),vc(6),sc(3),gc(2),vertex,rc(0:3,4)
      double complex rsq(0:3,4),rsqg(4),sq1(2,2),sq2(2,2),cImag,q(0:3)
      double complex gcc(2)
      parameter (cImag = ( 0.0d0, 1.0d0 ))
      integer i

#ifdef HELAS_CHECK
      double precision p0,p1,p2,p3,q0,q1,q2,q3,r0,r1,r2,r3,pm
      double precision epsi,o0,o1,o2,o3
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
      p0 = -dble( fic(5))
      p1 = -dble( fic(6))
      p2 = -dimag(fic(6))
      p3 = -dimag(fic(5))
      q0 = dble( roc(17))
      q1 = dble( roc(18))
      q2 = dimag(roc(18))
      q3 = dimag(roc(17))
      r0 = dble( sc(2))
      r1 = dble( sc(3))
      r2 = dimag(sc(3))
      r3 = dimag(sc(2))
      o0 = dble( vc(5))
      o1 = dble( vc(6))
      o2 = dimag(vc(6))
      o3 = dimag(vc(5))

      if ( abs(fic(1))+abs(fic(2))+abs(fic(3))+abs(fic(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : fic in iorvsc is zero spinor'
      endif
      if ( abs(fic(5))+abs(fic(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fic in iorvsc has zero momentum'
      endif
      if (abs(roc(1))+abs(roc(2))+abs(roc(3))+abs(roc(4))+abs(roc(5))
     & +abs(roc(6))+abs(roc(7))+abs(roc(8))abs(roc(9))+abs(roc(10))
     & +abs(roc(11))+abs(roc(12))+abs(roc(13))+abs(roc(14))+abs(roc(15))
     & +abs(roc(16)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : roc in iorvsc is zero spinor'
      endif
      if ( abs(roc(17))+abs(roc(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : roc in iorvsc has zero momentum'
      endif
      if ( sc(1).eq.cZero ) then
         write(stdo,*) ' helas-warn  : sc in iorvsc is zero scalar'
      endif
      if ( abs(sc(2))+abs(sc(3)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : sc in iorvsc has zero momentum'
      endif
	if ( abs(vc(5))+abs(vc(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : vc in iorvsc has zero momentum'
      endif
      pm = max( abs(p0),abs(q0),abs(r0),abs(o0),abs(p1),abs(q1),
     &          abs(r1),abs(o1),abs(p2),abs(q2),abs(r2),abs(o2),
     &          abs(p3),abs(q3),abs(r3),abs(o3))
      if ( abs(-fic(5)+roc(17)+sc(2)+vc(5))+abs(-fic(6)+roc(18)+sc(3)
     &                                     +vc(6)).ge.pm*epsi) then
         write(stdo,*)
     &        ' helas-error : fic,roc,vc,sc in iorvsc'
         write(stdo,*)
     &        '             :          have not balanced momenta'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*)
     &        ' helas-error : gc in iorvsc is zero coupling'
      endif
#endif

      gcc(1) = dconjg(gc(2))
      gcc(2) = dconjg(gc(1))

      q(0) = vc(1)
      q(1) = vc(2)
      q(2) = vc(3)
      q(3) = vc(4)

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

      vertex = sc(1)*(gcc(1)*(rsqg(1)*fic(1)+rsqg(2)*fic(2))
     &              + gcc(2)*(rsqg(3)*fic(3)+rsqg(4)*fic(4)))
      
      return
      end

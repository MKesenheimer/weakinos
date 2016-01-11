      subroutine irosxx(ri,fo,sc,gc , vertex)
c
c This subroutine computes an amplitude of the fermion-(Rarita-Schwinger
c fermion)-scalar coupling.
c
c input:
c       complex ri(18)         : flow-in  RS fermion               |ri>v
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex sc(3)          : input    scalar                    s(q)
c       complex gc(2)          : coupling constants                 GFRS  
c
c output:
c       complex vertex         : amplitude                   <fo|svq|ri>
c  
c- by Y.Takaesu - 2009/06/25
c    
      implicit none
      double complex fo(6),ri(18),sc(3),gc(2),vertex,rc(0:3,4)
      double complex sqr(0:3,4),gsqr(4),sq1(2,2),sq2(2,2),cImag
      parameter (cImag = ( 0.0d0, 1.0d0 ))
      double precision r(0:3)
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
      p0 = dble( fo(5))
      p1 = dble( fo(6))
      p2 = dimag(fo(6))
      p3 = dimag(fo(5))
      q0 = -dble( ri(17))
      q1 = -dble( ri(18))
      q2 = -dimag(ri(18))
      q3 = -dimag(ri(17))
      r(0) = dble( sc(2))
      r(1) = dble( sc(3))
      r(2) = dimag(sc(3))
      r(3) = dimag(sc(2))
      if ( abs(fo(1))+abs(fo(2))+abs(fo(3))+abs(fo(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : fo in irosxx is zero spinor'
      endif
      if ( abs(fo(5))+abs(fo(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fo in irosxx has zero momentum'
      endif
      if (abs(ri(1))+abs(ri(2))+abs(ri(3))+abs(ri(4))+abs(ri(5))
     & +abs(ri(6))+abs(ri(7))+abs(ri(8))abs(ri(9))+abs(ri(10))
     & +abs(ri(11))+abs(ri(12))+abs(ri(13))+abs(ri(14))+abs(ri(15))
     & +abs(ri(16)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : ri in irosxx is zero spinor'
      endif
      if ( abs(ri(17))+abs(ri(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : ri in irosxx has zero momentum'
      endif
      if ( sc(1).eq.cZero ) then
         write(stdo,*) ' helas-warn  : sc in irosxx is zero scalar'
      endif
      if ( abs(sc(2))+abs(sc(3)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : sc in irosxx has zero momentum'
      endif
      pm = max( abs(p0),abs(q0),abs(r0),abs(p1),abs(q1),abs(r1),
     &          abs(p2),abs(q2),abs(r2),abs(p3),abs(q3),abs(r3) )
      if ( abs(fo(5)-ri(17)+sc(2))+abs(fo(6)-ri(18)+sc(3))
     &                                               .ge.pm*epsi) then
         write(stdo,*)
     &        ' helas-error : fo,ri,sc in irosxx'
         write(stdo,*)
     &        '             :          have not balanced momenta'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*)
     &        ' helas-error : gc in irosxx is zero coupling'
      endif
#endif

      r(0) = dble( sc(2))
      r(1) = dble( sc(3))
      r(2) = dimag(sc(3))
      r(3) = dimag(sc(2))

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

      sq1(1,1) = r(0)-r(3)
      sq1(1,2) = -(r(1)-cImag*r(2))
      sq1(2,1) = -(r(1)+cImag*r(2))
      sq1(2,2) = r(0)+r(3)

      sq2(1,1) = r(0)+r(3)
      sq2(1,2) = r(1)-cImag*r(2)
      sq2(2,1) = r(1)+cImag*r(2)
      sq2(2,2) = r(0)-r(3)

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

      vertex = sc(1)*( dconjg(gc(1))*(gsqr(3)*fo(3)+gsqr(4)*fo(4))
     &                +dconjg(gc(2))*(gsqr(1)*fo(1)+gsqr(2)*fo(2)))      
      vertex = -vertex

      return
      end

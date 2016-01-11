      subroutine irovsx(ri,fo,vc,sc,gc , vertex)
c
c This subroutine computes an amplitude of the spin-2/3 Schwinger fermion
c -scalar-vector coupling.
c
c input:
c       complex ri(18)         : flow-in  spin-2/3 fermion          |ri>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex vc(6)          : input    vector                       v' 
c       complex sc(3)          : input    scalar                       s
c       complex gc(2)          : coupling constants                girovs
c
c output:
c       complex vertex         : amplitude                   <fo|v'|ri>s
c     
c- by Y.Takaesu - 2008/10/28
c
      implicit none
      double complex ri(18),fo(6),vc(6),sc(3),gc(2),vertex,q(0:3)
      double complex rc(0:3,4),svr(0:3,4),gsvr(4),sv1(2,2),sv2(2,2)
      double precision ps(4)

      double precision rZero, rOne
      parameter( rZero = 0.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )
      integer i

#ifdef HELAS_CHECK
      double precision p0,p1,p2,p3,r0,r1,r2,r3,pm
      double precision s0,s1,s2,s3
      double precision epsi
      parameter( epsi = 4.0d-5 )
      integer stdo
      parameter( stdo = 6 )
#endif

#ifdef HELAS_CHECK
      p0 = -dble( fo(5))
      p1 = -dble( fo(6))
      p2 = -dimag(fo(6))
      p3 = -dimag(fo(5))
      r0 = dble( ri(17))
      r1 = dble( ri(18))
      r2 = dimag(ri(18))
      r3 = dimag(ri(17))
      s0 = dble( sc(2))
      s1 = dble( sc(3))
      s2 = dimag(sc(3))
      s3 = dimag(sc(2))
      if ( abs(fo(1))+abs(fo(2))+abs(fo(3))+abs(fo(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : fo in irovsx is zero spinor'
      endif
      if ( abs(fo(5))+abs(fo(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fo in irovsx has zero momentum'
      endif
      if ( abs(ri(17))+abs(ri(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : ri in irovsx has zero momentum'
      endif
      if ( sc(1).eq.cZero ) then
         write(stdo,*) ' helas-warn  : sc in irovsx is zero scalar'
      endif
      if ( abs(sc(2))+abs(sc(3)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : sc in irovsx has zero momentum'
      endif
      pm = max( abs(p0),abs(s0),abs(r0),abs(p1),abs(s1),abs(r1),
     &          abs(p2),abs(s2),abs(r2),abs(p3),abs(s3),abs(r3) )
      if ( abs(-fo(5)+ri(17)+sc(2))+abs(-fo(6)+ri(18)+sc(3))
     &                                              .ge.pm*epsi ) then
         write(stdo,*)
     &        ' helas-error : fo,ri,sc in irovsx'
         write(stdo,*)
     &        '                        have not balanced momenta'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*) ' helas-error : gc in irovsx is zero coupling'
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

      q(0) = vc(1)
      q(1) = vc(2)
      q(2) = vc(3)
      q(3) = vc(4)

      sv1(1,1) = q(0)-q(3)
      sv1(1,2) = -(q(1)-cImag*q(2))
      sv1(2,1) = -(q(1)+cImag*q(2))
      sv1(2,2) = q(0)+q(3)

      sv2(1,1) = q(0)+q(3)
      sv2(1,2) = q(1)-cImag*q(2)
      sv2(2,1) = q(1)+cImag*q(2)
      sv2(2,2) = q(0)-q(3)

      do i=0,3
        svr(i,1) = rc(i,3)*sv1(1,1)+rc(i,4)*sv1(1,2)
        svr(i,2) = rc(i,3)*sv1(2,1)+rc(i,4)*sv1(2,2)
        svr(i,3) = rc(i,1)*sv2(1,1)+rc(i,2)*sv2(1,2)
        svr(i,4) = rc(i,1)*sv2(2,1)+rc(i,2)*sv2(2,2)
      enddo

      gsvr(1) = svr(0,3)-svr(1,4)+cImag*svr(2,4)-svr(3,3)
      gsvr(2) = svr(0,4)-svr(1,3)-cImag*svr(2,3)+svr(3,4)
      gsvr(3) = svr(0,1)+svr(1,2)-cImag*svr(2,2)+svr(3,1)
      gsvr(4) = svr(0,2)+svr(1,1)+cImag*svr(2,1)-svr(3,2)

      vertex=sc(1)*(dconjg(gc(1))*(fo(3)*gsvr(3)+fo(4)*gsvr(4))
     &            +dconjg(gc(2))*(fo(1)*gsvr(1)+fo(2)*gsvr(2)))

      return
      end

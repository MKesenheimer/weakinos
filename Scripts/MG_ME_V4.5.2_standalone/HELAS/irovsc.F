      subroutine irovsc(ric,foc,vc,sc,gc , vertex)
c
c This subroutine computes an amplitude of the antifermion-(antiRatita-Schwinger fermion)
c -scalar-vector coupling.
c
c input:
c       complex ric(18)        : flow-in  antiRSfermion            |ric>
c       complex foc(6)         : flow-out antifermion              <foc|
c       complex vc(6)          : input    vector                       v' 
c       complex sc(3)          : input    scalar                       s
c       complex gc(2)          : coupling constants                girovs
c
c output:
c       complex vertex         : amplitude                   <foc|v'|ric>s
c     
c- by Yoshitaro Takaesu - 2011/03/25
c
      implicit none
      double complex ric(18),foc(6),vc(6),sc(3),gc(2),vertex,q(0:3)
      double complex rc(0:3,4),svr(0:3,4),gsvr(4),sv1(2,2),sv2(2,2)
      double complex gcc(2)
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
      p0 = -dble( foc(5))
      p1 = -dble( foc(6))
      p2 = -dimag(foc(6))
      p3 = -dimag(foc(5))
      r0 = dble( ric(17))
      r1 = dble( ric(18))
      r2 = dimag(ric(18))
      r3 = dimag(ric(17))
      s0 = dble( sc(2))
      s1 = dble( sc(3))
      s2 = dimag(sc(3))
      s3 = dimag(sc(2))
      if ( abs(foc(1))+abs(foc(2))+abs(foc(3))+abs(foc(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : foc in irovsc is zero spinor'
      endif
      if ( abs(foc(5))+abs(foc(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : foc in irovsc has zero momentum'
      endif
      if ( abs(ric(17))+abs(ric(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : ric in irovsc has zero momentum'
      endif
      if ( sc(1).eq.cZero ) then
         write(stdo,*) ' helas-warn  : sc in irovsc is zero scalar'
      endif
      if ( abs(sc(2))+abs(sc(3)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : sc in irovsc has zero momentum'
      endif
      pm = max( abs(p0),abs(s0),abs(r0),abs(p1),abs(s1),abs(r1),
     &          abs(p2),abs(s2),abs(r2),abs(p3),abs(s3),abs(r3) )
      if ( abs(-foc(5)+ric(17)+sc(2))+abs(-foc(6)+ric(18)+sc(3))
     &                                              .ge.pm*epsi ) then
         write(stdo,*)
     &        ' helas-error : foc,ric,sc in irovsc'
         write(stdo,*)
     &        '                        have not balanced momenta'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*) ' helas-error : gc in irovsc is zero coupling'
      endif
#endif

      gcc(1) = dconjg(gc(2))
      gcc(2) = dconjg(gc(1)) 

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

      vertex=sc(1)*(dconjg(gcc(1))*(foc(3)*gsvr(3)+foc(4)*gsvr(4))
     &            +dconjg(gcc(2))*(foc(1)*gsvr(1)+foc(2)*gsvr(2)))

      return
      end

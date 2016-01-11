      subroutine iroscx(ric,foc,sc,gc , vertex)
c
c This subroutine computes an amplitude of the antifermion-(antiRarita-Schwinger
c fermion)-scalar coupling.
c
c input:
c       complex ric(18)        : flow-in  antiRSfermion            |ric>v
c       complex foc(6)         : flow-out antifermion               <foc|
c       complex sc(3)          : input    scalar                    s(q)
c       complex gc(2)          : coupling constants                 GFRS  
c
c output:
c       complex vertex         : amplitude                   <foc|svq|ric>
c  
c- by Yoshitaro Takaesu - 2011/03/25
c    
      implicit none
      double complex foc(6),ric(18),sc(3),gc(2),gcc(2),vertex,rc(0:3,4)
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
      p0 = dble( foc(5))
      p1 = dble( foc(6))
      p2 = dimag(foc(6))
      p3 = dimag(foc(5))
      q0 = -dble( ric(17))
      q1 = -dble( ric(18))
      q2 = -dimag(ric(18))
      q3 = -dimag(ric(17))
      r(0) = dble( sc(2))
      r(1) = dble( sc(3))
      r(2) = dimag(sc(3))
      r(3) = dimag(sc(2))
      if ( abs(foc(1))+abs(foc(2))+abs(foc(3))+abs(foc(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : foc in iroscx is zero spinor'
      endif
      if ( abs(foc(5))+abs(foc(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : foc in iroscx has zero momentum'
      endif
      if (abs(ric(1))+abs(ric(2))+abs(ric(3))+abs(ric(4))+abs(ric(5))
     & +abs(ric(6))+abs(ric(7))+abs(ric(8))abs(ric(9))+abs(ric(10))
     & +abs(ric(11))+abs(ric(12))+abs(ric(13))+abs(ric(14))+abs(ric(15))
     & +abs(ric(16)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : ric in iroscx is zero spinor'
      endif
      if ( abs(ric(17))+abs(ric(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : ric in iroscx has zero momentum'
      endif
      if ( sc(1).eq.cZero ) then
         write(stdo,*) ' helas-warn  : sc in iroscx is zero scalar'
      endif
      if ( abs(sc(2))+abs(sc(3)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : sc in iroscx has zero momentum'
      endif
      pm = max( abs(p0),abs(q0),abs(r0),abs(p1),abs(q1),abs(r1),
     &          abs(p2),abs(q2),abs(r2),abs(p3),abs(q3),abs(r3) )
      if ( abs(foc(5)-ric(17)+sc(2))+abs(foc(6)-ric(18)+sc(3))
     &                                               .ge.pm*epsi) then
         write(stdo,*)
     &        ' helas-error : foc,ric,sc in iroscx'
         write(stdo,*)
     &        '             :          have not balanced momenta'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*)
     &        ' helas-error : gc in iroscx is zero coupling'
      endif
#endif

      gcc(1) = dconjg(gc(2))
      gcc(2) = dconjg(gc(1)) 

      r(0) = dble( sc(2))
      r(1) = dble( sc(3))
      r(2) = dimag(sc(3))
      r(3) = dimag(sc(2))

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

      vertex = sc(1)*( dconjg(gcc(1))*(gsqr(3)*foc(3)+gsqr(4)*foc(4))
     &                +dconjg(gcc(2))*(gsqr(1)*foc(1)+gsqr(2)*foc(2)))      
 
      return
      end

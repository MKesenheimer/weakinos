      subroutine irovvc(ric,foc,vc1,vc2,gc , vertex)
c
c This subroutine computes an amplitude of the antifermion-(antiRatrita-Schwinger
c fermion)-vector-vector coupling.
c
c input:
c       complex ric(18)        : flow-in  antiRSfermion            |ric>
c       complex foc(6)         : flow-out antifermion              <foc|
c       complex vc1(6)         : input    vector                    v1
c       complex vc2(6)         : input    vector                    v2
c       complex gc(2)          : coupling constants                 gfrvv
c
c output:
c       complex vertex         : amplitude               <foc|[v1,v2]|ric>
c     
c- by Yoshitaro Takaesu - 2011/03/26
c
      implicit none
      double complex ric(18),foc(6),vc1(6),vc2(6),gc(2),vertex
      double complex rc(4,4),spvl(2,2),svcl(2,2),sssl(2,2),
     &                       spvr(2,2),svcr(2,2),sssr(2,2),gcc(2)
      double precision pv(4)

      double precision rZero, rOne
      parameter( rZero = 0.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )

#ifdef HELAS_CHECK
      double precision p0,p1,p2,p3,q0,q1,q2,q3,r0,r1,r2,r3,pm
      double precision epsi,s0,s1,s2,s3
      parameter( epsi = 4.0d-5 )
      integer stdo
      parameter( stdo = 6 )
#endif

#ifdef HELAS_CHECK
      p0 = -dble( foc(5))
      p1 = -dble( foc(6))
      p2 = -dimag(foc(6))
      p3 = -dimag(foc(5))
      q0 = dble( ric(17))
      q1 = dble( ric(18))
      q2 = dimag(ric(18))
      q3 = dimag(ric(17))
      r0 = dble( vc1(5))
      r1 = dble( vc1(6))
      r2 = dimag(vc1(6))
      r3 = dimag(vc1(5))
      s0 = dble( vc2(5))
      s1 = dble( vc2(6))
      s2 = dimag(vc2(6))
      s3 = dimag(vc2(5))
      if ( abs(foc(1))+abs(foc(2))+abs(foc(3))+abs(foc(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : foc in irovvc is zero spinor'
      endif
      if ( abs(foc(5))+abs(foc(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : foc in irovvc has zero momentum'
      endif
      if (abs(ric(1))+abs(ric(2))+abs(ric(3))+abs(ric(4))+abs(ric(5))
     & +abs(ric(6))+abs(ric(7))+abs(ric(8))abs(ric(9))+abs(ric(10))
     & +abs(ric(11))+abs(ric(12))+abs(ric(13))+abs(ric(14))+abs(ric(15))
     & +abs(ric(16)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : ric in irovvc is zero spinor'
      endif
      if ( abs(ric(17))+abs(ric(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : ro in irovvc has zero momentum'
      endif
      if ( abs(vc1(1))+abs(vc1(2))+abs(vc1(3))+abs(vc1(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : vc1 in irovvc is zero vector'
      endif
      if ( abs(vc1(5))+abs(vc1(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : vc1 in irovvc has zero momentum'
      endif
      if ( abs(vc2(1))+abs(vc2(2))+abs(vc2(3))+abs(vc2(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : vc2 in irovvc is zero vector'
      endif
      if ( abs(vc2(5))+abs(vc2(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : vc2 in irovvc has zero momentum'
      endif
      pm = max( abs(p0),abs(q0),abs(r0),abs(s0),abs(p1),abs(q1),abs(r1),
     &          abs(s1),abs(p2),abs(q2),abs(r2),abs(s2),abs(p3),abs(q3),
     &          abs(r3),abs(s3) )
      if ( abs(-ric(17)+foc(5)+vc1(5)*vc2(5))+abs(-ric(18)+foc(6)+vc(6)+vc2(6))
     &                                              .ge.pm*epsi ) then
         write(stdo,*)
     &        ' helas-error : ric,foc,vc1,vc2 in irovvc'
         write(stdo,*)
     &        '                        have not balanced momenta'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*) ' helas-error : gc in irovvc is zero coupling'
      endif
#endif

      gcc(1) = dconjg(gc(2))
      gcc(2) = dconjg(gc(1))

      rc(1,1) = ric(1)
      rc(1,2) = ric(2)
      rc(1,3) = ric(3)
      rc(1,4) = ric(4)
      rc(2,1) = ric(5)
      rc(2,2) = ric(6)
      rc(2,3) = ric(7)
      rc(2,4) = ric(8)
      rc(3,1) = ric(9)
      rc(3,2) = ric(10)
      rc(3,3) = ric(11)
      rc(3,4) = ric(12)
      rc(4,1) = ric(13)
      rc(4,2) = ric(14)
      rc(4,3) = ric(15)
      rc(4,4) = ric(16)

      vertex = cZero

      if ( gcc(2).ne.cZero ) then

        spvl(1,1) = vc1(1)+vc1(4)
        spvl(1,2) = vc1(2)-cImag*vc1(3)
        spvl(2,1) = vc1(2)+cImag*vc1(3)
        spvl(2,2) = vc1(1)-vc1(4)
        
        svcl(1,1) = vc2(1)-vc2(4)
        svcl(1,2) = -(vc2(2)-cImag*vc2(3))
        svcl(2,1) = -(vc2(2)+cImag*vc2(3))
        svcl(2,2) = vc2(1)+vc2(4)

        sssl(1,1) = spvl(1,1)*svcl(1,1)+spvl(1,2)*svcl(2,1)
        sssl(1,2) = spvl(1,1)*svcl(1,2)+spvl(1,2)*svcl(2,2) 
        sssl(2,1) = spvl(2,1)*svcl(1,1)+spvl(2,2)*svcl(2,1) 
        sssl(2,2) = spvl(2,1)*svcl(1,2)+spvl(2,2)*svcl(2,2)     

      vertex = dconjg(gcc(2)) 
     &          *( foc(1) *( sssl(1,1)*rc(1,3)+sssl(1,2)*rc(1,4) )
     &            +foc(2) *( sssl(2,1)*rc(1,3)+sssl(2,2)*rc(1,4) )
     &            -foc(1) *( sssl(2,1)*rc(2,3)+sssl(2,2)*rc(2,4) )
     &            -foc(2) *( sssl(1,1)*rc(2,3)+sssl(1,2)*rc(2,4) )
     &            -foc(1) *(-sssl(2,1)*rc(3,3)-sssl(2,2)*rc(3,4) )*cImag
     &            -foc(2) *( sssl(1,1)*rc(3,3)+sssl(1,2)*rc(3,4) )*cImag
     &            -foc(1) *( sssl(1,1)*rc(4,3)+sssl(1,2)*rc(4,4) )
     &            -foc(2) *(-sssl(2,1)*rc(4,3)-sssl(2,2)*rc(4,4) )
     &           )


        spvl(1,1) = vc2(1)+vc2(4)
        spvl(1,2) = vc2(2)-cImag*vc2(3)
        spvl(2,1) = vc2(2)+cImag*vc2(3)
        spvl(2,2) = vc2(1)-vc2(4)
        
        svcl(1,1) = vc1(1)-vc1(4)
        svcl(1,2) = -(vc1(2)-cImag*vc1(3))
        svcl(2,1) = -(vc1(2)+cImag*vc1(3))
        svcl(2,2) = vc1(1)+vc1(4)

        sssl(1,1) = spvl(1,1)*svcl(1,1)+spvl(1,2)*svcl(2,1)
        sssl(1,2) = spvl(1,1)*svcl(1,2)+spvl(1,2)*svcl(2,2) 
        sssl(2,1) = spvl(2,1)*svcl(1,1)+spvl(2,2)*svcl(2,1) 
        sssl(2,2) = spvl(2,1)*svcl(1,2)+spvl(2,2)*svcl(2,2)

      vertex = vertex
     &        -dconjg(gcc(2)) 
     &          *( foc(1) *( sssl(1,1)*rc(1,3)+sssl(1,2)*rc(1,4) )
     &            +foc(2) *( sssl(2,1)*rc(1,3)+sssl(2,2)*rc(1,4) )
     &            -foc(1) *( sssl(2,1)*rc(2,3)+sssl(2,2)*rc(2,4) )
     &            -foc(2) *( sssl(1,1)*rc(2,3)+sssl(1,2)*rc(2,4) )
     &            -foc(1) *(-sssl(2,1)*rc(3,3)-sssl(2,2)*rc(3,4) )*cImag
     &            -foc(2) *( sssl(1,1)*rc(3,3)+sssl(1,2)*rc(3,4) )*cImag
     &            -foc(1) *( sssl(1,1)*rc(4,3)+sssl(1,2)*rc(4,4) )
     &            -foc(2) *(-sssl(2,1)*rc(4,3)-sssl(2,2)*rc(4,4) )
     &           )
      
      endif

      if ( gcc(1).ne.cZero ) then

        spvr(1,1) = vc1(1)-vc1(4)
        spvr(1,2) = -(vc1(2)-cImag*vc1(3))
        spvr(2,1) = -(vc1(2)+cImag*vc1(3))
        spvr(2,2) = vc1(1)+vc1(4)
        
        svcr(1,1) = vc2(1)+vc2(4)
        svcr(1,2) = vc2(2)-cImag*vc2(3)
        svcr(2,1) = vc2(2)+cImag*vc2(3)
        svcr(2,2) = vc2(1)-vc2(4)

        sssr(1,1) = spvr(1,1)*svcr(1,1)+spvr(1,2)*svcr(2,1)
        sssr(1,2) = spvr(1,1)*svcr(1,2)+spvr(1,2)*svcr(2,2) 
        sssr(2,1) = spvr(2,1)*svcr(1,1)+spvr(2,2)*svcr(2,1) 
        sssr(2,2) = spvr(2,1)*svcr(1,2)+spvr(2,2)*svcr(2,2)

      vertex = vertex 
     &        +dconjg(gcc(1)) 
     &          *( foc(3) *( sssr(1,1)*rc(1,1)+sssr(1,2)*rc(1,2) )
     &            +foc(4) *( sssr(2,1)*rc(1,1)+sssr(2,2)*rc(1,2) )
     &            +foc(3) *( sssr(2,1)*rc(2,1)+sssr(2,2)*rc(2,2) )
     &            +foc(4) *( sssr(1,1)*rc(2,1)+sssr(1,2)*rc(2,2) )
     &            +foc(3) *(-sssr(2,1)*rc(3,1)-sssr(2,2)*rc(3,2) )*cImag
     &            +foc(4) *( sssr(1,1)*rc(3,1)+sssr(1,2)*rc(3,2) )*cImag
     &            +foc(3) *( sssr(1,1)*rc(4,1)+sssr(1,2)*rc(4,2) )
     &            +foc(4) *(-sssr(2,1)*rc(4,1)-sssr(2,2)*rc(4,2) )
     &           )

        spvr(1,1) = vc2(1)-vc2(4)
        spvr(1,2) = -(vc2(2)-cImag*vc2(3))
        spvr(2,1) = -(vc2(2)+cImag*vc2(3))
        spvr(2,2) = vc2(1)+vc2(4)
        
        svcr(1,1) = vc1(1)+vc1(4)
        svcr(1,2) = vc1(2)-cImag*vc1(3)
        svcr(2,1) = vc1(2)+cImag*vc1(3)
        svcr(2,2) = vc1(1)-vc1(4)

        sssr(1,1) = spvr(1,1)*svcr(1,1)+spvr(1,2)*svcr(2,1)
        sssr(1,2) = spvr(1,1)*svcr(1,2)+spvr(1,2)*svcr(2,2) 
        sssr(2,1) = spvr(2,1)*svcr(1,1)+spvr(2,2)*svcr(2,1) 
        sssr(2,2) = spvr(2,1)*svcr(1,2)+spvr(2,2)*svcr(2,2)

      vertex = vertex 
     &        -dconjg(gcc(1)) 
     &          *( foc(3) *( sssr(1,1)*rc(1,1)+sssr(1,2)*rc(1,2) )
     &            +foc(4) *( sssr(2,1)*rc(1,1)+sssr(2,2)*rc(1,2) )
     &            +foc(3) *( sssr(2,1)*rc(2,1)+sssr(2,2)*rc(2,2) )
     &            +foc(4) *( sssr(1,1)*rc(2,1)+sssr(1,2)*rc(2,2) )
     &            +foc(3) *(-sssr(2,1)*rc(3,1)-sssr(2,2)*rc(3,2) )*cImag
     &            +foc(4) *( sssr(1,1)*rc(3,1)+sssr(1,2)*rc(3,2) )*cImag
     &            +foc(3) *( sssr(1,1)*rc(4,1)+sssr(1,2)*rc(4,2) )
     &            +foc(4) *(-sssr(2,1)*rc(4,1)-sssr(2,2)*rc(4,2) )
     &           )

      endif

      vertex = -vertex

      return
      end

      subroutine irovvx(ri,fo,vc1,vc2,gc , vertex)
c
c This subroutine computes an amplitude of the fermion Ratrita-Schwinger
c fermion vector vector coupling.
c
c input:
c       complex ri(18)         : flow-in  Rarita-Schwinger fermion  |ri>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex vc1(6)         : input    vector                    v1
c       complex vc2(6)         : input    vector                    v2
c       complex gc(2)          : coupling constants                 gfrvv
c
c output:
c       complex vertex         : amplitude               <fo|[v1,v2]|ri>
c     
c- by Y.Takaesu - 2010/12/13
c
      implicit none
      double complex ri(18),fo(6),vc1(6),vc2(6),gc(2),vertex
      double complex rc(4,4),spvl(2,2),svcl(2,2),sssl(2,2),
     &                       spvr(2,2),svcr(2,2),sssr(2,2)
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
      p0 = -dble( fo(5))
      p1 = -dble( fo(6))
      p2 = -dimag(fo(6))
      p3 = -dimag(fo(5))
      q0 = dble( ri(17))
      q1 = dble( ri(18))
      q2 = dimag(ri(18))
      q3 = dimag(ri(17))
      r0 = dble( vc1(5))
      r1 = dble( vc1(6))
      r2 = dimag(vc1(6))
      r3 = dimag(vc1(5))
      s0 = dble( vc2(5))
      s1 = dble( vc2(6))
      s2 = dimag(vc2(6))
      s3 = dimag(vc2(5))
      if ( abs(fo(1))+abs(fo(2))+abs(fo(3))+abs(fo(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : fo in irovvx is zero spinor'
      endif
      if ( abs(fo(5))+abs(fo(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fo in irovvx has zero momentum'
      endif
      if (abs(ri(1))+abs(ri(2))+abs(ri(3))+abs(ri(4))+abs(ri(5))
     & +abs(ri(6))+abs(ri(7))+abs(ri(8))abs(ri(9))+abs(ri(10))
     & +abs(ri(11))+abs(ri(12))+abs(ri(13))+abs(ri(14))+abs(ri(15))
     & +abs(ri(16)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : ri in irovvx is zero spinor'
      endif
      if ( abs(ri(17))+abs(ri(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : ro in irovvx has zero momentum'
      endif
      if ( abs(vc1(1))+abs(vc1(2))+abs(vc1(3))+abs(vc1(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : vc1 in irovvx is zero vector'
      endif
      if ( abs(vc1(5))+abs(vc1(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : vc1 in irovvx has zero momentum'
      endif
      if ( abs(vc2(1))+abs(vc2(2))+abs(vc2(3))+abs(vc2(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : vc2 in irovvx is zero vector'
      endif
      if ( abs(vc2(5))+abs(vc2(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : vc2 in irovvx has zero momentum'
      endif
      pm = max( abs(p0),abs(q0),abs(r0),abs(s0),abs(p1),abs(q1),abs(r1),
     &          abs(s1),abs(p2),abs(q2),abs(r2),abs(s2),abs(p3),abs(q3),
     &          abs(r3),abs(s3) )
      if ( abs(-ri(17)+fo(5)+vc1(5)*vc2(5))+abs(-ri(18)+fo(6)+vc(6)+vc2(6))
     &                                              .ge.pm*epsi ) then
         write(stdo,*)
     &        ' helas-error : ri,fo,vc1,vc2 in irovvx'
         write(stdo,*)
     &        '                        have not balanced momenta'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*) ' helas-error : gc in irovvx is zero coupling'
      endif
#endif

      rc(1,1) = ri(1)
      rc(1,2) = ri(2)
      rc(1,3) = ri(3)
      rc(1,4) = ri(4)
      rc(2,1) = ri(5)
      rc(2,2) = ri(6)
      rc(2,3) = ri(7)
      rc(2,4) = ri(8)
      rc(3,1) = ri(9)
      rc(3,2) = ri(10)
      rc(3,3) = ri(11)
      rc(3,4) = ri(12)
      rc(4,1) = ri(13)
      rc(4,2) = ri(14)
      rc(4,3) = ri(15)
      rc(4,4) = ri(16)

      vertex = cZero

      if ( gc(2).ne.cZero ) then

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

      vertex = dconjg(gc(2)) 
     &          *( fo(1) *( sssl(1,1)*rc(1,3)+sssl(1,2)*rc(1,4) )
     &            +fo(2) *( sssl(2,1)*rc(1,3)+sssl(2,2)*rc(1,4) )
     &            -fo(1) *( sssl(2,1)*rc(2,3)+sssl(2,2)*rc(2,4) )
     &            -fo(2) *( sssl(1,1)*rc(2,3)+sssl(1,2)*rc(2,4) )
     &            -fo(1) *(-sssl(2,1)*rc(3,3)-sssl(2,2)*rc(3,4) )*cImag
     &            -fo(2) *( sssl(1,1)*rc(3,3)+sssl(1,2)*rc(3,4) )*cImag
     &            -fo(1) *( sssl(1,1)*rc(4,3)+sssl(1,2)*rc(4,4) )
     &            -fo(2) *(-sssl(2,1)*rc(4,3)-sssl(2,2)*rc(4,4) )
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
     &        -dconjg(gc(2)) 
     &          *( fo(1) *( sssl(1,1)*rc(1,3)+sssl(1,2)*rc(1,4) )
     &            +fo(2) *( sssl(2,1)*rc(1,3)+sssl(2,2)*rc(1,4) )
     &            -fo(1) *( sssl(2,1)*rc(2,3)+sssl(2,2)*rc(2,4) )
     &            -fo(2) *( sssl(1,1)*rc(2,3)+sssl(1,2)*rc(2,4) )
     &            -fo(1) *(-sssl(2,1)*rc(3,3)-sssl(2,2)*rc(3,4) )*cImag
     &            -fo(2) *( sssl(1,1)*rc(3,3)+sssl(1,2)*rc(3,4) )*cImag
     &            -fo(1) *( sssl(1,1)*rc(4,3)+sssl(1,2)*rc(4,4) )
     &            -fo(2) *(-sssl(2,1)*rc(4,3)-sssl(2,2)*rc(4,4) )
     &           )
      
      endif

      if ( gc(1).ne.cZero ) then

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
     &        +dconjg(gc(1)) 
     &          *( fo(3) *( sssr(1,1)*rc(1,1)+sssr(1,2)*rc(1,2) )
     &            +fo(4) *( sssr(2,1)*rc(1,1)+sssr(2,2)*rc(1,2) )
     &            +fo(3) *( sssr(2,1)*rc(2,1)+sssr(2,2)*rc(2,2) )
     &            +fo(4) *( sssr(1,1)*rc(2,1)+sssr(1,2)*rc(2,2) )
     &            +fo(3) *(-sssr(2,1)*rc(3,1)-sssr(2,2)*rc(3,2) )*cImag
     &            +fo(4) *( sssr(1,1)*rc(3,1)+sssr(1,2)*rc(3,2) )*cImag
     &            +fo(3) *( sssr(1,1)*rc(4,1)+sssr(1,2)*rc(4,2) )
     &            +fo(4) *(-sssr(2,1)*rc(4,1)-sssr(2,2)*rc(4,2) )
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
     &        -dconjg(gc(1)) 
     &          *( fo(3) *( sssr(1,1)*rc(1,1)+sssr(1,2)*rc(1,2) )
     &            +fo(4) *( sssr(2,1)*rc(1,1)+sssr(2,2)*rc(1,2) )
     &            +fo(3) *( sssr(2,1)*rc(2,1)+sssr(2,2)*rc(2,2) )
     &            +fo(4) *( sssr(1,1)*rc(2,1)+sssr(1,2)*rc(2,2) )
     &            +fo(3) *(-sssr(2,1)*rc(3,1)-sssr(2,2)*rc(3,2) )*cImag
     &            +fo(4) *( sssr(1,1)*rc(3,1)+sssr(1,2)*rc(3,2) )*cImag
     &            +fo(3) *( sssr(1,1)*rc(4,1)+sssr(1,2)*rc(4,2) )
     &            +fo(4) *(-sssr(2,1)*rc(4,1)-sssr(2,2)*rc(4,2) )
     &           )

      endif

      return
      end

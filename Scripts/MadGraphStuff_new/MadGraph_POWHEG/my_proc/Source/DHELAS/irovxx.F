      subroutine irovxx(ri,fo,vc,gc , vertex)
c
c This subroutine computes an amplitude of the fermion-(Rarita-Schwinger
c fermion)-vector coupling.
c
c input:
c       complex ri(18)         : flow-in  RS fermion               |ri>v
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex vc(6)          : input    vector                   v'(q)
c       complex gc(2)          : coupling constants                 GFRV
c
c output:
c       complex vertex         : amplitude               <fo|v[v',q]|ri>
c     
c- by Y.Takaesu - 2010/01/09
c
      implicit none
      double complex ri(18),fo(6),vc(6),gc(2),vertex
      double complex rc(4,4),spvl(2,2),svcl(2,2),sssl(2,2),
     &                       spvr(2,2),svcr(2,2),sssr(2,2)
      double precision pv(4)

      double precision rZero, rOne
      parameter( rZero = 0.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )

#ifdef HELAS_CHECK
      double precision p0,p1,p2,p3,q0,q1,q2,q3,r0,r1,r2,r3,pm
      double precision epsi
      parameter( epsi = 4.0d-5 )
      integer stdo
      parameter( stdo = 6 )
#endif
c
#ifdef HELAS_CHECK
      p0 = -dble( ri(17))
      p1 = -dble( ri(18))
      p2 = -dimag(ri(18))
      p3 = -dimag(ri(17))
      q0 = dble( fo(5))
      q1 = dble( fo(6))
      q2 = dimag(fo(6))
      q3 = dimag(fo(5))
      r0 = dble( vc(5))
      r1 = dble( vc(6))
      r2 = dimag(vc(6))
      r3 = dimag(vc(5))
      if ( abs(ri(17))+abs(ri(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : ri in irovxx has zero momentum'
      endif
      if ( abs(fo(5))+abs(fo(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fo in irovxx has zero momentum'
      endif
      if ( abs(vc(1))+abs(vc(2))+abs(vc(3))+abs(vc(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : vc in irovxx is zero vector'
      endif
      if ( abs(vc(5))+abs(vc(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : vc in irovxx has zero momentum'
      endif
      pm = max( abs(p0),abs(q0),abs(r0),abs(p1),abs(q1),abs(r1),
     &          abs(p2),abs(q2),abs(r2),abs(p3),abs(q3),abs(r3) )
      if ( abs(-ri(17)+fo(5)+vc(5))+abs(-ri(18)+fo(6)+vc(6))
     &                                              .ge.pm*epsi ) then
         write(stdo,*)
     &        ' helas-error : ri,fo,vc in irovxx'
         write(stdo,*)
     &        '                        have not balanced momenta'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*) ' helas-error : gc in irovxx is zero coupling'
      endif
#endif

      pv(1) = dble( vc(5))
      pv(2) = dble( vc(6))
      pv(3) = dimag(vc(6))
      pv(4) = dimag(vc(5))

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

        spvl(1,1) = pv(1)+pv(4)
        spvl(1,2) = pv(2)-cImag*pv(3)
        spvl(2,1) = pv(2)+cImag*pv(3)
        spvl(2,2) = pv(1)-pv(4)
        
        svcl(1,1) = vc(1)-vc(4)
        svcl(1,2) = -(vc(2)-cImag*vc(3))
        svcl(2,1) = -(vc(2)+cImag*vc(3))
        svcl(2,2) = vc(1)+vc(4)

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

        spvl(1,1) = vc(1)+vc(4)
        spvl(1,2) = vc(2)-cImag*vc(3)
        spvl(2,1) = vc(2)+cImag*vc(3)
        spvl(2,2) = vc(1)-vc(4)
        
        svcl(1,1) = pv(1)-pv(4)
        svcl(1,2) = -(pv(2)-cImag*pv(3))
        svcl(2,1) = -(pv(2)+cImag*pv(3))
        svcl(2,2) = pv(1)+pv(4)

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

        spvr(1,1) = pv(1)-pv(4)
        spvr(1,2) = -(pv(2)-cImag*pv(3))
        spvr(2,1) = -(pv(2)+cImag*pv(3))
        spvr(2,2) = pv(1)+pv(4)
        
        svcr(1,1) = vc(1)+vc(4)
        svcr(1,2) = vc(2)-cImag*vc(3)
        svcr(2,1) = vc(2)+cImag*vc(3)
        svcr(2,2) = vc(1)-vc(4)

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

        spvr(1,1) = vc(1)-vc(4)
        spvr(1,2) = -(vc(2)-cImag*vc(3))
        spvr(2,1) = -(vc(2)+cImag*vc(3))
        spvr(2,2) = vc(1)+vc(4)
        
        svcr(1,1) = pv(1)+pv(4)
        svcr(1,2) = pv(2)-cImag*pv(3)
        svcr(2,1) = pv(2)+cImag*pv(3)
        svcr(2,2) = pv(1)-pv(4)

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

      subroutine iorvxx(fi,ro,vc,gc , vertex)
c
c This subroutine computes an amplitude of the fermion-(Rarita-Schwinger 
c fermion)-vector coupling.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex ro(18)         : flow-out RS fermion               v<ro|
c       complex vc(6)          : input    vector                   v'(q)
c       complex gc(2)          : coupling constants                 GFRV
c
c output:
c       complex vertex         : amplitude               <ro|[q,v']v|fi>
c     
c- by K.Mawatari - 2008/02/27
c
      implicit none
      double complex fi(6),ro(18),vc(6),gc(2),vertex
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
      p0 = -dble( fi(5))
      p1 = -dble( fi(6))
      p2 = -dimag(fi(6))
      p3 = -dimag(fi(5))
      q0 = dble( ro(17))
      q1 = dble( ro(18))
      q2 = dimag(ro(18))
      q3 = dimag(ro(17))
      r0 = dble( vc(5))
      r1 = dble( vc(6))
      r2 = dimag(vc(6))
      r3 = dimag(vc(5))
      if ( abs(fi(1))+abs(fi(2))+abs(fi(3))+abs(fi(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : fi in iorvxx is zero spinor'
      endif
      if ( abs(fi(5))+abs(fi(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fi in iorvxx has zero momentum'
      endif
	if (abs(ro(1))+abs(ro(2))+abs(ro(3))+abs(ro(4))+abs(ro(5))
     & +abs(ro(6))+abs(ro(7))+abs(ro(8))abs(ro(9))+abs(ro(10))
     & +abs(ro(11))+abs(ro(12))+abs(ro(13))+abs(ro(14))+abs(ro(15))
     & +abs(ro(16)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : ro in iorvxx is zero spinor'
      endif
      if ( abs(ro(17))+abs(ro(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : ro in iorvxx has zero momentum'
      endif
      if ( abs(vc(1))+abs(vc(2))+abs(vc(3))+abs(vc(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : vc in iorvxx is zero vector'
      endif
      if ( abs(vc(5))+abs(vc(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : vc in iorvxx has zero momentum'
      endif
      pm = max( abs(p0),abs(q0),abs(r0),abs(p1),abs(q1),abs(r1),
     &          abs(p2),abs(q2),abs(r2),abs(p3),abs(q3),abs(r3) )
      if ( abs(-fi(5)+ro(17)+vc(5))+abs(-fi(6)+ro(18)+vc(6))
     &                                              .ge.pm*epsi ) then
         write(stdo,*)
     &        ' helas-error : fi,ro,vc in iorvxx'
         write(stdo,*)
     &        '                        have not balanced momenta'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*) ' helas-error : gc in iorvxx is zero coupling'
      endif
#endif

      pv(1) = dble( vc(5))
      pv(2) = dble( vc(6))
      pv(3) = dimag(vc(6))
      pv(4) = dimag(vc(5))

      rc(1,1) = ro(1)
      rc(1,2) = ro(2)
      rc(1,3) = ro(3)
      rc(1,4) = ro(4)
      rc(2,1) = ro(5)
      rc(2,2) = ro(6)
      rc(2,3) = ro(7)
      rc(2,4) = ro(8)
      rc(3,1) = ro(9)
      rc(3,2) = ro(10)
      rc(3,3) = ro(11)
      rc(3,4) = ro(12)
      rc(4,1) = ro(13)
      rc(4,2) = ro(14)
      rc(4,3) = ro(15)
      rc(4,4) = ro(16)

      vertex = cZero

      if ( gc(2).ne.cZero ) then

        spvl(1,1) = pv(1)-pv(4)
        spvl(1,2) = -(pv(2)-cImag*pv(3))
        spvl(2,1) = -(pv(2)+cImag*pv(3))
        spvl(2,2) = pv(1)+pv(4)
        
        svcl(1,1) = vc(1)+vc(4)
        svcl(1,2) = vc(2)-cImag*vc(3)
        svcl(2,1) = vc(2)+cImag*vc(3)
        svcl(2,2) = vc(1)-vc(4)

        sssl(1,1) = spvl(1,1)*svcl(1,1)+spvl(1,2)*svcl(2,1)
        sssl(1,2) = spvl(1,1)*svcl(1,2)+spvl(1,2)*svcl(2,2) 
        sssl(2,1) = spvl(2,1)*svcl(1,1)+spvl(2,2)*svcl(2,1) 
        sssl(2,2) = spvl(2,1)*svcl(1,2)+spvl(2,2)*svcl(2,2)

      vertex = gc(2) 
     &          *( rc(1,1) *( sssl(1,1)*fi(3)+sssl(1,2)*fi(4) )
     &            +rc(1,2) *( sssl(2,1)*fi(3)+sssl(2,2)*fi(4) )           
     &            -rc(2,1) *( sssl(1,2)*fi(3)+sssl(1,1)*fi(4) )
     &            -rc(2,2) *( sssl(2,2)*fi(3)+sssl(2,1)*fi(4) )
     &            -rc(3,1) *( sssl(1,2)*fi(3)-sssl(1,1)*fi(4) )*cImag
     &            -rc(3,2) *( sssl(2,2)*fi(3)-sssl(2,1)*fi(4) )*cImag
     &            -rc(4,1) *( sssl(1,1)*fi(3)-sssl(1,2)*fi(4) )
     &            -rc(4,2) *( sssl(2,1)*fi(3)-sssl(2,2)*fi(4) ) 
     &           )

        spvl(1,1) = vc(1)-vc(4)
        spvl(1,2) = -(vc(2)-cImag*vc(3))
        spvl(2,1) = -(vc(2)+cImag*vc(3))
        spvl(2,2) = vc(1)+vc(4)
        
        svcl(1,1) = pv(1)+pv(4)
        svcl(1,2) = pv(2)-cImag*pv(3)
        svcl(2,1) = pv(2)+cImag*pv(3)
        svcl(2,2) = pv(1)-pv(4)

        sssl(1,1) = spvl(1,1)*svcl(1,1)+spvl(1,2)*svcl(2,1)
        sssl(1,2) = spvl(1,1)*svcl(1,2)+spvl(1,2)*svcl(2,2) 
        sssl(2,1) = spvl(2,1)*svcl(1,1)+spvl(2,2)*svcl(2,1) 
        sssl(2,2) = spvl(2,1)*svcl(1,2)+spvl(2,2)*svcl(2,2)

      vertex = vertex
     &        -gc(2) 
     &          *( rc(1,1) *( sssl(1,1)*fi(3)+sssl(1,2)*fi(4) )
     &            +rc(1,2) *( sssl(2,1)*fi(3)+sssl(2,2)*fi(4) )           
     &            -rc(2,1) *( sssl(1,2)*fi(3)+sssl(1,1)*fi(4) )
     &            -rc(2,2) *( sssl(2,2)*fi(3)+sssl(2,1)*fi(4) )
     &            -rc(3,1) *( sssl(1,2)*fi(3)-sssl(1,1)*fi(4) )*cImag
     &            -rc(3,2) *( sssl(2,2)*fi(3)-sssl(2,1)*fi(4) )*cImag
     &            -rc(4,1) *( sssl(1,1)*fi(3)-sssl(1,2)*fi(4) )
     &            -rc(4,2) *( sssl(2,1)*fi(3)-sssl(2,2)*fi(4) ) 
     &           )
      
      endif

      if ( gc(1).ne.cZero ) then

        spvr(1,1) = pv(1)+pv(4)
        spvr(1,2) = pv(2)-cImag*pv(3)
        spvr(2,1) = pv(2)+cImag*pv(3)
        spvr(2,2) = pv(1)-pv(4)
        
        svcr(1,1) = vc(1)-vc(4)
        svcr(1,2) = -(vc(2)-cImag*vc(3))
        svcr(2,1) = -(vc(2)+cImag*vc(3))
        svcr(2,2) = vc(1)+vc(4)

        sssr(1,1) = spvr(1,1)*svcr(1,1)+spvr(1,2)*svcr(2,1)
        sssr(1,2) = spvr(1,1)*svcr(1,2)+spvr(1,2)*svcr(2,2) 
        sssr(2,1) = spvr(2,1)*svcr(1,1)+spvr(2,2)*svcr(2,1) 
        sssr(2,2) = spvr(2,1)*svcr(1,2)+spvr(2,2)*svcr(2,2)

      vertex = vertex 
     &        +gc(1) 
     &          *( rc(1,3) *( sssr(1,1)*fi(1)+sssr(1,2)*fi(2) )
     &            +rc(1,4) *( sssr(2,1)*fi(1)+sssr(2,2)*fi(2) )           
     &            +rc(2,3) *( sssr(1,2)*fi(1)+sssr(1,1)*fi(2) )
     &            +rc(2,4) *( sssr(2,2)*fi(1)+sssr(2,1)*fi(2) )
     &            +rc(3,3) *( sssr(1,2)*fi(1)-sssr(1,1)*fi(2) )*cImag
     &            +rc(3,4) *( sssr(2,2)*fi(1)-sssr(2,1)*fi(2) )*cImag
     &            +rc(4,3) *( sssr(1,1)*fi(1)-sssr(1,2)*fi(2) )
     &            +rc(4,4) *( sssr(2,1)*fi(1)-sssr(2,2)*fi(2) ) 
     &           )

        spvr(1,1) = vc(1)+vc(4)
        spvr(1,2) = vc(2)-cImag*vc(3)
        spvr(2,1) = vc(2)+cImag*vc(3)
        spvr(2,2) = vc(1)-vc(4)
        
        svcr(1,1) = pv(1)-pv(4)
        svcr(1,2) = -(pv(2)-cImag*pv(3))
        svcr(2,1) = -(pv(2)+cImag*pv(3))
        svcr(2,2) = pv(1)+pv(4)

        sssr(1,1) = spvr(1,1)*svcr(1,1)+spvr(1,2)*svcr(2,1)
        sssr(1,2) = spvr(1,1)*svcr(1,2)+spvr(1,2)*svcr(2,2) 
        sssr(2,1) = spvr(2,1)*svcr(1,1)+spvr(2,2)*svcr(2,1) 
        sssr(2,2) = spvr(2,1)*svcr(1,2)+spvr(2,2)*svcr(2,2)

      vertex = vertex 
     &        -gc(1) 
     &          *( rc(1,3) *( sssr(1,1)*fi(1)+sssr(1,2)*fi(2) )
     &            +rc(1,4) *( sssr(2,1)*fi(1)+sssr(2,2)*fi(2) )           
     &            +rc(2,3) *( sssr(1,2)*fi(1)+sssr(1,1)*fi(2) )
     &            +rc(2,4) *( sssr(2,2)*fi(1)+sssr(2,1)*fi(2) )
     &            +rc(3,3) *( sssr(1,2)*fi(1)-sssr(1,1)*fi(2) )*cImag
     &            +rc(3,4) *( sssr(2,2)*fi(1)-sssr(2,1)*fi(2) )*cImag
     &            +rc(4,3) *( sssr(1,1)*fi(1)-sssr(1,2)*fi(2) )
     &            +rc(4,4) *( sssr(2,1)*fi(1)-sssr(2,2)*fi(2) ) 
     &           )

      endif

      return
      end

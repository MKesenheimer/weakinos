      subroutine iorvvc(fic,roc,vc1,vc2,gc , vertex)
c
c This subroutine computes an amplitude of the antifermion-(antiRatrita-Schwinger
c fermion)-vector-vector coupling.
c
c input:
c       complex fic(6)         : flow-in  antifermion              |fic>
c       complex roc(18)        : flow-out antiRSfermion            <roc|
c       complex vc1(6)         : input    vector                    v1
c       complex vc2(6)         : input    vector                    v2
c       complex gc(2)          : coupling constants                 gfrvv
c
c output:
c       complex vertex         : amplitude               <roc|[v1,v2]|fic>
c     
c- by Yoshitaro Takaesu - 2011/03/26
c
      implicit none
      double complex fic(6),roc(18),vc1(6),vc2(6),gc(2),vertex
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
      p0 = -dble( fic(5))
      p1 = -dble( fic(6))
      p2 = -dimag(fic(6))
      p3 = -dimag(fic(5))
      q0 = dble( roc(17))
      q1 = dble( roc(18))
      q2 = dimag(roc(18))
      q3 = dimag(roc(17))
      r0 = dble( vc1(5))
      r1 = dble( vc1(6))
      r2 = dimag(vc1(6))
      r3 = dimag(vc1(5))
      s0 = dble( vc2(5))
      s1 = dble( vc2(6))
      s2 = dimag(vc2(6))
      s3 = dimag(vc2(5))
      if ( abs(fic(1))+abs(fic(2))+abs(fic(3))+abs(fic(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : fic in iorvvc is zero spinor'
      endif
      if ( abs(fic(5))+abs(fic(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fic in iorvvc has zero momentum'
      endif
      if (abs(roc(1))+abs(roc(2))+abs(roc(3))+abs(roc(4))+abs(roc(5))
     & +abs(roc(6))+abs(roc(7))+abs(roc(8))abs(roc(9))+abs(roc(10))
     & +abs(roc(11))+abs(roc(12))+abs(roc(13))+abs(roc(14))+abs(roc(15))
     & +abs(roc(16)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : roc in iorvvc is zero spinor'
      endif
      if ( abs(roc(17))+abs(roc(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : roc in iorvvc has zero momentum'
      endif
      if ( abs(vc1(1))+abs(vc1(2))+abs(vc1(3))+abs(vc1(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : vc1 in iorvvc is zero vector'
      endif
      if ( abs(vc1(5))+abs(vc1(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : vc1 in iorvvc has zero momentum'
      endif
      if ( abs(vc2(1))+abs(vc2(2))+abs(vc2(3))+abs(vc2(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : vc2 in iorvvc is zero vector'
      endif
      if ( abs(vc2(5))+abs(vc2(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : vc2 in iorvvc has zero momentum'
      endif
      pm = max( abs(p0),abs(q0),abs(r0),abs(s0),abs(p1),abs(q1),abs(r1),
     &          abs(s1),abs(p2),abs(q2),abs(r2),abs(s2),abs(p3),abs(q3),
     &          abs(r3),abs(s3) )
      if ( abs(-fic(5)+roc(17)+vc1(5)*vc2(5))+abs(-fic(6)+roc(18)+vc(6)+vc2(6))
     &                                              .ge.pm*epsi ) then
         write(stdo,*)
     &        ' helas-error : fic,roc,vc1,vc2 in iorvvc'
         write(stdo,*)
     &        '                        have not balanced momenta'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*) ' helas-error : gc in iorvvc is zero coupling'
      endif
#endif

      gcc(1) = dconjg(gc(2))
      gcc(2) = dconjg(gc(1))

      rc(1,1) = roc(1)
      rc(1,2) = roc(2)
      rc(1,3) = roc(3)
      rc(1,4) = roc(4)
      rc(2,1) = roc(5)
      rc(2,2) = roc(6)
      rc(2,3) = roc(7)
      rc(2,4) = roc(8)
      rc(3,1) = roc(9)
      rc(3,2) = roc(10)
      rc(3,3) = roc(11)
      rc(3,4) = roc(12)
      rc(4,1) = roc(13)
      rc(4,2) = roc(14)
      rc(4,3) = roc(15)
      rc(4,4) = roc(16)

      vertex = cZero

      if ( gcc(2).ne.cZero ) then

        spvl(1,1) = vc1(1)-vc1(4)
        spvl(1,2) = -(vc1(2)-cImag*vc1(3))
        spvl(2,1) = -(vc1(2)+cImag*vc1(3))
        spvl(2,2) = vc1(1)+vc1(4)
        
        svcl(1,1) = vc2(1)+vc2(4)
        svcl(1,2) = vc2(2)-cImag*vc2(3)
        svcl(2,1) = vc2(2)+cImag*vc2(3)
        svcl(2,2) = vc2(1)-vc2(4)

        sssl(1,1) = spvl(1,1)*svcl(1,1)+spvl(1,2)*svcl(2,1)
        sssl(1,2) = spvl(1,1)*svcl(1,2)+spvl(1,2)*svcl(2,2) 
        sssl(2,1) = spvl(2,1)*svcl(1,1)+spvl(2,2)*svcl(2,1) 
        sssl(2,2) = spvl(2,1)*svcl(1,2)+spvl(2,2)*svcl(2,2)

      vertex = gcc(2) 
     &          *( rc(1,1) *( sssl(1,1)*fic(3)+sssl(1,2)*fic(4) )
     &            +rc(1,2) *( sssl(2,1)*fic(3)+sssl(2,2)*fic(4) )           
     &            -rc(2,1) *( sssl(1,2)*fic(3)+sssl(1,1)*fic(4) )
     &            -rc(2,2) *( sssl(2,2)*fic(3)+sssl(2,1)*fic(4) )
     &            -rc(3,1) *( sssl(1,2)*fic(3)-sssl(1,1)*fic(4) )*cImag
     &            -rc(3,2) *( sssl(2,2)*fic(3)-sssl(2,1)*fic(4) )*cImag
     &            -rc(4,1) *( sssl(1,1)*fic(3)-sssl(1,2)*fic(4) )
     &            -rc(4,2) *( sssl(2,1)*fic(3)-sssl(2,2)*fic(4) ) 
     &           )

        spvl(1,1) = vc2(1)-vc2(4)
        spvl(1,2) = -(vc2(2)-cImag*vc2(3))
        spvl(2,1) = -(vc2(2)+cImag*vc2(3))
        spvl(2,2) = vc2(1)+vc2(4)
        
        svcl(1,1) = vc1(1)+vc1(4)
        svcl(1,2) = vc1(2)-cImag*vc1(3)
        svcl(2,1) = vc1(2)+cImag*vc1(3)
        svcl(2,2) = vc1(1)-vc1(4)

        sssl(1,1) = spvl(1,1)*svcl(1,1)+spvl(1,2)*svcl(2,1)
        sssl(1,2) = spvl(1,1)*svcl(1,2)+spvl(1,2)*svcl(2,2) 
        sssl(2,1) = spvl(2,1)*svcl(1,1)+spvl(2,2)*svcl(2,1) 
        sssl(2,2) = spvl(2,1)*svcl(1,2)+spvl(2,2)*svcl(2,2)

      vertex = vertex
     &        -gcc(2) 
     &          *( rc(1,1) *( sssl(1,1)*fic(3)+sssl(1,2)*fic(4) )
     &            +rc(1,2) *( sssl(2,1)*fic(3)+sssl(2,2)*fic(4) )           
     &            -rc(2,1) *( sssl(1,2)*fic(3)+sssl(1,1)*fic(4) )
     &            -rc(2,2) *( sssl(2,2)*fic(3)+sssl(2,1)*fic(4) )
     &            -rc(3,1) *( sssl(1,2)*fic(3)-sssl(1,1)*fic(4) )*cImag
     &            -rc(3,2) *( sssl(2,2)*fic(3)-sssl(2,1)*fic(4) )*cImag
     &            -rc(4,1) *( sssl(1,1)*fic(3)-sssl(1,2)*fic(4) )
     &            -rc(4,2) *( sssl(2,1)*fic(3)-sssl(2,2)*fic(4) ) 
     &           )
      
      endif

      if ( gcc(1).ne.cZero ) then

        spvr(1,1) = vc1(1)+vc1(4)
        spvr(1,2) = vc1(2)-cImag*vc1(3)
        spvr(2,1) = vc1(2)+cImag*vc1(3)
        spvr(2,2) = vc1(1)-vc1(4)
        
        svcr(1,1) = vc2(1)-vc2(4)
        svcr(1,2) = -(vc2(2)-cImag*vc2(3))
        svcr(2,1) = -(vc2(2)+cImag*vc2(3))
        svcr(2,2) = vc2(1)+vc2(4)

        sssr(1,1) = spvr(1,1)*svcr(1,1)+spvr(1,2)*svcr(2,1)
        sssr(1,2) = spvr(1,1)*svcr(1,2)+spvr(1,2)*svcr(2,2) 
        sssr(2,1) = spvr(2,1)*svcr(1,1)+spvr(2,2)*svcr(2,1) 
        sssr(2,2) = spvr(2,1)*svcr(1,2)+spvr(2,2)*svcr(2,2)

      vertex = vertex 
     &        +gcc(1) 
     &          *( rc(1,3) *( sssr(1,1)*fic(1)+sssr(1,2)*fic(2) )
     &            +rc(1,4) *( sssr(2,1)*fic(1)+sssr(2,2)*fic(2) )           
     &            +rc(2,3) *( sssr(1,2)*fic(1)+sssr(1,1)*fic(2) )
     &            +rc(2,4) *( sssr(2,2)*fic(1)+sssr(2,1)*fic(2) )
     &            +rc(3,3) *( sssr(1,2)*fic(1)-sssr(1,1)*fic(2) )*cImag
     &            +rc(3,4) *( sssr(2,2)*fic(1)-sssr(2,1)*fic(2) )*cImag
     &            +rc(4,3) *( sssr(1,1)*fic(1)-sssr(1,2)*fic(2) )
     &            +rc(4,4) *( sssr(2,1)*fic(1)-sssr(2,2)*fic(2) ) 
     &           )

        spvr(1,1) = vc2(1)+vc2(4)
        spvr(1,2) = vc2(2)-cImag*vc2(3)
        spvr(2,1) = vc2(2)+cImag*vc2(3)
        spvr(2,2) = vc2(1)-vc2(4)
        
        svcr(1,1) = vc1(1)-vc1(4)
        svcr(1,2) = -(vc1(2)-cImag*vc1(3))
        svcr(2,1) = -(vc1(2)+cImag*vc1(3))
        svcr(2,2) = vc1(1)+vc1(4)

        sssr(1,1) = spvr(1,1)*svcr(1,1)+spvr(1,2)*svcr(2,1)
        sssr(1,2) = spvr(1,1)*svcr(1,2)+spvr(1,2)*svcr(2,2) 
        sssr(2,1) = spvr(2,1)*svcr(1,1)+spvr(2,2)*svcr(2,1) 
        sssr(2,2) = spvr(2,1)*svcr(1,2)+spvr(2,2)*svcr(2,2)

      vertex = vertex 
     &        -gcc(1) 
     &          *( rc(1,3) *( sssr(1,1)*fic(1)+sssr(1,2)*fic(2) )
     &            +rc(1,4) *( sssr(2,1)*fic(1)+sssr(2,2)*fic(2) )           
     &            +rc(2,3) *( sssr(1,2)*fic(1)+sssr(1,1)*fic(2) )
     &            +rc(2,4) *( sssr(2,2)*fic(1)+sssr(2,1)*fic(2) )
     &            +rc(3,3) *( sssr(1,2)*fic(1)-sssr(1,1)*fic(2) )*cImag
     &            +rc(3,4) *( sssr(2,2)*fic(1)-sssr(2,1)*fic(2) )*cImag
     &            +rc(4,3) *( sssr(1,1)*fic(1)-sssr(1,2)*fic(2) )
     &            +rc(4,4) *( sssr(2,1)*fic(1)-sssr(2,2)*fic(2) ) 
     &           )

      endif

      vertex = -vertex

      return
      end

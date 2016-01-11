      subroutine jvirox(ri,fo,vc,gc,vmass,vwidth , jviro)
c
c This subroutine computes an off-shell vector current from an internal
c Rarita-Schwinger fermion, an external fermion, and a vector boson.
c The vector boson propagator is given in Feynman gauge for a massless 
c vector and in unitary gauge for a massive vector.
c
c input:
c       complex ri(6)          : flow-in  Rarita-Schwinger fermion  |ri>
c       complex fo(18)         : flow-out fermion                   <fo|
c       complex vc(6)          : vector boson                         v
c       complex gc(2)          : coupling constants                gfrvv
c       real    vmass          : mass  of OUTPUT vector v
c       real    vwidth         : width of OUTPUT vector v
c
c output:
c       complex jviro(6)       : vector current          j^mu(<fo|vv'|ri>)
c   
c- by Y.Takaesu - 2010/12/20
c  
      implicit none
      double complex ri(18),fo(6),gc(2),jviro(6),c0,c1,c2,c3,cs,d
      double complex vc(6)
      double precision q(0:3),vmass,vwidth,q2,vm2
      double complex cm2 ! mass**2- I Gamma mass (Fabio)


      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )

      double complex sq1(2,2),sq2(2,2),sqg1(0:3,2,2),sqg2(0:3,2,2)
      double complex gsq1(0:3,2,2),gsq2(0:3,2,2)
      double complex rc(0:3,4),sqgr(0:3,0:3,4),gsqgr(0:3,4)
      double complex spvl(2,2),svcl(2,2),sssl(2,2)
      double complex spvr(2,2),svcr(2,2),sssr(2,2),vertex1,vertex2
      integer i,j

#ifdef HELAS_CHECK
      integer stdo
      parameter( stdo = 6 )
#endif
c
#ifdef HELAS_CHECK
      if ( abs(fo(1))+abs(fo(2))+abs(fo(3))+abs(fo(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : fo in jvirox is zero spinor'
      endif
      if ( abs(fo(5))+abs(fo(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fo in jvirox has zero momentum'
      endif
      if (abs(ri(1))+abs(ri(2))+abs(ri(3))+abs(ri(4))+abs(ri(5))
     & +abs(ri(6))+abs(ri(7))+abs(ri(8))abs(ri(9))+abs(ri(10))
     & +abs(ri(11))+abs(ri(12))+abs(ri(13))+abs(ri(14))+abs(ri(15))
     & +abs(ri(16)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : ri in jvirox is zero spinor'
      endif
      if ( abs(ri(17))+abs(ri(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : ri in jvirox has zero momentum'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*) ' helas-error : gc in jvirox is zero coupling'
      endif
      if ( vmass.lt.rZero ) then
         write(stdo,*) ' helas-error : vmass in jvirox is negative'
         write(stdo,*) '             : vmass = ',vmass
      endif
      if ( vwidth.lt.rZero ) then
         write(stdo,*) ' helas-error : vwidth in jvirox is negative'
         write(stdo,*) '             : vwidth = ',vwidth
      endif
#endif

      jviro(5) = -ri(17)+fo(5)+vc(5)
      jviro(6) = -ri(18)+fo(6)+vc(6)

      q(0) = dble(jviro(5))
      q(1) = dble(jviro(6))
      q(2) = dimag(jviro(6))
      q(3) = dimag(jviro(5))
      
      
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)
      vm2 = vmass**2

#ifdef HELAS_CHECK
      if ( abs(jviro(5))+abs(jviro(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : jviro in jvirox has zero momentum'
      endif
      if ( vwidth.eq.rZero .and. q2.eq.vm2 ) then
         write(stdo,*)
     &        ' helas-error : jviro in jvirox is on vmass pole'
         write(stdo,*)
     &        '             : q     = ',q(0),q(1),q(2),q(3)
         write(stdo,*)
     &        '             : abs(q)= ',sqrt(abs(q2))
         jviro(1) = cZero
         jviro(2) = cZero
         jviro(3) = cZero
         jviro(4) = cZero
         return
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

      
      sq1(1,1) = vc(1)-vc(4)
      sq1(1,2) = -(vc(2)-cImag*vc(3))
      sq1(2,1) = -(vc(2)+cImag*vc(3))
      sq1(2,2) = vc(1)+vc(4)

      sq2(1,1) = vc(1)+vc(4)
      sq2(1,2) = vc(2)-cImag*vc(3)
      sq2(2,1) = vc(2)+cImag*vc(3)
      sq2(2,2) = vc(1)-vc(4)

      sqg1(0,1,1) = sq1(1,1)
      sqg1(0,1,2) = sq1(1,2)
      sqg1(0,2,1) = sq1(2,1)
      sqg1(0,2,2) = sq1(2,2)
      sqg1(1,1,1) = -sq1(1,2)
      sqg1(1,1,2) = -sq1(1,1)
      sqg1(1,2,1) = -sq1(2,2)
      sqg1(1,2,2) = -sq1(2,1)
      sqg1(2,1,1) = -cImag*sq1(1,2)
      sqg1(2,1,2) = cImag*sq1(1,1)
      sqg1(2,2,1) = -cImag*sq1(2,2)
      sqg1(2,2,2) = cImag*sq1(2,1)
      sqg1(3,1,1) = -sq1(1,1)
      sqg1(3,1,2) = sq1(1,2)
      sqg1(3,2,1) = -sq1(2,1)
      sqg1(3,2,2) = sq1(2,2)

      sqg2(0,1,1) = sq2(1,1)
      sqg2(0,1,2) = sq2(1,2)
      sqg2(0,2,1) = sq2(2,1)
      sqg2(0,2,2) = sq2(2,2)
      sqg2(1,1,1) = sq2(1,2)
      sqg2(1,1,2) = sq2(1,1)
      sqg2(1,2,1) = sq2(2,2)
      sqg2(1,2,2) = sq2(2,1)
      sqg2(2,1,1) = cImag*sq2(1,2)
      sqg2(2,1,2) = -cImag*sq2(1,1)
      sqg2(2,2,1) = cImag*sq2(2,2)
      sqg2(2,2,2) = -cImag*sq2(2,1)
      sqg2(3,1,1) = sq2(1,1)
      sqg2(3,1,2) = -sq2(1,2)
      sqg2(3,2,1) = sq2(2,1)
      sqg2(3,2,2) = -sq2(2,2)

      gsq1(0,1,1) = sq2(1,1)
      gsq1(0,1,2) = sq2(1,2)
      gsq1(0,2,1) = sq2(2,1)
      gsq1(0,2,2) = sq2(2,2)
      gsq1(1,1,1) = sq2(2,1)
      gsq1(1,1,2) = sq2(2,2)
      gsq1(1,2,1) = sq2(1,1)
      gsq1(1,2,2) = sq2(1,2)
      gsq1(2,1,1) = -cImag*sq2(2,1)
      gsq1(2,1,2) = -cImag*sq2(2,2)
      gsq1(2,2,1) = cImag*sq2(1,1)
      gsq1(2,2,2) = cImag*sq2(1,2)
      gsq1(3,1,1) = sq2(1,1)
      gsq1(3,1,2) = sq2(1,2)
      gsq1(3,2,1) = -sq2(2,1)
      gsq1(3,2,2) = -sq2(2,2)

      gsq2(0,1,1) = sq1(1,1)
      gsq2(0,1,2) = sq1(1,2)
      gsq2(0,2,1) = sq1(2,1)
      gsq2(0,2,2) = sq1(2,2)
      gsq2(1,1,1) = -sq1(2,1)
      gsq2(1,1,2) = -sq1(2,2)
      gsq2(1,2,1) = -sq1(1,1)
      gsq2(1,2,2) = -sq1(1,2)
      gsq2(2,1,1) = cImag*sq1(2,1)
      gsq2(2,1,2) = cImag*sq1(2,2)
      gsq2(2,2,1) = -cImag*sq1(1,1)
      gsq2(2,2,2) = -cImag*sq1(1,2)
      gsq2(3,1,1) = -sq1(1,1)
      gsq2(3,1,2) = -sq1(1,2)
      gsq2(3,2,1) = sq1(2,1)
      gsq2(3,2,2) = sq1(2,2)


      do i=0,3
       do j=0,3    
        sqgr(i,j,1) = rc(j,1)*(sqg1(i,1,1)-gsq1(i,1,1))+rc(j,2)
     &*(sqg1(i,1,2)-gsq1(i,1,2))
        sqgr(i,j,2) = rc(j,1)*(sqg1(i,2,1)-gsq1(i,2,1))+rc(j,2)
     &*(sqg1(i,2,2)-gsq1(i,2,2))   
        sqgr(i,j,3) = rc(j,3)*(sqg2(i,1,1)-gsq2(i,1,1))+rc(j,4)
     &*(sqg2(i,1,2)-gsq2(i,1,2))
        sqgr(i,j,4) = rc(j,3)*(sqg2(i,2,1)-gsq2(i,2,1))+rc(j,4)
     &*(sqg2(i,2,2)-gsq2(i,2,2))
       enddo
      enddo

      do i=0,3
       gsqgr(i,1) = sqgr(i,0,3)-sqgr(i,1,4)+cImag*sqgr(i,2,4)
     &              -sqgr(i,3,3) 
       gsqgr(i,2) = sqgr(i,0,4)-sqgr(i,1,3)-cImag*sqgr(i,2,3)
     &              +sqgr(i,3,4) 
       gsqgr(i,3) = sqgr(i,0,1)+sqgr(i,1,2)-cImag*sqgr(i,2,2)
     &              +sqgr(i,3,1)
       gsqgr(i,4) = sqgr(i,0,2)+sqgr(i,1,1)+cImag*sqgr(i,2,1)
     &              -sqgr(i,3,2) 
      enddo

       vertex1 = cZero
       vertex2 = cZero

      if ( gc(2).ne.cZero ) then

        spvl(1,1) = q(0)+q(3)
        spvl(1,2) = q(1)-cImag*q(2)
        spvl(2,1) = q(1)+cImag*q(2)
        spvl(2,2) = q(0)-q(3)
        
        svcl(1,1) = vc(1)-vc(4)
        svcl(1,2) = -(vc(2)-cImag*vc(3))
        svcl(2,1) = -(vc(2)+cImag*vc(3))
        svcl(2,2) = vc(1)+vc(4)

        sssl(1,1) = spvl(1,1)*svcl(1,1)+spvl(1,2)*svcl(2,1)
        sssl(1,2) = spvl(1,1)*svcl(1,2)+spvl(1,2)*svcl(2,2) 
        sssl(2,1) = spvl(2,1)*svcl(1,1)+spvl(2,2)*svcl(2,1) 
        sssl(2,2) = spvl(2,1)*svcl(1,2)+spvl(2,2)*svcl(2,2)     

      vertex2 = dconjg(gc(2)) 
     &          *( fo(1) *( sssl(1,1)*rc(0,3)+sssl(1,2)*rc(0,4) )
     &            +fo(2) *( sssl(2,1)*rc(0,3)+sssl(2,2)*rc(0,4) )
     &            -fo(1) *( sssl(2,1)*rc(1,3)+sssl(2,2)*rc(1,4) )
     &            -fo(2) *( sssl(1,1)*rc(1,3)+sssl(1,2)*rc(1,4) )
     &            -fo(1) *(-sssl(2,1)*rc(2,3)-sssl(2,2)*rc(2,4) )*cImag
     &            -fo(2) *( sssl(1,1)*rc(2,3)+sssl(1,2)*rc(2,4) )*cImag
     &            -fo(1) *( sssl(1,1)*rc(3,3)+sssl(1,2)*rc(3,4) )
     &            -fo(2) *(-sssl(2,1)*rc(3,3)-sssl(2,2)*rc(3,4) )
     &           )


        spvl(1,1) = vc(1)+vc(4)
        spvl(1,2) = vc(2)-cImag*vc(3)
        spvl(2,1) = vc(2)+cImag*vc(3)
        spvl(2,2) = vc(1)-vc(4)
        
        svcl(1,1) = q(0)-q(3)
        svcl(1,2) = -(q(1)-cImag*q(2))
        svcl(2,1) = -(q(1)+cImag*q(2))
        svcl(2,2) = q(0)+q(3)

        sssl(1,1) = spvl(1,1)*svcl(1,1)+spvl(1,2)*svcl(2,1)
        sssl(1,2) = spvl(1,1)*svcl(1,2)+spvl(1,2)*svcl(2,2) 
        sssl(2,1) = spvl(2,1)*svcl(1,1)+spvl(2,2)*svcl(2,1) 
        sssl(2,2) = spvl(2,1)*svcl(1,2)+spvl(2,2)*svcl(2,2)

      vertex2 = vertex2
     &        -dconjg(gc(2)) 
     &          *( fo(1) *( sssl(1,1)*rc(0,3)+sssl(1,2)*rc(0,4) )
     &            +fo(2) *( sssl(2,1)*rc(0,3)+sssl(2,2)*rc(0,4) )
     &            -fo(1) *( sssl(2,1)*rc(1,3)+sssl(2,2)*rc(1,4) )
     &            -fo(2) *( sssl(1,1)*rc(1,3)+sssl(1,2)*rc(1,4) )
     &            -fo(1) *(-sssl(2,1)*rc(2,3)-sssl(2,2)*rc(2,4) )*cImag
     &            -fo(2) *( sssl(1,1)*rc(2,3)+sssl(1,2)*rc(2,4) )*cImag
     &            -fo(1) *( sssl(1,1)*rc(3,3)+sssl(1,2)*rc(3,4) )
     &            -fo(2) *(-sssl(2,1)*rc(3,3)-sssl(2,2)*rc(3,4) )
     &           )
      
      endif

      if ( gc(1).ne.cZero ) then

        spvr(1,1) = q(0)-q(3)
        spvr(1,2) = -(q(1)-cImag*q(2))
        spvr(2,1) = -(q(1)+cImag*q(2))
        spvr(2,2) = q(0)+q(3)
        
        svcr(1,1) = vc(1)+vc(4)
        svcr(1,2) = vc(2)-cImag*vc(3)
        svcr(2,1) = vc(2)+cImag*vc(3)
        svcr(2,2) = vc(1)-vc(4)

        sssr(1,1) = spvr(1,1)*svcr(1,1)+spvr(1,2)*svcr(2,1)
        sssr(1,2) = spvr(1,1)*svcr(1,2)+spvr(1,2)*svcr(2,2) 
        sssr(2,1) = spvr(2,1)*svcr(1,1)+spvr(2,2)*svcr(2,1) 
        sssr(2,2) = spvr(2,1)*svcr(1,2)+spvr(2,2)*svcr(2,2)

      vertex1 =
     &        +dconjg(gc(1)) 
     &          *( fo(3) *( sssr(1,1)*rc(0,1)+sssr(1,2)*rc(0,2) )
     &            +fo(4) *( sssr(2,1)*rc(0,1)+sssr(2,2)*rc(0,2) )
     &            +fo(3) *( sssr(2,1)*rc(1,1)+sssr(2,2)*rc(1,2) )
     &            +fo(4) *( sssr(1,1)*rc(1,1)+sssr(1,2)*rc(1,2) )
     &            +fo(3) *(-sssr(2,1)*rc(2,1)-sssr(2,2)*rc(2,2) )*cImag
     &            +fo(4) *( sssr(1,1)*rc(2,1)+sssr(1,2)*rc(2,2) )*cImag
     &            +fo(3) *( sssr(1,1)*rc(3,1)+sssr(1,2)*rc(3,2) )
     &            +fo(4) *(-sssr(2,1)*rc(3,1)-sssr(2,2)*rc(3,2) )
     &           )

        spvr(1,1) = vc(1)-vc(4)
        spvr(1,2) = -(vc(2)-cImag*vc(3))
        spvr(2,1) = -(vc(2)+cImag*vc(3))
        spvr(2,2) = vc(1)+vc(4)
        
        svcr(1,1) = q(0)+q(3)
        svcr(1,2) = q(1)-cImag*q(2)
        svcr(2,1) = q(1)+cImag*q(2)
        svcr(2,2) = q(0)-q(3)

        sssr(1,1) = spvr(1,1)*svcr(1,1)+spvr(1,2)*svcr(2,1)
        sssr(1,2) = spvr(1,1)*svcr(1,2)+spvr(1,2)*svcr(2,2) 
        sssr(2,1) = spvr(2,1)*svcr(1,1)+spvr(2,2)*svcr(2,1) 
        sssr(2,2) = spvr(2,1)*svcr(1,2)+spvr(2,2)*svcr(2,2)

      vertex1 = vertex1 
     &        -dconjg(gc(1)) 
     &          *( fo(3) *( sssr(1,1)*rc(0,1)+sssr(1,2)*rc(0,2) )
     &            +fo(4) *( sssr(2,1)*rc(0,1)+sssr(2,2)*rc(0,2) )
     &            +fo(3) *( sssr(2,1)*rc(1,1)+sssr(2,2)*rc(1,2) )
     &            +fo(4) *( sssr(1,1)*rc(1,1)+sssr(1,2)*rc(1,2) )
     &            +fo(3) *(-sssr(2,1)*rc(2,1)-sssr(2,2)*rc(2,2) )*cImag
     &            +fo(4) *( sssr(1,1)*rc(2,1)+sssr(1,2)*rc(2,2) )*cImag
     &            +fo(3) *( sssr(1,1)*rc(3,1)+sssr(1,2)*rc(3,2) )
     &            +fo(4) *(-sssr(2,1)*rc(3,1)-sssr(2,2)*rc(3,2) )
     &           )

      endif


      jviro(1) = cZero
      jviro(2) = cZero
      jviro(3) = cZero
      jviro(4) = cZero

      if ( vmass.ne.rZero ) then

       d = rOne/dcmplx( q2-vm2, vmass*vwidth )

       if ( gc(2).ne.cZero ) then
        do i=1,4
         jviro(i) = jviro(i)-dconjg(gc(2))*d*(fo(1)*gsqgr(i-1,1)+fo(2)
     &             *gsqgr(i-1,2))-d*q(i-1)/vm2*vertex2
        enddo
       endif
      
       if ( gc(1).ne.cZero ) then
        do i=1,4
         jviro(i) = jviro(i)-dconjg(gc(1))*d*(fo(3)*gsqgr(i-1,3)+fo(4)
     &             *gsqgr(i-1,4))-d*q(i-1)/vm2*vertex1
        enddo
       endif

      else

       d = dcmplx( rOne/q2, rZero )

      if ( gc(2).ne.cZero ) then
        do i=1,4
         jviro(i) = jviro(i)-dconjg(gc(2))*d*(fo(1)*gsqgr(i-1,1)+fo(2)
     &             *gsqgr(i-1,2))
        enddo
       endif
    
       if ( gc(1).ne.cZero ) then
        do i=1,4
         jviro(i) = jviro(i)-dconjg(gc(1))*d*(fo(3)*gsqgr(i-1,3)+fo(4)
     &             *gsqgr(i-1,4))
        enddo
       endif
         

      end if
c
      return
      end

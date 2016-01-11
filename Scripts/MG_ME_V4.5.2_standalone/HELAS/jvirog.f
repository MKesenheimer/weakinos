c
c ----------------------------------------------------------------------
c
      subroutine jvirog(fi,fo,vc,gc,vmass,vwidth , jvirog1)
c
c This subroutine computes an off-shell vector current for the gaugino-
c goldstino-vector boson-vector boson vertex from an external fermion pair. The vector boson 
c propagator is given in feynman gauge for a massless vector and in 
c unitary gauge for a massive vector.
c
c input:
c       complex fi(6)          : flow-in goldstino
c       complex fo(6)          : flow-out gaugino
c       complex vc(6)          : input vector
c       complex gc(2)          : coupling constants
c       real    vmass          : mass  of output vector
c       real    vwidth         : width of output vector
c
c output:
c       complex jvirog(6)      : vector current
c
c- by Yoshitaro Takaesu - 2011/01/01
c
      implicit none
      double complex fi(6),fo(6),vc(6),gc(2),jvirog1(6)
      double complex sql(2,2),sqr(2,2),svl(2,2),svr(2,2),vgl(0:3,2,2)
      double complex vgr(0:3,2,2),gvl(0:3,2,2),gvr(0:3,2,2),qvl(2,2)
      double complex qvr(2,2),qvcl(2,2),qvcr(2,2),vql(2,2),vqr(2,2)
      double complex lol(0:3,2,2),gvcl(0:3,2,2),gvcr(0:3,2,2)
      double complex lor(0:3,2,2),flo(0:3,4),d
      double precision q(0:3),vmass,vwidth,q2,vm2
      integer i,j,k

      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
      double complex cZero,ci
      parameter( cZero = ( 0.0d0, 0.0d0 ), ci = ( 0.0d0, 1.0d0 ) )
c
      jvirog1(1) = cZero
      jvirog1(2) = cZero
      jvirog1(3) = cZero
      jvirog1(4) = cZero
      jvirog1(5) = vc(5) +fo(5) - fi(5)
      jvirog1(6) = vc(6) +fo(6) - fi(6)

      q(0) = dble( jvirog1(5))
      q(1) = dble( jvirog1(6))
      q(2) = dimag(jvirog1(6))
      q(3) = dimag(jvirog1(5))
      q2  = q(0)**2 - q(1)**2 - q(2)**2 - q(3)**2
      vm2 = vmass**2

c sq
      sql(1,1) = q(0)-q(3)
      sql(1,2) = -(q(1)-ci*q(2))
      sql(2,1) = -(q(1)+ci*q(2))
      sql(2,2) = q(0)+q(3)

      sqr(1,1) = q(0)+q(3)
      sqr(1,2) = q(1)-ci*q(2)
      sqr(2,1) = q(1)+ci*q(2)
      sqr(2,2) = q(0)-q(3)

c svc
      svl(1,1) = vc(1)-vc(4)
      svl(1,2) = -(vc(2)-ci*vc(3))
      svl(2,1) = -(vc(2)+ci*vc(3))
      svl(2,2) = vc(1)+vc(4)
      
      svr(1,1) = vc(1)+vc(4)
      svr(1,2) = vc(2)-ci*vc(3)
      svr(2,1) = vc(2)+ci*vc(3)
      svr(2,2) = vc(1)-vc(4)

c gamma^mu*svc
      gvl(0,1,1) = svr(1,1)
      gvl(0,1,2) = svr(1,2)
      gvl(0,2,1) = svr(2,1)
      gvl(0,2,2) = svr(2,2)
      gvl(1,1,1) = svr(2,1)
      gvl(1,1,2) = svr(2,2)
      gvl(1,2,1) = svr(1,1)
      gvl(1,2,2) = svr(1,2)
      gvl(2,1,1) = -ci*svr(2,1)
      gvl(2,1,2) = -ci*svr(2,2)
      gvl(2,2,1) = ci*svr(1,1)
      gvl(2,2,2) = ci*svr(1,2)
      gvl(3,1,1) = svr(1,1)
      gvl(3,1,2) = svr(1,2)
      gvl(3,2,1) = -svr(2,1)
      gvl(3,2,2) = -svr(2,2)

      gvr(0,1,1) = svl(1,1)
      gvr(0,1,2) = svl(1,2)
      gvr(0,2,1) = svl(2,1)
      gvr(0,2,2) = svl(2,2)
      gvr(1,1,1) = -svl(2,1)
      gvr(1,1,2) = -svl(2,2)
      gvr(1,2,1) = -svl(1,1)
      gvr(1,2,2) = -svl(1,2)
      gvr(2,1,1) = ci*svl(2,1)
      gvr(2,1,2) = ci*svl(2,2)
      gvr(2,2,1) = -ci*svl(1,1)
      gvr(2,2,2) = -ci*svl(1,2)
      gvr(3,1,1) = -svl(1,1)
      gvr(3,1,2) = -svl(1,2)
      gvr(3,2,1) = svl(2,1)
      gvr(3,2,2) = svl(2,2)

c svc*gamma^mu
      vgl(0,1,1) = svl(1,1)
      vgl(0,1,2) = svl(1,2)
      vgl(0,2,1) = svl(2,1)
      vgl(0,2,2) = svl(2,2)
      vgl(1,1,1) = -svl(1,2)
      vgl(1,1,2) = -svl(1,1)
      vgl(1,2,1) = -svl(2,2)
      vgl(1,2,2) = -svl(2,1)
      vgl(2,1,1) = -ci*svl(1,2)
      vgl(2,1,2) = ci*svl(1,1)
      vgl(2,2,1) = -ci*svl(2,2)
      vgl(2,2,2) = ci*svl(2,1)
      vgl(3,1,1) = -svl(1,1)
      vgl(3,1,2) = svl(1,2)
      vgl(3,2,1) = -svl(2,1)
      vgl(3,2,2) = svl(2,2)

      vgr(0,1,1) = svr(1,1)
      vgr(0,1,2) = svr(1,2)
      vgr(0,2,1) = svr(2,1)
      vgr(0,2,2) = svr(2,2)
      vgr(1,1,1) = svr(1,2)
      vgr(1,1,2) = svr(1,1)
      vgr(1,2,1) = svr(2,2)
      vgr(1,2,2) = svr(2,1)
      vgr(2,1,1) = ci*svr(1,2)
      vgr(2,1,2) = -ci*svr(1,1)
      vgr(2,2,1) = ci*svr(2,2)
      vgr(2,2,2) = -ci*svr(2,1)
      vgr(3,1,1) = svr(1,1)
      vgr(3,1,2) = -svr(1,2)
      vgr(3,2,1) = svr(2,1)
      vgr(3,2,2) = -svr(2,2)

c sq*svc
      qvl(1,1) = sql(1,1)*svr(1,1)+sql(1,2)*svr(2,1)
      qvl(1,2) = sql(1,1)*svr(1,2)+sql(1,2)*svr(2,2) 
      qvl(2,1) = sql(2,1)*svr(1,1)+sql(2,2)*svr(2,1) 
      qvl(2,2) = sql(2,1)*svr(1,2)+sql(2,2)*svr(2,2)
      
      qvr(1,1) = sqr(1,1)*svl(1,1)+sqr(1,2)*svl(2,1)
      qvr(1,2) = sqr(1,1)*svl(1,2)+sqr(1,2)*svl(2,2) 
      qvr(2,1) = sqr(2,1)*svl(1,1)+sqr(2,2)*svl(2,1) 
      qvr(2,2) = sqr(2,1)*svl(1,2)+sqr(2,2)*svl(2,2) 

c svc*sq
      vql(1,1) = svl(1,1)*sqr(1,1)+svl(1,2)*sqr(2,1)
      vql(1,2) = svl(1,1)*sqr(1,2)+svl(1,2)*sqr(2,2)
      vql(2,1) = svl(2,1)*sqr(1,1)+svl(2,2)*sqr(2,1) 
      vql(2,2) = svl(2,1)*sqr(1,2)+svl(2,2)*sqr(2,2)

      vqr(1,1) = svr(1,1)*sql(1,1)+svr(1,2)*sql(2,1)
      vqr(1,2) = svr(1,1)*sql(1,2)+svr(1,2)*sql(2,2) 
      vqr(2,1) = svr(2,1)*sql(1,1)+svr(2,2)*sql(2,1) 
      vqr(2,2) = svr(2,1)*sql(1,2)+svr(2,2)*sql(2,2)

c [gamma^mu,svc]
      do i = 0,3
         do j = 1,2
            do k = 1,2
               gvcl(i,j,k) = gvl(i,j,k) -vgl(i,j,k)
               gvcr(i,j,k) = gvr(i,j,k) -vgr(i,j,k)
            enddo
         enddo
      enddo

      if ( vmass.ne.rZero ) then
c [sq,svc]
         do i = 1,2
            do j = 1,2
               qvcl(i,j) = qvl(i,j) -vql(i,j)
               qvcr(i,j) = qvr(i,j) -vqr(i,j)
            enddo
         enddo

c ([gamma^mu,svc] -q^mu/vmass**2*[sq,svc])
         do i = 0,3
            do j = 1,2
               do k = 1,2
                  lol(i,j,k) = gvcl(i,j,k) -q(i)/vm2*qvcl(j,k)
                  lor(i,j,k) = gvcr(i,j,k) -q(i)/vm2*qvcr(j,k)
               enddo
            enddo
         enddo

c fo*([gamma^mu,svc] -q^mu/vmass**2*[sq,svc])
         do i = 0,3
            flo(i,1) = fo(1)*lol(i,1,1) +fo(2)*lol(i,2,1)
            flo(i,2) = fo(1)*lol(i,1,2) +fo(2)*lol(i,2,2)
            flo(i,3) = fo(3)*lor(i,1,1) +fo(4)*lor(i,2,1)
            flo(i,4) = fo(3)*lor(i,1,2) +fo(4)*lor(i,2,2)
         enddo

         d = rOne/dcmplx( q2-vm2, vmass*vwidth )
c  for the running width, use below instead of the above d.
c         d = rOne/dcmplx( q2-vm2, max(vwidth*q2/vmass,rZero) )

         jvirog1(1) = d*(dconjg(gc(2))*(flo(0,1)*fi(1) +flo(0,2)*fi(2))
     &              +dconjg(gc(1))*(flo(0,3)*fi(3) +flo(0,4)*fi(4)))
         jvirog1(2) = d*(dconjg(gc(2))*(flo(1,1)*fi(1) +flo(1,2)*fi(2))
     &              +dconjg(gc(1))*(flo(1,3)*fi(3) +flo(1,4)*fi(4)))
         jvirog1(3) = d*(dconjg(gc(2))*(flo(2,1)*fi(1) +flo(2,2)*fi(2))
     &              +dconjg(gc(1))*(flo(2,3)*fi(3) +flo(2,4)*fi(4)))
         jvirog1(4) = d*(dconjg(gc(2))*(flo(3,1)*fi(1) +flo(3,2)*fi(2))
     &              +dconjg(gc(1))*(flo(3,3)*fi(3) +flo(3,4)*fi(4)))          
        
      else
c ([gamma^mu,svc] -q^mu/vmass**2*[sq,svc])
         do i = 0,3
            do j = 1,2
               do k = 1,2
                  lol(i,j,k) = gvcl(i,j,k)
                  lor(i,j,k) = gvcr(i,j,k)
               enddo
            enddo
         enddo

c fo*([gamma^mu,svc] -q^mu/vmass**2*[sq,svc])
         do i = 0,3
            flo(i,1) = fo(1)*lol(i,1,1) +fo(2)*lol(i,2,1)
            flo(i,2) = fo(1)*lol(i,1,2) +fo(2)*lol(i,2,2)
            flo(i,3) = fo(3)*lor(i,1,1) +fo(4)*lor(i,2,1)
            flo(i,4) = fo(3)*lor(i,1,2) +fo(4)*lor(i,2,2)
         enddo
      
         d = rOne/q2

         jvirog1(1) = d*(dconjg(gc(2))*(flo(0,1)*fi(1) +flo(0,2)*fi(2))
     &              +dconjg(gc(1))*(flo(0,3)*fi(3) +flo(0,4)*fi(4)))
         jvirog1(2) = d*(dconjg(gc(2))*(flo(1,1)*fi(1) +flo(1,2)*fi(2))
     &              +dconjg(gc(1))*(flo(1,3)*fi(3) +flo(1,4)*fi(4)))
         jvirog1(3) = d*(dconjg(gc(2))*(flo(2,1)*fi(1) +flo(2,2)*fi(2))
     &              +dconjg(gc(1))*(flo(2,3)*fi(3) +flo(2,4)*fi(4)))
         jvirog1(4) = d*(dconjg(gc(2))*(flo(3,1)*fi(1) +flo(3,2)*fi(2))
     &              +dconjg(gc(1))*(flo(3,3)*fi(3) +flo(3,4)*fi(4)))          
         
      end if

      return
      end

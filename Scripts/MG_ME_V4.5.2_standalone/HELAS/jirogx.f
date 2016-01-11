c
c ----------------------------------------------------------------------
c
      subroutine jirogx(fi,fo,gc,vmass,vwidth, jirog)
c
c This subroutine computes an off-shell vector current for the fermion-
c Goldstino-vector boson vertex from an external fermion pair. The vector boson 
c propagator is given in feynman gauge for a massless vector and in 
c unitary gauge for a massive vector.
c
c input:
c       complex fi(6)          : flow-in  goldstino                 |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex gc(2)          : coupling constants                  gvf
c       real    vmass          : mass  of output vector v
c       real    vwidth         : width of output vector v
c
c output:
c       complex jirog(6)         : vector current          j^mu(<fo|v|fi>)
c
c- by Yoshitaro Takaesu - 2010/11/25
c
      implicit none
      double complex  fi(6), fo(6), gc(2), jirog(6),sql(2,2),sqr(2,2)
      double complex  d,sqgl(0:3,2,2),sqgr(0:3,2,2),gsql(0:3,2,2)
      double complex gsqr(0:3,2,2),fqg(0:3,4)
      double precision  q(0:3), vmass, vwidth, q2, vm2
      integer i
      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
      double complex ci, cZero
      parameter( ci = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )
c
      jirog(1) = cZero
      jirog(2) = cZero
      jirog(3) = cZero
      jirog(4) = cZero      
      jirog(5) = fo(5) - fi(5)
      jirog(6) = fo(6) - fi(6)

      q(0) = dble( jirog(5))
      q(1) = dble( jirog(6))
      q(2) = dimag(jirog(6))
      q(3) = dimag(jirog(5))
      q2  = q(0)**2 - q(1)**2 - q(2)**2 - q(3)**2
      vm2 = vmass**2

      sql(1,1) = q(0)-q(3)
      sql(1,2) = -(q(1)-ci*q(2))
      sql(2,1) = -(q(1)+ci*q(2))
      sql(2,2) = q(0)+q(3)

      sqr(1,1) = q(0)+q(3)
      sqr(1,2) = q(1)-ci*q(2)
      sqr(2,1) = q(1)+ci*q(2)
      sqr(2,2) = q(0)-q(3)

      sqgl(0,1,1) = sql(1,1)
      sqgl(0,1,2) = sql(1,2)
      sqgl(0,2,1) = sql(2,1)
      sqgl(0,2,2) = sql(2,2)
      sqgl(1,1,1) = -sql(1,2)
      sqgl(1,1,2) = -sql(1,1)
      sqgl(1,2,1) = -sql(2,2)
      sqgl(1,2,2) = -sql(2,1)
      sqgl(2,1,1) = -ci*sql(1,2)
      sqgl(2,1,2) = ci*sql(1,1)
      sqgl(2,2,1) = -ci*sql(2,2)
      sqgl(2,2,2) = ci*sql(2,1)
      sqgl(3,1,1) = -sql(1,1)
      sqgl(3,1,2) = sql(1,2)
      sqgl(3,2,1) = -sql(2,1)
      sqgl(3,2,2) = sql(2,2)

      sqgr(0,1,1) = sqr(1,1)
      sqgr(0,1,2) = sqr(1,2)
      sqgr(0,2,1) = sqr(2,1)
      sqgr(0,2,2) = sqr(2,2)
      sqgr(1,1,1) = sqr(1,2)
      sqgr(1,1,2) = sqr(1,1)
      sqgr(1,2,1) = sqr(2,2)
      sqgr(1,2,2) = sqr(2,1)
      sqgr(2,1,1) = ci*sqr(1,2)
      sqgr(2,1,2) = -ci*sqr(1,1)
      sqgr(2,2,1) = ci*sqr(2,2)
      sqgr(2,2,2) = -ci*sqr(2,1)
      sqgr(3,1,1) = sqr(1,1)
      sqgr(3,1,2) = -sqr(1,2)
      sqgr(3,2,1) = sqr(2,1)
      sqgr(3,2,2) = -sqr(2,2)

      gsql(0,1,1) = sqr(1,1)
      gsql(0,1,2) = sqr(1,2)
      gsql(0,2,1) = sqr(2,1)
      gsql(0,2,2) = sqr(2,2)
      gsql(1,1,1) = sqr(2,1)
      gsql(1,1,2) = sqr(2,2)
      gsql(1,2,1) = sqr(1,1)
      gsql(1,2,2) = sqr(1,2)
      gsql(2,1,1) = -ci*sqr(2,1)
      gsql(2,1,2) = -ci*sqr(2,2)
      gsql(2,2,1) = ci*sqr(1,1)
      gsql(2,2,2) = ci*sqr(1,2)
      gsql(3,1,1) = sqr(1,1)
      gsql(3,1,2) = sqr(1,2)
      gsql(3,2,1) = -sqr(2,1)
      gsql(3,2,2) = -sqr(2,2)

      gsqr(0,1,1) = sql(1,1)
      gsqr(0,1,2) = sql(1,2)
      gsqr(0,2,1) = sql(2,1)
      gsqr(0,2,2) = sql(2,2)
      gsqr(1,1,1) = -sql(2,1)
      gsqr(1,1,2) = -sql(2,2)
      gsqr(1,2,1) = -sql(1,1)
      gsqr(1,2,2) = -sql(1,2)
      gsqr(2,1,1) = ci*sql(2,1)
      gsqr(2,1,2) = ci*sql(2,2)
      gsqr(2,2,1) = -ci*sql(1,1)
      gsqr(2,2,2) = -ci*sql(1,2)
      gsqr(3,1,1) = -sql(1,1)
      gsqr(3,1,2) = -sql(1,2)
      gsqr(3,2,1) = sql(2,1)
      gsqr(3,2,2) = sql(2,2)

      do i = 0,3
         fqg(i,1) = fo(1)*(sqgl(i,1,1) -gsql(i,1,1))
     &             +fo(2)*(sqgl(i,2,1) -gsql(i,2,1))
         fqg(i,2) = fo(1)*(sqgl(i,1,2) -gsql(i,1,2))
     &             +fo(2)*(sqgl(i,2,2) -gsql(i,2,2))
         fqg(i,3) = fo(3)*(sqgr(i,1,1) -gsqr(i,1,1))
     &             +fo(4)*(sqgr(i,2,1) -gsqr(i,2,1))
         fqg(i,4) = fo(3)*(sqgr(i,1,2) -gsqr(i,1,2))
     &             +fo(4)*(sqgr(i,2,2) -gsqr(i,2,2))
      enddo

      if ( vmass.ne.rZero ) then

         d = -rOne/dcmplx( q2-vm2, vmass*vwidth )
c  for the running width, use below instead of the above d.
c         d = rOne/dcmplx( q2-vm2, max(vwidth*q2/vmass,rZero) )

         jirog(1) = d*(dconjg(gc(2))*(fqg(0,1)*fi(1) +fqg(0,2)*fi(2))
     &              +dconjg(gc(1))*(fqg(0,3)*fi(3) +fqg(0,4)*fi(4)))
         jirog(2) = d*(dconjg(gc(2))*(fqg(1,1)*fi(1) +fqg(1,2)*fi(2))
     &              +dconjg(gc(1))*(fqg(1,3)*fi(3) +fqg(1,4)*fi(4)))
         jirog(3) = d*(dconjg(gc(2))*(fqg(2,1)*fi(1) +fqg(2,2)*fi(2))
     &              +dconjg(gc(1))*(fqg(2,3)*fi(3) +fqg(2,4)*fi(4)))
         jirog(4) = d*(dconjg(gc(2))*(fqg(3,1)*fi(1) +fqg(3,2)*fi(2))
     &              +dconjg(gc(1))*(fqg(3,3)*fi(3) +fqg(3,4)*fi(4)))          
        
      else
         d = -rOne/q2

         jirog(1) = d*(dconjg(gc(2))*(fqg(0,1)*fi(1) +fqg(0,2)*fi(2))
     &              +dconjg(gc(1))*(fqg(0,3)*fi(3) +fqg(0,4)*fi(4)))
         jirog(2) = d*(dconjg(gc(2))*(fqg(1,1)*fi(1) +fqg(1,2)*fi(2))
     &              +dconjg(gc(1))*(fqg(1,3)*fi(3) +fqg(1,4)*fi(4)))
         jirog(3) = d*(dconjg(gc(2))*(fqg(2,1)*fi(1) +fqg(2,2)*fi(2))
     &              +dconjg(gc(1))*(fqg(2,3)*fi(3) +fqg(2,4)*fi(4)))
         jirog(4) = d*(dconjg(gc(2))*(fqg(3,1)*fi(1) +fqg(3,2)*fi(2))
     &              +dconjg(gc(1))*(fqg(3,3)*fi(3) +fqg(3,4)*fi(4)))  
         
         
        
      end if

      return
      end

   

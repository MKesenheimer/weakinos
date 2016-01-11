      subroutine orxxxx(p,rmass,nhel,nsr , ro)
c
c This subroutine computes a Rarita-Schwinger wavefunction of spin-3/2
c fermion with the flowing-IN fermion number.
c
c input:
c       real    p(0:3)           : four-momentum of RS fermion
c       real    rmass            : mass          of RS fermion
c       integer nhel = -3,-1,1,3 : helicity      of RS fermion
c                                  (1- and 1 is forbidden if rmass = 0)
c       integer nsr  = -1 or 1   : +1 for particle, -1 for anti-particle
c
c output:
c       complex ro(18)           : RS fermion wavefunction         |ro>v   
c     
c- by Y.Takaesu - 2011/01/11
c
      implicit none
      double precision p(0:3),rmass
      integer nhel,nsr
      double complex ro(18),fipp(4),fimm(4)

      double complex rc(6,4),ep(4),em(4),e0(4),fop(4),fom(4),chi(2)
      double precision pp,pt2,pt,pzpt,emp, sf(2),sfomeg(2),omega(2),pp3,
     &                 sqp0p3,sqm(0:1)      
      integer i,j,nsv,ip,im,nh

      double precision rZero, rHalf, rOne, rTwo, rThree, sqh,sq2,sq3
      parameter( rZero = 0.0d0, rHalf = 0.5d0 )
      parameter( rOne = 1.0d0, rTwo = 2.0d0, rThree = 3.0d0 )

#ifdef HELAS_CHECK
      double precision p2
      double precision epsi
      parameter( epsi = 2.0d-5 )
      integer stdo
      parameter( stdo = 6 )
#endif
c
#ifdef HELAS_CHECK
      pp = sqrt(p(1)**2+p(2)**2+p(3)**2)
      if ( abs(p(0))+pp.eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : p(0:3) in orxxxx is zero momentum'
      endif
      if ( p(0).le.rZero ) then
         write(stdo,*)
     &        ' helas-error : p(0:3) in orxxxx has non-positive energy'
         write(stdo,*)
     &        '             : p(0) = ',p(0)
      endif
      p2 = (p(0)-pp)*(p(0)+pp)
      if ( abs(p2-rmass**2).gt.p(0)**2*epsi ) then
         write(stdo,*)
     &        ' helas-error : p(0:3) in orxxxx has inappropriate mass'
         write(stdo,*)
     &        '             : p**2 = ',p2,' : rmass**2 = ',rmass**2
      endif
      if (abs(nhel).gt.3 .or. abs(nhel).eq.2 .or. abs(nhel).eq.0 ) then
         write(stdo,*) ' helas-error : nhel in orxxxx is not -3,-1,1,3'
         write(stdo,*) '             : nhel = ',nhel
      endif
      if (abs(nsr).ne.1) then
         write(stdo,*) ' helas-error : nsr in orxxxx is not -1,1'
         write(stdo,*) '             : nsr = ',nsr
      endif
#endif

      sqh = sqrt(rHalf)
      sq2 = sqrt(rTwo)
      sq3 = sqrt(rThree)

      pt2 = p(1)**2 + p(2)**2
      pp = min(p(0),sqrt(pt2+p(3)**2))
      pt = min(pp,sqrt(pt2))

      rc(5,1) = dcmplx(p(0),p(3))*nsr
      rc(6,1) = dcmplx(p(1),p(2))*nsr

      nsv = nsr ! nsv=+1 for final, -1 for initial

      if ( nhel.ge.1 ) then 
c construct eps+
         if ( pp.eq.rZero ) then
            ep(1) = dcmplx( rZero )
            ep(2) = dcmplx( -sqh )
            ep(3) = dcmplx( rZero , nsv*sqh )
            ep(4) = dcmplx( rZero )
         else
            ep(1) = dcmplx( rZero )
            ep(4) = dcmplx( pt/pp*sqh )
            if ( pt.ne.rZero ) then
               pzpt = p(3)/(pp*pt)*sqh
               ep(2) = dcmplx( -p(1)*pzpt , -nsv*p(2)/pt*sqh )
               ep(3) = dcmplx( -p(2)*pzpt ,  nsv*p(1)/pt*sqh )
            else
               ep(2) = dcmplx( -sqh )
               ep(3) = dcmplx( rZero , nsv*sign(sqh,p(3)) )
            endif
         endif
      end if

      if ( nhel.le.-1 ) then 
c construct eps-
         if ( pp.eq.rZero ) then
            em(1) = dcmplx( rZero )
            em(2) = dcmplx( sqh )
            em(3) = dcmplx( rZero , nsv*sqh )
            em(4) = dcmplx( rZero )
         else
            em(1) = dcmplx( rZero )
            em(4) = dcmplx( -pt/pp*sqh )
            if ( pt.ne.rZero ) then
               pzpt = -p(3)/(pp*pt)*sqh
               em(2) = dcmplx( -p(1)*pzpt , -nsv*p(2)/pt*sqh )
               em(3) = dcmplx( -p(2)*pzpt ,  nsv*p(1)/pt*sqh )
            else
               em(2) = dcmplx( sqh )
               em(3) = dcmplx( rZero , nsv*sign(sqh,p(3)) )
            endif
         endif
      end if

      if ( abs(nhel).le.1 ) then  
c construct eps0
         if ( pp.eq.rZero ) then
            e0(1) = dcmplx( rZero )
            e0(2) = dcmplx( rZero )
            e0(3) = dcmplx( rZero )
            e0(4) = dcmplx( rOne )
         else
            emp = p(0)/(rmass*pp)
            e0(1) = dcmplx( pp/rmass )
            e0(4) = dcmplx( p(3)*emp )
            if ( pt.ne.rZero ) then
               e0(2) = dcmplx( p(1)*emp )
               e0(3) = dcmplx( p(2)*emp )
            else
               e0(2) = dcmplx( rZero )
               e0(3) = dcmplx( rZero )
            endif
         end if
      end if

      if ( nhel.ge.-1 ) then
c constract spinor+ 
       nh = nsr

       if ( rmass.ne.rZero ) then

         pp = min(p(0),dsqrt(p(1)**2+p(2)**2+p(3)**2))

         if ( pp.eq.rZero ) then
            
            sqm(0) = dsqrt(abs(rmass)) ! possibility of negative fermion masses
            sqm(1) = sign(sqm(0),rmass) ! possibility of negative fermion masses
            ip = -((1+nh)/2)
            im =  (1-nh)/2
            
            fop(1) = im     * sqm(im)
            fop(2) = ip*nsr * sqm(im)
            fop(3) = im*nsr * sqm(-ip)
            fop(4) = ip     * sqm(-ip)
            
         else
            
            pp = min(p(0),dsqrt(p(1)**2+p(2)**2+p(3)**2))
            sf(1) = dble(1+nsr+(1-nsr)*nh)*rHalf
            sf(2) = dble(1+nsr-(1-nsr)*nh)*rHalf
            omega(1) = dsqrt(p(0)+pp)
            omega(2) = rmass/omega(1)
            ip = (3+nh)/2
            im = (3-nh)/2
            sfomeg(1) = sf(1)*omega(ip)
            sfomeg(2) = sf(2)*omega(im)
            pp3 = max(pp+p(3),rZero)
            chi(1) = dcmplx( dsqrt(pp3*rHalf/pp) )
            if ( pp3.eq.rZero ) then
               chi(2) = dcmplx(-nh )
            else
               chi(2) = dcmplx( nh*p(1) , -p(2) )/dsqrt(rTwo*pp*pp3)
            endif
            
            fop(1) = sfomeg(2)*chi(im)
            fop(2) = sfomeg(2)*chi(ip)
            fop(3) = sfomeg(1)*chi(im)
            fop(4) = sfomeg(1)*chi(ip)

         endif
         
      else
         
         if(p(1).eq.0d0.and.p(2).eq.0d0.and.p(3).lt.0d0) then
            sqp0p3 = 0d0
         else
            sqp0p3 = dsqrt(max(p(0)+p(3),rZero))*nsr
         end if
         chi(1) = dcmplx( sqp0p3 )
         if ( sqp0p3.eq.rZero ) then
            chi(2) = dcmplx(-nhel )*dsqrt(rTwo*p(0))
         else
            chi(2) = dcmplx( nh*p(1), -p(2) )/sqp0p3
         endif
         if ( nh.eq.1 ) then
            fop(1) = chi(1)
            fop(2) = chi(2)
            fop(3) = dcmplx( rZero )
            fop(4) = dcmplx( rZero )
         else
            fop(1) = dcmplx( rZero )
            fop(2) = dcmplx( rZero )
            fop(3) = chi(2)
            fop(4) = chi(1)
         endif
       endif
      endif

      if ( nhel.le.1 ) then
c constract spinor+ 
       nh = -nsr

      if ( rmass.ne.rZero ) then

         pp = min(p(0),dsqrt(p(1)**2+p(2)**2+p(3)**2))

         if ( pp.eq.rZero ) then
            
            sqm(0) = dsqrt(abs(rmass)) ! possibility of negative fermion masses
            sqm(1) = sign(sqm(0),rmass) ! possibility of negative fermion masses
            ip = -((1+nh)/2)
            im =  (1-nh)/2
            
            fom(1) = im     * sqm(im)
            fom(2) = ip*nsr * sqm(im)
            fom(3) = im*nsr * sqm(-ip)
            fom(4) = ip     * sqm(-ip)
            
         else
            
            pp = min(p(0),dsqrt(p(1)**2+p(2)**2+p(3)**2))
            sf(1) = dble(1+nsr+(1-nsr)*nh)*rHalf
            sf(2) = dble(1+nsr-(1-nsr)*nh)*rHalf
            omega(1) = dsqrt(p(0)+pp)
            omega(2) = rmass/omega(1)
            ip = (3+nh)/2
            im = (3-nh)/2
            sfomeg(1) = sf(1)*omega(ip)
            sfomeg(2) = sf(2)*omega(im)
            pp3 = max(pp+p(3),rZero)
            chi(1) = dcmplx( dsqrt(pp3*rHalf/pp) )
            if ( pp3.eq.rZero ) then
               chi(2) = dcmplx(-nh )
            else
               chi(2) = dcmplx( nh*p(1) , -p(2) )/dsqrt(rTwo*pp*pp3)
            endif
            
            fom(1) = sfomeg(2)*chi(im)
            fom(2) = sfomeg(2)*chi(ip)
            fom(3) = sfomeg(1)*chi(im)
            fom(4) = sfomeg(1)*chi(ip)

         endif
         
      else
         
         if(p(1).eq.0d0.and.p(2).eq.0d0.and.p(3).lt.0d0) then
            sqp0p3 = 0d0
         else
            sqp0p3 = dsqrt(max(p(0)+p(3),rZero))*nsr
         end if
         chi(1) = dcmplx( sqp0p3 )
         if ( sqp0p3.eq.rZero ) then
            chi(2) = dcmplx(-nhel )*dsqrt(rTwo*p(0))
         else
            chi(2) = dcmplx( nh*p(1), -p(2) )/sqp0p3
         endif
         if ( nh.eq.1 ) then
            fom(1) = chi(1)
            fom(2) = chi(2)
            fom(3) = dcmplx( rZero )
            fom(4) = dcmplx( rZero )
         else
            fom(1) = dcmplx( rZero )
            fom(2) = dcmplx( rZero )
            fom(3) = chi(2)
            fom(4) = chi(1)
         endif
       endif 
      endif
      
c spin-3/2 fermion wavefunction
      if ( nhel.eq.3 ) then
         do j = 1,4
            do i = 1,4
               rc(i,j) = ep(i)*fop(j)
            end do
         end do
      else if ( nhel.eq.1 ) then
         do j = 1,4
            do i = 1,4
              if     ( pt.eq.rZero .and. p(3).ge.0d0 ) then 
               rc(i,j) =  sq2/sq3*e0(i)*fop(j)
     &                    +rOne/sq3*ep(i)*fom(j)
              elseif ( pt.eq.rZero .and. p(3).lt.0d0 ) then 
               rc(i,j) =  sq2/sq3*e0(i)*fop(j)
     &                    -rOne/sq3*ep(i)*fom(j)
              else
               rc(i,j) =  sq2/sq3*e0(i)*fop(j) 
     &                  +rOne/sq3*ep(i)*fom(j)
     &                   *dcmplx(P(1),-nsr*P(2))/pt  
              endif
            end do
         end do
      else if ( nhel.eq.-1 ) then
         do j = 1,4
            do i = 1,4
              if     ( pt.eq.rZero .and.p(3).ge.0d0 ) then 
               rc(i,j) = rOne/sq3*em(i)*fop(j)
     &                   +sq2/sq3*e0(i)*fom(j)
              elseif ( pt.eq.rZero .and.p(3).lt.0d0 ) then 
               rc(i,j) = rOne/sq3*em(i)*fop(j)
     &                   -sq2/sq3*e0(i)*fom(j)
              else
               rc(i,j) = rOne/sq3*em(i)*fop(j) 
     &                  + sq2/sq3*e0(i)*fom(j)
     &                   *dcmplx(P(1),-nsr*P(2))/pt  
              endif
            end do
         end do
      else
         do j = 1,4
            do i = 1,4
              if     ( pt.eq.rZero .and. p(3).ge.0d0 ) then 
               rc(i,j) =  em(i)*fom(j)
              elseif ( pt.eq.rZero .and. p(3).lt.0d0 ) then 
               rc(i,j) = -em(i)*fom(j)
              else
               rc(i,j) =  em(i)*fom(j)*dcmplx(P(1),-nsr*P(2))/pt  
              endif
            end do
         end do
      end if

      ro(1) = rc(1,1)
      ro(2) = rc(1,2)
      ro(3) = rc(1,3)
      ro(4) = rc(1,4)
      ro(5) = rc(2,1)
      ro(6) = rc(2,2)
      ro(7) = rc(2,3)
      ro(8) = rc(2,4)
      ro(9) = rc(3,1)
      ro(10) = rc(3,2)
      ro(11) = rc(3,3)
      ro(12) = rc(3,4)
      ro(13) = rc(4,1)
      ro(14) = rc(4,2)
      ro(15) = rc(4,3)
      ro(16) = rc(4,4)
      ro(17) = rc(5,1)
      ro(18) = rc(6,1)

      return
      end

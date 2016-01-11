      subroutine hiorxx(fi,ro,gc,smass,swidth , hior)
c
c This subroutine computes an off-shell scalar current from a flowing-IN
c fermion and a flowing-OUT Rarita-Schwinger fermion.
c
c input:
c       complex fi(6)          : flow-in  fermion                   |fi>
c       complex ro(18)         : flow-out RS fermion               v<ro|
c       complex gc(2)          : coupling constants                 GFRS 
c       real    smass          : mass  of OUTPUT scalar s(q)
c       real    swidth         : width of OUTPUT scalar s(q)
c
c output:
c       complex hio(3)         : scalar current           j(<ro|sqv|fi>)
c
c- by Y.Takaesu - 2009/11/17
c
      implicit none
      double complex fi(6),ro(18),hior(3),gc(2),dn,rc(0:3,4),sq1(2,2) 
      double complex sq2(2,2),cImag,rsq(0:3,4),rsqg(4)
      double precision q(0:3),smass,swidth,q2
      integer i

      parameter (cImag = (0d0,1d0))

#ifdef HELAS_CHECK
      double precision rZero, cZero
      parameter( rZero = 0.0d0 )
      double complex cZero
      parameter( cZero = ( 0.0d0, 0.0d0 ) )
      integer stdo
      parameter( stdo = 6 )
#endif
c
#ifdef HELAS_CHECK
      if ( abs(fi(1))+abs(fi(2))+abs(fi(3))+abs(fi(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : fi in hiorxx is zero spinor'
      endif
      if ( abs(fi(5))+abs(fi(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fi in hiorxx has zero momentum'
      endif
       if (abs(ro(1))+abs(ro(2))+abs(ro(3))+abs(ro(4))+abs(ro(5))
     & +abs(ro(6))+abs(ro(7))+abs(ro(8))abs(ro(9))+abs(ro(10))
     & +abs(ro(11))+abs(ro(12))+abs(ro(13))+abs(ro(14))+abs(ro(15))
     & +abs(ro(16)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-warn  : ro in hiorxx is zero spinor'
      endif
      if ( abs(ro(17))+abs(ro(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : ro in hiorxx has zero momentum'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*)
     &        ' helas-error : gc in hiorxx is zero coupling'
      endif
      if ( smass.lt.rZero ) then
         write(stdo,*) ' helas-error : smass in hiorxx is negative'
         write(stdo,*) '             : smass = ',smass
      endif
      if ( swidth.lt.rZero ) then
         write(stdo,*) ' helas-error : swidth in hiorxx is negative'
         write(stdo,*) '             : swidth = ',swidth
      endif
#endif

      hior(2) = ro(17)-fi(5)
      hior(3) = ro(18)-fi(6)

      q(0) = dble( hior(2))
      q(1) = dble( hior(3))
      q(2) = dimag(hior(3))
      q(3) = dimag(hior(2))
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)

#ifdef HELAS_CHECK
      if ( abs(hior(2))+abs(hior(3)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : hior in hiorxx has zero momentum'
      endif
      if ( swidth.eq.rZero .and. q2.eq.smass**2 ) then
         write(stdo,*)
     &        ' helas-error : hior in hiorxx is on smass pole'
         write(stdo,*)
     &        '             : q     = ',q(0),q(1),q(2),q(3)
         write(stdo,*)
     &        '             : abs(q)= ',sqrt(abs(q2))
         hior(1) = cZero
         return
      endif
#endif

      rc(0,1) = ro(1)
      rc(0,2) = ro(2)
      rc(0,3) = ro(3)
      rc(0,4) = ro(4)
      rc(1,1) = ro(5)
      rc(1,2) = ro(6)
      rc(1,3) = ro(7)
      rc(1,4) = ro(8)
      rc(2,1) = ro(9)
      rc(2,2) = ro(10)
      rc(2,3) = ro(11)
      rc(2,4) = ro(12)
      rc(3,1) = ro(13)
      rc(3,2) = ro(14)
      rc(3,3) = ro(15)
      rc(3,4) = ro(16)

      sq1(1,1) = q(0)-q(3)
      sq1(1,2) = -(q(1)-cImag*q(2))
      sq1(2,1) = -(q(1)+cImag*q(2))
      sq1(2,2) = q(0)+q(3)

      sq2(1,1) = q(0)+q(3)
      sq2(1,2) = q(1)-cImag*q(2)
      sq2(2,1) = q(1)+cImag*q(2)
      sq2(2,2) = q(0)-q(3)

      do i=0,3
        rsq(i,1) = rc(i,3)*sq2(1,1)+rc(i,4)*sq2(2,1)
	rsq(i,2) = rc(i,3)*sq2(1,2)+rc(i,4)*sq2(2,2)
	rsq(i,3) = rc(i,1)*sq1(1,1)+rc(i,2)*sq1(2,1)
	rsq(i,4) = rc(i,1)*sq1(1,2)+rc(i,2)*sq1(2,2)
      enddo

      rsqg(1)=rsq(0,3)+rsq(1,4)+cImag*rsq(2,4)+rsq(3,3)
      rsqg(2)=rsq(0,4)+rsq(1,3)-cImag*rsq(2,3)-rsq(3,4)
      rsqg(3)=rsq(0,1)-rsq(1,2)-cImag*rsq(2,2)-rsq(3,1)
      rsqg(4)=rsq(0,2)-rsq(1,1)+cImag*rsq(2,1)+rsq(3,2)

      dn = 1d0/dcmplx( q2-smass**2, smass*swidth )

      hior(1) = dn*( gc(1)*(rsqg(1)*fi(1)+rsqg(2)*fi(2))
     &              +gc(2)*(rsqg(3)*fi(3)+rsqg(4)*fi(4)))
      
      return
      end

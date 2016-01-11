      subroutine hiroxx(ri,fo,gc,smass,swidth , hiro)
c
c This subroutine computes an off-shell scalar current from a flowing-IN
c Rarita-Schwinger fermion and a flowing-OUT fermion.
c
c input:
c       complex ri(18)         : flow-in  RS fermion               |ri>v
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex gc(2)          : coupling constants                 GFRS 
c       real    smass          : mass  of OUTPUT scalar s(q)
c       real    swidth         : width of OUTPUT scalar s(q)
c
c output:
c       complex hiro(3)         : scalar current          j(<fo|svq|ri>)
c  
c- by Y.Takaesu - 2009/10/26
c 

      implicit none
      double complex ri(18),fo(6),hiro(3),gc(2),dn,rc(0:3,4),sq1(2,2)
      double complex sq2(2,2),cImag,sqr(0:3,4),gsqr(4)
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
      if ( abs(fo(1))+abs(fo(2))+abs(fo(3))+abs(fo(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : fo in hiroxx is zero spinor'
      endif
      if ( abs(fo(5))+abs(fo(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : fo in hiroxx has zero momentum'
      endif
       if (abs(ri(1))+abs(ri(2))+abs(ri(3))+abs(ri(4))+abs(ri(5))
     & +abs(ri(6))+abs(ri(7))+abs(ri(8))abs(ri(9))+abs(ri(10))
     & +abs(ri(11))+abs(ri(12))+abs(ri(13))+abs(ri(14))+abs(ri(15))
     & +abs(ri(16)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-warn  : ri in hiroxx is zero spinor'
      endif
      if ( abs(ri(17))+abs(ri(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : ri in hiroxx has zero momentum'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*)
     &        ' helas-error : gc in hiroxx is zero coupling'
      endif
      if ( smass.lt.rZero ) then
         write(stdo,*) ' helas-error : smass in hiroxx is negative'
         write(stdo,*) '             : smass = ',smass
      endif
      if ( swidth.lt.rZero ) then
         write(stdo,*) ' helas-error : swidth in hiroxx is negative'
         write(stdo,*) '             : swidth = ',swidth
      endif
#endif

      hiro(2) = fo(5)-ri(17)
      hiro(3) = fo(6)-ri(18)

      q(0) = dble( hiro(2))
      q(1) = dble( hiro(3))
      q(2) = dimag(hiro(3))
      q(3) = dimag(hiro(2))
      q2 = q(0)**2-(q(1)**2+q(2)**2+q(3)**2)

#ifdef HELAS_CHECK
      if ( abs(hiro(2))+abs(hiro(3)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : hiro in hiroxx has zero momentum'
      endif
      if ( swidth.eq.rZero .and. q2.eq.smass**2 ) then
         write(stdo,*)
     &        ' helas-error : hiro in hiroxx is on smass pole'
         write(stdo,*)
     &        '             : q     = ',q(0),q(1),q(2),q(3)
         write(stdo,*)
     &        '             : abs(q)= ',sqrt(abs(q2))
         hiro(1) = cZero
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

      sq1(1,1) = q(0)-q(3)
      sq1(1,2) = -(q(1)-cImag*q(2))
      sq1(2,1) = -(q(1)+cImag*q(2))
      sq1(2,2) = q(0)+q(3)

      sq2(1,1) = q(0)+q(3)
      sq2(1,2) = q(1)-cImag*q(2)
      sq2(2,1) = q(1)+cImag*q(2)
      sq2(2,2) = q(0)-q(3)

      do i=0,3
        sqr(i,1) = rc(i,3)*sq1(1,1)+rc(i,4)*sq1(1,2)
        sqr(i,2) = rc(i,3)*sq1(2,1)+rc(i,4)*sq1(2,2)
	sqr(i,3) = rc(i,1)*sq2(1,1)+rc(i,2)*sq2(1,2)
	sqr(i,4) = rc(i,1)*sq2(2,1)+rc(i,2)*sq2(2,2)
      enddo

      gsqr(1)=sqr(0,3)-sqr(1,4)+cImag*sqr(2,4)-sqr(3,3)
      gsqr(2)=sqr(0,4)-sqr(1,3)-cImag*sqr(2,3)+sqr(3,4)
      gsqr(3)=sqr(0,1)+sqr(1,2)-cImag*sqr(2,2)+sqr(3,1)
      gsqr(4)=sqr(0,2)+sqr(1,1)+cImag*sqr(2,1)-sqr(3,2)

      dn = -1d0/dcmplx( q2-smass**2, smass*swidth )

      hiro(1) = dn*( dconjg(gc(1))*(gsqr(3)*fo(3)+gsqr(4)*fo(4))
     &              +dconjg(gc(2))*(gsqr(1)*fo(1)+gsqr(2)*fo(2)))
     &          
      
      return
      end

c
c ----------------------------------------------------------------------
c
      subroutine fvvirg(fi,vc1,vc2,gc,fmass,fwidth, fvvirg1)
c
c This subroutine computes an off-shell fermion wave function from the fermion-vector boson-vector boson
c -goldstino coupling.
c
c input:
c       complex fi(6)          : flow-in goldstino                   |fi>
c       complex vc1(6)         : input    vector1                     v1
c       complex vc2(6)         : input    vector2                     v2
c       complex gc(2)          : coupling constants                  gvf 
c
c output:
c       complex fvvirg1(6)      : off-shell fermion wavefunction   |f',v1,v2,fi>
c
c- by Yoshitaro Takaesu - 2011/01/01
c
      implicit none
      double complex fvvirg1(6), fi(6), gc(2), vc1(6),vc2(6)
      double complex svc1l(2,2),svc1r(2,2),svc2l(2,2),svc2r(2,2)
      double complex sssl(2,2),sssr(2,2),ssml(2,2),ssmr(2,2),ssr(2,2)
      double complex ssl(2,2),fvvl(2,2),fvvr(2,2), sfl(2,2),sfr(2,2),d
      integer i,j

      double precision rOne,f(0:3),fmass,fwidth,f2
      parameter( rOne = 1.0d0 )
      double complex ci, cZero
      parameter( ci = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )
c
      fvvirg1(1) = cZero
      fvvirg1(2) = cZero
      fvvirg1(3) = cZero
      fvvirg1(4) = cZero
      fvvirg1(5) = -vc1(5)-vc2(5)+fi(5)
      fvvirg1(6) = -vc1(6)-vc2(6)+fi(6)

      f(0) = dble(fvvirg1(5))
      f(1) = dble(fvvirg1(6))
      f(2) = dimag(fvvirg1(6))
      f(3) = dimag(fvvirg1(5))

      svc1r(1,1) = vc1(1)+vc1(4)
      svc1r(1,2) = vc1(2)-ci*vc1(3)
      svc1r(2,1) = vc1(2)+ci*vc1(3)
      svc1r(2,2) = vc1(1)-vc1(4)

      svc1l(1,1) = vc1(1)-vc1(4)
      svc1l(1,2) = -(vc1(2)-ci*vc1(3))
      svc1l(2,1) = -(vc1(2)+ci*vc1(3))
      svc1l(2,2) = vc1(1)+vc1(4)

      svc2r(1,1) = vc2(1)+vc2(4)
      svc2r(1,2) = vc2(2)-ci*vc2(3)
      svc2r(2,1) = vc2(2)+ci*vc2(3)
      svc2r(2,2) = vc2(1)-vc2(4)

      svc2l(1,1) = vc2(1)-vc2(4)
      svc2l(1,2) = -(vc2(2)-ci*vc2(3))
      svc2l(2,1) = -(vc2(2)+ci*vc2(3))
      svc2l(2,2) = vc2(1)+vc2(4)

      sfl(1,1) = f(0)-f(3)
      sfl(1,2) = -(f(1)-ci*f(2))
      sfl(2,1) = -(f(1)+ci*f(2))
      sfl(2,2) = f(0)+f(3)

      sfr(1,1) = f(0)+f(3)
      sfr(1,2) = f(1)-ci*f(2)
      sfr(2,1) = f(1)+ci*f(2)
      sfr(2,2) = f(0)-f(3)

      sssr(1,1) = svc1r(1,1)*svc2l(1,1)+svc1r(1,2)*svc2l(2,1)
      sssr(1,2) = svc1r(1,1)*svc2l(1,2)+svc1r(1,2)*svc2l(2,2) 
      sssr(2,1) = svc1r(2,1)*svc2l(1,1)+svc1r(2,2)*svc2l(2,1) 
      sssr(2,2) = svc1r(2,1)*svc2l(1,2)+svc1r(2,2)*svc2l(2,2) 
     
      sssl(1,1) = svc1l(1,1)*svc2r(1,1)+svc1l(1,2)*svc2r(2,1)
      sssl(1,2) = svc1l(1,1)*svc2r(1,2)+svc1l(1,2)*svc2r(2,2) 
      sssl(2,1) = svc1l(2,1)*svc2r(1,1)+svc1l(2,2)*svc2r(2,1) 
      sssl(2,2) = svc1l(2,1)*svc2r(1,2)+svc1l(2,2)*svc2r(2,2)

      ssmr(1,1) = svc2r(1,1)*svc1l(1,1)+svc2r(1,2)*svc1l(2,1)
      ssmr(1,2) = svc2r(1,1)*svc1l(1,2)+svc2r(1,2)*svc1l(2,2) 
      ssmr(2,1) = svc2r(2,1)*svc1l(1,1)+svc2r(2,2)*svc1l(2,1) 
      ssmr(2,2) = svc2r(2,1)*svc1l(1,2)+svc2r(2,2)*svc1l(2,2) 
     
      ssml(1,1) = svc2l(1,1)*svc1r(1,1)+svc2l(1,2)*svc1r(2,1)
      ssml(1,2) = svc2l(1,1)*svc1r(1,2)+svc2l(1,2)*svc1r(2,2) 
      ssml(2,1) = svc2l(2,1)*svc1r(1,1)+svc2l(2,2)*svc1r(2,1) 
      ssml(2,2) = svc2l(2,1)*svc1r(1,2)+svc2l(2,2)*svc1r(2,2)

       do i = 1,2
          do j=1,2
            ssr(i,j) = sssr(i,j) - ssmr(i,j)
          enddo
       enddo

        do i = 1,2
          do j=1,2
            ssl(i,j) = sssl(i,j) - ssml(i,j)
          enddo
       enddo

       fvvr(1,1) = sfr(1,1)*ssl(1,1)+sfr(1,2)*ssl(2,1)
       fvvr(1,2) = sfr(1,1)*ssl(1,2)+sfr(1,2)*ssl(2,2)
       fvvr(2,1) = sfr(2,1)*ssl(1,1)+sfr(2,2)*ssl(2,1)
       fvvr(2,2) = sfr(2,1)*ssl(1,2)+sfr(2,2)*ssl(2,2)

       fvvl(1,1) = sfl(1,1)*ssr(1,1)+sfl(1,2)*ssr(2,1)
       fvvl(1,2) = sfl(1,1)*ssr(1,2)+sfl(1,2)*ssr(2,2)
       fvvl(2,1) = sfl(2,1)*ssr(1,1)+sfl(2,2)*ssr(2,1)
       fvvl(2,2) = sfl(2,1)*ssr(1,2)+sfl(2,2)*ssr(2,2)
       
       f2 = f(0)**2-f(1)**2-f(2)**2-f(3)**2
       d = -1d0/(f2-fmass**2+ci*fmass*fwidth)

       fvvirg1(1) = dconjg(gc(1))*d*(fvvl(1,1)*fi(3) +fvvl(1,2)*fi(4))
       fvvirg1(2) = dconjg(gc(1))*d*(fvvl(2,1)*fi(3) +fvvl(2,2)*fi(4))
       fvvirg1(3) = dconjg(gc(1))*d*(fmass*ssr(1,1)*fi(3)
     &              +fmass*ssr(1,2)*fi(4))
       fvvirg1(4) = dconjg(gc(1))*d*(fmass*ssr(2,1)*fi(3)
     &              +fmass*ssr(2,2)*fi(4))

      if ( gc(2).ne.cZero ) then
          fvvirg1(1) = fvvirg1(1) +dconjg(gc(2))*d*(fmass*ssl(1,1)*fi(1) 
     &                 +fmass*ssl(1,2)*fi(2))
          fvvirg1(2) = fvvirg1(2) +dconjg(gc(2))*d*(
     &                 fmass*ssl(2,1)*fi(1) +fmass*ssl(2,2)*fi(2))
          fvvirg1(3) = fvvirg1(3) +dconjg(gc(2))*d*(fvvr(1,1)*fi(1)
     &                 +fvvr(1,2)*fi(2))
          fvvirg1(4) = fvvirg1(4) +dconjg(gc(2))*d*(fvvr(2,1)*fi(1)
     &                 +fvvr(2,2)*fi(2))     
      end if

     
      return
      end

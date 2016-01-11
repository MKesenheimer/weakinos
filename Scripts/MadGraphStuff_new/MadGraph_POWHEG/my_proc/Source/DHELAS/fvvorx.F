      subroutine fvvorx(ro,vc1,vc2,gc,fmass,fwidth , fvvor)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-OUT external Rarita-Schwinger fermion and a vector boson pair.
c
c
c input:
c       complex ro(18)         : flow-out Ratrita-Schwinger fermion  <ro|
c       complex vc1(6)         : input    vector                      v1
c       complex vc2(6)         : input    vector                      v2
c       complex gc(2)          : coupling constants                  gfrvv
c       real    fmass          : mass  of OUTPUT fermion f'
c       real    fwidth         : width of OUTPUT fermion f'
c
c output:
c       complex fvvor(6)        : off-shell fermion             <ro,v1,v2,f'|
c
c- by Y.Takaesu - 2009/11/04
c
      implicit none
      double complex ro(18),vc1(6),vc2(6),gc(2),fvvor(6),d
      double precision fmass,fwidth

      double precision rZero, rOne
      parameter( rZero = 0.0d0, rOne = 1.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )

      integer i
      double precision pf(0:3),pf2,sabsri,pv(0:3)
      double complex svcl(2,2),svcr(2,2),spvl(2,2),spvr(2,2),rc(0:3,4)
      double complex eprc(0:3,4),geprc(4),sssl(2,2),sssr(2,2)
      double complex spfl(2,2),spfr(2,2)

#ifdef HELAS_CHECK
      integer stdo
      parameter( stdo = 6 )
#endif

#ifdef HELAS_CHECK
      if ( abs(ro(1))+abs(ro(2))+abs(ro(3))+abs(ro(4))+abs(ro(5))
     & +abs(ro(6))+abs(ro(7))+abs(ro(8))abs(ro(9))+abs(ro(10))
     & +abs(ro(11))+abs(ro(12))+abs(ro(13))+abs(ro(14))+abs(ro(15))
     & +abs(ro(16)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : ro in fvvorx is zero spinor'
      endif
      if ( abs(ro(17))+abs(ro(18)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : ro in fvvorx has zero momentum'
      endif
      if ( abs(vc1(1))+abs(vc1(2)+abs(vc1(3))+abs(vc1(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : vc1 in fvvorx is zero vector'
      endif
      if ( abs(vc1(5))+abs(vc1(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : vc1 in fvvorx has zero momentum'
      endif
      if ( abs(vc2(1))+abs(vc2(2)+abs(vc2(3))+abs(vc2(4)).eq.rZero ) then
         write(stdo,*) ' helas-warn  : vc1 in fvvorx is zero vector'
      endif
      if ( abs(vc2(5))+abs(vc2(6)).eq.rZero ) then
         write(stdo,*)
     &        ' helas-error : vc2 in fvvorx has zero momentum'
      endif
      if ( gc(1).eq.cZero .and. gc(2).eq.cZero ) then
         write(stdo,*)
     &        ' helas-error : gc in fvvorx is zero coupling'
      endif
      if ( fmass.lt.rZero ) then
         write(stdo,*) ' helas-error : fmass in fvvorx is negative'
         write(stdo,*) '             : fmass = ',fmass
      endif
      if ( fwidth.lt.rZero ) then
         write(stdo,*) ' helas-error : fwidth in fvvorx is negative'
         write(stdo,*) '             : fwidth = ',fwidth
      endif
#endif

      fvvor(5) = ro(17)+vc1(5)+vc2(5)
      fvvor(6) = ro(18)+vc1(6)+vc2(6)

      pf(0) = dble( fvvor(5))
      pf(1) = dble( fvvor(6))
      pf(2) = dimag(fvvor(6))
      pf(3) = dimag(fvvor(5))
      pf2 = pf(0)**2-(pf(1)**2+pf(2)**2+pf(3)**2)

#ifdef HELAS_CHECK
      if ( abs(fvvor(5))+abs(fvvor(6)).eq.rZero ) then
          write(stdo,*)
     &        ' helas-error : fvvor in fvvorx has zero momentum'
       endif
       if ( fwidth.eq.rZero .and. pf2.eq.fmass**2 ) then
          write(stdo,*)
     &         ' helas-error : fvvor in fvvorx is on fmass pole'
          write(stdo,*)
     &         '             : p     = ',pf(0),pf(1),pf(2),pf(3)
          write(stdo,*)
     &         '             : abs(p)= ',sqrt(abs(pf2))
         fvvor(1) = cZero
         fvvor(2) = cZero
         fvvor(3) = cZero
         fvvor(4) = cZero
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

      fvvor(1) = cZero
      fvvor(2) = cZero
      fvvor(3) = cZero
      fvvor(4) = cZero

      svcl(1,1) = vc2(1)+vc2(4)
      svcl(1,2) = vc2(2)-cImag*vc2(3)
      svcl(2,1) = vc2(2)+cImag*vc2(3)
      svcl(2,2) = vc2(1)-vc2(4)

      svcr(1,1) = vc2(1)-vc2(4)
      svcr(1,2) = -(vc2(2)-cImag*vc2(3))
      svcr(2,1) = -(vc2(2)+cImag*vc2(3))
      svcr(2,2) = vc2(1)+vc2(4)
        
      spvl(1,1) = vc1(1)-vc1(4)
      spvl(1,2) = -(vc1(2)-cImag*vc1(3))
      spvl(2,1) = -(vc1(2)+cImag*vc1(3))
      spvl(2,2) = vc1(1)+vc1(4)

      spvr(1,1) = vc1(1)+vc1(4)
      spvr(1,2) = vc1(2)-cImag*vc1(3)
      spvr(2,1) = vc1(2)+cImag*vc1(3)
      spvr(2,2) = vc1(1)-vc1(4)

      spfl(1,1) = pf(0)+pf(3)
      spfl(1,2) = pf(1)-cImag*pf(2)
      spfl(2,1) = pf(1)+cImag*pf(2)
      spfl(2,2) = pf(0)-pf(3)

      spfr(1,1) = pf(0)-pf(3)
      spfr(1,2) = -(pf(1)-cImag*pf(2))
      spfr(2,1) = -(pf(1)+cImag*pf(2))
      spfr(2,2) = pf(0)+pf(3)

      d = -rOne/dcmplx(pf2-fmass**2,fmass*fwidth)

      if ( gc(1).ne.cZero ) then
        
       sssl(1,1) = svcl(1,1)*spvl(1,1)+svcl(1,2)*spvl(2,1)
       sssl(1,2) = svcl(1,1)*spvl(1,2)+svcl(1,2)*spvl(2,2) 
       sssl(2,1) = svcl(2,1)*spvl(1,1)+svcl(2,2)*spvl(2,1) 
       sssl(2,2) = svcl(2,1)*spvl(1,2)+svcl(2,2)*spvl(2,2) 
        
       do i=0,3
        eprc(i,1) = cZero
        eprc(i,2) = cZero
        eprc(i,3) = sssl(1,1)*rc(i,3)+sssl(2,1)*rc(i,4)
        eprc(i,4) = sssl(1,2)*rc(i,3)+sssl(2,2)*rc(i,4)
       enddo

       geprc(1) = eprc(0,3)+eprc(1,4)+cImag*eprc(2,4)+eprc(3,3)
       geprc(2) = eprc(0,4)+eprc(1,3)-cImag*eprc(2,3)-eprc(3,4)
       geprc(3) = cZero
       geprc(4) = cZero

       fvvor(1) = fvvor(1)-d*gc(1)*fmass*geprc(1)
       fvvor(2) = fvvor(2)-d*gc(1)*fmass*geprc(2)
       fvvor(3) = fvvor(3)-d*gc(1)*(spfr(1,1)*geprc(1)+spfr(2,1)
     &            *geprc(2))
       fvvor(4) = fvvor(4)-d*gc(1)*(spfr(1,2)*geprc(1)+spfr(2,2)
     &            *geprc(2))

       sssl(1,1) = spvr(1,1)*svcr(1,1)+spvr(1,2)*svcr(2,1)
       sssl(1,2) = spvr(1,1)*svcr(1,2)+spvr(1,2)*svcr(2,2) 
       sssl(2,1) = spvr(2,1)*svcr(1,1)+spvr(2,2)*svcr(2,1) 
       sssl(2,2) = spvr(2,1)*svcr(1,2)+spvr(2,2)*svcr(2,2) 
        
       do i=0,3
        eprc(i,1) = cZero
        eprc(i,2) = cZero
        eprc(i,3) = sssl(1,1)*rc(i,3)+sssl(2,1)*rc(i,4)
        eprc(i,4) = sssl(1,2)*rc(i,3)+sssl(2,2)*rc(i,4)
       enddo

       geprc(1) = eprc(0,3)+eprc(1,4)+cImag*eprc(2,4)+eprc(3,3)
       geprc(2) = eprc(0,4)+eprc(1,3)-cImag*eprc(2,3)-eprc(3,4)
       geprc(3) = cZero
       geprc(4) = cZero

       fvvor(1) = fvvor(1)+d*gc(1)*fmass*geprc(1)
       fvvor(2) = fvvor(2)+d*gc(1)*fmass*geprc(2)
       fvvor(3) = fvvor(3)+d*gc(1)*(spfr(1,1)*geprc(1)+spfr(2,1)
     &            *geprc(2))
       fvvor(4) = fvvor(4)+d*gc(1)*(spfr(1,2)*geprc(1)+spfr(2,2)
     &            *geprc(2))

     
      endif

      if ( gc(2).ne.cZero ) then

       sssr(1,1) = svcr(1,1)*spvr(1,1)+svcr(1,2)*spvr(2,1)
       sssr(1,2) = svcr(1,1)*spvr(1,2)+svcr(1,2)*spvr(2,2) 
       sssr(2,1) = svcr(2,1)*spvr(1,1)+svcr(2,2)*spvr(2,1) 
       sssr(2,2) = svcr(2,1)*spvr(1,2)+svcr(2,2)*spvr(2,2)

       do i=0,3
        eprc(i,1) = sssr(1,1)*rc(i,1)+sssr(2,1)*rc(i,2)
        eprc(i,2) = sssr(1,2)*rc(i,1)+sssr(2,2)*rc(i,2)
        eprc(i,3) = cZero
        eprc(i,4) = cZero
       enddo

       geprc(1) = cZero
       geprc(2) = cZero
       geprc(3) = eprc(0,1)-eprc(1,2)-cImag*eprc(2,2)-eprc(3,1)
       geprc(4) = eprc(0,2)-eprc(1,1)+cImag*eprc(2,1)+eprc(3,2)  

       fvvor(1) = fvvor(1)-d*gc(2)*(spfl(1,1)*geprc(3)+spfl(2,1)
     &            *geprc(4))
       fvvor(2) = fvvor(2)-d*gc(2)*(spfl(1,2)*geprc(3)+spfl(2,2)
     &            *geprc(4))
       fvvor(3) = fvvor(3)-d*gc(2)*fmass*geprc(3)
       fvvor(4) = fvvor(4)-d*gc(2)*fmass*geprc(4)

       sssr(1,1) = spvl(1,1)*svcl(1,1)+spvl(1,2)*svcl(2,1)
       sssr(1,2) = spvl(1,1)*svcl(1,2)+spvl(1,2)*svcl(2,2) 
       sssr(2,1) = spvl(2,1)*svcl(1,1)+spvl(2,2)*svcl(2,1) 
       sssr(2,2) = spvl(2,1)*svcl(1,2)+spvl(2,2)*svcl(2,2)

       do i=0,3
        eprc(i,1) = sssr(1,1)*rc(i,1)+sssr(2,1)*rc(i,2)
        eprc(i,2) = sssr(1,2)*rc(i,1)+sssr(2,2)*rc(i,2)
        eprc(i,3) = cZero
        eprc(i,4) = cZero
       enddo

       geprc(1) = cZero
       geprc(2) = cZero
       geprc(3) = eprc(0,1)-eprc(1,2)-cImag*eprc(2,2)-eprc(3,1)
       geprc(4) = eprc(0,2)-eprc(1,1)+cImag*eprc(2,1)+eprc(3,2)  

       fvvor(1) = fvvor(1)+d*gc(2)*(spfl(1,1)*geprc(3)+spfl(2,1)
     &            *geprc(4))
       fvvor(2) = fvvor(2)+d*gc(2)*(spfl(1,2)*geprc(3)+spfl(2,2)
     &            *geprc(4))
       fvvor(3) = fvvor(3)+d*gc(2)*fmass*geprc(3)
       fvvor(4) = fvvor(4)+d*gc(2)*fmass*geprc(4)


      endif   

      return
      end

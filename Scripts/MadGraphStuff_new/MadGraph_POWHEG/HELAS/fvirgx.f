c
c ----------------------------------------------------------------------
c
      subroutine fvirgx(fi,vc,gc,fmass,fwidth, fvirg)
c
c This subroutine computes an off-shell fermion wavefunction from a
c flowing-in external goldstino and a vector boson, for the fermion-boson-
c Goldstino vertex.
c
c input:
c       complex fi(6)          : flow-in goldstino                   <fi|
c       complex vc(6)          : input    vector                      v
c       complex gc(2)          : coupling constants                  gvf
c       real    fmass          : mass  of output fermion f'
c       real    fwidth         : width of output fermion f'
c
c output:
c       complex fvirg(6)         : off-shell fermion             <fo,v,f'|
c
c- by Yoshitaro Takaesu - 2010/11/17
c
      
      implicit none
      double complex  fi(6), vc(6), gc(2), fvirg(6), d, svcl(2,2)
      double complex  svcr(2,2), sql(2,2), sqr(2,2), spfl(2,2)
      double complex  spfr(2,2), qvcl(2,2), qvcr(2,2), vcql(2,2)
      double complex vcqr(2,2), pqvl(2,2), pqvr(2,2), pvql(2,2)
      double complex pvqr(2,2)
      double precision  fmass, fwidth
      double precision  pf(0:3), q(0:3), pf2

      double precision rOne
      parameter( rOne = 1.0d0 )
      double complex ci, cZero
      parameter( ci = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )
c
      fvirg(1) = cZero
      fvirg(2) = cZero
      fvirg(3) = cZero
      fvirg(4) = cZero
      fvirg(5) = fi(5) - vc(5)
      fvirg(6) = fi(6) - vc(6)

      q(0) = dble( vc(5))
      q(1) = dble( vc(6))
      q(2) = dimag(vc(6))
      q(3) = dimag(vc(5))

      pf(0) = dble( fvirg(5))
      pf(1) = dble( fvirg(6))
      pf(2) = dimag(fvirg(6))
      pf(3) = dimag(fvirg(5))
      pf2 = pf(0)**2 - pf(1)**2 - pf(2)**2 - pf(3)**2

      svcl(1,1) = vc(1)-vc(4)
      svcl(1,2) = -(vc(2)-ci*vc(3))
      svcl(2,1) = -(vc(2)+ci*vc(3))
      svcl(2,2) = vc(1)+vc(4)

      svcr(1,1) = vc(1)+vc(4)
      svcr(1,2) = vc(2)-ci*vc(3)
      svcr(2,1) = vc(2)+ci*vc(3)
      svcr(2,2) = vc(1)-vc(4)
        
      sql(1,1) = q(0)-q(3)
      sql(1,2) = -(q(1)-ci*q(2))
      sql(2,1) = -(q(1)+ci*q(2))
      sql(2,2) = q(0)+q(3)

      sqr(1,1) = q(0)+q(3)
      sqr(1,2) = q(1)-ci*q(2)
      sqr(2,1) = q(1)+ci*q(2)
      sqr(2,2) = q(0)-q(3)

      spfl(1,1) = pf(0)-pf(3)
      spfl(1,2) = -(pf(1)-ci*pf(2))
      spfl(2,1) = -(pf(1)+ci*pf(2))
      spfl(2,2) = pf(0)+pf(3)

      spfr(1,1) = pf(0)+pf(3)
      spfr(1,2) = pf(1)-ci*pf(2)
      spfr(2,1) = pf(1)+ci*pf(2)
      spfr(2,2) = pf(0)-pf(3)

      d = -rOne/dcmplx(pf2-fmass**2,fmass*fwidth)
      
      qvcl(1,1) = sql(1,1)*svcr(1,1)+sql(1,2)*svcr(2,1)
      qvcl(1,2) = sql(1,1)*svcr(1,2)+sql(1,2)*svcr(2,2) 
      qvcl(2,1) = sql(2,1)*svcr(1,1)+sql(2,2)*svcr(2,1) 
      qvcl(2,2) = sql(2,1)*svcr(1,2)+sql(2,2)*svcr(2,2)

      qvcr(1,1) = sqr(1,1)*svcl(1,1)+sqr(1,2)*svcl(2,1)
      qvcr(1,2) = sqr(1,1)*svcl(1,2)+sqr(1,2)*svcl(2,2) 
      qvcr(2,1) = sqr(2,1)*svcl(1,1)+sqr(2,2)*svcl(2,1) 
      qvcr(2,2) = sqr(2,1)*svcl(1,2)+sqr(2,2)*svcl(2,2)

      vcql(1,1) = svcl(1,1)*sqr(1,1)+svcl(1,2)*sqr(2,1)
      vcql(1,2) = svcl(1,1)*sqr(1,2)+svcl(1,2)*sqr(2,2) 
      vcql(2,1) = svcl(2,1)*sqr(1,1)+svcl(2,2)*sqr(2,1) 
      vcql(2,2) = svcl(2,1)*sqr(1,2)+svcl(2,2)*sqr(2,2)
      
      vcqr(1,1) = svcr(1,1)*sql(1,1)+svcr(1,2)*sql(2,1)
      vcqr(1,2) = svcr(1,1)*sql(1,2)+svcr(1,2)*sql(2,2) 
      vcqr(2,1) = svcr(2,1)*sql(1,1)+svcr(2,2)*sql(2,1) 
      vcqr(2,2) = svcr(2,1)*sql(1,2)+svcr(2,2)*sql(2,2)

      pqvl(1,1) = spfl(1,1)*qvcr(1,1)+spfl(1,2)*qvcr(2,1)
      pqvl(1,2) = spfl(1,1)*qvcr(1,2)+spfl(1,2)*qvcr(2,2) 
      pqvl(2,1) = spfl(2,1)*qvcr(1,1)+spfl(2,2)*qvcr(2,1) 
      pqvl(2,2) = spfl(2,1)*qvcr(1,2)+spfl(2,2)*qvcr(2,2)

      pqvr(1,1) = spfr(1,1)*qvcl(1,1)+spfr(1,2)*qvcl(2,1)
      pqvr(1,2) = spfr(1,1)*qvcl(1,2)+spfr(1,2)*qvcl(2,2) 
      pqvr(2,1) = spfr(2,1)*qvcl(1,1)+spfr(2,2)*qvcl(2,1) 
      pqvr(2,2) = spfr(2,1)*qvcl(1,2)+spfr(2,2)*qvcl(2,2)

      pvql(1,1) = spfl(1,1)*vcqr(1,1)+spfl(1,2)*vcqr(2,1)
      pvql(1,2) = spfl(1,1)*vcqr(1,2)+spfl(1,2)*vcqr(2,2) 
      pvql(2,1) = spfl(2,1)*vcqr(1,1)+spfl(2,2)*vcqr(2,1) 
      pvql(2,2) = spfl(2,1)*vcqr(1,2)+spfl(2,2)*vcqr(2,2)

      pvqr(1,1) = spfr(1,1)*vcql(1,1)+spfr(1,2)*vcql(2,1)
      pvqr(1,2) = spfr(1,1)*vcql(1,2)+spfr(1,2)*vcql(2,2) 
      pvqr(2,1) = spfr(2,1)*vcql(1,1)+spfr(2,2)*vcql(2,1) 
      pvqr(2,2) = spfr(2,1)*vcql(1,2)+spfr(2,2)*vcql(2,2)

      fvirg(1) = d*dconjg(gc(1))*( (pqvl(1,1)-pvql(1,1))*fi(3)
     &     +(pqvl(1,2)-pvql(1,2))*fi(4) )
      fvirg(2) = fvirg(2)+d*dconjg(gc(1))*( (pqvl(2,1)-pvql(2,1))*fi(3)
     &     +(pqvl(2,2)-pvql(2,2))*fi(4) )
      fvirg(3) = fvirg(3)+d*dconjg(gc(1))*( fmass*(qvcr(1,1)-vcqr(1,1))
     &     *fi(3) +fmass*(qvcr(1,2)-vcqr(1,2))*fi(4) )
      fvirg(4) = fvirg(4)+d*dconjg(gc(1))*( fmass*(qvcr(2,1)-vcqr(2,1))
     &     *fi(3) +fmass*(qvcr(2,2)-vcqr(2,2))*fi(4) )
      
     

      if (gc(2) .ne. 0d0) then
         fvirg(1) =fvirg(1)+d*dconjg(gc(2))*(fmass*(qvcl(1,1)-vcql(1,1))
     &        *fi(1) +fmass*(qvcl(1,2)-vcql(1,2))*fi(2) )
         fvirg(2) =fvirg(2)+d*dconjg(gc(2))*(fmass*(qvcl(2,1)-vcql(2,1))
     &        *fi(1) +fmass*(qvcl(2,2)-vcql(2,2))*fi(2) )
         fvirg(3) =fvirg(3)+d*dconjg(gc(2))*((pqvr(1,1)-pvqr(1,1))*fi(1)
     &        +(pqvr(1,2)-pvqr(1,2))*fi(2) )
         fvirg(4) =fvirg(4)+d*dconjg(gc(2))*((pqvr(2,1)-pvqr(2,1))*fi(1)
     &        +(pqvr(2,2)-pvqr(2,2))*fi(2) )
      endif
      
      return          
      end

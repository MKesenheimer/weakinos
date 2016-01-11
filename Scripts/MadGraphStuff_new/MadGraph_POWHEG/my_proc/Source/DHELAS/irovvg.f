c
c ----------------------------------------------------------------------
c
      subroutine irovvg(fi,fo,vc1,vc2,gc,  vertex)
c
c This subroutine computes an amplitude of the fermion-fermion-vector
c SUSY Goldstino coupling.
c
c input:
c       complex fi(6)          : flow-in  goldstino                   |fi>
c       complex fo(6)          : flow-out fermion                   <fo|
c       complex vc1(6)         : input    vector1                     v1
c       complex vc2(6)         : input    vector2                     v2
c       complex gc(2)          : coupling constants                  gvf 
c
c output:
c       complex vertex         : amplitude                     <fo|v1v2|fi>
c
c- by Yoshitaro Takaesu - 2011/01/01
c
      implicit none
      double complex  fi(6), fo(6), gc(2), vc1(6),vc2(6), vertex
      double complex svcl(2,2),svcr(2,2),spvl(2,2),spvr(2,2),sssl(2,2)
      double complex sssr(2,2),ssml(2,2),ssmr(2,2),ssr(2,2),ssl(2,2)
      integer i,j

      double precision rOne
      parameter( rOne = 1.0d0 )
      double complex cImag, cZero
      parameter( cImag = ( 0.0d0, 1.0d0 ), cZero = ( 0.0d0, 0.0d0 ) )
c    
      svcl(1,1) = vc1(1)+vc1(4)
      svcl(1,2) = vc1(2)-cImag*vc1(3)
      svcl(2,1) = vc1(2)+cImag*vc1(3)
      svcl(2,2) = vc1(1)-vc1(4)

      svcr(1,1) = vc1(1)-vc1(4)
      svcr(1,2) = -(vc1(2)-cImag*vc1(3))
      svcr(2,1) = -(vc1(2)+cImag*vc1(3))
      svcr(2,2) = vc1(1)+vc1(4)
        
      spvl(1,1) = vc2(1)-vc2(4)
      spvl(1,2) = -(vc2(2)-cImag*vc2(3))
      spvl(2,1) = -(vc2(2)+cImag*vc2(3))
      spvl(2,2) = vc2(1)+vc2(4)

      spvr(1,1) = vc2(1)+vc2(4)
      spvr(1,2) = vc2(2)-cImag*vc2(3)
      spvr(2,1) = vc2(2)+cImag*vc2(3)
      spvr(2,2) = vc2(1)-vc2(4)

       sssr(1,1) = svcl(1,1)*spvl(1,1)+svcl(1,2)*spvl(2,1)
       sssr(1,2) = svcl(1,1)*spvl(1,2)+svcl(1,2)*spvl(2,2) 
       sssr(2,1) = svcl(2,1)*spvl(1,1)+svcl(2,2)*spvl(2,1) 
       sssr(2,2) = svcl(2,1)*spvl(1,2)+svcl(2,2)*spvl(2,2) 
     
       sssl(1,1) = svcr(1,1)*spvr(1,1)+svcr(1,2)*spvr(2,1)
       sssl(1,2) = svcr(1,1)*spvr(1,2)+svcr(1,2)*spvr(2,2) 
       sssl(2,1) = svcr(2,1)*spvr(1,1)+svcr(2,2)*spvr(2,1) 
       sssl(2,2) = svcr(2,1)*spvr(1,2)+svcr(2,2)*spvr(2,2)

       ssmr(1,1) = spvr(1,1)*svcr(1,1)+spvr(1,2)*svcr(2,1)
       ssmr(1,2) = spvr(1,1)*svcr(1,2)+spvr(1,2)*svcr(2,2) 
       ssmr(2,1) = spvr(2,1)*svcr(1,1)+spvr(2,2)*svcr(2,1) 
       ssmr(2,2) = spvr(2,1)*svcr(1,2)+spvr(2,2)*svcr(2,2) 
       
       ssml(1,1) = spvl(1,1)*svcl(1,1)+spvl(1,2)*svcl(2,1)
       ssml(1,2) = spvl(1,1)*svcl(1,2)+spvl(1,2)*svcl(2,2)
       ssml(2,1) = spvl(2,1)*svcl(1,1)+spvl(2,2)*svcl(2,1) 
       ssml(2,2) = spvl(2,1)*svcl(1,2)+spvl(2,2)*svcl(2,2)

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

      vertex =  dconjg(gc(1))*((fo(3)*ssr(1,1)+fo(4)*ssr(2,1))*fi(3)
     &          + (fo(3)*ssr(1,2)+fo(4)*ssr(2,2))*fi(4))

      if ( gc(2).ne.cZero ) then
         vertex = vertex
     &        +dconjg(gc(2))*((fo(1)*ssl(1,1) + fo(2)*ssl(2,1))*fi(1)
     &        +(fo(1)*ssl(1,2) + fo(2)*ssl(2,2))*fi(2))
      end if

      return
      end

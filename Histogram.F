c############### generate histograms ###################################
      subroutine histogram(name,nbins,xm,xp,xval,weight)
        implicit none
#include "hist.h"
        character*(*) name
        integer k,nbins
        double precision xm,xp
        double precision xval,weight
        integer st,ind,iun,ncalls(nhistsmax)
        integer hits(nhistsmax,0:nbinsmax)
        double precision x(nhistsmax,0:nbinsmax)
        double precision ysum(nhistsmax,0:nbinsmax)
        double precision ymean(nhistsmax,0:nbinsmax)
        double precision yerr(nhistsmax,0:nbinsmax)
        double precision ysqsum(nhistsmax,0:nbinsmax)
        save x, ysum, ymean, yerr, ysqsum
        save hits
        data ncalls/nhistsmax * 0/
        save ncalls
        logical init(nhistsmax)
        data init/nhistsmax * .true./
        
        call index(name,ind,iun)
        if(init(ind)) then
          do k=0,nbins
            ysum(ind,k) = 0D0
            ymean(ind,k) = 0D0
            hits(ind,k) = 0
            ysqsum(ind,k) = 0D0
            yerr(ind,k) = 0D0
          enddo
          init(ind) = .false.
        endif
        
        do k=0,nbins
          x(ind,k) = xm + k*(xp-xm)/nbins
        enddo
        
        do k=0,(nbins-1)
          if(xval.ge.x(ind,k) .and. xval.lt.x(ind,k+1)) then
            hits(ind,k) = hits(ind,k) + 1
            ysum(ind,k) = ysum(ind,k) + weight
            ysqsum(ind,k) = ysqsum(ind,k) + weight**2D0
          endif
        enddo
        
        ncalls(ind) = ncalls(ind) + 1
        if(mod(ncalls(ind),1000).eq.0) then
          do k=0,(nbins-1)
            if(hits(ind,k).gt.0) then
              ymean(ind,k) = ysum(ind,k)/hits(ind,k)
              yerr(ind,k) = ysqsum(ind,k)-ymean(ind,k)**2D0
              yerr(ind,k) = dsqrt(yerr(ind,k))/hits(ind,k)
            endif
          enddo
          
          open(iun, file=trim(adjustl(name))//'.top',status='replace')
          write(iun,*) " # ",trim(name),"  index ",ind
          do k=0,(nbins-1)
            write(iun,"(4(ES14.7))") x(ind,k),x(ind,k+1),ymean(ind,k),yerr(ind,k)
          enddo
          close(iun)
        endif  
      end
      
      subroutine index(name,ind,iun)
        implicit none
#include "hist.h"
        character*(*) name
        integer ind,iun
        integer i
        do i=1,nhistsmax
          if(trim(histlist(i)).eq.trim(name)) then
            ind = i
            iun = unitlist(i)
            return
          endif  
        enddo
        nhists = nhists + 1
        histlist(nhists) = name
        !call newunit(iun)
        iun = ind + 50
        unitlist(nhists) = iun
        ind = nhists
      end
c#######################################################################
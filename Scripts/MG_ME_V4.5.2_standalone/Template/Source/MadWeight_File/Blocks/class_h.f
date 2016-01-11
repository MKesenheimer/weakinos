      subroutine class_h
c***************************************************************************
c     ECS in CLASS H
c  
c
c***************************************************************************
      implicit none
      include '../../genps.inc'
      include '../../nexternal.inc'
c
c     arguments
c      
      double precision x(20)
      integer p1,p2,n_var
c
c     local
c
      double precision normp1,normp2,jac_temp,det,sqrts
      double precision angles(2,2),px(2),py(2),jac_loc
      integer j
      double precision  pboost(0:3), measure1, measure2
c
c     global
c
      double precision momenta(0:3,-max_branch:2*nexternal)  ! records the momenta of external/intermediate legs     (MG order)
      double precision mvir2(-max_branch:2*nexternal)         ! records the sq invariant masses of intermediate particles (MG order)
      common /to_diagram_kin/ momenta, mvir2
      double precision pmass(nexternal)     ! records the pole mass of any particle of the diagram  (MG order)
      common / to_mass/pmass
      double precision Etot,pztot,misspx,misspy
      common /to_missingP/Etot,pztot,misspx,misspy
      double precision              S,X1,X2,PSWGT,JAC
      common /PHASESPACE/ S,X1,X2,PSWGT,JAC
      double precision c_point(1:nexternal,3,2)
      common/ph_sp_init/c_point
c
c     external
c
      double precision phi
      external phi
c---
c Begin code
c---
      jac_loc=1d0
      sqrts=dsqrt(s)
      if((Etot+dsqrt(misspx**2+misspy**2)).gt.sqrts) then
      jac=-1d0
      return
      endif

c
c     HERE I NEED TO BOOST EVERYTHING TO THE FRAME WHERE PT=0
c       

c      write(*,*) "misspx",misspx
c      write(*,*) "misspy",misspy
      pboost(0)=Etot
      pboost(1)=misspx
      pboost(2)=misspy
      pboost(3)=0d0


      measure1=1d0
       do j=3,nexternal
         measure1=measure1*(momenta(1,j)**2+momenta(2,j)**2+momenta(3,j)**2)
     & *dsqrt(momenta(1,j)**2+momenta(2,j)**2)/momenta(0,j)**2
       enddo

       do j=3,nexternal
c         write(*,*) "p",j,momenta(0,j), momenta(1,j),momenta(2,j),momenta(3,j)
         call boostx(momenta(0,j),pboost,momenta(0,j))
       enddo

      measure2=1d0
       do j=3,nexternal
         measure2=measure2*(momenta(1,j)**2+momenta(2,j)**2+momenta(3,j)**2)
     & *dsqrt(momenta(1,j)**2+momenta(2,j)**2)/momenta(0,j)**2
       enddo

c
c     NEED TO RE-EVALUATE ETOT
c
       Etot=0d0
       misspx=0d0
       misspy=0d0
       do j=3,nexternal
         misspx=misspx+momenta(1,j)
         misspy=misspy+momenta(2,j)
         Etot=Etot+momenta(0,j)
c         write(*,*) "p",j,momenta(0,j), momenta(1,j),momenta(2,j),momenta(3,j)
       enddo
c         write(*,*) "misspx", misspx
c         write(*,*) "misspy", misspy

c      write(*,*) "measure1", measure1
c      write(*,*) "measure2", measure2
c      write(*,*) "measure1/measure2", measure1/measure2

c      pause
c
c    fill initial momenta
c
      x1=((Etot)
     & +(pztot))/sqrts
      x2=((Etot)-
     & (pztot))/sqrts

        if(dabs(x1-0.5d0).gt.0.5d0.or.dabs(x2-0.5d0).gt.0.5d0) then
        jac=-1d0
        momenta(0,p1)=-1
        momenta(0,p2)=-1
        return
        endif

      momenta(0,1)=sqrts*x1/2d0
      momenta(1,1)=0d0
      momenta(2,1)=0d0
      momenta(3,1)=sqrts*x1/2d0
      momenta(0,2)=sqrts*x2/2d0
      momenta(1,2)=0d0
      momenta(2,2)=0d0
      momenta(3,2)=-sqrts*x2/2d0
c
c     flux factor
c
      jac_loc=jac_loc/(S**2*x1*x2) ! flux + jac x1,x2 -> Etot, Pztot
      jac=jac*jac_loc

c
c     also need to rescale the weight to compensate for the transformation of the TF under boost:
c
      jac=jac*measure2/measure1
      return
      end

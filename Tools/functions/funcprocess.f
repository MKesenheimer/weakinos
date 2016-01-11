c############### funcprocess.f #########################################
c last modified by MK, 19.12.2014
c adapted from dislepton
c various shared helpers
c 2012-07 AvM

c############### decode_pair subroutine ###############################
c decode two PDG ids for MSSM fermions from a single integer

      subroutine decode_pair(combination,ida,idb)
        implicit none

        integer combination,ida,idb,idaa,idbb

        idaa = mod(abs(combination),1000)
        idbb = abs(combination)/1000
        ida =  (idaa/100)*1000000 + mod(idaa,100)
        if( combination .gt. 0D0) then
          idb =  (idbb/100)*1000000 + mod(idbb,100)
        else if (combination .lt. 0D0) then
          idb =  -(idbb/100)*1000000 - mod(idbb,100)
        else 
          print*,"error in decode_pair."
          stop
        endif  

      end

c############### end decode_pair subroutine ############################

c############### encode_pair subroutine ################################
c encode two PDG ids for MSSM fermions into a single integer
c (encodes 1000022,1000022 to 122122 or 1000024,-1000024 to -124124 etc.)

      integer function encode_pair(ida,idb)
        implicit none

        integer ida,idb,idaa,idbb,combination,idar,idbr

        idaa = (ida / 1000000)*100 + mod(ida,100)
        idbb = (idb / 1000000)*100 + mod(idb,100)
        
        combination = abs(idaa) + 1000*abs(idbb)
       
        ! encode unequal charge of particles
        if( (idaa .lt. 0D0) .or. (idbb .lt. 0D0) ) then
          combination = -combination
        endif

        ! check if reconstruction works, then we are fine in any case
        call decode_pair(combination,idar,idbr)
        if ((ida.ne.idar).or.(idb.ne.idbr)) then
          print*,"invalid particle ID"
          print*,"combination",combination
          print*,"ida,idaa,idar",ida,idaa,idar
          print*,"idb,idbb,idbr",idb,idbb,idbr
          stop  
        endif
        
        encode_pair = combination

      end

c############### end encode_pair subroutine ############################

c############### check_4conservation subroutine ########################
c check if 4 momentum conservation is fulfilled

      subroutine check_4conservation(p,nleg)
        implicit none
        
#include "nlegborn.h"

        integer nleg,i
        double precision p(0:3,nleg) ! momentum vectors
        double precision pi(0:3) ! sum of incoming momenta
        double precision pf(0:3) ! sum of outgoing momenta
        double precision eps ! rel. err.
        
        pi(:) = 0d0
        pf(:) = 0d0
        
        do i=1,2
          pi(:) = pi(:) + p(:,i)
        enddo
        
        do i=3,nleg
          pf(:) = pf(:) + p(:,i)
        enddo
        
        ! check up to double precision of incoming energy
        !eps = 1d-15*pi(0)
        eps = 1d-10*pi(0)

#ifdef DEBUGQ
        print*,"p1 = ", p(:,1)
        print*,"p2 = ", p(:,2)
        print*,"p3 = ", p(:,3)
        print*,"p4 = ", p(:,4)
        if(nleg .eq. 5) print*,"p5 = ", p(:,5)
        if(nleg .eq. 6) print*,"p6 = ", p(:,6)
        print*, "Sum p in  = ", pi(:)
        print*, "Sum p out = ", pf(:)
#endif

#define CHECKMOM
#ifdef CHECKMOM
        if( (dabs(pi(0) - pf(0)) .gt. dabs(eps)) .or.
     &      (dabs(pi(1) - pf(1)) .gt. dabs(eps)) .or.
     &      (dabs(pi(2) - pf(2)) .gt. dabs(eps)) .or.
     &      (dabs(pi(3) - pf(3)) .gt. dabs(eps))) then
          print*, "Error: four momentum not conserved."
          print*,"p1 = ", p(:,1)
          print*,"p2 = ", p(:,2)
          print*,"p3 = ", p(:,3)
          print*,"p4 = ", p(:,4)
          if(nleg .eq. 5) print*,"p5 = ", p(:,5)
          if(nleg .eq. 6) print*,"p6 = ", p(:,6)
          print*, "Sum p in  = ", pi(:)
          print*, "Sum p out = ", pf(:)
          stop
        endif
#endif
      end

c############### end check_4conservation subroutine ####################

c############### subroutine set_squark_params ##########################
c pick SUSY masses relevant for specific initial state
      subroutine set_process(id, M1, M2, M3, M4)
        implicit none

#include "PhysPars.h"

        integer id(4)
        double precision M1, M2, M3, M4
        
        M3 = par_Fin1mass
        M4 = par_Fin2mass
        
        select case(abs(id(1)))
        case(1) ! d
          M1 = par_MD
        case(2) ! u
          M1 = par_MU
        case(3) ! s
          M1 = par_MS
        case(4) ! c
          M1 = par_MC
        case(5) ! b
          M1 = par_MB
        case(6) ! t
          print*, "top quarks not implemented yet."
          stop
        case default
          write(*,*) "encountered unhandled incoming quark ID ", id
          stop
        end select
        
        select case(abs(id(2)))
        case(1) ! d
          M2 = par_MD
        case(2) ! u
          M2 = par_MU
        case(3) ! s
          M2 = par_MS
        case(4) ! c
          M2 = par_MC
        case(5) ! b
          M2 = par_MB
        case(6) ! t
          print*, "top quarks not implemented yet."
          stop
        case default
          write(*,*) "encountered unhandled incoming quark ID ", id
          stop
        end select
        
#ifdef DEBUGQ
        print*,"id",id
        print*,"M1",M1
        print*,"M2",M2
        print*,"M3",M3
        print*,"M4",M4
#endif
        
      end

c############### end subroutine set_squark_params ######################

c############### subroutine switchmom ##################################
c Changes stuff for crossings (copied from MadGraph V4.5.2)
      subroutine switchmom(p1,p,ic,jc,nexternal)
        implicit none
        integer nexternal
        integer jc(nexternal),ic(nexternal)
        real*8 p1(0:3,nexternal),p(0:3,nexternal)
        integer i,j

        do i=1,nexternal
          do j=0,3
             p(j,ic(i))=p1(j,i)
          enddo
        enddo
        do i=1,nexternal
          jc(i)=1
        enddo
        jc(ic(1))=-1
        jc(ic(2))=-1
      end
c############### end subroutine switchmom ##############################

c############### function kaellen ######################################
c calculates the kaellen function and the sqrt of kaellen function
      double precision function kaellen(x, y, z)

        implicit none
        double precision x, y, z
        
        kaellen = x**2+y**2+z**2-2*(x*y+x*z+y*z)
        
      end function kaellen

      double precision function kaellenSqrt(x, y, z)

        implicit none
        double precision x, y, z
        
        kaellenSqrt = dsqrt(x**2+y**2+z**2-2*(x*y+x*z+y*z))
        
      end function kaellenSqrt
c############### end function kaellen ##################################


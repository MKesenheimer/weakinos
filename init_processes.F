c############### init_processes.f ######################################
c last modified by MK, 05.12.2015
c weakino pair production

c############### subroutine init_processes #############################
c     setup Born and real subprocesses
c     set process id using slepton types
      subroutine init_processes
        implicit none
        
#include "nlegborn.h"
#include "pwhg_flst.h"
#include "LesHouches.h"
#include "Flags.h"
#include "pwhg_st.h"
#include "pwhg_flst_add.h"
#include "osres.h"

        integer i1,i2,i3,i4,i5,k,ii(5)
        equivalence (i1,ii(1)),(i2,ii(2)),(i3,ii(3)),(i4,ii(4)),
     &              (i5,ii(5))
        logical debug
        parameter (debug=.true.)
        integer j
        real*8 powheginput
        external powheginput
        logical condition, mixed_channels
        integer max_flav,encode_pair

        ! set number of active flavors for incoming partons
        max_flav = int(powheginput('incomingflavors'))
        if (max_flav.eq.5)  print*,"WARNING: ignoring sbottom mixing"
        if (max_flav.gt.5)  stop "can't handle more than 5 flavours"

        ! whether to use fake virtuals or not
        if(powheginput("#fakevirtuals").gt.0) then
           if(int(powheginput("fakevirtuals")).gt.0) then
             flg_fakevirtuals = .true.
           endif
        else
          flg_fakevirtuals = .false.
        endif
        
        ! number of light quarks:
#ifdef NEGLECTBMASS
        st_nlight = 5
#else
        st_nlight = 4
#endif
        
        ! choose process
        i3 = int(powheginput("fin1"))
        i4 = int(powheginput("fin2"))

        ! channels with two neutralinos or two charginos in final state
        if (      ((i3.eq.1000022).and.(i4.eq.1000022))
     &       .or. ((i3.eq.1000022).and.(i4.eq.1000023))
     &       .or. ((i3.eq.1000022).and.(i4.eq.1000025))
     &       .or. ((i3.eq.1000022).and.(i4.eq.1000035))
     
     &       .or. ((i3.eq.1000023).and.(i4.eq.1000023))
     &       .or. ((i3.eq.1000023).and.(i4.eq.1000025))
     &       .or. ((i3.eq.1000023).and.(i4.eq.1000035))
     
     &       .or. ((i3.eq.1000025).and.(i4.eq.1000025))
     &       .or. ((i3.eq.1000025).and.(i4.eq.1000035))
     
     &       .or. ((i3.eq.1000035).and.(i4.eq.1000035))
     
     &       .or. ((i3.eq.1000024).and.(i4.eq.-1000024))
     &       .or. ((i3.eq.1000024).and.(i4.eq.-1000037))
     &       .or. ((i3.eq.1000037).and.(i4.eq.-1000037)) ) then
          mixed_channels = .false.

        ! channels with one neutralino and one chargino in final state
        else if ( ((i3.eq.1000022).and.(i4.eq.-1000024))
     &       .or. ((i3.eq.1000022).and.(i4.eq.1000024))
     &       .or. ((i3.eq.1000022).and.(i4.eq.-1000037))
     &       .or. ((i3.eq.1000022).and.(i4.eq.1000037))

     &       .or. ((i3.eq.1000023).and.(i4.eq.-1000024))
     &       .or. ((i3.eq.1000023).and.(i4.eq.1000024))
     &       .or. ((i3.eq.1000023).and.(i4.eq.-1000037))
     &       .or. ((i3.eq.1000023).and.(i4.eq.1000037))
     
     &       .or. ((i3.eq.1000025).and.(i4.eq.-1000024))
     &       .or. ((i3.eq.1000025).and.(i4.eq.1000024))
     &       .or. ((i3.eq.1000025).and.(i4.eq.-1000037))
     &       .or. ((i3.eq.1000025).and.(i4.eq.1000037))
     
     &       .or. ((i3.eq.1000035).and.(i4.eq.-1000024))
     &       .or. ((i3.eq.1000035).and.(i4.eq.1000024))
     &       .or. ((i3.eq.1000035).and.(i4.eq.-1000037))
     &       .or. ((i3.eq.1000035).and.(i4.eq.1000037)) ) then
          mixed_channels = .true.

        else
          print*, "requested final states not implemented"
          print*, "swap final states and try again"
          stop
        endif

        print*
        print*, "===================================================="
        print*, "        neutralino/chargino pair production         "
        print*, "      final states: ", i3, i4
        print*, "                                                    "
        if     (i3.eq.1000022) then
          print*, "             q                         n1           "
        elseif (i3.eq.1000023) then
          print*, "             q                         n2           "
        elseif (i3.eq.1000025) then
          print*, "             q                         n3           "
        elseif (i3.eq.1000035) then
          print*, "             q                         n4           "
        elseif (i3.eq.-1000024) then
          print*, "             q                         x1-          "
        elseif (i3.eq.1000024) then
          print*, "             q                         x1+          "
        elseif (i3.eq.-1000037) then
          print*, "             q                         x2-          "
        elseif (i3.eq.1000037) then
          print*, "             q                         x2+          "
        endif
        print*, "              \                       /             "
        print*, "               \                     /              "
        print*, "                \                   /               "
        print*, "                 \       Z/W       /                "
        print*, "                  ~~~~~~~~~~~~~~~~~                 "
        print*, "                 /                 \                "
        print*, "                /                   \               "
        print*, "               /                     \              "
        print*, "              /                       \             "
        if     (i4.eq.1000022) then
          print*, "             q*                        n1           "
        elseif (i4.eq.1000023) then
          print*, "             q*                        n2           "
        elseif (i4.eq.1000025) then
          print*, "             q*                        n3           "
        elseif (i4.eq.1000035) then
          print*, "             q*                        n4           "
        elseif (i4.eq.-1000024) then
          print*, "             q*                        x1-          "
        elseif (i4.eq.1000024) then
          print*, "             q*                        x1+          "
        elseif (i4.eq.-1000037) then
          print*, "             q*                        x2-          "
        elseif (i4.eq.1000037) then
          print*, "             q*                        x2+          "
        endif 
        print*, "                                                    "
        print*, "===================================================="
        print* 
        print*,' number of incoming flavors: ', max_flav

        ! set process id: encode types of slepton A and B
        lprup(1) = encode_pair(i3,i4)
        print*,' encoded final states to process ID: ', lprup(1)

        ! index of the first coloured particle in the final state
        ! (all subsequent particles are coloured)
        flst_lightpart=5

        ! Born subprocesses
        flst_nborn=0
        do i1=-max_flav,max_flav
          do i2=-max_flav,max_flav
            ! q qbar (with mixed flavor in the initial state, e.g. u dbar)
            ! d ubar channels
            if(mixed_channels.and.i4.lt.0) then
              if( (i1.ne.0) .and. (i2.ne.0) .and. (i1+i2.eq.-1)
     &           .and. (mod(i1,2).eq.0 .or. mod(i1,2).eq.1)
     &           .and. (mod(i2,2).eq.0 .or. mod(i2,2).eq.1) ) then
                flst_nborn=flst_nborn+1
                if(flst_nborn.gt.maxprocborn) goto 999
                do k=1,nlegborn
                  flst_born(k,flst_nborn)=ii(k)
                enddo
              endif
            ! dbar u channels
            elseif(mixed_channels .and. i4.gt.0) then
              if( (i1.ne.0) .and. (i2.ne.0) .and. (i1+i2.eq.1)
     &           .and. (mod(i1,2).eq.-1 .or. mod(i1,2).eq.0)
     &           .and. (mod(i2,2).eq.-1 .or. mod(i2,2).eq.0) ) then
                flst_nborn=flst_nborn+1
                if(flst_nborn.gt.maxprocborn) goto 999
                do k=1,nlegborn
                  flst_born(k,flst_nborn)=ii(k)
                enddo
              endif
            ! q qbar (same flavors in the initial state)
            else
              if(i1.ne.0.and.i1+i2.eq.0) then
                flst_nborn=flst_nborn+1
                if(flst_nborn.gt.maxprocborn) goto 999
                do k=1,nlegborn
                  flst_born(k,flst_nborn)=ii(k)
                enddo
              endif
            endif
          enddo
        enddo

        ! real subprocesses
        flst_nreal=0
        do i1=-max_flav,max_flav
          do i2=-max_flav,max_flav
            do i5=-max_flav,max_flav
              condition=.false.
              ! q qbar (with mixed flavor in the initial state, e.g. u dbar)
              ! d ubar channels
              if(mixed_channels.and.i4.lt.0) then
                if(.not.(i1.eq.0.and.i2.eq.0)
     &           .and. (mod(i1,2).eq.0 .or. mod(i1,2).eq.1)
     &           .and. (mod(i2,2).eq.0 .or. mod(i2,2).eq.1) ) then   ! exclude gg
                  if((i1.ne.0).and.(i1+i2.eq.-1).and.(i5.eq.0)) then ! q qbar -> g
                    condition=.true.
                  elseif((i1.eq.0).and.(i2.eq.i5-1)
     &               .and. (mod(i2,2).eq.0 .or. mod(i2,2).eq.1)
     &               .and. (mod(i5,2).eq.0 .or. mod(i5,2).eq.-1) ) then ! g q
                    condition=.true.
                  elseif((i2.eq.0).and.(i1.eq.i5-1)
     &               .and. (mod(i1,2).eq.0 .or. mod(i1,2).eq.1)
     &               .and. (mod(i5,2).eq.0 .or. mod(i5,2).eq.-1) ) then ! q g
                    condition=.true.
                  endif
                endif
              ! dbar u channels
              elseif(mixed_channels.and.i4.gt.0) then
                if(.not.(i1.eq.0.and.i2.eq.0)
     &           .and. (mod(i1,2).eq.-1 .or. mod(i1,2).eq.0)
     &           .and. (mod(i2,2).eq.-1 .or. mod(i2,2).eq.0) ) then ! exclude gg
                  if((i1.ne.0).and.(i1+i2.eq.1).and.(i5.eq.0)) then ! q qbar -> g
                    condition=.true.
                  elseif((i1.eq.0).and.(i2.eq.i5+1)
     &               .and. (mod(i2,2).eq.0 .or. mod(i2,2).eq.-1)
     &               .and. (mod(i5,2).eq.0 .or. mod(i5,2).eq.1) ) then ! g q
                    condition=.true.
                  elseif((i2.eq.0).and.(i1.eq.i5+1)
     &               .and. (mod(i1,2).eq.0 .or. mod(i1,2).eq.-1)
     &               .and. (mod(i5,2).eq.0 .or. mod(i5,2).eq.1) ) then ! q g
                    condition=.true.
                  endif
                endif
              ! q qbar (same flavors in the initial state)
              else
                if(.not.(i1.eq.0.and.i2.eq.0)) then                 ! exclude gg
                  if((i1.ne.0).and.(i1+i2.eq.0).and.(i5.eq.0)) then ! q qbar -> g
                    condition=.true.
                  elseif((i1.eq.0).and.(i2.eq.i5)) then ! g q
                    condition=.true.
                  elseif((i2.eq.0).and.(i1.eq.i5)) then ! q g
                    condition=.true.
                  endif
                endif
              endif
              if(condition) then
                flst_nreal=flst_nreal+1
                if(flst_nreal.gt.maxprocreal) goto 998
                do k=1,nlegreal
                  flst_real(k,flst_nreal)=ii(k)
                enddo
              endif
            enddo
          enddo
        enddo


#ifdef DEBUG
        flst_nborn = 1
        flst_nreal = 1
        flst_born(1,1) = -4
        flst_born(2,1) =  3
        flst_born(3,1) = 1000022
        flst_born(4,1) = -1000024
        flst_real(1,1) = -4
        flst_real(2,1) =  0
        flst_real(3,1) = 1000022
        flst_real(4,1) = -1000024
        flst_real(5,1) = -3
#endif
        
#ifdef DEBUG
           print*,' born processes',flst_nborn
           do j=1,flst_nborn
              print*,(flst_born(k,j),k=1,nlegborn)
           enddo
#endif

#ifdef DEBUG
        print*,' real processes',flst_nreal
        do j=1,flst_nreal
          print*,(flst_real(k,j),k=1,nlegreal)
        enddo
#endif

#if defined(DSUB_I) || defined(DSUB_II) || defined(DSUB_II_TEST)
        ! setup the flavor list of on-shell resonant diagrams
        call init_processes_osres
#ifdef DEBUG
        print*,' on-shell resonant real processes',flst_nosres
        do j=1,flst_nosres
          print*,(flst_osres(k,j),k=1,nlegreal)
        enddo
#endif
#endif

        return
 998    print*,'init_processes: increase maxprocreal'
        stop
 999    print*,'init_processes: increase maxprocborn'
        stop
      end
      
c############### end subroutine init_processes #########################


c############### subroutine init_osres #################################
c if the real-amplitudes are resonant split into a regular part and the 
c IR-divergent parts: this routine provides the corresponding flavor 
c lists and stores it in flst_osres
c this list should not be confused with flst_realres.
c osresID 1 = dl35: left handed up-type squark gets resonant with particles 3 & 5
c osresID 2 = dl45: left handed up-type squark gets resonant with particles 4 & 5
c osresID 3 = ul35: left handed down-type squark gets resonant with particles 3 & 5
c osresID 4 = ul45: left handed down-type squark gets resonant with particles 4 & 5
c osresID 5 = dr35: right handed up-type squark gets resonant with particles 3 & 5
c osresID 6 = dr45: right handed up-type squark gets resonant with particles 4 & 5
c osresID 7 = ur35: right handed down-type squark gets resonant with particles 3 & 5
c osresID 8 = ur45: right handed down-type squark gets resonant with particles 4 & 5
      subroutine init_processes_osres
        implicit none
        
#include "nlegborn.h"
#include "pwhg_flst.h"
#include "pwhg_flg.h"
#include "pwhg_flst_add.h"
#include "osres.h"

        ! local variables
        integer i1,i2
        integer k
        
        ! reset the number of on-shell resonant diagrams
        flst_nosres = 0
        ! fill the flavor list with on-shell resonant diagrams
        do i1=1,flst_nreal
          if((flst_real(1,i1).eq.0.or.flst_real(2,i1).eq.0)
     &       .and.flst_realres(nlegreal,i1).eq.0) then ! only for production-channels
            flst_nosres=flst_nosres+1 
            do i2=1,nlegreal
              flst_osres(i2,flst_nosres) = flst_real(i2,i1) ! copy existing real
            enddo
          endif
        enddo        

        ! fill the global variables for treating the on-shell resonances
        ! number of resonances
        nosres = 4 ! change cnosres in osres.h as well

        if(nosres.ne.cnosres) then
          print*,"error in init_processes.F: nosres and narraylngth not"
          print*,"equal: ", nosres, cnosres
          print*,"change nosres in init_processes or"//
     &           "cnosres in osres.h"
          stop
        endif

        ! the human readable identifiers
        osresID(1) = "dl35"
        osresID(2) = "dl45"
        osresID(3) = "ul35"
        osresID(4) = "ul45"
        !if(nosres.eq.8) then
        !  osresID(5) = "dr35"
        !  osresID(6) = "dr45"
        !  osresID(7) = "ur35"
        !  osresID(8) = "ur45"
        !endif

        do k=1,nosres
          if(len(osresID(k)).gt.4) then
            print*,"error in init_processes.F: array osresID not set"
            print*,"correctly. Increase field length."
            stop
          endif
        enddo

        ! save which final state particles get resonant in an 
        ! additional array. This will be helpful to generate the 
        ! on-shell real phase space and to calculate the matrix elements
        ! in subroutine setosresreal
        osreslegs(1,1) = 3
        osreslegs(1,2) = 5
        osreslegs(1,3) = 4
        
        osreslegs(2,1) = 4
        osreslegs(2,2) = 5
        osreslegs(2,3) = 3
        
        osreslegs(3,1) = 3
        osreslegs(3,2) = 5
        osreslegs(3,3) = 4
        
        osreslegs(4,1) = 4
        osreslegs(4,2) = 5
        osreslegs(4,3) = 3

        ! TODO: get rid of the compiler warnings
        !if(nosres.eq.8) then
        !  osreslegs(5,1) = 3
        !  osreslegs(5,2) = 5
        !  osreslegs(5,3) = 4
        !
        !  osreslegs(6,1) = 4
        !  osreslegs(6,2) = 5
        !  osreslegs(6,3) = 3
        !
        !  osreslegs(7,1) = 3
        !  osreslegs(7,2) = 5
        !  osreslegs(7,3) = 4
        !
        !  osreslegs(8,1) = 4
        !  osreslegs(8,2) = 5
        !  osreslegs(8,3) = 3
        !endif

        ! the rest of the variables defined in osres.h are set
        ! in init_couplings.f

#ifdef DEBUGQ
        osresID(1) = "dl35"
        osresID(2) = "dl35"
        osresID(3) = "dl35"
        osresID(4) = "dl35"
        
        osreslegs(1,1) = 3
        osreslegs(1,2) = 5
        osreslegs(1,3) = 4
        
        osreslegs(2,1) = 3
        osreslegs(2,2) = 5
        osreslegs(2,3) = 4
        
        osreslegs(3,1) = 3
        osreslegs(3,2) = 5
        osreslegs(3,3) = 4
        
        osreslegs(4,1) = 3
        osreslegs(4,2) = 5
        osreslegs(4,3) = 4
#endif
        
      end

c############### end subroutine init_osres #############################
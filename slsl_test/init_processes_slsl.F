c     setup Born and real subprocesses
c     set process id using slepton types
      subroutine init_processes
        implicit none
        
#include "nlegborn.h"
#include "pwhg_flst.h"
#include "LesHouches.h"
#include "Flags.h"
#include "pwhg_st.h"

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
        
        ! choose process ------------------------------------------------

        i3 = int(powheginput("fin1"))
        i4 = int(powheginput("fin2"))

        ! channels with two neutralinos or two charginos in final state
        if (      ((i3.eq.1000022).and.(i4.eq.1000022))
     &       .or. ((i3.eq.1000022).and.(i4.eq.1000023))
     &       .or. ((i3.eq.1000023).and.(i4.eq.1000023))
     &       .or. ((i3.eq.1000024).and.(i4.eq.-1000024))
     &       .or. ((i3.eq.1000024).and.(i4.eq.-1000037))
     &       .or. ((i3.eq.1000037).and.(i4.eq.-1000037))
     &       .or. ((i3.eq.2000011).and.(i4.eq.-2000011)) ) then
          mixed_channels = .false.

        ! channels with one neutralino and one chargino in final state
        else if ( ((i3.eq.1000022).and.(i4.eq.-1000024))
     &       .or. ((i3.eq.1000022).and.(i4.eq.1000024))
     &       .or. ((i3.eq.1000023).and.(i4.eq.-1000037))
     &       .or. ((i3.eq.1000023).and.(i4.eq.1000037))
     &       .or. ((i3.eq.1000023).and.(i4.eq.-1000024))
     &       .or. ((i3.eq.1000023).and.(i4.eq.1000024))
     &       .or. ((i3.eq.1000022).and.(i4.eq.-1000037))
     &       .or. ((i3.eq.1000022).and.(i4.eq.1000037)) ) then
          mixed_channels = .true.

        else
          print*, "requested final states not yet implemented"
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
        elseif (i3.eq.1000024) then
          print*, "             q                         x1+          "
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

#ifdef DEBUG
           print*,' born processes',flst_nborn
           do j=1,flst_nborn
              print*,(flst_born(k,j),k=1,nlegborn)
           enddo
#endif

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
     &           .and. (mod(i2,2).eq.0 .or. mod(i2,2).eq.1) ) then ! exclude gg
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
                if(.not.(i1.eq.0.and.i2.eq.0)) then ! exclude gg
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
        print*,' real processes',flst_nreal
        do j=1,flst_nreal
          print*,(flst_real(k,j),k=1,nlegreal)
        enddo
#endif

        return
 998    print*,'init_processes: increase maxprocreal'
        stop
 999    print*,'init_processes: increase maxprocborn'
        stop
      end

      subroutine sborn_proc(p_born,legs,wgt,wgtjk,wgtmunu)
      implicit none
      include "nexternal.inc"
      include "coupl.inc"
      double precision wgt
      double precision p_born(0:3,nexternal-1),wgt2(nexternal-1),
     &   wgtjk(nexternal-1,nexternal-1)
      integer legs(nexternal-1),lstr,i
      character*140 str
      logical calculatedBorn
      integer skip
      common/cBorn/calculatedBorn,skip
      double precision wgtmunu(0:3,0:3,nexternal-1)
      
 
      calculatedBorn=.false.
      
      call convert_to_string(nexternal-1,legs,str,lstr)
      
      if(str.eq."1-21000022-1000024") then ! d ubar -> ni + cha1-
         call sborn_cl_nIxJm_M_001(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_001(p_born,wgtjk)
         goto 20
      elseif(str.eq."-211000022-1000024") then ! ubar d -> ni + cha1-
         call sborn_cl_nIxJm_M_002(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_002(p_born,wgtjk)
         goto 20
      elseif(str.eq."-431000022-1000024") then ! cbar s -> ni + cha1-
         call sborn_cl_nIxJm_M_003(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_003(p_born,wgtjk)
         goto 20
      elseif(str.eq."3-41000022-1000024") then ! s cbar -> ni + cha1-
         call sborn_cl_nIxJm_M_004(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_004(p_born,wgtjk)
         goto 20
      elseif(str.eq."-1210000221000024") then
         call sborn_cl_nIxJ_p_001(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_001(p_born,wgtjk)
         goto 20
      elseif(str.eq."2-110000221000024") then
         call sborn_cl_nIxJ_p_002(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_002(p_born,wgtjk)
         goto 20
      elseif(str.eq."4-310000221000024") then
         call sborn_cl_nIxJ_p_003(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_003(p_born,wgtjk)
         goto 20
      elseif(str.eq."-3410000221000024") then
         call sborn_cl_nIxJ_p_004(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_004(p_born,wgtjk)
         goto 20
      
      elseif(str.eq."1-21000023-1000024") then ! d ubar -> ni + cha1-
         call sborn_cl_nIxJm_M_001(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_001(p_born,wgtjk)
         goto 20
      elseif(str.eq."-211000023-1000024") then ! ubar d -> ni + cha1-
         call sborn_cl_nIxJm_M_002(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_002(p_born,wgtjk)
         goto 20
      elseif(str.eq."-431000023-1000024") then ! cbar s -> ni + cha1-
         call sborn_cl_nIxJm_M_003(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_003(p_born,wgtjk)
         goto 20
      elseif(str.eq."3-41000023-1000024") then ! s cbar -> ni + cha1-
         call sborn_cl_nIxJm_M_004(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_004(p_born,wgtjk)
         goto 20
      elseif(str.eq."-1210000231000024") then
         call sborn_cl_nIxJ_p_001(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_001(p_born,wgtjk)
         goto 20
      elseif(str.eq."2-110000231000024") then
         call sborn_cl_nIxJ_p_002(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_002(p_born,wgtjk)
         goto 20
      elseif(str.eq."4-310000231000024") then
         call sborn_cl_nIxJ_p_003(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_003(p_born,wgtjk)
         goto 20
      elseif(str.eq."-3410000231000024") then
         call sborn_cl_nIxJ_p_004(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_004(p_born,wgtjk)
         goto 20
      
      elseif(str.eq."1-21000022-1000037") then ! d ubar -> ni + cha1-
         call sborn_cl_nIxJm_M_001(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_001(p_born,wgtjk)
         goto 20
      elseif(str.eq."-211000022-1000037") then ! ubar d -> ni + cha1-
         call sborn_cl_nIxJm_M_002(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_002(p_born,wgtjk)
         goto 20
      elseif(str.eq."-431000022-1000037") then ! cbar s -> ni + cha1-
         call sborn_cl_nIxJm_M_003(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_003(p_born,wgtjk)
         goto 20
      elseif(str.eq."3-41000022-1000037") then ! s cbar -> ni + cha1-
         call sborn_cl_nIxJm_M_004(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_004(p_born,wgtjk)
         goto 20
      elseif(str.eq."-1210000221000037") then
         call sborn_cl_nIxJ_p_001(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_001(p_born,wgtjk)
         goto 20
      elseif(str.eq."2-110000221000037") then
         call sborn_cl_nIxJ_p_002(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_002(p_born,wgtjk)
         goto 20
      elseif(str.eq."4-310000221000037") then
         call sborn_cl_nIxJ_p_003(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_003(p_born,wgtjk)
         goto 20
      elseif(str.eq."-3410000221000037") then
         call sborn_cl_nIxJ_p_004(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_004(p_born,wgtjk)
         goto 20
      
      elseif(str.eq."1-21000023-1000037") then ! d ubar -> ni + cha1-
         call sborn_cl_nIxJm_M_001(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_001(p_born,wgtjk)
         goto 20
      elseif(str.eq."-211000023-1000037") then ! ubar d -> ni + cha1-
         call sborn_cl_nIxJm_M_002(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_002(p_born,wgtjk)
         goto 20
      elseif(str.eq."-431000023-1000037") then ! cbar s -> ni + cha1-
         call sborn_cl_nIxJm_M_003(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_003(p_born,wgtjk)
         goto 20
      elseif(str.eq."3-41000023-1000037") then ! s cbar -> ni + cha1-
         call sborn_cl_nIxJm_M_004(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_004(p_born,wgtjk)
         goto 20
      elseif(str.eq."-1210000231000037") then
         call sborn_cl_nIxJ_p_001(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_001(p_born,wgtjk)
         goto 20
      elseif(str.eq."2-110000231000037") then
         call sborn_cl_nIxJ_p_002(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_002(p_born,wgtjk)
         goto 20
      elseif(str.eq."4-310000231000037") then
         call sborn_cl_nIxJ_p_003(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_003(p_born,wgtjk)
         goto 20
      elseif(str.eq."-3410000231000037") then
         call sborn_cl_nIxJ_p_004(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_004(p_born,wgtjk)
         goto 20
      
      elseif(str.eq."1-21000025-1000024") then ! d ubar -> ni + cha1-
         call sborn_cl_nIxJm_M_001(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_001(p_born,wgtjk)
         goto 20
      elseif(str.eq."-211000025-1000024") then ! ubar d -> ni + cha1-
         call sborn_cl_nIxJm_M_002(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_002(p_born,wgtjk)
         goto 20
      elseif(str.eq."-431000025-1000024") then ! cbar s -> ni + cha1-
         call sborn_cl_nIxJm_M_003(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_003(p_born,wgtjk)
         goto 20
      elseif(str.eq."3-41000025-1000024") then ! s cbar -> ni + cha1-
         call sborn_cl_nIxJm_M_004(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_004(p_born,wgtjk)
         goto 20
      elseif(str.eq."-1210000251000024") then
         call sborn_cl_nIxJ_p_001(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_001(p_born,wgtjk)
         goto 20
      elseif(str.eq."2-110000251000024") then
         call sborn_cl_nIxJ_p_002(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_002(p_born,wgtjk)
         goto 20
      elseif(str.eq."4-310000251000024") then
         call sborn_cl_nIxJ_p_003(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_003(p_born,wgtjk)
         goto 20
      elseif(str.eq."-3410000251000024") then
         call sborn_cl_nIxJ_p_004(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_004(p_born,wgtjk)
         goto 20
     
      elseif(str.eq."1-21000025-1000037") then ! d ubar -> ni + cha1-
         call sborn_cl_nIxJm_M_001(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_001(p_born,wgtjk)
         goto 20
      elseif(str.eq."-211000025-1000037") then ! ubar d -> ni + cha1-
         call sborn_cl_nIxJm_M_002(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_002(p_born,wgtjk)
         goto 20
      elseif(str.eq."-431000025-1000037") then ! cbar s -> ni + cha1-
         call sborn_cl_nIxJm_M_003(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_003(p_born,wgtjk)
         goto 20
      elseif(str.eq."3-41000025-1000037") then ! s cbar -> ni + cha1-
         call sborn_cl_nIxJm_M_004(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_004(p_born,wgtjk)
         goto 20
      elseif(str.eq."-1210000251000037") then
         call sborn_cl_nIxJ_p_001(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_001(p_born,wgtjk)
         goto 20
      elseif(str.eq."2-110000251000037") then
         call sborn_cl_nIxJ_p_002(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_002(p_born,wgtjk)
         goto 20
      elseif(str.eq."4-310000251000037") then
         call sborn_cl_nIxJ_p_003(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_003(p_born,wgtjk)
         goto 20
      elseif(str.eq."-3410000251000037") then
         call sborn_cl_nIxJ_p_004(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_004(p_born,wgtjk)
         goto 20
     
      elseif(str.eq."1-21000035-1000024") then ! d ubar -> ni + cha1-
         call sborn_cl_nIxJm_M_001(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_001(p_born,wgtjk)
         goto 20
      elseif(str.eq."-211000035-1000024") then ! ubar d -> ni + cha1-
         call sborn_cl_nIxJm_M_002(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_002(p_born,wgtjk)
         goto 20
      elseif(str.eq."-431000035-1000024") then ! cbar s -> ni + cha1-
         call sborn_cl_nIxJm_M_003(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_003(p_born,wgtjk)
         goto 20
      elseif(str.eq."3-41000035-1000024") then ! s cbar -> ni + cha1-
         call sborn_cl_nIxJm_M_004(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_004(p_born,wgtjk)
         goto 20
      elseif(str.eq."-1210000351000024") then
         call sborn_cl_nIxJ_p_001(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_001(p_born,wgtjk)
         goto 20
      elseif(str.eq."2-110000351000024") then
         call sborn_cl_nIxJ_p_002(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_002(p_born,wgtjk)
         goto 20
      elseif(str.eq."4-310000351000024") then
         call sborn_cl_nIxJ_p_003(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_003(p_born,wgtjk)
         goto 20
      elseif(str.eq."-3410000351000024") then
         call sborn_cl_nIxJ_p_004(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_004(p_born,wgtjk)
         goto 20
     
      elseif(str.eq."1-21000035-1000037") then ! d ubar -> ni + cha1-
         call sborn_cl_nIxJm_M_001(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_001(p_born,wgtjk)
         goto 20
      elseif(str.eq."-211000035-1000037") then ! ubar d -> ni + cha1-
         call sborn_cl_nIxJm_M_002(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_002(p_born,wgtjk)
         goto 20
      elseif(str.eq."-431000035-1000037") then ! cbar s -> ni + cha1-
         call sborn_cl_nIxJm_M_003(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_003(p_born,wgtjk)
         goto 20
      elseif(str.eq."3-41000035-1000037") then ! s cbar -> ni + cha1-
         call sborn_cl_nIxJm_M_004(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJm_M_004(p_born,wgtjk)
         goto 20
      elseif(str.eq."-1210000351000037") then
         call sborn_cl_nIxJ_p_001(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_001(p_born,wgtjk)
         goto 20
      elseif(str.eq."2-110000351000037") then
         call sborn_cl_nIxJ_p_002(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_002(p_born,wgtjk)
         goto 20
      elseif(str.eq."4-310000351000037") then
         call sborn_cl_nIxJ_p_003(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_003(p_born,wgtjk)
         goto 20
      elseif(str.eq."-3410000351000037") then
         call sborn_cl_nIxJ_p_004(p_born,wgtmunu,wgt2)
         call sborn_sf_nIxJ_p_004(p_born,wgtjk)
         goto 20
      endif
      
 20   wgt=0d0
      do i=1,nexternal-1
         if(wgt.eq.0d0 .and. wgt2(i).ne.0d0) then
            wgt=wgt2(i)
         elseif(wgt.ne.0d0 .and. wgt2(i).ne.0d0 .and.
     &           abs((wgt-wgt2(i))/wgt).gt.1d-7) then
            write(*,*) "Error #2 in sborn_proc ",i,wgt2
            stop
         endif
      enddo
      
      end
      
      
      
      
      
      subroutine born_color(legs,color)
      implicit none
      include "nexternal.inc"
      integer maxamps
      parameter(maxamps=6000)
      Double Precision amp2001(maxamps), jamp2001(0:maxamps)
      common/to_amps_001/amp2001,jamp2001
      Double Precision amp2002(maxamps), jamp2002(0:maxamps)
      common/to_amps_002/amp2002,jamp2002
      Double Precision amp2003(maxamps), jamp2003(0:maxamps)
      common/to_amps_003/amp2003,jamp2003
      Double Precision amp2004(maxamps), jamp2004(0:maxamps)
      common/to_amps_004/amp2004,jamp2004
      double precision jamp2cum(0:maxamps)
      integer ICOLUP(2,nexternal-1,maxamps)
      integer color(2,nexternal-1),color1(2,nexternal-1)
      double precision random,xtarget
      external random
      integer legs(nexternal-1),lstr,i,j
      character*140 str
      integer iflow,ifl

#ifdef DEBUGQ
      print*, "[DEBUG] str   = ", str
      print*, "[DEBUG] legs  = ", legs
      print*, "[DEBUG] color = ", color
      print*, "[DEBUG] amp2001 = ",amp2001(1)
      print*, "[DEBUG] jamp2001 = ",jamp2001(0)
      stop
#endif

      call convert_to_string(nexternal-1,legs,str,lstr)
      
      if(str.eq."1-21000022-1000024") then
         include "leshouches_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-211000022-1000024") then
         include "leshouches_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."-431000022-1000024") then
         include "leshouches_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."3-41000022-1000024") then
         include "leshouches_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-1210000221000024") then
         include "leshouches_005.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."2-110000221000024") then
         include "leshouches_006.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."4-310000221000024") then
         include "leshouches_007.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."-3410000221000024") then
         include "leshouches_008.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
         
      elseif(str.eq."1-21000023-1000024") then
         include "leshouches_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-211000023-1000024") then
         include "leshouches_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."-431000023-1000024") then
         include "leshouches_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."3-41000023-1000024") then
         include "leshouches_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-1210000231000024") then
         include "leshouches_005.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."2-110000231000024") then
         include "leshouches_006.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."4-310000231000024") then
         include "leshouches_007.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."-3410000231000024") then
         include "leshouches_008.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20

      elseif(str.eq."1-21000022-1000037") then
         include "leshouches_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-211000022-1000037") then
         include "leshouches_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."-431000022-1000037") then
         include "leshouches_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."3-41000022-1000037") then
         include "leshouches_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-1210000221000037") then
         include "leshouches_005.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."2-110000221000037") then
         include "leshouches_006.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."4-310000221000037") then
         include "leshouches_007.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."-3410000221000037") then
         include "leshouches_008.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      
      elseif(str.eq."1-21000023-1000037") then ! n2 + x2
         include "leshouches_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-211000023-1000037") then
         include "leshouches_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."-431000023-1000037") then
         include "leshouches_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."3-41000023-1000037") then
         include "leshouches_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-1210000231000037") then
         include "leshouches_005.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."2-110000231000037") then
         include "leshouches_006.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."4-310000231000037") then
         include "leshouches_007.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."-3410000231000037") then
         include "leshouches_008.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      
      elseif(str.eq."1-21000025-1000024") then ! n3 + x1
         include "leshouches_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-211000025-1000024") then
         include "leshouches_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."-431000025-1000024") then
         include "leshouches_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."3-41000025-1000024") then
         include "leshouches_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-1210000251000024") then
         include "leshouches_005.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."2-110000251000024") then
         include "leshouches_006.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."4-310000251000024") then
         include "leshouches_007.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."-3410000251000024") then
         include "leshouches_008.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
         
      elseif(str.eq."1-21000025-1000037") then ! n3 + x2
         include "leshouches_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-211000025-1000037") then
         include "leshouches_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."-431000025-1000037") then
         include "leshouches_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."3-41000025-1000037") then
         include "leshouches_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-1210000251000037") then
         include "leshouches_005.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."2-110000251000037") then
         include "leshouches_006.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."4-310000251000037") then
         include "leshouches_007.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."-3410000251000037") then
         include "leshouches_008.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
         
      elseif(str.eq."1-21000035-1000024") then ! n4 + x1
         include "leshouches_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-211000035-1000024") then
         include "leshouches_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."-431000035-1000024") then
         include "leshouches_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."3-41000035-1000024") then
         include "leshouches_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-1210000351000024") then
         include "leshouches_005.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."2-110000351000024") then
         include "leshouches_006.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."4-310000351000024") then
         include "leshouches_007.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."-3410000351000024") then
         include "leshouches_008.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
         
      elseif(str.eq."1-21000035-1000037") then ! n4 + x2
         include "leshouches_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-211000035-1000037") then
         include "leshouches_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."-431000035-1000037") then
         include "leshouches_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."3-41000035-1000037") then
         include "leshouches_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-1210000351000037") then
         include "leshouches_005.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."2-110000351000037") then
         include "leshouches_006.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."4-310000351000037") then
         include "leshouches_007.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."-3410000351000037") then
         include "leshouches_008.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
         
      endif
      
 20   continue
      xtarget=jamp2cum(iflow)*random()
      ifl=1
      do while(jamp2cum(ifl).lt.xtarget)
         ifl=ifl+1
      enddo
      do i=1,2
         do j=1,nexternal-1
            color(i,j)=ICOLUP(i,j,ifl)
         enddo
      enddo
      
      end
      
      
      
      
      subroutine convert_to_string(npart,id,string,lstring)
      implicit none
      integer npart,lstring,i
      integer id(npart)
      character*140 string
      character*8 s
      
      do i=1,140
         string(i:i)=' '
      enddo
      lstring=0
      do i=1,npart
         if(id(i).eq.21) id(i)=0
         if(abs(id(i)).le.9) then
            s=char(abs(id(i))+48)
         elseif(abs(id(i)).le.99)then
            s=char(abs(id(i))/10+48)
     &           //char(mod(abs(id(i)),10)+48)
               elseif(abs(id(i)).le.999) then
                  s=char(abs(id(i))/100+48)
     &           //char((abs(id(i))-(abs(id(i))/100)*100)/10+48)
     &           //char(mod(abs(id(i))-(abs(id(i))/100)*100,10)+48)
         elseif(abs(id(i)).le.9999) then
            s=char(abs(id(i))/1000+48)
     &        //char((abs(id(i))-(abs(id(i))/1000)*1000)/100+48)
     &        //char((abs(id(i))-(abs(id(i))/100)*100)/10+48)
     &        //char(mod(abs(id(i))-(abs(id(i))/100)*100,10)+48)
         elseif(abs(id(i)).le.99999) then
           s=char(abs(id(i))/10000+48)
     &        //char((abs(id(i))-(abs(id(i))/10000)*10000)/1000+48)
     &        //char((abs(id(i))-(abs(id(i))/1000)*1000)/100+48)
     &        //char((abs(id(i))-(abs(id(i))/100)*100)/10+48)
     &        //char(mod(abs(id(i))-(abs(id(i))/100)*100,10)+48)
         elseif(abs(id(i)).le.999999) then
            s=char(abs(id(i))/100000+48)
     &        //char((abs(id(i))-(abs(id(i))/100000)*100000)/10000+48)
     &        //char((abs(id(i))-(abs(id(i))/10000)*10000)/1000+48)
     &        //char((abs(id(i))-(abs(id(i))/1000)*1000)/100+48)
     &        //char((abs(id(i))-(abs(id(i))/100)*100)/10+48)
     &        //char(mod(abs(id(i))-(abs(id(i))/100)*100,10)+48)
         elseif(abs(id(i)).le.9999999) then
            s=char(abs(id(i))/1000000+48)
     &       //char((abs(id(i))-(abs(id(i))/1000000)*1000000)/100000+48)
     &       //char((abs(id(i))-(abs(id(i))/100000)*100000)/10000+48)
     &       //char((abs(id(i))-(abs(id(i))/10000)*10000)/1000+48)
     &       //char((abs(id(i))-(abs(id(i))/1000)*1000)/100+48)
     &       //char((abs(id(i))-(abs(id(i))/100)*100)/10+48)
     &       //char(mod(abs(id(i))-(abs(id(i))/100)*100,10)+48)
         else
            write(*,*) 'error, particle ID is too large',abs(id(i))
         endif
         if(id(i).ge.0) then
            if(id(i).le.9) then
               string(lstring+1:lstring+1)=s
               lstring=lstring+1
            elseif(id(i).le.99) then
               string(lstring+1:lstring+2)=s
               lstring=lstring+2
            elseif(id(i).le.999) then
               string(lstring+1:lstring+3)=s
               lstring=lstring+3
            elseif(id(i).le.9999) then
              string(lstring+1:lstring+4)=s
              lstring=lstring+4
            elseif(id(i).le.99999) then
              string(lstring+1:lstring+5)=s
              lstring=lstring+5
            elseif(id(i).le.999999) then
              string(lstring+1:lstring+6)=s
              lstring=lstring+6
            elseif(id(i).le.9999999) then
              string(lstring+1:lstring+7)=s
              lstring=lstring+7
            endif
         else
            if(abs(id(i)).le.9) then
               string(lstring+1:lstring+2)='-'//s
               lstring=lstring+2
            elseif(abs(id(i)).le.99) then
               string(lstring+1:lstring+3)='-'//s
               lstring=lstring+3
            elseif(abs(id(i)).le.999) then
               string(lstring+1:lstring+4)='-'//s
               lstring=lstring+4
            elseif(abs(id(i)).le.9999) then
               string(lstring+1:lstring+5)='-'//s
               lstring=lstring+5
            elseif(abs(id(i)).le.99999) then
               string(lstring+1:lstring+6)='-'//s
               lstring=lstring+6
            elseif(abs(id(i)).le.999999) then
               string(lstring+1:lstring+7)='-'//s
               lstring=lstring+7
            elseif(abs(id(i)).le.9999999) then
               string(lstring+1:lstring+8)='-'//s
               lstring=lstring+8
            endif
         endif
      enddo
      end

      subroutine sborn_proc(p_born,legs,wgt,wgtjk,wgtmunu)
      implicit none
      include "nexternal.inc"
      include "coupl.inc"
      double precision wgt
      double precision p_born(0:3,nexternal-1),wgt2(nexternal-1),
     &   wgtmunu(0:3,0:3,nexternal-1),wgtjk(nexternal-1,nexternal-1)
      integer legs(nexternal-1),lstr,i
      character*140 str
      logical calculatedBorn
      integer skip
      common/cBorn/calculatedBorn,skip
 
      calculatedBorn=.false.
      
      call convert_to_string(nexternal-1,legs,str,lstr)
      
      if (str.eq."-1110000221000022") then
         call sborn_cl_ninj_001(p_born,wgtmunu,wgt2)
         call sborn_sf_ninj_001(p_born,wgtjk)
         goto 20
      elseif (str.eq."1-110000221000022") then
         call sborn_cl_ninj_002(p_born,wgtmunu,wgt2)
         call sborn_sf_ninj_002(p_born,wgtjk)
         goto 20
      elseif (str.eq."-2210000221000022") then
         call sborn_cl_ninj_003(p_born,wgtmunu,wgt2)
         call sborn_sf_ninj_003(p_born,wgtjk)
         goto 20
      elseif (str.eq."2-210000221000022") then
         call sborn_cl_ninj_004(p_born,wgtmunu,wgt2)
         call sborn_sf_ninj_004(p_born,wgtjk)
         goto 20
      elseif (str.eq."-4410000221000022") then
         call sborn_cl_ninj_005(p_born,wgtmunu,wgt2)
         call sborn_sf_ninj_005(p_born,wgtjk)
         goto 20
      elseif (str.eq."4-410000221000022") then
         call sborn_cl_ninj_006(p_born,wgtmunu,wgt2)
         call sborn_sf_ninj_006(p_born,wgtjk)
         goto 20
      elseif (str.eq."-3310000221000022") then
         call sborn_cl_ninj_007(p_born,wgtmunu,wgt2)
         call sborn_sf_ninj_007(p_born,wgtjk)
         goto 20
      elseif (str.eq."3-310000221000022") then
         call sborn_cl_ninj_008(p_born,wgtmunu,wgt2)
         call sborn_sf_ninj_008(p_born,wgtjk)
         goto 20
      elseif (str.eq."-5510000221000022") then
         call sborn_cl_ninj_009(p_born,wgtmunu,wgt2)
         call sborn_sf_ninj_009(p_born,wgtjk)
         goto 20
      elseif (str.eq."5-510000221000022") then
         call sborn_cl_ninj_010(p_born,wgtmunu,wgt2)
         call sborn_sf_ninj_010(p_born,wgtjk)
         goto 20
      
      elseif(str.eq."-1110000221000023") then
         call sborn_cl_nInJ_001(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_001(p_born,wgtjk)
         goto 20
      elseif(str.eq."1-110000221000023") then
         call sborn_cl_nInJ_002(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_002(p_born,wgtjk)
         goto 20
      elseif(str.eq."-2210000221000023") then
         call sborn_cl_nInJ_003(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_003(p_born,wgtjk)
         goto 20
      elseif(str.eq."2-210000221000023") then
         call sborn_cl_nInJ_004(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_004(p_born,wgtjk)
         goto 20
      elseif(str.eq."-4410000221000023") then
         call sborn_cl_nInJ_005(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_005(p_born,wgtjk)
         goto 20
      elseif(str.eq."4-410000221000023") then
         call sborn_cl_nInJ_006(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_006(p_born,wgtjk)
         goto 20
      elseif(str.eq."-3310000221000023") then
         call sborn_cl_nInJ_007(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_007(p_born,wgtjk)
         goto 20
      elseif(str.eq."3-310000221000023") then
         call sborn_cl_nInJ_008(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_008(p_born,wgtjk)
         goto 20
      elseif(str.eq."-5510000221000023") then
         call sborn_cl_nInJ_009(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_009(p_born,wgtjk)
         goto 20
      elseif(str.eq."5-510000221000023") then
         call sborn_cl_nInJ_010(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_010(p_born,wgtjk)
         goto 20
         
      elseif(str.eq."-1110000231000023") then
         call sborn_cl_ninj_001(p_born,wgtmunu,wgt2)
         call sborn_sf_ninj_001(p_born,wgtjk)
         goto 20
      elseif(str.eq."1-110000231000023") then
         call sborn_cl_ninj_002(p_born,wgtmunu,wgt2)
         call sborn_sf_ninj_002(p_born,wgtjk)
         goto 20
      elseif(str.eq."-2210000231000023") then
         call sborn_cl_ninj_003(p_born,wgtmunu,wgt2)
         call sborn_sf_ninj_003(p_born,wgtjk)
         goto 20
      elseif(str.eq."2-210000231000023") then
         call sborn_cl_ninj_004(p_born,wgtmunu,wgt2)
         call sborn_sf_ninj_004(p_born,wgtjk)
         goto 20
      elseif(str.eq."-4410000231000023") then
         call sborn_cl_ninj_005(p_born,wgtmunu,wgt2)
         call sborn_sf_ninj_005(p_born,wgtjk)
         goto 20
      elseif(str.eq."4-410000231000023") then
         call sborn_cl_ninj_006(p_born,wgtmunu,wgt2)
         call sborn_sf_ninj_006(p_born,wgtjk)
         goto 20
      elseif(str.eq."-3310000231000023") then
         call sborn_cl_ninj_007(p_born,wgtmunu,wgt2)
         call sborn_sf_ninj_007(p_born,wgtjk)
         goto 20
      elseif(str.eq."3-310000231000023") then
         call sborn_cl_ninj_008(p_born,wgtmunu,wgt2)
         call sborn_sf_ninj_008(p_born,wgtjk)
         goto 20
      elseif(str.eq."-5510000231000023") then
         call sborn_cl_ninj_009(p_born,wgtmunu,wgt2)
         call sborn_sf_ninj_009(p_born,wgtjk)
         goto 20
      elseif(str.eq."5-510000231000023") then
         call sborn_cl_ninj_010(p_born,wgtmunu,wgt2)
         call sborn_sf_ninj_010(p_born,wgtjk)
         goto 20
         
      elseif(str.eq."-1110000221000025") then
         call sborn_cl_nInJ_001(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_001(p_born,wgtjk)
         goto 20
      elseif(str.eq."1-110000221000025") then
         call sborn_cl_nInJ_002(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_002(p_born,wgtjk)
         goto 20
      elseif(str.eq."-2210000221000025") then
         call sborn_cl_nInJ_003(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_003(p_born,wgtjk)
         goto 20
      elseif(str.eq."2-210000221000025") then
         call sborn_cl_nInJ_004(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_004(p_born,wgtjk)
         goto 20
      elseif(str.eq."-4410000221000025") then
         call sborn_cl_nInJ_005(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_005(p_born,wgtjk)
         goto 20
      elseif(str.eq."4-410000221000025") then
         call sborn_cl_nInJ_006(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_006(p_born,wgtjk)
         goto 20
      elseif(str.eq."-3310000221000025") then
         call sborn_cl_nInJ_007(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_007(p_born,wgtjk)
         goto 20
      elseif(str.eq."3-310000221000025") then
         call sborn_cl_nInJ_008(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_008(p_born,wgtjk)
         goto 20
      elseif(str.eq."-5510000221000025") then
         call sborn_cl_nInJ_009(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_009(p_born,wgtjk)
         goto 20
      elseif(str.eq."5-510000221000025") then
         call sborn_cl_nInJ_010(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_010(p_born,wgtjk)
         goto 20
         
      elseif(str.eq."-1110000221000035") then
         call sborn_cl_nInJ_001(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_001(p_born,wgtjk)
         goto 20
      elseif(str.eq."1-110000221000035") then
         call sborn_cl_nInJ_002(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_002(p_born,wgtjk)
         goto 20
      elseif(str.eq."-2210000221000035") then
         call sborn_cl_nInJ_003(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_003(p_born,wgtjk)
         goto 20
      elseif(str.eq."2-210000221000035") then
         call sborn_cl_nInJ_004(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_004(p_born,wgtjk)
         goto 20
      elseif(str.eq."-4410000221000035") then
         call sborn_cl_nInJ_005(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_005(p_born,wgtjk)
         goto 20
      elseif(str.eq."4-410000221000035") then
         call sborn_cl_nInJ_006(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_006(p_born,wgtjk)
         goto 20
      elseif(str.eq."-3310000221000035") then
         call sborn_cl_nInJ_007(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_007(p_born,wgtjk)
         goto 20
      elseif(str.eq."3-310000221000035") then
         call sborn_cl_nInJ_008(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_008(p_born,wgtjk)
         goto 20
      elseif(str.eq."-5510000221000035") then
         call sborn_cl_nInJ_009(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_009(p_born,wgtjk)
         goto 20
      elseif(str.eq."5-510000221000035") then
         call sborn_cl_nInJ_010(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_010(p_born,wgtjk)
         goto 20
         
      elseif(str.eq."-1110000231000025") then
         call sborn_cl_nInJ_001(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_001(p_born,wgtjk)
         goto 20
      elseif(str.eq."1-110000231000025") then
         call sborn_cl_nInJ_002(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_002(p_born,wgtjk)
         goto 20
      elseif(str.eq."-2210000231000025") then
         call sborn_cl_nInJ_003(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_003(p_born,wgtjk)
         goto 20
      elseif(str.eq."2-210000231000025") then
         call sborn_cl_nInJ_004(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_004(p_born,wgtjk)
         goto 20
      elseif(str.eq."-4410000231000025") then
         call sborn_cl_nInJ_005(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_005(p_born,wgtjk)
         goto 20
      elseif(str.eq."4-410000231000025") then
         call sborn_cl_nInJ_006(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_006(p_born,wgtjk)
         goto 20
      elseif(str.eq."-3310000231000025") then
         call sborn_cl_nInJ_007(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_007(p_born,wgtjk)
         goto 20
      elseif(str.eq."3-310000231000025") then
         call sborn_cl_nInJ_008(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_008(p_born,wgtjk)
         goto 20
      elseif(str.eq."-5510000231000025") then
         call sborn_cl_nInJ_009(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_009(p_born,wgtjk)
         goto 20
      elseif(str.eq."5-510000231000025") then
         call sborn_cl_nInJ_010(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_010(p_born,wgtjk)
         goto 20
         
      elseif(str.eq."-1110000231000035") then
         call sborn_cl_nInJ_001(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_001(p_born,wgtjk)
         goto 20
      elseif(str.eq."1-110000231000035") then
         call sborn_cl_nInJ_002(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_002(p_born,wgtjk)
         goto 20
      elseif(str.eq."-2210000231000035") then
         call sborn_cl_nInJ_003(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_003(p_born,wgtjk)
         goto 20
      elseif(str.eq."2-210000231000035") then
         call sborn_cl_nInJ_004(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_004(p_born,wgtjk)
         goto 20
      elseif(str.eq."-4410000231000035") then
         call sborn_cl_nInJ_005(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_005(p_born,wgtjk)
         goto 20
      elseif(str.eq."4-410000231000035") then
         call sborn_cl_nInJ_006(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_006(p_born,wgtjk)
         goto 20
      elseif(str.eq."-3310000231000035") then
         call sborn_cl_nInJ_007(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_007(p_born,wgtjk)
         goto 20
      elseif(str.eq."3-310000231000035") then
         call sborn_cl_nInJ_008(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_008(p_born,wgtjk)
         goto 20
      elseif(str.eq."-5510000231000035") then
         call sborn_cl_nInJ_009(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_009(p_born,wgtjk)
         goto 20
      elseif(str.eq."5-510000231000035") then
         call sborn_cl_nInJ_010(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_010(p_born,wgtjk)
         goto 20
         
      elseif(str.eq."-1110000251000035") then
         call sborn_cl_nInJ_001(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_001(p_born,wgtjk)
         goto 20
      elseif(str.eq."1-110000251000035") then
         call sborn_cl_nInJ_002(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_002(p_born,wgtjk)
         goto 20
      elseif(str.eq."-2210000251000035") then
         call sborn_cl_nInJ_003(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_003(p_born,wgtjk)
         goto 20
      elseif(str.eq."2-210000251000035") then
         call sborn_cl_nInJ_004(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_004(p_born,wgtjk)
         goto 20
      elseif(str.eq."-4410000251000035") then
         call sborn_cl_nInJ_005(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_005(p_born,wgtjk)
         goto 20
      elseif(str.eq."4-410000251000035") then
         call sborn_cl_nInJ_006(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_006(p_born,wgtjk)
         goto 20
      elseif(str.eq."-3310000251000035") then
         call sborn_cl_nInJ_007(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_007(p_born,wgtjk)
         goto 20
      elseif(str.eq."3-310000251000035") then
         call sborn_cl_nInJ_008(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_008(p_born,wgtjk)
         goto 20
      elseif(str.eq."-5510000251000035") then
         call sborn_cl_nInJ_009(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_009(p_born,wgtjk)
         goto 20
      elseif(str.eq."5-510000251000035") then
         call sborn_cl_nInJ_010(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_010(p_born,wgtjk)
         goto 20
         
      elseif(str.eq."-1110000251000025") then
         call sborn_cl_nInJ_001(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_001(p_born,wgtjk)
         goto 20
      elseif(str.eq."1-110000251000025") then
         call sborn_cl_nInJ_002(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_002(p_born,wgtjk)
         goto 20
      elseif(str.eq."-2210000251000025") then
         call sborn_cl_nInJ_003(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_003(p_born,wgtjk)
         goto 20
      elseif(str.eq."2-210000251000025") then
         call sborn_cl_nInJ_004(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_004(p_born,wgtjk)
         goto 20
      elseif(str.eq."-4410000251000025") then
         call sborn_cl_nInJ_005(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_005(p_born,wgtjk)
         goto 20
      elseif(str.eq."4-410000251000025") then
         call sborn_cl_nInJ_006(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_006(p_born,wgtjk)
         goto 20
      elseif(str.eq."-3310000251000025") then
         call sborn_cl_nInJ_007(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_007(p_born,wgtjk)
         goto 20
      elseif(str.eq."3-310000251000025") then
         call sborn_cl_nInJ_008(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_008(p_born,wgtjk)
         goto 20
      elseif(str.eq."-5510000251000025") then
         call sborn_cl_nInJ_009(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_009(p_born,wgtjk)
         goto 20
      elseif(str.eq."5-510000251000025") then
         call sborn_cl_nInJ_010(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_010(p_born,wgtjk)
         goto 20
         
      elseif(str.eq."-1110000351000035") then
         call sborn_cl_nInJ_001(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_001(p_born,wgtjk)
         goto 20
      elseif(str.eq."1-110000351000035") then
         call sborn_cl_nInJ_002(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_002(p_born,wgtjk)
         goto 20
      elseif(str.eq."-2210000351000035") then
         call sborn_cl_nInJ_003(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_003(p_born,wgtjk)
         goto 20
      elseif(str.eq."2-210000351000035") then
         call sborn_cl_nInJ_004(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_004(p_born,wgtjk)
         goto 20
      elseif(str.eq."-4410000351000035") then
         call sborn_cl_nInJ_005(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_005(p_born,wgtjk)
         goto 20
      elseif(str.eq."4-410000351000035") then
         call sborn_cl_nInJ_006(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_006(p_born,wgtjk)
         goto 20
      elseif(str.eq."-3310000351000035") then
         call sborn_cl_nInJ_007(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_007(p_born,wgtjk)
         goto 20
      elseif(str.eq."3-310000351000035") then
         call sborn_cl_nInJ_008(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_008(p_born,wgtjk)
         goto 20
      elseif(str.eq."-5510000351000035") then
         call sborn_cl_nInJ_009(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_009(p_born,wgtjk)
         goto 20
      elseif(str.eq."5-510000351000035") then
         call sborn_cl_nInJ_010(p_born,wgtmunu,wgt2)
         call sborn_sf_nInJ_010(p_born,wgtjk)
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
      Double Precision amp2005(maxamps), jamp2005(0:maxamps)
      common/to_amps_005/amp2005,jamp2005
      Double Precision amp2006(maxamps), jamp2006(0:maxamps)
      common/to_amps_006/amp2006,jamp2006
      Double Precision amp2007(maxamps), jamp2007(0:maxamps)
      common/to_amps_007/amp2007,jamp2007
      Double Precision amp2008(maxamps), jamp2008(0:maxamps)
      common/to_amps_008/amp2008,jamp2008
      Double Precision amp2009(maxamps), jamp2009(0:maxamps)
      common/to_amps_009/amp2009,jamp2009
      Double Precision amp2010(maxamps), jamp2010(0:maxamps)
      common/to_amps_010/amp2010,jamp2010
      double precision jamp2cum(0:maxamps)
      integer ICOLUP(2,nexternal-1,maxamps)
      integer color(2,nexternal-1),color1(2,nexternal-1)
      double precision random,xtarget
      external random
      integer legs(nexternal-1),lstr,i,j
      character*140 str
      integer iflow,ifl
      
      call convert_to_string(nexternal-1,legs,str,lstr)
      
      if (str.eq."-1110000221000022") then
         include "leshouches_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif (str.eq."1-110000221000022") then
         include "leshouches_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif (str.eq."-2210000221000022") then
         include "leshouches_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif (str.eq."2-210000221000022") then
         include "leshouches_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif (str.eq."-4410000221000022") then
         include "leshouches_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif (str.eq."4-410000221000022") then
         include "leshouches_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif (str.eq."-3310000221000022") then
         include "leshouches_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif (str.eq."3-310000221000022") then
         include "leshouches_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif (str.eq."-5510000221000022") then
         include "leshouches_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif (str.eq."5-510000221000022") then
         include "leshouches_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      
      elseif(str.eq."-1110000221000023") then
         include "leshouches_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."1-110000221000023") then
         include "leshouches_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."-2210000221000023") then
         include "leshouches_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."2-210000221000023") then
         include "leshouches_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-4410000221000023") then
         include "leshouches_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."4-410000221000023") then
         include "leshouches_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."-3310000221000023") then
         include "leshouches_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."3-310000221000023") then
         include "leshouches_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."-5510000221000023") then
         include "leshouches_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."5-510000221000023") then
         include "leshouches_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
         
      elseif(str.eq."-1110000231000023") then
         include "leshouches_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."1-110000231000023") then
         include "leshouches_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."-2210000231000023") then
         include "leshouches_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."2-210000231000023") then
         include "leshouches_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-4410000231000023") then
         include "leshouches_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."4-410000231000023") then
         include "leshouches_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."-3310000231000023") then
         include "leshouches_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."3-310000231000023") then
         include "leshouches_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."-5510000231000023") then
         include "leshouches_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."5-510000231000023") then
         include "leshouches_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
         
      elseif(str.eq."-1110000221000025") then
         include "leshouches_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."1-110000221000025") then
         include "leshouches_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."-2210000221000025") then
         include "leshouches_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."2-210000221000025") then
         include "leshouches_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-4410000221000025") then
         include "leshouches_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."4-410000221000025") then
         include "leshouches_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."-3310000221000025") then
         include "leshouches_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."3-310000221000025") then
         include "leshouches_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."-5510000221000025") then
         include "leshouches_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."5-510000221000025") then
         include "leshouches_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
         
      elseif(str.eq."-1110000221000035") then
         include "leshouches_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."1-110000221000035") then
         include "leshouches_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."-2210000221000035") then
         include "leshouches_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."2-210000221000035") then
         include "leshouches_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-4410000221000035") then
         include "leshouches_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."4-410000221000035") then
         include "leshouches_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."-3310000221000035") then
         include "leshouches_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."3-310000221000035") then
         include "leshouches_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."-5510000221000035") then
         include "leshouches_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."5-510000221000035") then
         include "leshouches_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
         
      elseif(str.eq."-1110000231000025") then
         include "leshouches_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."1-110000231000025") then
         include "leshouches_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."-2210000231000025") then
         include "leshouches_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."2-210000231000025") then
         include "leshouches_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-4410000231000025") then
         include "leshouches_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."4-410000231000025") then
         include "leshouches_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."-3310000231000025") then
         include "leshouches_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."3-310000231000025") then
         include "leshouches_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."-5510000231000025") then
         include "leshouches_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."5-510000231000025") then
         include "leshouches_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
         
      elseif(str.eq."-1110000231000035") then
         include "leshouches_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."1-110000231000035") then
         include "leshouches_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."-2210000231000035") then
         include "leshouches_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."2-210000231000035") then
         include "leshouches_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-4410000231000035") then
         include "leshouches_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."4-410000231000035") then
         include "leshouches_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."-3310000231000035") then
         include "leshouches_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."3-310000231000035") then
         include "leshouches_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."-5510000231000035") then
         include "leshouches_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."5-510000231000035") then
         include "leshouches_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
         
      elseif(str.eq."-1110000251000025") then
         include "leshouches_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."1-110000251000025") then
         include "leshouches_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."-2210000251000025") then
         include "leshouches_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."2-210000251000025") then
         include "leshouches_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-4410000251000025") then
         include "leshouches_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."4-410000251000025") then
         include "leshouches_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."-3310000251000025") then
         include "leshouches_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."3-310000251000025") then
         include "leshouches_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."-5510000251000025") then
         include "leshouches_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."5-510000251000025") then
         include "leshouches_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
         
      elseif(str.eq."-1110000251000035") then
         include "leshouches_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."1-110000251000035") then
         include "leshouches_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."-2210000251000035") then
         include "leshouches_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."2-210000251000035") then
         include "leshouches_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-4410000251000035") then
         include "leshouches_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."4-410000251000035") then
         include "leshouches_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."-3310000251000035") then
         include "leshouches_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."3-310000251000035") then
         include "leshouches_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."-5510000251000035") then
         include "leshouches_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."5-510000251000035") then
         include "leshouches_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
         
      elseif(str.eq."-1110000351000035") then
         include "leshouches_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."1-110000351000035") then
         include "leshouches_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."-2210000351000035") then
         include "leshouches_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."2-210000351000035") then
         include "leshouches_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-4410000351000035") then
         include "leshouches_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."4-410000351000035") then
         include "leshouches_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."-3310000351000035") then
         include "leshouches_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."3-310000351000035") then
         include "leshouches_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."-5510000351000035") then
         include "leshouches_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."5-510000351000035") then
         include "leshouches_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
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

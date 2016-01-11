      subroutine sreal_proc(p,legs,wgt)
      implicit none
      include "nexternal.inc"
      include "coupl.inc"
      double precision p(0:3,nexternal),wgt
      integer legs(nexternal),lstr
      character*140 str
      
      call convert_to_string(nexternal,legs,str,lstr)
      
      if (str.eq."-11100002210000230") then
         call srealmtrx_001(p,wgt)
         goto 20
      elseif (str.eq."-1010000221000023-1") then
         call srealmtrx_002(p,wgt)
         goto 20
      elseif (str.eq."1-1100002210000230") then
         call srealmtrx_003(p,wgt)
         goto 20
      elseif (str.eq."10100002210000231") then
         call srealmtrx_004(p,wgt)
         goto 20
      elseif (str.eq."-22100002210000230") then
         call srealmtrx_005(p,wgt)
         goto 20
      elseif (str.eq."-2010000221000023-2") then
         call srealmtrx_006(p,wgt)
         goto 20
      elseif (str.eq."2-2100002210000230") then
         call srealmtrx_007(p,wgt)
         goto 20
      elseif (str.eq."20100002210000232") then
         call srealmtrx_008(p,wgt)
         goto 20
      elseif (str.eq."-44100002210000230") then
         call srealmtrx_009(p,wgt)
         goto 20
      elseif (str.eq."-4010000221000023-4") then
         call srealmtrx_010(p,wgt)
         goto 20
      elseif (str.eq."4-4100002210000230") then
         call srealmtrx_011(p,wgt)
         goto 20
      elseif (str.eq."40100002210000234") then
         call srealmtrx_012(p,wgt)
         goto 20
      elseif (str.eq."-33100002210000230") then
         call srealmtrx_013(p,wgt)
         goto 20
      elseif (str.eq."-3010000221000023-3") then
         call srealmtrx_014(p,wgt)
         goto 20
      elseif (str.eq."3-3100002210000230") then
         call srealmtrx_015(p,wgt)
         goto 20
      elseif (str.eq."30100002210000233") then
         call srealmtrx_016(p,wgt)
         goto 20
      elseif (str.eq."-55100002210000230") then
         call srealmtrx_017(p,wgt)
         goto 20
      elseif (str.eq."-5010000221000023-5") then
         call srealmtrx_018(p,wgt)
         goto 20
      elseif (str.eq."5-5100002210000230") then
         call srealmtrx_019(p,wgt)
         goto 20
      elseif (str.eq."50100002210000235") then
         call srealmtrx_020(p,wgt)
         goto 20
      elseif (str.eq."0-110000221000023-1") then
         call srealmtrx_021(p,wgt)
         goto 20
      elseif (str.eq."01100002210000231") then
         call srealmtrx_022(p,wgt)
         goto 20
      elseif (str.eq."0-210000221000023-2") then
         call srealmtrx_023(p,wgt)
         goto 20
      elseif (str.eq."02100002210000232") then
         call srealmtrx_024(p,wgt)
         goto 20
      elseif (str.eq."0-410000221000023-4") then
         call srealmtrx_025(p,wgt)
         goto 20
      elseif (str.eq."04100002210000234") then
         call srealmtrx_026(p,wgt)
         goto 20
      elseif (str.eq."0-310000221000023-3") then
         call srealmtrx_027(p,wgt)
         goto 20
      elseif (str.eq."03100002210000233") then
         call srealmtrx_028(p,wgt)
         goto 20
      elseif (str.eq."0-510000221000023-5") then
         call srealmtrx_029(p,wgt)
         goto 20
      elseif (str.eq."05100002210000235") then
         call srealmtrx_030(p,wgt)
         goto 20
      endif
      
 20   continue
      return
      end
      
      
      subroutine real_color(legs,color)
      implicit none
      include "nexternal.inc"
      integer maxamps
      parameter (maxamps=6000)
      Double Precision amp2001(maxamps), jamp2001(0:maxamps)
      common/to_Ramps_001/amp2001,jamp2001
      Double Precision amp2002(maxamps), jamp2002(0:maxamps)
      common/to_Ramps_002/amp2002,jamp2002
      Double Precision amp2003(maxamps), jamp2003(0:maxamps)
      common/to_Ramps_003/amp2003,jamp2003
      Double Precision amp2004(maxamps), jamp2004(0:maxamps)
      common/to_Ramps_004/amp2004,jamp2004
      Double Precision amp2005(maxamps), jamp2005(0:maxamps)
      common/to_Ramps_005/amp2005,jamp2005
      Double Precision amp2006(maxamps), jamp2006(0:maxamps)
      common/to_Ramps_006/amp2006,jamp2006
      Double Precision amp2007(maxamps), jamp2007(0:maxamps)
      common/to_Ramps_007/amp2007,jamp2007
      Double Precision amp2008(maxamps), jamp2008(0:maxamps)
      common/to_Ramps_008/amp2008,jamp2008
      Double Precision amp2009(maxamps), jamp2009(0:maxamps)
      common/to_Ramps_009/amp2009,jamp2009
      Double Precision amp2010(maxamps), jamp2010(0:maxamps)
      common/to_Ramps_010/amp2010,jamp2010
      Double Precision amp2011(maxamps), jamp2011(0:maxamps)
      common/to_Ramps_011/amp2011,jamp2011
      Double Precision amp2012(maxamps), jamp2012(0:maxamps)
      common/to_Ramps_012/amp2012,jamp2012
      Double Precision amp2013(maxamps), jamp2013(0:maxamps)
      common/to_Ramps_013/amp2013,jamp2013
      Double Precision amp2014(maxamps), jamp2014(0:maxamps)
      common/to_Ramps_014/amp2014,jamp2014
      Double Precision amp2015(maxamps), jamp2015(0:maxamps)
      common/to_Ramps_015/amp2015,jamp2015
      Double Precision amp2016(maxamps), jamp2016(0:maxamps)
      common/to_Ramps_016/amp2016,jamp2016
      Double Precision amp2017(maxamps), jamp2017(0:maxamps)
      common/to_Ramps_017/amp2017,jamp2017
      Double Precision amp2018(maxamps), jamp2018(0:maxamps)
      common/to_Ramps_018/amp2018,jamp2018
      Double Precision amp2019(maxamps), jamp2019(0:maxamps)
      common/to_Ramps_019/amp2019,jamp2019
      Double Precision amp2020(maxamps), jamp2020(0:maxamps)
      common/to_Ramps_020/amp2020,jamp2020
      Double Precision amp2021(maxamps), jamp2021(0:maxamps)
      common/to_Ramps_021/amp2021,jamp2021
      Double Precision amp2022(maxamps), jamp2022(0:maxamps)
      common/to_Ramps_022/amp2022,jamp2022
      Double Precision amp2023(maxamps), jamp2023(0:maxamps)
      common/to_Ramps_023/amp2023,jamp2023
      Double Precision amp2024(maxamps), jamp2024(0:maxamps)
      common/to_Ramps_024/amp2024,jamp2024
      Double Precision amp2025(maxamps), jamp2025(0:maxamps)
      common/to_Ramps_025/amp2025,jamp2025
      Double Precision amp2026(maxamps), jamp2026(0:maxamps)
      common/to_Ramps_026/amp2026,jamp2026
      Double Precision amp2027(maxamps), jamp2027(0:maxamps)
      common/to_Ramps_027/amp2027,jamp2027
      Double Precision amp2028(maxamps), jamp2028(0:maxamps)
      common/to_Ramps_028/amp2028,jamp2028
      Double Precision amp2029(maxamps), jamp2029(0:maxamps)
      common/to_Ramps_029/amp2029,jamp2029
      Double Precision amp2030(maxamps), jamp2030(0:maxamps)
      common/to_Ramps_030/amp2030,jamp2030
      double precision jamp2cum(0:maxamps)
      integer ICOLUP(2,nexternal,maxamps)
      integer color(2,nexternal),color1(2,nexternal)
      double precision random,xtarget
      external random
      integer legs(nexternal),lstr,i,j
      character*140 str
      integer iflow,ifl
      
      call convert_to_string(nexternal,legs,str,lstr)
      
      if (str.eq."-11100002210000230") then
         include "leshouches_R_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif (str.eq."-1010000221000023-1") then
         include "leshouches_R_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif (str.eq."1-1100002210000230") then
         include "leshouches_R_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif (str.eq."10100002210000231") then
         include "leshouches_R_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif (str.eq."-22100002210000230") then
         include "leshouches_R_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif (str.eq."-2010000221000023-2") then
         include "leshouches_R_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif (str.eq."2-2100002210000230") then
         include "leshouches_R_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif (str.eq."20100002210000232") then
         include "leshouches_R_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif (str.eq."-44100002210000230") then
         include "leshouches_R_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif (str.eq."-4010000221000023-4") then
         include "leshouches_R_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif (str.eq."4-4100002210000230") then
         include "leshouches_R_011.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif (str.eq."40100002210000234") then
         include "leshouches_R_012.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
      elseif (str.eq."-33100002210000230") then
         include "leshouches_R_013.inc"
         iflow=nint(jamp2013(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2013(i)
         enddo
         goto 20
      elseif (str.eq."-3010000221000023-3") then
         include "leshouches_R_014.inc"
         iflow=nint(jamp2014(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2014(i)
         enddo
         goto 20
      elseif (str.eq."3-3100002210000230") then
         include "leshouches_R_015.inc"
         iflow=nint(jamp2015(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2015(i)
         enddo
         goto 20
      elseif (str.eq."30100002210000233") then
         include "leshouches_R_016.inc"
         iflow=nint(jamp2016(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2016(i)
         enddo
         goto 20
      elseif (str.eq."-55100002210000230") then
         include "leshouches_R_017.inc"
         iflow=nint(jamp2017(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2017(i)
         enddo
         goto 20
      elseif (str.eq."-5010000221000023-5") then
         include "leshouches_R_018.inc"
         iflow=nint(jamp2018(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2018(i)
         enddo
         goto 20
      elseif (str.eq."5-5100002210000230") then
         include "leshouches_R_019.inc"
         iflow=nint(jamp2019(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2019(i)
         enddo
         goto 20
      elseif (str.eq."50100002210000235") then
         include "leshouches_R_020.inc"
         iflow=nint(jamp2020(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2020(i)
         enddo
         goto 20
      elseif (str.eq."0-110000221000023-1") then
         include "leshouches_R_021.inc"
         iflow=nint(jamp2021(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2021(i)
         enddo
         goto 20
      elseif (str.eq."01100002210000231") then
         include "leshouches_R_022.inc"
         iflow=nint(jamp2022(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2022(i)
         enddo
         goto 20
      elseif (str.eq."0-210000221000023-2") then
         include "leshouches_R_023.inc"
         iflow=nint(jamp2023(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2023(i)
         enddo
         goto 20
      elseif (str.eq."02100002210000232") then
         include "leshouches_R_024.inc"
         iflow=nint(jamp2024(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2024(i)
         enddo
         goto 20
      elseif (str.eq."0-410000221000023-4") then
         include "leshouches_R_025.inc"
         iflow=nint(jamp2025(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2025(i)
         enddo
         goto 20
      elseif (str.eq."04100002210000234") then
         include "leshouches_R_026.inc"
         iflow=nint(jamp2026(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2026(i)
         enddo
         goto 20
      elseif (str.eq."0-310000221000023-3") then
         include "leshouches_R_027.inc"
         iflow=nint(jamp2027(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2027(i)
         enddo
         goto 20
      elseif (str.eq."03100002210000233") then
         include "leshouches_R_028.inc"
         iflow=nint(jamp2028(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2028(i)
         enddo
         goto 20
      elseif (str.eq."0-510000221000023-5") then
         include "leshouches_R_029.inc"
         iflow=nint(jamp2029(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2029(i)
         enddo
         goto 20
      elseif (str.eq."05100002210000235") then
         include "leshouches_R_030.inc"
         iflow=nint(jamp2030(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2030(i)
         enddo
         goto 20
      endif
      
 20   continue
      xtarget=jamp2cum(iflow)*random()
      ifl=1
      do while (jamp2cum(ifl).lt.xtarget)
         ifl=ifl+1
      enddo
      do i=1,2
         do j=1,nexternal
            color(i,j)=ICOLUP(i,j,ifl)
         enddo
      enddo
      
      return
      end
      
      
      
      

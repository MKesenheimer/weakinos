      subroutine sreal_proc(p,legs,wgt)
      implicit none
      include "nexternal.inc"
      include "coupl.inc"
      double precision p(0:3,nexternal),wgt
      integer legs(nexternal),lstr
      character*140 str
      
      call convert_to_string(nexternal,legs,str,lstr)
      
#ifdef DEBUGQ
      print*, "[DEBUG] str = ", str
#endif
      
      if (str.eq."-111000024-10000240") then
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DXD_XIXIG(p,wgt)
         goto 20
      elseif (str.eq."-101000024-1000024-1") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DXG_XIXIDX(p,wgt)
         goto 20
      elseif (str.eq."1-11000024-10000240") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_DDX_XIXIG(p,wgt)
         goto 20
      elseif (str.eq."101000024-10000241") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_DG_XIXID(p,wgt)
         goto 20
      elseif (str.eq."-221000024-10000240") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UXU_XIXIG(p,wgt)
         goto 20
      elseif (str.eq."-201000024-1000024-2") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UXG_XIXIUX(p,wgt)
         goto 20
      elseif (str.eq."2-21000024-10000240") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_UUX_XIXIG(p,wgt)
         goto 20
      elseif (str.eq."201000024-10000242") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_UG_XIXIU(p,wgt)
         goto 20
      elseif (str.eq."-441000024-10000240") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_UXU_XIXIG(p,wgt)
         goto 20
      elseif (str.eq."-401000024-1000024-4") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_UXG_XIXIUX(p,wgt)
         goto 20
      elseif (str.eq."4-41000024-10000240") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_UUX_XIXIG(p,wgt)
         goto 20
      elseif (str.eq."401000024-10000244") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_UG_XIXIU(p,wgt)
         goto 20
      elseif (str.eq."-331000024-10000240") then
         !call srealmtrx_013(p,wgt)
         call SMATRIX_DXD_XIXIG(p,wgt)
         goto 20
      elseif (str.eq."-301000024-1000024-3") then
         !call srealmtrx_014(p,wgt)
         call SMATRIX_DXG_XIXIDX(p,wgt)
         goto 20
      elseif (str.eq."3-31000024-10000240") then
         !call srealmtrx_015(p,wgt)
         call SMATRIX_DDX_XIXIG(p,wgt)
         goto 20
      elseif (str.eq."301000024-10000243") then
         !call srealmtrx_016(p,wgt)
         call SMATRIX_DG_XIXID(p,wgt) 
         goto 20
      elseif (str.eq."-551000024-10000240") then
         !call srealmtrx_017(p,wgt)
#ifdef NEGLECTBMASS
          call SMATRIX_DXD_XIXIG(p,wgt)
#else
          print*,"ERROR: MadGraph b-quark amplitudes not yet "
                 //"implemented."
          ! The same for -501000024-1000024-5, 5-51000024-10000240
          ! 501000024-10000245, 0-51000024-1000024-5, 051000024-10000245
          stop
#endif
         goto 20
      elseif (str.eq."-501000024-1000024-5") then
         !call srealmtrx_018(p,wgt)
         call SMATRIX_DXG_XIXIDX(p,wgt)
         goto 20
      elseif (str.eq."5-51000024-10000240") then
         !call srealmtrx_019(p,wgt)
         call SMATRIX_DDX_XIXIG(p,wgt)
         goto 20
      elseif (str.eq."501000024-10000245") then
         !call srealmtrx_020(p,wgt)
         call SMATRIX_DG_XIXID(p,wgt)
         goto 20
      elseif (str.eq."0-11000024-1000024-1") then
         !call srealmtrx_021(p,wgt)
         call SMATRIX_GDX_XIXIDX(p,wgt)
         goto 20
      elseif (str.eq."011000024-10000241") then
         !call srealmtrx_022(p,wgt)
         call SMATRIX_GD_XIXID(p,wgt)
         goto 20
      elseif (str.eq."0-21000024-1000024-2") then
         !call srealmtrx_023(p,wgt)
         call SMATRIX_GUX_XIXIUX(p,wgt)
         goto 20
      elseif (str.eq."021000024-10000242") then
         !call srealmtrx_024(p,wgt)
         call SMATRIX_GU_XIXIU(p,wgt)
         goto 20
      elseif (str.eq."0-41000024-1000024-4") then
         !call srealmtrx_025(p,wgt)
         call SMATRIX_GUX_XIXIUX(p,wgt)
         goto 20
      elseif (str.eq."041000024-10000244") then
         !call srealmtrx_026(p,wgt)
         call SMATRIX_GU_XIXIU(p,wgt)
         goto 20
      elseif (str.eq."0-31000024-1000024-3") then
         !call srealmtrx_027(p,wgt)
         call SMATRIX_GDX_XIXIDX(p,wgt)
         goto 20
      elseif (str.eq."031000024-10000243") then
         !call srealmtrx_028(p,wgt)
         call SMATRIX_GD_XIXID(p,wgt)
         goto 20
      elseif (str.eq."0-51000024-1000024-5") then
         !call srealmtrx_029(p,wgt)
         call SMATRIX_GDX_XIXIDX(p,wgt)
         goto 20
      elseif (str.eq."051000024-10000245") then
         !call srealmtrx_030(p,wgt)
         call SMATRIX_GD_XIXID(p,wgt)
         goto 20
         
      elseif(str.eq."-111000024-10000370") then
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DXD_XIXJG(p,wgt)
         goto 20
      elseif(str.eq."-101000024-1000037-1") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DXG_XIXJDX(p,wgt)
         goto 20
      elseif(str.eq."1-11000024-10000370") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_DDX_XIXJG(p,wgt)
         goto 20
      elseif(str.eq."101000024-10000371") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_DG_XIXJD(p,wgt)
         goto 20
      elseif(str.eq."-221000024-10000370") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UXU_XIXJG(p,wgt)
         goto 20
      elseif(str.eq."-201000024-1000037-2") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UXG_XIXJUX(p,wgt)
         goto 20
      elseif(str.eq."2-21000024-10000370") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_UUX_XIXJG(p,wgt)
         goto 20
      elseif(str.eq."201000024-10000372") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_UG_XIXJU(p,wgt)
         goto 20
      elseif(str.eq."-441000024-10000370") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_UXU_XIXJG(p,wgt)
         goto 20
      elseif(str.eq."-401000024-1000037-4") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_UXG_XIXJUX(p,wgt)
         goto 20
      elseif(str.eq."4-41000024-10000370") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_UUX_XIXJG(p,wgt)
         goto 20
      elseif(str.eq."401000024-10000374") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_UG_XIXJU(p,wgt)
         goto 20
      elseif(str.eq."-331000024-10000370") then
         !call srealmtrx_013(p,wgt)
         call SMATRIX_DXD_XIXJG(p,wgt)
         goto 20
      elseif(str.eq."-301000024-1000037-3") then
         !call srealmtrx_014(p,wgt)
         call SMATRIX_DXG_XIXJDX(p,wgt)
         goto 20
      elseif(str.eq."3-31000024-10000370") then
         !call srealmtrx_015(p,wgt)
         call SMATRIX_DDX_XIXJG(p,wgt)
         goto 20
      elseif(str.eq."301000024-10000373") then
         !call srealmtrx_016(p,wgt)
         call SMATRIX_DG_XIXJD(p,wgt) 
         goto 20
      elseif(str.eq."-551000024-10000370") then
         !call srealmtrx_017(p,wgt)
#ifdef NEGLECTBMASS
          call SMATRIX_DXD_XIXJG(p,wgt)
#else
          print*,"ERROR: MadGraph b-quark amplitudes not yet "
                 //"implemented."
          ! The same for -501000024-1000037-5, 5-51000024-10000370
          ! 501000024-10000375, 0-51000024-1000037-5, 051000024-10000375
          stop
#endif
         goto 20
      elseif(str.eq."-501000024-1000037-5") then
         !call srealmtrx_018(p,wgt)
         call SMATRIX_DXG_XIXJDX(p,wgt)
         goto 20
      elseif(str.eq."5-51000024-10000370") then
         !call srealmtrx_019(p,wgt)
         call SMATRIX_DDX_XIXJG(p,wgt)
         goto 20
      elseif(str.eq."501000024-10000375") then
         !call srealmtrx_020(p,wgt)
         call SMATRIX_DG_XIXJD(p,wgt)
         goto 20
      elseif(str.eq."0-11000024-1000037-1") then
         !call srealmtrx_021(p,wgt)
         call SMATRIX_GDX_XIXJDX(p,wgt)
         goto 20
      elseif(str.eq."011000024-10000371") then
         !call srealmtrx_022(p,wgt)
         call SMATRIX_GD_XIXJD(p,wgt)
         goto 20
      elseif(str.eq."0-21000024-1000037-2") then
         !call srealmtrx_023(p,wgt)
         call SMATRIX_GUX_XIXJUX(p,wgt)
         goto 20
      elseif(str.eq."021000024-10000372") then
         !call srealmtrx_024(p,wgt)
         call SMATRIX_GU_XIXJU(p,wgt)
         goto 20
      elseif(str.eq."0-41000024-1000037-4") then
         !call srealmtrx_025(p,wgt)
         call SMATRIX_GUX_XIXJUX(p,wgt)
         goto 20
      elseif(str.eq."041000024-10000374") then
         !call srealmtrx_026(p,wgt)
         call SMATRIX_GU_XIXJU(p,wgt)
         goto 20
      elseif(str.eq."0-31000024-1000037-3") then
         !call srealmtrx_027(p,wgt)
         call SMATRIX_GDX_XIXJDX(p,wgt)
         goto 20
      elseif(str.eq."031000024-10000373") then
         !call srealmtrx_028(p,wgt)
         call SMATRIX_GD_XIXJD(p,wgt)
         goto 20
      elseif(str.eq."0-51000024-1000037-5") then
         !call srealmtrx_029(p,wgt)
         call SMATRIX_GDX_XIXJDX(p,wgt)
         goto 20
      elseif(str.eq."051000024-10000375") then
         !call srealmtrx_030(p,wgt)
         call SMATRIX_GD_XIXJD(p,wgt)
         goto 20
         
      elseif(str.eq."-111000037-10000370") then
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DXD_XIXIG(p,wgt)
         goto 20
      elseif(str.eq."-101000037-1000037-1") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DXG_XIXIDX(p,wgt)
         goto 20
      elseif(str.eq."1-11000037-10000370") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_DDX_XIXIG(p,wgt)
         goto 20
      elseif(str.eq."101000037-10000371") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_DG_XIXID(p,wgt)
         goto 20
      elseif(str.eq."-221000037-10000370") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UXU_XIXIG(p,wgt)
         goto 20
      elseif(str.eq."-201000037-1000037-2") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UXG_XIXIUX(p,wgt)
         goto 20
      elseif(str.eq."2-21000037-10000370") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_UUX_XIXIG(p,wgt)
         goto 20
      elseif(str.eq."201000037-10000372") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_UG_XIXIU(p,wgt)
         goto 20
      elseif(str.eq."-441000037-10000370") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_UXU_XIXIG(p,wgt)
         goto 20
      elseif(str.eq."-401000037-1000037-4") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_UXG_XIXIUX(p,wgt)
         goto 20
      elseif(str.eq."4-41000037-10000370") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_UUX_XIXIG(p,wgt)
         goto 20
      elseif(str.eq."401000037-10000374") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_UG_XIXIU(p,wgt)
         goto 20
      elseif(str.eq."-331000037-10000370") then
         !call srealmtrx_013(p,wgt)
         call SMATRIX_DXD_XIXIG(p,wgt)
         goto 20
      elseif(str.eq."-301000037-1000037-3") then
         !call srealmtrx_014(p,wgt)
         call SMATRIX_DXG_XIXIDX(p,wgt)
         goto 20
      elseif(str.eq."3-31000037-10000370") then
         !call srealmtrx_015(p,wgt)
         call SMATRIX_DDX_XIXIG(p,wgt)
         goto 20
      elseif(str.eq."301000037-10000373") then
         !call srealmtrx_016(p,wgt)
         call SMATRIX_DG_XIXID(p,wgt) 
         goto 20
      elseif(str.eq."-551000037-10000370") then
         !call srealmtrx_017(p,wgt)
#ifdef NEGLECTBMASS
          call SMATRIX_DXD_XIXIG(p,wgt)
#else
          print*,"ERROR: MadGraph b-quark amplitudes not yet "
                 //"implemented."
          ! The same for -501000037-1000037-5, 5-51000037-10000370
          ! 501000037-10000375, 0-51000037-1000037-5, 051000037-10000375
          stop
#endif
         goto 20
      elseif(str.eq."-501000037-1000037-5") then
         !call srealmtrx_018(p,wgt)
         call SMATRIX_DXG_XIXIDX(p,wgt)
         goto 20
      elseif(str.eq."5-51000037-10000370") then
         !call srealmtrx_019(p,wgt)
         call SMATRIX_DDX_XIXIG(p,wgt)
         goto 20
      elseif(str.eq."501000037-10000375") then
         !call srealmtrx_020(p,wgt)
         call SMATRIX_DG_XIXID(p,wgt)
         goto 20
      elseif(str.eq."0-11000037-1000037-1") then
         !call srealmtrx_021(p,wgt)
         call SMATRIX_GDX_XIXIDX(p,wgt)
         goto 20
      elseif(str.eq."011000037-10000371") then
         !call srealmtrx_022(p,wgt)
         call SMATRIX_GD_XIXID(p,wgt)
         goto 20
      elseif(str.eq."0-21000037-1000037-2") then
         !call srealmtrx_023(p,wgt)
         call SMATRIX_GUX_XIXIUX(p,wgt)
         goto 20
      elseif(str.eq."021000037-10000372") then
         !call srealmtrx_024(p,wgt)
         call SMATRIX_GU_XIXIU(p,wgt)
         goto 20
      elseif(str.eq."0-41000037-1000037-4") then
         !call srealmtrx_025(p,wgt)
         call SMATRIX_GUX_XIXIUX(p,wgt)
         goto 20
      elseif(str.eq."041000037-10000374") then
         !call srealmtrx_026(p,wgt)
         call SMATRIX_GU_XIXIU(p,wgt)
         goto 20
      elseif(str.eq."0-31000037-1000037-3") then
         !call srealmtrx_027(p,wgt)
         call SMATRIX_GDX_XIXIDX(p,wgt)
         goto 20
      elseif(str.eq."031000037-10000373") then
         !call srealmtrx_028(p,wgt)
         call SMATRIX_GD_XIXID(p,wgt)
         goto 20
      elseif(str.eq."0-51000037-1000037-5") then
         !call srealmtrx_029(p,wgt)
         call SMATRIX_GDX_XIXIDX(p,wgt)
         goto 20
      elseif(str.eq."051000037-10000375") then
         !call srealmtrx_030(p,wgt)
         call SMATRIX_GD_XIXID(p,wgt)
         goto 20
      endif
      
 20   continue
      if(wgt .eq. 0D0) then
          print*,"WARNING: maybe error in real amplitudes: ", wgt
          stop
      endif

      return
      end
      
      
      subroutine real_color(legs,color)
      implicit none
      include "nexternal.inc"
      include "maxamps_xixj.inc"
      Double Precision amp2001(maxamps), jamp2001(0:maxflow)
      common/to_amps_dxd/amp2001,jamp2001
      Double Precision amp2002(maxamps), jamp2002(0:maxflow)
      common/to_amps_dxg/amp2002,jamp2002
      Double Precision amp2003(maxamps), jamp2003(0:maxflow)
      common/to_amps_ddx/amp2003,jamp2003
      Double Precision amp2004(maxamps), jamp2004(0:maxflow)
      common/to_amps_dg/amp2004,jamp2004
      Double Precision amp2005(maxamps), jamp2005(0:maxflow)
      common/to_amps_uxu/amp2005,jamp2005
      Double Precision amp2006(maxamps), jamp2006(0:maxflow)
      common/to_amps_uxg/amp2006,jamp2006
      Double Precision amp2007(maxamps), jamp2007(0:maxflow)
      common/to_amps_uux/amp2007,jamp2007
      Double Precision amp2008(maxamps), jamp2008(0:maxflow)
      common/to_amps_ug/amp2008,jamp2008
      Double Precision amp2009(maxamps), jamp2009(0:maxflow)
      !common/to_amps_uxu/amp2009,jamp2009
      Double Precision amp2010(maxamps), jamp2010(0:maxflow)
      !common/to_amps_uxg/amp2010,jamp2010
      Double Precision amp2011(maxamps), jamp2011(0:maxflow)
      !common/to_amps_uux/amp2011,jamp2011
      Double Precision amp2012(maxamps), jamp2012(0:maxflow)
      !common/to_amps_ug/amp2012,jamp2012
      Double Precision amp2013(maxamps), jamp2013(0:maxflow)
      !common/to_amps_dxd/amp2013,jamp2013
      Double Precision amp2014(maxamps), jamp2014(0:maxflow)
      !common/to_amps_dxg/amp2014,jamp2014
      Double Precision amp2015(maxamps), jamp2015(0:maxflow)
      !common/to_amps_ddx/amp2015,jamp2015
      Double Precision amp2016(maxamps), jamp2016(0:maxflow)
      !common/to_amps_dg/amp2016,jamp2016
      Double Precision amp2017(maxamps), jamp2017(0:maxflow)
      !common/to_amps_dxd/amp2017,jamp2017
      Double Precision amp2018(maxamps), jamp2018(0:maxflow)
      !common/to_amps_dxg/amp2018,jamp2018
      Double Precision amp2019(maxamps), jamp2019(0:maxflow)
      !common/to_amps_ddx/amp2019,jamp2019
      Double Precision amp2020(maxamps), jamp2020(0:maxflow)
      !common/to_amps_dg/amp2020,jamp2020
      Double Precision amp2021(maxamps), jamp2021(0:maxflow)
      common/to_amps_gdx/amp2021,jamp2021
      Double Precision amp2022(maxamps), jamp2022(0:maxflow)
      common/to_amps_gd/amp2022,jamp2022
      Double Precision amp2023(maxamps), jamp2023(0:maxflow)
      common/to_amps_gux/amp2023,jamp2023
      Double Precision amp2024(maxamps), jamp2024(0:maxflow)
      common/to_amps_gu/amp2024,jamp2024
      Double Precision amp2025(maxamps), jamp2025(0:maxflow)
      !common/to_amps_gux/amp2025,jamp2025
      Double Precision amp2026(maxamps), jamp2026(0:maxflow)
      !common/to_amps_gu/amp2026,jamp2026
      Double Precision amp2027(maxamps), jamp2027(0:maxflow)
      !common/to_amps_gdx/amp2027,jamp2027
      Double Precision amp2028(maxamps), jamp2028(0:maxflow)
      !common/to_amps_gd/amp2028,jamp2028
      Double Precision amp2029(maxamps), jamp2029(0:maxflow)
      !common/to_amps_gdx/amp2029,jamp2029
      Double Precision amp2030(maxamps), jamp2030(0:maxflow)
      !common/to_amps_gd/amp2030,jamp2030
      
      equivalence (amp2001,amp2013,amp2017)
      equivalence (jamp2001,jamp2013,jamp2017)
      
      equivalence (amp2002,amp2014,amp2018)
      equivalence (jamp2002,jamp2014,jamp2018)
      
      equivalence (amp2003,amp2015,amp2019)
      equivalence (jamp2003,jamp2015,jamp2019)
      
      equivalence (amp2004,amp2016,amp2020)
      equivalence (jamp2004,jamp2016,jamp2020)
      
      equivalence (amp2005,amp2009)
      equivalence (jamp2005,jamp2009)
      
      equivalence (amp2006,amp2010)
      equivalence (jamp2006,jamp2010)
      
      equivalence (amp2007,amp2011)
      equivalence (jamp2007,jamp2011)
      
      equivalence (amp2008,amp2012)
      equivalence (jamp2008,jamp2012)
      
      equivalence (amp2021,amp2027,amp2029)
      equivalence (jamp2021,jamp2027,jamp2029)
      
      equivalence (amp2022,amp2028,amp2030)
      equivalence (jamp2022,jamp2028,jamp2030)
      
      equivalence (amp2023,amp2025)
      equivalence (jamp2023,jamp2025)
      
      equivalence (amp2024,amp2026)
      equivalence (jamp2024,jamp2026)
      
      double precision jamp2cum(0:maxamps)
      integer ICOLUP(2,nexternal,maxamps)
      integer color(2,nexternal),color1(2,nexternal)
      double precision random,xtarget
      external random
      integer legs(nexternal),lstr,i,j
      character*140 str
      integer iflow,ifl
      
      call convert_to_string(nexternal,legs,str,lstr)

      if (str.eq."-111000024-10000240") then
         include "leshouches_R_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif (str.eq."-101000024-1000024-1") then
         include "leshouches_R_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif (str.eq."1-11000024-10000240") then
         include "leshouches_R_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif (str.eq."101000024-10000241") then
         include "leshouches_R_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif (str.eq."-221000024-10000240") then
         include "leshouches_R_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif (str.eq."-201000024-1000024-2") then
         include "leshouches_R_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif (str.eq."2-21000024-10000240") then
         include "leshouches_R_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif (str.eq."201000024-10000242") then
         include "leshouches_R_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif (str.eq."-441000024-10000240") then
         include "leshouches_R_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif (str.eq."-401000024-1000024-4") then
         include "leshouches_R_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif (str.eq."4-41000024-10000240") then
         include "leshouches_R_011.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif (str.eq."401000024-10000244") then
         include "leshouches_R_012.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
      elseif (str.eq."-331000024-10000240") then
         include "leshouches_R_013.inc"
         iflow=nint(jamp2013(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2013(i)
         enddo
         goto 20
      elseif (str.eq."-301000024-1000024-3") then
         include "leshouches_R_014.inc"
         iflow=nint(jamp2014(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2014(i)
         enddo
         goto 20
      elseif (str.eq."3-31000024-10000240") then
         include "leshouches_R_015.inc"
         iflow=nint(jamp2015(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2015(i)
         enddo
         goto 20
      elseif (str.eq."301000024-10000243") then
         include "leshouches_R_016.inc"
         iflow=nint(jamp2016(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2016(i)
         enddo
         goto 20
      elseif (str.eq."-551000024-10000240") then
         include "leshouches_R_017.inc"
         iflow=nint(jamp2017(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2017(i)
         enddo
         goto 20
      elseif (str.eq."-501000024-1000024-5") then
         include "leshouches_R_018.inc"
         iflow=nint(jamp2018(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2018(i)
         enddo
         goto 20
      elseif (str.eq."5-51000024-10000240") then
         include "leshouches_R_019.inc"
         iflow=nint(jamp2019(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2019(i)
         enddo
         goto 20
      elseif (str.eq."501000024-10000245") then
         include "leshouches_R_020.inc"
         iflow=nint(jamp2020(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2020(i)
         enddo
         goto 20
      elseif (str.eq."0-11000024-1000024-1") then
         include "leshouches_R_021.inc"
         iflow=nint(jamp2021(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2021(i)
         enddo
         goto 20
      elseif (str.eq."011000024-10000241") then
         include "leshouches_R_022.inc"
         iflow=nint(jamp2022(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2022(i)
         enddo
         goto 20
      elseif (str.eq."0-21000024-1000024-2") then
         include "leshouches_R_023.inc"
         iflow=nint(jamp2023(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2023(i)
         enddo
         goto 20
      elseif (str.eq."021000024-10000242") then
         include "leshouches_R_024.inc"
         iflow=nint(jamp2024(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2024(i)
         enddo
         goto 20
      elseif (str.eq."0-41000024-1000024-4") then
         include "leshouches_R_025.inc"
         iflow=nint(jamp2025(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2025(i)
         enddo
         goto 20
      elseif (str.eq."041000024-10000244") then
         include "leshouches_R_026.inc"
         iflow=nint(jamp2026(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2026(i)
         enddo
         goto 20
      elseif (str.eq."0-31000024-1000024-3") then
         include "leshouches_R_027.inc"
         iflow=nint(jamp2027(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2027(i)
         enddo
         goto 20
      elseif (str.eq."031000024-10000243") then
         include "leshouches_R_028.inc"
         iflow=nint(jamp2028(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2028(i)
         enddo
         goto 20
      elseif (str.eq."0-51000024-1000024-5") then
         include "leshouches_R_029.inc"
         iflow=nint(jamp2029(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2029(i)
         enddo
         goto 20
      elseif (str.eq."051000024-10000245") then
         include "leshouches_R_030.inc"
         iflow=nint(jamp2030(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2030(i)
         enddo
         goto 20
         
      elseif(str.eq."-111000024-10000370") then
         include "leshouches_R_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-101000024-1000037-1") then
         include "leshouches_R_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."1-11000024-10000370") then
         include "leshouches_R_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."101000024-10000371") then
         include "leshouches_R_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-221000024-10000370") then
         include "leshouches_R_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."-201000024-1000037-2") then
         include "leshouches_R_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."2-21000024-10000370") then
         include "leshouches_R_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."201000024-10000372") then
         include "leshouches_R_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."-441000024-10000370") then
         include "leshouches_R_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."-401000024-1000037-4") then
         include "leshouches_R_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."4-41000024-10000370") then
         include "leshouches_R_011.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."401000024-10000374") then
         include "leshouches_R_012.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
      elseif(str.eq."-331000024-10000370") then
         include "leshouches_R_013.inc"
         iflow=nint(jamp2013(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2013(i)
         enddo
         goto 20
      elseif(str.eq."-301000024-1000037-3") then
         include "leshouches_R_014.inc"
         iflow=nint(jamp2014(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2014(i)
         enddo
         goto 20
      elseif(str.eq."3-31000024-10000370") then
         include "leshouches_R_015.inc"
         iflow=nint(jamp2015(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2015(i)
         enddo
         goto 20
      elseif(str.eq."301000024-10000373") then
         include "leshouches_R_016.inc"
         iflow=nint(jamp2016(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2016(i)
         enddo
         goto 20
      elseif(str.eq."-551000024-10000370") then
         include "leshouches_R_017.inc"
         iflow=nint(jamp2017(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2017(i)
         enddo
         goto 20
      elseif(str.eq."-501000024-1000037-5") then
         include "leshouches_R_018.inc"
         iflow=nint(jamp2018(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2018(i)
         enddo
         goto 20
      elseif(str.eq."5-51000024-10000370") then
         include "leshouches_R_019.inc"
         iflow=nint(jamp2019(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2019(i)
         enddo
         goto 20
      elseif(str.eq."501000024-10000375") then
         include "leshouches_R_020.inc"
         iflow=nint(jamp2020(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2020(i)
         enddo
         goto 20
      elseif(str.eq."0-11000024-1000037-1") then
         include "leshouches_R_021.inc"
         iflow=nint(jamp2021(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2021(i)
         enddo
         goto 20
      elseif(str.eq."011000024-10000371") then
         include "leshouches_R_022.inc"
         iflow=nint(jamp2022(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2022(i)
         enddo
         goto 20
      elseif(str.eq."0-21000024-1000037-2") then
         include "leshouches_R_023.inc"
         iflow=nint(jamp2023(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2023(i)
         enddo
         goto 20
      elseif(str.eq."021000024-10000372") then
         include "leshouches_R_024.inc"
         iflow=nint(jamp2024(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2024(i)
         enddo
         goto 20
      elseif(str.eq."0-41000024-1000037-4") then
         include "leshouches_R_025.inc"
         iflow=nint(jamp2025(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2025(i)
         enddo
         goto 20
      elseif(str.eq."041000024-10000374") then
         include "leshouches_R_026.inc"
         iflow=nint(jamp2026(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2026(i)
         enddo
         goto 20
      elseif(str.eq."0-31000024-1000037-3") then
         include "leshouches_R_027.inc"
         iflow=nint(jamp2027(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2027(i)
         enddo
         goto 20
      elseif(str.eq."031000024-10000373") then
         include "leshouches_R_028.inc"
         iflow=nint(jamp2028(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2028(i)
         enddo
         goto 20
      elseif(str.eq."0-51000024-1000037-5") then
         include "leshouches_R_029.inc"
         iflow=nint(jamp2029(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2029(i)
         enddo
         goto 20
      elseif(str.eq."051000024-10000375") then
         include "leshouches_R_030.inc"
         iflow=nint(jamp2030(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2030(i)
         enddo
         goto 20
         
      elseif(str.eq."-111000037-10000370") then
         include "leshouches_R_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-101000037-1000037-1") then
         include "leshouches_R_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."1-11000037-10000370") then
         include "leshouches_R_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."101000037-10000371") then
         include "leshouches_R_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-221000037-10000370") then
         include "leshouches_R_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."-201000037-1000037-2") then
         include "leshouches_R_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."2-21000037-10000370") then
         include "leshouches_R_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."201000037-10000372") then
         include "leshouches_R_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."-441000037-10000370") then
         include "leshouches_R_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."-401000037-1000037-4") then
         include "leshouches_R_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."4-41000037-10000370") then
         include "leshouches_R_011.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."401000037-10000374") then
         include "leshouches_R_012.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
      elseif(str.eq."-331000037-10000370") then
         include "leshouches_R_013.inc"
         iflow=nint(jamp2013(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2013(i)
         enddo
         goto 20
      elseif(str.eq."-301000037-1000037-3") then
         include "leshouches_R_014.inc"
         iflow=nint(jamp2014(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2014(i)
         enddo
         goto 20
      elseif(str.eq."3-31000037-10000370") then
         include "leshouches_R_015.inc"
         iflow=nint(jamp2015(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2015(i)
         enddo
         goto 20
      elseif(str.eq."301000037-10000373") then
         include "leshouches_R_016.inc"
         iflow=nint(jamp2016(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2016(i)
         enddo
         goto 20
      elseif(str.eq."-551000037-10000370") then
         include "leshouches_R_017.inc"
         iflow=nint(jamp2017(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2017(i)
         enddo
         goto 20
      elseif(str.eq."-501000037-1000037-5") then
         include "leshouches_R_018.inc"
         iflow=nint(jamp2018(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2018(i)
         enddo
         goto 20
      elseif(str.eq."5-51000037-10000370") then
         include "leshouches_R_019.inc"
         iflow=nint(jamp2019(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2019(i)
         enddo
         goto 20
      elseif(str.eq."501000037-10000375") then
         include "leshouches_R_020.inc"
         iflow=nint(jamp2020(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2020(i)
         enddo
         goto 20
      elseif(str.eq."0-11000037-1000037-1") then
         include "leshouches_R_021.inc"
         iflow=nint(jamp2021(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2021(i)
         enddo
         goto 20
      elseif(str.eq."011000037-10000371") then
         include "leshouches_R_022.inc"
         iflow=nint(jamp2022(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2022(i)
         enddo
         goto 20
      elseif(str.eq."0-21000037-1000037-2") then
         include "leshouches_R_023.inc"
         iflow=nint(jamp2023(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2023(i)
         enddo
         goto 20
      elseif(str.eq."021000037-10000372") then
         include "leshouches_R_024.inc"
         iflow=nint(jamp2024(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2024(i)
         enddo
         goto 20
      elseif(str.eq."0-41000037-1000037-4") then
         include "leshouches_R_025.inc"
         iflow=nint(jamp2025(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2025(i)
         enddo
         goto 20
      elseif(str.eq."041000037-10000374") then
         include "leshouches_R_026.inc"
         iflow=nint(jamp2026(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2026(i)
         enddo
         goto 20
      elseif(str.eq."0-31000037-1000037-3") then
         include "leshouches_R_027.inc"
         iflow=nint(jamp2027(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2027(i)
         enddo
         goto 20
      elseif(str.eq."031000037-10000373") then
         include "leshouches_R_028.inc"
         iflow=nint(jamp2028(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2028(i)
         enddo
         goto 20
      elseif(str.eq."0-51000037-1000037-5") then
         include "leshouches_R_029.inc"
         iflow=nint(jamp2029(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2029(i)
         enddo
         goto 20
      elseif(str.eq."051000037-10000375") then
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


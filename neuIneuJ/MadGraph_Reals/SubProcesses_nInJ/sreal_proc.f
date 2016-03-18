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
      
      if (str.eq."-11100002210000220") then
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DXD_NINJG(p,wgt)
         goto 20
      elseif (str.eq."-1010000221000022-1") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif (str.eq."1-1100002210000220") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif (str.eq."10100002210000221") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt)
         goto 20
      elseif (str.eq."-22100002210000220") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UXU_NINJG(p,wgt)
         goto 20
      elseif (str.eq."-2010000221000022-2") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UXG_NINJUX(p,wgt)
         goto 20
      elseif (str.eq."2-2100002210000220") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_UUX_NINJG(p,wgt)
         goto 20
      elseif (str.eq."20100002210000222") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_UG_NINJU(p,wgt)
         goto 20
      elseif (str.eq."-44100002210000220") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_UXU_NINJG(p,wgt)
         goto 20
      elseif (str.eq."-4010000221000022-4") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_UXG_NINJUX(p,wgt)
         goto 20
      elseif (str.eq."4-4100002210000220") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_UUX_NINJG(p,wgt)
         goto 20
      elseif (str.eq."40100002210000224") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_UG_NINJU(p,wgt)
         goto 20
      elseif (str.eq."-33100002210000220") then
         !call srealmtrx_013(p,wgt)
         call SMATRIX_DXD_NINJG(p,wgt)
         goto 20
      elseif (str.eq."-3010000221000022-3") then
         !call srealmtrx_014(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif (str.eq."3-3100002210000220") then
         !call srealmtrx_015(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif (str.eq."30100002210000223") then
         !call srealmtrx_016(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt) 
         goto 20
      elseif (str.eq."-55100002210000220") then
         !call srealmtrx_017(p,wgt)
#ifdef NEGLECTBMASS
          call SMATRIX_DXD_NINJG(p,wgt)
#else
          print*,"ERROR: MadGraph b-quark amplitudes not yet "
                 //"implemented."
          ! The same for -5010000221000022-5, 5-5100002210000220
          ! 50100002210000225, 0-510000221000022-5, 05100002210000225
          stop
#endif
         goto 20
      elseif (str.eq."-5010000221000022-5") then
         !call srealmtrx_018(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif (str.eq."5-5100002210000220") then
         !call srealmtrx_019(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif (str.eq."50100002210000225") then
         !call srealmtrx_020(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt)
         goto 20
      elseif (str.eq."0-110000221000022-1") then
         !call srealmtrx_021(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif (str.eq."01100002210000221") then
         !call srealmtrx_022(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
      elseif (str.eq."0-210000221000022-2") then
         !call srealmtrx_023(p,wgt)
         call SMATRIX_GUX_NINJUX(p,wgt)
         goto 20
      elseif (str.eq."02100002210000222") then
         !call srealmtrx_024(p,wgt)
         call SMATRIX_GU_NINJU(p,wgt)
         goto 20
      elseif (str.eq."0-410000221000022-4") then
         !call srealmtrx_025(p,wgt)
         call SMATRIX_GUX_NINJUX(p,wgt)
         goto 20
      elseif (str.eq."04100002210000224") then
         !call srealmtrx_026(p,wgt)
         call SMATRIX_GU_NINJU(p,wgt)
         goto 20
      elseif (str.eq."0-310000221000022-3") then
         !call srealmtrx_027(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif (str.eq."03100002210000223") then
         !call srealmtrx_028(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
      elseif (str.eq."0-510000221000022-5") then
         !call srealmtrx_029(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif (str.eq."05100002210000225") then
         !call srealmtrx_030(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
         
      elseif(str.eq."-11100002210000230") then
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DXD_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-1010000221000023-1") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."1-1100002210000230") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."10100002210000231") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt)
         goto 20
      elseif(str.eq."-22100002210000230") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UXU_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-2010000221000023-2") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UXG_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."2-2100002210000230") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_UUX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."20100002210000232") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_UG_NINJU(p,wgt)
         goto 20
      elseif(str.eq."-44100002210000230") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_UXU_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-4010000221000023-4") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_UXG_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."4-4100002210000230") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_UUX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."40100002210000234") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_UG_NINJU(p,wgt)
         goto 20
      elseif(str.eq."-33100002210000230") then
         !call srealmtrx_013(p,wgt)
         call SMATRIX_DXD_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-3010000221000023-3") then
         !call srealmtrx_014(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."3-3100002210000230") then
         !call srealmtrx_015(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."30100002210000233") then
         !call srealmtrx_016(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt) 
         goto 20
      elseif(str.eq."-55100002210000230") then
         !call srealmtrx_017(p,wgt)
#ifdef NEGLECTBMASS
          call SMATRIX_DXD_NINJG(p,wgt)
#else
          print*,"ERROR: MadGraph b-quark amplitudes not yet "
                 //"implemented."
          ! The same for -5010000221000023-5, 5-5100002210000230
          ! 50100002210000235, 0-510000221000023-5, 05100002210000235
          stop
#endif
         goto 20
      elseif(str.eq."-5010000221000023-5") then
         !call srealmtrx_018(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."5-5100002210000230") then
         !call srealmtrx_019(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."50100002210000235") then
         !call srealmtrx_020(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-110000221000023-1") then
         !call srealmtrx_021(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."01100002210000231") then
         !call srealmtrx_022(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-210000221000023-2") then
         !call srealmtrx_023(p,wgt)
         call SMATRIX_GUX_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."02100002210000232") then
         !call srealmtrx_024(p,wgt)
         call SMATRIX_GU_NINJU(p,wgt)
         goto 20
      elseif(str.eq."0-410000221000023-4") then
         !call srealmtrx_025(p,wgt)
         call SMATRIX_GUX_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."04100002210000234") then
         !call srealmtrx_026(p,wgt)
         call SMATRIX_GU_NINJU(p,wgt)
         goto 20
      elseif(str.eq."0-310000221000023-3") then
         !call srealmtrx_027(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."03100002210000233") then
         !call srealmtrx_028(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-510000221000023-5") then
         !call srealmtrx_029(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."05100002210000235") then
         !call srealmtrx_030(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
         
      elseif(str.eq."-11100002310000230") then
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DXD_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-1010000231000023-1") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."1-1100002310000230") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."10100002310000231") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt)
         goto 20
      elseif(str.eq."-22100002310000230") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UXU_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-2010000231000023-2") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UXG_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."2-2100002310000230") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_UUX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."20100002310000232") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_UG_NINJU(p,wgt)
         goto 20
      elseif(str.eq."-44100002310000230") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_UXU_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-4010000231000023-4") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_UXG_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."4-4100002310000230") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_UUX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."40100002310000234") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_UG_NINJU(p,wgt)
         goto 20
      elseif(str.eq."-33100002310000230") then
         !call srealmtrx_013(p,wgt)
         call SMATRIX_DXD_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-3010000231000023-3") then
         !call srealmtrx_014(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."3-3100002310000230") then
         !call srealmtrx_015(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."30100002310000233") then
         !call srealmtrx_016(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt) 
         goto 20
      elseif(str.eq."-55100002310000230") then
         !call srealmtrx_017(p,wgt)
#ifdef NEGLECTBMASS
          call SMATRIX_DXD_NINJG(p,wgt)
#else
          print*,"ERROR: MadGraph b-quark amplitudes not yet "
                 //"implemented."
          ! The same for -5010000231000023-5, 5-5100002310000230
          ! 50100002310000235, 0-510000231000023-5, 05100002310000235
          stop
#endif
         goto 20
      elseif(str.eq."-5010000231000023-5") then
         !call srealmtrx_018(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."5-5100002310000230") then
         !call srealmtrx_019(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."50100002310000235") then
         !call srealmtrx_020(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-110000231000023-1") then
         !call srealmtrx_021(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."01100002310000231") then
         !call srealmtrx_022(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-210000231000023-2") then
         !call srealmtrx_023(p,wgt)
         call SMATRIX_GUX_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."02100002310000232") then
         !call srealmtrx_024(p,wgt)
         call SMATRIX_GU_NINJU(p,wgt)
         goto 20
      elseif(str.eq."0-410000231000023-4") then
         !call srealmtrx_025(p,wgt)
         call SMATRIX_GUX_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."04100002310000234") then
         !call srealmtrx_026(p,wgt)
         call SMATRIX_GU_NINJU(p,wgt)
         goto 20
      elseif(str.eq."0-310000231000023-3") then
         !call srealmtrx_027(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."03100002310000233") then
         !call srealmtrx_028(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-510000231000023-5") then
         !call srealmtrx_029(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."05100002310000235") then
         !call srealmtrx_030(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
         
      elseif(str.eq."-11100002210000250") then
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DXD_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-1010000221000025-1") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."1-1100002210000250") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."10100002210000251") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt)
         goto 20
      elseif(str.eq."-22100002210000250") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UXU_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-2010000221000025-2") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UXG_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."2-2100002210000250") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_UUX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."20100002210000252") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_UG_NINJU(p,wgt)
         goto 20
      elseif(str.eq."-44100002210000250") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_UXU_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-4010000221000025-4") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_UXG_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."4-4100002210000250") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_UUX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."40100002210000254") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_UG_NINJU(p,wgt)
         goto 20
      elseif(str.eq."-33100002210000250") then
         !call srealmtrx_013(p,wgt)
         call SMATRIX_DXD_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-3010000221000025-3") then
         !call srealmtrx_014(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."3-3100002210000250") then
         !call srealmtrx_015(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."30100002210000253") then
         !call srealmtrx_016(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt) 
         goto 20
      elseif(str.eq."-55100002210000250") then
         !call srealmtrx_017(p,wgt)
#ifdef NEGLECTBMASS
          call SMATRIX_DXD_NINJG(p,wgt)
#else
          print*,"ERROR: MadGraph b-quark amplitudes not yet "
                 //"implemented."
          ! The same for -5010000221000025-5, 5-5100002210000250
          ! 50100002210000255, 0-510000221000025-5, 05100002210000255
          stop
#endif
         goto 20
      elseif(str.eq."-5010000221000025-5") then
         !call srealmtrx_018(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."5-5100002210000250") then
         !call srealmtrx_019(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."50100002210000255") then
         !call srealmtrx_020(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-110000221000025-1") then
         !call srealmtrx_021(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."01100002210000251") then
         !call srealmtrx_022(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-210000221000025-2") then
         !call srealmtrx_023(p,wgt)
         call SMATRIX_GUX_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."02100002210000252") then
         !call srealmtrx_024(p,wgt)
         call SMATRIX_GU_NINJU(p,wgt)
         goto 20
      elseif(str.eq."0-410000221000025-4") then
         !call srealmtrx_025(p,wgt)
         call SMATRIX_GUX_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."04100002210000254") then
         !call srealmtrx_026(p,wgt)
         call SMATRIX_GU_NINJU(p,wgt)
         goto 20
      elseif(str.eq."0-310000221000025-3") then
         !call srealmtrx_027(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."03100002210000253") then
         !call srealmtrx_028(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-510000221000025-5") then
         !call srealmtrx_029(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."05100002210000255") then
         !call srealmtrx_030(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
         
      elseif(str.eq."-11100002210000350") then
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DXD_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-1010000221000035-1") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."1-1100002210000350") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."10100002210000351") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt)
         goto 20
      elseif(str.eq."-22100002210000350") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UXU_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-2010000221000035-2") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UXG_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."2-2100002210000350") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_UUX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."20100002210000352") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_UG_NINJU(p,wgt)
         goto 20
      elseif(str.eq."-44100002210000350") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_UXU_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-4010000221000035-4") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_UXG_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."4-4100002210000350") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_UUX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."40100002210000354") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_UG_NINJU(p,wgt)
         goto 20
      elseif(str.eq."-33100002210000350") then
         !call srealmtrx_013(p,wgt)
         call SMATRIX_DXD_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-3010000221000035-3") then
         !call srealmtrx_014(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."3-3100002210000350") then
         !call srealmtrx_015(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."30100002210000353") then
         !call srealmtrx_016(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt) 
         goto 20
      elseif(str.eq."-55100002210000350") then
         !call srealmtrx_017(p,wgt)
#ifdef NEGLECTBMASS
          call SMATRIX_DXD_NINJG(p,wgt)
#else
          print*,"ERROR: MadGraph b-quark amplitudes not yet "
                 //"implemented."
          ! The same for -5010000221000035-5, 5-5100002210000350
          ! 50100002210000355, 0-510000221000035-5, 05100002210000355
          stop
#endif
         goto 20
      elseif(str.eq."-5010000221000035-5") then
         !call srealmtrx_018(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."5-5100002210000350") then
         !call srealmtrx_019(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."50100002210000355") then
         !call srealmtrx_020(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-110000221000035-1") then
         !call srealmtrx_021(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."01100002210000351") then
         !call srealmtrx_022(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-210000221000035-2") then
         !call srealmtrx_023(p,wgt)
         call SMATRIX_GUX_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."02100002210000352") then
         !call srealmtrx_024(p,wgt)
         call SMATRIX_GU_NINJU(p,wgt)
         goto 20
      elseif(str.eq."0-410000221000035-4") then
         !call srealmtrx_025(p,wgt)
         call SMATRIX_GUX_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."04100002210000354") then
         !call srealmtrx_026(p,wgt)
         call SMATRIX_GU_NINJU(p,wgt)
         goto 20
      elseif(str.eq."0-310000221000035-3") then
         !call srealmtrx_027(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."03100002210000353") then
         !call srealmtrx_028(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-510000221000035-5") then
         !call srealmtrx_029(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."05100002210000355") then
         !call srealmtrx_030(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
         
      elseif(str.eq."-11100002310000250") then
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DXD_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-1010000231000025-1") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."1-1100002310000250") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."10100002310000251") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt)
         goto 20
      elseif(str.eq."-22100002310000250") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UXU_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-2010000231000025-2") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UXG_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."2-2100002310000250") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_UUX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."20100002310000252") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_UG_NINJU(p,wgt)
         goto 20
      elseif(str.eq."-44100002310000250") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_UXU_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-4010000231000025-4") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_UXG_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."4-4100002310000250") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_UUX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."40100002310000254") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_UG_NINJU(p,wgt)
         goto 20
      elseif(str.eq."-33100002310000250") then
         !call srealmtrx_013(p,wgt)
         call SMATRIX_DXD_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-3010000231000025-3") then
         !call srealmtrx_014(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."3-3100002310000250") then
         !call srealmtrx_015(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."30100002310000253") then
         !call srealmtrx_016(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt) 
         goto 20
      elseif(str.eq."-55100002310000250") then
         !call srealmtrx_017(p,wgt)
#ifdef NEGLECTBMASS
          call SMATRIX_DXD_NINJG(p,wgt)
#else
          print*,"ERROR: MadGraph b-quark amplitudes not yet "
                 //"implemented."
          ! The same for -5010000231000025-5, 5-5100002310000250
          ! 50100002310000255, 0-510000231000025-5, 05100002310000255
          stop
#endif
         goto 20
      elseif(str.eq."-5010000231000025-5") then
         !call srealmtrx_018(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."5-5100002310000250") then
         !call srealmtrx_019(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."50100002310000255") then
         !call srealmtrx_020(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-110000231000025-1") then
         !call srealmtrx_021(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."01100002310000251") then
         !call srealmtrx_022(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-210000231000025-2") then
         !call srealmtrx_023(p,wgt)
         call SMATRIX_GUX_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."02100002310000252") then
         !call srealmtrx_024(p,wgt)
         call SMATRIX_GU_NINJU(p,wgt)
         goto 20
      elseif(str.eq."0-410000231000025-4") then
         !call srealmtrx_025(p,wgt)
         call SMATRIX_GUX_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."04100002310000254") then
         !call srealmtrx_026(p,wgt)
         call SMATRIX_GU_NINJU(p,wgt)
         goto 20
      elseif(str.eq."0-310000231000025-3") then
         !call srealmtrx_027(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."03100002310000253") then
         !call srealmtrx_028(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-510000231000025-5") then
         !call srealmtrx_029(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."05100002310000255") then
         !call srealmtrx_030(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
         
      elseif(str.eq."-11100002310000350") then
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DXD_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-1010000231000035-1") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."1-1100002310000350") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."10100002310000351") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt)
         goto 20
      elseif(str.eq."-22100002310000350") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UXU_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-2010000231000035-2") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UXG_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."2-2100002310000350") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_UUX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."20100002310000352") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_UG_NINJU(p,wgt)
         goto 20
      elseif(str.eq."-44100002310000350") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_UXU_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-4010000231000035-4") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_UXG_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."4-4100002310000350") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_UUX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."40100002310000354") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_UG_NINJU(p,wgt)
         goto 20
      elseif(str.eq."-33100002310000350") then
         !call srealmtrx_013(p,wgt)
         call SMATRIX_DXD_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-3010000231000035-3") then
         !call srealmtrx_014(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."3-3100002310000350") then
         !call srealmtrx_015(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."30100002310000353") then
         !call srealmtrx_016(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt) 
         goto 20
      elseif(str.eq."-55100002310000350") then
         !call srealmtrx_017(p,wgt)
#ifdef NEGLECTBMASS
          call SMATRIX_DXD_NINJG(p,wgt)
#else
          print*,"ERROR: MadGraph b-quark amplitudes not yet "
                 //"implemented."
          ! The same for -5010000231000035-5, 5-5100002310000350
          ! 50100002310000355, 0-510000231000035-5, 05100002310000355
          stop
#endif
         goto 20
      elseif(str.eq."-5010000231000035-5") then
         !call srealmtrx_018(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."5-5100002310000350") then
         !call srealmtrx_019(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."50100002310000355") then
         !call srealmtrx_020(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-110000231000035-1") then
         !call srealmtrx_021(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."01100002310000351") then
         !call srealmtrx_022(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-210000231000035-2") then
         !call srealmtrx_023(p,wgt)
         call SMATRIX_GUX_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."02100002310000352") then
         !call srealmtrx_024(p,wgt)
         call SMATRIX_GU_NINJU(p,wgt)
         goto 20
      elseif(str.eq."0-410000231000035-4") then
         !call srealmtrx_025(p,wgt)
         call SMATRIX_GUX_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."04100002310000354") then
         !call srealmtrx_026(p,wgt)
         call SMATRIX_GU_NINJU(p,wgt)
         goto 20
      elseif(str.eq."0-310000231000035-3") then
         !call srealmtrx_027(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."03100002310000353") then
         !call srealmtrx_028(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-510000231000035-5") then
         !call srealmtrx_029(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."05100002310000355") then
         !call srealmtrx_030(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
         
      elseif(str.eq."-11100002510000250") then
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DXD_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-1010000251000025-1") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."1-1100002510000250") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."10100002510000251") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt)
         goto 20
      elseif(str.eq."-22100002510000250") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UXU_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-2010000251000025-2") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UXG_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."2-2100002510000250") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_UUX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."20100002510000252") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_UG_NINJU(p,wgt)
         goto 20
      elseif(str.eq."-44100002510000250") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_UXU_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-4010000251000025-4") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_UXG_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."4-4100002510000250") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_UUX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."40100002510000254") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_UG_NINJU(p,wgt)
         goto 20
      elseif(str.eq."-33100002510000250") then
         !call srealmtrx_013(p,wgt)
         call SMATRIX_DXD_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-3010000251000025-3") then
         !call srealmtrx_014(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."3-3100002510000250") then
         !call srealmtrx_015(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."30100002510000253") then
         !call srealmtrx_016(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt) 
         goto 20
      elseif(str.eq."-55100002510000250") then
         !call srealmtrx_017(p,wgt)
#ifdef NEGLECTBMASS
          call SMATRIX_DXD_NINJG(p,wgt)
#else
          print*,"ERROR: MadGraph b-quark amplitudes not yet "
                 //"implemented."
          ! The same for -5010000251000025-5, 5-5100002510000250
          ! 50100002510000255, 0-510000251000025-5, 05100002510000255
          stop
#endif
         goto 20
      elseif(str.eq."-5010000251000025-5") then
         !call srealmtrx_018(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."5-5100002510000250") then
         !call srealmtrx_019(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."50100002510000255") then
         !call srealmtrx_020(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-110000251000025-1") then
         !call srealmtrx_021(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."01100002510000251") then
         !call srealmtrx_022(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-210000251000025-2") then
         !call srealmtrx_023(p,wgt)
         call SMATRIX_GUX_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."02100002510000252") then
         !call srealmtrx_024(p,wgt)
         call SMATRIX_GU_NINJU(p,wgt)
         goto 20
      elseif(str.eq."0-410000251000025-4") then
         !call srealmtrx_025(p,wgt)
         call SMATRIX_GUX_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."04100002510000254") then
         !call srealmtrx_026(p,wgt)
         call SMATRIX_GU_NINJU(p,wgt)
         goto 20
      elseif(str.eq."0-310000251000025-3") then
         !call srealmtrx_027(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."03100002510000253") then
         !call srealmtrx_028(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-510000251000025-5") then
         !call srealmtrx_029(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."05100002510000255") then
         !call srealmtrx_030(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
         
      elseif(str.eq."-11100002510000350") then
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DXD_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-1010000251000035-1") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."1-1100002510000350") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."10100002510000351") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt)
         goto 20
      elseif(str.eq."-22100002510000350") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UXU_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-2010000251000035-2") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UXG_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."2-2100002510000350") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_UUX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."20100002510000352") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_UG_NINJU(p,wgt)
         goto 20
      elseif(str.eq."-44100002510000350") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_UXU_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-4010000251000035-4") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_UXG_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."4-4100002510000350") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_UUX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."40100002510000354") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_UG_NINJU(p,wgt)
         goto 20
      elseif(str.eq."-33100002510000350") then
         !call srealmtrx_013(p,wgt)
         call SMATRIX_DXD_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-3010000251000035-3") then
         !call srealmtrx_014(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."3-3100002510000350") then
         !call srealmtrx_015(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."30100002510000353") then
         !call srealmtrx_016(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt) 
         goto 20
      elseif(str.eq."-55100002510000350") then
         !call srealmtrx_017(p,wgt)
#ifdef NEGLECTBMASS
          call SMATRIX_DXD_NINJG(p,wgt)
#else
          print*,"ERROR: MadGraph b-quark amplitudes not yet "
                 //"implemented."
          ! The same for -5010000251000035-5, 5-5100002510000350
          ! 50100002510000355, 0-510000251000035-5, 05100002510000355
          stop
#endif
         goto 20
      elseif(str.eq."-5010000251000035-5") then
         !call srealmtrx_018(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."5-5100002510000350") then
         !call srealmtrx_019(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."50100002510000355") then
         !call srealmtrx_020(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-110000251000035-1") then
         !call srealmtrx_021(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."01100002510000351") then
         !call srealmtrx_022(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-210000251000035-2") then
         !call srealmtrx_023(p,wgt)
         call SMATRIX_GUX_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."02100002510000352") then
         !call srealmtrx_024(p,wgt)
         call SMATRIX_GU_NINJU(p,wgt)
         goto 20
      elseif(str.eq."0-410000251000035-4") then
         !call srealmtrx_025(p,wgt)
         call SMATRIX_GUX_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."04100002510000354") then
         !call srealmtrx_026(p,wgt)
         call SMATRIX_GU_NINJU(p,wgt)
         goto 20
      elseif(str.eq."0-310000251000035-3") then
         !call srealmtrx_027(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."03100002510000353") then
         !call srealmtrx_028(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-510000251000035-5") then
         !call srealmtrx_029(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."05100002510000355") then
         !call srealmtrx_030(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20        
         
      elseif(str.eq."-11100003510000350") then
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DXD_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-1010000351000035-1") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."1-1100003510000350") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."10100003510000351") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt)
         goto 20
      elseif(str.eq."-22100003510000350") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UXU_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-2010000351000035-2") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UXG_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."2-2100003510000350") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_UUX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."20100003510000352") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_UG_NINJU(p,wgt)
         goto 20
      elseif(str.eq."-44100003510000350") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_UXU_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-4010000351000035-4") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_UXG_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."4-4100003510000350") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_UUX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."40100003510000354") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_UG_NINJU(p,wgt)
         goto 20
      elseif(str.eq."-33100003510000350") then
         !call srealmtrx_013(p,wgt)
         call SMATRIX_DXD_NINJG(p,wgt)
         goto 20
      elseif(str.eq."-3010000351000035-3") then
         !call srealmtrx_014(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."3-3100003510000350") then
         !call srealmtrx_015(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."30100003510000353") then
         !call srealmtrx_016(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt) 
         goto 20
      elseif(str.eq."-55100003510000350") then
         !call srealmtrx_017(p,wgt)
#ifdef NEGLECTBMASS
          call SMATRIX_DXD_NINJG(p,wgt)
#else
          print*,"ERROR: MadGraph b-quark amplitudes not yet "
                 //"implemented."
          ! The same for -5010000351000035-5, 5-5100003510000350
          ! 50100003510000355, 0-510000351000035-5, 05100003510000355
          stop
#endif
         goto 20
      elseif(str.eq."-5010000351000035-5") then
         !call srealmtrx_018(p,wgt)
         call SMATRIX_DXG_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."5-5100003510000350") then
         !call srealmtrx_019(p,wgt)
         call SMATRIX_DDX_NINJG(p,wgt)
         goto 20
      elseif(str.eq."50100003510000355") then
         !call srealmtrx_020(p,wgt)
         call SMATRIX_DG_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-110000351000035-1") then
         !call srealmtrx_021(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."01100003510000351") then
         !call srealmtrx_022(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-210000351000035-2") then
         !call srealmtrx_023(p,wgt)
         call SMATRIX_GUX_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."02100003510000352") then
         !call srealmtrx_024(p,wgt)
         call SMATRIX_GU_NINJU(p,wgt)
         goto 20
      elseif(str.eq."0-410000351000035-4") then
         !call srealmtrx_025(p,wgt)
         call SMATRIX_GUX_NINJUX(p,wgt)
         goto 20
      elseif(str.eq."04100003510000354") then
         !call srealmtrx_026(p,wgt)
         call SMATRIX_GU_NINJU(p,wgt)
         goto 20
      elseif(str.eq."0-310000351000035-3") then
         !call srealmtrx_027(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."03100003510000353") then
         !call srealmtrx_028(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
      elseif(str.eq."0-510000351000035-5") then
         !call srealmtrx_029(p,wgt)
         call SMATRIX_GDX_NINJDX(p,wgt)
         goto 20
      elseif(str.eq."05100003510000355") then
         !call srealmtrx_030(p,wgt)
         call SMATRIX_GD_NINJD(p,wgt)
         goto 20
         
      endif
      
 20   continue
      if(wgt .eq. 0D0) then
          print*, "str = ", str
          print*,"WARNING: maybe error in real amplitudes: ", wgt
          stop
      endif

      return
      end
      
      
      subroutine real_color(legs,color)
      implicit none
      include "nexternal.inc"
      include "maxamps_ninj.inc"
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
      
      equivalence(amp2001,amp2013,amp2017)
      equivalence(jamp2001,jamp2013,jamp2017)
      
      equivalence(amp2002,amp2014,amp2018)
      equivalence(jamp2002,jamp2014,jamp2018)
      
      equivalence(amp2003,amp2015,amp2019)
      equivalence(jamp2003,jamp2015,jamp2019)
      
      equivalence(amp2004,amp2016,amp2020)
      equivalence(jamp2004,jamp2016,jamp2020)
      
      equivalence(amp2005,amp2009)
      equivalence(jamp2005,jamp2009)
      
      equivalence(amp2006,amp2010)
      equivalence(jamp2006,jamp2010)
      
      equivalence(amp2007,amp2011)
      equivalence(jamp2007,jamp2011)
      
      equivalence(amp2008,amp2012)
      equivalence(jamp2008,jamp2012)
      
      equivalence(amp2021,amp2027,amp2029)
      equivalence(jamp2021,jamp2027,jamp2029)
      
      equivalence(amp2022,amp2028,amp2030)
      equivalence(jamp2022,jamp2028,jamp2030)
      
      equivalence(amp2023,amp2025)
      equivalence(jamp2023,jamp2025)
      
      equivalence(amp2024,amp2026)
      equivalence(jamp2024,jamp2026)
      
      double precision jamp2cum(0:maxamps)
      integer ICOLUP(2,nexternal,maxamps)
      integer color(2,nexternal),color1(2,nexternal)
      double precision random,xtarget
      external random
      integer legs(nexternal),lstr,i,j
      character*140 str
      integer iflow,ifl
      
      call convert_to_string(nexternal,legs,str,lstr)

      if (str.eq."-11100002210000220") then
         include "leshouches_R_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif (str.eq."-1010000221000022-1") then
         include "leshouches_R_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif (str.eq."1-1100002210000220") then
         include "leshouches_R_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif (str.eq."10100002210000221") then
         include "leshouches_R_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif (str.eq."-22100002210000220") then
         include "leshouches_R_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif (str.eq."-2010000221000022-2") then
         include "leshouches_R_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif (str.eq."2-2100002210000220") then
         include "leshouches_R_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif (str.eq."20100002210000222") then
         include "leshouches_R_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif (str.eq."-44100002210000220") then
         include "leshouches_R_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif (str.eq."-4010000221000022-4") then
         include "leshouches_R_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif (str.eq."4-4100002210000220") then
         include "leshouches_R_011.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif (str.eq."40100002210000224") then
         include "leshouches_R_012.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
      elseif (str.eq."-33100002210000220") then
         include "leshouches_R_013.inc"
         iflow=nint(jamp2013(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2013(i)
         enddo
         goto 20
      elseif (str.eq."-3010000221000022-3") then
         include "leshouches_R_014.inc"
         iflow=nint(jamp2014(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2014(i)
         enddo
         goto 20
      elseif (str.eq."3-3100002210000220") then
         include "leshouches_R_015.inc"
         iflow=nint(jamp2015(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2015(i)
         enddo
         goto 20
      elseif (str.eq."30100002210000223") then
         include "leshouches_R_016.inc"
         iflow=nint(jamp2016(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2016(i)
         enddo
         goto 20
      elseif (str.eq."-55100002210000220") then
         include "leshouches_R_017.inc"
         iflow=nint(jamp2017(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2017(i)
         enddo
         goto 20
      elseif (str.eq."-5010000221000022-5") then
         include "leshouches_R_018.inc"
         iflow=nint(jamp2018(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2018(i)
         enddo
         goto 20
      elseif (str.eq."5-5100002210000220") then
         include "leshouches_R_019.inc"
         iflow=nint(jamp2019(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2019(i)
         enddo
         goto 20
      elseif (str.eq."50100002210000225") then
         include "leshouches_R_020.inc"
         iflow=nint(jamp2020(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2020(i)
         enddo
         goto 20
      elseif (str.eq."0-110000221000022-1") then
         include "leshouches_R_021.inc"
         iflow=nint(jamp2021(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2021(i)
         enddo
         goto 20
      elseif (str.eq."01100002210000221") then
         include "leshouches_R_022.inc"
         iflow=nint(jamp2022(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2022(i)
         enddo
         goto 20
      elseif (str.eq."0-210000221000022-2") then
         include "leshouches_R_023.inc"
         iflow=nint(jamp2023(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2023(i)
         enddo
         goto 20
      elseif (str.eq."02100002210000222") then
         include "leshouches_R_024.inc"
         iflow=nint(jamp2024(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2024(i)
         enddo
         goto 20
      elseif (str.eq."0-410000221000022-4") then
         include "leshouches_R_025.inc"
         iflow=nint(jamp2025(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2025(i)
         enddo
         goto 20
      elseif (str.eq."04100002210000224") then
         include "leshouches_R_026.inc"
         iflow=nint(jamp2026(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2026(i)
         enddo
         goto 20
      elseif (str.eq."0-310000221000022-3") then
         include "leshouches_R_027.inc"
         iflow=nint(jamp2027(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2027(i)
         enddo
         goto 20
      elseif (str.eq."03100002210000223") then
         include "leshouches_R_028.inc"
         iflow=nint(jamp2028(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2028(i)
         enddo
         goto 20
      elseif (str.eq."0-510000221000022-5") then
         include "leshouches_R_029.inc"
         iflow=nint(jamp2029(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2029(i)
         enddo
         goto 20
      elseif (str.eq."05100002210000225") then
         include "leshouches_R_030.inc"
         iflow=nint(jamp2030(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2030(i)
         enddo
         goto 20
         
      elseif(str.eq."-11100002210000230") then
         include "leshouches_R_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-1010000221000023-1") then
         include "leshouches_R_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."1-1100002210000230") then
         include "leshouches_R_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."10100002210000231") then
         include "leshouches_R_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-22100002210000230") then
         include "leshouches_R_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."-2010000221000023-2") then
         include "leshouches_R_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."2-2100002210000230") then
         include "leshouches_R_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."20100002210000232") then
         include "leshouches_R_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."-44100002210000230") then
         include "leshouches_R_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."-4010000221000023-4") then
         include "leshouches_R_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."4-4100002210000230") then
         include "leshouches_R_011.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."40100002210000234") then
         include "leshouches_R_012.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
      elseif(str.eq."-33100002210000230") then
         include "leshouches_R_013.inc"
         iflow=nint(jamp2013(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2013(i)
         enddo
         goto 20
      elseif(str.eq."-3010000221000023-3") then
         include "leshouches_R_014.inc"
         iflow=nint(jamp2014(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2014(i)
         enddo
         goto 20
      elseif(str.eq."3-3100002210000230") then
         include "leshouches_R_015.inc"
         iflow=nint(jamp2015(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2015(i)
         enddo
         goto 20
      elseif(str.eq."30100002210000233") then
         include "leshouches_R_016.inc"
         iflow=nint(jamp2016(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2016(i)
         enddo
         goto 20
      elseif(str.eq."-55100002210000230") then
         include "leshouches_R_017.inc"
         iflow=nint(jamp2017(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2017(i)
         enddo
         goto 20
      elseif(str.eq."-5010000221000023-5") then
         include "leshouches_R_018.inc"
         iflow=nint(jamp2018(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2018(i)
         enddo
         goto 20
      elseif(str.eq."5-5100002210000230") then
         include "leshouches_R_019.inc"
         iflow=nint(jamp2019(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2019(i)
         enddo
         goto 20
      elseif(str.eq."50100002210000235") then
         include "leshouches_R_020.inc"
         iflow=nint(jamp2020(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2020(i)
         enddo
         goto 20
      elseif(str.eq."0-110000221000023-1") then
         include "leshouches_R_021.inc"
         iflow=nint(jamp2021(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2021(i)
         enddo
         goto 20
      elseif(str.eq."01100002210000231") then
         include "leshouches_R_022.inc"
         iflow=nint(jamp2022(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2022(i)
         enddo
         goto 20
      elseif(str.eq."0-210000221000023-2") then
         include "leshouches_R_023.inc"
         iflow=nint(jamp2023(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2023(i)
         enddo
         goto 20
      elseif(str.eq."02100002210000232") then
         include "leshouches_R_024.inc"
         iflow=nint(jamp2024(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2024(i)
         enddo
         goto 20
      elseif(str.eq."0-410000221000023-4") then
         include "leshouches_R_025.inc"
         iflow=nint(jamp2025(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2025(i)
         enddo
         goto 20
      elseif(str.eq."04100002210000234") then
         include "leshouches_R_026.inc"
         iflow=nint(jamp2026(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2026(i)
         enddo
         goto 20
      elseif(str.eq."0-310000221000023-3") then
         include "leshouches_R_027.inc"
         iflow=nint(jamp2027(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2027(i)
         enddo
         goto 20
      elseif(str.eq."03100002210000233") then
         include "leshouches_R_028.inc"
         iflow=nint(jamp2028(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2028(i)
         enddo
         goto 20
      elseif(str.eq."0-510000221000023-5") then
         include "leshouches_R_029.inc"
         iflow=nint(jamp2029(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2029(i)
         enddo
         goto 20
      elseif(str.eq."05100002210000235") then
         include "leshouches_R_030.inc"
         iflow=nint(jamp2030(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2030(i)
         enddo
         goto 20
         
      elseif(str.eq."-11100002310000230") then
         include "leshouches_R_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-1010000231000023-1") then
         include "leshouches_R_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."1-1100002310000230") then
         include "leshouches_R_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."10100002310000231") then
         include "leshouches_R_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-22100002310000230") then
         include "leshouches_R_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."-2010000231000023-2") then
         include "leshouches_R_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."2-2100002310000230") then
         include "leshouches_R_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."20100002310000232") then
         include "leshouches_R_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."-44100002310000230") then
         include "leshouches_R_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."-4010000231000023-4") then
         include "leshouches_R_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."4-4100002310000230") then
         include "leshouches_R_011.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."40100002310000234") then
         include "leshouches_R_012.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
      elseif(str.eq."-33100002310000230") then
         include "leshouches_R_013.inc"
         iflow=nint(jamp2013(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2013(i)
         enddo
         goto 20
      elseif(str.eq."-3010000231000023-3") then
         include "leshouches_R_014.inc"
         iflow=nint(jamp2014(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2014(i)
         enddo
         goto 20
      elseif(str.eq."3-3100002310000230") then
         include "leshouches_R_015.inc"
         iflow=nint(jamp2015(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2015(i)
         enddo
         goto 20
      elseif(str.eq."30100002310000233") then
         include "leshouches_R_016.inc"
         iflow=nint(jamp2016(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2016(i)
         enddo
         goto 20
      elseif(str.eq."-55100002310000230") then
         include "leshouches_R_017.inc"
         iflow=nint(jamp2017(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2017(i)
         enddo
         goto 20
      elseif(str.eq."-5010000231000023-5") then
         include "leshouches_R_018.inc"
         iflow=nint(jamp2018(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2018(i)
         enddo
         goto 20
      elseif(str.eq."5-5100002310000230") then
         include "leshouches_R_019.inc"
         iflow=nint(jamp2019(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2019(i)
         enddo
         goto 20
      elseif(str.eq."50100002310000235") then
         include "leshouches_R_020.inc"
         iflow=nint(jamp2020(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2020(i)
         enddo
         goto 20
      elseif(str.eq."0-110000231000023-1") then
         include "leshouches_R_021.inc"
         iflow=nint(jamp2021(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2021(i)
         enddo
         goto 20
      elseif(str.eq."01100002310000231") then
         include "leshouches_R_022.inc"
         iflow=nint(jamp2022(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2022(i)
         enddo
         goto 20
      elseif(str.eq."0-210000231000023-2") then
         include "leshouches_R_023.inc"
         iflow=nint(jamp2023(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2023(i)
         enddo
         goto 20
      elseif(str.eq."02100002310000232") then
         include "leshouches_R_024.inc"
         iflow=nint(jamp2024(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2024(i)
         enddo
         goto 20
      elseif(str.eq."0-410000231000023-4") then
         include "leshouches_R_025.inc"
         iflow=nint(jamp2025(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2025(i)
         enddo
         goto 20
      elseif(str.eq."04100002310000234") then
         include "leshouches_R_026.inc"
         iflow=nint(jamp2026(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2026(i)
         enddo
         goto 20
      elseif(str.eq."0-310000231000023-3") then
         include "leshouches_R_027.inc"
         iflow=nint(jamp2027(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2027(i)
         enddo
         goto 20
      elseif(str.eq."03100002310000233") then
         include "leshouches_R_028.inc"
         iflow=nint(jamp2028(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2028(i)
         enddo
         goto 20
      elseif(str.eq."0-510000231000023-5") then
         include "leshouches_R_029.inc"
         iflow=nint(jamp2029(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2029(i)
         enddo
         goto 20
      elseif(str.eq."05100002310000235") then
         include "leshouches_R_030.inc"
         iflow=nint(jamp2030(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2030(i)
         enddo
         goto 20
         
      elseif(str.eq."-11100002210000250") then
         include "leshouches_R_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-1010000221000025-1") then
         include "leshouches_R_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."1-1100002210000250") then
         include "leshouches_R_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."10100002210000251") then
         include "leshouches_R_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-22100002210000250") then
         include "leshouches_R_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."-2010000221000025-2") then
         include "leshouches_R_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."2-2100002210000250") then
         include "leshouches_R_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."20100002210000252") then
         include "leshouches_R_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."-44100002210000250") then
         include "leshouches_R_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."-4010000221000025-4") then
         include "leshouches_R_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."4-4100002210000250") then
         include "leshouches_R_011.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."40100002210000254") then
         include "leshouches_R_012.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
      elseif(str.eq."-33100002210000250") then
         include "leshouches_R_013.inc"
         iflow=nint(jamp2013(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2013(i)
         enddo
         goto 20
      elseif(str.eq."-3010000221000025-3") then
         include "leshouches_R_014.inc"
         iflow=nint(jamp2014(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2014(i)
         enddo
         goto 20
      elseif(str.eq."3-3100002210000250") then
         include "leshouches_R_015.inc"
         iflow=nint(jamp2015(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2015(i)
         enddo
         goto 20
      elseif(str.eq."30100002210000253") then
         include "leshouches_R_016.inc"
         iflow=nint(jamp2016(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2016(i)
         enddo
         goto 20
      elseif(str.eq."-55100002210000250") then
         include "leshouches_R_017.inc"
         iflow=nint(jamp2017(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2017(i)
         enddo
         goto 20
      elseif(str.eq."-5010000221000025-5") then
         include "leshouches_R_018.inc"
         iflow=nint(jamp2018(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2018(i)
         enddo
         goto 20
      elseif(str.eq."5-5100002210000250") then
         include "leshouches_R_019.inc"
         iflow=nint(jamp2019(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2019(i)
         enddo
         goto 20
      elseif(str.eq."50100002210000255") then
         include "leshouches_R_020.inc"
         iflow=nint(jamp2020(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2020(i)
         enddo
         goto 20
      elseif(str.eq."0-110000221000025-1") then
         include "leshouches_R_021.inc"
         iflow=nint(jamp2021(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2021(i)
         enddo
         goto 20
      elseif(str.eq."01100002210000251") then
         include "leshouches_R_022.inc"
         iflow=nint(jamp2022(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2022(i)
         enddo
         goto 20
      elseif(str.eq."0-210000221000025-2") then
         include "leshouches_R_023.inc"
         iflow=nint(jamp2023(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2023(i)
         enddo
         goto 20
      elseif(str.eq."02100002210000252") then
         include "leshouches_R_024.inc"
         iflow=nint(jamp2024(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2024(i)
         enddo
         goto 20
      elseif(str.eq."0-410000221000025-4") then
         include "leshouches_R_025.inc"
         iflow=nint(jamp2025(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2025(i)
         enddo
         goto 20
      elseif(str.eq."04100002210000254") then
         include "leshouches_R_026.inc"
         iflow=nint(jamp2026(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2026(i)
         enddo
         goto 20
      elseif(str.eq."0-310000221000025-3") then
         include "leshouches_R_027.inc"
         iflow=nint(jamp2027(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2027(i)
         enddo
         goto 20
      elseif(str.eq."03100002210000253") then
         include "leshouches_R_028.inc"
         iflow=nint(jamp2028(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2028(i)
         enddo
         goto 20
      elseif(str.eq."0-510000221000025-5") then
         include "leshouches_R_029.inc"
         iflow=nint(jamp2029(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2029(i)
         enddo
         goto 20
      elseif(str.eq."05100002210000255") then
         include "leshouches_R_030.inc"
         iflow=nint(jamp2030(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2030(i)
         enddo
         goto 20
         
      elseif(str.eq."-11100002210000350") then
         include "leshouches_R_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-1010000221000035-1") then
         include "leshouches_R_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."1-1100002210000350") then
         include "leshouches_R_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."10100002210000351") then
         include "leshouches_R_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-22100002210000350") then
         include "leshouches_R_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."-2010000221000035-2") then
         include "leshouches_R_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."2-2100002210000350") then
         include "leshouches_R_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."20100002210000352") then
         include "leshouches_R_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."-44100002210000350") then
         include "leshouches_R_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."-4010000221000035-4") then
         include "leshouches_R_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."4-4100002210000350") then
         include "leshouches_R_011.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."40100002210000354") then
         include "leshouches_R_012.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
      elseif(str.eq."-33100002210000350") then
         include "leshouches_R_013.inc"
         iflow=nint(jamp2013(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2013(i)
         enddo
         goto 20
      elseif(str.eq."-3010000221000035-3") then
         include "leshouches_R_014.inc"
         iflow=nint(jamp2014(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2014(i)
         enddo
         goto 20
      elseif(str.eq."3-3100002210000350") then
         include "leshouches_R_015.inc"
         iflow=nint(jamp2015(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2015(i)
         enddo
         goto 20
      elseif(str.eq."30100002210000353") then
         include "leshouches_R_016.inc"
         iflow=nint(jamp2016(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2016(i)
         enddo
         goto 20
      elseif(str.eq."-55100002210000350") then
         include "leshouches_R_017.inc"
         iflow=nint(jamp2017(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2017(i)
         enddo
         goto 20
      elseif(str.eq."-5010000221000035-5") then
         include "leshouches_R_018.inc"
         iflow=nint(jamp2018(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2018(i)
         enddo
         goto 20
      elseif(str.eq."5-5100002210000350") then
         include "leshouches_R_019.inc"
         iflow=nint(jamp2019(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2019(i)
         enddo
         goto 20
      elseif(str.eq."50100002210000355") then
         include "leshouches_R_020.inc"
         iflow=nint(jamp2020(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2020(i)
         enddo
         goto 20
      elseif(str.eq."0-110000221000035-1") then
         include "leshouches_R_021.inc"
         iflow=nint(jamp2021(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2021(i)
         enddo
         goto 20
      elseif(str.eq."01100002210000351") then
         include "leshouches_R_022.inc"
         iflow=nint(jamp2022(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2022(i)
         enddo
         goto 20
      elseif(str.eq."0-210000221000035-2") then
         include "leshouches_R_023.inc"
         iflow=nint(jamp2023(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2023(i)
         enddo
         goto 20
      elseif(str.eq."02100002210000352") then
         include "leshouches_R_024.inc"
         iflow=nint(jamp2024(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2024(i)
         enddo
         goto 20
      elseif(str.eq."0-410000221000035-4") then
         include "leshouches_R_025.inc"
         iflow=nint(jamp2025(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2025(i)
         enddo
         goto 20
      elseif(str.eq."04100002210000354") then
         include "leshouches_R_026.inc"
         iflow=nint(jamp2026(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2026(i)
         enddo
         goto 20
      elseif(str.eq."0-310000221000035-3") then
         include "leshouches_R_027.inc"
         iflow=nint(jamp2027(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2027(i)
         enddo
         goto 20
      elseif(str.eq."03100002210000353") then
         include "leshouches_R_028.inc"
         iflow=nint(jamp2028(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2028(i)
         enddo
         goto 20
      elseif(str.eq."0-510000221000035-5") then
         include "leshouches_R_029.inc"
         iflow=nint(jamp2029(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2029(i)
         enddo
         goto 20
      elseif(str.eq."05100002210000355") then
         include "leshouches_R_030.inc"
         iflow=nint(jamp2030(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2030(i)
         enddo
         goto 20
         
      elseif(str.eq."-11100002310000250") then
         include "leshouches_R_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-1010000231000025-1") then
         include "leshouches_R_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."1-1100002310000250") then
         include "leshouches_R_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."10100002310000251") then
         include "leshouches_R_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-22100002310000250") then
         include "leshouches_R_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."-2010000231000025-2") then
         include "leshouches_R_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."2-2100002310000250") then
         include "leshouches_R_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."20100002310000252") then
         include "leshouches_R_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."-44100002310000250") then
         include "leshouches_R_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."-4010000231000025-4") then
         include "leshouches_R_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."4-4100002310000250") then
         include "leshouches_R_011.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."40100002310000254") then
         include "leshouches_R_012.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
      elseif(str.eq."-33100002310000250") then
         include "leshouches_R_013.inc"
         iflow=nint(jamp2013(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2013(i)
         enddo
         goto 20
      elseif(str.eq."-3010000231000025-3") then
         include "leshouches_R_014.inc"
         iflow=nint(jamp2014(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2014(i)
         enddo
         goto 20
      elseif(str.eq."3-3100002310000250") then
         include "leshouches_R_015.inc"
         iflow=nint(jamp2015(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2015(i)
         enddo
         goto 20
      elseif(str.eq."30100002310000253") then
         include "leshouches_R_016.inc"
         iflow=nint(jamp2016(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2016(i)
         enddo
         goto 20
      elseif(str.eq."-55100002310000250") then
         include "leshouches_R_017.inc"
         iflow=nint(jamp2017(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2017(i)
         enddo
         goto 20
      elseif(str.eq."-5010000231000025-5") then
         include "leshouches_R_018.inc"
         iflow=nint(jamp2018(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2018(i)
         enddo
         goto 20
      elseif(str.eq."5-5100002310000250") then
         include "leshouches_R_019.inc"
         iflow=nint(jamp2019(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2019(i)
         enddo
         goto 20
      elseif(str.eq."50100002310000255") then
         include "leshouches_R_020.inc"
         iflow=nint(jamp2020(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2020(i)
         enddo
         goto 20
      elseif(str.eq."0-110000231000025-1") then
         include "leshouches_R_021.inc"
         iflow=nint(jamp2021(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2021(i)
         enddo
         goto 20
      elseif(str.eq."01100002310000251") then
         include "leshouches_R_022.inc"
         iflow=nint(jamp2022(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2022(i)
         enddo
         goto 20
      elseif(str.eq."0-210000231000025-2") then
         include "leshouches_R_023.inc"
         iflow=nint(jamp2023(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2023(i)
         enddo
         goto 20
      elseif(str.eq."02100002310000252") then
         include "leshouches_R_024.inc"
         iflow=nint(jamp2024(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2024(i)
         enddo
         goto 20
      elseif(str.eq."0-410000231000025-4") then
         include "leshouches_R_025.inc"
         iflow=nint(jamp2025(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2025(i)
         enddo
         goto 20
      elseif(str.eq."04100002310000254") then
         include "leshouches_R_026.inc"
         iflow=nint(jamp2026(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2026(i)
         enddo
         goto 20
      elseif(str.eq."0-310000231000025-3") then
         include "leshouches_R_027.inc"
         iflow=nint(jamp2027(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2027(i)
         enddo
         goto 20
      elseif(str.eq."03100002310000253") then
         include "leshouches_R_028.inc"
         iflow=nint(jamp2028(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2028(i)
         enddo
         goto 20
      elseif(str.eq."0-510000231000025-5") then
         include "leshouches_R_029.inc"
         iflow=nint(jamp2029(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2029(i)
         enddo
         goto 20
      elseif(str.eq."05100002310000255") then
         include "leshouches_R_030.inc"
         iflow=nint(jamp2030(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2030(i)
         enddo
         goto 20
         
      elseif(str.eq."-11100002310000350") then
         include "leshouches_R_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-1010000231000035-1") then
         include "leshouches_R_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."1-1100002310000350") then
         include "leshouches_R_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."10100002310000351") then
         include "leshouches_R_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-22100002310000350") then
         include "leshouches_R_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."-2010000231000035-2") then
         include "leshouches_R_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."2-2100002310000350") then
         include "leshouches_R_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."20100002310000352") then
         include "leshouches_R_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."-44100002310000350") then
         include "leshouches_R_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."-4010000231000035-4") then
         include "leshouches_R_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."4-4100002310000350") then
         include "leshouches_R_011.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."40100002310000354") then
         include "leshouches_R_012.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
      elseif(str.eq."-33100002310000350") then
         include "leshouches_R_013.inc"
         iflow=nint(jamp2013(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2013(i)
         enddo
         goto 20
      elseif(str.eq."-3010000231000035-3") then
         include "leshouches_R_014.inc"
         iflow=nint(jamp2014(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2014(i)
         enddo
         goto 20
      elseif(str.eq."3-3100002310000350") then
         include "leshouches_R_015.inc"
         iflow=nint(jamp2015(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2015(i)
         enddo
         goto 20
      elseif(str.eq."30100002310000353") then
         include "leshouches_R_016.inc"
         iflow=nint(jamp2016(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2016(i)
         enddo
         goto 20
      elseif(str.eq."-55100002310000350") then
         include "leshouches_R_017.inc"
         iflow=nint(jamp2017(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2017(i)
         enddo
         goto 20
      elseif(str.eq."-5010000231000035-5") then
         include "leshouches_R_018.inc"
         iflow=nint(jamp2018(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2018(i)
         enddo
         goto 20
      elseif(str.eq."5-5100002310000350") then
         include "leshouches_R_019.inc"
         iflow=nint(jamp2019(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2019(i)
         enddo
         goto 20
      elseif(str.eq."50100002310000355") then
         include "leshouches_R_020.inc"
         iflow=nint(jamp2020(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2020(i)
         enddo
         goto 20
      elseif(str.eq."0-110000231000035-1") then
         include "leshouches_R_021.inc"
         iflow=nint(jamp2021(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2021(i)
         enddo
         goto 20
      elseif(str.eq."01100002310000351") then
         include "leshouches_R_022.inc"
         iflow=nint(jamp2022(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2022(i)
         enddo
         goto 20
      elseif(str.eq."0-210000231000035-2") then
         include "leshouches_R_023.inc"
         iflow=nint(jamp2023(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2023(i)
         enddo
         goto 20
      elseif(str.eq."02100002310000352") then
         include "leshouches_R_024.inc"
         iflow=nint(jamp2024(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2024(i)
         enddo
         goto 20
      elseif(str.eq."0-410000231000035-4") then
         include "leshouches_R_025.inc"
         iflow=nint(jamp2025(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2025(i)
         enddo
         goto 20
      elseif(str.eq."04100002310000354") then
         include "leshouches_R_026.inc"
         iflow=nint(jamp2026(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2026(i)
         enddo
         goto 20
      elseif(str.eq."0-310000231000035-3") then
         include "leshouches_R_027.inc"
         iflow=nint(jamp2027(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2027(i)
         enddo
         goto 20
      elseif(str.eq."03100002310000353") then
         include "leshouches_R_028.inc"
         iflow=nint(jamp2028(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2028(i)
         enddo
         goto 20
      elseif(str.eq."0-510000231000035-5") then
         include "leshouches_R_029.inc"
         iflow=nint(jamp2029(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2029(i)
         enddo
         goto 20
      elseif(str.eq."05100002310000355") then
         include "leshouches_R_030.inc"
         iflow=nint(jamp2030(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2030(i)
         enddo
         goto 20
         
      elseif(str.eq."-11100002510000250") then
         include "leshouches_R_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-1010000251000025-1") then
         include "leshouches_R_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."1-1100002510000250") then
         include "leshouches_R_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."10100002510000251") then
         include "leshouches_R_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-22100002510000250") then
         include "leshouches_R_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."-2010000251000025-2") then
         include "leshouches_R_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."2-2100002510000250") then
         include "leshouches_R_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."20100002510000252") then
         include "leshouches_R_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."-44100002510000250") then
         include "leshouches_R_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."-4010000251000025-4") then
         include "leshouches_R_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."4-4100002510000250") then
         include "leshouches_R_011.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."40100002510000254") then
         include "leshouches_R_012.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
      elseif(str.eq."-33100002510000250") then
         include "leshouches_R_013.inc"
         iflow=nint(jamp2013(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2013(i)
         enddo
         goto 20
      elseif(str.eq."-3010000251000025-3") then
         include "leshouches_R_014.inc"
         iflow=nint(jamp2014(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2014(i)
         enddo
         goto 20
      elseif(str.eq."3-3100002510000250") then
         include "leshouches_R_015.inc"
         iflow=nint(jamp2015(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2015(i)
         enddo
         goto 20
      elseif(str.eq."30100002510000253") then
         include "leshouches_R_016.inc"
         iflow=nint(jamp2016(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2016(i)
         enddo
         goto 20
      elseif(str.eq."-55100002510000250") then
         include "leshouches_R_017.inc"
         iflow=nint(jamp2017(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2017(i)
         enddo
         goto 20
      elseif(str.eq."-5010000251000025-5") then
         include "leshouches_R_018.inc"
         iflow=nint(jamp2018(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2018(i)
         enddo
         goto 20
      elseif(str.eq."5-5100002510000250") then
         include "leshouches_R_019.inc"
         iflow=nint(jamp2019(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2019(i)
         enddo
         goto 20
      elseif(str.eq."50100002510000255") then
         include "leshouches_R_020.inc"
         iflow=nint(jamp2020(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2020(i)
         enddo
         goto 20
      elseif(str.eq."0-110000251000025-1") then
         include "leshouches_R_021.inc"
         iflow=nint(jamp2021(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2021(i)
         enddo
         goto 20
      elseif(str.eq."01100002510000251") then
         include "leshouches_R_022.inc"
         iflow=nint(jamp2022(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2022(i)
         enddo
         goto 20
      elseif(str.eq."0-210000251000025-2") then
         include "leshouches_R_023.inc"
         iflow=nint(jamp2023(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2023(i)
         enddo
         goto 20
      elseif(str.eq."02100002510000252") then
         include "leshouches_R_024.inc"
         iflow=nint(jamp2024(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2024(i)
         enddo
         goto 20
      elseif(str.eq."0-410000251000025-4") then
         include "leshouches_R_025.inc"
         iflow=nint(jamp2025(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2025(i)
         enddo
         goto 20
      elseif(str.eq."04100002510000254") then
         include "leshouches_R_026.inc"
         iflow=nint(jamp2026(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2026(i)
         enddo
         goto 20
      elseif(str.eq."0-310000251000025-3") then
         include "leshouches_R_027.inc"
         iflow=nint(jamp2027(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2027(i)
         enddo
         goto 20
      elseif(str.eq."03100002510000253") then
         include "leshouches_R_028.inc"
         iflow=nint(jamp2028(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2028(i)
         enddo
         goto 20
      elseif(str.eq."0-510000251000025-5") then
         include "leshouches_R_029.inc"
         iflow=nint(jamp2029(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2029(i)
         enddo
         goto 20
      elseif(str.eq."05100002510000255") then
         include "leshouches_R_030.inc"
         iflow=nint(jamp2030(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2030(i)
         enddo
         goto 20
         
      elseif(str.eq."-11100002510000350") then
         include "leshouches_R_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-1010000251000035-1") then
         include "leshouches_R_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."1-1100002510000350") then
         include "leshouches_R_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."10100002510000351") then
         include "leshouches_R_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-22100002510000350") then
         include "leshouches_R_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."-2010000251000035-2") then
         include "leshouches_R_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."2-2100002510000350") then
         include "leshouches_R_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."20100002510000352") then
         include "leshouches_R_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."-44100002510000350") then
         include "leshouches_R_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."-4010000251000035-4") then
         include "leshouches_R_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."4-4100002510000350") then
         include "leshouches_R_011.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."40100002510000354") then
         include "leshouches_R_012.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
      elseif(str.eq."-33100002510000350") then
         include "leshouches_R_013.inc"
         iflow=nint(jamp2013(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2013(i)
         enddo
         goto 20
      elseif(str.eq."-3010000251000035-3") then
         include "leshouches_R_014.inc"
         iflow=nint(jamp2014(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2014(i)
         enddo
         goto 20
      elseif(str.eq."3-3100002510000350") then
         include "leshouches_R_015.inc"
         iflow=nint(jamp2015(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2015(i)
         enddo
         goto 20
      elseif(str.eq."30100002510000353") then
         include "leshouches_R_016.inc"
         iflow=nint(jamp2016(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2016(i)
         enddo
         goto 20
      elseif(str.eq."-55100002510000350") then
         include "leshouches_R_017.inc"
         iflow=nint(jamp2017(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2017(i)
         enddo
         goto 20
      elseif(str.eq."-5010000251000035-5") then
         include "leshouches_R_018.inc"
         iflow=nint(jamp2018(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2018(i)
         enddo
         goto 20
      elseif(str.eq."5-5100002510000350") then
         include "leshouches_R_019.inc"
         iflow=nint(jamp2019(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2019(i)
         enddo
         goto 20
      elseif(str.eq."50100002510000355") then
         include "leshouches_R_020.inc"
         iflow=nint(jamp2020(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2020(i)
         enddo
         goto 20
      elseif(str.eq."0-110000251000035-1") then
         include "leshouches_R_021.inc"
         iflow=nint(jamp2021(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2021(i)
         enddo
         goto 20
      elseif(str.eq."01100002510000351") then
         include "leshouches_R_022.inc"
         iflow=nint(jamp2022(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2022(i)
         enddo
         goto 20
      elseif(str.eq."0-210000251000035-2") then
         include "leshouches_R_023.inc"
         iflow=nint(jamp2023(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2023(i)
         enddo
         goto 20
      elseif(str.eq."02100002510000352") then
         include "leshouches_R_024.inc"
         iflow=nint(jamp2024(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2024(i)
         enddo
         goto 20
      elseif(str.eq."0-410000251000035-4") then
         include "leshouches_R_025.inc"
         iflow=nint(jamp2025(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2025(i)
         enddo
         goto 20
      elseif(str.eq."04100002510000354") then
         include "leshouches_R_026.inc"
         iflow=nint(jamp2026(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2026(i)
         enddo
         goto 20
      elseif(str.eq."0-310000251000035-3") then
         include "leshouches_R_027.inc"
         iflow=nint(jamp2027(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2027(i)
         enddo
         goto 20
      elseif(str.eq."03100002510000353") then
         include "leshouches_R_028.inc"
         iflow=nint(jamp2028(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2028(i)
         enddo
         goto 20
      elseif(str.eq."0-510000251000035-5") then
         include "leshouches_R_029.inc"
         iflow=nint(jamp2029(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2029(i)
         enddo
         goto 20
      elseif(str.eq."05100002510000355") then
         include "leshouches_R_030.inc"
         iflow=nint(jamp2030(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2030(i)
         enddo
         goto 20
         
      elseif(str.eq."-11100003510000350") then
         include "leshouches_R_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-1010000351000035-1") then
         include "leshouches_R_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."1-1100003510000350") then
         include "leshouches_R_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."10100003510000351") then
         include "leshouches_R_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-22100003510000350") then
         include "leshouches_R_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."-2010000351000035-2") then
         include "leshouches_R_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."2-2100003510000350") then
         include "leshouches_R_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."20100003510000352") then
         include "leshouches_R_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."-44100003510000350") then
         include "leshouches_R_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."-4010000351000035-4") then
         include "leshouches_R_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."4-4100003510000350") then
         include "leshouches_R_011.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."40100003510000354") then
         include "leshouches_R_012.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
      elseif(str.eq."-33100003510000350") then
         include "leshouches_R_013.inc"
         iflow=nint(jamp2013(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2013(i)
         enddo
         goto 20
      elseif(str.eq."-3010000351000035-3") then
         include "leshouches_R_014.inc"
         iflow=nint(jamp2014(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2014(i)
         enddo
         goto 20
      elseif(str.eq."3-3100003510000350") then
         include "leshouches_R_015.inc"
         iflow=nint(jamp2015(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2015(i)
         enddo
         goto 20
      elseif(str.eq."30100003510000353") then
         include "leshouches_R_016.inc"
         iflow=nint(jamp2016(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2016(i)
         enddo
         goto 20
      elseif(str.eq."-55100003510000350") then
         include "leshouches_R_017.inc"
         iflow=nint(jamp2017(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2017(i)
         enddo
         goto 20
      elseif(str.eq."-5010000351000035-5") then
         include "leshouches_R_018.inc"
         iflow=nint(jamp2018(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2018(i)
         enddo
         goto 20
      elseif(str.eq."5-5100003510000350") then
         include "leshouches_R_019.inc"
         iflow=nint(jamp2019(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2019(i)
         enddo
         goto 20
      elseif(str.eq."50100003510000355") then
         include "leshouches_R_020.inc"
         iflow=nint(jamp2020(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2020(i)
         enddo
         goto 20
      elseif(str.eq."0-110000351000035-1") then
         include "leshouches_R_021.inc"
         iflow=nint(jamp2021(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2021(i)
         enddo
         goto 20
      elseif(str.eq."01100003510000351") then
         include "leshouches_R_022.inc"
         iflow=nint(jamp2022(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2022(i)
         enddo
         goto 20
      elseif(str.eq."0-210000351000035-2") then
         include "leshouches_R_023.inc"
         iflow=nint(jamp2023(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2023(i)
         enddo
         goto 20
      elseif(str.eq."02100003510000352") then
         include "leshouches_R_024.inc"
         iflow=nint(jamp2024(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2024(i)
         enddo
         goto 20
      elseif(str.eq."0-410000351000035-4") then
         include "leshouches_R_025.inc"
         iflow=nint(jamp2025(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2025(i)
         enddo
         goto 20
      elseif(str.eq."04100003510000354") then
         include "leshouches_R_026.inc"
         iflow=nint(jamp2026(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2026(i)
         enddo
         goto 20
      elseif(str.eq."0-310000351000035-3") then
         include "leshouches_R_027.inc"
         iflow=nint(jamp2027(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2027(i)
         enddo
         goto 20
      elseif(str.eq."03100003510000353") then
         include "leshouches_R_028.inc"
         iflow=nint(jamp2028(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2028(i)
         enddo
         goto 20
      elseif(str.eq."0-510000351000035-5") then
         include "leshouches_R_029.inc"
         iflow=nint(jamp2029(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2029(i)
         enddo
         goto 20
      elseif(str.eq."05100003510000355") then
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
      do while(jamp2cum(ifl).lt.xtarget)
         ifl=ifl+1
      enddo
      do i=1,2
         do j=1,nexternal
            color(i,j)=ICOLUP(i,j,ifl)
         enddo
      enddo
      
      return
      end


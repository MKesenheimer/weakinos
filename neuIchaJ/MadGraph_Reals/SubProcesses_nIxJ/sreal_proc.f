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
      
      if(str.eq."1-21000022-10000240") then ! d ubar -> ni + cha1- g
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DUX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."101000022-10000242") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DG_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."-211000022-10000240") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_UXD_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-201000022-1000024-1") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_UXG_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."-431000022-10000240") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UXD_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-401000022-1000024-3") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UXG_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."3-41000022-10000240") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_DUX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."301000022-10000244") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_DG_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."011000022-10000242") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_GD_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."0-21000022-1000024-1") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_GUX_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."0-41000022-1000024-3") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_GUX_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."031000022-10000244") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_GD_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."-12100002210000240") then
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DXU_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-1010000221000024-2") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DXG_NIXJUX(p,wgt)
         goto 20
      elseif(str.eq."2-1100002210000240") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_UDX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."20100002210000241") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_UG_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."4-3100002210000240") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UDX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."40100002210000243") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UG_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."-34100002210000240") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_DXU_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-3010000221000024-4") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_DXG_NIXJUX(p,wgt)
         goto 20
      elseif(str.eq."0-110000221000024-2") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_GDX_NIXJUX(p,wgt)
         goto 20
      elseif(str.eq."02100002210000241") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_GU_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."04100002210000243") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_GU_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."0-310000221000024-4") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_GDX_NIXJUX(p,wgt)
         goto 20
         
      elseif(str.eq."1-21000023-10000240") then ! d ubar -> ni + cha1- g
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DUX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."101000023-10000242") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DG_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."-211000023-10000240") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_UXD_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-201000023-1000024-1") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_UXG_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."-431000023-10000240") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UXD_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-401000023-1000024-3") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UXG_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."3-41000023-10000240") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_DUX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."301000023-10000244") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_DG_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."011000023-10000242") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_GD_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."0-21000023-1000024-1") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_GUX_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."0-41000023-1000024-3") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_GUX_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."031000023-10000244") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_GD_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."-12100002310000240") then
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DXU_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-1010000231000024-2") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DXG_NIXJUX(p,wgt)
         goto 20
      elseif(str.eq."2-1100002310000240") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_UDX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."20100002310000241") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_UG_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."4-3100002310000240") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UDX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."40100002310000243") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UG_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."-34100002310000240") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_DXU_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-3010000231000024-4") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_DXG_NIXJUX(p,wgt)
         goto 20
      elseif(str.eq."0-110000231000024-2") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_GDX_NIXJUX(p,wgt)
         goto 20
      elseif(str.eq."02100002310000241") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_GU_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."04100002310000243") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_GU_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."0-310000231000024-4") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_GDX_NIXJUX(p,wgt)
         goto 20
         
      elseif(str.eq."1-21000022-10000370") then ! d ubar -> ni + cha1- g
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DUX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."101000022-10000372") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DG_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."-211000022-10000370") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_UXD_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-201000022-1000037-1") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_UXG_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."-431000022-10000370") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UXD_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-401000022-1000037-3") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UXG_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."3-41000022-10000370") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_DUX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."301000022-10000374") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_DG_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."011000022-10000372") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_GD_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."0-21000022-1000037-1") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_GUX_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."0-41000022-1000037-3") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_GUX_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."031000022-10000374") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_GD_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."-12100002210000370") then
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DXU_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-1010000221000037-2") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DXG_NIXJUX(p,wgt)
         goto 20
      elseif(str.eq."2-1100002210000370") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_UDX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."20100002210000371") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_UG_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."4-3100002210000370") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UDX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."40100002210000373") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UG_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."-34100002210000370") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_DXU_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-3010000221000037-4") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_DXG_NIXJUX(p,wgt)
         goto 20
      elseif(str.eq."0-110000221000037-2") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_GDX_NIXJUX(p,wgt)
         goto 20
      elseif(str.eq."02100002210000371") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_GU_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."04100002210000373") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_GU_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."0-310000221000037-4") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_GDX_NIXJUX(p,wgt)
         goto 20
         
      elseif(str.eq."1-21000023-10000370") then ! d ubar -> ni + cha1- g
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DUX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."101000023-10000372") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DG_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."-211000023-10000370") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_UXD_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-201000023-1000037-1") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_UXG_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."-431000023-10000370") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UXD_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-401000023-1000037-3") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UXG_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."3-41000023-10000370") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_DUX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."301000023-10000374") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_DG_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."011000023-10000372") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_GD_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."0-21000023-1000037-1") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_GUX_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."0-41000023-1000037-3") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_GUX_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."031000023-10000374") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_GD_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."-12100002310000370") then
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DXU_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-1010000231000037-2") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DXG_NIXJUX(p,wgt)
         goto 20
      elseif(str.eq."2-1100002310000370") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_UDX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."20100002310000371") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_UG_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."4-3100002310000370") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UDX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."40100002310000373") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UG_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."-34100002310000370") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_DXU_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-3010000231000037-4") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_DXG_NIXJUX(p,wgt)
         goto 20
      elseif(str.eq."0-110000231000037-2") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_GDX_NIXJUX(p,wgt)
         goto 20
      elseif(str.eq."02100002310000371") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_GU_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."04100002310000373") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_GU_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."0-310000231000037-4") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_GDX_NIXJUX(p,wgt)
         goto 20
      
      elseif(str.eq."1-21000025-10000240") then ! d ubar -> ni + cha1- g
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DUX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."101000025-10000242") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DG_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."-211000025-10000240") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_UXD_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-201000025-1000024-1") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_UXG_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."-431000025-10000240") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UXD_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-401000025-1000024-3") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UXG_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."3-41000025-10000240") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_DUX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."301000025-10000244") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_DG_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."011000025-10000242") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_GD_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."0-21000025-1000024-1") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_GUX_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."0-41000025-1000024-3") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_GUX_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."031000025-10000244") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_GD_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."-12100002510000240") then
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DXU_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-1010000251000024-2") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DXG_NIXJUX(p,wgt)
         goto 20
      elseif(str.eq."2-1100002510000240") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_UDX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."20100002510000241") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_UG_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."4-3100002510000240") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UDX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."40100002510000243") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UG_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."-34100002510000240") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_DXU_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-3010000251000024-4") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_DXG_NIXJUX(p,wgt)
         goto 20
      elseif(str.eq."0-110000251000024-2") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_GDX_NIXJUX(p,wgt)
         goto 20
      elseif(str.eq."02100002510000241") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_GU_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."04100002510000243") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_GU_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."0-310000251000024-4") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_GDX_NIXJUX(p,wgt)
         goto 20
         
      elseif(str.eq."1-21000025-10000370") then ! d ubar -> ni + cha2- g
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DUX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."101000025-10000372") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DG_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."-211000025-10000370") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_UXD_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-201000025-1000037-1") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_UXG_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."-431000025-10000370") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UXD_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-401000025-1000037-3") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UXG_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."3-41000025-10000370") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_DUX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."301000025-10000374") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_DG_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."011000025-10000372") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_GD_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."0-21000025-1000037-1") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_GUX_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."0-41000025-1000037-3") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_GUX_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."031000025-10000374") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_GD_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."-12100002510000370") then
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DXU_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-1010000251000037-2") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DXG_NIXJUX(p,wgt)
         goto 20
      elseif(str.eq."2-1100002510000370") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_UDX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."20100002510000371") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_UG_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."4-3100002510000370") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UDX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."40100002510000373") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UG_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."-34100002510000370") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_DXU_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-3010000251000037-4") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_DXG_NIXJUX(p,wgt)
         goto 20
      elseif(str.eq."0-110000251000037-2") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_GDX_NIXJUX(p,wgt)
         goto 20
      elseif(str.eq."02100002510000371") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_GU_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."04100002510000373") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_GU_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."0-310000251000037-4") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_GDX_NIXJUX(p,wgt)
         goto 20
         
      elseif(str.eq."1-21000035-10000240") then ! d ubar -> ni + cha1- g
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DUX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."101000035-10000242") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DG_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."-211000035-10000240") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_UXD_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-201000035-1000024-1") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_UXG_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."-431000035-10000240") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UXD_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-401000035-1000024-3") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UXG_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."3-41000035-10000240") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_DUX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."301000035-10000244") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_DG_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."011000035-10000242") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_GD_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."0-21000035-1000024-1") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_GUX_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."0-41000035-1000024-3") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_GUX_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."031000035-10000244") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_GD_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."-12100003510000240") then
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DXU_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-1010000351000024-2") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DXG_NIXJUX(p,wgt)
         goto 20
      elseif(str.eq."2-1100003510000240") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_UDX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."20100003510000241") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_UG_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."4-3100003510000240") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UDX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."40100003510000243") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UG_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."-34100003510000240") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_DXU_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-3010000351000024-4") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_DXG_NIXJUX(p,wgt)
         goto 20
      elseif(str.eq."0-110000351000024-2") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_GDX_NIXJUX(p,wgt)
         goto 20
      elseif(str.eq."02100003510000241") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_GU_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."04100003510000243") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_GU_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."0-310000351000024-4") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_GDX_NIXJUX(p,wgt)
         goto 20
         
      elseif(str.eq."1-21000035-10000370") then ! d ubar -> ni + cha2- g
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DUX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."101000035-10000372") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DG_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."-211000035-10000370") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_UXD_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-201000035-1000037-1") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_UXG_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."-431000035-10000370") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UXD_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-401000035-1000037-3") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UXG_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."3-41000035-10000370") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_DUX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."301000035-10000374") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_DG_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."011000035-10000372") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_GD_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."0-21000035-1000037-1") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_GUX_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."0-41000035-1000037-3") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_GUX_NIXJDX(p,wgt)
         goto 20
      elseif(str.eq."031000035-10000374") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_GD_NIXJU(p,wgt)
         goto 20
      elseif(str.eq."-12100003510000370") then
         !call srealmtrx_001(p,wgt)
         call SMATRIX_DXU_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-1010000351000037-2") then
         !call srealmtrx_002(p,wgt)
         call SMATRIX_DXG_NIXJUX(p,wgt)
         goto 20
      elseif(str.eq."2-1100003510000370") then
         !call srealmtrx_003(p,wgt)
         call SMATRIX_UDX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."20100003510000371") then
         !call srealmtrx_004(p,wgt)
         call SMATRIX_UG_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."4-3100003510000370") then
         !call srealmtrx_005(p,wgt)
         call SMATRIX_UDX_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."40100003510000373") then
         !call srealmtrx_006(p,wgt)
         call SMATRIX_UG_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."-34100003510000370") then
         !call srealmtrx_007(p,wgt)
         call SMATRIX_DXU_NIXJG(p,wgt)
         goto 20
      elseif(str.eq."-3010000351000037-4") then
         !call srealmtrx_008(p,wgt)
         call SMATRIX_DXG_NIXJUX(p,wgt)
         goto 20
      elseif(str.eq."0-110000351000037-2") then
         !call srealmtrx_009(p,wgt)
         call SMATRIX_GDX_NIXJUX(p,wgt)
         goto 20
      elseif(str.eq."02100003510000371") then
         !call srealmtrx_010(p,wgt)
         call SMATRIX_GU_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."04100003510000373") then
         !call srealmtrx_011(p,wgt)
         call SMATRIX_GU_NIXJD(p,wgt)
         goto 20
      elseif(str.eq."0-310000351000037-4") then
         !call srealmtrx_012(p,wgt)
         call SMATRIX_GDX_NIXJUX(p,wgt)
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
      include "maxamps_nixj.inc"
      Double Precision amp2001(maxamps), jamp2001(0:maxflow)
      common/to_amps_dux/amp2001,jamp2001
      Double Precision amp2002(maxamps), jamp2002(0:maxflow)
      common/to_amps_dg/amp2002,jamp2002
      Double Precision amp2003(maxamps), jamp2003(0:maxflow)
      common/to_amps_uxd/amp2003,jamp2003
      Double Precision amp2004(maxamps), jamp2004(0:maxflow)
      common/to_amps_uxg/amp2004,jamp2004
      Double Precision amp2005(maxamps), jamp2005(0:maxflow)
      !common/to_amps_uxd/amp2005,jamp2005
      Double Precision amp2006(maxamps), jamp2006(0:maxflow)
      !common/to_amps_uxg/amp2006,jamp2006
      Double Precision amp2007(maxamps), jamp2007(0:maxflow)
      !common/to_amps_dux/amp2007,jamp2007
      Double Precision amp2008(maxamps), jamp2008(0:maxflow)
      !common/to_amps_dg/amp2008,jamp2008
      Double Precision amp2009(maxamps), jamp2009(0:maxflow)
      common/to_amps_gd/amp2009,jamp2009
      Double Precision amp2010(maxamps), jamp2010(0:maxflow)
      common/to_amps_gux/amp2010,jamp2010
      Double Precision amp2011(maxamps), jamp2011(0:maxflow)
      !common/to_amps_gux/amp2011,jamp2011
      Double Precision amp2012(maxamps), jamp2012(0:maxflow)
      !common/to_amps_gd/amp2012,jamp2012
      double precision jamp2cum(0:maxamps)
      integer ICOLUP(2,nexternal,maxamps)
      integer color(2,nexternal),color1(2,nexternal)
      double precision random,xtarget
      external random
      integer legs(nexternal),lstr,i,j
      character*140 str
      integer iflow,ifl
      
      equivalence(amp2001,amp2007)
      equivalence(jamp2001,jamp2007)
      
      equivalence(amp2002,amp2008)
      equivalence(jamp2002,jamp2008)
      
      equivalence(amp2003,amp2005)
      equivalence(jamp2003,jamp2005)
      
      equivalence(amp2004,amp2006)
      equivalence(jamp2004,jamp2006)
      
      equivalence(amp2009,amp2012)
      equivalence(jamp2009,jamp2012)
      
      equivalence(amp2010,amp2011)
      equivalence(jamp2010,jamp2011)
      
      call convert_to_string(nexternal,legs,str,lstr)
      
      if(str.eq."1-21000022-10000240") then
         include "leshouches_R_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."101000022-10000242") then
         include "leshouches_R_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."-211000022-10000240") then
         include "leshouches_R_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."-201000022-1000024-1") then
         include "leshouches_R_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-431000022-10000240") then
         include "leshouches_R_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."-401000022-1000024-3") then
         include "leshouches_R_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."3-41000022-10000240") then
         include "leshouches_R_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."301000022-10000244") then
         include "leshouches_R_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."011000022-10000242") then
         include "leshouches_R_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."0-21000022-1000024-1") then
         include "leshouches_R_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."0-41000022-1000024-3") then
         include "leshouches_R_011.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."031000022-10000244") then
         include "leshouches_R_012.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
      elseif(str.eq."-12100002210000240") then
         include "leshouches_R_013.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-1010000221000024-2") then
         include "leshouches_R_014.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."2-1100002210000240") then
         include "leshouches_R_015.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."20100002210000241") then
         include "leshouches_R_016.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."4-3100002210000240") then
         include "leshouches_R_017.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."40100002210000243") then
         include "leshouches_R_018.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."-34100002210000240") then
         include "leshouches_R_019.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."-3010000221000024-4") then
         include "leshouches_R_020.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."0-110000221000024-2") then
         include "leshouches_R_021.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."02100002210000241") then
         include "leshouches_R_022.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."04100002210000243") then
         include "leshouches_R_023.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."0-310000221000024-4") then
         include "leshouches_R_024.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
         
      elseif(str.eq."1-21000023-10000240") then
         include "leshouches_R_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."101000023-10000242") then
         include "leshouches_R_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."-211000023-10000240") then
         include "leshouches_R_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."-201000023-1000024-1") then
         include "leshouches_R_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-431000023-10000240") then
         include "leshouches_R_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."-401000023-1000024-3") then
         include "leshouches_R_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."3-41000023-10000240") then
         include "leshouches_R_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."301000023-10000244") then
         include "leshouches_R_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."011000023-10000242") then
         include "leshouches_R_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."0-21000023-1000024-1") then
         include "leshouches_R_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."0-41000023-1000024-3") then
         include "leshouches_R_011.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."031000023-10000244") then
         include "leshouches_R_012.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
      elseif(str.eq."-12100002310000240") then
         include "leshouches_R_013.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-1010000231000024-2") then
         include "leshouches_R_014.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."2-1100002310000240") then
         include "leshouches_R_015.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."20100002310000241") then
         include "leshouches_R_016.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."4-3100002310000240") then
         include "leshouches_R_017.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."40100002310000243") then
         include "leshouches_R_018.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."-34100002310000240") then
         include "leshouches_R_019.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."-3010000231000024-4") then
         include "leshouches_R_020.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."0-110000231000024-2") then
         include "leshouches_R_021.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."02100002310000241") then
         include "leshouches_R_022.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."04100002310000243") then
         include "leshouches_R_023.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."0-310000231000024-4") then
         include "leshouches_R_024.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
         
      elseif(str.eq."1-21000022-10000370") then
         include "leshouches_R_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."101000022-10000372") then
         include "leshouches_R_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."-211000022-10000370") then
         include "leshouches_R_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."-201000022-1000037-1") then
         include "leshouches_R_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-431000022-10000370") then
         include "leshouches_R_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."-401000022-1000037-3") then
         include "leshouches_R_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."3-41000022-10000370") then
         include "leshouches_R_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."301000022-10000374") then
         include "leshouches_R_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."011000022-10000372") then
         include "leshouches_R_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."0-21000022-1000037-1") then
         include "leshouches_R_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."0-41000022-1000037-3") then
         include "leshouches_R_011.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."031000022-10000374") then
         include "leshouches_R_012.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
      elseif(str.eq."-12100002210000370") then
         include "leshouches_R_013.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-1010000221000037-2") then
         include "leshouches_R_014.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."2-1100002210000370") then
         include "leshouches_R_015.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."20100002210000371") then
         include "leshouches_R_016.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."4-3100002210000370") then
         include "leshouches_R_017.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."40100002210000373") then
         include "leshouches_R_018.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."-34100002210000370") then
         include "leshouches_R_019.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."-3010000221000037-4") then
         include "leshouches_R_020.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."0-110000221000037-2") then
         include "leshouches_R_021.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."02100002210000371") then
         include "leshouches_R_022.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."04100002210000373") then
         include "leshouches_R_023.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."0-310000221000037-4") then
         include "leshouches_R_024.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
      
      elseif(str.eq."1-21000023-10000370") then
         include "leshouches_R_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."101000023-10000372") then
         include "leshouches_R_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."-211000023-10000370") then
         include "leshouches_R_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."-201000023-1000037-1") then
         include "leshouches_R_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-431000023-10000370") then
         include "leshouches_R_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."-401000023-1000037-3") then
         include "leshouches_R_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."3-41000023-10000370") then
         include "leshouches_R_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."301000023-10000374") then
         include "leshouches_R_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."011000023-10000372") then
         include "leshouches_R_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."0-21000023-1000037-1") then
         include "leshouches_R_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."0-41000023-1000037-3") then
         include "leshouches_R_011.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."031000023-10000374") then
         include "leshouches_R_012.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
      elseif(str.eq."-12100002310000370") then
         include "leshouches_R_013.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-1010000231000037-2") then
         include "leshouches_R_014.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."2-1100002310000370") then
         include "leshouches_R_015.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."20100002310000371") then
         include "leshouches_R_016.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."4-3100002310000370") then
         include "leshouches_R_017.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."40100002310000373") then
         include "leshouches_R_018.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."-34100002310000370") then
         include "leshouches_R_019.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."-3010000231000037-4") then
         include "leshouches_R_020.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."0-110000231000037-2") then
         include "leshouches_R_021.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."02100002310000371") then
         include "leshouches_R_022.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."04100002310000373") then
         include "leshouches_R_023.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."0-310000231000037-4") then
         include "leshouches_R_024.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
      
      elseif(str.eq."1-21000025-10000240") then  ! n3 + x1
         include "leshouches_R_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."101000025-10000242") then
         include "leshouches_R_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."-211000025-10000240") then
         include "leshouches_R_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."-201000025-1000024-1") then
         include "leshouches_R_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-431000025-10000240") then
         include "leshouches_R_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."-401000025-1000024-3") then
         include "leshouches_R_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."3-41000025-10000240") then
         include "leshouches_R_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."301000025-10000244") then
         include "leshouches_R_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."011000025-10000242") then
         include "leshouches_R_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."0-21000025-1000024-1") then
         include "leshouches_R_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."0-41000025-1000024-3") then
         include "leshouches_R_011.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."031000025-10000244") then
         include "leshouches_R_012.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
      elseif(str.eq."-12100002510000240") then
         include "leshouches_R_013.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-1010000251000024-2") then
         include "leshouches_R_014.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."2-1100002510000240") then
         include "leshouches_R_015.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."20100002510000241") then
         include "leshouches_R_016.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."4-3100002510000240") then
         include "leshouches_R_017.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."40100002510000243") then
         include "leshouches_R_018.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."-34100002510000240") then
         include "leshouches_R_019.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."-3010000251000024-4") then
         include "leshouches_R_020.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."0-110000251000024-2") then
         include "leshouches_R_021.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."02100002510000241") then
         include "leshouches_R_022.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."04100002510000243") then
         include "leshouches_R_023.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."0-310000251000024-4") then
         include "leshouches_R_024.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
         
      elseif(str.eq."1-21000025-10000370") then  ! n3 + x2
         include "leshouches_R_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."101000025-10000372") then
         include "leshouches_R_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."-211000025-10000370") then
         include "leshouches_R_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."-201000025-1000037-1") then
         include "leshouches_R_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-431000025-10000370") then
         include "leshouches_R_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."-401000025-1000037-3") then
         include "leshouches_R_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."3-41000025-10000370") then
         include "leshouches_R_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."301000025-10000374") then
         include "leshouches_R_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."011000025-10000372") then
         include "leshouches_R_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."0-21000025-1000037-1") then
         include "leshouches_R_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."0-41000025-1000037-3") then
         include "leshouches_R_011.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."031000025-10000374") then
         include "leshouches_R_012.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
      elseif(str.eq."-12100002510000370") then
         include "leshouches_R_013.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-1010000251000037-2") then
         include "leshouches_R_014.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."2-1100002510000370") then
         include "leshouches_R_015.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."20100002510000371") then
         include "leshouches_R_016.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."4-3100002510000370") then
         include "leshouches_R_017.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."40100002510000373") then
         include "leshouches_R_018.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."-34100002510000370") then
         include "leshouches_R_019.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."-3010000251000037-4") then
         include "leshouches_R_020.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."0-110000251000037-2") then
         include "leshouches_R_021.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."02100002510000371") then
         include "leshouches_R_022.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."04100002510000373") then
         include "leshouches_R_023.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."0-310000251000037-4") then
         include "leshouches_R_024.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
         
      elseif(str.eq."1-21000035-10000240") then  ! n3 + x1
         include "leshouches_R_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."101000035-10000242") then
         include "leshouches_R_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."-211000035-10000240") then
         include "leshouches_R_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."-201000035-1000024-1") then
         include "leshouches_R_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-431000035-10000240") then
         include "leshouches_R_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."-401000035-1000024-3") then
         include "leshouches_R_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."3-41000035-10000240") then
         include "leshouches_R_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."301000035-10000244") then
         include "leshouches_R_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."011000035-10000242") then
         include "leshouches_R_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."0-21000035-1000024-1") then
         include "leshouches_R_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."0-41000035-1000024-3") then
         include "leshouches_R_011.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."031000035-10000244") then
         include "leshouches_R_012.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
      elseif(str.eq."-12100003510000240") then
         include "leshouches_R_013.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-1010000351000024-2") then
         include "leshouches_R_014.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."2-1100003510000240") then
         include "leshouches_R_015.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."20100003510000241") then
         include "leshouches_R_016.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."4-3100003510000240") then
         include "leshouches_R_017.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."40100003510000243") then
         include "leshouches_R_018.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."-34100003510000240") then
         include "leshouches_R_019.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."-3010000351000024-4") then
         include "leshouches_R_020.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."0-110000351000024-2") then
         include "leshouches_R_021.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."02100003510000241") then
         include "leshouches_R_022.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."04100003510000243") then
         include "leshouches_R_023.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."0-310000351000024-4") then
         include "leshouches_R_024.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
         
      elseif(str.eq."1-21000035-10000370") then  ! n3 + x2
         include "leshouches_R_001.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."101000035-10000372") then
         include "leshouches_R_002.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."-211000035-10000370") then
         include "leshouches_R_003.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."-201000035-1000037-1") then
         include "leshouches_R_004.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."-431000035-10000370") then
         include "leshouches_R_005.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."-401000035-1000037-3") then
         include "leshouches_R_006.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."3-41000035-10000370") then
         include "leshouches_R_007.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."301000035-10000374") then
         include "leshouches_R_008.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."011000035-10000372") then
         include "leshouches_R_009.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."0-21000035-1000037-1") then
         include "leshouches_R_010.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."0-41000035-1000037-3") then
         include "leshouches_R_011.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."031000035-10000374") then
         include "leshouches_R_012.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
         enddo
         goto 20
      elseif(str.eq."-12100003510000370") then
         include "leshouches_R_013.inc"
         iflow=nint(jamp2001(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2001(i)
         enddo
         goto 20
      elseif(str.eq."-1010000351000037-2") then
         include "leshouches_R_014.inc"
         iflow=nint(jamp2002(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2002(i)
         enddo
         goto 20
      elseif(str.eq."2-1100003510000370") then
         include "leshouches_R_015.inc"
         iflow=nint(jamp2003(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2003(i)
         enddo
         goto 20
      elseif(str.eq."20100003510000371") then
         include "leshouches_R_016.inc"
         iflow=nint(jamp2004(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2004(i)
         enddo
         goto 20
      elseif(str.eq."4-3100003510000370") then
         include "leshouches_R_017.inc"
         iflow=nint(jamp2005(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2005(i)
         enddo
         goto 20
      elseif(str.eq."40100003510000373") then
         include "leshouches_R_018.inc"
         iflow=nint(jamp2006(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2006(i)
         enddo
         goto 20
      elseif(str.eq."-34100003510000370") then
         include "leshouches_R_019.inc"
         iflow=nint(jamp2007(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2007(i)
         enddo
         goto 20
      elseif(str.eq."-3010000351000037-4") then
         include "leshouches_R_020.inc"
         iflow=nint(jamp2008(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2008(i)
         enddo
         goto 20
      elseif(str.eq."0-110000351000037-2") then
         include "leshouches_R_021.inc"
         iflow=nint(jamp2009(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2009(i)
         enddo
         goto 20
      elseif(str.eq."02100003510000371") then
         include "leshouches_R_022.inc"
         iflow=nint(jamp2010(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2010(i)
         enddo
         goto 20
      elseif(str.eq."04100003510000373") then
         include "leshouches_R_023.inc"
         iflow=nint(jamp2011(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2011(i)
         enddo
         goto 20
      elseif(str.eq."0-310000351000037-4") then
         include "leshouches_R_024.inc"
         iflow=nint(jamp2012(0))
         jamp2cum(0)=0d0
         do i=1,iflow
            jamp2cum(i)=jamp2cum(i-1)+jamp2012(i)
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
      
      
      
      

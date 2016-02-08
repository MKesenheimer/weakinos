      subroutine sreal_proc_res(p,legs,chan,wgt)
      implicit none
      include "nexternal.inc"
      include "coupl.inc"
      double precision p(0:3,nexternal),wgt
      integer legs(nexternal),lstr
      character*140 str
      character*4 chan
      call convert_to_string(nexternal,legs,str,lstr)
      
#ifdef DEBUGQ
      print*
      print*, "[DEBUG] str = ", str
      !print*, "p1",p(0:3,1)
      !print*, "p2",p(0:3,2)
      !print*, "p3",p(0:3,3)
      !print*, "p4",p(0:3,4)
      !print*, "p5",p(0:3,5)
      !print*, "legs",legs
      print*, chan,wgt
#endif

      ! here: only left handed squarks can become resonant
      if( chan.ne."ul35" .and. chan.ne."ul45" .and. 
     &    chan.ne."dl35" .and. chan.ne."dl45" .and.
     &    chan.ne."allr ") then
        print*,"wrong channel number: ", chan
        print*,"choose ul35, ul45, dl35, dl45 or allr."
        stop
      endif

      wgt = 0D0      
      if (str.eq."-101000024-1000024-1") then
         !call srealmtrx_002_RES(p,chan,wgt)
         call SMATRIX_DXG_XIXIDX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif (str.eq."101000024-10000241") then
         !call srealmtrx_004_RES(p,chan,wgt)
         call SMATRIX_DG_XIXID_RES(p,chan,wgt) ! >> Fehler in srealmtrx_004
         goto 20
      elseif (str.eq."-201000024-1000024-2") then
         !call srealmtrx_006_RES(p,chan,wgt)
         call SMATRIX_UXG_XIXIUX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif (str.eq."201000024-10000242") then
         !call srealmtrx_008_RES(p,chan,wgt)
         call SMATRIX_UG_XIXIU_RES(p,chan,wgt) ! >> Fehler in srealmtrx_008
         goto 20
      elseif (str.eq."-401000024-1000024-4") then
         !call srealmtrx_010_RES(p,chan,wgt)
         !call SMATRIX_CXG_XIXICX_RES(p,chan,wgt) ! >> OK!
         call SMATRIX_UXG_XIXIUX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif (str.eq."401000024-10000244") then
         !call srealmtrx_012_RES(p,chan,wgt)
         !call SMATRIX_CG_XIXIC_RES(p,chan,wgt) ! >> Fehler in srealmtrx_012
         call SMATRIX_UG_XIXIU_RES(p,chan,wgt)
         goto 20
      elseif (str.eq."-301000024-1000024-3") then
         !call srealmtrx_014_RES(p,chan,wgt)
         !call SMATRIX_SXG_XIXISX_RES(p,chan,wgt) ! >> OK!
         call SMATRIX_DXG_XIXIDX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif (str.eq."301000024-10000243") then
         !call srealmtrx_016_RES(p,chan,wgt)
         !call SMATRIX_SG_XIXIS_RES(p,chan,wgt) ! >> Fehler in srealmtrx_016
         call SMATRIX_DG_XIXID_RES(p,chan,wgt) 
         goto 20
      elseif (str.eq."-501000024-1000024-5") then
         !call srealmtrx_018_RES(p,chan,wgt)
         !call SMATRIX_BXG_XIXIBX_RES(p,chan,wgt) ! >> OK!
         call SMATRIX_DXG_XIXIDX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif (str.eq."501000024-10000245") then
         !call srealmtrx_020_RES(p,chan,wgt)
         !call SMATRIX_BG_XIXIB_RES(p,chan,wgt) ! >> Fehler in srealmtrx_020
         call SMATRIX_DG_XIXID_RES(p,chan,wgt)
         goto 20
      elseif (str.eq."0-11000024-1000024-1") then
         !call srealmtrx_021_RES(p,chan,wgt)
         call SMATRIX_GDX_XIXIDX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif (str.eq."011000024-10000241") then
         !call srealmtrx_022_RES(p,chan,wgt)
         call SMATRIX_GD_XIXID_RES(p,chan,wgt) ! >> Fehler in srealmtrx_022
         goto 20
      elseif (str.eq."0-21000024-1000024-2") then
         !call srealmtrx_023_RES(p,chan,wgt)
         call SMATRIX_GUX_XIXIUX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif (str.eq."021000024-10000242") then
         !call srealmtrx_024_RES(p,chan,wgt)
         call SMATRIX_GU_XIXIU_RES(p,chan,wgt) ! >> Fehler in srealmtrx_024
         goto 20
      elseif (str.eq."0-41000024-1000024-4") then
         !call srealmtrx_025_RES(p,chan,wgt)
         !call SMATRIX_GCX_XIXICX_RES(p,chan,wgt) ! >> OK!
         call SMATRIX_GUX_XIXIUX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif (str.eq."041000024-10000244") then
         !call srealmtrx_026_RES(p,chan,wgt)
         !call SMATRIX_GC_XIXIC_RES(p,chan,wgt) ! >> Fehler in srealmtrx_026
         call SMATRIX_GU_XIXIU_RES(p,chan,wgt)
         goto 20
      elseif (str.eq."0-31000024-1000024-3") then
         !call srealmtrx_027_RES(p,chan,wgt)
         !call SMATRIX_GSX_XIXISX_RES(p,chan,wgt) ! >> OK!
         call SMATRIX_GDX_XIXIDX_RES(p,chan,wgt)
         goto 20
      elseif (str.eq."031000024-10000243") then
         !call srealmtrx_028_RES(p,chan,wgt)
         !call SMATRIX_GS_XIXIS_RES(p,chan,wgt) ! >> Fehler in srealmtrx_028
         call SMATRIX_GD_XIXID_RES(p,chan,wgt)
         goto 20
      elseif (str.eq."0-51000024-1000024-5") then
         !call srealmtrx_029_RES(p,chan,wgt)
         !call SMATRIX_GBX_XIXIBX_RES(p,chan,wgt) ! >> OK!
         call SMATRIX_GDX_XIXIDX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif (str.eq."051000024-10000245") then
         !call srealmtrx_030_RES(p,chan,wgt)
         !call SMATRIX_GB_XIXIB_RES(p,chan,wgt) ! >> Fehler in srealmtrx_030
         call SMATRIX_GD_XIXID_RES(p,chan,wgt)
         goto 20
         
      elseif(str.eq."-101000024-1000037-1") then
         !call srealmtrx_002_RES(p,chan,wgt)
         call SMATRIX_DXG_XIXJDX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif(str.eq."101000024-10000371") then
         !call srealmtrx_004_RES(p,chan,wgt)
         call SMATRIX_DG_XIXJD_RES(p,chan,wgt) ! >> Fehler in srealmtrx_004
         goto 20
      elseif(str.eq."-201000024-1000037-2") then
         !call srealmtrx_006_RES(p,chan,wgt)
         call SMATRIX_UXG_XIXJUX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif(str.eq."201000024-10000372") then
         !call srealmtrx_008_RES(p,chan,wgt)
         call SMATRIX_UG_XIXJU_RES(p,chan,wgt) ! >> Fehler in srealmtrx_008
         goto 20
      elseif(str.eq."-401000024-1000037-4") then
         !call srealmtrx_010_RES(p,chan,wgt)
         !call SMATRIX_CXG_XIXJCX_RES(p,chan,wgt) ! >> OK!
         call SMATRIX_UXG_XIXJUX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif(str.eq."401000024-10000374") then
         !call srealmtrx_012_RES(p,chan,wgt)
         !call SMATRIX_CG_XIXJC_RES(p,chan,wgt) ! >> Fehler in srealmtrx_012
         call SMATRIX_UG_XIXJU_RES(p,chan,wgt)
         goto 20
      elseif(str.eq."-301000024-1000037-3") then
         !call srealmtrx_014_RES(p,chan,wgt)
         !call SMATRIX_SXG_XIXJSX_RES(p,chan,wgt) ! >> OK!
         call SMATRIX_DXG_XIXJDX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif(str.eq."301000024-10000373") then
         !call srealmtrx_016_RES(p,chan,wgt)
         !call SMATRIX_SG_XIXJS_RES(p,chan,wgt) ! >> Fehler in srealmtrx_016
         call SMATRIX_DG_XIXJD_RES(p,chan,wgt) 
         goto 20
      elseif(str.eq."-501000024-1000037-5") then
         !call srealmtrx_018_RES(p,chan,wgt)
         !call SMATRIX_BXG_XIXJBX_RES(p,chan,wgt) ! >> OK!
         call SMATRIX_DXG_XIXJDX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif(str.eq."501000024-10000375") then
         !call srealmtrx_020_RES(p,chan,wgt)
         !call SMATRIX_BG_XIXJB_RES(p,chan,wgt) ! >> Fehler in srealmtrx_020
         call SMATRIX_DG_XIXJD_RES(p,chan,wgt)
         goto 20
      elseif(str.eq."0-11000024-1000037-1") then
         !call srealmtrx_021_RES(p,chan,wgt)
         call SMATRIX_GDX_XIXJDX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif(str.eq."011000024-10000371") then
         !call srealmtrx_022_RES(p,chan,wgt)
         call SMATRIX_GD_XIXJD_RES(p,chan,wgt) ! >> Fehler in srealmtrx_022
         goto 20
      elseif(str.eq."0-21000024-1000037-2") then
         !call srealmtrx_023_RES(p,chan,wgt)
         call SMATRIX_GUX_XIXJUX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif(str.eq."021000024-10000372") then
         !call srealmtrx_024_RES(p,chan,wgt)
         call SMATRIX_GU_XIXJU_RES(p,chan,wgt) ! >> Fehler in srealmtrx_024
         goto 20
      elseif(str.eq."0-41000024-1000037-4") then
         !call srealmtrx_025_RES(p,chan,wgt)
         !call SMATRIX_GCX_XIXJCX_RES(p,chan,wgt) ! >> OK!
         call SMATRIX_GUX_XIXJUX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif(str.eq."041000024-10000374") then
         !call srealmtrx_026_RES(p,chan,wgt)
         !call SMATRIX_GC_XIXJC_RES(p,chan,wgt) ! >> Fehler in srealmtrx_026
         call SMATRIX_GU_XIXJU_RES(p,chan,wgt)
         goto 20
      elseif(str.eq."0-31000024-1000037-3") then
         !call srealmtrx_027_RES(p,chan,wgt)
         !call SMATRIX_GSX_XIXJSX_RES(p,chan,wgt) ! >> OK!
         call SMATRIX_GDX_XIXJDX_RES(p,chan,wgt)
         goto 20
      elseif(str.eq."031000024-10000373") then
         !call srealmtrx_028_RES(p,chan,wgt)
         !call SMATRIX_GS_XIXJS_RES(p,chan,wgt) ! >> Fehler in srealmtrx_028
         call SMATRIX_GD_XIXJD_RES(p,chan,wgt)
         goto 20
      elseif(str.eq."0-51000024-1000037-5") then
         !call srealmtrx_029_RES(p,chan,wgt)
         !call SMATRIX_GBX_XIXJBX_RES(p,chan,wgt) ! >> OK!
         call SMATRIX_GDX_XIXJDX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif(str.eq."051000024-10000375") then
         !call srealmtrx_030_RES(p,chan,wgt)
         !call SMATRIX_GB_XIXJB_RES(p,chan,wgt) ! >> Fehler in srealmtrx_030
         call SMATRIX_GD_XIXJD_RES(p,chan,wgt)
         goto 20
         
      elseif(str.eq."-101000037-1000037-1") then
         !call srealmtrx_002_RES(p,chan,wgt)
         call SMATRIX_DXG_XIXIDX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif(str.eq."101000037-10000371") then
         !call srealmtrx_004_RES(p,chan,wgt)
         call SMATRIX_DG_XIXID_RES(p,chan,wgt) ! >> Fehler in srealmtrx_004
         goto 20
      elseif(str.eq."-201000037-1000037-2") then
         !call srealmtrx_006_RES(p,chan,wgt)
         call SMATRIX_UXG_XIXIUX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif(str.eq."201000037-10000372") then
         !call srealmtrx_008_RES(p,chan,wgt)
         call SMATRIX_UG_XIXIU_RES(p,chan,wgt) ! >> Fehler in srealmtrx_008
         goto 20
      elseif(str.eq."-401000037-1000037-4") then
         !call srealmtrx_010_RES(p,chan,wgt)
         !call SMATRIX_CXG_XIXICX_RES(p,chan,wgt) ! >> OK!
         call SMATRIX_UXG_XIXIUX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif(str.eq."401000037-10000374") then
         !call srealmtrx_012_RES(p,chan,wgt)
         !call SMATRIX_CG_XIXIC_RES(p,chan,wgt) ! >> Fehler in srealmtrx_012
         call SMATRIX_UG_XIXIU_RES(p,chan,wgt)
         goto 20
      elseif(str.eq."-301000037-1000037-3") then
         !call srealmtrx_014_RES(p,chan,wgt)
         !call SMATRIX_SXG_XIXISX_RES(p,chan,wgt) ! >> OK!
         call SMATRIX_DXG_XIXIDX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif(str.eq."301000037-10000373") then
         !call srealmtrx_016_RES(p,chan,wgt)
         !call SMATRIX_SG_XIXIS_RES(p,chan,wgt) ! >> Fehler in srealmtrx_016
         call SMATRIX_DG_XIXID_RES(p,chan,wgt) 
         goto 20
      elseif(str.eq."-501000037-1000037-5") then
         !call srealmtrx_018_RES(p,chan,wgt)
         !call SMATRIX_BXG_XIXIBX_RES(p,chan,wgt) ! >> OK!
         call SMATRIX_DXG_XIXIDX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif(str.eq."501000037-10000375") then
         !call srealmtrx_020_RES(p,chan,wgt)
         !call SMATRIX_BG_XIXIB_RES(p,chan,wgt) ! >> Fehler in srealmtrx_020
         call SMATRIX_DG_XIXID_RES(p,chan,wgt)
         goto 20
      elseif(str.eq."0-11000037-1000037-1") then
         !call srealmtrx_021_RES(p,chan,wgt)
         call SMATRIX_GDX_XIXIDX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif(str.eq."011000037-10000371") then
         !call srealmtrx_022_RES(p,chan,wgt)
         call SMATRIX_GD_XIXID_RES(p,chan,wgt) ! >> Fehler in srealmtrx_022
         goto 20
      elseif(str.eq."0-21000037-1000037-2") then
         !call srealmtrx_023_RES(p,chan,wgt)
         call SMATRIX_GUX_XIXIUX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif(str.eq."021000037-10000372") then
         !call srealmtrx_024_RES(p,chan,wgt)
         call SMATRIX_GU_XIXIU_RES(p,chan,wgt) ! >> Fehler in srealmtrx_024
         goto 20
      elseif(str.eq."0-41000037-1000037-4") then
         !call srealmtrx_025_RES(p,chan,wgt)
         !call SMATRIX_GCX_XIXICX_RES(p,chan,wgt) ! >> OK!
         call SMATRIX_GUX_XIXIUX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif(str.eq."041000037-10000374") then
         !call srealmtrx_026_RES(p,chan,wgt)
         !call SMATRIX_GC_XIXIC_RES(p,chan,wgt) ! >> Fehler in srealmtrx_026
         call SMATRIX_GU_XIXIU_RES(p,chan,wgt)
         goto 20
      elseif(str.eq."0-31000037-1000037-3") then
         !call srealmtrx_027_RES(p,chan,wgt)
         !call SMATRIX_GSX_XIXISX_RES(p,chan,wgt) ! >> OK!
         call SMATRIX_GDX_XIXIDX_RES(p,chan,wgt)
         goto 20
      elseif(str.eq."031000037-10000373") then
         !call srealmtrx_028_RES(p,chan,wgt)
         !call SMATRIX_GS_XIXIS_RES(p,chan,wgt) ! >> Fehler in srealmtrx_028
         call SMATRIX_GD_XIXID_RES(p,chan,wgt)
         goto 20
      elseif(str.eq."0-51000037-1000037-5") then
         !call srealmtrx_029_RES(p,chan,wgt)
         !call SMATRIX_GBX_XIXIBX_RES(p,chan,wgt) ! >> OK!
         call SMATRIX_GDX_XIXIDX_RES(p,chan,wgt) ! >> OK!
         goto 20
      elseif(str.eq."051000037-10000375") then
         !call srealmtrx_030_RES(p,chan,wgt)
         !call SMATRIX_GB_XIXIB_RES(p,chan,wgt) ! >> Fehler in srealmtrx_030
         call SMATRIX_GD_XIXID_RES(p,chan,wgt)
         goto 20
      endif   
      
 20   continue
 
      return
      end

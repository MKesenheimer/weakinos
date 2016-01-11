      subroutine coup9(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 9
      if(readlha) then


      MGVX88(1) = (2*EE*PvL0hTo2)/3. + (2*EE*PvL1hTo2)/3.
      MGVX88(2) = (2*EE*PvR1hTo2)/3. + (2*EE*PvR2hTo2)/3.
      MGVX89(1) = (2*EE*vL0hTo*vL0to)/3. + (2*EE*vL1hTo*vL1to)/3.
      MGVX89(2) = (2*EE*vR1hTo*vR1to)/3. + (2*EE*vR2hTo*vR2to)/3.
      MGVX91(1) = (2*EE*PvL0up2)/3. + (2*EE*PvL1up2)/3.
      MGVX91(2) = (2*EE)/3.
      MGVX93(1) = (2*EE*PvL0ch2)/3. + (2*EE*PvL1ch2)/3.
      MGVX93(2) = (2*EE*PvR1ch2)/3. + (2*EE*PvR2ch2)/3.
      MGVX95(1) = (2*EE*PvL0to2)/3. + (2*EE*PvL1to2)/3.
      MGVX95(2) = (2*EE*PvR1to2)/3. + (2*EE*PvR2to2)/3.
      MGVX108(1) = (g0*v0W*vL0do*vL0hUp)/Sqrt2 + (gt*v1W*vL1do*vL1h
     +Up)/Sqrt2
      MGVX108(2) = 0
      MGVX109(1) = (g0*v0W*vL0hDo*vL0hUp)/Sqrt2 + (gt*v1W*vL1hDo*vL
     +1hUp)/Sqrt2
      MGVX109(2) = (gt*v1W)/Sqrt2
      MGVX110(1) = (g0*v0W*vL0hCh*vL0st)/Sqrt2 + (gt*v1W*vL1hCh*vL1
     +st)/Sqrt2
      MGVX110(2) = (gt*v1W*vR1hCh*vR1st)/Sqrt2
      MGVX111(1) = (g0*v0W*vL0hCh*vL0hSt)/Sqrt2 + (gt*v1W*vL1hCh*vL
     +1hSt)/Sqrt2
      MGVX111(2) = (gt*v1W*vR1hCh*vR1hSt)/Sqrt2
      MGVX112(1) = (g0*v0W*vL0bo*vL0hTo)/Sqrt2 + (gt*v1W*vL1bo*vL1h
     +To)/Sqrt2
      MGVX112(2) = (gt*v1W*vR1bo*vR1hTo)/Sqrt2
      MGVX113(1) = (g0*v0W*vL0hBo*vL0hTo)/Sqrt2 + (gt*v1W*vL1hBo*vL
     +1hTo)/Sqrt2
      MGVX113(2) = (gt*v1W*vR1hBo*vR1hTo)/Sqrt2
      MGVX114(1) = (g0*v0W*vL0do*vL0up)/Sqrt2 + (gt*v1W*vL1do*vL1up
     +)/Sqrt2
      MGVX114(2) = 0
      MGVX115(1) = (g0*v0W*vL0hDo*vL0up)/Sqrt2 + (gt*v1W*vL1hDo*vL1
     +up)/Sqrt2
      MGVX115(2) = 0
      MGVX116(1) = (g0*v0W*vL0ch*vL0st)/Sqrt2 + (gt*v1W*vL1ch*vL1st
     +)/Sqrt2
      MGVX116(2) = (gt*v1W*vR1ch*vR1st)/Sqrt2
      MGVX117(1) = (g0*v0W*vL0ch*vL0hSt)/Sqrt2 + (gt*v1W*vL1ch*vL1h
     +St)/Sqrt2
      MGVX117(2) = (gt*v1W*vR1ch*vR1hSt)/Sqrt2
      endif

      return
      end

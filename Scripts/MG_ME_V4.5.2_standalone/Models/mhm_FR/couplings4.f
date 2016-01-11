      subroutine coup4(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 4
      if(readlha) then
      MGVX46 = Pg02*Pv0ZP2*v0W*v0WP + Pgt2*Pv1ZP2*v1W*v1WP
      MGVX27 = -(Pg02*Pv0WP4) - Pgt2*Pv1WP4
      MGVX35 = Pg02*Pv0WP2*Pv0Z2 + Pgt2*Pv1WP2*Pv1Z2
      MGVX43 = Pg02*Pv0WP2*v0Z*v0ZP + Pgt2*Pv1WP2*v1Z*v1ZP
      MGVX47 = Pg02*Pv0WP2*Pv0ZP2 + Pgt2*Pv1WP2*Pv1ZP2


      MGVX72(1) = -(EE*PvL0el2) - EE*PvL1el2
      MGVX72(2) = -EE
      MGVX73(1) = -(EE*vL0el*vL0hEl) - EE*vL1el*vL1hEl
      MGVX73(2) = 0
      MGVX74(1) = -(EE*PvL0mu2) - EE*PvL1mu2
      MGVX74(2) = -(EE*PvR1mu2) - EE*PvR2mu2
      MGVX75(1) = -(EE*vL0hMu*vL0mu) - EE*vL1hMu*vL1mu
      MGVX75(2) = -(EE*vR1hMu*vR1mu) - EE*vR2hMu*vR2mu
      MGVX76(1) = -(EE*PvL0ta2) - EE*PvL1ta2
      MGVX76(2) = -(EE*PvR1ta2) - EE*PvR2ta2
      MGVX77(1) = -(EE*vL0hTa*vL0ta) - EE*vL1hTa*vL1ta
      MGVX77(2) = -(EE*vR1hTa*vR1ta) - EE*vR2hTa*vR2ta
      MGVX79(1) = -(EE*PvL0hEl2) - EE*PvL1hEl2
      MGVX79(2) = -EE
      MGVX81(1) = -(EE*PvL0hMu2) - EE*PvL1hMu2
      MGVX81(2) = -(EE*PvR1hMu2) - EE*PvR2hMu2
      MGVX83(1) = -(EE*PvL0hTa2) - EE*PvL1hTa2
      MGVX83(2) = -(EE*PvR1hTa2) - EE*PvR2hTa2
      MGVX96(1) = (g0*v0W*vL0el*vL0hNu)/Sqrt2 + (gt*v1W*vL1el*vL1hN
     +u)/Sqrt2
      MGVX96(2) = 0
      endif

      return
      end

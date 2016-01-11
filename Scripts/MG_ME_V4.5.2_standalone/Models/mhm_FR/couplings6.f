      subroutine coup6(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 6
      if(readlha) then


      MGVX148(1) = (g0*v0WP*vL0hNu*vL0ta)/Sqrt2 + (gt*v1WP*vL1hNu*v
     +L1ta)/Sqrt2
      MGVX148(2) = -((gt*v1WP*vR1ta)/Sqrt2)
      MGVX149(1) = (g0*v0WP*vL0hNu*vL0hTa)/Sqrt2 + (gt*v1WP*vL1hNu*
     +vL1hTa)/Sqrt2
      MGVX149(2) = -((gt*v1WP*vR1hTa)/Sqrt2)
      MGVX150(1) = (g0*v0WP*vL0el*vL0nu)/Sqrt2 + (gt*v1WP*vL1el*vL1
     +nu)/Sqrt2
      MGVX150(2) = 0
      MGVX151(1) = (g0*v0WP*vL0hEl*vL0nu)/Sqrt2 + (gt*v1WP*vL1hEl*v
     +L1nu)/Sqrt2
      MGVX151(2) = 0
      MGVX152(1) = (g0*v0WP*vL0mu*vL0nu)/Sqrt2 + (gt*v1WP*vL1mu*vL1
     +nu)/Sqrt2
      MGVX152(2) = 0
      MGVX153(1) = (g0*v0WP*vL0hMu*vL0nu)/Sqrt2 + (gt*v1WP*vL1hMu*v
     +L1nu)/Sqrt2
      MGVX153(2) = 0
      MGVX154(1) = (g0*v0WP*vL0nu*vL0ta)/Sqrt2 + (gt*v1WP*vL1nu*vL1
     +ta)/Sqrt2
      MGVX154(2) = 0
      MGVX155(1) = (g0*v0WP*vL0hTa*vL0nu)/Sqrt2 + (gt*v1WP*vL1hTa*v
     +L1nu)/Sqrt2
      MGVX155(2) = 0
      MGVX204(1) = -(g0*PvL0el2*v0Z)/2. - (gt*PvL1el2*v1Z)/2. - (gp
     +*PvL0el2*v2Z)/2. - (gp*PvL1el2*v2Z)/2.
      MGVX204(2) = -(gp*v2Z)
      MGVX205(1) = -(g0*v0Z*vL0el*vL0hEl)/2. - (gp*v2Z*vL0el*vL0hEl
     +)/2. - (gt*v1Z*vL1el*vL1hEl)/2. - (gp*v2Z*vL1el*vL1hEl)/2.
      MGVX205(2) = 0
      MGVX206(1) = -(g0*PvL0mu2*v0Z)/2. - (gt*PvL1mu2*v1Z)/2. - (gp
     +*PvL0mu2*v2Z)/2. - (gp*PvL1mu2*v2Z)/2.
      MGVX206(2) = -(gt*PvR1mu2*v1Z)/2. - (gp*PvR1mu2*v2Z)/2. - gp*
     +PvR2mu2*v2Z
      MGVX207(1) = -(g0*v0Z*vL0hMu*vL0mu)/2. - (gp*v2Z*vL0hMu*vL0mu
     +)/2. - (gt*v1Z*vL1hMu*vL1mu)/2. - (gp*v2Z*vL1hMu*vL1mu)/2.
      MGVX207(2) = -(gt*v1Z*vR1hMu*vR1mu)/2. - (gp*v2Z*vR1hMu*vR1mu
     +)/2. - gp*v2Z*vR2hMu*vR2mu
      MGVX208(1) = -(g0*PvL0ta2*v0Z)/2. - (gt*PvL1ta2*v1Z)/2. - (gp
     +*PvL0ta2*v2Z)/2. - (gp*PvL1ta2*v2Z)/2.
      MGVX208(2) = -(gt*PvR1ta2*v1Z)/2. - (gp*PvR1ta2*v2Z)/2. - gp*
     +PvR2ta2*v2Z
      MGVX209(1) = -(g0*v0Z*vL0hTa*vL0ta)/2. - (gp*v2Z*vL0hTa*vL0ta
     +)/2. - (gt*v1Z*vL1hTa*vL1ta)/2. - (gp*v2Z*vL1hTa*vL1ta)/2.
      MGVX209(2) = -(gt*v1Z*vR1hTa*vR1ta)/2. - (gp*v2Z*vR1hTa*vR1ta
     +)/2. - gp*v2Z*vR2hTa*vR2ta
      MGVX211(1) = -(g0*PvL0hEl2*v0Z)/2. - (gt*PvL1hEl2*v1Z)/2. - (
     +gp*PvL0hEl2*v2Z)/2. - (gp*PvL1hEl2*v2Z)/2.
      MGVX211(2) = -(gt*v1Z)/2. - (gp*v2Z)/2.
      endif

      return
      end

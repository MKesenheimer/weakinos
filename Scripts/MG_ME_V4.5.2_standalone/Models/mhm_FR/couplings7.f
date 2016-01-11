      subroutine coup7(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 7
      if(readlha) then


      MGVX213(1) = -(g0*PvL0hMu2*v0Z)/2. - (gt*PvL1hMu2*v1Z)/2. - (
     +gp*PvL0hMu2*v2Z)/2. - (gp*PvL1hMu2*v2Z)/2.
      MGVX213(2) = -(gt*PvR1hMu2*v1Z)/2. - (gp*PvR1hMu2*v2Z)/2. - g
     +p*PvR2hMu2*v2Z
      MGVX215(1) = -(g0*PvL0hTa2*v0Z)/2. - (gt*PvL1hTa2*v1Z)/2. - (
     +gp*PvL0hTa2*v2Z)/2. - (gp*PvL1hTa2*v2Z)/2.
      MGVX215(2) = -(gt*PvR1hTa2*v1Z)/2. - (gp*PvR1hTa2*v2Z)/2. - g
     +p*PvR2hTa2*v2Z
      MGVX216(1) = (g0*PvL0hNu2*v0Z)/2. + (gt*PvL1hNu2*v1Z)/2. - (g
     +p*PvL0hNu2*v2Z)/2. - (gp*PvL1hNu2*v2Z)/2.
      MGVX216(2) = (gt*v1Z)/2. - (gp*v2Z)/2.
      MGVX217(1) = (g0*v0Z*vL0hNu*vL0nu)/2. - (gp*v2Z*vL0hNu*vL0nu)
     +/2. + (gt*v1Z*vL1hNu*vL1nu)/2. - (gp*v2Z*vL1hNu*vL1nu)/2.
      MGVX217(2) = 0
      MGVX223(1) = (g0*PvL0nu2*v0Z)/2. + (gt*PvL1nu2*v1Z)/2. - (gp*
     +PvL0nu2*v2Z)/2. - (gp*PvL1nu2*v2Z)/2.
      MGVX223(2) = 0
      MGVX252(1) = -(g0*PvL0el2*v0ZP)/2. - (gt*PvL1el2*v1ZP)/2. - (
     +gp*PvL0el2*v2ZP)/2. - (gp*PvL1el2*v2ZP)/2.
      MGVX252(2) = -(gp*v2ZP)
      MGVX253(1) = -(g0*v0ZP*vL0el*vL0hEl)/2. - (gp*v2ZP*vL0el*vL0h
     +El)/2. - (gt*v1ZP*vL1el*vL1hEl)/2. - (gp*v2ZP*vL1el*vL1hEl)/
     +2.
      MGVX253(2) = 0
      MGVX254(1) = -(g0*PvL0mu2*v0ZP)/2. - (gt*PvL1mu2*v1ZP)/2. - (
     +gp*PvL0mu2*v2ZP)/2. - (gp*PvL1mu2*v2ZP)/2.
      MGVX254(2) = -(gt*PvR1mu2*v1ZP)/2. - (gp*PvR1mu2*v2ZP)/2. - g
     +p*PvR2mu2*v2ZP
      MGVX255(1) = -(g0*v0ZP*vL0hMu*vL0mu)/2. - (gp*v2ZP*vL0hMu*vL0
     +mu)/2. - (gt*v1ZP*vL1hMu*vL1mu)/2. - (gp*v2ZP*vL1hMu*vL1mu)/
     +2.
      MGVX255(2) = -(gt*v1ZP*vR1hMu*vR1mu)/2. - (gp*v2ZP*vR1hMu*vR1
     +mu)/2. - gp*v2ZP*vR2hMu*vR2mu
      MGVX256(1) = -(g0*PvL0ta2*v0ZP)/2. - (gt*PvL1ta2*v1ZP)/2. - (
     +gp*PvL0ta2*v2ZP)/2. - (gp*PvL1ta2*v2ZP)/2.
      MGVX256(2) = -(gt*PvR1ta2*v1ZP)/2. - (gp*PvR1ta2*v2ZP)/2. - g
     +p*PvR2ta2*v2ZP
      MGVX257(1) = -(g0*v0ZP*vL0hTa*vL0ta)/2. - (gp*v2ZP*vL0hTa*vL0
     +ta)/2. - (gt*v1ZP*vL1hTa*vL1ta)/2. - (gp*v2ZP*vL1hTa*vL1ta)/
     +2.
      MGVX257(2) = -(gt*v1ZP*vR1hTa*vR1ta)/2. - (gp*v2ZP*vR1hTa*vR1
     +ta)/2. - gp*v2ZP*vR2hTa*vR2ta
      MGVX259(1) = -(g0*PvL0hEl2*v0ZP)/2. - (gt*PvL1hEl2*v1ZP)/2. -
     + (gp*PvL0hEl2*v2ZP)/2. - (gp*PvL1hEl2*v2ZP)/2.
      MGVX259(2) = -(gt*v1ZP)/2. - (gp*v2ZP)/2.
      MGVX261(1) = -(g0*PvL0hMu2*v0ZP)/2. - (gt*PvL1hMu2*v1ZP)/2. -
     + (gp*PvL0hMu2*v2ZP)/2. - (gp*PvL1hMu2*v2ZP)/2.
      MGVX261(2) = -(gt*PvR1hMu2*v1ZP)/2. - (gp*PvR1hMu2*v2ZP)/2. -
     + gp*PvR2hMu2*v2ZP
      MGVX263(1) = -(g0*PvL0hTa2*v0ZP)/2. - (gt*PvL1hTa2*v1ZP)/2. -
     + (gp*PvL0hTa2*v2ZP)/2. - (gp*PvL1hTa2*v2ZP)/2.
      MGVX263(2) = -(gt*PvR1hTa2*v1ZP)/2. - (gp*PvR1hTa2*v2ZP)/2. -
     + gp*PvR2hTa2*v2ZP
      MGVX264(1) = (g0*PvL0hNu2*v0ZP)/2. + (gt*PvL1hNu2*v1ZP)/2. - 
     +(gp*PvL0hNu2*v2ZP)/2. - (gp*PvL1hNu2*v2ZP)/2.
      MGVX264(2) = (gt*v1ZP)/2. - (gp*v2ZP)/2.
      endif

      return
      end

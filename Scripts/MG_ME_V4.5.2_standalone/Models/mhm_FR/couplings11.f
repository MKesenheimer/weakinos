      subroutine coup11(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 11
      if(readlha) then


      MGVX193(1) = -(g0*v0Z*vL0do*vL0hDo)/2. + (gp*v2Z*vL0do*vL0hDo
     +)/6. - (gt*v1Z*vL1do*vL1hDo)/2. + (gp*v2Z*vL1do*vL1hDo)/6.
      MGVX193(2) = 0
      MGVX194(1) = -(g0*PvL0st2*v0Z)/2. - (gt*PvL1st2*v1Z)/2. + (gp
     +*PvL0st2*v2Z)/6. + (gp*PvL1st2*v2Z)/6.
      MGVX194(2) = -(gt*PvR1st2*v1Z)/2. + (gp*PvR1st2*v2Z)/6. - (gp
     +*PvR2st2*v2Z)/3.
      MGVX195(1) = -(g0*v0Z*vL0hSt*vL0st)/2. + (gp*v2Z*vL0hSt*vL0st
     +)/6. - (gt*v1Z*vL1hSt*vL1st)/2. + (gp*v2Z*vL1hSt*vL1st)/6.
      MGVX195(2) = -(gt*v1Z*vR1hSt*vR1st)/2. + (gp*v2Z*vR1hSt*vR1st
     +)/6. - (gp*v2Z*vR2hSt*vR2st)/3.
      MGVX196(1) = -(g0*PvL0bo2*v0Z)/2. - (gt*PvL1bo2*v1Z)/2. + (gp
     +*PvL0bo2*v2Z)/6. + (gp*PvL1bo2*v2Z)/6.
      MGVX196(2) = -(gt*PvR1bo2*v1Z)/2. + (gp*PvR1bo2*v2Z)/6. - (gp
     +*PvR2bo2*v2Z)/3.
      MGVX197(1) = -(g0*v0Z*vL0bo*vL0hBo)/2. + (gp*v2Z*vL0bo*vL0hBo
     +)/6. - (gt*v1Z*vL1bo*vL1hBo)/2. + (gp*v2Z*vL1bo*vL1hBo)/6.
      MGVX197(2) = -(gt*v1Z*vR1bo*vR1hBo)/2. + (gp*v2Z*vR1bo*vR1hBo
     +)/6. - (gp*v2Z*vR2bo*vR2hBo)/3.
      MGVX199(1) = -(g0*PvL0hDo2*v0Z)/2. - (gt*PvL1hDo2*v1Z)/2. + (
     +gp*PvL0hDo2*v2Z)/6. + (gp*PvL1hDo2*v2Z)/6.
      MGVX199(2) = -(gt*v1Z)/2. + (gp*v2Z)/6.
      MGVX201(1) = -(g0*PvL0hSt2*v0Z)/2. - (gt*PvL1hSt2*v1Z)/2. + (
     +gp*PvL0hSt2*v2Z)/6. + (gp*PvL1hSt2*v2Z)/6.
      MGVX201(2) = -(gt*PvR1hSt2*v1Z)/2. + (gp*PvR1hSt2*v2Z)/6. - (
     +gp*PvR2hSt2*v2Z)/3.
      MGVX203(1) = -(g0*PvL0hBo2*v0Z)/2. - (gt*PvL1hBo2*v1Z)/2. + (
     +gp*PvL0hBo2*v2Z)/6. + (gp*PvL1hBo2*v2Z)/6.
      MGVX203(2) = -(gt*PvR1hBo2*v1Z)/2. + (gp*PvR1hBo2*v2Z)/6. - (
     +gp*PvR2hBo2*v2Z)/3.
      MGVX228(1) = (g0*PvL0hUp2*v0Z)/2. + (gt*PvL1hUp2*v1Z)/2. + (g
     +p*PvL0hUp2*v2Z)/6. + (gp*PvL1hUp2*v2Z)/6.
      MGVX228(2) = (gt*v1Z)/2. + (gp*v2Z)/6.
      MGVX229(1) = (g0*v0Z*vL0hUp*vL0up)/2. + (gp*v2Z*vL0hUp*vL0up)
     +/6. + (gt*v1Z*vL1hUp*vL1up)/2. + (gp*v2Z*vL1hUp*vL1up)/6.
      MGVX229(2) = 0
      MGVX230(1) = (g0*PvL0hCh2*v0Z)/2. + (gt*PvL1hCh2*v1Z)/2. + (g
     +p*PvL0hCh2*v2Z)/6. + (gp*PvL1hCh2*v2Z)/6.
      MGVX230(2) = (gt*PvR1hCh2*v1Z)/2. + (gp*PvR1hCh2*v2Z)/6. + (2
     +*gp*PvR2hCh2*v2Z)/3.
      MGVX231(1) = (g0*v0Z*vL0ch*vL0hCh)/2. + (gp*v2Z*vL0ch*vL0hCh)
     +/6. + (gt*v1Z*vL1ch*vL1hCh)/2. + (gp*v2Z*vL1ch*vL1hCh)/6.
      MGVX231(2) = (gt*v1Z*vR1ch*vR1hCh)/2. + (gp*v2Z*vR1ch*vR1hCh)
     +/6. + (2*gp*v2Z*vR2ch*vR2hCh)/3.
      MGVX232(1) = (g0*PvL0hTo2*v0Z)/2. + (gt*PvL1hTo2*v1Z)/2. + (g
     +p*PvL0hTo2*v2Z)/6. + (gp*PvL1hTo2*v2Z)/6.
      MGVX232(2) = (gt*PvR1hTo2*v1Z)/2. + (gp*PvR1hTo2*v2Z)/6. + (2
     +*gp*PvR2hTo2*v2Z)/3.
      MGVX233(1) = (g0*v0Z*vL0hTo*vL0to)/2. + (gp*v2Z*vL0hTo*vL0to)
     +/6. + (gt*v1Z*vL1hTo*vL1to)/2. + (gp*v2Z*vL1hTo*vL1to)/6.
      MGVX233(2) = (gt*v1Z*vR1hTo*vR1to)/2. + (gp*v2Z*vR1hTo*vR1to)
     +/6. + (2*gp*v2Z*vR2hTo*vR2to)/3.
      MGVX235(1) = (g0*PvL0up2*v0Z)/2. + (gt*PvL1up2*v1Z)/2. + (gp*
     +PvL0up2*v2Z)/6. + (gp*PvL1up2*v2Z)/6.
      MGVX235(2) = (2*gp*v2Z)/3.
      endif

      return
      end

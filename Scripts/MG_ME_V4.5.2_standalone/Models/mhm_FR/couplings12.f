      subroutine coup12(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 12
      if(readlha) then


      MGVX237(1) = (g0*PvL0ch2*v0Z)/2. + (gt*PvL1ch2*v1Z)/2. + (gp*
     +PvL0ch2*v2Z)/6. + (gp*PvL1ch2*v2Z)/6.
      MGVX237(2) = (gt*PvR1ch2*v1Z)/2. + (gp*PvR1ch2*v2Z)/6. + (2*g
     +p*PvR2ch2*v2Z)/3.
      MGVX239(1) = (g0*PvL0to2*v0Z)/2. + (gt*PvL1to2*v1Z)/2. + (gp*
     +PvL0to2*v2Z)/6. + (gp*PvL1to2*v2Z)/6.
      MGVX239(2) = (gt*PvR1to2*v1Z)/2. + (gp*PvR1to2*v2Z)/6. + (2*g
     +p*PvR2to2*v2Z)/3.
      MGVX240(1) = -(g0*PvL0do2*v0ZP)/2. - (gt*PvL1do2*v1ZP)/2. + (
     +gp*PvL0do2*v2ZP)/6. + (gp*PvL1do2*v2ZP)/6.
      MGVX240(2) = -(gp*v2ZP)/3.
      MGVX241(1) = -(g0*v0ZP*vL0do*vL0hDo)/2. + (gp*v2ZP*vL0do*vL0h
     +Do)/6. - (gt*v1ZP*vL1do*vL1hDo)/2. + (gp*v2ZP*vL1do*vL1hDo)/
     +6.
      MGVX241(2) = 0
      MGVX242(1) = -(g0*PvL0st2*v0ZP)/2. - (gt*PvL1st2*v1ZP)/2. + (
     +gp*PvL0st2*v2ZP)/6. + (gp*PvL1st2*v2ZP)/6.
      MGVX242(2) = -(gt*PvR1st2*v1ZP)/2. + (gp*PvR1st2*v2ZP)/6. - (
     +gp*PvR2st2*v2ZP)/3.
      MGVX243(1) = -(g0*v0ZP*vL0hSt*vL0st)/2. + (gp*v2ZP*vL0hSt*vL0
     +st)/6. - (gt*v1ZP*vL1hSt*vL1st)/2. + (gp*v2ZP*vL1hSt*vL1st)/
     +6.
      MGVX243(2) = -(gt*v1ZP*vR1hSt*vR1st)/2. + (gp*v2ZP*vR1hSt*vR1
     +st)/6. - (gp*v2ZP*vR2hSt*vR2st)/3.
      MGVX244(1) = -(g0*PvL0bo2*v0ZP)/2. - (gt*PvL1bo2*v1ZP)/2. + (
     +gp*PvL0bo2*v2ZP)/6. + (gp*PvL1bo2*v2ZP)/6.
      MGVX244(2) = -(gt*PvR1bo2*v1ZP)/2. + (gp*PvR1bo2*v2ZP)/6. - (
     +gp*PvR2bo2*v2ZP)/3.
      MGVX245(1) = -(g0*v0ZP*vL0bo*vL0hBo)/2. + (gp*v2ZP*vL0bo*vL0h
     +Bo)/6. - (gt*v1ZP*vL1bo*vL1hBo)/2. + (gp*v2ZP*vL1bo*vL1hBo)/
     +6.
      MGVX245(2) = -(gt*v1ZP*vR1bo*vR1hBo)/2. + (gp*v2ZP*vR1bo*vR1h
     +Bo)/6. - (gp*v2ZP*vR2bo*vR2hBo)/3.
      MGVX247(1) = -(g0*PvL0hDo2*v0ZP)/2. - (gt*PvL1hDo2*v1ZP)/2. +
     + (gp*PvL0hDo2*v2ZP)/6. + (gp*PvL1hDo2*v2ZP)/6.
      MGVX247(2) = -(gt*v1ZP)/2. + (gp*v2ZP)/6.
      MGVX249(1) = -(g0*PvL0hSt2*v0ZP)/2. - (gt*PvL1hSt2*v1ZP)/2. +
     + (gp*PvL0hSt2*v2ZP)/6. + (gp*PvL1hSt2*v2ZP)/6.
      MGVX249(2) = -(gt*PvR1hSt2*v1ZP)/2. + (gp*PvR1hSt2*v2ZP)/6. -
     + (gp*PvR2hSt2*v2ZP)/3.
      MGVX251(1) = -(g0*PvL0hBo2*v0ZP)/2. - (gt*PvL1hBo2*v1ZP)/2. +
     + (gp*PvL0hBo2*v2ZP)/6. + (gp*PvL1hBo2*v2ZP)/6.
      MGVX251(2) = -(gt*PvR1hBo2*v1ZP)/2. + (gp*PvR1hBo2*v2ZP)/6. -
     + (gp*PvR2hBo2*v2ZP)/3.
      MGVX276(1) = (g0*PvL0hUp2*v0ZP)/2. + (gt*PvL1hUp2*v1ZP)/2. + 
     +(gp*PvL0hUp2*v2ZP)/6. + (gp*PvL1hUp2*v2ZP)/6.
      MGVX276(2) = (gt*v1ZP)/2. + (gp*v2ZP)/6.
      MGVX277(1) = (g0*v0ZP*vL0hUp*vL0up)/2. + (gp*v2ZP*vL0hUp*vL0u
     +p)/6. + (gt*v1ZP*vL1hUp*vL1up)/2. + (gp*v2ZP*vL1hUp*vL1up)/6
     +.
      MGVX277(2) = 0
      MGVX278(1) = (g0*PvL0hCh2*v0ZP)/2. + (gt*PvL1hCh2*v1ZP)/2. + 
     +(gp*PvL0hCh2*v2ZP)/6. + (gp*PvL1hCh2*v2ZP)/6.
      MGVX278(2) = (gt*PvR1hCh2*v1ZP)/2. + (gp*PvR1hCh2*v2ZP)/6. + 
     +(2*gp*PvR2hCh2*v2ZP)/3.
      MGVX279(1) = (g0*v0ZP*vL0ch*vL0hCh)/2. + (gp*v2ZP*vL0ch*vL0hC
     +h)/6. + (gt*v1ZP*vL1ch*vL1hCh)/2. + (gp*v2ZP*vL1ch*vL1hCh)/6
     +.
      MGVX279(2) = (gt*v1ZP*vR1ch*vR1hCh)/2. + (gp*v2ZP*vR1ch*vR1hC
     +h)/6. + (2*gp*v2ZP*vR2ch*vR2hCh)/3.
      endif

      return
      end

      subroutine coup13(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 13
      if(readlha) then


      MGVX280(1) = (g0*PvL0hTo2*v0ZP)/2. + (gt*PvL1hTo2*v1ZP)/2. + 
     +(gp*PvL0hTo2*v2ZP)/6. + (gp*PvL1hTo2*v2ZP)/6.
      MGVX280(2) = (gt*PvR1hTo2*v1ZP)/2. + (gp*PvR1hTo2*v2ZP)/6. + 
     +(2*gp*PvR2hTo2*v2ZP)/3.
      MGVX281(1) = (g0*v0ZP*vL0hTo*vL0to)/2. + (gp*v2ZP*vL0hTo*vL0t
     +o)/6. + (gt*v1ZP*vL1hTo*vL1to)/2. + (gp*v2ZP*vL1hTo*vL1to)/6
     +.
      MGVX281(2) = (gt*v1ZP*vR1hTo*vR1to)/2. + (gp*v2ZP*vR1hTo*vR1t
     +o)/6. + (2*gp*v2ZP*vR2hTo*vR2to)/3.
      MGVX283(1) = (g0*PvL0up2*v0ZP)/2. + (gt*PvL1up2*v1ZP)/2. + (g
     +p*PvL0up2*v2ZP)/6. + (gp*PvL1up2*v2ZP)/6.
      MGVX283(2) = (2*gp*v2ZP)/3.
      MGVX285(1) = (g0*PvL0ch2*v0ZP)/2. + (gt*PvL1ch2*v1ZP)/2. + (g
     +p*PvL0ch2*v2ZP)/6. + (gp*PvL1ch2*v2ZP)/6.
      MGVX285(2) = (gt*PvR1ch2*v1ZP)/2. + (gp*PvR1ch2*v2ZP)/6. + (2
     +*gp*PvR2ch2*v2ZP)/3.
      MGVX287(1) = (g0*PvL0to2*v0ZP)/2. + (gt*PvL1to2*v1ZP)/2. + (g
     +p*PvL0to2*v2ZP)/6. + (gp*PvL1to2*v2ZP)/6.
      MGVX287(2) = (gt*PvR1to2*v1ZP)/2. + (gp*PvR1to2*v2ZP)/6. + (2
     +*gp*PvR2to2*v2ZP)/3.
      endif

      return
      end

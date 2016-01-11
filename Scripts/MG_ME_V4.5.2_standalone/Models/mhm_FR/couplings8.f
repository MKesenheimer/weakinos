      subroutine coup8(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 8
      if(readlha) then


      MGVX265(1) = (g0*v0ZP*vL0hNu*vL0nu)/2. - (gp*v2ZP*vL0hNu*vL0n
     +u)/2. + (gt*v1ZP*vL1hNu*vL1nu)/2. - (gp*v2ZP*vL1hNu*vL1nu)/2
     +.
      MGVX265(2) = 0
      MGVX271(1) = (g0*PvL0nu2*v0ZP)/2. + (gt*PvL1nu2*v1ZP)/2. - (g
     +p*PvL0nu2*v2ZP)/2. - (gp*PvL1nu2*v2ZP)/2.
      MGVX271(2) = 0
      MGVX60(1) = -(EE*PvL0do2)/3. - (EE*PvL1do2)/3.
      MGVX60(2) = -EE/3.
      MGVX61(1) = -(EE*vL0do*vL0hDo)/3. - (EE*vL1do*vL1hDo)/3.
      MGVX61(2) = 0
      MGVX62(1) = -(EE*PvL0st2)/3. - (EE*PvL1st2)/3.
      MGVX62(2) = -(EE*PvR1st2)/3. - (EE*PvR2st2)/3.
      MGVX63(1) = -(EE*vL0hSt*vL0st)/3. - (EE*vL1hSt*vL1st)/3.
      MGVX63(2) = -(EE*vR1hSt*vR1st)/3. - (EE*vR2hSt*vR2st)/3.
      MGVX64(1) = -(EE*PvL0bo2)/3. - (EE*PvL1bo2)/3.
      MGVX64(2) = -(EE*PvR1bo2)/3. - (EE*PvR2bo2)/3.
      MGVX65(1) = -(EE*vL0bo*vL0hBo)/3. - (EE*vL1bo*vL1hBo)/3.
      MGVX65(2) = -(EE*vR1bo*vR1hBo)/3. - (EE*vR2bo*vR2hBo)/3.
      MGVX67(1) = -(EE*PvL0hDo2)/3. - (EE*PvL1hDo2)/3.
      MGVX67(2) = -EE/3.
      MGVX69(1) = -(EE*PvL0hSt2)/3. - (EE*PvL1hSt2)/3.
      MGVX69(2) = -(EE*PvR1hSt2)/3. - (EE*PvR2hSt2)/3.
      MGVX71(1) = -(EE*PvL0hBo2)/3. - (EE*PvL1hBo2)/3.
      MGVX71(2) = -(EE*PvR1hBo2)/3. - (EE*PvR2hBo2)/3.
      MGVX84(1) = (2*EE*PvL0hUp2)/3. + (2*EE*PvL1hUp2)/3.
      MGVX84(2) = (2*EE)/3.
      MGVX85(1) = (2*EE*vL0hUp*vL0up)/3. + (2*EE*vL1hUp*vL1up)/3.
      MGVX85(2) = 0
      MGVX86(1) = (2*EE*PvL0hCh2)/3. + (2*EE*PvL1hCh2)/3.
      MGVX86(2) = (2*EE*PvR1hCh2)/3. + (2*EE*PvR2hCh2)/3.
      MGVX87(1) = (2*EE*vL0ch*vL0hCh)/3. + (2*EE*vL1ch*vL1hCh)/3.
      MGVX87(2) = (2*EE*vR1ch*vR1hCh)/3. + (2*EE*vR2ch*vR2hCh)/3.
      endif

      return
      end

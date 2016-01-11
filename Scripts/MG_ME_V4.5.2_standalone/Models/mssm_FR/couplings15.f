      subroutine coup15(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 15
      if(readlha) then


      MGVX682(1) = (ee*OppL24)/(2.*cw*sw) - (ee*OppR42)/(2.*cw*sw)
      MGVX682(2) = -(ee*OppL42)/(2.*cw*sw) + (ee*OppR24)/(2.*cw*sw)
      MGVX683(1) = (ee*OppL33)/(2.*cw*sw) - (ee*OppR33)/(2.*cw*sw)
      MGVX683(2) = -(ee*OppL33)/(2.*cw*sw) + (ee*OppR33)/(2.*cw*sw)
      MGVX684(1) = (ee*OppL34)/(2.*cw*sw) - (ee*OppR43)/(2.*cw*sw)
      MGVX684(2) = -(ee*OppL43)/(2.*cw*sw) + (ee*OppR34)/(2.*cw*sw)
      MGVX685(1) = (ee*OppL44)/(2.*cw*sw) - (ee*OppR44)/(2.*cw*sw)
      MGVX685(2) = -(ee*OppL44)/(2.*cw*sw) + (ee*OppR44)/(2.*cw*sw)
      MGVX711(1) = -((Cosalp*yukl3)/Sqrt2)
      MGVX711(2) = -((Cosalp*yukl3)/Sqrt2)
      MGVX724(1) = (Sinalp*yukl3)/Sqrt2
      MGVX724(2) = (Sinalp*yukl3)/Sqrt2
      MGVX729(1) = ((0,-1)*Sinbet*yukl3)/Sqrt2
      MGVX729(2) = ((0,1)*Sinbet*yukl3)/Sqrt2
      MGVX730(1) = Sinbet*yukl3
      MGVX730(2) = 0
      MGVX731(1) = 0
      MGVX731(2) = Sinbet*yukl3
      MGVX9(1) = ee/3.
      MGVX9(2) = ee/3.
      MGVX15(1) = (-2*ee)/3.
      MGVX15(2) = (-2*ee)/3.
      MGVX700(1) = ((0,-1)*Sinbet*yukd3)/Sqrt2
      MGVX700(2) = ((0,1)*Sinbet*yukd3)/Sqrt2
      MGVX733(1) = ((0,-1)*Cosbet*yuku3)/Sqrt2
      MGVX733(2) = ((0,1)*Cosbet*yuku3)/Sqrt2
      MGVX404(1) = -(CONJNN11*CONJRd16*ee*Sqrt2)/(3.*cw) - CONJNN13
     +*CONJRd13*yukd3
      MGVX404(2) = -(CONJRd13*ee*NN11)/(3.*cw*Sqrt2) + (CONJRd13*ee
     +*NN12)/(Sqrt2*sw) - CONJRd16*NN13*yukd3
      MGVX405(1) = -(CONJNN21*CONJRd16*ee*Sqrt2)/(3.*cw) - CONJNN23
     +*CONJRd13*yukd3
      MGVX405(2) = -(CONJRd13*ee*NN21)/(3.*cw*Sqrt2) + (CONJRd13*ee
     +*NN22)/(Sqrt2*sw) - CONJRd16*NN23*yukd3
      endif

      return
      end

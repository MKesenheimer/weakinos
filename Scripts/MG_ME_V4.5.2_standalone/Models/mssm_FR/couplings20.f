      subroutine coup20(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 20
      if(readlha) then


      MGVX486(1) = -(CONJNN31*ee*Rd62)/(3.*cw*Sqrt2) + (CONJNN32*ee
     +*Rd62)/(Sqrt2*sw)
      MGVX486(2) = 0
      MGVX487(1) = -(CONJNN41*ee*Rd62)/(3.*cw*Sqrt2) + (CONJNN42*ee
     +*Rd62)/(Sqrt2*sw)
      MGVX487(2) = 0
      MGVX702(1) = Cosbet*yuku3
      MGVX702(2) = Sinbet*yukd3
      MGVX695(1) = (Sinalp*yukd3)/Sqrt2
      MGVX695(2) = (Sinalp*yukd3)/Sqrt2
      MGVX732(1) = -((Cosalp*yuku3)/Sqrt2)
      MGVX732(2) = -((Cosalp*yuku3)/Sqrt2)
      MGVX686(1) = -((Cosalp*yukd3)/Sqrt2)
      MGVX686(2) = -((Cosalp*yukd3)/Sqrt2)
      MGVX734(1) = -((Sinalp*yuku3)/Sqrt2)
      MGVX734(2) = -((Sinalp*yuku3)/Sqrt2)
      MGVX701(1) = Sinbet*yukd3
      MGVX701(2) = Cosbet*yuku3
      MGVX645(1) = CONJRu13*CONJUU12*yukd3
      MGVX645(2) = -((CONJRu13*ee*VV11)/sw) + CONJRu16*VV12*yuku3
      MGVX652(1) = CONJRu13*CONJUU22*yukd3
      MGVX652(2) = -((CONJRu13*ee*VV21)/sw) + CONJRu16*VV22*yuku3
      MGVX558(1) = (2*CONJNN11*CONJRu16*ee*Sqrt2)/(3.*cw) - CONJNN1
     +4*CONJRu13*yuku3
      MGVX558(2) = -(CONJRu13*ee*NN11)/(3.*cw*Sqrt2) - (CONJRu13*ee
     +*NN12)/(Sqrt2*sw) - CONJRu16*NN14*yuku3
      MGVX559(1) = (2*CONJNN21*CONJRu16*ee*Sqrt2)/(3.*cw) - CONJNN2
     +4*CONJRu13*yuku3
      MGVX559(2) = -(CONJRu13*ee*NN21)/(3.*cw*Sqrt2) - (CONJRu13*ee
     +*NN22)/(Sqrt2*sw) - CONJRu16*NN24*yuku3
      MGVX560(1) = (2*CONJNN31*CONJRu16*ee*Sqrt2)/(3.*cw) - CONJNN3
     +4*CONJRu13*yuku3
      MGVX560(2) = -(CONJRu13*ee*NN31)/(3.*cw*Sqrt2) - (CONJRu13*ee
     +*NN32)/(Sqrt2*sw) - CONJRu16*NN34*yuku3
      MGVX561(1) = (2*CONJNN41*CONJRu16*ee*Sqrt2)/(3.*cw) - CONJNN4
     +4*CONJRu13*yuku3
      MGVX561(2) = -(CONJRu13*ee*NN41)/(3.*cw*Sqrt2) - (CONJRu13*ee
     +*NN42)/(Sqrt2*sw) - CONJRu16*NN44*yuku3
      MGVX565(1) = -(CONJNN11*ee*Ru13)/(3.*cw*Sqrt2) - (CONJNN12*ee
     +*Ru13)/(Sqrt2*sw) - CONJNN14*Ru16*yuku3
      MGVX565(2) = (2*ee*NN11*Ru16*Sqrt2)/(3.*cw) - NN14*Ru13*yuku3
      endif

      return
      end

      subroutine coup16(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 16
      if(readlha) then


      MGVX406(1) = -(CONJNN31*CONJRd16*ee*Sqrt2)/(3.*cw) - CONJNN33
     +*CONJRd13*yukd3
      MGVX406(2) = -(CONJRd13*ee*NN31)/(3.*cw*Sqrt2) + (CONJRd13*ee
     +*NN32)/(Sqrt2*sw) - CONJRd16*NN33*yukd3
      MGVX407(1) = -(CONJNN41*CONJRd16*ee*Sqrt2)/(3.*cw) - CONJNN43
     +*CONJRd13*yukd3
      MGVX407(2) = -(CONJRd13*ee*NN41)/(3.*cw*Sqrt2) + (CONJRd13*ee
     +*NN42)/(Sqrt2*sw) - CONJRd16*NN43*yukd3
      MGVX630(1) = CONJRd13*CONJVV12*yuku3
      MGVX630(2) = -((CONJRd13*ee*UU11)/sw) + CONJRd16*UU12*yukd3
      MGVX638(1) = CONJRd13*CONJVV22*yuku3
      MGVX638(2) = -((CONJRd13*ee*UU21)/sw) + CONJRd16*UU22*yukd3
      MGVX456(1) = -((CONJUU11*ee*Rd13)/sw) + CONJUU12*Rd16*yukd3
      MGVX456(2) = Rd13*VV12*yuku3
      MGVX457(1) = -((CONJUU21*ee*Rd13)/sw) + CONJUU22*Rd16*yukd3
      MGVX457(2) = Rd13*VV22*yuku3
      MGVX458(1) = -(CONJNN11*ee*Rd13)/(3.*cw*Sqrt2) + (CONJNN12*ee
     +*Rd13)/(Sqrt2*sw) - CONJNN13*Rd16*yukd3
      MGVX458(2) = -(ee*NN11*Rd16*Sqrt2)/(3.*cw) - NN13*Rd13*yukd3
      MGVX459(1) = -(CONJNN21*ee*Rd13)/(3.*cw*Sqrt2) + (CONJNN22*ee
     +*Rd13)/(Sqrt2*sw) - CONJNN23*Rd16*yukd3
      MGVX459(2) = -(ee*NN21*Rd16*Sqrt2)/(3.*cw) - NN23*Rd13*yukd3
      MGVX460(1) = -(CONJNN31*ee*Rd13)/(3.*cw*Sqrt2) + (CONJNN32*ee
     +*Rd13)/(Sqrt2*sw) - CONJNN33*Rd16*yukd3
      MGVX460(2) = -(ee*NN31*Rd16*Sqrt2)/(3.*cw) - NN33*Rd13*yukd3
      MGVX461(1) = -(CONJNN41*ee*Rd13)/(3.*cw*Sqrt2) + (CONJNN42*ee
     +*Rd13)/(Sqrt2*sw) - CONJNN43*Rd16*yukd3
      MGVX461(2) = -(ee*NN41*Rd16*Sqrt2)/(3.*cw) - NN43*Rd13*yukd3
      MGVX408(1) = -(CONJNN11*CONJRd26*ee*Sqrt2)/(3.*cw) - CONJNN13
     +*CONJRd23*yukd3
      MGVX408(2) = -(CONJRd23*ee*NN11)/(3.*cw*Sqrt2) + (CONJRd23*ee
     +*NN12)/(Sqrt2*sw) - CONJRd26*NN13*yukd3
      MGVX409(1) = -(CONJNN21*CONJRd26*ee*Sqrt2)/(3.*cw) - CONJNN23
     +*CONJRd23*yukd3
      MGVX409(2) = -(CONJRd23*ee*NN21)/(3.*cw*Sqrt2) + (CONJRd23*ee
     +*NN22)/(Sqrt2*sw) - CONJRd26*NN23*yukd3
      MGVX410(1) = -(CONJNN31*CONJRd26*ee*Sqrt2)/(3.*cw) - CONJNN33
     +*CONJRd23*yukd3
      MGVX410(2) = -(CONJRd23*ee*NN31)/(3.*cw*Sqrt2) + (CONJRd23*ee
     +*NN32)/(Sqrt2*sw) - CONJRd26*NN33*yukd3
      MGVX411(1) = -(CONJNN41*CONJRd26*ee*Sqrt2)/(3.*cw) - CONJNN43
     +*CONJRd23*yukd3
      MGVX411(2) = -(CONJRd23*ee*NN41)/(3.*cw*Sqrt2) + (CONJRd23*ee
     +*NN42)/(Sqrt2*sw) - CONJRd26*NN43*yukd3
      MGVX631(1) = CONJRd23*CONJVV12*yuku3
      MGVX631(2) = -((CONJRd23*ee*UU11)/sw) + CONJRd26*UU12*yukd3
      endif

      return
      end

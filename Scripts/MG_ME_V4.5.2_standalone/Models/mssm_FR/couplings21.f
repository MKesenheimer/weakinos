      subroutine coup21(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 21
      if(readlha) then


      MGVX566(1) = -(CONJNN21*ee*Ru13)/(3.*cw*Sqrt2) - (CONJNN22*ee
     +*Ru13)/(Sqrt2*sw) - CONJNN24*Ru16*yuku3
      MGVX566(2) = (2*ee*NN21*Ru16*Sqrt2)/(3.*cw) - NN24*Ru13*yuku3
      MGVX567(1) = -(CONJNN31*ee*Ru13)/(3.*cw*Sqrt2) - (CONJNN32*ee
     +*Ru13)/(Sqrt2*sw) - CONJNN34*Ru16*yuku3
      MGVX567(2) = (2*ee*NN31*Ru16*Sqrt2)/(3.*cw) - NN34*Ru13*yuku3
      MGVX568(1) = -(CONJNN41*ee*Ru13)/(3.*cw*Sqrt2) - (CONJNN42*ee
     +*Ru13)/(Sqrt2*sw) - CONJNN44*Ru16*yuku3
      MGVX568(2) = (2*ee*NN41*Ru16*Sqrt2)/(3.*cw) - NN44*Ru13*yuku3
      MGVX563(1) = -((CONJVV11*ee*Ru13)/sw) + CONJVV12*Ru16*yuku3
      MGVX563(2) = Ru13*UU12*yukd3
      MGVX564(1) = -((CONJVV21*ee*Ru13)/sw) + CONJVV22*Ru16*yuku3
      MGVX564(2) = Ru13*UU22*yukd3
      MGVX570(1) = (2*CONJNN11*CONJRu24*ee*Sqrt2)/(3.*cw)
      MGVX570(2) = 0
      MGVX571(1) = (2*CONJNN21*CONJRu24*ee*Sqrt2)/(3.*cw)
      MGVX571(2) = 0
      MGVX572(1) = (2*CONJNN31*CONJRu24*ee*Sqrt2)/(3.*cw)
      MGVX572(2) = 0
      MGVX573(1) = (2*CONJNN41*CONJRu24*ee*Sqrt2)/(3.*cw)
      MGVX573(2) = 0
      MGVX575(1) = 0
      MGVX575(2) = (2*ee*NN11*Ru24*Sqrt2)/(3.*cw)
      MGVX576(1) = 0
      MGVX576(2) = (2*ee*NN21*Ru24*Sqrt2)/(3.*cw)
      MGVX577(1) = 0
      MGVX577(2) = (2*ee*NN31*Ru24*Sqrt2)/(3.*cw)
      MGVX578(1) = 0
      MGVX578(2) = (2*ee*NN41*Ru24*Sqrt2)/(3.*cw)
      MGVX580(1) = (2*CONJNN11*CONJRu35*ee*Sqrt2)/(3.*cw)
      MGVX580(2) = 0
      MGVX581(1) = (2*CONJNN21*CONJRu35*ee*Sqrt2)/(3.*cw)
      MGVX581(2) = 0
      endif

      return
      end

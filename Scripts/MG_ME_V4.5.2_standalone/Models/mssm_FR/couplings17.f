      subroutine coup17(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 17
      if(readlha) then


      MGVX639(1) = CONJRd23*CONJVV22*yuku3
      MGVX639(2) = -((CONJRd23*ee*UU21)/sw) + CONJRd26*UU22*yukd3
      MGVX462(1) = -((CONJUU11*ee*Rd23)/sw) + CONJUU12*Rd26*yukd3
      MGVX462(2) = Rd23*VV12*yuku3
      MGVX463(1) = -((CONJUU21*ee*Rd23)/sw) + CONJUU22*Rd26*yukd3
      MGVX463(2) = Rd23*VV22*yuku3
      MGVX464(1) = -(CONJNN11*ee*Rd23)/(3.*cw*Sqrt2) + (CONJNN12*ee
     +*Rd23)/(Sqrt2*sw) - CONJNN13*Rd26*yukd3
      MGVX464(2) = -(ee*NN11*Rd26*Sqrt2)/(3.*cw) - NN13*Rd23*yukd3
      MGVX465(1) = -(CONJNN21*ee*Rd23)/(3.*cw*Sqrt2) + (CONJNN22*ee
     +*Rd23)/(Sqrt2*sw) - CONJNN23*Rd26*yukd3
      MGVX465(2) = -(ee*NN21*Rd26*Sqrt2)/(3.*cw) - NN23*Rd23*yukd3
      MGVX466(1) = -(CONJNN31*ee*Rd23)/(3.*cw*Sqrt2) + (CONJNN32*ee
     +*Rd23)/(Sqrt2*sw) - CONJNN33*Rd26*yukd3
      MGVX466(2) = -(ee*NN31*Rd26*Sqrt2)/(3.*cw) - NN33*Rd23*yukd3
      MGVX467(1) = -(CONJNN41*ee*Rd23)/(3.*cw*Sqrt2) + (CONJNN42*ee
     +*Rd23)/(Sqrt2*sw) - CONJNN43*Rd26*yukd3
      MGVX467(2) = -(ee*NN41*Rd26*Sqrt2)/(3.*cw) - NN43*Rd23*yukd3
      MGVX412(1) = -(CONJNN11*CONJRd35*ee*Sqrt2)/(3.*cw)
      MGVX412(2) = 0
      MGVX413(1) = -(CONJNN21*CONJRd35*ee*Sqrt2)/(3.*cw)
      MGVX413(2) = 0
      MGVX414(1) = -(CONJNN31*CONJRd35*ee*Sqrt2)/(3.*cw)
      MGVX414(2) = 0
      MGVX415(1) = -(CONJNN41*CONJRd35*ee*Sqrt2)/(3.*cw)
      MGVX415(2) = 0
      MGVX468(1) = 0
      MGVX468(2) = -(ee*NN11*Rd35*Sqrt2)/(3.*cw)
      MGVX469(1) = 0
      MGVX469(2) = -(ee*NN21*Rd35*Sqrt2)/(3.*cw)
      MGVX470(1) = 0
      MGVX470(2) = -(ee*NN31*Rd35*Sqrt2)/(3.*cw)
      MGVX471(1) = 0
      MGVX471(2) = -(ee*NN41*Rd35*Sqrt2)/(3.*cw)
      endif

      return
      end

      subroutine coup2(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 2
      if(readlha) then
      MGVX3 = -(gw*sw)
      MGVX5 = cw*gw
      MGVX4 = Pgw2*Psw2
      MGVX7 = cw*Pgw2*sw
      MGVX6 = -Pgw2
      MGVX8 = Pcw2*Pgw2


      MGVX18(1) = ee
      MGVX18(2) = ee
      MGVX24(1) = -(ee/(Sqrt2*sw))
      MGVX24(2) = 0
      MGVX43(1) = (cw*ee)/(2.*sw) - (ee*sw)/(2.*cw)
      MGVX43(2) = -((ee*sw)/cw)
      MGVX49(1) = -(cw*ee)/(2.*sw) - (ee*sw)/(2.*cw)
      MGVX49(2) = 0
      MGVX141(1) = -((CONJGL11*Sqrt2)/v)
      MGVX141(2) = 0
      MGVX142(1) = -((CONJGL21*Sqrt2)/v)
      MGVX142(2) = 0
      MGVX143(1) = -((CONJGL31*Sqrt2)/v)
      MGVX143(2) = 0
      MGVX144(1) = -((CONJGL12*Sqrt2)/v)
      MGVX144(2) = 0
      MGVX145(1) = -((CONJGL22*Sqrt2)/v)
      MGVX145(2) = 0
      MGVX146(1) = -((CONJGL32*Sqrt2)/v)
      MGVX146(2) = 0
      MGVX147(1) = -((CONJGL13*Sqrt2)/v)
      MGVX147(2) = 0
      MGVX148(1) = -((CONJGL23*Sqrt2)/v)
      MGVX148(2) = 0
      MGVX149(1) = -((CONJGL33*Sqrt2)/v)
      MGVX149(2) = 0
      MGVX159(1) = 0
      MGVX159(2) = -((GL11*Sqrt2)/v)
      MGVX160(1) = 0
      MGVX160(2) = -((GL12*Sqrt2)/v)
      MGVX161(1) = 0
      MGVX161(2) = -((GL13*Sqrt2)/v)
      MGVX162(1) = 0
      MGVX162(2) = -((GL21*Sqrt2)/v)
      MGVX163(1) = 0
      MGVX163(2) = -((GL22*Sqrt2)/v)
      MGVX164(1) = 0
      MGVX164(2) = -((GL23*Sqrt2)/v)
      MGVX165(1) = 0
      MGVX165(2) = -((GL31*Sqrt2)/v)
      MGVX166(1) = 0
      MGVX166(2) = -((GL32*Sqrt2)/v)
      MGVX167(1) = 0
      MGVX167(2) = -((GL33*Sqrt2)/v)
      MGVX177(1) = -((CONJGL11*TH21)/v) + ((0,1)*CONJGL11*TH31)/v -
     + (TH11*yukl1)/v
      MGVX177(2) = -((GL11*TH21)/v) - ((0,1)*GL11*TH31)/v - (TH11*y
     +ukl1)/v
      MGVX178(1) = -((CONJGL21*TH21)/v) + ((0,1)*CONJGL21*TH31)/v
      MGVX178(2) = -((GL12*TH21)/v) - ((0,1)*GL12*TH31)/v
      MGVX179(1) = -((CONJGL31*TH21)/v) + ((0,1)*CONJGL31*TH31)/v
      MGVX179(2) = -((GL13*TH21)/v) - ((0,1)*GL13*TH31)/v
      MGVX180(1) = -((CONJGL12*TH21)/v) + ((0,1)*CONJGL12*TH31)/v
      MGVX180(2) = -((GL21*TH21)/v) - ((0,1)*GL21*TH31)/v
      MGVX181(1) = -((CONJGL22*TH21)/v) + ((0,1)*CONJGL22*TH31)/v -
     + (TH11*yukl2)/v
      MGVX181(2) = -((GL22*TH21)/v) - ((0,1)*GL22*TH31)/v - (TH11*y
     +ukl2)/v
      MGVX182(1) = -((CONJGL32*TH21)/v) + ((0,1)*CONJGL32*TH31)/v
      MGVX182(2) = -((GL23*TH21)/v) - ((0,1)*GL23*TH31)/v
      MGVX183(1) = -((CONJGL13*TH21)/v) + ((0,1)*CONJGL13*TH31)/v
      MGVX183(2) = -((GL31*TH21)/v) - ((0,1)*GL31*TH31)/v
      MGVX184(1) = -((CONJGL23*TH21)/v) + ((0,1)*CONJGL23*TH31)/v
      MGVX184(2) = -((GL32*TH21)/v) - ((0,1)*GL32*TH31)/v
      MGVX185(1) = -((CONJGL33*TH21)/v) + ((0,1)*CONJGL33*TH31)/v -
     + (TH11*yukl3)/v
      MGVX185(2) = -((GL33*TH21)/v) - ((0,1)*GL33*TH31)/v - (TH11*y
     +ukl3)/v
      MGVX204(1) = -((CONJGL11*TH22)/v) + ((0,1)*CONJGL11*TH32)/v -
     + (TH12*yukl1)/v
      MGVX204(2) = -((GL11*TH22)/v) - ((0,1)*GL11*TH32)/v - (TH12*y
     +ukl1)/v
      MGVX205(1) = -((CONJGL21*TH22)/v) + ((0,1)*CONJGL21*TH32)/v
      MGVX205(2) = -((GL12*TH22)/v) - ((0,1)*GL12*TH32)/v
      MGVX206(1) = -((CONJGL31*TH22)/v) + ((0,1)*CONJGL31*TH32)/v
      MGVX206(2) = -((GL13*TH22)/v) - ((0,1)*GL13*TH32)/v
      MGVX207(1) = -((CONJGL12*TH22)/v) + ((0,1)*CONJGL12*TH32)/v
      MGVX207(2) = -((GL21*TH22)/v) - ((0,1)*GL21*TH32)/v
      MGVX208(1) = -((CONJGL22*TH22)/v) + ((0,1)*CONJGL22*TH32)/v -
     + (TH12*yukl2)/v
      MGVX208(2) = -((GL22*TH22)/v) - ((0,1)*GL22*TH32)/v - (TH12*y
     +ukl2)/v
      MGVX209(1) = -((CONJGL32*TH22)/v) + ((0,1)*CONJGL32*TH32)/v
      MGVX209(2) = -((GL23*TH22)/v) - ((0,1)*GL23*TH32)/v
      MGVX210(1) = -((CONJGL13*TH22)/v) + ((0,1)*CONJGL13*TH32)/v
      MGVX210(2) = -((GL31*TH22)/v) - ((0,1)*GL31*TH32)/v
      MGVX211(1) = -((CONJGL23*TH22)/v) + ((0,1)*CONJGL23*TH32)/v
      MGVX211(2) = -((GL32*TH22)/v) - ((0,1)*GL32*TH32)/v
      MGVX212(1) = -((CONJGL33*TH22)/v) + ((0,1)*CONJGL33*TH32)/v -
     + (TH12*yukl3)/v
      MGVX212(2) = -((GL33*TH22)/v) - ((0,1)*GL33*TH32)/v - (TH12*y
     +ukl3)/v
      MGVX231(1) = -((CONJGL11*TH23)/v) + ((0,1)*CONJGL11*TH33)/v -
     + (TH13*yukl1)/v
      MGVX231(2) = -((GL11*TH23)/v) - ((0,1)*GL11*TH33)/v - (TH13*y
     +ukl1)/v
      MGVX232(1) = -((CONJGL21*TH23)/v) + ((0,1)*CONJGL21*TH33)/v
      MGVX232(2) = -((GL12*TH23)/v) - ((0,1)*GL12*TH33)/v
      MGVX233(1) = -((CONJGL31*TH23)/v) + ((0,1)*CONJGL31*TH33)/v
      MGVX233(2) = -((GL13*TH23)/v) - ((0,1)*GL13*TH33)/v
      MGVX234(1) = -((CONJGL12*TH23)/v) + ((0,1)*CONJGL12*TH33)/v
      MGVX234(2) = -((GL21*TH23)/v) - ((0,1)*GL21*TH33)/v
      endif

      return
      end

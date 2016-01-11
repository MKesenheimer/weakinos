      subroutine coup4(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 4
      if(readlha) then

      MGVX57 = -6*l1*PTH113*v - (3*CONJl7*PTH213*v)/2. - (3*l7*PTH2
     +13*v)/2. + (0,1.5)*CONJl7*PTH313*v - (0,1.5)*l7*PTH313*v - 3
     +*l3*PTH212*TH11*v - 3*l4*PTH212*TH11*v - 6*l5*PTH212*TH11*v 
     +- 3*l3*PTH312*TH11*v - 3*l4*PTH312*TH11*v + 6*l5*PTH312*TH11
     +*v - (9*CONJl6*PTH112*TH21*v)/2. - (9*l6*PTH112*TH21*v)/2. -
     + (3*CONJl7*PTH312*TH21*v)/2. - (3*l7*PTH312*TH21*v)/2. + (0,
     +4.5)*CONJl6*PTH112*TH31*v - (0,4.5)*l6*PTH112*TH31*v + (0,1.
     +5)*CONJl7*PTH212*TH31*v - (0,1.5)*l7*PTH212*TH31*v
      MGVX61 = -6*l1*PTH112*TH12*v - l3*PTH212*TH12*v - l4*PTH212*T
     +H12*v - 2*l5*PTH212*TH12*v - l3*PTH312*TH12*v - l4*PTH312*TH
     +12*v + 2*l5*PTH312*TH12*v - 3*CONJl6*TH11*TH12*TH21*v - 3*l6
     +*TH11*TH12*TH21*v - (3*CONJl6*PTH112*TH22*v)/2. - (3*l6*PTH1
     +12*TH22*v)/2. - (3*CONJl7*PTH212*TH22*v)/2. - (3*l7*PTH212*T
     +H22*v)/2. - (CONJl7*PTH312*TH22*v)/2. - (l7*PTH312*TH22*v)/2
     +. - 2*l3*TH11*TH21*TH22*v - 2*l4*TH11*TH21*TH22*v - 4*l5*TH1
     +1*TH21*TH22*v + (0,3)*CONJl6*TH11*TH12*TH31*v - (0,3)*l6*TH1
     +1*TH12*TH31*v + (0,1)*CONJl7*TH21*TH22*TH31*v - (0,1)*l7*TH2
     +1*TH22*TH31*v + (0,1.5)*CONJl6*PTH112*TH32*v - (0,1.5)*l6*PT
     +H112*TH32*v + (0,0.5)*CONJl7*PTH212*TH32*v - (0,0.5)*l7*PTH2
     +12*TH32*v + (0,1.5)*CONJl7*PTH312*TH32*v - (0,1.5)*l7*PTH312
     +*TH32*v - 2*l3*TH11*TH31*TH32*v - 2*l4*TH11*TH31*TH32*v + 4*
     +l5*TH11*TH31*TH32*v - CONJl7*TH21*TH31*TH32*v - l7*TH21*TH31
     +*TH32*v
      MGVX71 = -6*l1*PTH112*TH13*v - l3*PTH212*TH13*v - l4*PTH212*T
     +H13*v - 2*l5*PTH212*TH13*v - l3*PTH312*TH13*v - l4*PTH312*TH
     +13*v + 2*l5*PTH312*TH13*v - 3*CONJl6*TH11*TH13*TH21*v - 3*l6
     +*TH11*TH13*TH21*v - (3*CONJl6*PTH112*TH23*v)/2. - (3*l6*PTH1
     +12*TH23*v)/2. - (3*CONJl7*PTH212*TH23*v)/2. - (3*l7*PTH212*T
     +H23*v)/2. - (CONJl7*PTH312*TH23*v)/2. - (l7*PTH312*TH23*v)/2
     +. - 2*l3*TH11*TH21*TH23*v - 2*l4*TH11*TH21*TH23*v - 4*l5*TH1
     +1*TH21*TH23*v + (0,3)*CONJl6*TH11*TH13*TH31*v - (0,3)*l6*TH1
     +1*TH13*TH31*v + (0,1)*CONJl7*TH21*TH23*TH31*v - (0,1)*l7*TH2
     +1*TH23*TH31*v + (0,1.5)*CONJl6*PTH112*TH33*v - (0,1.5)*l6*PT
     +H112*TH33*v + (0,0.5)*CONJl7*PTH212*TH33*v - (0,0.5)*l7*PTH2
     +12*TH33*v + (0,1.5)*CONJl7*PTH312*TH33*v - (0,1.5)*l7*PTH312
     +*TH33*v - 2*l3*TH11*TH31*TH33*v - 2*l4*TH11*TH31*TH33*v + 4*
     +l5*TH11*TH31*TH33*v - CONJl7*TH21*TH31*TH33*v - l7*TH21*TH31
     +*TH33*v
      MGVX64 = -6*l1*PTH122*TH11*v - l3*PTH222*TH11*v - l4*PTH222*T
     +H11*v - 2*l5*PTH222*TH11*v - l3*PTH322*TH11*v - l4*PTH322*TH
     +11*v + 2*l5*PTH322*TH11*v - (3*CONJl6*PTH122*TH21*v)/2. - (3
     +*l6*PTH122*TH21*v)/2. - (3*CONJl7*PTH222*TH21*v)/2. - (3*l7*
     +PTH222*TH21*v)/2. - (CONJl7*PTH322*TH21*v)/2. - (l7*PTH322*T
     +H21*v)/2. - 3*CONJl6*TH11*TH12*TH22*v - 3*l6*TH11*TH12*TH22*
     +v - 2*l3*TH12*TH21*TH22*v - 2*l4*TH12*TH21*TH22*v - 4*l5*TH1
     +2*TH21*TH22*v + (0,1.5)*CONJl6*PTH122*TH31*v - (0,1.5)*l6*PT
     +H122*TH31*v + (0,0.5)*CONJl7*PTH222*TH31*v - (0,0.5)*l7*PTH2
     +22*TH31*v + (0,1.5)*CONJl7*PTH322*TH31*v - (0,1.5)*l7*PTH322
     +*TH31*v + (0,3)*CONJl6*TH11*TH12*TH32*v - (0,3)*l6*TH11*TH12
     +*TH32*v + (0,1)*CONJl7*TH21*TH22*TH32*v - (0,1)*l7*TH21*TH22
     +*TH32*v - 2*l3*TH12*TH31*TH32*v - 2*l4*TH12*TH31*TH32*v + 4*
     +l5*TH12*TH31*TH32*v - CONJl7*TH22*TH31*TH32*v - l7*TH22*TH31
     +*TH32*v
      MGVX74 = -6*l1*TH11*TH12*TH13*v - (3*CONJl6*TH12*TH13*TH21*v)
     +/2. - (3*l6*TH12*TH13*TH21*v)/2. - (3*CONJl6*TH11*TH13*TH22*
     +v)/2. - (3*l6*TH11*TH13*TH22*v)/2. - l3*TH13*TH21*TH22*v - l
     +4*TH13*TH21*TH22*v - 2*l5*TH13*TH21*TH22*v - (3*CONJl6*TH11*
     +TH12*TH23*v)/2. - (3*l6*TH11*TH12*TH23*v)/2. - l3*TH12*TH21*
     +TH23*v - l4*TH12*TH21*TH23*v - 2*l5*TH12*TH21*TH23*v - l3*TH
     +11*TH22*TH23*v - l4*TH11*TH22*TH23*v - 2*l5*TH11*TH22*TH23*v
     + - (3*CONJl7*TH21*TH22*TH23*v)/2. - (3*l7*TH21*TH22*TH23*v)/
     +2. + (0,1.5)*CONJl6*TH12*TH13*TH31*v - (0,1.5)*l6*TH12*TH13*
     +TH31*v + (0,0.5)*CONJl7*TH22*TH23*TH31*v - (0,0.5)*l7*TH22*T
     +H23*TH31*v + (0,1.5)*CONJl6*TH11*TH13*TH32*v - (0,1.5)*l6*TH
     +11*TH13*TH32*v + (0,0.5)*CONJl7*TH21*TH23*TH32*v - (0,0.5)*l
     +7*TH21*TH23*TH32*v - l3*TH13*TH31*TH32*v - l4*TH13*TH31*TH32
     +*v + 2*l5*TH13*TH31*TH32*v - (CONJl7*TH23*TH31*TH32*v)/2. - 
     +(l7*TH23*TH31*TH32*v)/2. + (0,1.5)*CONJl6*TH11*TH12*TH33*v -
     + (0,1.5)*l6*TH11*TH12*TH33*v + (0,0.5)*CONJl7*TH21*TH22*TH33
     +*v - (0,0.5)*l7*TH21*TH22*TH33*v - l3*TH12*TH31*TH33*v - l4*
     +TH12*TH31*TH33*v + 2*l5*TH12*TH31*TH33*v - (CONJl7*TH22*TH31
     +*TH33*v)/2. - (l7*TH22*TH31*TH33*v)/2. - l3*TH11*TH32*TH33*v
     + - l4*TH11*TH32*TH33*v + 2*l5*TH11*TH32*TH33*v - (CONJl7*TH2
     +1*TH32*TH33*v)/2. - (l7*TH21*TH32*TH33*v)/2. + (0,1.5)*CONJl
     +7*TH31*TH32*TH33*v - (0,1.5)*l7*TH31*TH32*TH33*v
      MGVX80 = -6*l1*PTH132*TH11*v - l3*PTH232*TH11*v - l4*PTH232*T
     +H11*v - 2*l5*PTH232*TH11*v - l3*PTH332*TH11*v - l4*PTH332*TH
     +11*v + 2*l5*PTH332*TH11*v - (3*CONJl6*PTH132*TH21*v)/2. - (3
     +*l6*PTH132*TH21*v)/2. - (3*CONJl7*PTH232*TH21*v)/2. - (3*l7*
     +PTH232*TH21*v)/2. - (CONJl7*PTH332*TH21*v)/2. - (l7*PTH332*T
     +H21*v)/2. - 3*CONJl6*TH11*TH13*TH23*v - 3*l6*TH11*TH13*TH23*
     +v - 2*l3*TH13*TH21*TH23*v - 2*l4*TH13*TH21*TH23*v - 4*l5*TH1
     +3*TH21*TH23*v + (0,1.5)*CONJl6*PTH132*TH31*v - (0,1.5)*l6*PT
     +H132*TH31*v + (0,0.5)*CONJl7*PTH232*TH31*v - (0,0.5)*l7*PTH2
     +32*TH31*v + (0,1.5)*CONJl7*PTH332*TH31*v - (0,1.5)*l7*PTH332
     +*TH31*v + (0,3)*CONJl6*TH11*TH13*TH33*v - (0,3)*l6*TH11*TH13
     +*TH33*v + (0,1)*CONJl7*TH21*TH23*TH33*v - (0,1)*l7*TH21*TH23
     +*TH33*v - 2*l3*TH13*TH31*TH33*v - 2*l4*TH13*TH31*TH33*v + 4*
     +l5*TH13*TH31*TH33*v - CONJl7*TH23*TH31*TH33*v - l7*TH23*TH31
     +*TH33*v
      MGVX55 = -(l3*TH11*v) - (CONJl7*TH21*v)/2. - (l7*TH21*v)/2. +
     + (0,0.5)*CONJl7*TH31*v - (0,0.5)*l7*TH31*v
      MGVX66 = -6*l1*PTH123*v - (3*CONJl7*PTH223*v)/2. - (3*l7*PTH2
     +23*v)/2. + (0,1.5)*CONJl7*PTH323*v - (0,1.5)*l7*PTH323*v - 3
     +*l3*PTH222*TH12*v - 3*l4*PTH222*TH12*v - 6*l5*PTH222*TH12*v 
     +- 3*l3*PTH322*TH12*v - 3*l4*PTH322*TH12*v + 6*l5*PTH322*TH12
     +*v - (9*CONJl6*PTH122*TH22*v)/2. - (9*l6*PTH122*TH22*v)/2. -
     + (3*CONJl7*PTH322*TH22*v)/2. - (3*l7*PTH322*TH22*v)/2. + (0,
     +4.5)*CONJl6*PTH122*TH32*v - (0,4.5)*l6*PTH122*TH32*v + (0,1.
     +5)*CONJl7*PTH222*TH32*v - (0,1.5)*l7*PTH222*TH32*v
      MGVX76 = -6*l1*PTH122*TH13*v - l3*PTH222*TH13*v - l4*PTH222*T
     +H13*v - 2*l5*PTH222*TH13*v - l3*PTH322*TH13*v - l4*PTH322*TH
     +13*v + 2*l5*PTH322*TH13*v - 3*CONJl6*TH12*TH13*TH22*v - 3*l6
     +*TH12*TH13*TH22*v - (3*CONJl6*PTH122*TH23*v)/2. - (3*l6*PTH1
     +22*TH23*v)/2. - (3*CONJl7*PTH222*TH23*v)/2. - (3*l7*PTH222*T
     +H23*v)/2. - (CONJl7*PTH322*TH23*v)/2. - (l7*PTH322*TH23*v)/2
     +. - 2*l3*TH12*TH22*TH23*v - 2*l4*TH12*TH22*TH23*v - 4*l5*TH1
     +2*TH22*TH23*v + (0,3)*CONJl6*TH12*TH13*TH32*v - (0,3)*l6*TH1
     +2*TH13*TH32*v + (0,1)*CONJl7*TH22*TH23*TH32*v - (0,1)*l7*TH2
     +2*TH23*TH32*v + (0,1.5)*CONJl6*PTH122*TH33*v - (0,1.5)*l6*PT
     +H122*TH33*v + (0,0.5)*CONJl7*PTH222*TH33*v - (0,0.5)*l7*PTH2
     +22*TH33*v + (0,1.5)*CONJl7*PTH322*TH33*v - (0,1.5)*l7*PTH322
     +*TH33*v - 2*l3*TH12*TH32*TH33*v - 2*l4*TH12*TH32*TH33*v + 4*
     +l5*TH12*TH32*TH33*v - CONJl7*TH22*TH32*TH33*v - l7*TH22*TH32
     +*TH33*v
      MGVX82 = -6*l1*PTH132*TH12*v - l3*PTH232*TH12*v - l4*PTH232*T
     +H12*v - 2*l5*PTH232*TH12*v - l3*PTH332*TH12*v - l4*PTH332*TH
     +12*v + 2*l5*PTH332*TH12*v - (3*CONJl6*PTH132*TH22*v)/2. - (3
     +*l6*PTH132*TH22*v)/2. - (3*CONJl7*PTH232*TH22*v)/2. - (3*l7*
     +PTH232*TH22*v)/2. - (CONJl7*PTH332*TH22*v)/2. - (l7*PTH332*T
     +H22*v)/2. - 3*CONJl6*TH12*TH13*TH23*v - 3*l6*TH12*TH13*TH23*
     +v - 2*l3*TH13*TH22*TH23*v - 2*l4*TH13*TH22*TH23*v - 4*l5*TH1
     +3*TH22*TH23*v + (0,1.5)*CONJl6*PTH132*TH32*v - (0,1.5)*l6*PT
     +H132*TH32*v + (0,0.5)*CONJl7*PTH232*TH32*v - (0,0.5)*l7*PTH2
     +32*TH32*v + (0,1.5)*CONJl7*PTH332*TH32*v - (0,1.5)*l7*PTH332
     +*TH32*v + (0,3)*CONJl6*TH12*TH13*TH33*v - (0,3)*l6*TH12*TH13
     +*TH33*v + (0,1)*CONJl7*TH22*TH23*TH33*v - (0,1)*l7*TH22*TH23
     +*TH33*v - 2*l3*TH13*TH32*TH33*v - 2*l4*TH13*TH32*TH33*v + 4*
     +l5*TH13*TH32*TH33*v - CONJl7*TH23*TH32*TH33*v - l7*TH23*TH32
     +*TH33*v
      MGVX59 = -(l3*TH12*v) - (CONJl7*TH22*v)/2. - (l7*TH22*v)/2. +
     + (0,0.5)*CONJl7*TH32*v - (0,0.5)*l7*TH32*v

      MGVX226(1) = -((CONJGD11*TH23)/v) + ((0,1)*CONJGD11*TH33)/v -
     + (TH13*yukd1)/v
      MGVX226(2) = -((GD11*TH23)/v) - ((0,1)*GD11*TH33)/v - (TH13*y
     +ukd1)/v
      MGVX227(1) = -((CONJGD21*TH23)/v) + ((0,1)*CONJGD21*TH33)/v
      MGVX227(2) = -((GD12*TH23)/v) - ((0,1)*GD12*TH33)/v
      MGVX228(1) = -((CONJGD32*TH23)/v) + ((0,1)*CONJGD32*TH33)/v
      MGVX228(2) = -((GD23*TH23)/v) - ((0,1)*GD23*TH33)/v
      MGVX229(1) = -((CONJGD12*TH23)/v) + ((0,1)*CONJGD12*TH33)/v
      MGVX229(2) = -((GD21*TH23)/v) - ((0,1)*GD21*TH33)/v
      MGVX230(1) = -((CONJGD22*TH23)/v) + ((0,1)*CONJGD22*TH33)/v -
     + (TH13*yukd2)/v
      MGVX230(2) = -((GD22*TH23)/v) - ((0,1)*GD22*TH33)/v - (TH13*y
     +ukd2)/v
      MGVX243(1) = -((CONJGU23*TH23)/v) - ((0,1)*CONJGU23*TH33)/v
      MGVX243(2) = -((GU32*TH23)/v) + ((0,1)*GU32*TH33)/v
      MGVX244(1) = -((CONJGU33*TH23)/v) - ((0,1)*CONJGU33*TH33)/v -
     + (TH13*yuku3)/v
      MGVX244(2) = -((GU33*TH23)/v) + ((0,1)*GU33*TH33)/v - (TH13*y
     +uku3)/v
      MGVX245(1) = -((CONJGU13*TH23)/v) - ((0,1)*CONJGU13*TH33)/v
      MGVX245(2) = -((GU31*TH23)/v) + ((0,1)*GU31*TH33)/v
      MGVX246(1) = -((CONJGU21*TH23)/v) - ((0,1)*CONJGU21*TH33)/v
      MGVX246(2) = -((GU12*TH23)/v) + ((0,1)*GU12*TH33)/v
      MGVX247(1) = -((CONJGU31*TH23)/v) - ((0,1)*CONJGU31*TH33)/v
      MGVX247(2) = -((GU13*TH23)/v) + ((0,1)*GU13*TH33)/v
      MGVX248(1) = -((CONJGU11*TH23)/v) - ((0,1)*CONJGU11*TH33)/v -
     + (TH13*yuku1)/v
      MGVX248(2) = -((GU11*TH23)/v) + ((0,1)*GU11*TH33)/v - (TH13*y
     +uku1)/v
      MGVX150(1) = (CONJGU32*Sqrt2)/v
      MGVX150(2) = -((GD23*Sqrt2)/v)
      MGVX151(1) = (CONJGU12*Sqrt2)/v
      MGVX151(2) = -((GD21*Sqrt2)/v)
      MGVX152(1) = (CONJGU22*Sqrt2)/v
      MGVX152(2) = -((GD22*Sqrt2)/v)
      MGVX153(1) = (CONJGU33*Sqrt2)/v
      MGVX153(2) = -((GD33*Sqrt2)/v)
      MGVX154(1) = (CONJGU13*Sqrt2)/v
      MGVX154(2) = -((GD31*Sqrt2)/v)
      MGVX155(1) = (CONJGU23*Sqrt2)/v
      MGVX155(2) = -((GD32*Sqrt2)/v)
      MGVX156(1) = (CONJGU31*Sqrt2)/v
      MGVX156(2) = -((GD13*Sqrt2)/v)
      MGVX157(1) = (CONJGU11*Sqrt2)/v
      MGVX157(2) = -((GD11*Sqrt2)/v)
      MGVX158(1) = (CONJGU21*Sqrt2)/v
      MGVX158(2) = -((GD12*Sqrt2)/v)
      MGVX132(1) = -((CONJGD23*Sqrt2)/v)
      MGVX132(2) = (GU32*Sqrt2)/v
      MGVX133(1) = -((CONJGD33*Sqrt2)/v)
      MGVX133(2) = (GU33*Sqrt2)/v
      MGVX134(1) = -((CONJGD13*Sqrt2)/v)
      MGVX134(2) = (GU31*Sqrt2)/v
      MGVX135(1) = -((CONJGD21*Sqrt2)/v)
      MGVX135(2) = (GU12*Sqrt2)/v
      MGVX136(1) = -((CONJGD31*Sqrt2)/v)
      MGVX136(2) = (GU13*Sqrt2)/v
      MGVX137(1) = -((CONJGD11*Sqrt2)/v)
      MGVX137(2) = (GU11*Sqrt2)/v
      MGVX138(1) = -((CONJGD22*Sqrt2)/v)
      MGVX138(2) = (GU22*Sqrt2)/v
      MGVX139(1) = -((CONJGD32*Sqrt2)/v)
      MGVX139(2) = (GU23*Sqrt2)/v
      MGVX140(1) = -((CONJGD12*Sqrt2)/v)
      MGVX140(2) = (GU21*Sqrt2)/v
      MGVX27(1) = -((CKM21*ee)/(Sqrt2*sw))
      MGVX27(2) = 0
      MGVX28(1) = -((CKM22*ee)/(Sqrt2*sw))
      MGVX28(2) = 0
      MGVX30(1) = -((CKM11*ee)/(Sqrt2*sw))
      MGVX30(2) = 0
      MGVX31(1) = -((CKM12*ee)/(Sqrt2*sw))
      MGVX31(2) = 0
      MGVX36(1) = -((CONJCKM21*ee)/(Sqrt2*sw))
      MGVX36(2) = 0
      MGVX37(1) = -((CONJCKM11*ee)/(Sqrt2*sw))
      MGVX37(2) = 0
      MGVX38(1) = -((CONJCKM22*ee)/(Sqrt2*sw))
      MGVX38(2) = 0
      MGVX39(1) = -((CONJCKM12*ee)/(Sqrt2*sw))
      MGVX39(2) = 0
      MGVX40(1) = (cw*ee)/(2.*sw) + (ee*sw)/(6.*cw)
      MGVX40(2) = -(ee*sw)/(3.*cw)
      MGVX46(1) = -(cw*ee)/(2.*sw) + (ee*sw)/(6.*cw)
      MGVX46(2) = (2*ee*sw)/(3.*cw)
      endif

      return
      end

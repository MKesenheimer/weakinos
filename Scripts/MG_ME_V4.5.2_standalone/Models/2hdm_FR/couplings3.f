      subroutine coup3(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 3
      if(readlha) then


      MGVX235(1) = -((CONJGL22*TH23)/v) + ((0,1)*CONJGL22*TH33)/v -
     + (TH13*yukl2)/v
      MGVX235(2) = -((GL22*TH23)/v) - ((0,1)*GL22*TH33)/v - (TH13*y
     +ukl2)/v
      MGVX236(1) = -((CONJGL32*TH23)/v) + ((0,1)*CONJGL32*TH33)/v
      MGVX236(2) = -((GL23*TH23)/v) - ((0,1)*GL23*TH33)/v
      MGVX237(1) = -((CONJGL13*TH23)/v) + ((0,1)*CONJGL13*TH33)/v
      MGVX237(2) = -((GL31*TH23)/v) - ((0,1)*GL31*TH33)/v
      MGVX238(1) = -((CONJGL23*TH23)/v) + ((0,1)*CONJGL23*TH33)/v
      MGVX238(2) = -((GL32*TH23)/v) - ((0,1)*GL32*TH33)/v
      MGVX239(1) = -((CONJGL33*TH23)/v) + ((0,1)*CONJGL33*TH33)/v -
     + (TH13*yukl3)/v
      MGVX239(2) = -((GL33*TH23)/v) - ((0,1)*GL33*TH33)/v - (TH13*y
     +ukl3)/v
      MGVX15(1) = ee/3.
      MGVX15(2) = ee/3.
      MGVX21(1) = (-2*ee)/3.
      MGVX21(2) = (-2*ee)/3.
      MGVX168(1) = -((CONJGD33*TH21)/v) + ((0,1)*CONJGD33*TH31)/v -
     + (TH11*yukd3)/v
      MGVX168(2) = -((GD33*TH21)/v) - ((0,1)*GD33*TH31)/v - (TH11*y
     +ukd3)/v
      MGVX169(1) = -((CONJGD13*TH21)/v) + ((0,1)*CONJGD13*TH31)/v
      MGVX169(2) = -((GD31*TH21)/v) - ((0,1)*GD31*TH31)/v
      MGVX170(1) = -((CONJGD23*TH21)/v) + ((0,1)*CONJGD23*TH31)/v
      MGVX170(2) = -((GD32*TH21)/v) - ((0,1)*GD32*TH31)/v
      MGVX186(1) = -((CONJGU22*TH21)/v) - ((0,1)*CONJGU22*TH31)/v -
     + (TH11*yuku2)/v
      MGVX186(2) = -((GU22*TH21)/v) + ((0,1)*GU22*TH31)/v - (TH11*y
     +uku2)/v
      MGVX187(1) = -((CONJGU32*TH21)/v) - ((0,1)*CONJGU32*TH31)/v
      MGVX187(2) = -((GU23*TH21)/v) + ((0,1)*GU23*TH31)/v
      MGVX188(1) = -((CONJGU12*TH21)/v) - ((0,1)*CONJGU12*TH31)/v
      MGVX188(2) = -((GU21*TH21)/v) + ((0,1)*GU21*TH31)/v
      MGVX171(1) = -((CONJGD31*TH21)/v) + ((0,1)*CONJGD31*TH31)/v
      MGVX171(2) = -((GD13*TH21)/v) - ((0,1)*GD13*TH31)/v
      MGVX172(1) = -((CONJGD11*TH21)/v) + ((0,1)*CONJGD11*TH31)/v -
     + (TH11*yukd1)/v
      MGVX172(2) = -((GD11*TH21)/v) - ((0,1)*GD11*TH31)/v - (TH11*y
     +ukd1)/v
      MGVX173(1) = -((CONJGD21*TH21)/v) + ((0,1)*CONJGD21*TH31)/v
      MGVX173(2) = -((GD12*TH21)/v) - ((0,1)*GD12*TH31)/v
      MGVX174(1) = -((CONJGD32*TH21)/v) + ((0,1)*CONJGD32*TH31)/v
      MGVX174(2) = -((GD23*TH21)/v) - ((0,1)*GD23*TH31)/v
      MGVX175(1) = -((CONJGD12*TH21)/v) + ((0,1)*CONJGD12*TH31)/v
      MGVX175(2) = -((GD21*TH21)/v) - ((0,1)*GD21*TH31)/v
      MGVX176(1) = -((CONJGD22*TH21)/v) + ((0,1)*CONJGD22*TH31)/v -
     + (TH11*yukd2)/v
      MGVX176(2) = -((GD22*TH21)/v) - ((0,1)*GD22*TH31)/v - (TH11*y
     +ukd2)/v
      MGVX189(1) = -((CONJGU23*TH21)/v) - ((0,1)*CONJGU23*TH31)/v
      MGVX189(2) = -((GU32*TH21)/v) + ((0,1)*GU32*TH31)/v
      MGVX190(1) = -((CONJGU33*TH21)/v) - ((0,1)*CONJGU33*TH31)/v -
     + (TH11*yuku3)/v
      MGVX190(2) = -((GU33*TH21)/v) + ((0,1)*GU33*TH31)/v - (TH11*y
     +uku3)/v
      MGVX191(1) = -((CONJGU13*TH21)/v) - ((0,1)*CONJGU13*TH31)/v
      MGVX191(2) = -((GU31*TH21)/v) + ((0,1)*GU31*TH31)/v
      MGVX192(1) = -((CONJGU21*TH21)/v) - ((0,1)*CONJGU21*TH31)/v
      MGVX192(2) = -((GU12*TH21)/v) + ((0,1)*GU12*TH31)/v
      MGVX193(1) = -((CONJGU31*TH21)/v) - ((0,1)*CONJGU31*TH31)/v
      MGVX193(2) = -((GU13*TH21)/v) + ((0,1)*GU13*TH31)/v
      MGVX194(1) = -((CONJGU11*TH21)/v) - ((0,1)*CONJGU11*TH31)/v -
     + (TH11*yuku1)/v
      MGVX194(2) = -((GU11*TH21)/v) + ((0,1)*GU11*TH31)/v - (TH11*y
     +uku1)/v
      MGVX195(1) = -((CONJGD33*TH22)/v) + ((0,1)*CONJGD33*TH32)/v -
     + (TH12*yukd3)/v
      MGVX195(2) = -((GD33*TH22)/v) - ((0,1)*GD33*TH32)/v - (TH12*y
     +ukd3)/v
      MGVX196(1) = -((CONJGD13*TH22)/v) + ((0,1)*CONJGD13*TH32)/v
      MGVX196(2) = -((GD31*TH22)/v) - ((0,1)*GD31*TH32)/v
      MGVX197(1) = -((CONJGD23*TH22)/v) + ((0,1)*CONJGD23*TH32)/v
      MGVX197(2) = -((GD32*TH22)/v) - ((0,1)*GD32*TH32)/v
      MGVX213(1) = -((CONJGU22*TH22)/v) - ((0,1)*CONJGU22*TH32)/v -
     + (TH12*yuku2)/v
      MGVX213(2) = -((GU22*TH22)/v) + ((0,1)*GU22*TH32)/v - (TH12*y
     +uku2)/v
      MGVX214(1) = -((CONJGU32*TH22)/v) - ((0,1)*CONJGU32*TH32)/v
      MGVX214(2) = -((GU23*TH22)/v) + ((0,1)*GU23*TH32)/v
      MGVX215(1) = -((CONJGU12*TH22)/v) - ((0,1)*CONJGU12*TH32)/v
      MGVX215(2) = -((GU21*TH22)/v) + ((0,1)*GU21*TH32)/v
      MGVX198(1) = -((CONJGD31*TH22)/v) + ((0,1)*CONJGD31*TH32)/v
      MGVX198(2) = -((GD13*TH22)/v) - ((0,1)*GD13*TH32)/v
      MGVX199(1) = -((CONJGD11*TH22)/v) + ((0,1)*CONJGD11*TH32)/v -
     + (TH12*yukd1)/v
      MGVX199(2) = -((GD11*TH22)/v) - ((0,1)*GD11*TH32)/v - (TH12*y
     +ukd1)/v
      MGVX200(1) = -((CONJGD21*TH22)/v) + ((0,1)*CONJGD21*TH32)/v
      MGVX200(2) = -((GD12*TH22)/v) - ((0,1)*GD12*TH32)/v
      MGVX201(1) = -((CONJGD32*TH22)/v) + ((0,1)*CONJGD32*TH32)/v
      MGVX201(2) = -((GD23*TH22)/v) - ((0,1)*GD23*TH32)/v
      MGVX202(1) = -((CONJGD12*TH22)/v) + ((0,1)*CONJGD12*TH32)/v
      MGVX202(2) = -((GD21*TH22)/v) - ((0,1)*GD21*TH32)/v
      MGVX203(1) = -((CONJGD22*TH22)/v) + ((0,1)*CONJGD22*TH32)/v -
     + (TH12*yukd2)/v
      MGVX203(2) = -((GD22*TH22)/v) - ((0,1)*GD22*TH32)/v - (TH12*y
     +ukd2)/v
      MGVX216(1) = -((CONJGU23*TH22)/v) - ((0,1)*CONJGU23*TH32)/v
      MGVX216(2) = -((GU32*TH22)/v) + ((0,1)*GU32*TH32)/v
      MGVX217(1) = -((CONJGU33*TH22)/v) - ((0,1)*CONJGU33*TH32)/v -
     + (TH12*yuku3)/v
      MGVX217(2) = -((GU33*TH22)/v) + ((0,1)*GU33*TH32)/v - (TH12*y
     +uku3)/v
      MGVX218(1) = -((CONJGU13*TH22)/v) - ((0,1)*CONJGU13*TH32)/v
      MGVX218(2) = -((GU31*TH22)/v) + ((0,1)*GU31*TH32)/v
      MGVX219(1) = -((CONJGU21*TH22)/v) - ((0,1)*CONJGU21*TH32)/v
      MGVX219(2) = -((GU12*TH22)/v) + ((0,1)*GU12*TH32)/v
      MGVX220(1) = -((CONJGU31*TH22)/v) - ((0,1)*CONJGU31*TH32)/v
      MGVX220(2) = -((GU13*TH22)/v) + ((0,1)*GU13*TH32)/v
      MGVX221(1) = -((CONJGU11*TH22)/v) - ((0,1)*CONJGU11*TH32)/v -
     + (TH12*yuku1)/v
      MGVX221(2) = -((GU11*TH22)/v) + ((0,1)*GU11*TH32)/v - (TH12*y
     +uku1)/v
      MGVX222(1) = -((CONJGD33*TH23)/v) + ((0,1)*CONJGD33*TH33)/v -
     + (TH13*yukd3)/v
      MGVX222(2) = -((GD33*TH23)/v) - ((0,1)*GD33*TH33)/v - (TH13*y
     +ukd3)/v
      MGVX223(1) = -((CONJGD13*TH23)/v) + ((0,1)*CONJGD13*TH33)/v
      MGVX223(2) = -((GD31*TH23)/v) - ((0,1)*GD31*TH33)/v
      MGVX224(1) = -((CONJGD23*TH23)/v) + ((0,1)*CONJGD23*TH33)/v
      MGVX224(2) = -((GD32*TH23)/v) - ((0,1)*GD32*TH33)/v
      MGVX240(1) = -((CONJGU22*TH23)/v) - ((0,1)*CONJGU22*TH33)/v -
     + (TH13*yuku2)/v
      MGVX240(2) = -((GU22*TH23)/v) + ((0,1)*GU22*TH33)/v - (TH13*y
     +uku2)/v
      MGVX241(1) = -((CONJGU32*TH23)/v) - ((0,1)*CONJGU32*TH33)/v
      MGVX241(2) = -((GU23*TH23)/v) + ((0,1)*GU23*TH33)/v
      MGVX242(1) = -((CONJGU12*TH23)/v) - ((0,1)*CONJGU12*TH33)/v
      MGVX242(2) = -((GU21*TH23)/v) + ((0,1)*GU21*TH33)/v
      MGVX225(1) = -((CONJGD31*TH23)/v) + ((0,1)*CONJGD31*TH33)/v
      MGVX225(2) = -((GD13*TH23)/v) - ((0,1)*GD13*TH33)/v
      endif

      return
      end

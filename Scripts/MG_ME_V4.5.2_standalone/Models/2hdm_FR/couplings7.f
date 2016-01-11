      subroutine coup7(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 7
      if(readlha) then


      MGVX203(1) = -((GDR33*TH22)/v) + ((0,1)*GDR33*TH32)/v - (TH12
     +*yukd3)/v
      MGVX203(2) = -((GDR33*TH22)/v) - ((0,1)*GDR33*TH32)/v - (TH12
     +*yukd3)/v
      MGVX197(1) = -((GDR13*TH22)/v) + ((0,1)*GDR13*TH32)/v
      MGVX197(2) = -((GDR31*TH22)/v) - ((0,1)*GDR31*TH32)/v
      MGVX200(1) = -((GDR23*TH22)/v) + ((0,1)*GDR23*TH32)/v
      MGVX200(2) = -((GDR32*TH22)/v) - ((0,1)*GDR32*TH32)/v
      MGVX217(1) = -((GUR22*TH22)/v) - ((0,1)*GUR22*TH32)/v - (TH12
     +*yuku2)/v
      MGVX217(2) = -((GUR22*TH22)/v) + ((0,1)*GUR22*TH32)/v - (TH12
     +*yuku2)/v
      MGVX220(1) = -((GUR32*TH22)/v) - ((0,1)*GUR32*TH32)/v
      MGVX220(2) = -((GUR23*TH22)/v) + ((0,1)*GUR23*TH32)/v
      MGVX214(1) = -((GUR12*TH22)/v) - ((0,1)*GUR12*TH32)/v
      MGVX214(2) = -((GUR21*TH22)/v) + ((0,1)*GUR21*TH32)/v
      MGVX201(1) = -((GDR31*TH22)/v) + ((0,1)*GDR31*TH32)/v
      MGVX201(2) = -((GDR13*TH22)/v) - ((0,1)*GDR13*TH32)/v
      MGVX195(1) = -((GDR11*TH22)/v) + ((0,1)*GDR11*TH32)/v - (TH12
     +*yukd1)/v
      MGVX195(2) = -((GDR11*TH22)/v) - ((0,1)*GDR11*TH32)/v - (TH12
     +*yukd1)/v
      MGVX198(1) = -((GDR21*TH22)/v) + ((0,1)*GDR21*TH32)/v
      MGVX198(2) = -((GDR12*TH22)/v) - ((0,1)*GDR12*TH32)/v
      MGVX202(1) = -((GDR32*TH22)/v) + ((0,1)*GDR32*TH32)/v
      MGVX202(2) = -((GDR23*TH22)/v) - ((0,1)*GDR23*TH32)/v
      MGVX196(1) = -((GDR12*TH22)/v) + ((0,1)*GDR12*TH32)/v
      MGVX196(2) = -((GDR21*TH22)/v) - ((0,1)*GDR21*TH32)/v
      MGVX199(1) = -((GDR22*TH22)/v) + ((0,1)*GDR22*TH32)/v - (TH12
     +*yukd2)/v
      MGVX199(2) = -((GDR22*TH22)/v) - ((0,1)*GDR22*TH32)/v - (TH12
     +*yukd2)/v
      MGVX218(1) = -((GUR23*TH22)/v) - ((0,1)*GUR23*TH32)/v
      MGVX218(2) = -((GUR32*TH22)/v) + ((0,1)*GUR32*TH32)/v
      MGVX221(1) = -((GUR33*TH22)/v) - ((0,1)*GUR33*TH32)/v - (TH12
     +*yuku3)/v
      MGVX221(2) = -((GUR33*TH22)/v) + ((0,1)*GUR33*TH32)/v - (TH12
     +*yuku3)/v
      MGVX215(1) = -((GUR13*TH22)/v) - ((0,1)*GUR13*TH32)/v
      MGVX215(2) = -((GUR31*TH22)/v) + ((0,1)*GUR31*TH32)/v
      endif

      return
      end

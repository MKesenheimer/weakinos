      subroutine coup8(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 8
      if(readlha) then


      MGVX216(1) = -((GUR21*TH22)/v) - ((0,1)*GUR21*TH32)/v
      MGVX216(2) = -((GUR12*TH22)/v) + ((0,1)*GUR12*TH32)/v
      MGVX219(1) = -((GUR31*TH22)/v) - ((0,1)*GUR31*TH32)/v
      MGVX219(2) = -((GUR13*TH22)/v) + ((0,1)*GUR13*TH32)/v
      MGVX213(1) = -((GUR11*TH22)/v) - ((0,1)*GUR11*TH32)/v - (TH12
     +*yuku1)/v
      MGVX213(2) = -((GUR11*TH22)/v) + ((0,1)*GUR11*TH32)/v - (TH12
     +*yuku1)/v
      MGVX230(1) = -((GDR33*TH23)/v) + ((0,1)*GDR33*TH33)/v - (TH13
     +*yukd3)/v
      MGVX230(2) = -((GDR33*TH23)/v) - ((0,1)*GDR33*TH33)/v - (TH13
     +*yukd3)/v
      MGVX224(1) = -((GDR13*TH23)/v) + ((0,1)*GDR13*TH33)/v
      MGVX224(2) = -((GDR31*TH23)/v) - ((0,1)*GDR31*TH33)/v
      MGVX227(1) = -((GDR23*TH23)/v) + ((0,1)*GDR23*TH33)/v
      MGVX227(2) = -((GDR32*TH23)/v) - ((0,1)*GDR32*TH33)/v
      MGVX244(1) = -((GUR22*TH23)/v) - ((0,1)*GUR22*TH33)/v - (TH13
     +*yuku2)/v
      MGVX244(2) = -((GUR22*TH23)/v) + ((0,1)*GUR22*TH33)/v - (TH13
     +*yuku2)/v
      MGVX247(1) = -((GUR32*TH23)/v) - ((0,1)*GUR32*TH33)/v
      MGVX247(2) = -((GUR23*TH23)/v) + ((0,1)*GUR23*TH33)/v
      MGVX241(1) = -((GUR12*TH23)/v) - ((0,1)*GUR12*TH33)/v
      MGVX241(2) = -((GUR21*TH23)/v) + ((0,1)*GUR21*TH33)/v
      MGVX228(1) = -((GDR31*TH23)/v) + ((0,1)*GDR31*TH33)/v
      MGVX228(2) = -((GDR13*TH23)/v) - ((0,1)*GDR13*TH33)/v
      MGVX222(1) = -((GDR11*TH23)/v) + ((0,1)*GDR11*TH33)/v - (TH13
     +*yukd1)/v
      MGVX222(2) = -((GDR11*TH23)/v) - ((0,1)*GDR11*TH33)/v - (TH13
     +*yukd1)/v
      MGVX225(1) = -((GDR21*TH23)/v) + ((0,1)*GDR21*TH33)/v
      MGVX225(2) = -((GDR12*TH23)/v) - ((0,1)*GDR12*TH33)/v
      MGVX229(1) = -((GDR32*TH23)/v) + ((0,1)*GDR32*TH33)/v
      MGVX229(2) = -((GDR23*TH23)/v) - ((0,1)*GDR23*TH33)/v
      MGVX223(1) = -((GDR12*TH23)/v) + ((0,1)*GDR12*TH33)/v
      MGVX223(2) = -((GDR21*TH23)/v) - ((0,1)*GDR21*TH33)/v
      MGVX226(1) = -((GDR22*TH23)/v) + ((0,1)*GDR22*TH33)/v - (TH13
     +*yukd2)/v
      MGVX226(2) = -((GDR22*TH23)/v) - ((0,1)*GDR22*TH33)/v - (TH13
     +*yukd2)/v
      endif

      return
      end

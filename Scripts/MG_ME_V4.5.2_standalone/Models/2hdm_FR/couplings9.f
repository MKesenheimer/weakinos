      subroutine coup9(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 9
      if(readlha) then


      MGVX245(1) = -((GUR23*TH23)/v) - ((0,1)*GUR23*TH33)/v
      MGVX245(2) = -((GUR32*TH23)/v) + ((0,1)*GUR32*TH33)/v
      MGVX248(1) = -((GUR33*TH23)/v) - ((0,1)*GUR33*TH33)/v - (TH13
     +*yuku3)/v
      MGVX248(2) = -((GUR33*TH23)/v) + ((0,1)*GUR33*TH33)/v - (TH13
     +*yuku3)/v
      MGVX242(1) = -((GUR13*TH23)/v) - ((0,1)*GUR13*TH33)/v
      MGVX242(2) = -((GUR31*TH23)/v) + ((0,1)*GUR31*TH33)/v
      MGVX243(1) = -((GUR21*TH23)/v) - ((0,1)*GUR21*TH33)/v
      MGVX243(2) = -((GUR12*TH23)/v) + ((0,1)*GUR12*TH33)/v
      MGVX246(1) = -((GUR31*TH23)/v) - ((0,1)*GUR31*TH33)/v
      MGVX246(2) = -((GUR13*TH23)/v) + ((0,1)*GUR13*TH33)/v
      MGVX240(1) = -((GUR11*TH23)/v) - ((0,1)*GUR11*TH33)/v - (TH13
     +*yuku1)/v
      MGVX240(2) = -((GUR11*TH23)/v) + ((0,1)*GUR11*TH33)/v - (TH13
     +*yuku1)/v
      MGVX157(1) = (GUR32*Sqrt2)/v
      MGVX157(2) = -((GDR23*Sqrt2)/v)
      MGVX151(1) = (GUR12*Sqrt2)/v
      MGVX151(2) = -((GDR21*Sqrt2)/v)
      MGVX154(1) = (GUR22*Sqrt2)/v
      MGVX154(2) = -((GDR22*Sqrt2)/v)
      MGVX158(1) = (GUR33*Sqrt2)/v
      MGVX158(2) = -((GDR33*Sqrt2)/v)
      MGVX152(1) = (GUR13*Sqrt2)/v
      MGVX152(2) = -((GDR31*Sqrt2)/v)
      MGVX155(1) = (GUR23*Sqrt2)/v
      MGVX155(2) = -((GDR32*Sqrt2)/v)
      MGVX156(1) = (GUR31*Sqrt2)/v
      MGVX156(2) = -((GDR13*Sqrt2)/v)
      MGVX150(1) = (GUR11*Sqrt2)/v
      MGVX150(2) = -((GDR11*Sqrt2)/v)
      MGVX153(1) = (GUR21*Sqrt2)/v
      MGVX153(2) = -((GDR12*Sqrt2)/v)
      endif

      return
      end

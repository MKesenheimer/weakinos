      subroutine coup10(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 10
      if(readlha) then


      MGVX137(1) = -((GDR23*Sqrt2)/v)
      MGVX137(2) = (GUR32*Sqrt2)/v
      MGVX140(1) = -((GDR33*Sqrt2)/v)
      MGVX140(2) = (GUR33*Sqrt2)/v
      MGVX134(1) = -((GDR13*Sqrt2)/v)
      MGVX134(2) = (GUR31*Sqrt2)/v
      MGVX135(1) = -((GDR21*Sqrt2)/v)
      MGVX135(2) = (GUR12*Sqrt2)/v
      MGVX138(1) = -((GDR31*Sqrt2)/v)
      MGVX138(2) = (GUR13*Sqrt2)/v
      MGVX132(1) = -((GDR11*Sqrt2)/v)
      MGVX132(2) = (GUR11*Sqrt2)/v
      MGVX136(1) = -((GDR22*Sqrt2)/v)
      MGVX136(2) = (GUR22*Sqrt2)/v
      MGVX139(1) = -((GDR32*Sqrt2)/v)
      MGVX139(2) = (GUR23*Sqrt2)/v
      MGVX133(1) = -((GDR12*Sqrt2)/v)
      MGVX133(2) = (GUR21*Sqrt2)/v
      MGVX22(1) = -((CKMR21*ee)/(Sqrt2*sw))
      MGVX22(2) = 0
      MGVX23(1) = -((CKMR22*ee)/(Sqrt2*sw))
      MGVX23(2) = 0
      MGVX25(1) = -((CKMR11*ee)/(Sqrt2*sw))
      MGVX25(2) = 0
      MGVX26(1) = -((CKMR12*ee)/(Sqrt2*sw))
      MGVX26(2) = 0
      MGVX35(1) = (cw*ee)/(2.*sw) + (ee*sw)/(6.*cw)
      MGVX35(2) = -(ee*sw)/(3.*cw)
      MGVX41(1) = -(cw*ee)/(2.*sw) + (ee*sw)/(6.*cw)
      MGVX41(2) = (2*ee*sw)/(3.*cw)
      endif

      return
      end

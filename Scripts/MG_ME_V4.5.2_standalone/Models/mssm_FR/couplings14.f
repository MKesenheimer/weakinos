      subroutine coup14(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 14
      if(readlha) then


      MGVX667(1) = (ee*OL22)/sw
      MGVX667(2) = (ee*OR22)/sw
      MGVX668(1) = (ee*OL31)/sw
      MGVX668(2) = (ee*OR31)/sw
      MGVX669(1) = (ee*OL32)/sw
      MGVX669(2) = (ee*OR32)/sw
      MGVX670(1) = (ee*OL41)/sw
      MGVX670(2) = (ee*OR41)/sw
      MGVX671(1) = (ee*OL42)/sw
      MGVX671(2) = (ee*OR42)/sw
      MGVX672(1) = (ee*OpL11)/(cw*sw)
      MGVX672(2) = (ee*OpR11)/(cw*sw)
      MGVX673(1) = (ee*OpL12)/(cw*sw)
      MGVX673(2) = (ee*OpR12)/(cw*sw)
      MGVX674(1) = (ee*OpL21)/(cw*sw)
      MGVX674(2) = (ee*OpR21)/(cw*sw)
      MGVX675(1) = (ee*OpL22)/(cw*sw)
      MGVX675(2) = (ee*OpR22)/(cw*sw)
      MGVX676(1) = (ee*OppL11)/(2.*cw*sw) - (ee*OppR11)/(2.*cw*sw)
      MGVX676(2) = -(ee*OppL11)/(2.*cw*sw) + (ee*OppR11)/(2.*cw*sw)
      MGVX677(1) = (ee*OppL12)/(2.*cw*sw) - (ee*OppR21)/(2.*cw*sw)
      MGVX677(2) = -(ee*OppL21)/(2.*cw*sw) + (ee*OppR12)/(2.*cw*sw)
      MGVX678(1) = (ee*OppL13)/(2.*cw*sw) - (ee*OppR31)/(2.*cw*sw)
      MGVX678(2) = -(ee*OppL31)/(2.*cw*sw) + (ee*OppR13)/(2.*cw*sw)
      MGVX679(1) = (ee*OppL14)/(2.*cw*sw) - (ee*OppR41)/(2.*cw*sw)
      MGVX679(2) = -(ee*OppL41)/(2.*cw*sw) + (ee*OppR14)/(2.*cw*sw)
      MGVX680(1) = (ee*OppL22)/(2.*cw*sw) - (ee*OppR22)/(2.*cw*sw)
      MGVX680(2) = -(ee*OppL22)/(2.*cw*sw) + (ee*OppR22)/(2.*cw*sw)
      MGVX681(1) = (ee*OppL23)/(2.*cw*sw) - (ee*OppR32)/(2.*cw*sw)
      MGVX681(2) = -(ee*OppL32)/(2.*cw*sw) + (ee*OppR23)/(2.*cw*sw)
      endif

      return
      end

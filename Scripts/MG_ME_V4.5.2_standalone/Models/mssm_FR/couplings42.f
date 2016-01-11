      subroutine coup42(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 42
      if(readlha) then

      MGVX250 = -(CONJRu63*ee*Ru13)/(2.*cw*sw)
      MGVX126 = (CONJRd13*ee*Ru63)/(Sqrt2*sw)
      MGVX127 = (CONJRd23*ee*Ru63)/(Sqrt2*sw)
      MGVX249 = (CONJRu13*ee*Ru63)/(2.*cw*sw)
      MGVX265 = (CONJRu63*ee*Ru63)/(2.*cw*sw) - (2*ee*sw)/(3.*cw)

      endif

      return
      end

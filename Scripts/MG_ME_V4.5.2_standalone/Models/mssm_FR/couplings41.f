      subroutine coup41(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 41
      if(readlha) then

      MGVX121 = (CONJRl51*ee*RLn31)/(Sqrt2*sw)
      MGVX146 = -((CONJRu13*ee*Rd13)/(Sqrt2*sw))
      MGVX148 = -((CONJRu13*ee*Rd23)/(Sqrt2*sw))
      MGVX122 = (CONJRd13*ee*Ru13)/(Sqrt2*sw)
      MGVX123 = (CONJRd23*ee*Ru13)/(Sqrt2*sw)
      MGVX248 = (CONJRu13*ee*Ru13)/(2.*cw*sw) - (2*ee*sw)/(3.*cw)
      MGVX253 = (-2*ee*sw)/(3.*cw)
      MGVX151 = -((CONJRu42*ee*Rd62)/(Sqrt2*sw))
      MGVX124 = (CONJRd62*ee*Ru42)/(Sqrt2*sw)
      MGVX259 = (CONJRu42*ee*Ru42)/(2.*cw*sw) - (2*ee*sw)/(3.*cw)
      MGVX150 = -((CONJRu51*ee*Rd51)/(Sqrt2*sw))
      MGVX125 = (CONJRd51*ee*Ru51)/(Sqrt2*sw)
      MGVX262 = (CONJRu51*ee*Ru51)/(2.*cw*sw) - (2*ee*sw)/(3.*cw)
      MGVX147 = -((CONJRu63*ee*Rd13)/(Sqrt2*sw))
      MGVX149 = -((CONJRu63*ee*Rd23)/(Sqrt2*sw))

      endif

      return
      end

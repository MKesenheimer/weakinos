      subroutine coup40(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 40
      if(readlha) then

      MGVX210 = -(CONJRl13*ee*Rl13)/(2.*cw*sw) + (ee*sw)/cw
      MGVX211 = (ee*sw)/cw
      MGVX213 = -(CONJRl42*ee*Rl42)/(2.*cw*sw) + (ee*sw)/cw
      MGVX214 = -(CONJRl51*ee*Rl51)/(2.*cw*sw) + (ee*sw)/cw
      MGVX236 = (CONJRl63*ee*Rl13)/(2.*cw*sw)
      MGVX238 = -(CONJRl13*ee*Rl63)/(2.*cw*sw)
      MGVX215 = -(CONJRl63*ee*Rl63)/(2.*cw*sw) + (ee*sw)/cw
      MGVX153 = -((CONJRLn13*ee*Rl13)/(Sqrt2*sw))
      MGVX159 = -((CONJRLn13*ee*Rl63)/(Sqrt2*sw))
      MGVX116 = (CONJRl13*ee*RLn13)/(Sqrt2*sw)
      MGVX117 = (CONJRl63*ee*RLn13)/(Sqrt2*sw)
      MGVX216 = ee/(2.*cw*sw)
      MGVX155 = -((CONJRLn22*ee*Rl42)/(Sqrt2*sw))
      MGVX119 = (CONJRl42*ee*RLn22)/(Sqrt2*sw)
      MGVX157 = -((CONJRLn31*ee*Rl51)/(Sqrt2*sw))

      endif

      return
      end

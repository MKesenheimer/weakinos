      subroutine coup32(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 32
      if(readlha) then

      MGVX163 = (CONJRu13*ee**2*Rd23)/(3.*Sqrt2*sw)
      MGVX171 = (CONJRu63*ee**2*Rd23)/(3.*Sqrt2*sw)
      MGVX224 = (-2*ee**2*sw)/(9.*cw)
      MGVX228 = (CONJRd51*ee**2*Rd51)/(3.*cw*sw) - (2*ee**2*sw)/(9.
     +*cw)
      MGVX136 = (CONJRd51*ee**2*Ru51)/(3.*Sqrt2*sw)
      MGVX168 = (CONJRu51*ee**2*Rd51)/(3.*Sqrt2*sw)
      MGVX230 = (CONJRd62*ee**2*Rd62)/(3.*cw*sw) - (2*ee**2*sw)/(9.
     +*cw)
      MGVX134 = (CONJRd62*ee**2*Ru42)/(3.*Sqrt2*sw)
      MGVX166 = (CONJRu42*ee**2*Rd62)/(3.*Sqrt2*sw)
      MGVX251 = (2*CONJRu13*ee**2*Ru13)/(3.*cw*sw) - (8*ee**2*sw)/(
     +9.*cw)
      MGVX268 = (2*CONJRu13*ee**2*Ru63)/(3.*cw*sw)
      MGVX266 = (2*CONJRu63*ee**2*Ru13)/(3.*cw*sw)
      MGVX254 = (-8*ee**2*sw)/(9.*cw)
      MGVX260 = (2*CONJRu42*ee**2*Ru42)/(3.*cw*sw) - (8*ee**2*sw)/(
     +9.*cw)
      MGVX263 = (2*CONJRu51*ee**2*Ru51)/(3.*cw*sw) - (8*ee**2*sw)/(
     +9.*cw)

      endif

      return
      end

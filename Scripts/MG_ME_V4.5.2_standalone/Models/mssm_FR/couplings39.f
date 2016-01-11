      subroutine coup39(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 39
      if(readlha) then

      MGVX219 = -(CONJRd13*ee*Rd13)/(2.*cw*sw) + (ee*sw)/(3.*cw)
      MGVX231 = (CONJRd23*ee*Rd13)/(2.*cw*sw)
      MGVX233 = -(CONJRd13*ee*Rd23)/(2.*cw*sw)
      MGVX221 = -(CONJRd23*ee*Rd23)/(2.*cw*sw) + (ee*sw)/(3.*cw)
      MGVX223 = (ee*sw)/(3.*cw)
      MGVX227 = -(CONJRd51*ee*Rd51)/(2.*cw*sw) + (ee*sw)/(3.*cw)
      MGVX229 = -(CONJRd62*ee*Rd62)/(2.*cw*sw) + (ee*sw)/(3.*cw)
      MGVX144 = ((0,0.5)*ee*Cos(bet)**2)/sw + ((0,0.5)*ee*Sin(bet)*
     +*2)/sw
      MGVX208 = ((0,-0.5)*Cosalp*Cosbet*ee)/(cw*sw) - ((0,0.5)*ee*S
     +inalp*Sinbet)/(cw*sw)
      MGVX145 = (Cosalp*Cosbet*ee)/(2.*sw) + (ee*Sinalp*Sinbet)/(2.
     +*sw)
      MGVX239 = ((0,-0.5)*Cosbet*ee*Sinalp)/(cw*sw) + ((0,0.5)*Cosa
     +lp*ee*Sinbet)/(cw*sw)
      MGVX161 = (Cosbet*ee*Sinalp)/(2.*sw) - (Cosalp*ee*Sinbet)/(2.
     +*sw)
      MGVX55 = ee*Cos(bet)**2 + ee*Sin(bet)**2
      MGVX209 = (ee*Cos(bet)**2)/(2.*cw*sw) - (ee*sw*Cos(bet)**2)/c
     +w + (ee*Sin(bet)**2)/(2.*cw*sw) - (ee*sw*Sin(bet)**2)/cw
      MGVX56 = -ee

      endif

      return
      end

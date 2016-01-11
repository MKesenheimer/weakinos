      subroutine coup10(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 10
      if(readlha) then


      MGVX118(1) = (g0*v0W*vL0bo*vL0to)/Sqrt2 + (gt*v1W*vL1bo*vL1to
     +)/Sqrt2
      MGVX118(2) = (gt*v1W*vR1bo*vR1to)/Sqrt2
      MGVX119(1) = (g0*v0W*vL0hBo*vL0to)/Sqrt2 + (gt*v1W*vL1hBo*vL1
     +to)/Sqrt2
      MGVX119(2) = (gt*v1W*vR1hBo*vR1to)/Sqrt2
      MGVX156(1) = (g0*v0WP*vL0do*vL0hUp)/Sqrt2 + (gt*v1WP*vL1do*vL
     +1hUp)/Sqrt2
      MGVX156(2) = 0
      MGVX157(1) = (g0*v0WP*vL0hDo*vL0hUp)/Sqrt2 + (gt*v1WP*vL1hDo*
     +vL1hUp)/Sqrt2
      MGVX157(2) = (gt*v1WP)/Sqrt2
      MGVX158(1) = (g0*v0WP*vL0hCh*vL0st)/Sqrt2 + (gt*v1WP*vL1hCh*v
     +L1st)/Sqrt2
      MGVX158(2) = (gt*v1WP*vR1hCh*vR1st)/Sqrt2
      MGVX159(1) = (g0*v0WP*vL0hCh*vL0hSt)/Sqrt2 + (gt*v1WP*vL1hCh*
     +vL1hSt)/Sqrt2
      MGVX159(2) = (gt*v1WP*vR1hCh*vR1hSt)/Sqrt2
      MGVX160(1) = (g0*v0WP*vL0bo*vL0hTo)/Sqrt2 + (gt*v1WP*vL1bo*vL
     +1hTo)/Sqrt2
      MGVX160(2) = (gt*v1WP*vR1bo*vR1hTo)/Sqrt2
      MGVX161(1) = (g0*v0WP*vL0hBo*vL0hTo)/Sqrt2 + (gt*v1WP*vL1hBo*
     +vL1hTo)/Sqrt2
      MGVX161(2) = (gt*v1WP*vR1hBo*vR1hTo)/Sqrt2
      MGVX162(1) = (g0*v0WP*vL0do*vL0up)/Sqrt2 + (gt*v1WP*vL1do*vL1
     +up)/Sqrt2
      MGVX162(2) = 0
      MGVX163(1) = (g0*v0WP*vL0hDo*vL0up)/Sqrt2 + (gt*v1WP*vL1hDo*v
     +L1up)/Sqrt2
      MGVX163(2) = 0
      MGVX164(1) = (g0*v0WP*vL0ch*vL0st)/Sqrt2 + (gt*v1WP*vL1ch*vL1
     +st)/Sqrt2
      MGVX164(2) = (gt*v1WP*vR1ch*vR1st)/Sqrt2
      MGVX165(1) = (g0*v0WP*vL0ch*vL0hSt)/Sqrt2 + (gt*v1WP*vL1ch*vL
     +1hSt)/Sqrt2
      MGVX165(2) = (gt*v1WP*vR1ch*vR1hSt)/Sqrt2
      MGVX166(1) = (g0*v0WP*vL0bo*vL0to)/Sqrt2 + (gt*v1WP*vL1bo*vL1
     +to)/Sqrt2
      MGVX166(2) = (gt*v1WP*vR1bo*vR1to)/Sqrt2
      MGVX167(1) = (g0*v0WP*vL0hBo*vL0to)/Sqrt2 + (gt*v1WP*vL1hBo*v
     +L1to)/Sqrt2
      MGVX167(2) = (gt*v1WP*vR1hBo*vR1to)/Sqrt2
      MGVX192(1) = -(g0*PvL0do2*v0Z)/2. - (gt*PvL1do2*v1Z)/2. + (gp
     +*PvL0do2*v2Z)/6. + (gp*PvL1do2*v2Z)/6.
      MGVX192(2) = -(gp*v2Z)/3.
      endif

      return
      end

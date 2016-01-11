      subroutine coup5(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 5
      if(readlha) then


      MGVX97(1) = (g0*v0W*vL0hEl*vL0hNu)/Sqrt2 + (gt*v1W*vL1hEl*vL1
     +hNu)/Sqrt2
      MGVX97(2) = (gt*v1W)/Sqrt2
      MGVX98(1) = (g0*v0W*vL0hNu*vL0mu)/Sqrt2 + (gt*v1W*vL1hNu*vL1m
     +u)/Sqrt2
      MGVX98(2) = -((gt*v1W*vR1mu)/Sqrt2)
      MGVX99(1) = (g0*v0W*vL0hMu*vL0hNu)/Sqrt2 + (gt*v1W*vL1hMu*vL1
     +hNu)/Sqrt2
      MGVX99(2) = -((gt*v1W*vR1hMu)/Sqrt2)
      MGVX100(1) = (g0*v0W*vL0hNu*vL0ta)/Sqrt2 + (gt*v1W*vL1hNu*vL1
     +ta)/Sqrt2
      MGVX100(2) = -((gt*v1W*vR1ta)/Sqrt2)
      MGVX101(1) = (g0*v0W*vL0hNu*vL0hTa)/Sqrt2 + (gt*v1W*vL1hNu*vL
     +1hTa)/Sqrt2
      MGVX101(2) = -((gt*v1W*vR1hTa)/Sqrt2)
      MGVX102(1) = (g0*v0W*vL0el*vL0nu)/Sqrt2 + (gt*v1W*vL1el*vL1nu
     +)/Sqrt2
      MGVX102(2) = 0
      MGVX103(1) = (g0*v0W*vL0hEl*vL0nu)/Sqrt2 + (gt*v1W*vL1hEl*vL1
     +nu)/Sqrt2
      MGVX103(2) = 0
      MGVX104(1) = (g0*v0W*vL0mu*vL0nu)/Sqrt2 + (gt*v1W*vL1mu*vL1nu
     +)/Sqrt2
      MGVX104(2) = 0
      MGVX105(1) = (g0*v0W*vL0hMu*vL0nu)/Sqrt2 + (gt*v1W*vL1hMu*vL1
     +nu)/Sqrt2
      MGVX105(2) = 0
      MGVX106(1) = (g0*v0W*vL0nu*vL0ta)/Sqrt2 + (gt*v1W*vL1nu*vL1ta
     +)/Sqrt2
      MGVX106(2) = 0
      MGVX107(1) = (g0*v0W*vL0hTa*vL0nu)/Sqrt2 + (gt*v1W*vL1hTa*vL1
     +nu)/Sqrt2
      MGVX107(2) = 0
      MGVX144(1) = (g0*v0WP*vL0el*vL0hNu)/Sqrt2 + (gt*v1WP*vL1el*vL
     +1hNu)/Sqrt2
      MGVX144(2) = 0
      MGVX145(1) = (g0*v0WP*vL0hEl*vL0hNu)/Sqrt2 + (gt*v1WP*vL1hEl*
     +vL1hNu)/Sqrt2
      MGVX145(2) = (gt*v1WP)/Sqrt2
      MGVX146(1) = (g0*v0WP*vL0hNu*vL0mu)/Sqrt2 + (gt*v1WP*vL1hNu*v
     +L1mu)/Sqrt2
      MGVX146(2) = -((gt*v1WP*vR1mu)/Sqrt2)
      MGVX147(1) = (g0*v0WP*vL0hMu*vL0hNu)/Sqrt2 + (gt*v1WP*vL1hMu*
     +vL1hNu)/Sqrt2
      MGVX147(2) = -((gt*v1WP*vR1hMu)/Sqrt2)
      endif

      return
      end

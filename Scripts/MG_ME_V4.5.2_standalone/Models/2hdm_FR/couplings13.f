      subroutine coup13(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 13
      if(readlha) then

      MGVX87 = -6*l1*PTH133*TH12 - 3*lR7*PTH233*TH12 - 3*l3*PTH232*
     +TH12*TH13 - 3*l4*PTH232*TH12*TH13 - 6*l5*PTH232*TH12*TH13 - 
     +3*l3*PTH332*TH12*TH13 - 3*l4*PTH332*TH12*TH13 + 6*l5*PTH332*
     +TH12*TH13 - 3*lR6*PTH133*TH22 - 6*l2*PTH233*TH22 - 9*lR7*PTH
     +232*TH13*TH22 - 3*lR7*PTH332*TH13*TH22 - 9*lR6*PTH132*TH12*T
     +H23 - 3*lR7*PTH332*TH12*TH23 - 3*l3*PTH132*TH22*TH23 - 3*l4*
     +PTH132*TH22*TH23 - 6*l5*PTH132*TH22*TH23 - 6*l2*PTH332*TH22*
     +TH23 - 6*l2*PTH333*TH32 - 3*l3*PTH132*TH32*TH33 - 3*l4*PTH13
     +2*TH32*TH33 + 6*l5*PTH132*TH32*TH33 - 6*l2*PTH232*TH32*TH33 
     +- 6*lR7*TH13*TH23*TH32*TH33
      MGVX73 = -(l3*TH12*TH13) - lR7*TH13*TH22 - lR7*TH12*TH23 - 2*
     +l2*TH22*TH23 - 2*l2*TH32*TH33
      MGVX88 = -6*l1*PTH134 - 6*l3*PTH132*PTH232 - 6*l4*PTH132*PTH2
     +32 - 12*l5*PTH132*PTH232 - 6*l2*PTH234 - 6*l3*PTH132*PTH332 
     +- 6*l4*PTH132*PTH332 + 12*l5*PTH132*PTH332 - 12*l2*PTH232*PT
     +H332 - 6*l2*PTH334 - 12*lR7*PTH233*TH13 - 12*lR6*PTH133*TH23
     + - 12*lR7*PTH332*TH13*TH23
      MGVX79 = -(l3*PTH132) - 2*l2*PTH232 - 2*l2*PTH332 - 2*lR7*TH1
     +3*TH23
      MGVX52 = -4*l2
      MGVX102 = (ee**2*TH11*v)/(2.*Psw2)
      MGVX123 = ee**2*TH11*v + (ee**2*Pcw2*TH11*v)/(2.*Psw2) + (ee*
     +*2*Psw2*TH11*v)/(2.*Pcw2)
      MGVX104 = (ee**2*TH12*v)/(2.*Psw2)
      MGVX125 = ee**2*TH12*v + (ee**2*Pcw2*TH12*v)/(2.*Psw2) + (ee*
     +*2*Psw2*TH12*v)/(2.*Pcw2)
      MGVX107 = (ee**2*TH13*v)/(2.*Psw2)
      MGVX128 = ee**2*TH13*v + (ee**2*Pcw2*TH13*v)/(2.*Psw2) + (ee*
     +*2*Psw2*TH13*v)/(2.*Pcw2)
      MGVX103 = (ee**2*PTH112)/(2.*Psw2) + (ee**2*PTH212)/(2.*Psw2)
     + + (ee**2*PTH312)/(2.*Psw2)
      MGVX124 = ee**2*PTH112 + (ee**2*Pcw2*PTH112)/(2.*Psw2) + (ee*
     +*2*Psw2*PTH112)/(2.*Pcw2) + ee**2*PTH212 + (ee**2*Pcw2*PTH21
     +2)/(2.*Psw2) + (ee**2*Psw2*PTH212)/(2.*Pcw2) + ee**2*PTH312 
     ++ (ee**2*Pcw2*PTH312)/(2.*Psw2) + (ee**2*Psw2*PTH312)/(2.*Pc
     +w2)
      MGVX105 = (ee**2*TH11*TH12)/(2.*Psw2) + (ee**2*TH21*TH22)/(2.
     +*Psw2) + (ee**2*TH31*TH32)/(2.*Psw2)
      MGVX126 = ee**2*TH11*TH12 + (ee**2*Pcw2*TH11*TH12)/(2.*Psw2) 
     ++ (ee**2*Psw2*TH11*TH12)/(2.*Pcw2) + ee**2*TH21*TH22 + (ee**
     +2*Pcw2*TH21*TH22)/(2.*Psw2) + (ee**2*Psw2*TH21*TH22)/(2.*Pcw
     +2) + ee**2*TH31*TH32 + (ee**2*Pcw2*TH31*TH32)/(2.*Psw2) + (e
     +e**2*Psw2*TH31*TH32)/(2.*Pcw2)

      endif

      return
      end

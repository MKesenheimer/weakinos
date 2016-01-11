      subroutine coup14(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 14
      if(readlha) then

      MGVX108 = (ee**2*TH11*TH13)/(2.*Psw2) + (ee**2*TH21*TH23)/(2.
     +*Psw2) + (ee**2*TH31*TH33)/(2.*Psw2)
      MGVX129 = ee**2*TH11*TH13 + (ee**2*Pcw2*TH11*TH13)/(2.*Psw2) 
     ++ (ee**2*Psw2*TH11*TH13)/(2.*Pcw2) + ee**2*TH21*TH23 + (ee**
     +2*Pcw2*TH21*TH23)/(2.*Psw2) + (ee**2*Psw2*TH21*TH23)/(2.*Pcw
     +2) + ee**2*TH31*TH33 + (ee**2*Pcw2*TH31*TH33)/(2.*Psw2) + (e
     +e**2*Psw2*TH31*TH33)/(2.*Pcw2)
      MGVX95 = (ee**2*TH21)/(2.*sw) - ((0,0.5)*ee**2*TH31)/sw
      MGVX119 = -(ee**2*TH21)/(2.*cw) + ((0,0.5)*ee**2*TH31)/cw
      MGVX89 = (ee**2*TH21)/(2.*sw) + ((0,0.5)*ee**2*TH31)/sw
      MGVX116 = -(ee**2*TH21)/(2.*cw) - ((0,0.5)*ee**2*TH31)/cw
      MGVX106 = (ee**2*PTH122)/(2.*Psw2) + (ee**2*PTH222)/(2.*Psw2)
     + + (ee**2*PTH322)/(2.*Psw2)
      MGVX127 = ee**2*PTH122 + (ee**2*Pcw2*PTH122)/(2.*Psw2) + (ee*
     +*2*Psw2*PTH122)/(2.*Pcw2) + ee**2*PTH222 + (ee**2*Pcw2*PTH22
     +2)/(2.*Psw2) + (ee**2*Psw2*PTH222)/(2.*Pcw2) + ee**2*PTH322 
     ++ (ee**2*Pcw2*PTH322)/(2.*Psw2) + (ee**2*Psw2*PTH322)/(2.*Pc
     +w2)
      MGVX109 = (ee**2*TH12*TH13)/(2.*Psw2) + (ee**2*TH22*TH23)/(2.
     +*Psw2) + (ee**2*TH32*TH33)/(2.*Psw2)
      MGVX130 = ee**2*TH12*TH13 + (ee**2*Pcw2*TH12*TH13)/(2.*Psw2) 
     ++ (ee**2*Psw2*TH12*TH13)/(2.*Pcw2) + ee**2*TH22*TH23 + (ee**
     +2*Pcw2*TH22*TH23)/(2.*Psw2) + (ee**2*Psw2*TH22*TH23)/(2.*Pcw
     +2) + ee**2*TH32*TH33 + (ee**2*Pcw2*TH32*TH33)/(2.*Psw2) + (e
     +e**2*Psw2*TH32*TH33)/(2.*Pcw2)
      MGVX97 = (ee**2*TH22)/(2.*sw) - ((0,0.5)*ee**2*TH32)/sw
      MGVX120 = -(ee**2*TH22)/(2.*cw) + ((0,0.5)*ee**2*TH32)/cw
      MGVX91 = (ee**2*TH22)/(2.*sw) + ((0,0.5)*ee**2*TH32)/sw
      MGVX117 = -(ee**2*TH22)/(2.*cw) - ((0,0.5)*ee**2*TH32)/cw
      MGVX110 = (ee**2*PTH132)/(2.*Psw2) + (ee**2*PTH232)/(2.*Psw2)
     + + (ee**2*PTH332)/(2.*Psw2)

      endif

      return
      end

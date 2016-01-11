      subroutine coup1(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 1
      MGVX1 = G
      MGVX2 = ggT1

      MGVX75 = (-2*ee*G)/3.
      MGVX89 = (4*ee*G)/3.
      MGVX331 = G**2
      MGVX240 = -((CONJRd13*ee*G*Rd13)/(cw*sw)) + (2*ee*G*sw)/(3.*c
     +w)
      MGVX247 = -((CONJRd13*ee*G*Rd23)/(cw*sw))
      MGVX132 = (CONJRd13*ee*G*Ru13*Sqrt2)/sw
      MGVX140 = (CONJRd13*ee*G*Ru63*Sqrt2)/sw
      MGVX246 = -((CONJRd23*ee*G*Rd13)/(cw*sw))
      MGVX164 = (CONJRu13*ee*G*Rd13*Sqrt2)/sw
      MGVX172 = (CONJRu63*ee*G*Rd13*Sqrt2)/sw
      MGVX241 = -((CONJRd23*ee*G*Rd23)/(cw*sw)) + (2*ee*G*sw)/(3.*c
     +w)
      MGVX133 = (CONJRd23*ee*G*Ru13*Sqrt2)/sw
      MGVX141 = (CONJRd23*ee*G*Ru63*Sqrt2)/sw
      MGVX165 = (CONJRu13*ee*G*Rd23*Sqrt2)/sw
      MGVX173 = (CONJRu63*ee*G*Rd23*Sqrt2)/sw
      MGVX242 = (2*ee*G*sw)/(3.*cw)
      MGVX244 = -((CONJRd51*ee*G*Rd51)/(cw*sw)) + (2*ee*G*sw)/(3.*c
     +w)
      MGVX137 = (CONJRd51*ee*G*Ru51*Sqrt2)/sw
      MGVX169 = (CONJRu51*ee*G*Rd51*Sqrt2)/sw
      MGVX245 = -((CONJRd62*ee*G*Rd62)/(cw*sw)) + (2*ee*G*sw)/(3.*c
     +w)
      MGVX135 = (CONJRd62*ee*G*Ru42*Sqrt2)/sw
      MGVX167 = (CONJRu42*ee*G*Rd62*Sqrt2)/sw
      MGVX252 = (CONJRu13*ee*G*Ru13)/(cw*sw) - (4*ee*G*sw)/(3.*cw)
      MGVX269 = (CONJRu13*ee*G*Ru63)/(cw*sw)
      MGVX267 = (CONJRu63*ee*G*Ru13)/(cw*sw)
      MGVX255 = (-4*ee*G*sw)/(3.*cw)
      MGVX261 = (CONJRu42*ee*G*Ru42)/(cw*sw) - (4*ee*G*sw)/(3.*cw)
      MGVX264 = (CONJRu51*ee*G*Ru51)/(cw*sw) - (4*ee*G*sw)/(3.*cw)
      MGVX271 = (CONJRu63*ee*G*Ru63)/(cw*sw) - (4*ee*G*sw)/(3.*cw)
      MGVX74 = G

      MGVX552(1) = CONJRd16*G*Sqrt2
      MGVX552(2) = -(CONJRd13*G*Sqrt2)
      MGVX546(1) = -(G*Rd13*Sqrt2)
      MGVX546(2) = G*Rd16*Sqrt2
      MGVX553(1) = CONJRd26*G*Sqrt2
      MGVX553(2) = -(CONJRd23*G*Sqrt2)
      MGVX547(1) = -(G*Rd23*Sqrt2)
      MGVX547(2) = G*Rd26*Sqrt2
      MGVX554(1) = CONJRd35*G*Sqrt2
      MGVX554(2) = 0
      MGVX548(1) = 0
      MGVX548(2) = G*Rd35*Sqrt2
      MGVX555(1) = CONJRd44*G*Sqrt2
      MGVX555(2) = 0
      MGVX549(1) = 0
      MGVX549(2) = G*Rd44*Sqrt2
      MGVX556(1) = 0
      MGVX556(2) = -(CONJRd51*G*Sqrt2)
      MGVX550(1) = -(G*Rd51*Sqrt2)
      MGVX550(2) = 0
      MGVX557(1) = 0
      MGVX557(2) = -(CONJRd62*G*Sqrt2)
      MGVX551(1) = -(G*Rd62*Sqrt2)
      MGVX551(2) = 0
      GG(1) = -G
      GG(2) = -G
      MGVX562(1) = CONJRu16*G*Sqrt2
      MGVX562(2) = -(CONJRu13*G*Sqrt2)
      MGVX569(1) = -(G*Ru13*Sqrt2)
      MGVX569(2) = G*Ru16*Sqrt2
      MGVX574(1) = CONJRu24*G*Sqrt2
      MGVX574(2) = 0
      MGVX579(1) = 0
      MGVX579(2) = G*Ru24*Sqrt2
      MGVX584(1) = CONJRu35*G*Sqrt2
      MGVX584(2) = 0
      MGVX589(1) = 0
      MGVX589(2) = G*Ru35*Sqrt2
      MGVX594(1) = 0
      MGVX594(2) = -(CONJRu42*G*Sqrt2)
      MGVX601(1) = -(G*Ru42*Sqrt2)
      MGVX601(2) = 0
      MGVX606(1) = 0
      MGVX606(2) = -(CONJRu51*G*Sqrt2)
      MGVX613(1) = -(G*Ru51*Sqrt2)
      MGVX613(2) = 0
      MGVX618(1) = CONJRu66*G*Sqrt2
      MGVX618(2) = -(CONJRu63*G*Sqrt2)
      MGVX625(1) = -(G*Ru63*Sqrt2)
      MGVX625(2) = G*Ru66*Sqrt2

      return
      end

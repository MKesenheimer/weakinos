      subroutine coup25(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 25
      if(readlha) then

      MGVX753 = -(Cosbet*ee**2*Sinalp*Sinbet*vd)/(4.*Pcw2*Psw2) + (
     +Cosalp*ee**2*vd*Cos(bet)**2)/(4.*Pcw2*Psw2) - (Cosalp*ee**2*
     +vd*Sin(bet)**2)/(4.*Pcw2*Psw2) + (ee**2*Sinalp*Tanbet*vd*Sin
     +(bet)**2)/(4.*Pcw2*Psw2)
      MGVX688 = ((0,1)*CONJRd16*CONJtd33*Rd13*Sinbet)/Sqrt2 - ((0,1
     +)*CONJRd13*Rd16*Sinbet*td33)/Sqrt2 + ((0,1)*CONJRd16*Cosbet*
     +MMMU*Rd13*yukd3)/Sqrt2 - ((0,1)*CONJMMMU*CONJRd13*Cosbet*Rd1
     +6*yukd3)/Sqrt2
      MGVX692 = ((0,1)*CONJRd16*CONJtd33*Rd23*Sinbet)/Sqrt2 - ((0,1
     +)*CONJRd13*Rd26*Sinbet*td33)/Sqrt2 + ((0,1)*CONJRd16*Cosbet*
     +MMMU*Rd23*yukd3)/Sqrt2 - ((0,1)*CONJMMMU*CONJRd13*Cosbet*Rd2
     +6*yukd3)/Sqrt2
      MGVX690 = ((0,1)*CONJRd26*CONJtd33*Rd13*Sinbet)/Sqrt2 - ((0,1
     +)*CONJRd23*Rd16*Sinbet*td33)/Sqrt2 + ((0,1)*CONJRd26*Cosbet*
     +MMMU*Rd13*yukd3)/Sqrt2 - ((0,1)*CONJMMMU*CONJRd23*Cosbet*Rd1
     +6*yukd3)/Sqrt2
      MGVX694 = ((0,1)*CONJRd26*CONJtd33*Rd23*Sinbet)/Sqrt2 - ((0,1
     +)*CONJRd23*Rd26*Sinbet*td33)/Sqrt2 + ((0,1)*CONJRd26*Cosbet*
     +MMMU*Rd23*yukd3)/Sqrt2 - ((0,1)*CONJMMMU*CONJRd23*Cosbet*Rd2
     +6*yukd3)/Sqrt2
      MGVX714 = ((0,1)*CONJRl16*CONJte33*Rl13*Sinbet)/Sqrt2 - ((0,1
     +)*CONJRl13*Rl16*Sinbet*te33)/Sqrt2 + ((0,1)*CONJRl16*Cosbet*
     +MMMU*Rl13*yukl3)/Sqrt2 - ((0,1)*CONJMMMU*CONJRl13*Cosbet*Rl1
     +6*yukl3)/Sqrt2
      MGVX719 = ((0,1)*CONJRl16*CONJte33*Rl63*Sinbet)/Sqrt2 - ((0,1
     +)*CONJRl13*Rl66*Sinbet*te33)/Sqrt2 + ((0,1)*CONJRl16*Cosbet*
     +MMMU*Rl63*yukl3)/Sqrt2 - ((0,1)*CONJMMMU*CONJRl13*Cosbet*Rl6
     +6*yukl3)/Sqrt2
      MGVX715 = ((0,1)*CONJRl66*CONJte33*Rl13*Sinbet)/Sqrt2 - ((0,1
     +)*CONJRl63*Rl16*Sinbet*te33)/Sqrt2 + ((0,1)*CONJRl66*Cosbet*
     +MMMU*Rl13*yukl3)/Sqrt2 - ((0,1)*CONJMMMU*CONJRl63*Cosbet*Rl1
     +6*yukl3)/Sqrt2
      MGVX720 = ((0,1)*CONJRl66*CONJte33*Rl63*Sinbet)/Sqrt2 - ((0,1
     +)*CONJRl63*Rl66*Sinbet*te33)/Sqrt2 + ((0,1)*CONJRl66*Cosbet*
     +MMMU*Rl63*yukl3)/Sqrt2 - ((0,1)*CONJMMMU*CONJRl63*Cosbet*Rl6
     +6*yukl3)/Sqrt2
      MGVX737 = ((0,1)*CONJRu16*CONJtu33*Cosbet*Ru13)/Sqrt2 - ((0,1
     +)*CONJRu13*Cosbet*Ru16*tu33)/Sqrt2 + ((0,1)*CONJRu16*MMMU*Ru
     +13*Sinbet*yuku3)/Sqrt2 - ((0,1)*CONJMMMU*CONJRu13*Ru16*Sinbe
     +t*yuku3)/Sqrt2
      MGVX743 = ((0,1)*CONJRu16*CONJtu33*Cosbet*Ru63)/Sqrt2 - ((0,1
     +)*CONJRu13*Cosbet*Ru66*tu33)/Sqrt2 + ((0,1)*CONJRu16*MMMU*Ru
     +63*Sinbet*yuku3)/Sqrt2 - ((0,1)*CONJMMMU*CONJRu13*Ru66*Sinbe
     +t*yuku3)/Sqrt2
      MGVX740 = ((0,1)*CONJRu66*CONJtu33*Cosbet*Ru13)/Sqrt2 - ((0,1
     +)*CONJRu63*Cosbet*Ru16*tu33)/Sqrt2 + ((0,1)*CONJRu66*MMMU*Ru
     +13*Sinbet*yuku3)/Sqrt2 - ((0,1)*CONJMMMU*CONJRu63*Ru16*Sinbe
     +t*yuku3)/Sqrt2
      MGVX746 = ((0,1)*CONJRu66*CONJtu33*Cosbet*Ru63)/Sqrt2 - ((0,1
     +)*CONJRu63*Cosbet*Ru66*tu33)/Sqrt2 + ((0,1)*CONJRu66*MMMU*Ru
     +63*Sinbet*yuku3)/Sqrt2 - ((0,1)*CONJMMMU*CONJRu63*Ru66*Sinbe
     +t*yuku3)/Sqrt2
      MGVX705 = CONJRd16*CONJtd33*Ru13*Sinbet + CONJRd13*Cosbet*Ru1
     +6*tu33 - (CONJRd13*ee**2*Ru13*Sinbet*vd)/(Psw2*Sqrt2) + (CON
     +JRd13*Pyukd32*Ru13*Sinbet*vd)/Sqrt2 + (CONJRd13*Pyuku32*Ru13
     +*Sinbet*vd)/Sqrt2 + CONJRd16*Cosbet*MMMU*Ru13*yukd3 + CONJMM
     +MU*CONJRd13*Ru16*Sinbet*yuku3 + (CONJRd16*Cosbet*Ru16*vd*yuk
     +d3*yuku3)/Sqrt2 + (CONJRd16*Ru16*Sinbet*Tanbet*vd*yukd3*yuku
     +3)/Sqrt2
      MGVX709 = CONJRd16*CONJtd33*Ru63*Sinbet + CONJRd13*Cosbet*Ru6
     +6*tu33 - (CONJRd13*ee**2*Ru63*Sinbet*vd)/(Psw2*Sqrt2) + (CON
     +JRd13*Pyukd32*Ru63*Sinbet*vd)/Sqrt2 + (CONJRd13*Pyuku32*Ru63
     +*Sinbet*vd)/Sqrt2 + CONJRd16*Cosbet*MMMU*Ru63*yukd3 + CONJMM
     +MU*CONJRd13*Ru66*Sinbet*yuku3 + (CONJRd16*Cosbet*Ru66*vd*yuk
     +d3*yuku3)/Sqrt2 + (CONJRd16*Ru66*Sinbet*Tanbet*vd*yukd3*yuku
     +3)/Sqrt2

      endif

      return
      end

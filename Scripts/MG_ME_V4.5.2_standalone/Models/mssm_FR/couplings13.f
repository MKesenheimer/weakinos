      subroutine coup13(readlha)

      implicit none
      logical readlha

      include 'input.inc'
      include 'coupl.inc'
      include 'intparam_definition.inc'


c Interactions associated with 13
      if(readlha) then


      MGVX644(1) = CONJRLn13*CONJUU12*yukl3
      MGVX644(2) = -((CONJRLn13*ee*VV11)/sw)
      MGVX649(1) = 0
      MGVX649(2) = -((CONJRLn31*ee*VV21)/sw)
      MGVX650(1) = 0
      MGVX650(2) = -((CONJRLn22*ee*VV21)/sw)
      MGVX651(1) = CONJRLn13*CONJUU22*yukl3
      MGVX651(2) = -((CONJRLn13*ee*VV21)/sw)
      MGVX656(1) = (CONJOL11*ee)/sw
      MGVX656(2) = (CONJOR11*ee)/sw
      MGVX657(1) = (CONJOL21*ee)/sw
      MGVX657(2) = (CONJOR21*ee)/sw
      MGVX658(1) = (CONJOL31*ee)/sw
      MGVX658(2) = (CONJOR31*ee)/sw
      MGVX659(1) = (CONJOL41*ee)/sw
      MGVX659(2) = (CONJOR41*ee)/sw
      MGVX660(1) = (CONJOL12*ee)/sw
      MGVX660(2) = (CONJOR12*ee)/sw
      MGVX661(1) = (CONJOL22*ee)/sw
      MGVX661(2) = (CONJOR22*ee)/sw
      MGVX662(1) = (CONJOL32*ee)/sw
      MGVX662(2) = (CONJOR32*ee)/sw
      MGVX663(1) = (CONJOL42*ee)/sw
      MGVX663(2) = (CONJOR42*ee)/sw
      MGVX664(1) = (ee*OL11)/sw
      MGVX664(2) = (ee*OR11)/sw
      MGVX665(1) = (ee*OL12)/sw
      MGVX665(2) = (ee*OR12)/sw
      MGVX666(1) = (ee*OL21)/sw
      MGVX666(2) = (ee*OR21)/sw
      endif

      return
      end

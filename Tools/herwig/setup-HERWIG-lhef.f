c     opens input file and counts number of events, setting maxev
      subroutine getmaxev(maxev)
      implicit none
      integer maxev
      call opencount(maxev)
      end
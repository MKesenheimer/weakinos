      subroutine newunit(iun)
      implicit none
      integer iun
      logical ok
      integer j
c units 97 and 99 are used for lhe files;
c keep reserved.
      do j=20,96 ! MK: changed 10 -> 20. SLHAlib uses unit=10 to read from SLHA files
         inquire(unit=j,opened=ok)
         if(.not.ok) then
            iun=j
            return
         endif
      enddo
      end

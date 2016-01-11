* timing.F
* this program reads and writes an SLHA file 1000 times
*   assuming that actual disc I/O happens only the first time,
*   the timing of this program may be used to obtain an estimate
*   for the actual performance of the SLHALib
* last modified 17 Jan 09 th


	program SLHATiming
	implicit none

#include "SLHA.h"

#define Strip(s) s(1:lnblnk(s))

	integer error, repeat
	double complex slhadata(nslhadata)

	character*80 argv

#if U77EXT
	integer iargc, lnblnk
	external iargc
#endif

	if( iargc() .ne. 1 ) then
	  call getarg(0, argv)
	  print *, "USAGE:"
	  print *, Strip(argv), " FILE"
	  print *, "reads FILE and writes it out again as FILE.out"
	  stop
	endif

	call getarg(1, argv)

	do repeat = 1, 1000

	call SLHARead(error, slhadata, argv, 1)
	if( error .ne. 0 ) stop "Read error"

	call SLHAWrite(error, slhadata, Strip(argv)//".out")
	if( error .ne. 0 ) stop "Write error"

	enddo
	end


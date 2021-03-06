*#define CHECK

#if 0
* uubar_x2x2_SquaredME.F
* assembly of squared matrix element
* generated by FormCalc 8.4 on 18-Jun-2015 10:15
#endif

#include "uubar_x2x2_vars.h"
#include "contains.h"


************************************************************************

	subroutine uubar_x2x2_SquaredMEHel(res, flags)
	implicit none
	ResType res(HelDim(*))
	integer flags

#include "uubar_x2x2_vars.h"

        ! MK:
        HelType amp, ampctr
        ResType ampsq, ampsqctr
        
        ! MK: t- and s-channel amplitudes without a factor 4*Pi*Alfa
        ! needed for susy restoring
        ComplexType tchannel, schannel

* BEGIN ABBR_HEL
	call uubar_x2x2_abbr0h
	TEST(flags, BIT_LOOP)
	call uubar_x2x2_abbr1h
	ENDTEST(flags, BIT_LOOP)
* END ABBR_HEL

        ! MK:
        ampsq = 0
        ampsqctr = 0
* BEGIN FF_TREE
	Ctree(HelInd(1)) = 0

	call uubar_x2x2_born
	
        ! MK: copy the s-channel diagrams but factor out 4*Pi*Alfa
        schannel = Ctree(HelInd(1))/(4*Pi*Alfa)

	LOOP(Sfe5, 1,2,1)
	call uubar_x2x2_born_Sfe5
	ENDLOOP(Sfe5)
	
        ! MK: extract the t-channel diagrams (without 4*Pi*Alfa)
        tchannel = Ctree(HelInd(1))/(4*Pi*Alfa) - schannel
* END FF_TREE

* MK: calculate the SUSY-restoring counter term        
        ampctr = MatSUN(1,1)*(schannel+tchannel)
        ! put back the factors and square the matrix element
        ! just keep terms proportional two Alfa*AlfaS
        ampsqctr = 16D0*Pi*Alfa*El*dZe1y*Re(Conjugate(ampctr)*tchannel)
* MK: end calculation of SUSY-restoring counter term

* BEGIN M2_TREE
	amp = 0
	amp = amp + MatSUN(1,1)*Ctree(HelInd(1))
	ampsq = ampsq + Re(Conjugate(Ctree(HelInd(1)))*amp)
* END M2_TREE
	res(HelInd(1)) = ampsq

	ampsq = 0
	TEST(flags, BIT_LOOP)
* BEGIN FF_LOOP
	Cloop(HelInd(1)) = 0

	call uubar_x2x2_vert

	LOOP(Sfe5, 1,2,1)
	call uubar_x2x2_self_Sfe5
	call uubar_x2x2_vert_Sfe5
	call uubar_x2x2_box_Sfe5
	ENDLOOP(Sfe5)

	LOOP(Sfe6, 1,2,1)
	LOOP(Sfe5, 1,2,1)
	call uubar_x2x2_self_Sfe5_Sfe6
	call uubar_x2x2_vert_Sfe5_Sfe6
	call uubar_x2x2_box_Sfe5_Sfe6
	ENDLOOP(Sfe5)
	ENDLOOP(Sfe6)

	LOOP(Sfe7, 1,2,1)
	LOOP(Sfe6, 1,2,1)
	LOOP(Sfe5, 1,2,1)
	call uubar_x2x2_self_Sfe5_Sfe6_Sfe7
	ENDLOOP(Sfe5)
	ENDLOOP(Sfe6)
	ENDLOOP(Sfe7)
* END FF_LOOP

* BEGIN M2_LOOP
	amp = 0
	amp = amp + MatSUN(1,1)*Cloop(HelInd(1))
	ampsq = ampsq + Re(Conjugate(Ctree(HelInd(1)))*amp)
	ampsq = ampsq + ampsq
* END M2_LOOP
	ENDTEST(flags, BIT_LOOP)
        ! MK:
        res(HelInd(2)) = ampsq + ampsqctr

	end

************************************************************************

	subroutine uubar_x2x2_SquaredME(result, helicities, flags)
	implicit none
	RealType result(*)
	integer*8 helicities
	integer flags

#include "uubar_x2x2_vars.h"

* BEGIN VAR_DECL
	SIMD_ONLY(integer v)
	integer i, h, hmax, hsimd
	parameter (hmax = 16)
	parameter (hsimd = SIMD_CEIL(hmax))
	integer*8 hbits
	ResType res(HelDim(2),hsimd)
	RealType rtree, rloop
	external uubar_x2x2_SquaredMEHel
* END VAR_DECL

* BEGIN HSEL_DECL
	RealType norm
	ResType hseltest_v(HelDim(hsimd))
	RealType hseltest_s(0:hmax-1)
	equivalence (hseltest_v, hseltest_s)
	RealType hselmin
	integer hseli
	save hseltest_v, hselmin, hseli
* END HSEL_DECL

	data MatSUN /1*bogus/


	PAR_PREP(res,res(HelInd(1),2),Sub1,ends,Sub111,enda,seq,endhel)

* BEGIN INVARIANTS
	S = SInvariant(1,2)
	T = TInvariant(1,3)
	U = TInvariant(2,3)
* END INVARIANTS

	TEST(flags, BIT_RESET)
* BEGIN ABBR_S
	seq(1) = seq(1) + 1
	INI_S(seq)
	call uubar_x2x2_abbr0s
	TEST(flags, BIT_LOOP)
	call uubar_x2x2_abbr1s
	ENDTEST(flags, BIT_LOOP)
* END ABBR_S
	ENDTEST(flags, BIT_RESET)

* BEGIN ABBR_ANGLE
	seq(2) = seq(2) + 1
	INI_ANGLE(seq)
	TEST(flags, BIT_LOOP)
	call uubar_x2x2_abbr1a
	ENDTEST(flags, BIT_LOOP)
* END ABBR_ANGLE

* BEGIN HEL_LOOP
	SIMD_ONLY(v = 1)
	h = 0

	do i = 0, hmax - 1
* BEGIN HSEL_IF
	  if( hseltest_s(i) .lt. hselmin ) cycle
* END HSEL_IF

	  Hel(1) = -1 + 2*Mod(i/8,2)
	  Hel(2) = -1 + 2*Mod(i/4,2)
	  Hel(3) = -1 + 2*Mod(i/2,2)
	  Hel(4) = -1 + 2*Mod(i/1,2)
	  hbits =
     &      + MASK_HEL(1)
     &      + MASK_HEL(2)
     &      + MASK_HEL(3)
     &      + MASK_HEL(4)
	  if( iand(helicities, hbits) .ne. hbits ) cycle

	  SIMD_ONLY(call VecCopy(v, LEGS, Hel(1)))
	  SIMD_MULT(v = mod(v, SIMD) + 1)
	  SIMD_MULT(if( v .eq. 1 ) then)
	  h = h + 1
	  PAR_EXEC(uubar_x2x2_SquaredMEHel, res(HelInd(1),h), flags)
	  SIMD_MULT(endif)
	enddo

	SIMD_MULT(if( v .ne. 1 ) then)
	SIMD_MULT(h = h + 1)
	SIMD_MULT(PAR_EXEC(uubar_x2x2_SquaredMEHel,
     &                     res(HelInd(1),h), flags))
	SIMD_MULT(endif)
* END HEL_LOOP

	PAR_SYNC()
	DEINI(seq)

	SIMD_MULT(if( v .ne. 1 ) res(v:SIMD,:,h) = 0)
* BEGIN RESULT
	rtree = 0
	rloop = 0
	do i = 1, h
	  rtree = rtree + HelSum(res(HelInd(1),i))
	  rloop = rloop + HelSum(res(HelInd(2),i))
	enddo
	result(1) = rtree
	TEST(flags, BIT_LOOP)
#ifdef WF_RENORMALIZATION
	rloop = rloop + (WF_RENORMALIZATION)*rtree
#endif
	result(2) = rloop
	ENDTEST(flags, BIT_LOOP)
* END RESULT

* BEGIN HSEL_SET
	if( hseli .lt. hseln ) then
	  norm = 1/(rtree + rloop)
	  do i = 1, hsimd
	    hseltest_v(HelInd(i)) = hseltest_v(HelInd(i)) +
     &        abs(norm*(res(HelInd(1),i) + res(HelInd(2),i)))
	  enddo
	  hseli = hseli + 1
	  if( hseli .eq. hseln ) then
	    hselmin = 0
	    do i = 0, hmax - 1
	      hselmin = max(hselmin, hseltest_s(i))
	    enddo
	    hselmin = hselmin*hseleps
	    do i = 0, hmax - 1
	      if( hseltest_s(i) .ge. hselmin ) cycle
	      Hel(1) = -1 + 2*Mod(i/8,2)
	      Hel(2) = -1 + 2*Mod(i/4,2)
	      Hel(3) = -1 + 2*Mod(i/2,2)
	      Hel(4) = -1 + 2*Mod(i/1,2)
	      print 123, Hel
123	      format(' neglecting ', LEGS I3)
	    enddo
	  endif
	endif
* END HSEL_SET

#ifdef CHECK
	print *, 'S =', S
	print *, 'T =', T
	print *, 'U =', U
	print *, 'tree =', rtree
	print *, 'loop =', rloop
	stop
#endif

#include "contains.h"
* END SQUAREDME
	end


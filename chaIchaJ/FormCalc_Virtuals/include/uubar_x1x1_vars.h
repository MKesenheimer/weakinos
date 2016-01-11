#if 0
* vars.h
* variable declarations
* generated by FormCalc 8.4 on 12-Jun-2015 10:14
#endif

#ifndef VARS_H
#define VARS_H

#define LEGS 4

#include "decl.h"

#else

#include "decl.h"

	ComplexType Sub1(2), Sub2(2), Sub4, Sub5
	ComplexType Sub31, Sub89, Sub90, Sub49, Sub67, Sub74, Sub80
	ComplexType Sub30, Sub109, Sub75, Sub81, Sub50, Sub68, Sub66
	ComplexType Sub63, Sub79, Sub48, Sub53, Sub73, Sub70, Sub64
	ComplexType Sub69, Sub55, Sub56, Sub77, Sub76, Sub83
	ComplexType Sub72(2), Sub62(2), Sub58(2), Sub59(2), Sub52(2)
	ComplexType Sub47(2), Sub54(2), Sub51(2), Sub57(2), Sub60(2)
	ComplexType Sub61(2), Sub82(2), Sub65(2), Sub71(2), Sub78(2)
	ComplexType Sub99(2), Sub84(2), Sub9(2), Sub42(2), Sub43(2)
	ComplexType Sub8(2,2), Sub37(2,2), Sub20(2,2), Sub45(2,2)
	ComplexType Sub21(2,2), Sub22(2,2), Sub28(2,2,2), Sub35
	ComplexType Sub32, Sub110, Sub108, Sub34, Sub33, Sub36
	ComplexType Sub115(2), Sub7(2,2), Sub29, Sub92, pave2(2)
	ComplexType pave3, pave7(Nbb), pave8(Ncc), pave10(Ncc,2)
	ComplexType pave12(Ncc,2,2)
	common /varXs/ Sub1, Sub2, Sub4, Sub5, Sub31, Sub89, Sub90
	common /varXs/ Sub49, Sub67, Sub74, Sub80, Sub30, Sub109
	common /varXs/ Sub75, Sub81, Sub50, Sub68, Sub66, Sub63
	common /varXs/ Sub79, Sub48, Sub53, Sub73, Sub70, Sub64
	common /varXs/ Sub69, Sub55, Sub56, Sub77, Sub76, Sub83
	common /varXs/ Sub72, Sub62, Sub58, Sub59, Sub52, Sub47
	common /varXs/ Sub54, Sub51, Sub57, Sub60, Sub61, Sub82
	common /varXs/ Sub65, Sub71, Sub78, Sub99, Sub84, Sub9
	common /varXs/ Sub42, Sub43, Sub8, Sub37, Sub20, Sub45
	common /varXs/ Sub21, Sub22, Sub28, Sub35, Sub32, Sub110
	common /varXs/ Sub108, Sub34, Sub33, Sub36, Sub115, Sub7
	common /varXs/ Sub29, Sub92, pave2, pave3, pave7, pave8
	common /varXs/ pave10, pave12

	ComplexType Sub111(2), Sub112(2), Sub113(2), Sub123(2)
	ComplexType Sub103(2), Sub102(2), Sub122(2), pave1(Nbb,2)
	ComplexType pave4(Nbb), pave5(Ndd,2,2), pave6(Ndd,2)
	ComplexType pave9(Ncc,2), pave11(Ncc,2)
	RealType S, T, U
	common /varXa/ Sub111, Sub112, Sub113, Sub123, Sub103, Sub102
	common /varXa/ Sub122, pave1, pave4, pave5, pave6, pave9
	common /varXa/ pave11, S, T, U

	HelType F9, F7, F1, F10, F3, F4, F8, F2, F6, F5, Sub6, Sub3
	HelType F11, F18, F15, F13, F14, F16, F17, F21, F12, F22, F23
	HelType F19, F25, F20, F24, Sub88, Sub87, Sub91, Sub97
	HelType Sub95, Sub85(HelDim(2)), Sub86(HelDim(2))
	HelType Sub114(HelDim(2)), Sub116(HelDim(2))
	HelType Sub96(HelDim(2)), Sub41(HelDim(2))
	HelType Sub105(HelDim(2),2), Sub119(HelDim(2),2)
	HelType Sub106(HelDim(2),2), Sub120(HelDim(2),2)
	HelType Sub121(HelDim(2),2), Sub40(HelDim(2),2)
	HelType Sub46(HelDim(2),2), Sub104(HelDim(2),2)
	HelType Sub107(HelDim(2),2), Sub94(HelDim(2),2)
	HelType Sub100(HelDim(2),2), Sub38(HelDim(2),2)
	HelType Sub39(HelDim(2),2), Sub44(HelDim(2),2), Sub93, Sub98
	HelType Sub124(HelDim(2)), Sub117(HelDim(2),2)
	HelType Sub101(HelDim(2),2), Sub118(HelDim(2),2)
	common /varXh/ F9, F7, F1, F10, F3, F4, F8, F2, F6, F5, Sub6
	common /varXh/ Sub3, F11, F18, F15, F13, F14, F16, F17, F21
	common /varXh/ F12, F22, F23, F19, F25, F20, F24, Sub88
	common /varXh/ Sub87, Sub91, Sub97, Sub95, Sub85, Sub86
	common /varXh/ Sub114, Sub116, Sub96, Sub41, Sub105, Sub119
	common /varXh/ Sub106, Sub120, Sub121, Sub40, Sub46, Sub104
	common /varXh/ Sub107, Sub94, Sub100, Sub38, Sub39, Sub44
	common /varXh/ Sub93, Sub98, Sub124, Sub117, Sub101, Sub118

	integer seq(2), Hel(4)
	common /helind/ seq, Hel

	integer Sfe5, Sfe6, Sfe7
	common /indices/ Sfe5, Sfe6, Sfe7

	HelType Ctree(HelDim(1)), Cloop(HelDim(1))
	ComplexType MatSUN(1,1)
	common /uubar_x1x1_formfactors/ Ctree, Cloop, MatSUN

#if PARALLEL
	marker ends, enda, endhel
	common /varXs/ ends
	common /varXa/ enda
	common /helind/ endhel
#endif

#endif

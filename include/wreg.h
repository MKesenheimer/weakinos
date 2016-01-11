c############### wreg.h ################################################
c last modified by MK, date

#ifndef WREG_H
#define WREG_H

#if defined(DSUB_I) || defined(DSUB_II) || defined(DSUB_II_TEST)
      double precision WREG  ! regulator to suppress the resonant diagrams
      parameter (WREG=1D0)   ! in a gauge-invariant way (default 1D0)      
#define WDLR WREG
#define WDRR WREG
#define WULR WREG
#define WURR WREG

#else
c if no subtraction scheme is defined, reset to physical widths (default zero)
      double precision WREG  ! define for Reals_osres.f here regulator too,
      parameter (WREG=0D0)   ! but set it to zero (this should not be modified)
#define WDLR WDL
#define WDRR WDR
#define WULR WUL
#define WURR WUR
#endif

#endif

c############### end wreg.h ############################################
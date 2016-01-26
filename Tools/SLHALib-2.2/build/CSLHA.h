#ifndef FTYPES_H
#define FTYPES_H

#if 1
#define FORTRAN(s) s##_
#else
#define FORTRAN(s) s
#endif

typedef int INTEGER;
typedef const INTEGER CINTEGER;
typedef double DOUBLE_PRECISION;
typedef const DOUBLE_PRECISION CDOUBLE_PRECISION;
typedef struct { DOUBLE_PRECISION re, im; } DOUBLE_COMPLEX;
typedef const DOUBLE_COMPLEX CDOUBLE_COMPLEX;
typedef char CHARACTER;
typedef const CHARACTER CCHARACTER;

#ifdef __cplusplus

#include <complex>
typedef std::complex<double> double_complex;
#define ToComplex(c) double_complex(c.re, c.im)
#define ToComplex2(r,i) double_complex(r, i)
#define Re(x) std::real(x)
#define Im(x) std::imag(x)

#else

typedef DOUBLE_COMPLEX double_complex;
#define ToComplex(c) c
#define ToComplex2(r,i) (double_complex){r, i}
#define Re(x) (x).re
#define Im(x) (x).im

#endif

#endif

/*
	CSLHA.h
		C/C++ wrapper functions for the SLHALib subroutines
		last modified 27 Oct 09 th
*/


#ifndef CSLHA_H
#define CSLHA_H

#include <string.h>
#include <stdarg.h>
#include "SLHADefs.h"
#include "PDG.h"

#define Slhadata(i) slhadata[i-1]
#define SlhaData(i) slhadata[i-1]


#ifdef __cplusplus
extern "C" {
#endif

extern void FORTRAN(slhaclear)(DOUBLE_COMPLEX *slhadata);

extern void FORTRAN(slharead)(INTEGER *error,
  DOUBLE_COMPLEX *slhadata, const char *filename,
  CINTEGER *abort, const int filename_len);

extern void FORTRAN(slhawrite)(INTEGER *error,
  CDOUBLE_COMPLEX *slhadata, const char *filename,
  const int filename_len);

extern void FORTRAN(slhaputinfo)(DOUBLE_COMPLEX *slhablock,
  CINTEGER *code, const char *text, const int text_len);

extern void FORTRAN(slhagetinfo)(DOUBLE_COMPLEX *slhablock,
  CINTEGER *i, char *text, const int text_len);

extern INTEGER FORTRAN(slhanewdecay)(DOUBLE_COMPLEX *slhadata,
  CDOUBLE_PRECISION *width, CINTEGER *parent_id);

extern INTEGER FORTRAN(slhafinddecay)(DOUBLE_COMPLEX *slhadata,
  CINTEGER *parent_id);

extern void FORTRAN(slhaadddecay)(DOUBLE_COMPLEX *slhadata,
  CDOUBLE_PRECISION *br, CINTEGER *decay, CINTEGER *nchildren,
  CINTEGER *child1_id, CINTEGER *child2_id,
  CINTEGER *child3_id, CINTEGER *child4_id);

extern DOUBLE_PRECISION FORTRAN(slhagetdecay)(
  CDOUBLE_COMPLEX *slhadata,
  CINTEGER *parent_id, CINTEGER *nchildren,
  CINTEGER *child1_id, CINTEGER *child2_id,
  CINTEGER *child3_id, CINTEGER *child4_id);

extern INTEGER FORTRAN(slhadecaytable)(CDOUBLE_COMPLEX *slhadata,
  CINTEGER *parent_id, DOUBLE_PRECISION *width, INTEGER *id,
  CINTEGER *maxparticles, CINTEGER *maxchannels);

extern INTEGER FORTRAN(slhaexist)(CDOUBLE_COMPLEX *slhablock,
  CINTEGER *length);

extern INTEGER FORTRAN(slhavalid)(CDOUBLE_COMPLEX *slhablock,
  CINTEGER *length);

extern void FORTRAN(slhapdgname)(CINTEGER *code,
  const char *name, const int name_len);

#ifdef __cplusplus
}
#endif


static inline void SLHAClear(double_complex *slhadata)
{
  FORTRAN(slhaclear)((DOUBLE_COMPLEX *)slhadata);
}


static inline void SLHARead(int *error,
  double_complex *slhadata, const char *filename, const int abort)
{
  FORTRAN(slharead)(error, (DOUBLE_COMPLEX *)slhadata,
    filename, &abort, strlen(filename));
}


static inline void SLHAWrite(int *error,
  const double_complex *slhadata, const char *filename)
{
  FORTRAN(slhawrite)(error, (CDOUBLE_COMPLEX *)slhadata,
    filename, strlen(filename));
}


static inline void SLHAPutInfo(double_complex *slhablock,
  const int code, const char *text)
{
  FORTRAN(slhaputinfo)((DOUBLE_COMPLEX *)slhablock, &code,
    text, strlen(text));
}


static inline void SLHAGetInfo(double_complex *slhablock,
  const int i, char *text, const int textsize)
{
  int off;

  FORTRAN(slhagetinfo)((DOUBLE_COMPLEX *)slhablock, &i,
    text, textsize);

  off = textsize - 1;
  while( off && text[off - 1] == ' ' ) --off;
  text[off] = 0;
}


static inline int SLHANewDecay(double_complex *slhadata,
  const double width, const int parent_id)
{
  return FORTRAN(slhanewdecay)((DOUBLE_COMPLEX *)slhadata,
    &width, &parent_id);
}


static inline int SLHAFindDecay(double_complex *slhadata,
  const int parent_id)
{
  return FORTRAN(slhafinddecay)((DOUBLE_COMPLEX *)slhadata,
    &parent_id);
}


static inline void SLHAAddDecay(double_complex *slhadata,
  const double br, const int decay,
  const int nchildren, const int child1_id, ...)
{
  va_list child_ids;
  int i, child_id[4];

  va_start(child_ids, child1_id);
  for( i = 0; i < nchildren; ++i )
    child_id[i] = va_arg(child_ids, int);
  va_end(child_ids);

  FORTRAN(slhaadddecay)((DOUBLE_COMPLEX *)slhadata,
    &br, &decay, &nchildren,
    &child_id[0], &child_id[1], &child_id[2], &child_id[3]);
}


static inline void SLHAGetDecay(double_complex *slhadata,
  const int parent_id,
  const int nchildren, const int child1_id, ...)
{
  va_list child_ids;
  int i, child_id[4];

  va_start(child_ids, child1_id);
  for( i = 0; i < nchildren; ++i )
    child_id[i] = va_arg(child_ids, int);
  va_end(child_ids);

  FORTRAN(slhagetdecay)((DOUBLE_COMPLEX *)slhadata,
    &parent_id, &nchildren,
    &child_id[0], &child_id[1], &child_id[2], &child_id[3]);
}


static inline int SLHADecayTable(const double_complex *slhadata,
  const int parent_id, double *width, int *id,
  const int maxparticles, const int maxchannels)
{
  return FORTRAN(slhadecaytable)((CDOUBLE_COMPLEX *)slhadata,
    &parent_id, width, id, &maxparticles, &maxchannels);
}


static inline int SLHAExist(const double_complex *slhablock,
  const int length)
{
  return FORTRAN(slhaexist)((CDOUBLE_COMPLEX *)slhablock, &length);
}


static inline int SLHAValid(const double_complex *slhablock,
  const int length)
{
  return FORTRAN(slhavalid)((CDOUBLE_COMPLEX *)slhablock, &length);
}


static inline void SLHAPDGName(const int code,
  char name[PDGLen+1])
{
  char *n;

  FORTRAN(slhapdgname)(&code, name, PDGLen);
  n = (char *)memchr(name, ' ', PDGLen);
  *(n ? n : name + PDGLen) = 0;
}

#endif


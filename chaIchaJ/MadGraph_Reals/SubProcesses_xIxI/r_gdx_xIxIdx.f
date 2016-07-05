      SUBROUTINE SMATRIX_GDX_XIXIDX(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : g d~ -> xI+ xI- d~  
C  
C Crossing   1 is g d~ -> xI+ xI- d~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps.inc"
      Include "nexternal.inc"
      Include "maxamps_xixj.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER(             NCOMB=  32, NCROSS=  1)
      INTEGER    THEL
      PARAMETER(THEL=NCOMB*NCROSS)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 MATRIX_GDX_XIXIDX
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I
      LOGICAL GOODHEL(NCOMB,NCROSS)
      INTEGER NGRAPHS
      REAL*8 hwgt, xtot, xtry, xrej, xr, yfrac(0:ncomb)
      INTEGER idum, ngood, igood(ncomb), jhel, j, jj
      LOGICAL warned
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2(maxamps), jamp2(0:maxflow)
      common/to_amps_gdx_xixi/  amp2,       jamp2

      character*79         hel_buff
      !common/to_helicity/  hel_buff

      REAL*8 POL(2)
      !common/to_polarization/ POL

      integer          isum_hel
      logical                    multi_channel
      !common/to_matrix/isum_hel, multi_channel
      INTEGER MAPCONFIG(0:LMAXCONFIGS), ICONFIG
      !common/to_mconfigs/mapconfig, iconfig
      DATA NTRY,IDUM /0,-1/
      DATA xtry, xrej, ngood /0,0,0/
      DATA warned, isum_hel/.false.,0/
      DATA multi_channel/.true./
      SAVE yfrac, igood, jhel
      DATA NGRAPHS /    7/          
      DATA jamp2(0) /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA(NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA(NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA(NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA(NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA(NHEL(IHEL,   5),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA(NHEL(IHEL,   6),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA(NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA(NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA(NHEL(IHEL,   9),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA(NHEL(IHEL,  10),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA(NHEL(IHEL,  11),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA(NHEL(IHEL,  12),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA(NHEL(IHEL,  13),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA(NHEL(IHEL,  14),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA(NHEL(IHEL,  15),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA(NHEL(IHEL,  16),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA(NHEL(IHEL,  17),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA(NHEL(IHEL,  18),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA(NHEL(IHEL,  19),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA(NHEL(IHEL,  20),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA(NHEL(IHEL,  21),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA(NHEL(IHEL,  22),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA(NHEL(IHEL,  23),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA(NHEL(IHEL,  24),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA(NHEL(IHEL,  25),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA(NHEL(IHEL,  26),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA(NHEL(IHEL,  27),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA(NHEL(IHEL,  28),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA(NHEL(IHEL,  29),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA(NHEL(IHEL,  30),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA(NHEL(IHEL,  31),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA(NHEL(IHEL,  32),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA(  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA(IDEN(IHEL),IHEL=  1,  1) /  96/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      CALL SWITCHMOM(P1,P,IC(1,IPROC),JC,NEXTERNAL)
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO

          DO IHEL=1,NGRAPHS
              amp2(ihel)=0d0
          ENDDO
          DO IHEL=1,int(jamp2(0))
              jamp2(ihel)=0d0
          ENDDO
      ANS(IPROC) = 0D0
          DO IHEL=1,NCOMB
             IF(GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=MATRIX_GDX_XIXIDX(P,NHEL(1,IHEL),JC(1))            
               ANS(IPROC)=ANS(IPROC)+T
               IF(T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION MATRIX_GDX_XIXIDX(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : g d~ -> xI+ xI- d~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER(NGRAPHS=   7,NEIGEN=  1) 
      include "genps.inc"
      include "nexternal.inc"
      include "maxamps_xixj.inc"
      INTEGER    NWAVEFUNCS, NCOLOR
      PARAMETER(NWAVEFUNCS=  13, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER(ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL)
C  
C LOCAL VARIABLES 
C  
      INTEGER I,J
      COMPLEX*16 ZTEMP
      REAL*8 DENOM(NCOLOR), CF(NCOLOR,NCOLOR)
      COMPLEX*16 AMP(NGRAPHS), JAMP(NCOLOR)
      COMPLEX*16 W(18,NWAVEFUNCS)
C  
C GLOBAL VARIABLES
C  
      Double Precision amp2(maxamps), jamp2(0:maxflow)
      common/to_amps_gdx_xixi/  amp2,       jamp2
      include "coupl.inc"
C
C VARIABLES TO APPLY THE DIAGRAM SUBTRACTION SCHEMES 
C
#include "dsubtraction_xixi.h"
C  
C COLOR DATA
C  
      DATA Denom(1)/            1/                                       
      DATA(CF(i,1),i=1,1) /     4/                                  
C               T[ 2, 5, 1]                                                
C ----------
C BEGIN CODE
C ----------
      CALL VXXXXX(P(0,1),ZERO,NHEL(1),-1*IC(1),W(1,1))        
      CALL OXXXXX(P(0,2),ZERO,NHEL(2),-1*IC(2),W(1,2))        
      CALL IXXXXX(P(0,3),MXI,NHEL(3),-1*IC(3),W(1,3))         
      CALL OXXXXX(P(0,4),MXI,NHEL(4),+1*IC(4),W(1,4))         
      CALL IXXXXX(P(0,5),ZERO,NHEL(5),-1*IC(5),W(1,5))    
      
      ! diagram 1
      CALL HIOXXX(W(1,3),W(1,2),GULXIP,MUL,WUL,W(1,6))                                                          
      CALL FSOXXX(W(1,4),W(1,6),GULXIM,ZERO,ZERO,W(1,7))                                                          
      CALL IOVXXX(W(1,5),W(1,7),W(1,1),GG,AMP(1))        
      ! diagram 2 <- resonant
      CALL HVSXXX(W(1,1),W(1,6),GC,MUL,WULR,W(1,8))     
      CALL IOSXXX(W(1,5),W(1,4),W(1,8),GULXIM,AMP(2))     
      ! diagram 3 <- resonant
      CALL FVOXXX(W(1,2),W(1,1),GG,ZERO,ZERO,W(1,9))     
      CALL HIOXXX(W(1,3),W(1,9),GULXIP,MUL,WULR,W(1,10))                                                          
      CALL IOSXXX(W(1,5),W(1,4),W(1,10),GULXIM,AMP(3))      
      ! diagram 4
      CALL JIOXXX(W(1,3),W(1,4),GAX,ZERO,AWIDTH,W(1,11))    
      CALL FVIXXX(W(1,5),W(1,1),GG,ZERO,ZERO,W(1,12))     
      CALL IOVXXX(W(1,12),W(1,2),W(1,11),GAD,AMP(4))     
      ! diagram 5
      CALL JIOXXX(W(1,3),W(1,4),GZXII,ZMASS,ZWIDTH,W(1,13))                                                          
      CALL IOVXXX(W(1,12),W(1,2),W(1,13),GZD,AMP(5))      
      ! diagram 6
      CALL IOVXXX(W(1,5),W(1,9),W(1,11),GAD,AMP(6))    
      ! diagram 7
      CALL IOVXXX(W(1,5),W(1,9),W(1,13),GZD,AMP(7))      
      
c delete resonant diagrams (diagram removal type I)
#ifdef DR_I
      AMP(2)  = 0 ! diagram 2,  s45, L, XI-
      AMP(3)  = 0 ! diagram 3,  s45, L, XI-
#endif
      
      JAMP(   1) = +AMP(   1)+AMP(   2)+AMP(   3)-AMP(   4)-AMP(   5)
     &             -AMP(   6)-AMP(   7)
      MATRIX_GDX_XIXIDX = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP =(0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          MATRIX_GDX_XIXIDX =MATRIX_GDX_XIXIDX
     &                       +ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      
c delete the resonant squared amplitude but keep the interference term
c (diagram removal type II)
#ifdef DR_II
      JAMPR(1)  = +AMP(   2)+AMP(   3)
     
      MATRIX_RESONANT = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMPR(J)
          ENDDO
          MATRIX_RESONANT =MATRIX_RESONANT
     &                      +ZTEMP*DCONJG(JAMPR(I))/DENOM(I)   
      ENDDO

      MATRIX_GDX_XIXIDX = MATRIX_GDX_XIXIDX - MATRIX_RESONANT
#endif

c delete the on-shell contributions of the resonant diagrams but keep
c the interference term and the off-shell contributions
#ifdef DSUB_I
      S   = momsum2sq(p(0:3,1), p(0:3,2))
      S35 = momsum2sq(p(0:3,3), p(0:3,5))
      S45 = momsum2sq(p(0:3,4), p(0:3,5))
      
      RATIO45L   = 0D0
      COUNTER45L = 0D0
      
      ! apply on shell conditions (don't subtract anything if the squark
      ! is off-shell, but if the quark is on-shell subtract the divergent
      ! part) -> Prospino scheme 1211.0286
      ! S45 = MUL^2, m4 = MXI, m5 = 0D0
      if( (S.ge.(MUL+MXI)**2) .and. (MUL.ge.MXI)) then
        call off_to_on(p,"ul45",p_OS)
        RATIO45L = (MUL*WREG)**2/((S45-MUL**2)**2+(MUL*WREG)**2)
        COUNTER45L = RATIO45L*MATRIX_GDX_XIXIDX_RES(p_OS,NHEL,IC,"ul45")
      endif
      
      MATRIX_GDX_XIXIDX = MATRIX_GDX_XIXIDX - COUNTER45L
#endif
      
      ! amp2 and jamp2 are not used
      Do I = 1, NGRAPHS
          amp2(i)=amp2(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2(i)=Jamp2(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       

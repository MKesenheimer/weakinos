      REAL*8 FUNCTION MATRIX_UXG_XIXJUX_RES(P,NHEL,IC,CHAN)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : u~ g -> xi+ xj- u~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER(NGRAPHS=   5,NEIGEN=  1) 
      include "genps.inc"
      include "nexternal.inc"
      include "maxamps_xixj.inc"
      INTEGER    NWAVEFUNCS, NCOLOR
      PARAMETER(NWAVEFUNCS=  14, NCOLOR=   1) 
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
      include "coupl.inc"
C
C VARIABLES TO APPLY THE DIAGRAM SUBTRACTION SCHEMES 
C
#define _EXTFUNCT_
#include "dsubtraction_xixj.h"
C  
C COLOR DATA
C  
      DATA Denom(1)/            1/                                       
      DATA(CF(i,1),i=1,1) /     4/                                  
C               T[ 1, 5, 2]                                                
C ----------
C BEGIN CODE
C ----------
      CALL OXXXXX(P(0,1),ZERO,NHEL(1),-1*IC(1),W(1,1))        
      CALL VXXXXX(P(0,2),ZERO,NHEL(2),-1*IC(2),W(1,2))        
      CALL IXXXXX(P(0,3),MXI,NHEL(3),-1*IC(3),W(1,3))         
      CALL OXXXXX(P(0,4),MXJ,NHEL(4),+1*IC(4),W(1,4))         
      CALL IXXXXX(P(0,5),ZERO,NHEL(5),-1*IC(5),W(1,5))       
      
      ! diagram 1
      CALL JIOXXX(W(1,3),W(1,4),GZXIJ,ZMASS,ZWIDTH,W(1,6))                                                          
      CALL FVIXXX(W(1,5),W(1,2),GG,ZERO,ZERO,W(1,7))     
      CALL IOVXXX(W(1,7),W(1,1),W(1,6),GZU,AMP(1))           
      ! diagram 2
      CALL FVOXXX(W(1,1),W(1,2),GG,ZERO,ZERO,W(1,8))     
      CALL IOVXXX(W(1,5),W(1,8),W(1,6),GZU,AMP(2))       
      ! diagram 3 
      CALL IXXXXX(P(0,4),MXJ,NHEL(4),-1*IC(4),W(1,9))         
      CALL OXXXXX(P(0,5),ZERO,NHEL(5),+1*IC(5),W(1,10))        
      CALL HIOXXX(W(1,9),W(1,1),GDLXJP,MDL,WDL,W(1,11))                                                          
      CALL FVOCXX(W(1,10),W(1,2),GG,ZERO,ZERO,W(1,12))    
      CALL IOSCXX(W(1,3),W(1,12),W(1,11),GDLXIM,AMP(3))        
      ! diagram 4 <- resonant
      CALL HIOCXX(W(1,3),W(1,10),GDLXIM,MDL,WDLR,W(1,13))                                                          
      CALL VSSXXX(W(1,2),W(1,11),W(1,13),GC,AMP(4))             
      ! diagram 5 <- resonant
      CALL HIOXXX(W(1,9),W(1,8),GDLXJP,MDL,WDLR,W(1,14))                                                          
      CALL IOSCXX(W(1,3),W(1,10),W(1,14),GDLXIM,AMP(5))      
      
      ! S35 = MDL^2
      if (CHAN.eq.0) then
        JAMP(1) = +AMP(4)+AMP(5)
      else
        print*,"error in r_gu_xixju_res.f"
        stop
      endif
      
      MATRIX_UXG_XIXJUX_RES = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP =(0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          MATRIX_UXG_XIXJUX_RES =MATRIX_UXG_XIXJUX_RES
     &                       +ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO

      END
       
       

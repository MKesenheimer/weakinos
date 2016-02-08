      SUBROUTINE SMATRIX_DXG_NINJDX(P1,ANS)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : d~ g -> ni nj d~  
C  
C Crossing   1 is d~ g -> ni nj d~  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "genps.inc"
      Include "nexternal.inc"
      Include "maxamps_ninj.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER(             NCOMB=  32, NCROSS=  1)
      INTEGER    THEL
      PARAMETER(THEL=NCOMB*NCROSS)
C
C FINAL STATES
#include "finalstate.h"
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T, P(0:3,NEXTERNAL)
      REAL*8 MATRIX_DXG_NINJDX
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
      common/to_amps_dxg/  amp2,       jamp2

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
      DATA NGRAPHS /   14/          
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
                 T=MATRIX_DXG_NINJDX(P,NHEL(1,IHEL),JC(1))            
               ANS(IPROC)=ANS(IPROC)+T
               IF(T .NE. 0D0 .AND. .NOT.    GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
      if(final1.eq.final2) then ! equal final states
        ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC)*2)
      else
        ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      endif  
      ENDDO
      END
       
       
      REAL*8 FUNCTION MATRIX_DXG_NINJDX(P,NHEL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : d~ g -> ni nj d~  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER(NGRAPHS=  14,NEIGEN=  1) 
      include "genps.inc"
      include "nexternal.inc"
      include "maxamps_ninj.inc"
      INTEGER    NWAVEFUNCS, NCOLOR
      PARAMETER(NWAVEFUNCS=  24, NCOLOR=   1) 
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
      common/to_amps_dxg/  amp2,       jamp2
      include "coupl.inc"
C
C VARIABLES TO APPLY THE DIAGRAM SUBTRACTION SCHEMES 
C
#include "dsubtraction_ninj.h"
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
      CALL OXXXXX(P(0,3),MNI,NHEL(3),+1*IC(3),W(1,3))         
      CALL IXXXXX(P(0,4),MNJ,NHEL(4),-1*IC(4),W(1,4))         
      CALL IXXXXX(P(0,5),ZERO,NHEL(5),-1*IC(5),W(1,5))   
      
      ! diagram 1
      CALL HIOXXX(W(1,4),W(1,1),GDLNJP,MDL,WDL,W(1, 6))                                                          
      CALL FVIXXX(W(1,5),W(1,2),GG,ZERO,ZERO,W(1,7))     
      CALL IOSXXX(W(1,7),W(1,3),W(1,6),GDLNIM,AMP(1))   
      ! diagram 2
      CALL HIOXXX(W(1,4),W(1,1),GDRNJP,MDR,WDR,W(1, 8))                                                          
      CALL IOSXXX(W(1,7),W(1,3),W(1,8),GDRNIM,AMP(2))   
      
      ! diagram 3 <- resonant
      CALL HIOXXX(W(1,5),W(1,3),GDLNIM,MDL,WDLR,W(1, 9))                                                          
      CALL VSSXXX(W(1,2),W(1,6),W(1,9),GC,AMP(3))     
      ! diagram 4 <- resonant
      CALL HIOXXX(W(1,5),W(1,3),GDRNIM,MDR,WDRR,W(1, 10))                                                          
      CALL VSSXXX(W(1,2),W(1,8),W(1,10),GC,AMP(4))     
      
      ! diagram 5 <- resonant
      CALL FVOXXX(W(1,1),W(1,2),GG,ZERO,ZERO,W(1,11))     
      CALL HIOXXX(W(1,4),W(1,11),GDLNJP,MDL,WDLR,W(1, 12))                                                          
      CALL IOSXXX(W(1,5),W(1,3),W(1,12),GDLNIM,AMP(5))       
      ! diagram 6 <- resonant
      CALL HIOXXX(W(1,4),W(1,11),GDRNJP,MDR,WDRR,W(1, 13))                                                          
      CALL IOSXXX(W(1,5),W(1,3),W(1,13),GDRNIM,AMP(6))    
      
      ! diagram 7
      CALL JIOCXX(W(1,4),W(1,3),GZNIJ,ZMASS,ZWIDTH,W(1, 14))                                                          
      CALL IOVXXX(W(1,7),W(1,1),W(1,14),GZD,AMP(7))       
      ! diagram 8
      CALL IOVXXX(W(1,5),W(1,11),W(1,14),GZD,AMP(8))        
      
      ! diagram 9
      CALL IXXXXX(P(0,3),MNI,NHEL(3),-1*IC(3),W(1,15))         
      CALL OXXXXX(P(0,4),MNJ,NHEL(4),+1*IC(4),W(1,16))         
      CALL HIOXXX(W(1,15),W(1,1),GDLNIP,MDL,WDL,W(1, 17))                                                          
      CALL FSOXXX(W(1,16),W(1,17),GDLNJM,ZERO,ZERO,W(1, 18))                                                          
      CALL IOVXXX(W(1,5),W(1,18),W(1,2),GG,AMP(9))          
      ! diagram 10
      CALL HIOXXX(W(1,15),W(1,1),GDRNIP,MDR,WDR,W(1, 19))                                                          
      CALL FSOXXX(W(1,16),W(1,19),GDRNJM,ZERO,ZERO,W(1, 20))                                                          
      CALL IOVXXX(W(1,5),W(1,20),W(1,2),GG,AMP(10))          
      
      ! diagram 11 <- resonant
      CALL HVSXXX(W(1,2),W(1,17),GC,MDL,WDLR,W(1,21))     
      CALL IOSXXX(W(1,5),W(1,16),W(1,21),GDLNJM,AMP(11))         
      ! diagram 12 <- resonant
      CALL HVSXXX(W(1,2),W(1,19),GC,MDR,WDRR,W(1,22))     
      CALL IOSXXX(W(1,5),W(1,16),W(1,22),GDRNJM,AMP(12))       
      
      ! diagram 13 <- resonant
      CALL HIOXXX(W(1,15),W(1,11),GDLNIP,MDL,WDLR,W(1, 23))                                                          
      CALL IOSXXX(W(1,5),W(1,16),W(1,23),GDLNJM,AMP(13))      
      ! diagram 14 <- resonant
      CALL HIOXXX(W(1,15),W(1,11),GDRNIP,MDR,WDRR,W(1, 24))                                                          
      CALL IOSXXX(W(1,5),W(1,16),W(1,24),GDRNJM,AMP(14))    
      
c delete resonant diagrams (diagram removal type I)
#ifdef DR_I
      AMP(3) = 0
      AMP(4) = 0
      AMP(5) = 0
      AMP(6) = 0
      AMP(11) = 0
      AMP(12) = 0
      AMP(13) = 0
      AMP(14) = 0
#endif
      
      JAMP(   1) = -AMP(   1)-AMP(   2)-AMP(   3)-AMP(   4)-AMP(   5)
     &             -AMP(   6)+AMP(   7)+AMP(   8)+AMP(   9)+AMP(  10)
     &             +AMP(  11)+AMP(  12)+AMP(  13)+AMP(  14)
      MATRIX_DXG_NINJDX = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP =(0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          MATRIX_DXG_NINJDX =MATRIX_DXG_NINJDX
     &                       +ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      
c delete the resonant squared amplitude but keep the interference term
c (diagram removal type II)
#ifdef DR_II
      JAMPR(1)  = -AMP(3)-AMP(4)-AMP(5)-AMP(6)+AMP(11)+AMP(12)
     &            +AMP(13)+AMP(14)
     
      MATRIX_RESONANT = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMPR(J)
          ENDDO
          MATRIX_RESONANT =MATRIX_RESONANT
     &                      +ZTEMP*DCONJG(JAMPR(I))/DENOM(I)   
      ENDDO

      MATRIX_DXG_NINJDX = MATRIX_DXG_NINJDX - MATRIX_RESONANT
#endif

c if DSUB_II is used, the resonant matrix element |M_R|**2 will be added
c back in subroutine setosresreal
#if defined(DSUB_II) || defined(DSUB_II_TEST)
      MATRIX_DXG_NINJDX = MATRIX_DXG_NINJDX
     &                - MATRIX_DXG_NINJDX_RES(p,NHEL,IC,"dl35")
     &                - MATRIX_DXG_NINJDX_RES(p,NHEL,IC,"dr35")
     &                - MATRIX_DXG_NINJDX_RES(p,NHEL,IC,"dl45")
     &                - MATRIX_DXG_NINJDX_RES(p,NHEL,IC,"dr45")
#endif

c delete the on-shell contributions of the resonant diagrams but keep
c the interference term and the off-shell contributions
#ifdef DSUB_I
      S   = momsum2sq(p(0:3,1), p(0:3,2))
      S35 = momsum2sq(p(0:3,3), p(0:3,5))
      S45 = momsum2sq(p(0:3,4), p(0:3,5))
      
      RATIO35L   = 0D0
      RATIO35R   = 0D0
      RATIO45L   = 0D0
      RATIO45R   = 0D0
      COUNTER35L = 0D0
      COUNTER35R = 0D0
      COUNTER45L = 0D0
      COUNTER45R = 0D0
      
      ! apply on shell conditions (don't subtract anything if the squark
      ! is off-shell, but if the quark is on-shell subtract the divergent
      ! part) -> Prospino scheme 1211.0286
      ! S35 = MDL^2, m3 = MNI, m5 = 0D0
      if( (S.ge.(MDL+dabs(MNJ))**2) .and. (MDL.ge.dabs(MNI))) then
        !off_to_on(p,i,j,k,mij,mi,mj,mk,p_OS)
        call off_to_on(p,"dl35",p_OS)              ! off_to_on the momenta p to on-shell momenta p_OS
        RATIO35L = (MDL*WREG)**2/((S35-MDL**2)**2+(MDL*WREG)**2)  ! calculate the ratio of the breit wigner functions
        COUNTER35L = RATIO35L*MATRIX_DXG_NINJDX_RES(p_OS,NHEL,IC,"dl35") ! generate the counter term
      endif
      ! S35 = MDR^2, m3 = MNI, m5 = 0D0
      if( (S.ge.(MDR+dabs(MNJ))**2) .and. (MDR.ge.dabs(MNI))) then
        call off_to_on(p,"dr35",p_OS)
        RATIO35R = (MDR*WREG)**2/((S35-MDR**2)**2+(MDR*WREG)**2)
        COUNTER35R = RATIO35R*MATRIX_DXG_NINJDX_RES(p_OS,NHEL,IC,"dr35")
      endif
      ! S45 = MDL^2, m4 = MNJ, m5 = 0D0
      if( (S.ge.(MDL+dabs(MNI))**2) .and. (MDL.ge.dabs(MNJ))) then
        call off_to_on(p,"dl45",p_OS)
        RATIO45L = (MDL*WREG)**2/((S45-MDL**2)**2+(MDL*WREG)**2)
        COUNTER45L = RATIO45L*MATRIX_DXG_NINJDX_RES(p_OS,NHEL,IC,"dl45")
      endif
      ! S45 = MDR^2, m4 = MNJ, m5 = 0D0
      if( (S.ge.(MDR+dabs(MNI))**2) .and. (MDR.ge.dabs(MNJ))) then
        call off_to_on(p,"dr45",p_OS)
        RATIO45R = (MDR*WREG)**2/((S45-MDR**2)**2+(MDR*WREG)**2)
        COUNTER45R = RATIO45R*MATRIX_DXG_NINJDX_RES(p_OS,NHEL,IC,"dr45")
      endif
      
      MATRIX_DXG_NINJDX = MATRIX_DXG_NINJDX - COUNTER35L - COUNTER35R
     &                                      - COUNTER45L - COUNTER45R
#endif
      
      ! amp2 and jamp2 isn't used
      Do I = 1, NGRAPHS
          amp2(i)=amp2(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2(i)=Jamp2(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
C      CALL GAUGECHECK(JAMP,ZTEMP,EIGEN_VEC,EIGEN_VAL,NCOLOR,NEIGEN) 
      END
       
       

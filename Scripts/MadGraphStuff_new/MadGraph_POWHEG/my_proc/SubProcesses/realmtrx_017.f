      SUBROUTINE SREALMTRX_017(P1,ANS)
C  
C Generated by MadGraph II                                              
C MadGraph StandAlone Version
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C AND HELICITIES
C FOR THE POINT IN PHASE SPACE P(0:3,NEXTERNAL)
C  
C FOR PROCESS : b~ b -> n1 n1 g  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      Include "nexternal.inc"
      INTEGER                 NCOMB,     NCROSS         
      PARAMETER (             NCOMB=  32, NCROSS=  1)
      INTEGER    THEL
      PARAMETER (THEL=NCOMB*NCROSS)
      INTEGER NGRAPHS
      PARAMETER (NGRAPHS=  20)
C  
C ARGUMENTS 
C  
      REAL*8 P1(0:3,NEXTERNAL),ANS(NCROSS)
C  
C LOCAL VARIABLES 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL,NCOMB),NTRY
      REAL*8 T
      REAL*8 REALMTRX_017
      REAL*8 ZERO
      PARAMETER(ZERO=0d0)
      INTEGER IHEL,IDEN(NCROSS),IC(NEXTERNAL,NCROSS)
      INTEGER IPROC,JC(NEXTERNAL), I,L,K
      LOGICAL GOODHEL(NCOMB,NCROSS)
      DATA NTRY/0/
      INTEGER NGOOD,igood(ncomb),jhel
      data ngood /0/
      save igood,jhel
      REAL*8 hwgt
      integer maxamps
      parameter (maxamps=6000)
      Double Precision amp2(maxamps), jamp2(0:maxamps)
      common/to_Ramps_017/  amp2,       jamp2

      integer j,jj
      integer max_bhel
      parameter ( max_bhel =          32 )

      INTEGER NCOLOR
      DATA NCOLOR   /   1/          
      DATA GOODHEL/THEL*.FALSE./
      DATA (NHEL(IHEL,   1),IHEL=1, 5) /-1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,   2),IHEL=1, 5) /-1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,   3),IHEL=1, 5) /-1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,   4),IHEL=1, 5) /-1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,   5),IHEL=1, 5) /-1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,   6),IHEL=1, 5) /-1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,   7),IHEL=1, 5) /-1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,   8),IHEL=1, 5) /-1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,   9),IHEL=1, 5) /-1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  10),IHEL=1, 5) /-1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  11),IHEL=1, 5) /-1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  12),IHEL=1, 5) /-1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  13),IHEL=1, 5) /-1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  14),IHEL=1, 5) /-1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  15),IHEL=1, 5) /-1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  16),IHEL=1, 5) /-1, 1, 1, 1, 1/
      DATA (NHEL(IHEL,  17),IHEL=1, 5) / 1,-1,-1,-1,-1/
      DATA (NHEL(IHEL,  18),IHEL=1, 5) / 1,-1,-1,-1, 1/
      DATA (NHEL(IHEL,  19),IHEL=1, 5) / 1,-1,-1, 1,-1/
      DATA (NHEL(IHEL,  20),IHEL=1, 5) / 1,-1,-1, 1, 1/
      DATA (NHEL(IHEL,  21),IHEL=1, 5) / 1,-1, 1,-1,-1/
      DATA (NHEL(IHEL,  22),IHEL=1, 5) / 1,-1, 1,-1, 1/
      DATA (NHEL(IHEL,  23),IHEL=1, 5) / 1,-1, 1, 1,-1/
      DATA (NHEL(IHEL,  24),IHEL=1, 5) / 1,-1, 1, 1, 1/
      DATA (NHEL(IHEL,  25),IHEL=1, 5) / 1, 1,-1,-1,-1/
      DATA (NHEL(IHEL,  26),IHEL=1, 5) / 1, 1,-1,-1, 1/
      DATA (NHEL(IHEL,  27),IHEL=1, 5) / 1, 1,-1, 1,-1/
      DATA (NHEL(IHEL,  28),IHEL=1, 5) / 1, 1,-1, 1, 1/
      DATA (NHEL(IHEL,  29),IHEL=1, 5) / 1, 1, 1,-1,-1/
      DATA (NHEL(IHEL,  30),IHEL=1, 5) / 1, 1, 1,-1, 1/
      DATA (NHEL(IHEL,  31),IHEL=1, 5) / 1, 1, 1, 1,-1/
      DATA (NHEL(IHEL,  32),IHEL=1, 5) / 1, 1, 1, 1, 1/
      DATA (  IC(IHEL,  1),IHEL=1, 5) / 1, 2, 3, 4, 5/
      DATA (IDEN(IHEL),IHEL=  1,  1) /  72/
C ----------
C BEGIN CODE
C ----------
      NTRY=NTRY+1
      DO IPROC=1,NCROSS
      DO IHEL=1,NEXTERNAL
         JC(IHEL) = +1
      ENDDO
      DO IHEL=1,NGRAPHS
          amp2(ihel)=0d0
      ENDDO
      jamp2(0)=dble(NCOLOR)
      DO IHEL=1,int(jamp2(0))
          jamp2(ihel)=0d0
      ENDDO
      ANS(IPROC) = 0D0
          DO IHEL=1,NCOMB
             IF (GOODHEL(IHEL,IPROC) .OR. NTRY .LT. 2) THEN
                 T=REALMTRX_017(P1,NHEL(1,IHEL),IHEL,JC(1))              
               ANS(IPROC)=ANS(IPROC)+T
               IF (T .GT. 0D0 .AND. .NOT. GOODHEL(IHEL,IPROC)) THEN
                   GOODHEL(IHEL,IPROC)=.TRUE.
                   NGOOD = NGOOD +1
                   IGOOD(NGOOD) = IHEL
               ENDIF
             ENDIF
          ENDDO
      ANS(IPROC)=ANS(IPROC)/DBLE(IDEN(IPROC))
      ENDDO
      END
       
       
      REAL*8 FUNCTION REALMTRX_017(P,NHEL,HELL,IC)
C  
C Generated by MadGraph II                                              
C RETURNS AMPLITUDE SQUARED SUMMED/AVG OVER COLORS
C FOR THE POINT WITH EXTERNAL LINES W(0:6,NEXTERNAL)
C  
C FOR PROCESS : b~ b -> n1 n1 g  
C  
      IMPLICIT NONE
C  
C CONSTANTS
C  
      INTEGER    NGRAPHS,    NEIGEN 
      PARAMETER (NGRAPHS=  20,NEIGEN=  1) 
      include "nexternal.inc"
      INTEGER    NWAVEFUNCS     , NCOLOR
      PARAMETER (NWAVEFUNCS=  25, NCOLOR=   1) 
      REAL*8     ZERO
      PARAMETER (ZERO=0D0)
C  
C ARGUMENTS 
C  
      REAL*8 P(0:3,NEXTERNAL)
      INTEGER NHEL(NEXTERNAL), IC(NEXTERNAL), HELL
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
      integer maxamps
      parameter (maxamps=6000)
      Double Precision amp2(maxamps), jamp2(0:maxamps)
      common/to_Ramps_017/  amp2,       jamp2
      integer max_bhel
      parameter ( max_bhel =          32 )
      include "coupl.inc"
C  
C COLOR DATA
C  
      DATA Denom(1  )/            1/                                       
      DATA (CF(i,1  ),i=1  ,1  ) /     4/                                  
C               T[ 1, 2, 5]                                                
C ----------
C BEGIN CODE
C ----------
      CALL OXXXXX(P(0,1   ),BMASS ,NHEL(1   ),-1*IC(1   ),W(1,1   ))       
      CALL IXXXXX(P(0,2   ),BMASS ,NHEL(2   ),+1*IC(2   ),W(1,2   ))       
      CALL IXXXXX(P(0,3   ),MN1 ,NHEL(3   ),-1*IC(3   ),W(1,3   ))         
      CALL OXXXXX(P(0,4   ),MN1 ,NHEL(4   ),+1*IC(4   ),W(1,4   ))         
      CALL VXXXXX(P(0,5   ),ZERO ,NHEL(5   ),+1*IC(5   ),W(1,5   ))        
      CALL HIOXXX(W(1,2   ),W(1,4   ),GB1N1M ,MB1     ,WB1     ,W(1,       
     &     6   ))                                                          
      CALL FVOXXX(W(1,1   ),W(1,5   ),GG ,BMASS   ,ZERO    ,W(1,7   ))     
      CALL IOSXXX(W(1,3   ),W(1,7   ),W(1,6   ),GB1N1P ,AMP(1   ))         
      CALL HIOXXX(W(1,2   ),W(1,4   ),GB2N1M ,MB2     ,WB2     ,W(1,       
     &     8   ))                                                          
      CALL IOSXXX(W(1,3   ),W(1,7   ),W(1,8   ),GB2N1P ,AMP(2   ))         
      CALL HIOXXX(W(1,3   ),W(1,1   ),GB1N1P ,MB1     ,WB1     ,W(1,       
     &     9   ))                                                          
      CALL FSOXXX(W(1,4   ),W(1,9   ),GB1N1M ,BMASS   ,ZERO    ,W(1,       
     &     10  ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,10  ),W(1,5   ),GG ,AMP(3   ))             
      CALL HIOXXX(W(1,3   ),W(1,1   ),GB2N1P ,MB2     ,WB2     ,W(1,       
     &     11  ))                                                          
      CALL FSOXXX(W(1,4   ),W(1,11  ),GB2N1M ,BMASS   ,ZERO    ,W(1,       
     &     12  ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,12  ),W(1,5   ),GG ,AMP(4   ))             
      CALL VSSXXX(W(1,5   ),W(1,9   ),W(1,6   ),GC ,AMP(5   ))             
      CALL VSSXXX(W(1,5   ),W(1,11  ),W(1,8   ),GC ,AMP(6   ))             
      CALL JIOXXX(W(1,3   ),W(1,4   ),GZN11 ,ZMASS   ,ZWIDTH  ,W(1,        
     &     13  ))                                                          
      CALL IOVXXX(W(1,2   ),W(1,7   ),W(1,13  ),GZD ,AMP(7   ))            
      CALL HIOXXX(W(1,3   ),W(1,4   ),GH1N11 ,MH1     ,WH1     ,W(1,       
     &     14  ))                                                          
      CALL IOSXXX(W(1,2   ),W(1,7   ),W(1,14  ),GH1BB ,AMP(8   ))          
      CALL HIOXXX(W(1,3   ),W(1,4   ),GH2N11 ,MH2     ,WH2     ,W(1,       
     &     15  ))                                                          
      CALL IOSXXX(W(1,2   ),W(1,7   ),W(1,15  ),GH2BB ,AMP(9   ))          
      CALL HIOXXX(W(1,3   ),W(1,4   ),GH3N11 ,MH3     ,WH3     ,W(1,       
     &     16  ))                                                          
      CALL IOSXXX(W(1,2   ),W(1,7   ),W(1,16  ),GH3BB ,AMP(10  ))          
      CALL FVIXXX(W(1,2   ),W(1,5   ),GG ,BMASS   ,ZERO    ,W(1,17  ))     
      CALL IOVXXX(W(1,17  ),W(1,1   ),W(1,13  ),GZD ,AMP(11  ))            
      CALL IOSXXX(W(1,17  ),W(1,1   ),W(1,14  ),GH1BB ,AMP(12  ))          
      CALL IOSXXX(W(1,17  ),W(1,1   ),W(1,15  ),GH2BB ,AMP(13  ))          
      CALL IOSXXX(W(1,17  ),W(1,1   ),W(1,16  ),GH3BB ,AMP(14  ))          
      CALL IXXXXX(P(0,1   ),BMASS ,NHEL(1   ),+1*IC(1   ),W(1,18  ))       
      CALL OXXXXX(P(0,3   ),MN1 ,NHEL(3   ),+1*IC(3   ),W(1,19  ))         
      CALL HIOXXX(W(1,2   ),W(1,19  ),GB1N1M ,MB1     ,WB1     ,W(1,       
     &     20  ))                                                          
      CALL FSOCXX(W(1,4   ),W(1,20  ),GB1N1P  ,BMASS   ,ZERO    ,W(1,      
     &     21  ))                                                          
      CALL IOVCXX(W(1,18  ),W(1,21  ),W(1,5   ),GG  ,AMP(15  ))            
      CALL HIOXXX(W(1,2   ),W(1,19  ),GB2N1M ,MB2     ,WB2     ,W(1,       
     &     22  ))                                                          
      CALL FSOCXX(W(1,4   ),W(1,22  ),GB2N1P  ,BMASS   ,ZERO    ,W(1,      
     &     23  ))                                                          
      CALL IOVCXX(W(1,18  ),W(1,23  ),W(1,5   ),GG  ,AMP(16  ))            
      CALL HIOCXX(W(1,18  ),W(1,4   ),GB1N1P  ,MB1     ,WB1     ,W(1,      
     &     24  ))                                                          
      CALL IOSXXX(W(1,17  ),W(1,19  ),W(1,24  ),GB1N1M ,AMP(17  ))         
      CALL HIOCXX(W(1,18  ),W(1,4   ),GB2N1P  ,MB2     ,WB2     ,W(1,      
     &     25  ))                                                          
      CALL IOSXXX(W(1,17  ),W(1,19  ),W(1,25  ),GB2N1M ,AMP(18  ))         
      CALL VSSXXX(W(1,5   ),W(1,24  ),W(1,20  ),GC ,AMP(19  ))             
      CALL VSSXXX(W(1,5   ),W(1,25  ),W(1,22  ),GC ,AMP(20  ))             
      JAMP(   1) = +AMP(   1)+AMP(   2)+AMP(   3)+AMP(   4)+AMP(   5)
     &             +AMP(   6)-AMP(   7)-AMP(   8)-AMP(   9)-AMP(  10)
     &             -AMP(  11)-AMP(  12)-AMP(  13)-AMP(  14)+AMP(  15)
     &             +AMP(  16)+AMP(  17)+AMP(  18)+AMP(  19)+AMP(  20)
      REALMTRX_017 = 0.D0 
      DO I = 1, NCOLOR
          ZTEMP = (0.D0,0.D0)
          DO J = 1, NCOLOR
              ZTEMP = ZTEMP + CF(J,I)*JAMP(J)
          ENDDO
          REALMTRX_017 =REALMTRX_017+ZTEMP*DCONJG(JAMP(I))/DENOM(I)   
      ENDDO
      Do I = 1, NGRAPHS
          amp2(i)=amp2(i)+amp(i)*dconjg(amp(i))
      Enddo
      Do I = 1, NCOLOR
          Jamp2(i)=Jamp2(i)+Jamp(i)*dconjg(Jamp(i))
      Enddo
      END
       
       

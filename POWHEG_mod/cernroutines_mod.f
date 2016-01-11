c MK: copied and modified version of cernroutines.f, revision 3154
c changes marked with "! MK:"

c# 1 "dilog64.F"
c# 1 "<built-in>"
c# 1 "<command line>"
c# 1 "dilog64.F"
*
* $Id: dilog64.F,v 1.1.1.1 1996/04/01 15:02:05 mclareni Exp $
*
* $Log: dilog64.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:05  mclareni
* Mathlib gen
*
*
c# 1 "/usr/local/home/video/cernlib/2005/include/gen/pilot.h" 1
























c# 40 "/usr/local/home/video/cernlib/2005/include/gen/pilot.h"

c# 57 "/usr/local/home/video/cernlib/2005/include/gen/pilot.h"



























































c# 10 "dilog64.F" 2

      FUNCTION DDILOG(X)
c# 1 "/usr/local/home/video/cernlib/2005/include/gen/imp64.inc" 1
*
* $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
*
* $Log: imp64.inc,v $
* Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
* Mathlib gen
*
*
* imp64.inc
*







      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c# 13 "dilog64.F" 2




      DIMENSION C(0:19)

      PARAMETER (Z1 = 1, HF = Z1/2)
      PARAMETER (PI = 3.14159 26535 89793 24D0)
      PARAMETER (PI3 = PI**2/3, PI6 = PI**2/6, PI12 = PI**2/12)

      DATA C( 0) / 0.42996 69356 08136 97D0/
      DATA C( 1) / 0.40975 98753 30771 05D0/
      DATA C( 2) /-0.01858 84366 50145 92D0/
      DATA C( 3) / 0.00145 75108 40622 68D0/
      DATA C( 4) /-0.00014 30418 44423 40D0/
      DATA C( 5) / 0.00001 58841 55418 80D0/
      DATA C( 6) /-0.00000 19078 49593 87D0/
      DATA C( 7) / 0.00000 02419 51808 54D0/
      DATA C( 8) /-0.00000 00319 33412 74D0/
      DATA C( 9) / 0.00000 00043 45450 63D0/
      DATA C(10) /-0.00000 00006 05784 80D0/
      DATA C(11) / 0.00000 00000 86120 98D0/
      DATA C(12) /-0.00000 00000 12443 32D0/
      DATA C(13) / 0.00000 00000 01822 56D0/
      DATA C(14) /-0.00000 00000 00270 07D0/
      DATA C(15) / 0.00000 00000 00040 42D0/
      DATA C(16) /-0.00000 00000 00006 10D0/
      DATA C(17) / 0.00000 00000 00000 93D0/
      DATA C(18) /-0.00000 00000 00000 14D0/
      DATA C(19) /+0.00000 00000 00000 02D0/

      IF(X .EQ. 1) THEN
       H=PI6
      ELSEIF(X .EQ. -1) THEN
       H=-PI12
      ELSE
       T=-X
       IF(T .LE. -2) THEN
        Y=-1/(1+T)
        S=1
        A=-PI3+HF*(LOG(-T)**2-LOG(1+1/T)**2)
       ELSEIF(T .LT. -1) THEN
        Y=-1-T
        S=-1
        A=LOG(-T)
        A=-PI6+A*(A+LOG(1+1/T))
       ELSE IF(T .LE. -HF) THEN
        Y=-(1+T)/T
        S=1
        A=LOG(-T)
        A=-PI6+A*(-HF*A+LOG(1+T))
       ELSE IF(T .LT. 0) THEN
        Y=-T/(1+T)
        S=-1
        A=HF*LOG(1+T)**2
       ELSE IF(T .LE. 1) THEN
        Y=T
        S=1
        A=0
       ELSE
        Y=1/T
        S=-1
        A=PI6+HF*LOG(T)**2
       ENDIF
       H=Y+Y-1
       ALFA=H+H
       B1=0
       B2=0
       DO 1 I = 19,0,-1
       B0=C(I)+ALFA*B1-B2
       B2=B1
    1  B1=B0
       H=-(S*(B0-H*B2)+A)
      ENDIF

      DDILOG=H




      RETURN
      END


*
* $Id: dzero.F,v 1.1.1.1 1996/02/15 17:49:07 mclareni Exp $
*
* $Log: dzero.F,v $
* Revision 1.1.1.1  1996/02/15 17:49:07  mclareni
* Kernlib
*










































      SUBROUTINE DZERO(A,B,X0,R,EPS,MXF,F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
* $Id: c205body.inc,v 1.1.1.1 1996/02/15 17:49:07 mclareni Exp $
*
* $Log: c205body.inc,v $
* Revision 1.1.1.1  1996/02/15 17:49:07  mclareni
* Kernlib
*
*


*
*
* c205body.inc
*
      LOGICAL MFLAG,RFLAG
      EXTERNAL F
 
      PARAMETER (ONE = 1, HALF = ONE/2)
 
      XA=MIN(A,B)
      XB=MAX(A,B)
      FA=F(XA,1)
      FB=F(XB,2)
      IF(FA*FB .GT. 0) GO TO 5
      MC=0
 
    1 X0=HALF*(XA+XB)
      R=X0-XA
      EE=EPS*(ABS(X0)+1)
      IF(R .LE. EE) GO TO 4
      F1=FA
      X1=XA
      F2=FB
      X2=XB
 
    2 FX=F(X0,2)
      MC=MC+1
      IF(MC .GT. MXF) GO TO 6
      IF(FX*FA .GT. 0) THEN
       XA=X0
       FA=FX
      ELSE
       XB=X0
       FB=FX
      END IF
 
    3 U1=F1-F2
      U2=X1-X2
      U3=F2-FX
      U4=X2-X0
      IF(U2 .EQ. 0 .OR. U4 .EQ. 0) GO TO 1
      F3=FX
      X3=X0
      U1=U1/U2
      U2=U3/U4
      CA=U1-U2
      CB=(X1+X2)*U2-(X2+X0)*U1
      CC=(X1-X0)*F1-X1*(CA*X1+CB)
      IF(CA .EQ. 0) THEN
       IF(CB .EQ. 0) GO TO 1
       X0=-CC/CB
      ELSE
       U3=CB/(2*CA)
       U4=U3*U3-CC/CA
       IF(U4 .LT. 0) GO TO 1
       X0=-U3+SIGN(SQRT(U4),X0+U3)
      END IF
      IF(X0 .LT. XA .OR. X0 .GT. XB) GO TO 1
 
      R=MIN(ABS(X0-X3),ABS(X0-X2))
      EE=EPS*(ABS(X0)+1)
      IF(R .GT. EE) THEN
       F1=F2
       X1=X2
       F2=F3
       X2=X3
       GO TO 2
      END IF
 
      FX=F(X0,2)
      IF(FX .EQ. 0) GO TO 4
      IF(FX*FA .LT. 0) THEN
       XX=X0-EE
       IF(XX .LE. XA) GO TO 4
       FF=F(XX,2)
       FB=FF
       XB=XX
      ELSE
       XX=X0+EE
       IF(XX .GE. XB) GO TO 4
       FF=F(XX,2)
       FA=FF
       XA=XX
      END IF
      IF(FX*FF .GT. 0) THEN
       MC=MC+2
       IF(MC .GT. MXF) GO TO 6
       F1=F3
       X1=X3
       F2=FX
       X2=X0
       X0=XX
       FX=FF
       GO TO 3
      END IF
 
    4 R=EE
      FF=F(X0,3)
      RETURN
    5 CALL KERMTR('C205.1',LGFILE,MFLAG,RFLAG)
      IF(MFLAG) THEN
       IF(LGFILE .EQ. 0) WRITE(*,100)
       IF(LGFILE .NE. 0) WRITE(LGFILE,100)
      END IF
      IF(.NOT.RFLAG) CALL ABEND
      R=-2*(XB-XA)
      X0=0
      RETURN
    6 CALL KERMTR('C205.2',LGFILE,MFLAG,RFLAG)
      IF(MFLAG) THEN
       IF(LGFILE .EQ. 0) WRITE(*,101)
       IF(LGFILE .NE. 0) WRITE(LGFILE,101)
      END IF
      IF(.NOT.RFLAG) CALL ABEND
      R=-HALF*ABS(XB-XA)
      X0=0
      RETURN


  100 FORMAT(1X,'***** CERN C205 DZERO ... F(A) AND F(B)',
     1          ' HAVE THE SAME SIGN')
  101 FORMAT(1X,'***** CERN C205 DZERO ... TOO MANY FUNCTION CALLS')
      END


c# 1 "kerset.F"
c# 1 "<built-in>"
c# 1 "<command line>"
c# 1 "kerset.F"
*
* $Id: kerset.F,v 1.1.1.1 1996/02/15 17:48:35 mclareni Exp $
*
* $Log: kerset.F,v $
* Revision 1.1.1.1  1996/02/15 17:48:35  mclareni
* Kernlib
*
*
c# 1 "/usr/local/home/video/cernlib/2005/include/kernnum/pilot.h" 1
c# 21 "/usr/local/home/video/cernlib/2005/include/kernnum/pilot.h"

c# 33 "/usr/local/home/video/cernlib/2005/include/kernnum/pilot.h"

c# 10 "kerset.F" 2
          SUBROUTINE KERSET(ERCODE,LGFILE,LIMITM,LIMITR)
                    PARAMETER(KOUNTE  =  27)
          CHARACTER*6         ERCODE,   CODE(KOUNTE)
          LOGICAL             MFLAG,    RFLAG
          INTEGER             KNTM(KOUNTE),       KNTR(KOUNTE)
          DATA      LOGF      /  0  /
          DATA      CODE(1), KNTM(1), KNTR(1)  / 'C204.1', 255, 255 /
          DATA      CODE(2), KNTM(2), KNTR(2)  / 'C204.2', 255, 255 /
          DATA      CODE(3), KNTM(3), KNTR(3)  / 'C204.3', 255, 255 /
          DATA      CODE(4), KNTM(4), KNTR(4)  / 'C205.1', 255, 255 /
          DATA      CODE(5), KNTM(5), KNTR(5)  / 'C205.2', 255, 255 /
          DATA      CODE(6), KNTM(6), KNTR(6)  / 'C305.1', 255, 255 /
          DATA      CODE(7), KNTM(7), KNTR(7)  / 'C308.1', 255, 255 /
          DATA      CODE(8), KNTM(8), KNTR(8)  / 'C312.1', 255, 255 /
          DATA      CODE(9), KNTM(9), KNTR(9)  / 'C313.1', 255, 255 /
          DATA      CODE(10),KNTM(10),KNTR(10) / 'C336.1', 255, 255 /
          DATA      CODE(11),KNTM(11),KNTR(11) / 'C337.1', 255, 255 /
          DATA      CODE(12),KNTM(12),KNTR(12) / 'C341.1', 255, 255 /
          DATA      CODE(13),KNTM(13),KNTR(13) / 'D103.1', 255, 255 /
          DATA      CODE(14),KNTM(14),KNTR(14) / 'D106.1', 255, 255 /
          DATA      CODE(15),KNTM(15),KNTR(15) / 'D209.1', 255, 255 /
          DATA      CODE(16),KNTM(16),KNTR(16) / 'D509.1', 255, 255 /
          DATA      CODE(17),KNTM(17),KNTR(17) / 'E100.1', 255, 255 /
          DATA      CODE(18),KNTM(18),KNTR(18) / 'E104.1', 255, 255 /
          DATA      CODE(19),KNTM(19),KNTR(19) / 'E105.1', 255, 255 /
          DATA      CODE(20),KNTM(20),KNTR(20) / 'E208.1', 255, 255 /
          DATA      CODE(21),KNTM(21),KNTR(21) / 'E208.2', 255, 255 /
          DATA      CODE(22),KNTM(22),KNTR(22) / 'F010.1', 255,   0 /
          DATA      CODE(23),KNTM(23),KNTR(23) / 'F011.1', 255,   0 /
          DATA      CODE(24),KNTM(24),KNTR(24) / 'F012.1', 255,   0 /
          DATA      CODE(25),KNTM(25),KNTR(25) / 'F406.1', 255,   0 /
          DATA      CODE(26),KNTM(26),KNTR(26) / 'G100.1', 255, 255 /
          DATA      CODE(27),KNTM(27),KNTR(27) / 'G100.2', 255, 255 /
          LOGF  =  LGFILE
             L  =  0
          IF(ERCODE .NE. ' ')  THEN
             DO 10  L = 1, 6
                IF(ERCODE(1:L) .EQ. ERCODE)  GOTO 12
  10            CONTINUE
  12         CONTINUE
          ENDIF
          DO 14     I  =  1, KOUNTE
             IF(L .EQ. 0)  GOTO 13
             IF(CODE(I)(1:L) .NE. ERCODE(1:L))  GOTO 14
  13         IF(LIMITM.GE.0) KNTM(I)  =  LIMITM
             IF(LIMITR.GE.0) KNTR(I)  =  LIMITR
  14         CONTINUE
          RETURN
          ENTRY KERMTR(ERCODE,LOG,MFLAG,RFLAG)
          LOG  =  LOGF
          DO 20     I  =  1, KOUNTE
             IF(ERCODE .EQ. CODE(I))  GOTO 21
  20         CONTINUE
          WRITE(*,1000)  ERCODE
          CALL ABEND
          RETURN
  21      RFLAG  =  KNTR(I) .GE. 1
          IF(RFLAG  .AND.  (KNTR(I) .LT. 255))  KNTR(I)  =  KNTR(I) - 1
          MFLAG  =  KNTM(I) .GE. 1
          IF(MFLAG  .AND.  (KNTM(I) .LT. 255))  KNTM(I)  =  KNTM(I) - 1
          IF(.NOT. RFLAG)  THEN
             IF(LOGF .LT. 1)  THEN
                WRITE(*,1001)  CODE(I)
             ELSE
                WRITE(LOGF,1001)  CODE(I)
             ENDIF
          ENDIF
          IF(MFLAG .AND. RFLAG)  THEN
             IF(LOGF .LT. 1)  THEN
                WRITE(*,1002)  CODE(I)
             ELSE
                WRITE(LOGF,1002)  CODE(I)
             ENDIF
          ENDIF
          RETURN
1000      FORMAT(' KERNLIB LIBRARY ERROR. ' /
     +           ' ERROR CODE ',A6,' NOT RECOGNIZED BY KERMTR',
     +           ' ERROR MONITOR. RUN ABORTED.')
1001      FORMAT(/' ***** RUN TERMINATED BY CERN LIBRARY ERROR ',
     +           'CONDITION ',A6)
1002      FORMAT(/' ***** CERN LIBRARY ERROR CONDITION ',A6)
          END

c# 1 "rm48.F"
c# 1 "<built-in>"
c# 1 "<command line>"
c# 1 "rm48.F"
*
* $Id: rm48.F,v 1.2 1996/12/12 16:32:06 cernlib Exp $
*
* $Log: rm48.F,v $
* Revision 1.2  1996/12/12 16:32:06  cernlib
* Variables ONE and ZERO added to SAVE statement, courtesy R.Veenhof
*
* Revision 1.1.1.1  1996/04/01 15:02:55  mclareni
* Mathlib gen
*
*
c# 1 "/usr/local/home/video/cernlib/2005/include/gen/pilot.h" 1
























c# 40 "/usr/local/home/video/cernlib/2005/include/gen/pilot.h"

c# 57 "/usr/local/home/video/cernlib/2005/include/gen/pilot.h"



























































c# 13 "rm48.F" 2
      SUBROUTINE RM48(RVEC,LENV)
C     Double-precision version of
C Universal random number generator proposed by Marsaglia and Zaman
C in report FSU-SCRI-87-50
C        based on RANMAR, modified by F. James, to generate vectors
C        of pseudorandom numbers RVEC of length LENV, where the numbers
C        in RVEC are numbers with at least 48-bit mantissas.
C   Input and output entry points: RM48IN, RM48UT.
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C!!!  Calling sequences for RM48:                                    ++
C!!!      CALL RM48 (RVEC, LEN)     returns a vector RVEC of LEN     ++
C!!!                   64-bit random floating point numbers between  ++
C!!!                   zero and one.                                 ++
C!!!      CALL RM48IN(I1,N1,N2)   initializes the generator from one ++
C!!!                   64-bit integer I1, and number counts N1,N2    ++
C!!!                  (for initializing, set N1=N2=0, but to restart ++
C!!!                    a previously generated sequence, use values  ++ 
C!!!                    output by RM48UT)                            ++ 
C!!!      CALL RM48UT(I1,N1,N2)   outputs the value of the original  ++
C!!!                  seed and the two number counts, to be used     ++
C!!!                  for restarting by initializing to I1 and       ++  
C!!!                  skipping N2*100000000+N1 numbers.              ++
C!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C for 32-bit machines, use IMPLICIT DOUBLE PRECISION
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION RVEC(*)
      COMMON/R48ST1/U(97),C,I97,J97
      PARAMETER (MODCNS=1000000000)
      SAVE CD, CM, TWOM24, NTOT, NTOT2, IJKL,TWOM49, ONE, ZERO
      save /R48ST1/
      DATA NTOT,NTOT2,IJKL/-1,0,0/
C
      IF (NTOT .GE. 0)  GO TO 50
C
C        Default initialization. User has called RM48 without RM48IN.
      IJKL = 54217137
      NTOT = 0
      NTOT2 = 0
      KALLED = 0
      GO TO 1
C
      ENTRY      RM48IN(IJKLIN, NTOTIN,NTOT2N)
C         Initializing routine for RM48, may be called before
C         generating pseudorandom numbers with RM48.   The input
C         values should be in the ranges:  0<=IJKLIN<=900 OOO OOO
C                                          0<=NTOTIN<=999 999 999
C                                          0<=NTOT2N<<999 999 999!
C To get the standard values in Marsaglia's paper, IJKLIN=54217137
C                                            NTOTIN,NTOT2N=0
      IJKL = IJKLIN
      NTOT = MAX(NTOTIN,0)
      NTOT2= MAX(NTOT2N,0)
      KALLED = 1
C          always come here to initialize
    1 CONTINUE
      IJ = IJKL/30082
      KL = IJKL - 30082*IJ
      I = MOD(IJ/177, 177) + 2
      J = MOD(IJ, 177)     + 2
      K = MOD(KL/169, 178) + 1
      L = MOD(KL, 169)
      WRITE(6,'(A,I10,2X,2I10)') ' RM48 INITIALIZED:',IJKL,NTOT,NTOT2
CCC      PRINT '(A,4I10)', '   I,J,K,L= ',I,J,K,L
      ONE = 1.
      HALF = 0.5
      ZERO = 0.
      DO 2 II= 1, 97
      S = 0.
      T = HALF
      DO 3 JJ= 1, 48
         M = MOD(MOD(I*J,179)*K, 179)
         I = J
         J = K
         K = M
         L = MOD(53*L+1, 169)
         IF (MOD(L*M,64) .GE. 32)  S = S+T
    3    T = HALF*T
    2 U(II) = S
      TWOM49 = T
      TWOM24 = ONE
      DO 4 I24= 1, 24
    4 TWOM24 = HALF*TWOM24
      C  =   362436.*TWOM24
      CD =  7654321.*TWOM24
      CM = 16777213.*TWOM24
      I97 = 97
      J97 = 33
C       Complete initialization by skipping
C            (NTOT2*MODCNS + NTOT) random numbers
      DO 45 LOOP2= 1, NTOT2+1
      NOW = MODCNS
      IF (LOOP2 .EQ. NTOT2+1)  NOW=NTOT
      IF (NOW .GT. 0)  THEN
      WRITE(6,'(A,I15)') ' RM48IN SKIPPING OVER ',NOW
c Bug! fixed by P. Nason, 13/12/2012
c          DO 40 IDUM = 1, NTOT
          DO 40 IDUM = 1, NOW
          UNI = U(I97)-U(J97)
          IF (UNI .LT. ZERO)  UNI=UNI+ONE
          U(I97) = UNI
          I97 = I97-1
          IF (I97 .EQ. 0)  I97=97
          J97 = J97-1
          IF (J97 .EQ. 0)  J97=97
          C = C - CD
          IF (C .LT. ZERO)  C=C+CM
   40     CONTINUE
      ENDIF
   45 CONTINUE
      IF (KALLED .EQ. 1)  RETURN
C
C          Normal entry to generate LENV random numbers
   50 CONTINUE
      DO 100 IVEC= 1, LENV
      UNI = U(I97)-U(J97)
      IF (UNI .LT. ZERO)  UNI=UNI+ONE
      U(I97) = UNI
      I97 = I97-1
      IF (I97 .EQ. 0)  I97=97
      J97 = J97-1
      IF (J97 .EQ. 0)  J97=97
      C = C - CD
      IF (C .LT. ZERO)  C=C+CM
      UNI = UNI-C
      IF (UNI .LT. ZERO) UNI=UNI+ONE
      RVEC(IVEC) = UNI
C             Replace exact zeros by 2**-49
         IF (UNI .EQ. ZERO)  THEN
            RVEC(IVEC) = TWOM49
         ENDIF
  100 CONTINUE
      NTOT = NTOT + LENV
         IF (NTOT .GE. MODCNS)  THEN
         NTOT2 = NTOT2 + 1
         NTOT = NTOT - MODCNS
         ENDIF
      RETURN
C           Entry to output current status
      ENTRY RM48UT(IJKLUT,NTOTUT,NTOT2T)
      IJKLUT = IJKL
      NTOTUT = NTOT
      NTOT2T = NTOT2
      RETURN
      END


c# 1 "abend.F"
c# 1 "<built-in>"
c# 1 "<command line>"
c# 1 "abend.F"
*
* $Id: abend.F,v 1.1.1.1 1996/02/15 17:50:37 mclareni Exp $
*
* $Log: abend.F,v $
* Revision 1.1.1.1  1996/02/15 17:50:37  mclareni
* Kernlib
*
*
c# 1 "/usr/local/home/video/cernlib/2005/include/kerngen/pilot.h" 1
c# 94 "/usr/local/home/video/cernlib/2005/include/kerngen/pilot.h"





c# 10 "abend.F" 2





      SUBROUTINE ABEND
C
C CERN PROGLIB# Z035    ABEND           .VERSION KERNFOR  4.31  911111
C ORIG.  8/02/88  JZ
C

      STOP  7
      END
c# 1 "lenocc.F"
c# 1 "<built-in>"
c# 1 "<command line>"
c# 1 "lenocc.F"
*
* $Id: lenocc.F,v 1.1.1.1 1996/02/15 17:49:49 mclareni Exp $
*
* $Log: lenocc.F,v $
* Revision 1.1.1.1  1996/02/15 17:49:49  mclareni
* Kernlib
*
*
c# 1 "/usr/local/home/video/cernlib/2005/include/kerngen/pilot.h" 1
c# 94 "/usr/local/home/video/cernlib/2005/include/kerngen/pilot.h"





c# 10 "lenocc.F" 2
      FUNCTION LENOCC (CHV)
C
C CERN PROGLIB# M507    LENOCC          .VERSION KERNFOR  4.21  890323
C ORIG. March 85, A.Petrilli, re-write 21/02/89, JZ
C
C-    Find last non-blank character in CHV

      CHARACTER    CHV*(*)

      N = LEN(CHV)

      DO 17  JJ= N,1,-1
      IF (CHV(JJ:JJ).NE.' ') GO TO 99
   17 CONTINUE
      JJ = 0

   99 LENOCC = JJ
      RETURN
      END
c# 1 "mtlset.F"
c# 1 "<built-in>"
c# 1 "<command line>"
c# 1 "mtlset.F"
*
* $Id: mtlset.F,v 1.1.1.1 1996/04/01 15:02:53 mclareni Exp $
*
* $Log: mtlset.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:53  mclareni
* Mathlib gen
*
*
c# 1 "/usr/local/home/video/cernlib/2005/include/gen/pilot.h" 1
























c# 40 "/usr/local/home/video/cernlib/2005/include/gen/pilot.h"

c# 57 "/usr/local/home/video/cernlib/2005/include/gen/pilot.h"



























































c# 10 "mtlset.F" 2
      SUBROUTINE MTLSET(ERC,NLG,MXM,MXR)

      PARAMETER (KTE = 132)
      CHARACTER*6 ERC,CODE(KTE)
      LOGICAL LMF,LRF
      DIMENSION KNTM(KTE),KNTR(KTE)

      DATA ILG /0/

C     renumber the data statements after putting new codes in Unix with:
C     awk -F'[()]' '{ printf"%s(%s)%s(%s)%s(%s)%s\n",$1,NR,$3,NR,$5,NR,$7 }'
C     and modify KTE to the number of lines below

      DATA CODE(1),KNTM(1),KNTR(1) / 'B100.1', 255, 255 /
      DATA CODE(2),KNTM(2),KNTR(2) / 'B300.1', 255, 255 /
      DATA CODE(3),KNTM(3),KNTR(3) / 'B300.2', 255, 255 /
      DATA CODE(4),KNTM(4),KNTR(4) / 'C200.0', 255, 255 /
      DATA CODE(5),KNTM(5),KNTR(5) / 'C200.1', 255, 255 /
      DATA CODE(6),KNTM(6),KNTR(6) / 'C200.2', 255, 255 /
      DATA CODE(7),KNTM(7),KNTR(7) / 'C200.3', 255, 255 /
      DATA CODE(8),KNTM(8),KNTR(8) / 'C201.0', 255, 255 /
      DATA CODE(9),KNTM(9),KNTR(9) / 'C202.0', 255, 255 /
      DATA CODE(10),KNTM(10),KNTR(10) / 'C202.1', 255, 255 /
      DATA CODE(11),KNTM(11),KNTR(11) / 'C202.2', 255, 255 /
      DATA CODE(12),KNTM(12),KNTR(12) / 'C205.1', 255, 255 /
      DATA CODE(13),KNTM(13),KNTR(13) / 'C205.2', 255, 255 /
      DATA CODE(14),KNTM(14),KNTR(14) / 'C207.0', 255, 255 /
      DATA CODE(15),KNTM(15),KNTR(15) / 'C208.0', 255, 255 /
      DATA CODE(16),KNTM(16),KNTR(16) / 'C209.0', 255, 255 /
      DATA CODE(17),KNTM(17),KNTR(17) / 'C209.1', 255, 255 /
      DATA CODE(18),KNTM(18),KNTR(18) / 'C209.2', 255, 255 /
      DATA CODE(19),KNTM(19),KNTR(19) / 'C209.3', 255, 255 /
      DATA CODE(20),KNTM(20),KNTR(20) / 'C210.1', 255, 255 /
      DATA CODE(21),KNTM(21),KNTR(21) / 'C302.1', 255, 255 /
      DATA CODE(22),KNTM(22),KNTR(22) / 'C303.1', 255, 255 /
      DATA CODE(23),KNTM(23),KNTR(23) / 'C304.1', 255, 255 /
      DATA CODE(24),KNTM(24),KNTR(24) / 'C305.1', 255, 255 /
      DATA CODE(25),KNTM(25),KNTR(25) / 'C306.1', 255, 255 /
      DATA CODE(26),KNTM(26),KNTR(26) / 'C307.1', 255, 255 /
      DATA CODE(27),KNTM(27),KNTR(27) / 'C312.1', 255, 255 /
      DATA CODE(28),KNTM(28),KNTR(28) / 'C313.1', 255, 255 /
      DATA CODE(29),KNTM(29),KNTR(29) / 'C315.1', 255, 255 /
      DATA CODE(30),KNTM(30),KNTR(30) / 'C316.1', 255, 255 /
      DATA CODE(31),KNTM(31),KNTR(31) / 'C316.2', 255, 255 /
      DATA CODE(32),KNTM(32),KNTR(32) / 'C320.1', 255, 255 /
      DATA CODE(33),KNTM(33),KNTR(33) / 'C321.1', 255, 255 /
      DATA CODE(34),KNTM(34),KNTR(34) / 'C323.1', 255, 255 /
      DATA CODE(35),KNTM(35),KNTR(35) / 'C327.1', 255, 255 /
      DATA CODE(36),KNTM(36),KNTR(36) / 'C328.1', 255, 255 /
      DATA CODE(37),KNTM(37),KNTR(37) / 'C328.2', 255, 255 /
      DATA CODE(38),KNTM(38),KNTR(38) / 'C328.3', 255, 255 /
      DATA CODE(39),KNTM(39),KNTR(39) / 'C330.1', 255, 255 /
      DATA CODE(40),KNTM(40),KNTR(40) / 'C330.2', 255, 255 /
      DATA CODE(41),KNTM(41),KNTR(41) / 'C330.3', 255, 255 /
      DATA CODE(42),KNTM(42),KNTR(42) / 'C331.1', 255, 255 /
      DATA CODE(43),KNTM(43),KNTR(43) / 'C331.2', 255, 255 /
      DATA CODE(44),KNTM(44),KNTR(44) / 'C334.1', 255, 255 /
      DATA CODE(45),KNTM(45),KNTR(45) / 'C334.2', 255, 255 /
      DATA CODE(46),KNTM(46),KNTR(46) / 'C334.3', 255, 255 /
      DATA CODE(47),KNTM(47),KNTR(47) / 'C334.4', 255, 255 /
      DATA CODE(48),KNTM(48),KNTR(48) / 'C334.5', 255, 255 /
      DATA CODE(49),KNTM(49),KNTR(49) / 'C334.6', 255, 255 /
      DATA CODE(50),KNTM(50),KNTR(50) / 'C336.1', 255, 255 /
      DATA CODE(51),KNTM(51),KNTR(51) / 'C337.1', 255, 255 /
      DATA CODE(52),KNTM(52),KNTR(52) / 'C338.1', 255, 255 /
      DATA CODE(53),KNTM(53),KNTR(53) / 'C340.1', 255, 255 /
      DATA CODE(54),KNTM(54),KNTR(54) / 'C343.1', 255, 255 /
      DATA CODE(55),KNTM(55),KNTR(55) / 'C343.2', 255, 255 /
      DATA CODE(56),KNTM(56),KNTR(56) / 'C343.3', 255, 255 /
      DATA CODE(57),KNTM(57),KNTR(57) / 'C343.4', 255, 255 /
      DATA CODE(58),KNTM(58),KNTR(58) / 'C344.1', 255, 255 /
      DATA CODE(59),KNTM(59),KNTR(59) / 'C344.2', 255, 255 /
      DATA CODE(60),KNTM(60),KNTR(60) / 'C344.3', 255, 255 /
      DATA CODE(61),KNTM(61),KNTR(61) / 'C344.4', 255, 255 /
      DATA CODE(62),KNTM(62),KNTR(62) / 'C345.1', 255, 255 /
      DATA CODE(63),KNTM(63),KNTR(63) / 'C346.1', 255, 255 /
      DATA CODE(64),KNTM(64),KNTR(64) / 'C346.2', 255, 255 /
      DATA CODE(65),KNTM(65),KNTR(65) / 'C346.3', 255, 255 /
      DATA CODE(66),KNTM(66),KNTR(66) / 'C347.1', 255, 255 /
      DATA CODE(67),KNTM(67),KNTR(67) / 'C347.2', 255, 255 /
      DATA CODE(68),KNTM(68),KNTR(68) / 'C347.3', 255, 255 /
      DATA CODE(69),KNTM(69),KNTR(69) / 'C347.4', 255, 255 /
      DATA CODE(70),KNTM(70),KNTR(70) / 'C347.5', 255, 255 /
      DATA CODE(71),KNTM(71),KNTR(71) / 'C347.6', 255, 255 /
      DATA CODE(72),KNTM(72),KNTR(72) / 'C348.1', 255, 255 /
      DATA CODE(73),KNTM(73),KNTR(73) / 'C349.1', 255, 255 /
      DATA CODE(74),KNTM(74),KNTR(74) / 'C349.2', 255, 255 /
      DATA CODE(75),KNTM(75),KNTR(75) / 'C349.3', 255, 255 /
      DATA CODE(76),KNTM(76),KNTR(76) / 'D101.1', 255, 255 /
      DATA CODE(77),KNTM(77),KNTR(77) / 'D103.1', 255, 255 /
      DATA CODE(78),KNTM(78),KNTR(78) / 'D104.1', 255, 255 /
      DATA CODE(79),KNTM(79),KNTR(79) / 'D104.2', 255, 255 /
      DATA CODE(80),KNTM(80),KNTR(80) / 'D105.1', 255, 255 /
      DATA CODE(81),KNTM(81),KNTR(81) / 'D105.2', 255, 255 /
      DATA CODE(82),KNTM(82),KNTR(82) / 'D107.1', 255, 255 /
      DATA CODE(83),KNTM(83),KNTR(83) / 'D110.0', 255, 255 /
      DATA CODE(84),KNTM(84),KNTR(84) / 'D110.1', 255, 255 /
      DATA CODE(85),KNTM(85),KNTR(85) / 'D110.2', 255, 255 /
      DATA CODE(86),KNTM(86),KNTR(86) / 'D110.3', 255, 255 /
      DATA CODE(87),KNTM(87),KNTR(87) / 'D110.4', 255, 255 /
      DATA CODE(88),KNTM(88),KNTR(88) / 'D110.5', 255, 255 /
      DATA CODE(89),KNTM(89),KNTR(89) / 'D110.6', 255, 255 /
      DATA CODE(90),KNTM(90),KNTR(90) / 'D113.1', 255, 255 /
      DATA CODE(91),KNTM(91),KNTR(91) / 'D201.1', 255, 255 /
      DATA CODE(92),KNTM(92),KNTR(92) / 'D202.1', 255, 255 /
      DATA CODE(93),KNTM(93),KNTR(93) / 'D401.1', 255, 255 /
      DATA CODE(94),KNTM(94),KNTR(94) / 'D601.1', 255, 255 /
      DATA CODE(95),KNTM(95),KNTR(95) / 'E210.1', 255, 255 /
      DATA CODE(96),KNTM(96),KNTR(96) / 'E210.2', 255, 255 /
      DATA CODE(97),KNTM(97),KNTR(97) / 'E210.3', 255, 255 /
      DATA CODE(98),KNTM(98),KNTR(98) / 'E210.4', 255, 255 /
      DATA CODE(99),KNTM(99),KNTR(99) / 'E210.5', 255, 255 /
      DATA CODE(100),KNTM(100),KNTR(100) / 'E210.6', 255, 255 /
      DATA CODE(101),KNTM(101),KNTR(101) / 'E210.7', 255, 255 /
      DATA CODE(102),KNTM(102),KNTR(102) / 'E211.0', 255, 255 /
      DATA CODE(103),KNTM(103),KNTR(103) / 'E211.1', 255, 255 /
      DATA CODE(104),KNTM(104),KNTR(104) / 'E211.2', 255, 255 /
      DATA CODE(105),KNTM(105),KNTR(105) / 'E211.3', 255, 255 /
      DATA CODE(106),KNTM(106),KNTR(106) / 'E211.4', 255, 255 /
      DATA CODE(107),KNTM(107),KNTR(107) / 'E406.0', 255, 255 /
      DATA CODE(108),KNTM(108),KNTR(108) / 'E406.1', 255, 255 /
      DATA CODE(109),KNTM(109),KNTR(109) / 'E407.0', 255, 255 /
      DATA CODE(110),KNTM(110),KNTR(110) / 'E408.0', 255, 255 /
      DATA CODE(111),KNTM(111),KNTR(111) / 'E408.1', 255, 255 /
      DATA CODE(112),KNTM(112),KNTR(112) / 'F500.0', 255, 255 /
      DATA CODE(113),KNTM(113),KNTR(113) / 'F500.1', 255, 255 /
      DATA CODE(114),KNTM(114),KNTR(114) / 'F500.2', 255, 255 /
      DATA CODE(115),KNTM(115),KNTR(115) / 'F500.3', 255, 255 /
      DATA CODE(116),KNTM(116),KNTR(116) / 'G100.1', 255, 255 /
      DATA CODE(117),KNTM(117),KNTR(117) / 'G100.2', 255, 255 /
      DATA CODE(118),KNTM(118),KNTR(118) / 'G101.1', 255, 255 /
      DATA CODE(119),KNTM(119),KNTR(119) / 'G101.2', 255, 255 /
      DATA CODE(120),KNTM(120),KNTR(120) / 'G105.1', 255, 255 /
      DATA CODE(121),KNTM(121),KNTR(121) / 'G106.1', 255, 255 /
      DATA CODE(122),KNTM(122),KNTR(122) / 'G106.2', 255, 255 /
      DATA CODE(123),KNTM(123),KNTR(123) / 'G116.1', 255, 255 /
      DATA CODE(124),KNTM(124),KNTR(124) / 'G116.2', 255, 255 /
      DATA CODE(125),KNTM(125),KNTR(125) / 'H101.0', 255, 255 /
      DATA CODE(126),KNTM(126),KNTR(126) / 'H101.1', 255, 255 /
      DATA CODE(127),KNTM(127),KNTR(127) / 'H101.2', 255, 255 /
      DATA CODE(128),KNTM(128),KNTR(128) / 'H301.1', 255, 255 /
      DATA CODE(129),KNTM(129),KNTR(129) / 'U501.1', 255, 255 /
      DATA CODE(130),KNTM(130),KNTR(130) / 'V202.1', 255, 255 /
      DATA CODE(131),KNTM(131),KNTR(131) / 'V202.2', 255, 255 /
      DATA CODE(132),KNTM(132),KNTR(132) / 'V202.3', 255, 255 /

c# 175 "mtlset.F"

      ILG=NLG
      L=0
      IF(ERC .NE. ' ') THEN
       DO 10 L = 1,6
       IF(ERC(1:L) .EQ. ERC) GOTO 12
   10  CONTINUE
   12  CONTINUE
      ENDIF
      DO 14 I = 1,KTE
      IF(L .EQ. 0 .OR. CODE(I)(1:L) .EQ. ERC(1:L)) THEN
       IF(MXM .GE. 0) KNTM(I)=MXM
       IF(MXR .GE. 0) KNTR(I)=MXR
      ENDIF
   14 CONTINUE
      RETURN

      ENTRY MTLMTR(ERC,MLG,LMF,LRF)

      MLG=ILG
      DO 20 I = 1,KTE
      IF(ERC .EQ. CODE(I))  GOTO 21
   20 CONTINUE
      WRITE(*,100) ERC
      CALL ABEND
      RETURN

   21 LMF=KNTM(I) .GE. 1
      LRF=KNTR(I) .GE. 1
      IF(LMF .AND. KNTM(I) .LT. 255)  KNTM(I)=KNTM(I)-1
      IF(LRF .AND. KNTR(I) .LT. 255)  KNTR(I)=KNTR(I)-1
      IF(.NOT.LRF) THEN
       IF(ILG .LT. 1) WRITE(  *,101) CODE(I)
       IF(ILG .GE. 1) WRITE(ILG,101) CODE(I)
      ENDIF
      RETURN
  100 FORMAT(7X,'***** CERN N002 MTLSET ... ERROR N002: ',
     1'ERROR CODE ',A6,' NOT RECOGNIZED BY ERROR MONITOR. RUN ABORTED.')
  101 FORMAT(7X,'***** CERN N002 MTLSET ... ERROR NOO2.1: ',
     1'RUN TERMINATED BY LIBRARY ERROR CONDITION ',A6)
      END
c# 1 "mtlprt.F"
c# 1 "<built-in>"
c# 1 "<command line>"
c# 1 "mtlprt.F"
*
* $Id: mtlprt.F,v 1.1.1.1 1996/04/01 15:02:52 mclareni Exp $
*
* $Log: mtlprt.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:52  mclareni
* Mathlib gen
*
*
c# 1 "/usr/local/home/video/cernlib/2005/include/gen/pilot.h" 1
























c# 40 "/usr/local/home/video/cernlib/2005/include/gen/pilot.h"

c# 57 "/usr/local/home/video/cernlib/2005/include/gen/pilot.h"



























































c# 10 "mtlprt.F" 2
      SUBROUTINE MTLPRT(NAME,ERC,TEXT)
      CHARACTER*(*) NAME,ERC,TEXT
      LOGICAL LMF,LRF

      IF(ERC(5:6).NE.'.0') THEN
        CALL MTLMTR(ERC,MLG,LMF,LRF)
      ELSE
        LMF=.TRUE.
        LRF=.FALSE.
      ENDIF
      IF(LMF) THEN
        LT=LENOCC(TEXT)
        IF(MLG .LT. 1) WRITE(  *,100) ERC(1:4),NAME,ERC,TEXT(1:LT)
        IF(MLG .GE. 1) WRITE(MLG,100) ERC(1:4),NAME,ERC,TEXT(1:LT)
      ENDIF
      IF(.NOT.LRF) CALL ABEND
      RETURN
100   FORMAT(7X,'***** CERN ',A,1X,A,' ERROR ',A,': ',A)
      END
c# 1 "permu.F"
c# 1 "<built-in>"
c# 1 "<command line>"
c# 1 "permu.F"
*
* $Id: permu.F,v 1.1.1.1 1996/04/01 15:02:57 mclareni Exp $
*
* $Log: permu.F,v $
* Revision 1.1.1.1  1996/04/01 15:02:57  mclareni
* Mathlib gen
*
*
c# 1 "/usr/local/home/video/cernlib/2005/include/gen/pilot.h" 1
























c# 40 "/usr/local/home/video/cernlib/2005/include/gen/pilot.h"

c# 57 "/usr/local/home/video/cernlib/2005/include/gen/pilot.h"



























































c# 10 "permu.F" 2
      SUBROUTINE PERMU(IA,N)
C
      CHARACTER*(*) NAME
      PARAMETER(NAME='PERMUT')
C
      DIMENSION IA(*)
      PARAMETER (IFD = 12)
      DIMENSION IFCT(0:IFD),IV(IFD+1)
      CHARACTER*80 ERRTXT

      DATA (IFCT(I),I=0,IFD) /1,1,2,6,24,120,720,5040,40320,362880,
     1                        3628800,39916800,479001600/

      IF(N .LE. 0) RETURN
      IF(IA(1) .EQ. 0) THEN
       DO 11 I = 1,N
   11  IA(I)=I
       IF(N .EQ. 1) IA(1)=0
      ELSE
       DO 12 K1 = N,2,-1
       K=K1
       IF(IA(K-1) .LT. IA(K)) GO TO 14
   12  CONTINUE
       IA(1)=0
       RETURN
   14  KN=K+N
       DO 15 L = K,KN/2
       IB=IA(KN-L)
       IA(KN-L)=IA(L)
   15  IA(L)=IB
       DO 16 L1 = K,N
       L=L1
       IF(IA(L) .GT. IA(K-1)) GO TO 17
   16  CONTINUE
   17  IB=IA(K-1)
       IA(K-1)=IA(L)
       IA(L)=IB
      ENDIF
      RETURN

      ENTRY PERMUT(NRP,N,IA)

      IF(N .LE. 0) RETURN
      IF(N .GT. IFD) THEN
       WRITE(ERRTXT,101) N
       CALL MTLPRT(NAME,'V202.1',ERRTXT)
      ELSEIF(NRP .GT. IFCT(N)) THEN
       IA(1)=0
       CALL MTLPRT(NAME,'V202.2',
     +             'PERMUTATION OUTSIDE LEXICON REQUESTED')
      ELSE
       DO 21 I = 1,N
   21  IV(I)=I
       IO=NRP-1
       DO 22 M = N-1,1,-1
       IN=IO/IFCT(M)+1
       IO=MOD(IO,IFCT(M))
       IA(N-M)=IV(IN)
       DO 23 I = IN,M
   23  IV(I)=IV(I+1)
   22  CONTINUE
       IA(N)=IV(1)
      END IF
      RETURN

      ENTRY COMBI(IA,N,J)

      IF(N .LE. 0 .OR. J .LE. 0) RETURN
      IF(J .GT. N) THEN
       WRITE(ERRTXT,103) J,N
       CALL MTLPRT(NAME,'V202.3',ERRTXT)
      ELSEIF(IA(1) .EQ. 0) THEN
       DO 31 I = 1,J
   31  IA(I)=I
       IA(J+1)=0
      ELSE
       DO 32 I1 = 1,N
       I=I1
       IF(IA(I+1) .NE. IA(I)+1) GO TO 33
   32  IA(I)=I
   33  IA(I)=IA(I)+1
       IF(IA(J) .EQ. N+1) IA(1)=0
      ENDIF
      RETURN
  101 FORMAT('N = ',I20,' TOO BIG')
  103 FORMAT('J = ',I5,' > ',I5,' = N IS NOT PERMITTED')
      END




c# 1 "sortzv.F"
c# 1 "<built-in>"
c# 1 "<command line>"
c# 1 "sortzv.F"
*
* $Id: sortzv.F,v 1.1.1.1 1996/02/15 17:49:50 mclareni Exp $
*
* $Log: sortzv.F,v $
* Revision 1.1.1.1  1996/02/15 17:49:50  mclareni
* Kernlib
*
*
c# 1 "/usr/local/home/video/cernlib/2005/include/kerngen/pilot.h" 1
c# 94 "/usr/local/home/video/cernlib/2005/include/kerngen/pilot.h"





c# 10 "sortzv.F" 2
      SUBROUTINE SORTZV (A,INDEX,N1,MODE,NWAY,NSORT)
C
C CERN PROGLIB# M101    SORTZV          .VERSION KERNFOR  3.15  820113
C ORIG. 02/10/75
C
      ! MK: added this line to adapt to newer compilers
      INTEGER A
      DIMENSION A(N1),INDEX(N1)
C
C
      N = N1
      IF (N.LE.0)            RETURN
      IF (NSORT.NE.0) GO TO 2
      DO 1 I=1,N
    1 INDEX(I)=I
C
    2 IF (N.EQ.1)            RETURN
      IF (MODE)    10,20,30
   10 CALL SORTTI(A,INDEX,N)
      GO TO 40
C
   20 CALL SORTTC(A,INDEX,N)
      GO TO 40
C
   30 CALL SORTTF (A,INDEX,N)
C
   40 IF (NWAY.EQ.0) GO TO 50
      N2 = N/2
      DO 41 I=1,N2
      ISWAP = INDEX(I)
      K = N+1-I
      INDEX(I) = INDEX(K)
   41 INDEX(K) = ISWAP
   50 RETURN
      END
*     ========================================
      SUBROUTINE SORTTF (A,INDEX,N1)
C
      ! MK: added this line to adapt to newer compilers 
      INTEGER A
      DIMENSION A(N1),INDEX(N1)
C
      N = N1
      DO 3 I1=2,N
      I3 = I1
      I33 = INDEX(I3)
      AI = A(I33)
    1 I2 = I3/2
      IF (I2) 3,3,2
    2 I22 = INDEX(I2)
      IF (AI.LE.A (I22)) GO TO 3
      INDEX (I3) = I22
      I3 = I2
      GO TO 1
    3 INDEX (I3) = I33
    4 I3 = INDEX (N)
      INDEX (N) = INDEX (1)
      AI = A(I3)
      N = N-1
      IF (N-1) 12,12,5
    5 I1 = 1
    6 I2 = I1 + I1
      IF (I2.LE.N) I22= INDEX(I2)
      IF (I2-N) 7,9,11
    7 I222 = INDEX (I2+1)
      IF (A(I22)-A(I222)) 8,9,9
    8 I2 = I2+1
      I22 = I222
    9 IF (AI-A(I22)) 10,11,11
   10 INDEX(I1) = I22
      I1 = I2
      GO TO 6
   11 INDEX (I1) = I3
      GO TO 4
   12 INDEX (1) = I3
      RETURN
      END
*     ========================================
      SUBROUTINE SORTTI (A,INDEX,N1)
C
      INTEGER A,AI
      DIMENSION A(N1),INDEX(N1)
C
      N = N1
      DO 3 I1=2,N
      I3 = I1
      I33 = INDEX(I3)
      AI = A(I33)
    1 I2 = I3/2
      IF (I2) 3,3,2
    2 I22 = INDEX(I2)
      IF (AI.LE.A (I22)) GO TO 3
      INDEX (I3) = I22
      I3 = I2
      GO TO 1
    3 INDEX (I3) = I33
    4 I3 = INDEX (N)
      INDEX (N) = INDEX (1)
      AI = A(I3)
      N = N-1
      IF (N-1) 12,12,5
    5 I1 = 1
    6 I2 = I1 + I1
      IF (I2.LE.N) I22= INDEX(I2)
      IF (I2-N) 7,9,11
    7 I222 = INDEX (I2+1)
      IF (A(I22)-A(I222)) 8,9,9
    8 I2 = I2+1
      I22 = I222
    9 IF (AI-A(I22)) 10,11,11
   10 INDEX(I1) = I22
      I1 = I2
      GO TO 6
   11 INDEX (I1) = I3
      GO TO 4
   12 INDEX (1) = I3
      RETURN
      END
*     ========================================
      SUBROUTINE SORTTC (A,INDEX,N1)
C
      INTEGER A,AI
      DIMENSION A(N1),INDEX(N1)
C
      N = N1
      DO 3 I1=2,N
      I3 = I1
      I33 = INDEX(I3)
      AI = A(I33)
    1 I2 = I3/2
      IF (I2) 3,3,2
    2 I22 = INDEX(I2)
      IF(ICMPCH(AI,A(I22)))3,3,21
   21 INDEX (I3) = I22
      I3 = I2
      GO TO 1
    3 INDEX (I3) = I33
    4 I3 = INDEX (N)
      INDEX (N) = INDEX (1)
      AI = A(I3)
      N = N-1
      IF (N-1) 12,12,5
    5 I1 = 1
    6 I2 = I1 + I1
      IF (I2.LE.N) I22= INDEX(I2)
      IF (I2-N) 7,9,11
    7 I222 = INDEX (I2+1)
      IF (ICMPCH(A(I22),A(I222))) 8,9,9
    8 I2 = I2+1
      I22 = I222
    9 IF (ICMPCH(AI,A(I22))) 10,11,11
   10 INDEX(I1) = I22
      I1 = I2
      GO TO 6
   11 INDEX (I1) = I3
      GO TO 4
   12 INDEX (1) = I3
      RETURN
      END
*     ========================================
      FUNCTION ICMPCH(IC1,IC2)
C     FUNCTION TO COMPARE TWO 4 CHARACTER EBCDIC STRINGS - IC1,IC2
C     ICMPCH=-1 IF HEX VALUE OF IC1 IS LESS THAN IC2
C     ICMPCH=0  IF HEX VALUES OF IC1 AND IC2 ARE THE SAME
C     ICMPCH=+1 IF HEX VALUES OF IC1 IS GREATER THAN IC2
      I1=IC1
      I2=IC2
      IF(I1.GE.0.AND.I2.GE.0)GOTO 40
      IF(I1.GE.0)GOTO 60
      IF(I2.GE.0)GOTO 80
      I1=-I1
      I2=-I2
      IF(I1-I2)80,70,60
 40   IF(I1-I2)60,70,80
 60   ICMPCH=-1
      RETURN
 70   ICMPCH=0
      RETURN
 80   ICMPCH=1
      RETURN
      END





c **********  THIS IS NOT A CERN ROUTINE!!  ***************
      logical function pwhg_isfinite(x)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      real * 8 x
c According to ieee standards, a NaN is the only real not
c satisfying x.eq.x.
      if (.not.(x.eq.x)) then
         pwhg_isfinite = .false.
         if(kn_jacborn.ne.0) then
            call increasecnt('NaN exception')
         else
            call increasecnt('NaN exception with kn_jacborn=0')
         endif
         return
      endif
c Put constraint to avoid denormals
      if(x.gt.1.or.x.lt.-1) then
         if (1/x.eq.0) then
            pwhg_isfinite = .false.
            call increasecnt('Inf exception')
            return
         endif
      endif
      pwhg_isfinite = .true.
      end

      subroutine initcnt
      implicit none
      integer maxnum
      parameter (maxnum=100)
      logical ini
      character * 100 keywords(maxnum)
      real * 8 counters(maxnum)
      integer ncounters
      common/ccounters/keywords,counters,ncounters
      data ini/.true./
      save ini
      if(ini) then
         ncounters=0
         ini=.false.
      endif
      end

      subroutine increasecnt(string)
      implicit none
      character *(*) string
      integer maxnum
      parameter (maxnum=100)
      character * 100 keywords(maxnum)
      real * 8 counters(maxnum)
      integer ncounters
      common/ccounters/keywords,counters,ncounters
      integer ini,j
      call initcnt
      do j=1,ncounters
         if(string.eq.keywords(j)) then
            counters(j)=counters(j)+1
            return
         endif
      enddo
c not found
      if(ncounters.eq.maxnum) then
         write(*,*) 'ERROR: increasecnt too many counters requested'
         stop
      endif
      ncounters=ncounters+1
      keywords(ncounters)=string
      counters(ncounters)=1
      end

      subroutine addtocnt(string,value)
      implicit none
      character *(*) string
      real * 8 value
      integer maxnum
      parameter (maxnum=100)
      character * 100 keywords(maxnum)
      real * 8 counters(maxnum)
      integer ncounters
      common/ccounters/keywords,counters,ncounters
      integer j
      call initcnt
      do j=1,ncounters
         if(string.eq.keywords(j)) then
            counters(j)=counters(j)+value
            return
         endif
      enddo
c not found
      if(ncounters.eq.maxnum) then
         write(*,*) 'ERROR: increasecnt too many counters requested'
         stop
      endif
      ncounters=ncounters+1
      keywords(ncounters)=string
      counters(ncounters)=value
      end

      subroutine resetcnt(string)
      implicit none
      character *(*) string
      integer maxnum
      parameter (maxnum=100)
      character * 100 keywords(maxnum)
      real * 8 counters(maxnum)
      integer ncounters
      common/ccounters/keywords,counters,ncounters
      integer j
      do j=1,ncounters
         if(string.eq.keywords(j)) then
            counters(j)=0
            return
         endif
      enddo
c not found
      if(ncounters.eq.maxnum) then
         write(*,*) 'ERROR: increasecnt too many counters requested'
         stop
      endif
      ncounters=ncounters+1
      keywords(ncounters)=string
      counters(ncounters)=0
      end

      subroutine incrcntrs(string,k)
      implicit none
      integer k
      character *(*) string
      character * 20 line
      line=string
      line(17:)=' '
      write(line(17:),'(i4)') k
      call increasecnt(line)
      end

      subroutine incrcntrs2(string,j,k)
      implicit none
      integer j,k
      character *(*) string
      character * 20 line
      line=string
      line(17:)=' '
      write(line(13:),'(i4)') j
      write(line(17:),'(i4)') k
      call increasecnt(line)
      end

      subroutine write_counters
      implicit none
      include 'pwhg_rnd.h'
      integer iun
      character * 20 pwgprefix
      integer lprefix
      common/cpwgprefix/pwgprefix,lprefix
      real * 8 powheginput
      integer stage
      character * 3 cst
      call newunit(iun)
      if(rnd_cwhichseed.eq.'none') then
         open(unit=iun,file=pwgprefix(1:lprefix)//'counters.dat'
     1     ,status='unknown')
      else
         stage=powheginput("#parallelstage")
         if(stage.gt.0) then
            write(cst(3:3),'(i1)') stage
            cst(1:2)='st'
            open(unit=iun,file=pwgprefix(1:lprefix)//'counters'//'-'
     1       //cst//'-'//rnd_cwhichseed//'.dat',status='unknown')
         else
            open(unit=iun,file=pwgprefix(1:lprefix)//'counters'//
     1           rnd_cwhichseed//'.dat',status='unknown')
         endif
      endif
      call printcnt(iun)
      close(iun)
      end




      subroutine printcnt(iun)
      implicit none
      integer iun
      integer maxnum
      parameter (maxnum=100)
      character * 100 keywords(maxnum)
      real * 8 counters(maxnum)
      integer ncounters
      common/ccounters/keywords,counters,ncounters
      integer j,k
      if(ncounters.eq.0) return
      write(iun,*)
      write(iun,*) 'Value of counters at end of run:'
      do j=1,ncounters
         k=100
 1       if(keywords(j)(k:k).eq.' ') then
            k=k-1
            goto 1
         endif
         write(iun,*) keywords(j)(1:k),' = ',counters(j)
      enddo
      write(iun,*)
      end

      function getcnt(string)
      implicit none
      real * 8 getcnt
      character *(*) string
      integer maxnum
      parameter (maxnum=100)
      character * 100 keywords(maxnum)
      real * 8 counters(maxnum)
      integer ncounters
      common/ccounters/keywords,counters,ncounters
      integer j
      do j=1,ncounters
         if(string.eq.keywords(j)) then
            getcnt=counters(j)
            return
         endif
      enddo
      getcnt=0
      end


      subroutine reset_timer
      implicit none
      integer nmax
      parameter (nmax=10)
      integer n,i1(nmax),i2(nmax),i3(nmax)
      common/pwhg_cmytimer/n,i1,i2,i3
      save /pwhg_cmytimer/
      logical ini
      data ini/.true./
      save ini
      if(ini) then
         n=0
         ini=.false.
      endif
      if(n.ge.nmax) then
         write(*,*) ' more than ',nmax,'reentrant calls of reset_timer'
         write(*,*) ' increase the parameter nmax in reset_timer'
         call pwhg_exit(-1)
      endif
      n=n+1
      call system_clock(i1(n),i2(n),i3(n))
      end

      subroutine get_timer(seconds)
      implicit none
      real * 8 seconds
      integer j1,j2,j3
      integer nmax
      parameter (nmax=10)
      integer n,i1(nmax),i2(nmax),i3(nmax)
      common/pwhg_cmytimer/n,i1,i2,i3
      save /pwhg_cmytimer/
      if(n.le.0) then
         write(*,*) ' ************ ERROR ******************'
         write(*,*) ' get_timer called with no previous set_timer call'
         call pwhg_exit(-1)
      endif
c j1: seconds*j2;
c j3: max value of j1; after that j1 restart from zero
c With gfortran, j3 corresponds roughly to 24 days.
      call system_clock(j1,j2,j3)
c the whole point is to make sure we are not wrapping around
      if(j1.ge.i1(n)) then
c No wrap around
         seconds=dble(j1-i1(n))/i2(n)
      else
c Assume we wrapped around
         seconds=dble(j1)/i2(n)+dble(i3(n)-i1(n))/i2(n)
      endif
      n=n-1
      end


*
* $Id: dgauss.f,v 1.2 2007/01/07 00:14:07 madgraph Exp $
*
* $Log: dgauss.f,v $
* Revision 1.2  2007/01/07 00:14:07  madgraph
* Merged version 4.1 to main branch
*
* Revision 1.1.2.1  2006/09/29 23:35:23  madgraph
* Introducing MLM matching
*
* Revision 1.1.1.1  1996/04/01 15:02:13  mclareni
* Mathlib gen
*
*
      FUNCTION DGAUSS(F,A,B,EPS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      CHARACTER NAME*(*)
      PARAMETER (NAME = 'DGAUSS')

*
* $Id: dgauss.f,v 1.2 2007/01/07 00:14:07 madgraph Exp $
*
* $Log: dgauss.f,v $
* Revision 1.2  2007/01/07 00:14:07  madgraph
* Merged version 4.1 to main branch
*
* Revision 1.1.2.1  2006/09/29 23:35:23  madgraph
* Introducing MLM matching
*
* Revision 1.1.1.1  1996/04/01 15:02:13  mclareni
* Mathlib gen
*
*
*
* gausscod.inc
*
      DIMENSION W(12),X(12)

      PARAMETER (Z1 = 1, HF = Z1/2, CST = 5*Z1/1000)

      DATA X( 1) /9.6028985649753623D-1/, W( 1) /1.0122853629037626D-1/
      DATA X( 2) /7.9666647741362674D-1/, W( 2) /2.2238103445337447D-1/
      DATA X( 3) /5.2553240991632899D-1/, W( 3) /3.1370664587788729D-1/
      DATA X( 4) /1.8343464249564980D-1/, W( 4) /3.6268378337836198D-1/
      DATA X( 5) /9.8940093499164993D-1/, W( 5) /2.7152459411754095D-2/
      DATA X( 6) /9.4457502307323258D-1/, W( 6) /6.2253523938647893D-2/
      DATA X( 7) /8.6563120238783174D-1/, W( 7) /9.5158511682492785D-2/
      DATA X( 8) /7.5540440835500303D-1/, W( 8) /1.2462897125553387D-1/
      DATA X( 9) /6.1787624440264375D-1/, W( 9) /1.4959598881657673D-1/
      DATA X(10) /4.5801677765722739D-1/, W(10) /1.6915651939500254D-1/
      DATA X(11) /2.8160355077925891D-1/, W(11) /1.8260341504492359D-1/
      DATA X(12) /9.5012509837637440D-2/, W(12) /1.8945061045506850D-1/

      H=0
      IF(B .EQ. A) GO TO 99
      CONST=CST/ABS(B-A)
      BB=A
    1 AA=BB
      BB=B
    2 C1=HF*(BB+AA)
      C2=HF*(BB-AA)
      S8=0
      DO 3 I = 1,4
      U=C2*X(I)
    3 S8=S8+W(I)*(F(C1+U)+F(C1-U))
      S16=0
      DO 4 I = 5,12
      U=C2*X(I)
    4 S16=S16+W(I)*(F(C1+U)+F(C1-U))
      S16=C2*S16
      IF(ABS(S16-C2*S8) .LE. EPS*(1+ABS(S16))) THEN
       H=H+S16
       IF(BB .NE. B) GO TO 1
      ELSE
       BB=C1
       IF(1+CONST*ABS(C2) .NE. 1) GO TO 2
       H=0
       WRITE(*,*) NAME,'ERROR: TOO HIGH ACCURACY REQUIRED'
       GO TO 99
      END IF
      

   99 DGAUSS=H
      RETURN
      END



      subroutine pwhg_exit(iret)
      integer iret
      real * 8 tmp,powheginput
      tmp=powheginput('print unused tokens')
      call write_counters
      call exit(iret)
      end

      SUBROUTINE DGSET(A,B,N,X,W)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL DGQUAD
      DIMENSION X(*),W(*)
      CALL D107D1(2,DGQUAD,A,B,N,X,W)
      RETURN
      END

      FUNCTION DGQUAD(F,A,B,N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION X(1),W(1)
      EXTERNAL F

      CALL D107D1(1,F,A,B,N,X,W)
      DGQUAD=X(1)
      RETURN
      END


      SUBROUTINE D107D1(MODE,F,A,B,N,X,W)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      CHARACTER*6 NAME(2)
      DATA NAME /'DGQUAD','DGSET'/
      DIMENSION X(*),W(*),KTBA(97),U(273),V(273)
      PARAMETER (Z1 = 1, HF = Z1/2)

      CHARACTER*80 ERRTXT

      DATA KTBA
     1/0,1,2,4,6,9,12,16,20,25,30,36,42,49,56,64,3*0,72,3*0,82,7*0,94,
     2 7*0,110,7*0,130,15*0,154,15*0,186,15*0,226,0/

C  N=2.
      DATA U(1)  /5.7735026918962576D-1/, V(1)  /1/
C  N=3.
      DATA U(2)  /7.7459666924148338D-1/, V(2)  /5.5555555555555556D-1/
      DATA U(3)  /0/                    , V(3)  /8.8888888888888889D-1/
C  N=4.
      DATA U(4)  /8.6113631159405258D-1/, V(4)  /3.4785484513745386D-1/
      DATA U(5)  /3.3998104358485626D-1/, V(5)  /6.5214515486254614D-1/
C  N=5.
      DATA U(6)  /9.0617984593866399D-1/, V(6)  /2.3692688505618909D-1/
      DATA U(7)  /5.3846931010568309D-1/, V(7)  /4.7862867049936647D-1/
      DATA U(8)  /0/,                     V(8)  /5.6888888888888889D-1/
C  N=6.
      DATA U(9)  /9.3246951420315203D-1/, V(9)  /1.7132449237917035D-1/
      DATA U(10) /6.6120938646626451D-1/, V(10) /3.6076157304813861D-1/
      DATA U(11) /2.3861918608319691D-1/, V(11) /4.6791393457269105D-1/
C  N=7.
      DATA U(12) /9.4910791234275852D-1/, V(12) /1.2948496616886969D-1/
      DATA U(13) /7.4153118559939444D-1/, V(13) /2.7970539148927667D-1/
      DATA U(14) /4.0584515137739717D-1/, V(14) /3.8183005050511894D-1/
      DATA U(15) /0/                    , V(15) /4.1795918367346939D-1/
C  N=8.
      DATA U(16) /9.6028985649753623D-1/, V(16) /1.0122853629037626D-1/
      DATA U(17) /7.9666647741362674D-1/, V(17) /2.2238103445337447D-1/
      DATA U(18) /5.2553240991632899D-1/, V(18) /3.1370664587788729D-1/
      DATA U(19) /1.8343464249564980D-1/, V(19) /3.6268378337836198D-1/
C  N=9.
      DATA U(20) /9.6816023950762609D-1/, V(20) /8.1274388361574412D-2/
      DATA U(21) /8.3603110732663579D-1/, V(21) /1.8064816069485740D-1/
      DATA U(22) /6.1337143270059040D-1/, V(22) /2.6061069640293546D-1/
      DATA U(23) /3.2425342340380893D-1/, V(23) /3.1234707704000284D-1/
      DATA U(24) /0/                    , V(24) /3.3023935500125976D-1/
C  N=10.
      DATA U(25) /9.7390652851717172D-1/, V(25) /6.6671344308688138D-2/
      DATA U(26) /8.6506336668898451D-1/, V(26) /1.4945134915058059D-1/
      DATA U(27) /6.7940956829902441D-1/, V(27) /2.1908636251598204D-1/
      DATA U(28) /4.3339539412924719D-1/, V(28) /2.6926671930999636D-1/
      DATA U(29) /1.4887433898163121D-1/, V(29) /2.9552422471475287D-1/
C  N=11.
      DATA U(30) /9.7822865814605699D-1/, V(30) /5.5668567116173666D-2/
      DATA U(31) /8.8706259976809530D-1/, V(31) /1.2558036946490462D-1/
      DATA U(32) /7.3015200557404932D-1/, V(32) /1.8629021092773425D-1/
      DATA U(33) /5.1909612920681182D-1/, V(33) /2.3319376459199048D-1/
      DATA U(34) /2.6954315595234497D-1/, V(34) /2.6280454451024666D-1/
      DATA U(35) /0/                    , V(35) /2.7292508677790063D-1/
C  N=12.
      DATA U(36) /9.8156063424671925D-1/, V(36) /4.7175336386511827D-2/
      DATA U(37) /9.0411725637047486D-1/, V(37) /1.0693932599531843D-1/
      DATA U(38) /7.6990267419430469D-1/, V(38) /1.6007832854334623D-1/
      DATA U(39) /5.8731795428661745D-1/, V(39) /2.0316742672306592D-1/
      DATA U(40) /3.6783149899818019D-1/, V(40) /2.3349253653835481D-1/
      DATA U(41) /1.2523340851146892D-1/, V(41) /2.4914704581340279D-1/
C  N=13.
      DATA U(42) /9.8418305471858815D-1/, V(42) /4.0484004765315880D-2/
      DATA U(43) /9.1759839922297797D-1/, V(43) /9.2121499837728448D-2/
      DATA U(44) /8.0157809073330991D-1/, V(44) /1.3887351021978724D-1/
      DATA U(45) /6.4234933944034022D-1/, V(45) /1.7814598076194574D-1/
      DATA U(46) /4.4849275103644685D-1/, V(46) /2.0781604753688850D-1/
      DATA U(47) /2.3045831595513479D-1/, V(47) /2.2628318026289724D-1/
      DATA U(48) /0/                    , V(48) /2.3255155323087391D-1/
C  N=14.
      DATA U(49) /9.8628380869681234D-1/, V(49) /3.5119460331751863D-2/
      DATA U(50) /9.2843488366357352D-1/, V(50) /8.0158087159760210D-2/
      DATA U(51) /8.2720131506976499D-1/, V(51) /1.2151857068790318D-1/
      DATA U(52) /6.8729290481168547D-1/, V(52) /1.5720316715819353D-1/
      DATA U(53) /5.1524863635815409D-1/, V(53) /1.8553839747793781D-1/
      DATA U(54) /3.1911236892788976D-1/, V(54) /2.0519846372129560D-1/
      DATA U(55) /1.0805494870734366D-1/, V(55) /2.1526385346315779D-1/
C  N=15.
      DATA U(56) /9.8799251802048543D-1/, V(56) /3.0753241996117268D-2/
      DATA U(57) /9.3727339240070590D-1/, V(57) /7.0366047488108125D-2/
      DATA U(58) /8.4820658341042722D-1/, V(58) /1.0715922046717194D-1/
      DATA U(59) /7.2441773136017005D-1/, V(59) /1.3957067792615431D-1/
      DATA U(60) /5.7097217260853885D-1/, V(60) /1.6626920581699393D-1/
      DATA U(61) /3.9415134707756337D-1/, V(61) /1.8616100001556221D-1/
      DATA U(62) /2.0119409399743452D-1/, V(62) /1.9843148532711158D-1/
      DATA U(63) /0/                    , V(63) /2.0257824192556127D-1/
C  N=16.
      DATA U(64) /9.8940093499164993D-1/, V(64) /2.7152459411754095D-2/
      DATA U(65) /9.4457502307323258D-1/, V(65) /6.2253523938647893D-2/
      DATA U(66) /8.6563120238783174D-1/, V(66) /9.5158511682492785D-2/
      DATA U(67) /7.5540440835500303D-1/, V(67) /1.2462897125553387D-1/
      DATA U(68) /6.1787624440264375D-1/, V(68) /1.4959598881657673D-1/
      DATA U(69) /4.5801677765722739D-1/, V(69) /1.6915651939500254D-1/
      DATA U(70) /2.8160355077925891D-1/, V(70) /1.8260341504492359D-1/
      DATA U(71) /9.5012509837637440D-2/, V(71) /1.8945061045506850D-1/
C  N=20.
      DATA U(72) /9.9312859918509492D-1/, V(72) /1.7614007139152118D-2/
      DATA U(73) /9.6397192727791379D-1/, V(73) /4.0601429800386941D-2/
      DATA U(74) /9.1223442825132591D-1/, V(74) /6.2672048334109064D-2/
      DATA U(75) /8.3911697182221882D-1/, V(75) /8.3276741576704749D-2/
      DATA U(76) /7.4633190646015079D-1/, V(76) /1.0193011981724044D-1/
      DATA U(77) /6.3605368072651503D-1/, V(77) /1.1819453196151842D-1/
      DATA U(78) /5.1086700195082710D-1/, V(78) /1.3168863844917663D-1/
      DATA U(79) /3.7370608871541956D-1/, V(79) /1.4209610931838205D-1/
      DATA U(80) /2.2778585114164508D-1/, V(80) /1.4917298647260374D-1/
      DATA U(81) /7.6526521133497334D-2/, V(81) /1.5275338713072585D-1/
C  N=24.
      DATA U(82) /9.9518721999702136D-1/, V(82) /1.2341229799987200D-2/
      DATA U(83) /9.7472855597130950D-1/, V(83) /2.8531388628933663D-2/
      DATA U(84) /9.3827455200273276D-1/, V(84) /4.4277438817419806D-2/
      DATA U(85) /8.8641552700440103D-1/, V(85) /5.9298584915436781D-2/
      DATA U(86) /8.2000198597390292D-1/, V(86) /7.3346481411080306D-2/
      DATA U(87) /7.4012419157855436D-1/, V(87) /8.6190161531953276D-2/
      DATA U(88) /6.4809365193697557D-1/, V(88) /9.7618652104113888D-2/
      DATA U(89) /5.4542147138883954D-1/, V(89) /1.0744427011596563D-1/
      DATA U(90) /4.3379350762604514D-1/, V(90) /1.1550566805372560D-1/
      DATA U(91) /3.1504267969616337D-1/, V(91) /1.2167047292780339D-1/
      DATA U(92) /1.9111886747361631D-1/, V(92) /1.2583745634682830D-1/
      DATA U(93) /6.4056892862605626D-2/, V(93) /1.2793819534675216D-1/
C  N=32.
      DATA U(94) /9.9726386184948156D-1/, V(94) /7.0186100094700966D-3/
      DATA U(95) /9.8561151154526834D-1/, V(95) /1.6274394730905671D-2/
      DATA U(96) /9.6476225558750643D-1/, V(96) /2.5392065309262059D-2/
      DATA U(97) /9.3490607593773969D-1/, V(97) /3.4273862913021433D-2/
      DATA U(98) /8.9632115576605212D-1/, V(98) /4.2835898022226681D-2/
      DATA U(99) /8.4936761373256997D-1/, V(99) /5.0998059262376176D-2/
      DATA U(100)/7.9448379596794241D-1/, V(100)/5.8684093478535547D-2/
      DATA U(101)/7.3218211874028968D-1/, V(101)/6.5822222776361847D-2/
      DATA U(102)/6.6304426693021520D-1/, V(102)/7.2345794108848506D-2/
      DATA U(103)/5.8771575724076233D-1/, V(103)/7.8193895787070306D-2/
      DATA U(104)/5.0689990893222939D-1/, V(104)/8.3311924226946755D-2/
      DATA U(105)/4.2135127613063535D-1/, V(105)/8.7652093004403811D-2/
      DATA U(106)/3.3186860228212765D-1/, V(106)/9.1173878695763885D-2/
      DATA U(107)/2.3928736225213707D-1/, V(107)/9.3844399080804566D-2/
      DATA U(108)/1.4447196158279649D-1/, V(108)/9.5638720079274859D-2/
      DATA U(109)/4.8307665687738316D-2/, V(109)/9.6540088514727801D-2/
C  N=40.
      DATA U(110)/9.9823770971055920D-1/, V(110)/4.5212770985331913D-3/
      DATA U(111)/9.9072623869945701D-1/, V(111)/1.0498284531152814D-2/
      DATA U(112)/9.7725994998377426D-1/, V(112)/1.6421058381907889D-2/
      DATA U(113)/9.5791681921379166D-1/, V(113)/2.2245849194166957D-2/
      DATA U(114)/9.3281280827867653D-1/, V(114)/2.7937006980023401D-2/
      DATA U(115)/9.0209880696887430D-1/, V(115)/3.3460195282547847D-2/
      DATA U(116)/8.6595950321225950D-1/, V(116)/3.8782167974472018D-2/
      DATA U(117)/8.2461223083331166D-1/, V(117)/4.3870908185673272D-2/
      DATA U(118)/7.7830565142651939D-1/, V(118)/4.8695807635072232D-2/
      DATA U(119)/7.2731825518992710D-1/, V(119)/5.3227846983936824D-2/
      DATA U(120)/6.7195668461417955D-1/, V(120)/5.7439769099391551D-2/
      DATA U(121)/6.1255388966798024D-1/, V(121)/6.1306242492928939D-2/
      DATA U(122)/5.4946712509512820D-1/, V(122)/6.4804013456601038D-2/
      DATA U(123)/4.8307580168617871D-1/, V(123)/6.7912045815233904D-2/
      DATA U(124)/4.1377920437160500D-1/, V(124)/7.0611647391286780D-2/
      DATA U(125)/3.4199409082575847D-1/, V(125)/7.2886582395804059D-2/
      DATA U(126)/2.6815218500725368D-1/, V(126)/7.4723169057968264D-2/
      DATA U(127)/1.9269758070137110D-1/, V(127)/7.6110361900626242D-2/
      DATA U(128)/1.1608407067525521D-1/, V(128)/7.7039818164247966D-2/
      DATA U(129)/3.8772417506050822D-2/, V(129)/7.7505947978424811D-2/
C  N=48.
      DATA U(130)/9.9877100725242612D-1/, V(130)/3.1533460523058386D-3/
      DATA U(131)/9.9353017226635076D-1/, V(131)/7.3275539012762621D-3/
      DATA U(132)/9.8412458372282686D-1/, V(132)/1.1477234579234539D-2/
      DATA U(133)/9.7059159254624725D-1/, V(133)/1.5579315722943849D-2/
      DATA U(134)/9.5298770316043086D-1/, V(134)/1.9616160457355528D-2/
      DATA U(135)/9.3138669070655433D-1/, V(135)/2.3570760839324379D-2/
      DATA U(136)/9.0587913671556967D-1/, V(136)/2.7426509708356948D-2/
      DATA U(137)/8.7657202027424789D-1/, V(137)/3.1167227832798089D-2/
      DATA U(138)/8.4358826162439353D-1/, V(138)/3.4777222564770439D-2/
      DATA U(139)/8.0706620402944263D-1/, V(139)/3.8241351065830706D-2/
      DATA U(140)/7.6715903251574034D-1/, V(140)/4.1545082943464749D-2/
      DATA U(141)/7.2403413092381465D-1/, V(141)/4.4674560856694280D-2/
      DATA U(142)/6.7787237963266391D-1/, V(142)/4.7616658492490475D-2/
      DATA U(143)/6.2886739677651362D-1/, V(143)/5.0359035553854475D-2/
      DATA U(144)/5.7722472608397270D-1/, V(144)/5.2890189485193667D-2/
      DATA U(145)/5.2316097472223303D-1/, V(145)/5.5199503699984163D-2/
      DATA U(146)/4.6690290475095840D-1/, V(146)/5.7277292100403216D-2/
      DATA U(147)/4.0868648199071673D-1/, V(147)/5.9114839698395636D-2/
      DATA U(148)/3.4875588629216074D-1/, V(148)/6.0704439165893880D-2/
      DATA U(149)/2.8736248735545558D-1/, V(149)/6.2039423159892664D-2/
      DATA U(150)/2.2476379039468906D-1/, V(150)/6.3114192286254026D-2/
      DATA U(151)/1.6122235606889172D-1/, V(151)/6.3924238584648187D-2/
      DATA U(152)/9.7004699209462699D-2/, V(152)/6.4466164435950082D-2/
      DATA U(153)/3.2380170962869362D-2/, V(153)/6.4737696812683923D-2/
C  N=64.
      DATA U(154)/9.9930504173577214D-1/, V(154)/1.7832807216964329D-3/
      DATA U(155)/9.9634011677195528D-1/, V(155)/4.1470332605624676D-3/
      DATA U(156)/9.9101337147674432D-1/, V(156)/6.5044579689783629D-3/
      DATA U(157)/9.8333625388462596D-1/, V(157)/8.8467598263639477D-3/
      DATA U(158)/9.7332682778991096D-1/, V(158)/1.1168139460131129D-2/
      DATA U(159)/9.6100879965205372D-1/, V(159)/1.3463047896718643D-2/
      DATA U(160)/9.4641137485840282D-1/, V(160)/1.5726030476024719D-2/
      DATA U(161)/9.2956917213193958D-1/, V(161)/1.7951715775697343D-2/
      DATA U(162)/9.1052213707850281D-1/, V(162)/2.0134823153530209D-2/
      DATA U(163)/8.8931544599511412D-1/, V(163)/2.2270173808383254D-2/
      DATA U(164)/8.6599939815409282D-1/, V(164)/2.4352702568710873D-2/
      DATA U(165)/8.4062929625258036D-1/, V(165)/2.6377469715054659D-2/
      DATA U(166)/8.1326531512279756D-1/, V(166)/2.8339672614259483D-2/
      DATA U(167)/7.8397235894334141D-1/, V(167)/3.0234657072402479D-2/
      DATA U(168)/7.5281990726053190D-1/, V(168)/3.2057928354851554D-2/
      DATA U(169)/7.1988185017161083D-1/, V(169)/3.3805161837141609D-2/
      DATA U(170)/6.8523631305423324D-1/, V(170)/3.5472213256882384D-2/
      DATA U(171)/6.4896547125465734D-1/, V(171)/3.7055128540240046D-2/
      DATA U(172)/6.1115535517239325D-1/, V(172)/3.8550153178615629D-2/
      DATA U(173)/5.7189564620263403D-1/, V(173)/3.9953741132720341D-2/
      DATA U(174)/5.3127946401989455D-1/, V(174)/4.1262563242623529D-2/
      DATA U(175)/4.8940314570705296D-1/, V(175)/4.2473515123653589D-2/
      DATA U(176)/4.4636601725346409D-1/, V(176)/4.3583724529323453D-2/
      DATA U(177)/4.0227015796399160D-1/, V(177)/4.4590558163756563D-2/
      DATA U(178)/3.5722015833766812D-1/, V(178)/4.5491627927418144D-2/
      DATA U(179)/3.1132287199021096D-1/, V(179)/4.6284796581314417D-2/
      DATA U(180)/2.6468716220876742D-1/, V(180)/4.6968182816210017D-2/
      DATA U(181)/2.1742364374000708D-1/, V(181)/4.7540165714830309D-2/
      DATA U(182)/1.6964442042399282D-1/, V(182)/4.7999388596458308D-2/
      DATA U(183)/1.2146281929612055D-1/, V(183)/4.8344762234802957D-2/
      DATA U(184)/7.2993121787799039D-2/, V(184)/4.8575467441503427D-2/
      DATA U(185)/2.4350292663424433D-2/, V(185)/4.8690957009139720D-2/
C  N=80.
      DATA U(186)/9.9955382265163063D-1/, V(186)/1.1449500031869415D-3/
      DATA U(187)/9.9764986439823769D-1/, V(187)/2.6635335895126817D-3/
      DATA U(188)/9.9422754096568828D-1/, V(188)/4.1803131246948952D-3/
      DATA U(189)/9.8929130249975553D-1/, V(189)/5.6909224514031986D-3/
      DATA U(190)/9.8284857273862907D-1/, V(190)/7.1929047681173128D-3/
      DATA U(191)/9.7490914058572779D-1/, V(191)/8.6839452692608584D-3/
      DATA U(192)/9.6548508904379925D-1/, V(192)/1.0161766041103065D-2/
      DATA U(193)/9.5459076634363491D-1/, V(193)/1.1624114120797827D-2/
      DATA U(194)/9.4224276130987267D-1/, V(194)/1.3068761592401339D-2/
      DATA U(195)/9.2845987717244580D-1/, V(195)/1.4493508040509076D-2/
      DATA U(196)/9.1326310257175765D-1/, V(196)/1.5896183583725688D-2/
      DATA U(197)/8.9667557943877068D-1/, V(197)/1.7274652056269306D-2/
      DATA U(198)/8.7872256767821383D-1/, V(198)/1.8626814208299031D-2/
      DATA U(199)/8.5943140666311110D-1/, V(199)/1.9950610878141999D-2/
      DATA U(200)/8.3883147358025528D-1/, V(200)/2.1244026115782006D-2/
      DATA U(201)/8.1695413868146347D-1/, V(201)/2.2505090246332462D-2/
      DATA U(202)/7.9383271750460545D-1/, V(202)/2.3731882865930101D-2/
      DATA U(203)/7.6950242013504137D-1/, V(203)/2.4922535764115491D-2/
      DATA U(204)/7.4400029758359727D-1/, V(204)/2.6075235767565118D-2/
      DATA U(205)/7.1736518536209988D-1/, V(205)/2.7188227500486381D-2/
      DATA U(206)/6.8963764434202760D-1/, V(206)/2.8259816057276862D-2/
      DATA U(207)/6.6085989898611980D-1/, V(207)/2.9288369583267848D-2/
      DATA U(208)/6.3107577304687197D-1/, V(208)/3.0272321759557981D-2/
      DATA U(209)/6.0033062282975174D-1/, V(209)/3.1210174188114702D-2/
      DATA U(210)/5.6867126812270978D-1/, V(210)/3.2100498673487773D-2/
      DATA U(211)/5.3614592089713193D-1/, V(211)/3.2941939397645401D-2/
      DATA U(212)/5.0280411188878499D-1/, V(212)/3.3733214984611523D-2/
      DATA U(213)/4.6869661517054448D-1/, V(213)/3.4473120451753929D-2/
      DATA U(214)/4.3387537083175609D-1/, V(214)/3.5160529044747593D-2/
      DATA U(215)/3.9839340588196923D-1/, V(215)/3.5794393953416055D-2/
      DATA U(216)/3.6230475349948732D-1/, V(216)/3.6373749905835978D-2/
      DATA U(217)/3.2566437074770191D-1/, V(217)/3.6897714638276009D-2/
      DATA U(218)/2.8852805488451185D-1/, V(218)/3.7365490238730490D-2/
      DATA U(219)/2.5095235839227212D-1/, V(219)/3.7776364362001397D-2/
      DATA U(220)/2.1299450285766613D-1/, V(220)/3.8129711314477638D-2/
      DATA U(221)/1.7471229183264681D-1/, V(221)/3.8424993006959423D-2/
      DATA U(222)/1.3616402280914389D-1/, V(222)/3.8661759774076463D-2/
      DATA U(223)/9.7408398441584599D-2/, V(223)/3.8839651059051969D-2/
      DATA U(224)/5.8504437152420669D-2/, V(224)/3.8958395962769531D-2/
      DATA U(225)/1.9511383256793998D-2/, V(225)/3.9017813656306655D-2/
C  N=96.
      DATA U(226)/9.9968950388323077D-1/, V(226)/7.9679206555201243D-4/
      DATA U(227)/9.9836437586318168D-1/, V(227)/1.8539607889469217D-3/
      DATA U(228)/9.9598184298720929D-1/, V(228)/2.9107318179349464D-3/
      DATA U(229)/9.9254390032376262D-1/, V(229)/3.9645543384446867D-3/
      DATA U(230)/9.8805412632962380D-1/, V(230)/5.0142027429275177D-3/
      DATA U(231)/9.8251726356301468D-1/, V(231)/6.0585455042359617D-3/
      DATA U(232)/9.7593917458513647D-1/, V(232)/7.0964707911538653D-3/
      DATA U(233)/9.6832682846326421D-1/, V(233)/8.1268769256987592D-3/
      DATA U(234)/9.5968829144874254D-1/, V(234)/9.1486712307833866D-3/
      DATA U(235)/9.5003271778443764D-1/, V(235)/1.0160770535008416D-2/
      DATA U(236)/9.3937033975275522D-1/, V(236)/1.1162102099838499D-2/
      DATA U(237)/9.2771245672230869D-1/, V(237)/1.2151604671088320D-2/
      DATA U(238)/9.1507142312089807D-1/, V(238)/1.3128229566961573D-2/
      DATA U(239)/9.0146063531585234D-1/, V(239)/1.4090941772314861D-2/
      DATA U(240)/8.8689451740242042D-1/, V(240)/1.5038721026994938D-2/
      DATA U(241)/8.7138850590929650D-1/, V(241)/1.5970562902562291D-2/
      DATA U(242)/8.5495903343460146D-1/, V(242)/1.6885479864245172D-2/
      DATA U(243)/8.3762351122818712D-1/, V(243)/1.7782502316045261D-2/
      DATA U(244)/8.1940031073793168D-1/, V(244)/1.8660679627411467D-2/
      DATA U(245)/8.0030874413914082D-1/, V(245)/1.9519081140145022D-2/
      DATA U(246)/7.8036904386743322D-1/, V(246)/2.0356797154333325D-2/
      DATA U(247)/7.5960234117664750D-1/, V(247)/2.1172939892191299D-2/
      DATA U(248)/7.3803064374440013D-1/, V(248)/2.1966644438744349D-2/
      DATA U(249)/7.1567681234896763D-1/, V(249)/2.2737069658329374D-2/
      DATA U(250)/6.9256453664217156D-1/, V(250)/2.3483399085926220D-2/
      DATA U(251)/6.6871831004391615D-1/, V(251)/2.4204841792364691D-2/
      DATA U(252)/6.4416340378496712D-1/, V(252)/2.4900633222483610D-2/
      DATA U(253)/6.1892584012546857D-1/, V(253)/2.5570036005349361D-2/
      DATA U(254)/5.9303236477757208D-1/, V(254)/2.6212340735672414D-2/
      DATA U(255)/5.6651041856139717D-1/, V(255)/2.6826866725591762D-2/
      DATA U(256)/5.3938810832435744D-1/, V(256)/2.7412962726029243D-2/
      DATA U(257)/5.1169417715466767D-1/, V(257)/2.7970007616848334D-2/
      DATA U(258)/4.8345797392059636D-1/, V(258)/2.8497411065085386D-2/
      DATA U(259)/4.5470942216774301D-1/, V(259)/2.8994614150555237D-2/
      DATA U(260)/4.2547898840730055D-1/, V(260)/2.9461089958167906D-2/
      DATA U(261)/3.9579764982890860D-1/, V(261)/2.9896344136328386D-2/
      DATA U(262)/3.6569686147231364D-1/, V(262)/3.0299915420827594D-2/
      DATA U(263)/3.3520852289262542D-1/, V(263)/3.0671376123669149D-2/
      DATA U(264)/3.0436494435449635D-1/, V(264)/3.1010332586313837D-2/
      DATA U(265)/2.7319881259104914D-1/, V(265)/3.1316425596861356D-2/
      DATA U(266)/2.4174315616384001D-1/, V(266)/3.1589330770727167D-2/
      DATA U(267)/2.1003131046056720D-1/, V(267)/3.1828758894411006D-2/
      DATA U(268)/1.7809688236761860D-1/, V(268)/3.2034456231992663D-2/
      DATA U(269)/1.4597371465489694D-1/, V(269)/3.2206204794030251D-2/
      DATA U(270)/1.1369585011066592D-1/, V(270)/3.2343822568575928D-2/
      DATA U(271)/8.1297495464425559D-2/, V(271)/3.2447163714064269D-2/
      DATA U(272)/4.8812985136049731D-2/, V(272)/3.2516118713868836D-2/
      DATA U(273)/1.6276744849602970D-2/, V(273)/3.2550614492363166D-2/

      IF(KTBA(MIN(MAX(1,N),97)) .EQ. 0) THEN
       X(1)=0
       WRITE(ERRTXT,101) N
       CALL MTLPRT(NAME(MODE),'D107.1',ERRTXT)
       RETURN
      ENDIF
      ALFA=HF*(B+A)
      BETA=HF*(B-A)
      IF(MODE .EQ. 1) THEN
       SUM=0
       J1=MOD(N,2)
       J2=KTBA(N)+(N-1)/2
       DO 1 J = KTBA(N),J2-J1
       DELTA=BETA*U(J)
       SUM=SUM+V(J)*(F(ALFA+DELTA)+F(ALFA-DELTA))
    1  CONTINUE
       IF(J1 .EQ. 1) SUM=SUM+V(J2)*F(ALFA)
       X(1)=BETA*SUM
      ELSE
       J1=KTBA(N)-1
       J2=N+1
       DO 2 J=1,J2/2
       WTEMP=BETA*V(J1+J)
       DELTA=BETA*U(J1+J)
       X(J)=ALFA-DELTA
       W(J)=WTEMP
       X(J2-J)=ALFA+DELTA
       W(J2-J)=WTEMP
    2  CONTINUE
      ENDIF
      RETURN
  101 FORMAT('N = ',I5,' IS NON-PERMISSIBLE')
      END


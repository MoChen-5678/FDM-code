!*****************************************************************************!
!
MODULE DBESIK
!
!*****************************************************************************!
!
USE ERRORS
USE DImach
USE LNgamma
USE DBSK01
implicit none
!--- Contains DASYIK, DBSKNU, DBESI, DBESK

Contains
!--- DECK DASYIK
!
!=============================================================================!
!
SUBROUTINE DASYIK (X, FNU, KODE, FLGIK, RA, ARG, IN, Y)
!
!=============================================================================!
!
!***BEGIN PROLOGUE  DASYIK
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DBESI and DBESK
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (ASYIK-S, DASYIK-D)
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!                    DASYIK computes Bessel functions I and K
!                  for arguments X.GT.0.0 and orders FNU.GE.35
!                  on FLGIK = 1 and FLGIK = -1 respectively.
!
!                                    INPUT
!
!      X    - Argument, X.GT.0.0D0
!      FNU  - Order of first Bessel function
!      KODE - A parameter to indicate the scaling option
!             KODE=1 returns Y(I)=        I/SUB(FNU+I-1)/(X), I=1,IN
!                    or      Y(I)=        K/SUB(FNU+I-1)/(X), I=1,IN
!                    on FLGIK = 1.0D0 or FLGIK = -1.0D0
!             KODE=2 returns Y(I)=EXP(-X)*I/SUB(FNU+I-1)/(X), I=1,IN
!                    or      Y(I)=EXP( X)*K/SUB(FNU+I-1)/(X), I=1,IN
!                    on FLGIK = 1.0D0 or FLGIK = -1.0D0
!     FLGIK - Selection parameter for I or K FUNCTION
!             FLGIK =  1.0D0 gives the I function
!             FLGIK = -1.0D0 gives the K function
!        RA - SQRT(1.+Z*Z), Z=X/FNU
!       ARG - Argument of the leading exponential
!        IN - Number of functions desired, IN=1 or 2
!
!                                    OUTPUT
!
!         Y - A vector whose first IN components contain the sequence
!
!     Abstract  **** A double precision routine ****
!         DASYIK implements the uniform asymptotic expansion of
!         the I and K Bessel functions for FNU.GE.35 and real
!         X.GT.0.0D0. The forms are identical except for a change
!         in sign of some of the terms. This change in sign is
!         accomplished by means of the FLAG FLGIK = 1 or -1.
!
!***SEE ALSO  DBESI, DBESK
!***ROUTINES CALLED  D1MACH
!***REVISION HISTORY  (YYMMDD)
!   750101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910408  Updated the AUTHOR section.  (WRB)
!***END PROLOGUE  DASYIK
!
      INTEGER IN, J, JN, K, KK, KODE, L
      DOUBLE PRECISION AK,AP,ARG,C,COEF,CON,ETX,FLGIK,FN,FNU,GLN,RA
      DOUBLE PRECISION S1, S2, T, TOL, T2, X, Y, Z

      DIMENSION Y(*), C(65), CON(2)
      SAVE CON, C
      DATA CON(1), CON(2)  /3.98942280401432678D-01,    1.25331413731550025D+00/
      DATA C(1), C(2), C(3), C(4), C(5), C(6), C(7), C(8), C(9), C(10), &
     &     C(11), C(12), C(13), C(14), C(15), C(16), C(17), C(18), &
     &     C(19), C(20), C(21), C(22), C(23), C(24)/ &
     &       -2.08333333333333D-01,        1.25000000000000D-01, &
     &        3.34201388888889D-01,       -4.01041666666667D-01, &
     &        7.03125000000000D-02,       -1.02581259645062D+00, &
     &        1.84646267361111D+00,       -8.91210937500000D-01, &
     &        7.32421875000000D-02,        4.66958442342625D+00, &
     &       -1.12070026162230D+01,        8.78912353515625D+00, &
     &       -2.36408691406250D+00,        1.12152099609375D-01, &
     &       -2.82120725582002D+01,        8.46362176746007D+01, &
     &       -9.18182415432400D+01,        4.25349987453885D+01, &
     &       -7.36879435947963D+00,        2.27108001708984D-01, &
     &        2.12570130039217D+02,       -7.65252468141182D+02, &
     &        1.05999045252800D+03,       -6.99579627376133D+02/
      DATA C(25), C(26), C(27), C(28), C(29), C(30), C(31), C(32), &
     &     C(33), C(34), C(35), C(36), C(37), C(38), C(39), C(40), &
     &     C(41), C(42), C(43), C(44), C(45), C(46), C(47), C(48)/ &
     &        2.18190511744212D+02,       -2.64914304869516D+01, &
     &        5.72501420974731D-01,       -1.91945766231841D+03, &
     &        8.06172218173731D+03,       -1.35865500064341D+04, &
     &        1.16553933368645D+04,       -5.30564697861340D+03, &
     &        1.20090291321635D+03,       -1.08090919788395D+02, &
     &        1.72772750258446D+00,        2.02042913309661D+04, &
     &       -9.69805983886375D+04,        1.92547001232532D+05, &
     &       -2.03400177280416D+05,        1.22200464983017D+05, &
     &       -4.11926549688976D+04,        7.10951430248936D+03, &
     &       -4.93915304773088D+02,        6.07404200127348D+00, &
     &       -2.42919187900551D+05,        1.31176361466298D+06, &
     &       -2.99801591853811D+06,        3.76327129765640D+06/
      DATA C(49), C(50), C(51), C(52), C(53), C(54), C(55), C(56), &
     &     C(57), C(58), C(59), C(60), C(61), C(62), C(63), C(64), &
     &     C(65)/ &
     &       -2.81356322658653D+06,        1.26836527332162D+06, &
     &       -3.31645172484564D+05,        4.52187689813627D+04, &
     &       -2.49983048181121D+03,        2.43805296995561D+01, &
     &        3.28446985307204D+06,       -1.97068191184322D+07, &
     &        5.09526024926646D+07,       -7.41051482115327D+07, &
     &        6.63445122747290D+07,       -3.75671766607634D+07, &
     &        1.32887671664218D+07,       -2.78561812808645D+06, &
     &        3.08186404612662D+05,       -1.38860897537170D+04, &
     &        1.10017140269247D+02/
!***FIRST EXECUTABLE STATEMENT  DASYIK
      TOL = D1MACH(3)
      TOL = MAX(TOL,1.0D-15)
      FN = FNU
      Z  = (3.0D0-FLGIK)/2.0D0
      KK = INT(Z)
      DO 50 JN=1,IN
        IF (JN.EQ.1) GO TO 10
        FN = FN - FLGIK
        Z = X/FN
        RA = SQRT(1.0D0+Z*Z)
        GLN = LOG((1.0D0+RA)/Z)
        ETX = KODE - 1
        T = RA*(1.0D0-ETX) + ETX/(Z+RA)
        ARG = FN*(T-GLN)*FLGIK
   10   COEF = EXP(ARG)
        T = 1.0D0/RA
        T2 = T*T
        T = T/FN
        T = SIGN(T,FLGIK)
        S2 = 1.0D0
        AP = 1.0D0
        L = 0
        DO 30 K=2,11
          L = L + 1
          S1 = C(L)
          DO 20 J=2,K
            L = L + 1
            S1 = S1*T2 + C(L)
   20     CONTINUE
          AP = AP*T
          AK = AP*S1
          S2 = S2 + AK
          IF (MAX(ABS(AK),ABS(AP)) .LT.TOL) GO TO 40
   30   CONTINUE
   40   CONTINUE
      T = ABS(T)
      Y(JN) = S2*COEF*SQRT(T)*CON(KK)
   50 CONTINUE
      RETURN
END SUBROUTINE DASYIK
!--- DECK DBSKNU
!
!=============================================================================!
!
SUBROUTINE DBSKNU (X, FNU, KODE, N, Y, NZ)
!
!=============================================================================!
!
!***BEGIN PROLOGUE  DBSKNU
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DBESK
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (BESKNU-S, DBSKNU-D)
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Abstract  **** A DOUBLE PRECISION routine ****
!         DBSKNU computes N member sequences of K Bessel functions
!         K/SUB(FNU+I-1)/(X), I=1,N for non-negative orders FNU and
!         positive X. Equations of the references are implemented on
!         small orders DNU for K/SUB(DNU)/(X) and K/SUB(DNU+1)/(X).
!         Forward recursion with the three term recursion relation
!         generates higher orders FNU+I-1, I=1,...,N. The parameter
!         KODE permits K/SUB(FNU+I-1)/(X) values or scaled values
!         EXP(X)*K/SUB(FNU+I-1)/(X), I=1,N to be returned.
!
!         To start the recursion FNU is normalized to the interval
!         -0.5.LE.DNU.LT.0.5. A special form of the power series is
!         implemented on 0.LT.X.LE.X1 while the Miller algorithm for the
!         K Bessel function in terms of the confluent hypergeometric
!         function U(FNU+0.5,2*FNU+1,X) is implemented on X1.LT.X.LE.X2.
!         For X.GT.X2, the asymptotic expansion for large X is used.
!         When FNU is a half odd integer, a special formula for
!         DNU=-0.5 and DNU+1.0=0.5 is used to start the recursion.
!
!         The maximum number of significant digits obtainable
!         is the smaller of 14 and the number of digits carried in
!         DOUBLE PRECISION arithmetic.
!
!         DBSKNU assumes that a significant digit SINH function is
!         available.
!
!     Description of Arguments
!
!         INPUT      X,FNU are DOUBLE PRECISION
!           X      - X.GT.0.0D0
!           FNU    - Order of initial K function, FNU.GE.0.0D0
!           N      - Number of members of the sequence, N.GE.1
!           KODE   - A parameter to indicate the scaling option
!                    KODE= 1  returns
!                             Y(I)=       K/SUB(FNU+I-1)/(X)
!                                  I=1,...,N
!                        = 2  returns
!                             Y(I)=EXP(X)*K/SUB(FNU+I-1)/(X)
!                                  I=1,...,N
!
!         OUTPUT     Y is DOUBLE PRECISION
!           Y      - A vector whose first N components contain values
!                    for the sequence
!                    Y(I)=       K/SUB(FNU+I-1)/(X), I=1,...,N or
!                    Y(I)=EXP(X)*K/SUB(FNU+I-1)/(X), I=1,...,N
!                    depending on KODE
!           NZ     - Number of components set to zero due to
!                    underflow,
!                    NZ= 0   , normal return
!                    NZ.NE.0 , first NZ components of Y set to zero
!                              due to underflow, Y(I)=0.0D0,I=1,...,NZ
!
!     Error Conditions
!         Improper input arguments - a fatal error
!         Overflow - a fatal error
!         Underflow with KODE=1 - a non-fatal error (NZ.NE.0)
!
!***SEE ALSO  DBESK
!***REFERENCES  N. M. Temme, On the numerical evaluation of the modified
!                 Bessel function of the third kind, Journal of
!                 Computational Physics 19, (1975), pp. 324-337.
!***ROUTINES CALLED  D1MACH, DGAMMA, I1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790201  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!   900727  Added EXTERNAL statement.  (WRB)
!   910408  Updated the AUTHOR and REFERENCES sections.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DBSKNU
!
      INTEGER I, IFLAG, INU, J, K, KK, KODE, KODED, N, NN, NZ
      DOUBLE PRECISION A,AK,A1,A2,B,BK,CC,CK,COEF,CX,DK,DNU,DNU2,ELIM
      DOUBLE PRECISION ETEST, EX, F, FC, FHS, FK, FKS, FLRX, FMU, FNU, G1, G2, P, PI
      DOUBLE PRECISION PT, P1, P2, Q, RTHPI, RX, S, SMU, SQK, ST, S1, S2, TM, TOL, T1
      DOUBLE PRECISION T2, X, X1, X2, Y
      DIMENSION A(160), B(160), Y(*), CC(8)
      SAVE X1, X2, PI, RTHPI, CC
      DATA X1, X2 / 2.0D0, 17.0D0 /
      DATA PI,RTHPI        / 3.14159265358979D+00, 1.25331413731550D+00/
      DATA CC(1), CC(2), CC(3), CC(4), CC(5), CC(6), CC(7), CC(8)&
     &                     / 5.77215664901533D-01,-4.20026350340952D-02,&
     &-4.21977345555443D-02, 7.21894324666300D-03,-2.15241674114900D-04,&
     &-2.01348547807000D-05, 1.13302723200000D-06, 6.11609500000000D-09/
!***FIRST EXECUTABLE STATEMENT  DBSKNU
      KK = -I1MACH(15)
      ELIM = 2.303D0*(KK*D1MACH(5)-3.0D0)
      AK = D1MACH(3)
      TOL = MAX(AK,1.0D-15)
      IF (X.LE.0.0D0) GO TO 350
      IF (FNU.LT.0.0D0) GO TO 360
      IF (KODE.LT.1 .OR. KODE.GT.2) GO TO 370
      IF (N.LT.1) GO TO 380
      NZ = 0
      IFLAG = 0
      KODED = KODE
      RX = 2.0D0/X
      INU = INT(FNU+0.5D0)
      DNU = FNU - INU
      IF (ABS(DNU).EQ.0.5D0) GO TO 120
      DNU2 = 0.0D0
      IF (ABS(DNU).LT.TOL) GO TO 10
      DNU2 = DNU*DNU
   10 CONTINUE
      IF (X.GT.X1) GO TO 120
!
!     SERIES FOR X.LE.X1
!
      A1 = 1.0D0 - DNU
      A2 = 1.0D0 + DNU
      T1 = 1.0D0/DGAMMA(A1)
      T2 = 1.0D0/DGAMMA(A2)
      IF (ABS(DNU).GT.0.1D0) GO TO 40
!     SERIES FOR F0 TO RESOLVE INDETERMINACY FOR SMALL ABS(DNU)
      S = CC(1)
      AK = 1.0D0
      DO 20 K=2,8
        AK = AK*DNU2
        TM = CC(K)*AK
        S = S + TM
        IF (ABS(TM).LT.TOL) GO TO 30
   20 CONTINUE
   30 G1 = -S
      GO TO 50
   40 CONTINUE
      G1 = (T1-T2)/(DNU+DNU)
   50 CONTINUE
      G2 = (T1+T2)*0.5D0
      SMU = 1.0D0
      FC = 1.0D0
      FLRX = LOG(RX)
      FMU = DNU*FLRX
      IF (DNU.EQ.0.0D0) GO TO 60
      FC = DNU*PI
      FC = FC/SIN(FC)
      IF (FMU.NE.0.0D0) SMU = SINH(FMU)/FMU
   60 CONTINUE
      F = FC*(G1*COSH(FMU)+G2*FLRX*SMU)
      FC = EXP(FMU)
      P = 0.5D0*FC/T2
      Q = 0.5D0/(FC*T1)
      AK = 1.0D0
      CK = 1.0D0
      BK = 1.0D0
      S1 = F
      S2 = P
      IF (INU.GT.0 .OR. N.GT.1) GO TO 90
      IF (X.LT.TOL) GO TO 80
      CX = X*X*0.25D0
   70 CONTINUE
      F = (AK*F+P+Q)/(BK-DNU2)
      P = P/(AK-DNU)
      Q = Q/(AK+DNU)
      CK = CK*CX/AK
      T1 = CK*F
      S1 = S1 + T1
      BK = BK + AK + AK + 1.0D0
      AK = AK + 1.0D0
      S = ABS(T1)/(1.0D0+ABS(S1))
      IF (S.GT.TOL) GO TO 70
   80 CONTINUE
      Y(1) = S1
      IF (KODED.EQ.1) RETURN
      Y(1) = S1*EXP(X)
      RETURN
   90 CONTINUE
      IF (X.LT.TOL) GO TO 110
      CX = X*X*0.25D0
  100 CONTINUE
      F = (AK*F+P+Q)/(BK-DNU2)
      P = P/(AK-DNU)
      Q = Q/(AK+DNU)
      CK = CK*CX/AK
      T1 = CK*F
      S1 = S1 + T1
      T2 = CK*(P-AK*F)
      S2 = S2 + T2
      BK = BK + AK + AK + 1.0D0
      AK = AK + 1.0D0
      S = ABS(T1)/(1.0D0+ABS(S1)) + ABS(T2)/(1.0D0+ABS(S2))
      IF (S.GT.TOL) GO TO 100
  110 CONTINUE
      S2 = S2*RX
      IF (KODED.EQ.1) GO TO 170
      F = EXP(X)
      S1 = S1*F
      S2 = S2*F
      GO TO 170
  120 CONTINUE
      COEF = RTHPI/SQRT(X)
      IF (KODED.EQ.2) GO TO 130
      IF (X.GT.ELIM) GO TO 330
      COEF = COEF*EXP(-X)
  130 CONTINUE
      IF (ABS(DNU).EQ.0.5D0) GO TO 340
      IF (X.GT.X2) GO TO 280
!
!     MILLER ALGORITHM FOR X1.LT.X.LE.X2
!
      ETEST = COS(PI*DNU)/(PI*X*TOL)
      FKS = 1.0D0
      FHS = 0.25D0
      FK = 0.0D0
      CK = X + X + 2.0D0
      P1 = 0.0D0
      P2 = 1.0D0
      K = 0
  140 CONTINUE
      K = K + 1
      FK = FK + 1.0D0
      AK = (FHS-DNU2)/(FKS+FK)
      BK = CK/(FK+1.0D0)
      PT = P2
      P2 = BK*P2 - AK*P1
      P1 = PT
      A(K) = AK
      B(K) = BK
      CK = CK + 2.0D0
      FKS = FKS + FK + FK + 1.0D0
      FHS = FHS + FK + FK
      IF (ETEST.GT.FK*P1) GO TO 140
      KK = K
      S = 1.0D0
      P1 = 0.0D0
      P2 = 1.0D0
      DO 150 I=1,K
        PT = P2
        P2 = (B(KK)*P2-P1)/A(KK)
        P1 = PT
        S = S + P2
        KK = KK - 1
  150 CONTINUE
      S1 = COEF*(P2/S)
      IF (INU.GT.0 .OR. N.GT.1) GO TO 160
      GO TO 200
  160 CONTINUE
      S2 = S1*(X+DNU+0.5D0-P1/P2)/X
!
!     FORWARD RECURSION ON THE THREE TERM RECURSION RELATION
!
  170 CONTINUE
      CK = (DNU+DNU+2.0D0)/X
      IF (N.EQ.1) INU = INU - 1
      IF (INU.GT.0) GO TO 180
      IF (N.GT.1) GO TO 200
      S1 = S2
      GO TO 200
  180 CONTINUE
      DO 190 I=1,INU
        ST = S2
        S2 = CK*S2 + S1
        S1 = ST
        CK = CK + RX
  190 CONTINUE
      IF (N.EQ.1) S1 = S2
  200 CONTINUE
      IF (IFLAG.EQ.1) GO TO 220
      Y(1) = S1
      IF (N.EQ.1) RETURN
      Y(2) = S2
      IF (N.EQ.2) RETURN
      DO 210 I=3,N
        Y(I) = CK*Y(I-1) + Y(I-2)
        CK = CK + RX
  210 CONTINUE
      RETURN
!     IFLAG=1 CASES
  220 CONTINUE
      S = -X + LOG(S1)
      Y(1) = 0.0D0
      NZ = 1
      IF (S.LT.-ELIM) GO TO 230
      Y(1) = EXP(S)
      NZ = 0
  230 CONTINUE
      IF (N.EQ.1) RETURN
      S = -X + LOG(S2)
      Y(2) = 0.0D0
      NZ = NZ + 1
      IF (S.LT.-ELIM) GO TO 240
      NZ = NZ - 1
      Y(2) = EXP(S)
  240 CONTINUE
      IF (N.EQ.2) RETURN
      KK = 2
      IF (NZ.LT.2) GO TO 260
      DO 250 I=3,N
        KK = I
        ST = S2
        S2 = CK*S2 + S1
        S1 = ST
        CK = CK + RX
        S = -X + LOG(S2)
        NZ = NZ + 1
        Y(I) = 0.0D0
        IF (S.LT.-ELIM) GO TO 250
        Y(I) = EXP(S)
        NZ = NZ - 1
        GO TO 260
  250 CONTINUE
      RETURN
  260 CONTINUE
      IF (KK.EQ.N) RETURN
      S2 = S2*CK + S1
      CK = CK + RX
      KK = KK + 1
      Y(KK) = EXP(-X+LOG(S2))
      IF (KK.EQ.N) RETURN
      KK = KK + 1
      DO 270 I=KK,N
        Y(I) = CK*Y(I-1) + Y(I-2)
        CK = CK + RX
  270 CONTINUE
      RETURN
!
!     ASYMPTOTIC EXPANSION FOR LARGE X, X.GT.X2
!
!     IFLAG=0 MEANS NO UNDERFLOW OCCURRED
!     IFLAG=1 MEANS AN UNDERFLOW OCCURRED- COMPUTATION PROCEEDS WITH
!     KODED=2 AND A TEST FOR ON SCALE VALUES IS MADE DURING FORWARD
!     RECURSION
  280 CONTINUE
      NN = 2
      IF (INU.EQ.0 .AND. N.EQ.1) NN = 1
      DNU2 = DNU + DNU
      FMU = 0.0D0
      IF (ABS(DNU2).LT.TOL) GO TO 290
      FMU = DNU2*DNU2
  290 CONTINUE
      EX = X*8.0D0
      S2 = 0.0D0
      DO 320 K=1,NN
        S1 = S2
        S = 1.0D0
        AK = 0.0D0
        CK = 1.0D0
        SQK = 1.0D0
        DK = EX
        DO 300 J=1,30
          CK = CK*(FMU-SQK)/DK
          S = S + CK
          DK = DK + EX
          AK = AK + 8.0D0
          SQK = SQK + AK
          IF (ABS(CK).LT.TOL) GO TO 310
  300   CONTINUE
  310   S2 = S*COEF
        FMU = FMU + 8.0D0*DNU + 4.0D0
  320 CONTINUE
      IF (NN.GT.1) GO TO 170
      S1 = S2
      GO TO 200
  330 CONTINUE
      KODED = 2
      IFLAG = 1
      GO TO 120
!
!     FNU=HALF ODD INTEGER CASE
!
  340 CONTINUE
      S1 = COEF
      S2 = COEF
      GO TO 170
!
!
  350 CALL XERMSG ('SLATEC', 'DBSKNU', 'X NOT GREATER THAN ZERO', 2, 1)
      RETURN
  360 CALL XERMSG ('SLATEC', 'DBSKNU', 'FNU NOT ZERO OR POSITIVE', 2, 1)
      RETURN
  370 CALL XERMSG ('SLATEC', 'DBSKNU', 'KODE NOT 1 OR 2', 2, 1)
      RETURN
  380 CALL XERMSG ('SLATEC', 'DBSKNU', 'N NOT GREATER THAN 0', 2, 1)
      RETURN
END SUBROUTINE DBSKNU
!--- DECK DBESI
!
!=============================================================================!
!
SUBROUTINE DBESI (X, ALPHA, KODE, N, Y, NZ)
!
!=============================================================================!
!
!***BEGIN PROLOGUE  DBESI
!***PURPOSE  Compute an N member sequence of I Bessel functions
!            I/SUB(ALPHA+K-1)/(X), K=1,...,N or scaled Bessel functions
!            EXP(-X)*I/SUB(ALPHA+K-1)/(X), K=1,...,N for nonnegative
!            ALPHA and X.
!***LIBRARY   SLATEC
!***CATEGORY  C10B3
!***TYPE      DOUBLE PRECISION (BESI-S, DBESI-D)
!***KEYWORDS  I BESSEL FUNCTION, SPECIAL FUNCTIONS
!***AUTHOR  Amos, D. E., (SNLA)
!           Daniel, S. L., (SNLA)
!***DESCRIPTION
!
!     Abstract  **** a double precision routine ****
!         DBESI computes an N member sequence of I Bessel functions
!         I/sub(ALPHA+K-1)/(X), K=1,...,N or scaled Bessel functions
!         EXP(-X)*I/sub(ALPHA+K-1)/(X), K=1,...,N for nonnegative ALPHA
!         and X.  A combination of the power series, the asymptotic
!         expansion for X to infinity, and the uniform asymptotic
!         expansion for NU to infinity are applied over subdivisions of
!         the (NU,X) plane.  For values not covered by one of these
!         formulae, the order is incremented by an integer so that one
!         of these formulae apply.  Backward recursion is used to reduce
!         orders by integer values.  The asymptotic expansion for X to
!         infinity is used only when the entire sequence (specifically
!         the last member) lies within the region covered by the
!         expansion.  Leading terms of these expansions are used to test
!         for over or underflow where appropriate.  If a sequence is
!         requested and the last member would underflow, the result is
!         set to zero and the next lower order tried, etc., until a
!         member comes on scale or all are set to zero.  An overflow
!         cannot occur with scaling.
!
!         The maximum number of significant digits obtainable
!         is the smaller of 14 and the number of digits carried in
!         double precision arithmetic.
!
!     Description of Arguments
!
!         Input      X,ALPHA are double precision
!           X      - X .GE. 0.0D0
!           ALPHA  - order of first member of the sequence,
!                    ALPHA .GE. 0.0D0
!           KODE   - a parameter to indicate the scaling option
!                    KODE=1 returns
!                           Y(K)=        I/sub(ALPHA+K-1)/(X),
!                                K=1,...,N
!                    KODE=2 returns
!                           Y(K)=EXP(-X)*I/sub(ALPHA+K-1)/(X),
!                                K=1,...,N
!           N      - number of members in the sequence, N .GE. 1
!
!         Output     Y is double precision
!           Y      - a vector whose first N components contain
!                    values for I/sub(ALPHA+K-1)/(X) or scaled
!                    values for EXP(-X)*I/sub(ALPHA+K-1)/(X),
!                    K=1,...,N depending on KODE
!           NZ     - number of components of Y set to zero due to
!                    underflow,
!                    NZ=0   , normal return, computation completed
!                    NZ .NE. 0, last NZ components of Y set to zero,
!                             Y(K)=0.0D0, K=N-NZ+1,...,N.
!
!     Error Conditions
!         Improper input arguments - a fatal error
!         Overflow with KODE=1 - a fatal error
!         Underflow - a non-fatal error(NZ .NE. 0)
!
!***REFERENCES  D. E. Amos, S. L. Daniel and M. K. Weston, CDC 6600
!                 subroutines IBESS and JBESS for Bessel functions
!                 I(NU,X) and J(NU,X), X .GE. 0, NU .GE. 0, ACM
!                 Transactions on Mathematical Software 3, (1977),
!                 pp. 76-92.
!               F. W. J. Olver, Tables of Bessel Functions of Moderate
!                 or Large Orders, NPL Mathematical Tables 6, Her
!                 Majesty's Stationery Office, London, 1962.
!***ROUTINES CALLED  D1MACH, DASYIK, DLNGAM, I1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   750101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DBESI
!
      INTEGER I, IALP, IN, INLIM, IS, I1, K, KK, KM, KODE, KT, N, NN, NS, NZ 
      DOUBLE PRECISION AIN,AK,AKM,ALPHA,ANS,AP,ARG,ATOL,TOLLN,DFN
      DOUBLE PRECISION DTM, DX, EARG, ELIM, ETX, FLGIK,FN, FNF, FNI,FNP1,FNU,GLN,RA
      DOUBLE PRECISION RTTPI, S, SX, SXO2, S1, S2, T, TA, TB, TEMP, TFN, TM, TOL
      DOUBLE PRECISION TRX, T2, X, XO2, XO2L, Y, Z

      DIMENSION Y(*), TEMP(3)
      SAVE RTTPI, INLIM
      DATA RTTPI           / 3.98942280401433D-01/
      DATA INLIM           /          80         /
!***FIRST EXECUTABLE STATEMENT  DBESI
      NZ = 0
      KT = 1
!     I1MACH(15) REPLACES I1MACH(12) IN A DOUBLE PRECISION CODE
!     I1MACH(14) REPLACES I1MACH(11) IN A DOUBLE PRECISION CODE
      RA = D1MACH(3)
      TOL = MAX(RA,1.0D-15)
      I1 = -I1MACH(15)
      GLN = D1MACH(5)
      ELIM = 2.303D0*(I1*GLN-3.0D0)
!     TOLLN = -LN(TOL)
      I1 = I1MACH(14)+1
      TOLLN = 2.303D0*GLN*I1
      TOLLN = MIN(TOLLN,34.5388D0)
!      IF (N-1) 590, 10, 20

      IF (N-1.lt.0) then
            go to  590
      else if(N-1.eq.0) then
            go to 10
      else
            go to 20
      end if

   10 KT = 2
   20 NN = N
      IF (KODE.LT.1 .OR. KODE.GT.2) GO TO 570
!      IF (X) 600, 30, 80

      IF (X.lt.0.0d0) then
            go to 600
      else if (X.eq.0.0d0) then
            go to 30
      else
            go to 80
      end if

!   30 IF (ALPHA) 580, 40, 50

   30 IF (ALPHA.lt.0.0d0) then
            go to 580
      else if (ALPHA.eq. 0.0d0) then
            go to 40
      else
            go to 50
      end if

   40 Y(1) = 1.0D0
      IF (N.EQ.1) RETURN
      I1 = 2
      GO TO 60
   50 I1 = 1
   60 DO 70 I=I1,N
        Y(I) = 0.0D0
   70 CONTINUE
      RETURN
   80 CONTINUE
      IF (ALPHA.LT.0.0D0) GO TO 580
!
      IALP = INT(ALPHA)
      FNI = IALP + N - 1
      FNF = ALPHA - IALP
      DFN = FNI + FNF
      FNU = DFN
      IN = 0
      XO2 = X*0.5D0
      SXO2 = XO2*XO2
      ETX = KODE - 1
      SX = ETX*X
!
!     DECISION TREE FOR REGION WHERE SERIES, ASYMPTOTIC EXPANSION FOR X
!     TO INFINITY AND ASYMPTOTIC EXPANSION FOR NU TO INFINITY ARE
!     APPLIED.
!
      IF (SXO2.LE.(FNU+1.0D0)) GO TO 90
      IF (X.LE.12.0D0) GO TO 110
      FN = 0.55D0*FNU*FNU
      FN = MAX(17.0D0,FN)
      IF (X.GE.FN) GO TO 430
      ANS = MAX(36.0D0-FNU,0.0D0)
      NS = INT(ANS)
      FNI = FNI + NS
      DFN = FNI + FNF
      FN = DFN
      IS = KT
      KM = N - 1 + NS
      IF (KM.GT.0) IS = 3
      GO TO 120
   90 FN = FNU
      FNP1 = FN + 1.0D0
      XO2L = LOG(XO2)
      IS = KT
      IF (X.LE.0.5D0) GO TO 230
      NS = 0
  100 FNI = FNI + NS
      DFN = FNI + FNF
      FN = DFN
      FNP1 = FN + 1.0D0
      IS = KT
      IF (N-1+NS.GT.0) IS = 3
      GO TO 230
  110 XO2L = LOG(XO2)
      NS = INT(SXO2-FNU)
      GO TO 100
  120 CONTINUE
!
!     OVERFLOW TEST ON UNIFORM ASYMPTOTIC EXPANSION
!
      IF (KODE.EQ.2) GO TO 130
      IF (ALPHA.LT.1.0D0) GO TO 150
      Z = X/ALPHA
      RA = SQRT(1.0D0+Z*Z)
      GLN = LOG((1.0D0+RA)/Z)
      T = RA*(1.0D0-ETX) + ETX/(Z+RA)
      ARG = ALPHA*(T-GLN)
      IF (ARG.GT.ELIM) GO TO 610
      IF (KM.EQ.0) GO TO 140
  130 CONTINUE
!
!     UNDERFLOW TEST ON UNIFORM ASYMPTOTIC EXPANSION
!
      Z = X/FN
      RA = SQRT(1.0D0+Z*Z)
      GLN = LOG((1.0D0+RA)/Z)
      T = RA*(1.0D0-ETX) + ETX/(Z+RA)
      ARG = FN*(T-GLN)
  140 IF (ARG.LT.(-ELIM)) GO TO 280
      GO TO 190
  150 IF (X.GT.ELIM) GO TO 610
      GO TO 130
!
!     UNIFORM ASYMPTOTIC EXPANSION FOR NU TO INFINITY
!
  160 IF (KM.NE.0) GO TO 170
      Y(1) = TEMP(3)
      RETURN
  170 TEMP(1) = TEMP(3)
      IN = NS
      KT = 1
      I1 = 0
  180 CONTINUE
      IS = 2
      FNI = FNI - 1.0D0
      DFN = FNI + FNF
      FN = DFN
      IF(I1.EQ.2) GO TO 350
      Z = X/FN
      RA = SQRT(1.0D0+Z*Z)
      GLN = LOG((1.0D0+RA)/Z)
      T = RA*(1.0D0-ETX) + ETX/(Z+RA)
      ARG = FN*(T-GLN)
  190 CONTINUE
      I1 = ABS(3-IS)
      I1 = MAX(I1,1)
      FLGIK = 1.0D0
      CALL DASYIK(X,FN,KODE,FLGIK,RA,ARG,I1,TEMP(IS))
      GO TO (180, 350, 510), IS
!
!     SERIES FOR (X/2)**2.LE.NU+1
!
  230 CONTINUE
      GLN = DLNGAM(FNP1)
      ARG = FN*XO2L - GLN - SX
      IF (ARG.LT.(-ELIM)) GO TO 300
      EARG = EXP(ARG)
  240 CONTINUE
      S = 1.0D0
      IF (X.LT.TOL) GO TO 260
      AK = 3.0D0
      T2 = 1.0D0
      T = 1.0D0
      S1 = FN
      DO 250 K=1,17
        S2 = T2 + S1
        T = T*SXO2/S2
        S = S + T
        IF (ABS(T).LT.TOL) GO TO 260
        T2 = T2 + AK
        AK = AK + 2.0D0
        S1 = S1 + FN
  250 CONTINUE
  260 CONTINUE
      TEMP(IS) = S*EARG
      GO TO (270, 350, 500), IS
  270 EARG = EARG*FN/XO2
      FNI = FNI - 1.0D0
      DFN = FNI + FNF
      FN = DFN
      IS = 2
      GO TO 240
!
!     SET UNDERFLOW VALUE AND UPDATE PARAMETERS
!
  280 Y(NN) = 0.0D0
      NN = NN - 1
      FNI = FNI - 1.0D0
      DFN = FNI + FNF
      FN = DFN
!      IF (NN-1) 340, 290, 130
      IF (NN-1.lt.0) then
            go to 340
      else if (NN-1.eq.0) then
            go to 290
      else
            go to 130
      end if

  290 KT = 2
      IS = 2
      GO TO 130
  300 Y(NN) = 0.0D0
      NN = NN - 1
      FNP1 = FN
      FNI = FNI - 1.0D0
      DFN = FNI + FNF
      FN = DFN
!      IF (NN-1) 340, 310, 320
      IF (NN-1.lt.0) then
            go to 340
      else if (NN-1.eq.0) then
            go to 310
      else
            go to 320
      end if

  310 KT = 2
      IS = 2
  320 IF (SXO2.LE.FNP1) GO TO 330
      GO TO 130
  330 ARG = ARG - XO2L + LOG(FNP1)
      IF (ARG.LT.(-ELIM)) GO TO 300
      GO TO 230
  340 NZ = N - NN
      RETURN
!
!     BACKWARD RECURSION SECTION
!
  350 CONTINUE
      NZ = N - NN
  360 CONTINUE
      IF(KT.EQ.2) GO TO 420
      S1 = TEMP(1)
      S2 = TEMP(2)
      TRX = 2.0D0/X
      DTM = FNI
      TM = (DTM+FNF)*TRX
      IF (IN.EQ.0) GO TO 390
!     BACKWARD RECUR TO INDEX ALPHA+NN-1
      DO 380 I=1,IN
        S = S2
        S2 = TM*S2 + S1
        S1 = S
        DTM = DTM - 1.0D0
        TM = (DTM+FNF)*TRX
  380 CONTINUE
      Y(NN) = S1
      IF (NN.EQ.1) RETURN
      Y(NN-1) = S2
      IF (NN.EQ.2) RETURN
      GO TO 400
  390 CONTINUE
!     BACKWARD RECUR FROM INDEX ALPHA+NN-1 TO ALPHA
      Y(NN) = S1
      Y(NN-1) = S2
      IF (NN.EQ.2) RETURN
  400 K = NN + 1
      DO 410 I=3,NN
        K = K - 1
        Y(K-2) = TM*Y(K-1) + Y(K)
        DTM = DTM - 1.0D0
        TM = (DTM+FNF)*TRX
  410 CONTINUE
      RETURN
  420 Y(1) = TEMP(2)
      RETURN
!
!     ASYMPTOTIC EXPANSION FOR X TO INFINITY
!
  430 CONTINUE
      EARG = RTTPI/SQRT(X)
      IF (KODE.EQ.2) GO TO 440
      IF (X.GT.ELIM) GO TO 610
      EARG = EARG*EXP(X)
  440 ETX = 8.0D0*X
      IS = KT
      IN = 0
      FN = FNU
  450 DX = FNI + FNI
      TM = 0.0D0
      IF (FNI.EQ.0.0D0 .AND. ABS(FNF).LT.TOL) GO TO 460
      TM = 4.0D0*FNF*(FNI+FNI+FNF)
  460 CONTINUE
      DTM = DX*DX
      S1 = ETX
      TRX = DTM - 1.0D0
      DX = -(TRX+TM)/ETX
      T = DX
      S = 1.0D0 + DX
      ATOL = TOL*ABS(S)
      S2 = 1.0D0
      AK = 8.0D0
      DO 470 K=1,25
        S1 = S1 + ETX
        S2 = S2 + AK
        DX = DTM - S2
        AP = DX + TM
        T = -T*AP/S1
        S = S + T
        IF (ABS(T).LE.ATOL) GO TO 480
        AK = AK + 8.0D0
  470 CONTINUE
  480 TEMP(IS) = S*EARG
      IF(IS.EQ.2) GO TO 360
      IS = 2
      FNI = FNI - 1.0D0
      DFN = FNI + FNF
      FN = DFN
      GO TO 450
!
!     BACKWARD RECURSION WITH NORMALIZATION BY
!     ASYMPTOTIC EXPANSION FOR NU TO INFINITY OR POWER SERIES.
!
  500 CONTINUE
!     COMPUTATION OF LAST ORDER FOR SERIES NORMALIZATION
      AKM = MAX(3.0D0-FN,0.0D0)
      KM = INT(AKM)
      TFN = FN + KM
      TA = (GLN+TFN-0.9189385332D0-0.0833333333D0/TFN)/(TFN+0.5D0)
      TA = XO2L - TA
      TB = -(1.0D0-1.0D0/TFN)/TFN
      AIN = TOLLN/(-TA+SQRT(TA*TA-TOLLN*TB)) + 1.5D0
      IN = INT(AIN)
      IN = IN + KM
      GO TO 520
  510 CONTINUE
!     COMPUTATION OF LAST ORDER FOR ASYMPTOTIC EXPANSION NORMALIZATION
      T = 1.0D0/(FN*RA)
      AIN = TOLLN/(GLN+SQRT(GLN*GLN+T*TOLLN)) + 1.5D0
      IN = INT(AIN)
      IF (IN.GT.INLIM) GO TO 160
  520 CONTINUE
      TRX = 2.0D0/X
      DTM = FNI + IN
      TM = (DTM+FNF)*TRX
      TA = 0.0D0
      TB = TOL
      KK = 1
  530 CONTINUE
!
!     BACKWARD RECUR UNINDEXED
!
      DO 540 I=1,IN
        S = TB
        TB = TM*TB + TA
        TA = S
        DTM = DTM - 1.0D0
        TM = (DTM+FNF)*TRX
  540 CONTINUE
!     NORMALIZATION
      IF (KK.NE.1) GO TO 550
      TA = (TA/TB)*TEMP(3)
      TB = TEMP(3)
      KK = 2
      IN = NS
      IF (NS.NE.0) GO TO 530
  550 Y(NN) = TB
      NZ = N - NN
      IF (NN.EQ.1) RETURN
      TB = TM*TB + TA
      K = NN - 1
      Y(K) = TB
      IF (NN.EQ.2) RETURN
      DTM = DTM - 1.0D0
      TM = (DTM+FNF)*TRX
      KM = K - 1
!
!     BACKWARD RECUR INDEXED
!
      DO 560 I=1,KM
        Y(K-1) = TM*Y(K) + Y(K+1)
        DTM = DTM - 1.0D0
        TM = (DTM+FNF)*TRX
        K = K - 1
  560 CONTINUE
      RETURN
!
!
!
  570 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESI', 'SCALING OPTION, KODE, NOT 1 OR 2.', 2, 1)
      RETURN
  580 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESI', 'ORDER, ALPHA, LESS THAN ZERO.', 2, 1)
      RETURN
  590 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESI', 'N LESS THAN ONE.', 2, 1)
      RETURN
  600 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESI', 'X LESS THAN ZERO.', 2, 1)
      RETURN
  610 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESI', 'OVERFLOW, X TOO LARGE FOR KODE = 1.', 6, 1)
      RETURN
END SUBROUTINE DBESI

!--- DECK DBESK
!
!=============================================================================!
!
SUBROUTINE DBESK (X, FNU, KODE, N, Y, NZ)
!
!=============================================================================!
!
!***BEGIN PROLOGUE  DBESK
!***PURPOSE  Implement forward recursion on the three term recursion
!            relation for a sequence of non-negative order Bessel
!            functions K/SUB(FNU+I-1)/(X), or scaled Bessel functions
!            EXP(X)*K/SUB(FNU+I-1)/(X), I=1,...,N for real, positive
!            X and non-negative orders FNU.
!***LIBRARY   SLATEC
!***CATEGORY  C10B3
!***TYPE      DOUBLE PRECISION (BESK-S, DBESK-D)
!***KEYWORDS  K BESSEL FUNCTION, SPECIAL FUNCTIONS
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Abstract  **** a double precision routine ****
!         DBESK implements forward recursion on the three term
!         recursion relation for a sequence of non-negative order Bessel
!         functions K/sub(FNU+I-1)/(X), or scaled Bessel functions
!         EXP(X)*K/sub(FNU+I-1)/(X), I=1,..,N for real X .GT. 0.0D0 and
!         non-negative orders FNU.  If FNU .LT. NULIM, orders FNU and
!         FNU+1 are obtained from DBSKNU to start the recursion.  If
!         FNU .GE. NULIM, the uniform asymptotic expansion is used for
!         orders FNU and FNU+1 to start the recursion.  NULIM is 35 or
!         70 depending on whether N=1 or N .GE. 2.  Under and overflow
!         tests are made on the leading term of the asymptotic expansion
!         before any extensive computation is done.
!
!         The maximum number of significant digits obtainable
!         is the smaller of 14 and the number of digits carried in
!         double precision arithmetic.
!
!     Description of Arguments
!
!         Input      X,FNU are double precision
!           X      - X .GT. 0.0D0
!           FNU    - order of the initial K function, FNU .GE. 0.0D0
!           KODE   - a parameter to indicate the scaling option
!                    KODE=1 returns Y(I)=       K/sub(FNU+I-1)/(X),
!                                        I=1,...,N
!                    KODE=2 returns Y(I)=EXP(X)*K/sub(FNU+I-1)/(X),
!                                        I=1,...,N
!           N      - number of members in the sequence, N .GE. 1
!
!         Output     Y is double precision
!           Y      - a vector whose first N components contain values
!                    for the sequence
!                    Y(I)=       k/sub(FNU+I-1)/(X), I=1,...,N  or
!                    Y(I)=EXP(X)*K/sub(FNU+I-1)/(X), I=1,...,N
!                    depending on KODE
!           NZ     - number of components of Y set to zero due to
!                    underflow with KODE=1,
!                    NZ=0   , normal return, computation completed
!                    NZ .NE. 0, first NZ components of Y set to zero
!                             due to underflow, Y(I)=0.0D0, I=1,...,NZ
!
!     Error Conditions
!         Improper input arguments - a fatal error
!         Overflow - a fatal error
!         Underflow with KODE=1 -  a non-fatal error (NZ .NE. 0)
!
!***REFERENCES  F. W. J. Olver, Tables of Bessel Functions of Moderate
!                 or Large Orders, NPL Mathematical Tables 6, Her
!                 Majesty's Stationery Office, London, 1962.
!               N. M. Temme, On the numerical evaluation of the modified
!                 Bessel function of the third kind, Journal of
!                 Computational Physics 19, (1975), pp. 324-337.
!***ROUTINES CALLED  D1MACH, DASYIK, DBESK0, DBESK1, DBSK0E, DBSK1E,
!                    DBSKNU, I1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790201  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DBESK
!
      INTEGER I, J, K, KODE, MZ, N, NB, ND, NN, NUD, NULIM, NZ
      DOUBLE PRECISION CN,DNU,ELIM,ETX,FLGIK,FN,FNN,FNU,GLN,GNU,RTZ
      DOUBLE PRECISION S, S1, S2, T, TM, TRX, W, X, XLIM, Y, ZN
      DIMENSION W(2), NULIM(2), Y(*)
      SAVE NULIM
      DATA NULIM(1),NULIM(2) / 35 , 70 /
!***FIRST EXECUTABLE STATEMENT  DBESK
      NN = -I1MACH(15)
      ELIM = 2.303D0*(NN*D1MACH(5)-3.0D0)
      XLIM = D1MACH(1)*1.0D+3
      IF (KODE.LT.1 .OR. KODE.GT.2) GO TO 280
      IF (FNU.LT.0.0D0) GO TO 290
      IF (X.LE.0.0D0) GO TO 300
      IF (X.LT.XLIM) GO TO 320
      IF (N.LT.1) GO TO 310
      ETX = KODE - 1
!
!     ND IS A DUMMY VARIABLE FOR N
!     GNU IS A DUMMY VARIABLE FOR FNU
!     NZ = NUMBER OF UNDERFLOWS ON KODE=1
!
      ND = N
      NZ = 0
      NUD = INT(FNU)
      DNU = FNU - NUD
      GNU = FNU
      NN = MIN(2,ND)
      FN = FNU + N - 1
      FNN = FN
      IF (FN.LT.2.0D0) GO TO 150
!
!     OVERFLOW TEST  (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION)
!     FOR THE LAST ORDER, FNU+N-1.GE.NULIM
!
      ZN = X/FN
      IF (ZN.EQ.0.0D0) GO TO 320
      RTZ = SQRT(1.0D0+ZN*ZN)
      GLN = LOG((1.0D0+RTZ)/ZN)
      T = RTZ*(1.0D0-ETX) + ETX/(ZN+RTZ)
      CN = -FN*(T-GLN)
      IF (CN.GT.ELIM) GO TO 320
      IF (NUD.LT.NULIM(NN)) GO TO 30
      IF (NN.EQ.1) GO TO 20
   10 CONTINUE
!
!     UNDERFLOW TEST (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION)
!     FOR THE FIRST ORDER, FNU.GE.NULIM
!
      FN = GNU
      ZN = X/FN
      RTZ = SQRT(1.0D0+ZN*ZN)
      GLN = LOG((1.0D0+RTZ)/ZN)
      T = RTZ*(1.0D0-ETX) + ETX/(ZN+RTZ)
      CN = -FN*(T-GLN)
   20 CONTINUE
      IF (CN.LT.-ELIM) GO TO 230
!
!     ASYMPTOTIC EXPANSION FOR ORDERS FNU AND FNU+1.GE.NULIM
!
      FLGIK = -1.0D0
      CALL DASYIK(X,GNU,KODE,FLGIK,RTZ,CN,NN,Y)
      IF (NN.EQ.1) GO TO 240
      TRX = 2.0D0/X
      TM = (GNU+GNU+2.0D0)/X
      GO TO 130
!
   30 CONTINUE
      IF (KODE.EQ.2) GO TO 40
!
!     UNDERFLOW TEST (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION IN X)
!     FOR ORDER DNU
!
      IF (X.GT.ELIM) GO TO 230
   40 CONTINUE
      IF (DNU.NE.0.0D0) GO TO 80
      IF (KODE.EQ.2) GO TO 50
      S1 = DBESK0(X)
      GO TO 60
   50 S1 = DBSK0E(X)
   60 CONTINUE
      IF (NUD.EQ.0 .AND. ND.EQ.1) GO TO 120
      IF (KODE.EQ.2) GO TO 70
      S2 = DBESK1(X)
      GO TO 90
   70 S2 = DBSK1E(X)
      GO TO 90
   80 CONTINUE
      NB = 2
      IF (NUD.EQ.0 .AND. ND.EQ.1) NB = 1
      CALL DBSKNU(X, DNU, KODE, NB, W, NZ)
      S1 = W(1)
      IF (NB.EQ.1) GO TO 120
      S2 = W(2)
   90 CONTINUE
      TRX = 2.0D0/X
      TM = (DNU+DNU+2.0D0)/X
!     FORWARD RECUR FROM DNU TO FNU+1 TO GET Y(1) AND Y(2)
      IF (ND.EQ.1) NUD = NUD - 1
      IF (NUD.GT.0) GO TO 100
      IF (ND.GT.1) GO TO 120
      S1 = S2
      GO TO 120
  100 CONTINUE
      DO 110 I=1,NUD
        S = S2
        S2 = TM*S2 + S1
        S1 = S
        TM = TM + TRX
  110 CONTINUE
      IF (ND.EQ.1) S1 = S2
  120 CONTINUE
      Y(1) = S1
      IF (ND.EQ.1) GO TO 240
      Y(2) = S2
  130 CONTINUE
      IF (ND.EQ.2) GO TO 240
!     FORWARD RECUR FROM FNU+2 TO FNU+N-1
      DO 140 I=3,ND
        Y(I) = TM*Y(I-1) + Y(I-2)
        TM = TM + TRX
  140 CONTINUE
      GO TO 240
!
  150 CONTINUE
!     UNDERFLOW TEST FOR KODE=1
      IF (KODE.EQ.2) GO TO 160
      IF (X.GT.ELIM) GO TO 230
  160 CONTINUE
!     OVERFLOW TEST
      IF (FN.LE.1.0D0) GO TO 170
      IF (-FN*(LOG(X)-0.693D0).GT.ELIM) GO TO 320
  170 CONTINUE
      IF (DNU.EQ.0.0D0) GO TO 180
      CALL DBSKNU(X, FNU, KODE, ND, Y, MZ)
      GO TO 240
  180 CONTINUE
      J = NUD
      IF (J.EQ.1) GO TO 210
      J = J + 1
      IF (KODE.EQ.2) GO TO 190
      Y(J) = DBESK0(X)
      GO TO 200
  190 Y(J) = DBSK0E(X)
  200 IF (ND.EQ.1) GO TO 240
      J = J + 1
  210 IF (KODE.EQ.2) GO TO 220
      Y(J) = DBESK1(X)
      GO TO 240
  220 Y(J) = DBSK1E(X)
      GO TO 240
!
!     UPDATE PARAMETERS ON UNDERFLOW
!
  230 CONTINUE
      NUD = NUD + 1
      ND = ND - 1
      IF (ND.EQ.0) GO TO 240
      NN = MIN(2,ND)
      GNU = GNU + 1.0D0
      IF (FNN.LT.2.0D0) GO TO 230
      IF (NUD.LT.NULIM(NN)) GO TO 230
      GO TO 10
  240 CONTINUE
      NZ = N - ND
      IF (NZ.EQ.0) RETURN
      IF (ND.EQ.0) GO TO 260
      DO 250 I=1,ND
        J = N - I + 1
        K = ND - I + 1
        Y(J) = Y(K)
  250 CONTINUE
  260 CONTINUE
      DO 270 I=1,NZ
        Y(I) = 0.0D0
  270 CONTINUE
      RETURN
!
!
!
  280 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESK', 'SCALING OPTION, KODE, NOT 1 OR 2', 2, 1)
      RETURN
  290 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESK', 'ORDER, FNU, LESS THAN ZERO', 2, 1)
      RETURN
  300 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESK', 'X LESS THAN OR EQUAL TO ZERO', 2, 1)
      RETURN
  310 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESK', 'N LESS THAN ONE', 2, 1)
      RETURN
  320 CONTINUE
      CALL XERMSG ('SLATEC', 'DBESK', 'OVERFLOW, FNU OR N TOO LARGE OR X TOO SMALL', 6, 1)
      RETURN
END SUBROUTINE DBESK

END MODULE DBESIK
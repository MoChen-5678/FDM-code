!*****************************************************************************!
!
Module DEFQ2
!
!*****************************************************************************!
!
Use DImach
Use ERRORS
Implicit none
!

Contains
!DECK DENORM
!
!=============================================================================!
!
FUNCTION DENORM (N, X)
!
!=============================================================================!
!
!***BEGIN PROLOGUE  DENORM
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DNSQ and DNSQE
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (ENORM-S, DENORM-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     Given an N-vector X, this function calculates the
!     Euclidean norm of X.
!
!     The Euclidean norm is computed by accumulating the sum of
!     squares in three different sums. The sums of squares for the
!     small and large components are scaled so that no overflows
!     occur. Non-destructive underflows are permitted. Underflows
!     and overflows do not occur in the computation of the unscaled
!     sum of squares for the intermediate components.
!     The definitions of small, intermediate and large components
!     depend on two constants, RDWARF and RGIANT. The main
!     restrictions on these constants are that RDWARF**2 not
!     underflow and RGIANT**2 not overflow. The constants
!     given here are suitable for every known computer.
!
!     The function statement is
!
!       DOUBLE PRECISION FUNCTION DENORM(N,X)
!
!     where
!
!       N is a positive integer input variable.
!
!       X is an input array of length N.
!
!***SEE ALSO  DNSQ, DNSQE
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DENORM
      INTEGER I, N
      DOUBLE PRECISION AGIANT, FLOATN, ONE, RDWARF, RGIANT, S1, S2, S3, X(*), X1MAX, X3MAX, XABS, ZERO
      double precision DENORM
      SAVE ONE, ZERO, RDWARF, RGIANT
      DATA ONE/1.0D0/,ZERO/0.0D0/,RDWARF/3.834D-20/,RGIANT/1.304D19/
!***FIRST EXECUTABLE STATEMENT  DENORM
      S1 = ZERO
      S2 = ZERO
      S3 = ZERO
      X1MAX = ZERO
      X3MAX = ZERO
      FLOATN = N
      AGIANT = RGIANT/FLOATN
      DO 90 I = 1, N
         XABS = ABS(X(I))
         IF (XABS .GT. RDWARF .AND. XABS .LT. AGIANT) GO TO 70
            IF (XABS .LE. RDWARF) GO TO 30
!
!              SUM FOR LARGE COMPONENTS.
!
               IF (XABS .LE. X1MAX) GO TO 10
                  S1 = ONE + S1*(X1MAX/XABS)**2
                  X1MAX = XABS
                  GO TO 20
   10          CONTINUE
                  S1 = S1 + (XABS/X1MAX)**2
   20          CONTINUE
               GO TO 60
   30       CONTINUE
!
!              SUM FOR SMALL COMPONENTS.
!
               IF (XABS .LE. X3MAX) GO TO 40
                  S3 = ONE + S3*(X3MAX/XABS)**2
                  X3MAX = XABS
                  GO TO 50
   40          CONTINUE
                  IF (XABS .NE. ZERO) S3 = S3 + (XABS/X3MAX)**2
   50          CONTINUE
   60       CONTINUE
            GO TO 80
   70    CONTINUE
!
!           SUM FOR INTERMEDIATE COMPONENTS.
!
            S2 = S2 + XABS**2
   80    CONTINUE
   90    CONTINUE
!
!     CALCULATION OF NORM.
!
      IF (S1 .EQ. ZERO) GO TO 100
         DENORM = X1MAX*SQRT(S1+(S2/X1MAX)/X1MAX)
         GO TO 130
  100 CONTINUE
         IF (S2 .EQ. ZERO) GO TO 110
            IF (S2 .GE. X3MAX) DENORM = SQRT(S2*(ONE+(X3MAX/S2)*(X3MAX*S3)))
            IF (S2 .LT. X3MAX) DENORM = SQRT(X3MAX*((S2/X3MAX)+(X3MAX*S3)))
            GO TO 120
  110    CONTINUE
            DENORM = X3MAX*SQRT(S3)
  120    CONTINUE
  130 CONTINUE
      RETURN
!
!     LAST CARD OF FUNCTION DENORM.
!
END FUNCTION DENORM

!DECK DFDJC1
!=============================================================================!
!
SUBROUTINE DFDJC1 (FCN, N, X, FVEC, FJAC, LDFJAC, IFLAG, ML, MU, EPSFCN, WA1, WA2)
!
!=============================================================================!
!
!***BEGIN PROLOGUE  DFDJC1
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DNSQ and DNSQE
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (FDJAC1-S, DFDJC1-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     This subroutine computes a forward-difference approximation
!     to the N by N Jacobian matrix associated with a specified
!     problem of N functions in N variables. If the Jacobian has
!     a banded form, then function evaluations are saved by only
!     approximating the nonzero terms.
!
!     The subroutine statement is
!
!       SUBROUTINE DFDJC1(FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU,EPSFCN,
!                         WA1,WA2)
!
!     where
!
!       FCN is the name of the user-supplied subroutine which
!         calculates the functions. FCN must be declared
!         in an EXTERNAL statement in the user calling
!         program, and should be written as follows.
!
!         SUBROUTINE FCN(N,X,FVEC,IFLAG)
!         INTEGER N,IFLAG
!         DOUBLE PRECISION X(N),FVEC(N)
!         ----------
!         Calculate the functions at X and
!         return this vector in FVEC.
!         ----------
!         RETURN
!
!         The value of IFLAG should not be changed by FCN unless
!         the user wants to terminate execution of DFDJC1.
!         In this case set IFLAG to a negative integer.
!
!       N is a positive integer input variable set to the number
!         of functions and variables.
!
!       X is an input array of length N.
!
!       FVEC is an input array of length N which must contain the
!         functions evaluated at X.
!
!       FJAC is an output N by N array which contains the
!         approximation to the Jacobian matrix evaluated at X.
!
!       LDFJAC is a positive integer input variable not less than N
!         which specifies the leading dimension of the array FJAC.
!
!       IFLAG is an integer variable which can be used to terminate
!         the execution of DFDJC1. See description of FCN.
!
!       ML is a nonnegative integer input variable which specifies
!         the number of subdiagonals within the band of the
!         Jacobian matrix. If the Jacobian is not banded, set
!         ML to at least N - 1.
!
!       EPSFCN is an input variable used in determining a suitable
!         step length for the forward-difference approximation. This
!         approximation assumes that the relative errors in the
!         functions are of the order of EPSFCN. If EPSFCN is less
!         than the machine precision, it is assumed that the relative
!         errors in the functions are of the order of the machine
!         precision.
!
!       MU is a nonnegative integer input variable which specifies
!         the number of superdiagonals within the band of the
!         Jacobian matrix. If the Jacobian is not banded, set
!         MU to at least N - 1.
!
!       WA1 and WA2 are work arrays of length N. If ML + MU + 1 is at
!         least N, then the Jacobian is considered dense, and WA2 is
!         not referenced.
!
!***SEE ALSO  DNSQ, DNSQE
!***ROUTINES CALLED  D1MACH
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DFDJC1
      INTEGER I, IFLAG, J, K, LDFJAC, ML, MSUM, MU, N
      DOUBLE PRECISION EPS, EPSFCN, EPSMCH, FJAC(LDFJAC,*), FVEC(*), H, TEMP, WA1(*), WA2(*), X(*), ZERO
      External FCN
      SAVE ZERO
      DATA ZERO /0.0D0/
!
!     EPSMCH IS THE MACHINE PRECISION.
!
!***FIRST EXECUTABLE STATEMENT  DFDJC1
      EPSMCH = D1MACH(4)
!
      EPS = SQRT(MAX(EPSFCN,EPSMCH))
      MSUM = ML + MU + 1
      IF (MSUM .LT. N) GO TO 40
!
!        COMPUTATION OF DENSE APPROXIMATE JACOBIAN.
!
         DO 20 J = 1, N
            TEMP = X(J)
            H = EPS*ABS(TEMP)
            IF (H .EQ. ZERO) H = EPS
            X(J) = TEMP + H
            CALL FCN(N,X,WA1,IFLAG)
            IF (IFLAG .LT. 0) GO TO 30
            X(J) = TEMP
            DO 10 I = 1, N
               FJAC(I,J) = (WA1(I) - FVEC(I))/H
   10          CONTINUE
   20       CONTINUE
   30    CONTINUE
         GO TO 110
   40 CONTINUE
!
!        COMPUTATION OF BANDED APPROXIMATE JACOBIAN.
!
         DO 90 K = 1, MSUM
            DO 60 J = K, N, MSUM
               WA2(J) = X(J)
               H = EPS*ABS(WA2(J))
               IF (H .EQ. ZERO) H = EPS
               X(J) = WA2(J) + H
   60          CONTINUE
            CALL FCN(N,X,WA1,IFLAG)
            IF (IFLAG .LT. 0) GO TO 100
            DO 80 J = K, N, MSUM
               X(J) = WA2(J)
               H = EPS*ABS(WA2(J))
               IF (H .EQ. ZERO) H = EPS
               DO 70 I = 1, N
                  FJAC(I,J) = ZERO
                  IF (I .GE. J - MU .AND. I .LE. J + ML) FJAC(I,J) = (WA1(I) - FVEC(I))/H
   70             CONTINUE
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
  110 CONTINUE
      RETURN
!
!     LAST CARD OF SUBROUTINE DFDJC1.
!
END SUBROUTINE DFDJC1

!DECK DQFORM
!=============================================================================!
!
SUBROUTINE DQFORM (M, N, Q, LDQ, WA)
!
!=============================================================================!
!***BEGIN PROLOGUE  DQFORM
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DNSQ and DNSQE
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (QFORM-S, DQFORM-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     This subroutine proceeds from the computed QR factorization of
!     an M by N matrix A to accumulate the M by M orthogonal matrix
!     Q from its factored form.
!
!     The subroutine statement is
!
!       SUBROUTINE DQFORM(M,N,Q,LDQ,WA)
!
!     where
!
!       M is a positive integer input variable set to the number
!         of rows of A and the order of Q.
!
!       N is a positive integer input variable set to the number
!         of columns of A.
!
!       Q is an M by M array. On input the full lower trapezoid in
!         the first MIN(M,N) columns of Q contains the factored form.
!         On output Q has been accumulated into a square matrix.
!
!       LDQ is a positive integer input variable not less than M
!         which specifies the leading dimension of the array Q.
!
!       WA is a work array of length M.
!
!***SEE ALSO  DNSQ, DNSQE
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DQFORM
      INTEGER I, J, JM1, K, L, LDQ, M, MINMN, N, NP1
      DOUBLE PRECISION ONE, Q(LDQ,*), SUM, TEMP, WA(*), ZERO
      SAVE ONE, ZERO
      DATA ONE,ZERO /1.0D0,0.0D0/
!
!     ZERO OUT UPPER TRIANGLE OF Q IN THE FIRST MIN(M,N) COLUMNS.
!
!***FIRST EXECUTABLE STATEMENT  DQFORM
      MINMN = MIN(M,N)
      IF (MINMN .LT. 2) GO TO 30
      DO 20 J = 2, MINMN
         JM1 = J - 1
         DO 10 I = 1, JM1
            Q(I,J) = ZERO
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
!
!     INITIALIZE REMAINING COLUMNS TO THOSE OF THE IDENTITY MATRIX.
!
      NP1 = N + 1
      IF (M .LT. NP1) GO TO 60
      DO 50 J = NP1, M
         DO 40 I = 1, M
            Q(I,J) = ZERO
   40       CONTINUE
         Q(J,J) = ONE
   50    CONTINUE
   60 CONTINUE
!
!     ACCUMULATE Q FROM ITS FACTORED FORM.
!
      DO 120 L = 1, MINMN
         K = MINMN - L + 1
         DO 70 I = K, M
            WA(I) = Q(I,K)
            Q(I,K) = ZERO
   70       CONTINUE
         Q(K,K) = ONE
         IF (WA(K) .EQ. ZERO) GO TO 110
         DO 100 J = K, M
            SUM = ZERO
            DO 80 I = K, M
               SUM = SUM + Q(I,J)*WA(I)
   80          CONTINUE
            TEMP = SUM/WA(K)
            DO 90 I = K, M
               Q(I,J) = Q(I,J) - TEMP*WA(I)
   90          CONTINUE
  100       CONTINUE
  110    CONTINUE
  120    CONTINUE
      RETURN
!
!     LAST CARD OF SUBROUTINE DQFORM.
!
END SUBROUTINE DQFORM

!DECK DQRFAC
!=============================================================================!
!
SUBROUTINE DQRFAC (M, N, A, LDA, PIVOT, IPVT, LIPVT, SIGMA, ACNORM, WA)
!
!=============================================================================!
!
!***BEGIN PROLOGUE  DQRFAC
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DNLS1, DNLS1E, DNSQ and DNSQE
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (QRFAC-S, DQRFAC-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!   **** Double Precision version of QRFAC ****
!
!     This subroutine uses Householder transformations with column
!     pivoting (optional) to compute a QR factorization of the
!     M by N matrix A. That is, DQRFAC determines an orthogonal
!     matrix Q, a permutation matrix P, and an upper trapezoidal
!     matrix R with diagonal elements of nonincreasing magnitude,
!     such that A*P = Q*R. The Householder transformation for
!     column K, K = 1,2,...,MIN(M,N), is of the form
!
!                           T
!           I - (1/U(K))*U*U
!
!     where U has zeros in the first K-1 positions. The form of
!     this transformation and the method of pivoting first
!     appeared in the corresponding LINPACK subroutine.
!
!     The subroutine statement is
!
!       SUBROUTINE DQRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,SIGMA,ACNORM,WA)
!
!     where
!
!       M is a positive integer input variable set to the number
!         of rows of A.
!
!       N is a positive integer input variable set to the number
!         of columns of A.
!
!       A is an M by N array. On input A contains the matrix for
!         which the QR factorization is to be computed. On output
!         the strict upper trapezoidal part of A contains the strict
!         upper trapezoidal part of R, and the lower trapezoidal
!         part of A contains a factored form of Q (the non-trivial
!         elements of the U vectors described above).
!
!       LDA is a positive integer input variable not less than M
!         which specifies the leading dimension of the array A.
!
!       PIVOT is a logical input variable. If pivot is set .TRUE.,
!         then column pivoting is enforced. If pivot is set .FALSE.,
!         then no column pivoting is done.
!
!       IPVT is an integer output array of length LIPVT. IPVT
!         defines the permutation matrix P such that A*P = Q*R.
!         Column J of P is column IPVT(J) of the identity matrix.
!         If pivot is .FALSE., IPVT is not referenced.
!
!       LIPVT is a positive integer input variable. If PIVOT is
!             .FALSE., then LIPVT may be as small as 1. If PIVOT is
!             .TRUE., then LIPVT must be at least N.
!
!       SIGMA is an output array of length N which contains the
!         diagonal elements of R.
!
!       ACNORM is an output array of length N which contains the
!         norms of the corresponding columns of the input matrix A.
!         If this information is not needed, then ACNORM can coincide
!         with SIGMA.
!
!       WA is a work array of length N. If pivot is .FALSE., then WA
!         can coincide with SIGMA.
!
!***SEE ALSO  DNLS1, DNLS1E, DNSQ, DNSQE
!***ROUTINES CALLED  D1MACH, DENORM
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DQRFAC
      INTEGER M,N,LDA,LIPVT
      INTEGER IPVT(*)
      LOGICAL PIVOT
      SAVE ONE, P05, ZERO
      DOUBLE PRECISION A(LDA,*),SIGMA(*),ACNORM(*),WA(*)
      INTEGER I,J,JP1,K,KMAX,MINMN
      DOUBLE PRECISION AJNORM,EPSMCH,ONE,P05,SUM,TEMP,ZERO
      DATA ONE,P05,ZERO /1.0D0,5.0D-2,0.0D0/
!***FIRST EXECUTABLE STATEMENT  DQRFAC
      EPSMCH = D1MACH(4)
!
!     COMPUTE THE INITIAL COLUMN NORMS AND INITIALIZE SEVERAL ARRAYS.
!
      DO 10 J = 1, N
         ACNORM(J) = DENORM(M,A(1,J))
         SIGMA(J) = ACNORM(J)
         WA(J) = SIGMA(J)
         IF (PIVOT) IPVT(J) = J
   10    CONTINUE
!
!     REDUCE A TO R WITH HOUSEHOLDER TRANSFORMATIONS.
!
      MINMN = MIN(M,N)
      DO 110 J = 1, MINMN
         IF (.NOT.PIVOT) GO TO 40
!
!        BRING THE COLUMN OF LARGEST NORM INTO THE PIVOT POSITION.
!
         KMAX = J
         DO 20 K = J, N
            IF (SIGMA(K) .GT. SIGMA(KMAX)) KMAX = K
   20       CONTINUE
         IF (KMAX .EQ. J) GO TO 40
         DO 30 I = 1, M
            TEMP = A(I,J)
            A(I,J) = A(I,KMAX)
            A(I,KMAX) = TEMP
   30       CONTINUE
         SIGMA(KMAX) = SIGMA(J)
         WA(KMAX) = WA(J)
         K = IPVT(J)
         IPVT(J) = IPVT(KMAX)
         IPVT(KMAX) = K
   40    CONTINUE
!
!        COMPUTE THE HOUSEHOLDER TRANSFORMATION TO REDUCE THE
!        J-TH COLUMN OF A TO A MULTIPLE OF THE J-TH UNIT VECTOR.
!
         AJNORM = DENORM(M-J+1,A(J,J))
         IF (AJNORM .EQ. ZERO) GO TO 100
         IF (A(J,J) .LT. ZERO) AJNORM = -AJNORM
         DO 50 I = J, M
            A(I,J) = A(I,J)/AJNORM
   50       CONTINUE
         A(J,J) = A(J,J) + ONE
!
!        APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS
!        AND UPDATE THE NORMS.
!
         JP1 = J + 1
         IF (N .LT. JP1) GO TO 100
         DO 90 K = JP1, N
            SUM = ZERO
            DO 60 I = J, M
               SUM = SUM + A(I,J)*A(I,K)
   60          CONTINUE
            TEMP = SUM/A(J,J)
            DO 70 I = J, M
               A(I,K) = A(I,K) - TEMP*A(I,J)
   70          CONTINUE
            IF (.NOT.PIVOT .OR. SIGMA(K) .EQ. ZERO) GO TO 80
            TEMP = A(J,K)/SIGMA(K)
            SIGMA(K) = SIGMA(K)*SQRT(MAX(ZERO,ONE-TEMP**2))
            IF (P05*(SIGMA(K)/WA(K))**2 .GT. EPSMCH) GO TO 80
            SIGMA(K) = DENORM(M-J,A(JP1,K))
            WA(K) = SIGMA(K)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
         SIGMA(J) = -AJNORM
  110    CONTINUE
      RETURN
!
!     LAST CARD OF SUBROUTINE DQRFAC.
!
END SUBROUTINE DQRFAC

!DECK DDOGLG
!=============================================================================!
!
SUBROUTINE DDOGLG (N, R, LR, DIAG, QTB, DELTA, X, WA1, WA2)
!
!=============================================================================!
!
!***BEGIN PROLOGUE  DDOGLG
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DNSQ and DNSQE
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (DOGLEG-S, DDOGLG-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     Given an M by N matrix A, an N by N nonsingular diagonal
!     matrix D, an M-vector B, and a positive number DELTA, the
!     problem is to determine the convex combination X of the
!     Gauss-Newton and scaled gradient directions that minimizes
!     (A*X - B) in the least squares sense, subject to the
!     restriction that the Euclidean norm of D*X be at most DELTA.
!
!     This subroutine completes the solution of the problem
!     if it is provided with the necessary information from the
!     QR factorization of A. That is, if A = Q*R, where Q has
!     orthogonal columns and R is an upper triangular matrix,
!     then DDOGLG expects the full upper triangle of R and
!     the first N components of (Q transpose)*B.
!
!     The subroutine statement is
!
!       SUBROUTINE DDOGLG(N,R,LR,DIAG,QTB,DELTA,X,WA1,WA2)
!
!     where
!
!       N is a positive integer input variable set to the order of R.
!
!       R is an input array of length LR which must contain the upper
!         triangular matrix R stored by rows.
!
!       LR is a positive integer input variable not less than
!         (N*(N+1))/2.
!
!       DIAG is an input array of length N which must contain the
!         diagonal elements of the matrix D.
!
!       QTB is an input array of length N which must contain the first
!         N elements of the vector (Q transpose)*B.
!
!       DELTA is a positive input variable which specifies an upper
!         bound on the Euclidean norm of D*X.
!
!       X is an output array of length N which contains the desired
!         convex combination of the Gauss-Newton direction and the
!         scaled gradient direction.
!
!       WA1 and WA2 are work arrays of length N.
!
!***SEE ALSO  DNSQ, DNSQE
!***ROUTINES CALLED  D1MACH, DENORM
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DDOGLG
      INTEGER I, J, JJ, JP1, K, L, LR, N
      DOUBLE PRECISION ALPHA, BNORM, DELTA, DIAG(*), EPSMCH, GNORM, ONE, QNORM, QTB(*), &
            & R(*), SGNORM, SUM, TEMP, WA1(*), WA2(*), X(*), ZERO
      SAVE ONE, ZERO
      DATA ONE/1.0D0/, ZERO/0.0D0/
!
!     EPSMCH IS THE MACHINE PRECISION.
!
!***FIRST EXECUTABLE STATEMENT  DDOGLG
      EPSMCH = D1MACH(4)
!
!     FIRST, CALCULATE THE GAUSS-NEWTON DIRECTION.
!
      JJ = (N*(N + 1))/2 + 1
      DO 50 K = 1, N
         J = N - K + 1
         JP1 = J + 1
         JJ = JJ - K
         L = JJ + 1
         SUM = ZERO
         IF (N .LT. JP1) GO TO 20
         DO 10 I = JP1, N
            SUM = SUM + R(L)*X(I)
            L = L + 1
   10       CONTINUE
   20    CONTINUE
         TEMP = R(JJ)
         IF (TEMP .NE. ZERO) GO TO 40
         L = J
         DO 30 I = 1, J
            TEMP = MAX(TEMP,ABS(R(L)))
            L = L + N - I
   30       CONTINUE
         TEMP = EPSMCH*TEMP
         IF (TEMP .EQ. ZERO) TEMP = EPSMCH
   40    CONTINUE
         X(J) = (QTB(J) - SUM)/TEMP
   50    CONTINUE
!
!     TEST WHETHER THE GAUSS-NEWTON DIRECTION IS ACCEPTABLE.
!
      DO 60 J = 1, N
         WA1(J) = ZERO
         WA2(J) = DIAG(J)*X(J)
   60    CONTINUE
      QNORM = DENORM(N,WA2)
      IF (QNORM .LE. DELTA) GO TO 140
!
!     THE GAUSS-NEWTON DIRECTION IS NOT ACCEPTABLE.
!     NEXT, CALCULATE THE SCALED GRADIENT DIRECTION.
!
      L = 1
      DO 80 J = 1, N
         TEMP = QTB(J)
         DO 70 I = J, N
            WA1(I) = WA1(I) + R(L)*TEMP
            L = L + 1
   70       CONTINUE
         WA1(J) = WA1(J)/DIAG(J)
   80    CONTINUE
!
!     CALCULATE THE NORM OF THE SCALED GRADIENT AND TEST FOR
!     THE SPECIAL CASE IN WHICH THE SCALED GRADIENT IS ZERO.
!
      GNORM = DENORM(N,WA1)
      SGNORM = ZERO
      ALPHA = DELTA/QNORM
      IF (GNORM .EQ. ZERO) GO TO 120
!
!     CALCULATE THE POINT ALONG THE SCALED GRADIENT
!     AT WHICH THE QUADRATIC IS MINIMIZED.
!
      DO 90 J = 1, N
         WA1(J) = (WA1(J)/GNORM)/DIAG(J)
   90    CONTINUE
      L = 1
      DO 110 J = 1, N
         SUM = ZERO
         DO 100 I = J, N
            SUM = SUM + R(L)*WA1(I)
            L = L + 1
  100       CONTINUE
         WA2(J) = SUM
  110    CONTINUE
      TEMP = DENORM(N,WA2)
      SGNORM = (GNORM/TEMP)/TEMP
!
!     TEST WHETHER THE SCALED GRADIENT DIRECTION IS ACCEPTABLE.
!
      ALPHA = ZERO
      IF (SGNORM .GE. DELTA) GO TO 120
!
!     THE SCALED GRADIENT DIRECTION IS NOT ACCEPTABLE.
!     FINALLY, CALCULATE THE POINT ALONG THE DOGLEG
!     AT WHICH THE QUADRATIC IS MINIMIZED.
!
      BNORM = DENORM(N,QTB)
      TEMP = (BNORM/GNORM)*(BNORM/QNORM)*(SGNORM/DELTA)
      TEMP = TEMP - (DELTA/QNORM)*(SGNORM/DELTA)**2 + SQRT((TEMP-(DELTA/QNORM))**2 &
            & +(ONE-(DELTA/QNORM)**2)*(ONE-(SGNORM/DELTA)**2))
      ALPHA = ((DELTA/QNORM)*(ONE - (SGNORM/DELTA)**2))/TEMP
  120 CONTINUE
!
!     FORM APPROPRIATE CONVEX COMBINATION OF THE GAUSS-NEWTON
!     DIRECTION AND THE SCALED GRADIENT DIRECTION.
!
      TEMP = (ONE - ALPHA)*MIN(SGNORM,DELTA)
      DO 130 J = 1, N
         X(J) = TEMP*WA1(J) + ALPHA*X(J)
  130    CONTINUE
  140 CONTINUE
      RETURN
!
!     LAST CARD OF SUBROUTINE DDOGLG.
!
END SUBROUTINE DDOGLG

End Module DEFQ2
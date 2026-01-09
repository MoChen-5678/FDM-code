!*************************************************************************************************!
!
Module Eigenlib
!
!*************************************************************************************************!
!
Use LPKlib
IMPLICIT NONE

!--- DORMTR DORMQL DORM2L DSTEDC DSTERF

Contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
SUBROUTINE DORMTR( SIDE, UPLO, TRANS, M, N, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS, UPLO
      INTEGER :: INFO, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORMTR overwrites the general real M-by-N matrix C with
!
!                  SIDE = 'L'     SIDE = 'R'
!  TRANS = 'N':      Q * C          C * Q
!  TRANS = 'T':      Q**T * C       C * Q**T
!
!  where Q is a real orthogonal matrix of order nq, with nq = m if
!  SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
!  nq-1 elementary reflectors, as returned by DSYTRD:
!
!  if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1);
!
!  if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1).
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q**T from the Left;
!          = 'R': apply Q or Q**T from the Right.
!
!  UPLO    (input) CHARACTER*1
!          = 'U': Upper triangle of A contains elementary reflectors
!                 from DSYTRD;
!          = 'L': Lower triangle of A contains elementary reflectors
!                 from DSYTRD.
!
!  TRANS   (input) CHARACTER*1
!          = 'N':  No transpose, apply Q;
!          = 'T':  Transpose, apply Q**T.
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension
!                               (LDA,M) if SIDE = 'L'
!                               (LDA,N) if SIDE = 'R'
!          The vectors which define the elementary reflectors, as
!          returned by DSYTRD.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!          LDA >= max(1,M) if SIDE = 'L'; LDA >= max(1,N) if SIDE = 'R'.
!
!  TAU     (input) DOUBLE PRECISION array, dimension
!                               (M-1) if SIDE = 'L'
!                               (N-1) if SIDE = 'R'
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DSYTRD.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.
!          If SIDE = 'L', LWORK >= max(1,N);
!          if SIDE = 'R', LWORK >= max(1,M).
!          For optimum performance LWORK >= N*NB if SIDE = 'L', and
!          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
!          blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, UPPER
      INTEGER :: I1, I2, IINFO, LWKOPT, MI, NB, NI, NQ, NW, nik8, nk8
!     ..
!     .. External Functions ..
      !LOGICAL            LSAME
      !INTEGER            ILAENV
      !EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
      !EXTERNAL           DORMQL, DORMQR, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
!
!     NQ is the order of Q and NW is the minimum dimension of WORK
!
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.LSAME( TRANS, 'N' ) .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
!
      IF( INFO.EQ.0 ) THEN
         IF( UPPER ) THEN
            IF( LEFT ) THEN
               nik8=-1; nk8 = 1; NB = ILAENV( nk8, 'DORMQL', SIDE // TRANS, M-1, N, M-1, nik8 )
            ELSE
               nik8=-1; nk8 = 1; NB = ILAENV( nk8, 'DORMQL', SIDE // TRANS, M, N-1, N-1, nik8 )
            END IF
         ELSE
            IF( LEFT ) THEN
               nik8=-1; nk8 = 1; NB = ILAENV( nk8, 'DORMQR', SIDE // TRANS, M-1, N, M-1, nik8 )
            ELSE
               nik8=-1; nk8 = 1; NB = ILAENV( nk8, 'DORMQR', SIDE // TRANS, M, N-1, N-1, nik8 )
            END IF
         END IF
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORMTR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. NQ.EQ.1 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      IF( LEFT ) THEN
         MI = M - 1
         NI = N
      ELSE
         MI = M
         NI = N - 1
      END IF
!
      IF( UPPER ) THEN
!
!        Q was determined by a call to DSYTRD with UPLO = 'U'
!
         CALL DORMQL( SIDE, TRANS, MI, NI, NQ-1, A( 1, 2 ), LDA, TAU, C, LDC, WORK, LWORK, IINFO )
      ELSE
!
!        Q was determined by a call to DSYTRD with UPLO = 'L'
!
         IF( LEFT ) THEN
            I1 = 2
            I2 = 1
         ELSE
            I1 = 1
            I2 = 2
         END IF
         CALL DORMQR( SIDE, TRANS, MI, NI, NQ-1, A( 2, 1 ), LDA, TAU, &
     &                C( I1, I2 ), LDC, WORK, LWORK, IINFO )
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of DORMTR
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
END Subroutine DORMTR
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
SUBROUTINE DORMQL( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER:: INFO, K, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORMQL overwrites the general real M-by-N matrix C with
!
!                  SIDE = 'L'     SIDE = 'R'
!  TRANS = 'N':      Q * C          C * Q
!  TRANS = 'T':      Q**T * C       C * Q**T
!
!  where Q is a real orthogonal matrix defined as the product of k
!  elementary reflectors
!
!        Q = H(k) . . . H(2) H(1)
!
!  as returned by DGEQLF. Q is of order M if SIDE = 'L' and of order N
!  if SIDE = 'R'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q**T from the Left;
!          = 'R': apply Q or Q**T from the Right.
!
!  TRANS   (input) CHARACTER*1
!          = 'N':  No transpose, apply Q;
!          = 'T':  Transpose, apply Q**T.
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
!          The i-th column must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          DGEQLF in the last k columns of its array argument A.
!          A is modified by the routine but restored on exit.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!          If SIDE = 'L', LDA >= max(1,M);
!          if SIDE = 'R', LDA >= max(1,N).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQLF.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.
!          If SIDE = 'L', LWORK >= max(1,N);
!          if SIDE = 'R', LWORK >= max(1,M).
!          For optimum performance LWORK >= N*NB if SIDE = 'L', and
!          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
!          blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER :: NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      INTEGER :: I, I1, I2, I3, IB, IINFO, IWS, LDWORK, LWKOPT, MI, NB, NBMIN, NI, NQ, NW, nik8, nk8
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   T( LDT, NBMAX )
!     ..
!     .. External Functions ..
      !LOGICAL            LSAME
      !INTEGER            ILAENV
      !EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DLARFB, DLARFT, DORM2L, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
!
!     NQ is the order of Q and NW is the minimum dimension of WORK
!
      IF( LEFT ) THEN
         NQ = M
         NW = MAX( 1, N )
      ELSE
         NQ = N
         NW = MAX( 1, M )
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
!
      IF( INFO.EQ.0 ) THEN
         IF( M.EQ.0 .OR. N.EQ.0 ) THEN
            LWKOPT = 1
         ELSE
!
!           Determine the block size.  NB may be at most NBMAX, where
!           NBMAX is used to define the local array T.
!
            nik8=1; nk8 = -1; NB = MIN( NBMAX, ILAENV( nik8, 'DORMQL', SIDE // TRANS, M, N, K, nk8 ) )
            LWKOPT = NW*NB
         END IF
         WORK( 1 ) = LWKOPT
!
         IF( LWORK.LT.NW .AND. .NOT.LQUERY ) THEN
            INFO = -12
         END IF
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORMQL', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
         RETURN
      END IF
!
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IWS = NW*NB
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            nik8=2; nk8 = -1; NBMIN = MAX( 2, ILAENV( nik8, 'DORMQL', SIDE // TRANS, M, N, K, nk8 ) )
         END IF
      ELSE
         IWS = NW
      END IF
!
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
!
!        Use unblocked code
!
         CALL DORM2L( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, IINFO )
      ELSE
!
!        Use blocked code
!
         IF( ( LEFT .AND. NOTRAN ) .OR. ( .NOT.LEFT .AND. .NOT.NOTRAN ) ) THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         END IF
!
         IF( LEFT ) THEN
            NI = N
         ELSE
            MI = M
         END IF
!
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
!
!           Form the triangular factor of the block reflector
!           H = H(i+ib-1) . . . H(i+1) H(i)
!
            CALL DLARFT( 'Backward', 'Columnwise', NQ-K+I+IB-1, IB, &
     &                   A( 1, I ), LDA, TAU( I ), T, LDT )
            IF( LEFT ) THEN
!
!              H or H' is applied to C(1:m-k+i+ib-1,1:n)
!
               MI = M - K + I + IB - 1
            ELSE
!
!              H or H' is applied to C(1:m,1:n-k+i+ib-1)
!
               NI = N - K + I + IB - 1
            END IF
!
!           Apply H or H'
!
            CALL DLARFB( SIDE, TRANS, 'Backward', 'Columnwise', MI, NI, &
     &                   IB, A( 1, I ), LDA, T, LDT, C, LDC, WORK,      &
     &                   LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of DORMQL
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
END Subroutine DORMQL
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
SUBROUTINE DORM2L( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, INFO )
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER :: INFO, K, LDA, LDC, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORM2L overwrites the general real m by n matrix C with
!
!        Q * C  if SIDE = 'L' and TRANS = 'N', or
!
!        Q'* C  if SIDE = 'L' and TRANS = 'T', or
!
!        C * Q  if SIDE = 'R' and TRANS = 'N', or
!
!        C * Q' if SIDE = 'R' and TRANS = 'T',
!
!  where Q is a real orthogonal matrix defined as the product of k
!  elementary reflectors
!
!        Q = H(k) . . . H(2) H(1)
!
!  as returned by DGEQLF. Q is of order m if SIDE = 'L' and of order n
!  if SIDE = 'R'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q' from the Left
!          = 'R': apply Q or Q' from the Right
!
!  TRANS   (input) CHARACTER*1
!          = 'N': apply Q  (No transpose)
!          = 'T': apply Q' (Transpose)
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
!          The i-th column must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          DGEQLF in the last k columns of its array argument A.
!          A is modified by the routine but restored on exit.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!          If SIDE = 'L', LDA >= max(1,M);
!          if SIDE = 'R', LDA >= max(1,N).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQLF.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension
!                                   (N) if SIDE = 'L',
!                                   (M) if SIDE = 'R'
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER :: I, I1, I2, I3, MI, NI, NQ, nik8
      DOUBLE PRECISION   AII
!     ..
!     .. External Functions ..
!      LOGICAL            LSAME
!      EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
!      EXTERNAL           DLARF, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
!
!     NQ is the order of Q
!
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORM2L', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) RETURN
!
      IF( ( LEFT .AND. NOTRAN ) .OR. ( .NOT.LEFT .AND. .NOT.NOTRAN ) ) THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
!
      IF( LEFT ) THEN
         NI = N
      ELSE
         MI = M
      END IF
!
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
!
!           H(i) is applied to C(1:m-k+i,1:n)
!
            MI = M - K + I
         ELSE
!
!           H(i) is applied to C(1:m,1:n-k+i)
!
            NI = N - K + I
         END IF
!
!        Apply H(i)
!
         AII = A( NQ-K+I, I )
         A( NQ-K+I, I ) = ONE
         nik8 = 1; CALL DLARF( SIDE, MI, NI, A( 1, I ), nik8, TAU( I ), C, LDC, WORK )
         A( NQ-K+I, I ) = AII
   10 CONTINUE
      RETURN
!
!     End of DORM2L
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
END Subroutine DORM2L
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
SUBROUTINE DSTEDC( COMPZ, N, D, E, Z, LDZ, WORK, LWORK, IWORK, LIWORK, INFO )
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!  -- LAPACK driver routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER :: INFO, LDZ, LIWORK, LWORK, N
!     ..
!     .. Array Arguments ..
      INTEGER :: IWORK( * )
      DOUBLE PRECISION   D( * ), E( * ), WORK( * ), Z( LDZ, * )
!     ..
!
!  Purpose
!  =======
!
!  DSTEDC computes all eigenvalues and, optionally, eigenvectors of a
!  symmetric tridiagonal matrix using the divide and conquer method.
!  The eigenvectors of a full or band real symmetric matrix can also be
!  found if DSYTRD or DSPTRD or DSBTRD has been used to reduce this
!  matrix to tridiagonal form.
!
!  This code makes very mild assumptions about floating point
!  arithmetic. It will work on machines with a guard digit in
!  add/subtract, or on those binary machines without guard digits
!  which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
!  It could conceivably fail on hexadecimal or decimal machines
!  without guard digits, but we know of none.  See DLAED3 for details.
!
!  Arguments
!  =========
!
!  COMPZ   (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only.
!          = 'I':  Compute eigenvectors of tridiagonal matrix also.
!          = 'V':  Compute eigenvectors of original dense symmetric
!                  matrix also.  On entry, Z contains the orthogonal
!                  matrix used to reduce the original matrix to
!                  tridiagonal form.
!
!  N       (input) INTEGER
!          The dimension of the symmetric tridiagonal matrix.  N >= 0.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the diagonal elements of the tridiagonal matrix.
!          On exit, if INFO = 0, the eigenvalues in ascending order.
!
!  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, the subdiagonal elements of the tridiagonal matrix.
!          On exit, E has been destroyed.
!
!  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
!          On entry, if COMPZ = 'V', then Z contains the orthogonal
!          matrix used in the reduction to tridiagonal form.
!          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the
!          orthonormal eigenvectors of the original symmetric matrix,
!          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
!          of the symmetric tridiagonal matrix.
!          If  COMPZ = 'N', then Z is not referenced.
!
!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.  LDZ >= 1.
!          If eigenvectors are desired, then LDZ >= max(1,N).
!
!  WORK    (workspace/output) DOUBLE PRECISION array,
!                                         dimension (LWORK)
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.
!          If COMPZ = 'N' or N <= 1 then LWORK must be at least 1.
!          If COMPZ = 'V' and N > 1 then LWORK must be at least
!                         ( 1 + 3*N + 2*N*lg N + 3*N**2 ),
!                         where lg( N ) = smallest integer k such
!                         that 2**k >= N.
!          If COMPZ = 'I' and N > 1 then LWORK must be at least
!                         ( 1 + 4*N + N**2 ).
!          Note that for COMPZ = 'I' or 'V', then if N is less than or
!          equal to the minimum divide size, usually 25, then LWORK need
!          only be max(1,2*(N-1)).
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))
!          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
!
!  LIWORK  (input) INTEGER
!          The dimension of the array IWORK.
!          If COMPZ = 'N' or N <= 1 then LIWORK must be at least 1.
!          If COMPZ = 'V' and N > 1 then LIWORK must be at least
!                         ( 6 + 6*N + 5*N*lg N ).
!          If COMPZ = 'I' and N > 1 then LIWORK must be at least
!                         ( 3 + 5*N ).
!          Note that for COMPZ = 'I' or 'V', then if N is less than or
!          equal to the minimum divide size, usually 25, then LIWORK
!          need only be 1.
!
!          If LIWORK = -1, then a workspace query is assumed; the
!          routine only calculates the optimal size of the IWORK array,
!          returns this value as the first entry of the IWORK array, and
!          no error message related to LIWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  The algorithm failed to compute an eigenvalue while
!                working on the submatrix lying in rows and columns
!                INFO/(N+1) through mod(INFO,N+1).
!
!  Further Details
!  ===============
!
!  Based on contributions by
!     Jeff Rutter, Computer Science Division, University of California
!     at Berkeley, USA
!  Modified by Francoise Tisseur, University of Tennessee.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER :: FINISH, I, ICOMPZ, II, J, K, LGN, LIWMIN, &
     &                   LWMIN, M, SMLSIZ, START, STOREZ, STRTRW, nik8, nk8
      DOUBLE PRECISION   EPS, ORGNRM, P, TINY
!     ..
!     .. External Functions ..
!     LOGICAL            LSAME
!     INTEGER            ILAENV
!     DOUBLE PRECISION   DLAMCH, DLANST
!     EXTERNAL           LSAME, ILAENV, DLAMCH, DLANST
!     ..
!     .. External Subroutines ..
!     EXTERNAL           DGEMM, DLACPY, DLAED0, DLASCL, DLASET, DLASRT, &
!    &                   DSTEQR, DSTERF, DSWAP, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, INT, LOG, MAX, MOD, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 .OR. LIWORK.EQ.-1 )
!
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1, N ) ) ) THEN
         INFO = -6
      END IF
!
      IF( INFO.EQ.0 ) THEN
!
!        Compute the workspace requirements
!
         nik8 = 9; nk8 = 0; SMLSIZ = ILAENV( nik8, 'DSTEDC', ' ', nk8, nk8, nk8, nk8 )
         IF( N.LE.1 .OR. ICOMPZ.EQ.0 ) THEN
            LIWMIN = 1
            LWMIN = 1
         ELSE IF( N.LE.SMLSIZ ) THEN
            LIWMIN = 1
            LWMIN = 2*( N - 1 )
         ELSE
            LGN = INT( LOG( DBLE( N ) )/LOG( TWO ) )
            IF( 2**LGN.LT.N ) LGN = LGN + 1
            IF( 2**LGN.LT.N ) LGN = LGN + 1
            IF( ICOMPZ.EQ.1 ) THEN
               LWMIN = 1 + 3*N + 2*N*LGN + 3*N**2
               LIWMIN = 6 + 6*N + 5*N*LGN
            ELSE IF( ICOMPZ.EQ.2 ) THEN
               LWMIN = 1 + 4*N + N**2
               LIWMIN = 3 + 5*N
            END IF
         END IF
         WORK( 1 ) = LWMIN
         IWORK( 1 ) = LIWMIN
!
         IF( LWORK.LT.LWMIN .AND. .NOT. LQUERY ) THEN
            INFO = -8
         ELSE IF( LIWORK.LT.LIWMIN .AND. .NOT. LQUERY ) THEN
            INFO = -10
         END IF
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSTEDC', -INFO )
         RETURN
      ELSE IF (LQUERY) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) RETURN
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.NE.0 ) Z( 1, 1 ) = ONE
         RETURN
      END IF
!
!     If the following conditional clause is removed, then the routine
!     will use the Divide and Conquer routine to compute only the
!     eigenvalues, which requires (3N + 3N**2) real workspace and
!     (2 + 5N + 2N lg(N)) integer workspace.
!     Since on many architectures DSTERF is much faster than any other
!     algorithm for finding eigenvalues only, it is used here
!     as the default. If the conditional clause is removed, then
!     information on the size of workspace needs to be changed.
!
!     If COMPZ = 'N', use DSTERF to compute the eigenvalues.
!
      IF( ICOMPZ.EQ.0 ) THEN
         CALL DSTERF( N, D, E, INFO )
         GO TO 50
      END IF
!
!     If N is smaller than the minimum divide size (SMLSIZ+1), then
!     solve the problem with another solver.
!
      IF( N.LE.SMLSIZ ) THEN
!
         CALL DSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
!
      ELSE
!
!        If COMPZ = 'V', the Z matrix must be stored elsewhere for later
!        use.
!
         IF( ICOMPZ.EQ.1 ) THEN
            STOREZ = 1 + N*N
         ELSE
            STOREZ = 1
         END IF
!
         IF( ICOMPZ.EQ.2 ) THEN
            CALL DLASET( 'Full', N, N, ZERO, ONE, Z, LDZ )
         END IF
!
!        Scale.
!
         ORGNRM = DLANST( 'M', N, D, E )
         IF( ORGNRM.EQ.ZERO ) GO TO 50
!
         EPS = DLAMCH( 'Epsilon' )
!
         START = 1
!
!        while ( START <= N )
!
   10    CONTINUE
         IF( START.LE.N ) THEN
!
!           Let FINISH be the position of the next subdiagonal entry
!           such that E( FINISH ) <= TINY or FINISH = N if no such
!           subdiagonal exists.  The matrix identified by the elements
!           between START and FINISH constitutes an independent
!           sub-problem.
!
            FINISH = START
   20       CONTINUE
            IF( FINISH.LT.N ) THEN
               TINY = EPS*SQRT( ABS( D( FINISH ) ) )* SQRT( ABS( D( FINISH+1 ) ) )
               IF( ABS( E( FINISH ) ).GT.TINY ) THEN
                  FINISH = FINISH + 1
                  GO TO 20
               END IF
            END IF
!
!           (Sub) Problem determined.  Compute its size and solve it.
!
            M = FINISH - START + 1
            IF( M.EQ.1 ) THEN
               START = FINISH + 1
               GO TO 10
            END IF
            IF( M.GT.SMLSIZ ) THEN
!
!              Scale.
!
               ORGNRM = DLANST( 'M', M, D( START ), E( START ) )
               nik8 = 1; nk8 = 0; CALL DLASCL( 'G', nk8, nk8, ORGNRM, ONE, M, nik8, D( START ), M, INFO )
               CALL DLASCL( 'G', nk8, nk8, ORGNRM, ONE, M-1, nik8, E( START ), M-1, INFO )
!
               IF( ICOMPZ.EQ.1 ) THEN
                  STRTRW = 1
               ELSE
                  STRTRW = START
               END IF
               CALL DLAED0( ICOMPZ, N, M, D( START ), E( START ),  &
     &                      Z( STRTRW, START ), LDZ, WORK( 1 ), N, &
     &                      WORK( STOREZ ), IWORK, INFO )
               IF( INFO.NE.0 ) THEN
                  INFO = ( INFO / ( M+1 )+START-1 )*( N+1 ) + &
     &                   MOD( INFO, ( M+1 ) ) + START - 1
                  GO TO 50
               END IF
!
!              Scale back.
!
               nik8 = 1; nk8 = 0; CALL DLASCL( 'G', nk8, nk8, ONE, ORGNRM, M, nik8, D( START ), M, INFO )
!
            ELSE
               IF( ICOMPZ.EQ.1 ) THEN
!
!                 Since QR won't update a Z matrix which is larger than
!                 the length of D, we must solve the sub-problem in a
!                 workspace and then multiply back into Z.
!
                  CALL DSTEQR( 'I', M, D( START ), E( START ), WORK, M, WORK( M*M+1 ), INFO )
                  CALL DLACPY( 'A', N, M, Z( 1, START ), LDZ, WORK( STOREZ ), N )
                  CALL DGEMM( 'N', 'N', N, M, M, ONE, WORK( STOREZ ), N, WORK, M, ZERO, Z( 1, START ), LDZ )
               ELSE IF( ICOMPZ.EQ.2 ) THEN
                  CALL DSTEQR( 'I', M, D( START ), E( START ), Z( START, START ), LDZ, WORK, INFO )
               ELSE
                  CALL DSTERF( M, D( START ), E( START ), INFO )
               END IF
               IF( INFO.NE.0 ) THEN
                  INFO = START*( N+1 ) + FINISH
                  GO TO 50
               END IF
            END IF
!
            START = FINISH + 1
            GO TO 10
         END IF
!
!        endwhile
!
!        If the problem split any number of times, then the eigenvalues
!        will not be properly ordered.  Here we permute the eigenvalues
!        (and the associated eigenvectors) into ascending order.
!
         IF( M.NE.N ) THEN
            IF( ICOMPZ.EQ.0 ) THEN
!
!              Use Quick Sort
!
               CALL DLASRT( 'I', N, D, INFO )
!
            ELSE
!
!              Use Selection Sort to minimize swaps of eigenvectors
!
               DO 40 II = 2, N
                  I = II - 1
                  K = I
                  P = D( I )
                  DO 30 J = II, N
                     IF( D( J ).LT.P ) THEN
                        K = J
                        P = D( J )
                     END IF
   30             CONTINUE
                  IF( K.NE.I ) THEN
                     D( K ) = D( I )
                     D( I ) = P
                     nk8 = 1; CALL DSWAP( N, Z( 1, I ), nk8, Z( 1, K ), nk8 )
                  END IF
   40          CONTINUE
            END IF
         END IF
      END IF
!
   50 CONTINUE
      WORK( 1 ) = LWMIN
      IWORK( 1 ) = LIWMIN
!
      RETURN
!
!     End of DSTEDC
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
END Subroutine DSTEDC
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
SUBROUTINE DSTERF( N, D, E, INFO )
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      INTEGER :: INFO, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
!     ..
!
!  Purpose
!  =======
!
!  DSTERF computes all eigenvalues of a symmetric tridiagonal matrix
!  using the Pal-Walker-Kahan variant of the QL or QR algorithm.
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix.  N >= 0.
!
!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the n diagonal elements of the tridiagonal matrix.
!          On exit, if INFO = 0, the eigenvalues in ascending order.
!
!  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, the (n-1) subdiagonal elements of the tridiagonal
!          matrix.
!          On exit, E has been destroyed.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  the algorithm failed to find all of the eigenvalues in
!                a total of 30*N iterations; if INFO = i, then i
!                elements of E have not converged to zero.
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, THREE = 3.0D0 )
      INTEGER :: MAXIT
      PARAMETER          ( MAXIT = 30 )
!     ..
!     .. Local Scalars ..
      INTEGER :: I, ISCALE, JTOT, L, L1, LEND, LENDSV, LSV, M, NMAXIT, nik8, nk8 
      DOUBLE PRECISION   ALPHA, ANORM, BB, C, EPS, EPS2, GAMMA, OLDC,    &
     &                   OLDGAM, P, R, RT1, RT2, RTE, S, SAFMAX, SAFMIN, &
     &                   SIGMA, SSFMAX, SSFMIN
!     ..
!     .. External Functions ..
!     DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
!     EXTERNAL           DLAMCH, DLANST, DLAPY2
!     ..
!     .. External Subroutines ..
!     EXTERNAL           DLAE2, DLASCL, DLASRT, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      INFO = 0
!
!     Quick return if possible
!
      IF( N.LT.0 ) THEN
         INFO = -1
         CALL XERBLA( 'DSTERF', -INFO )
         RETURN
      END IF
      IF( N.LE.1 ) RETURN
!
!     Determine the unit roundoff for this environment.
!
      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
!
!     Compute the eigenvalues of the tridiagonal matrix.
!
      NMAXIT = N*MAXIT
      SIGMA = ZERO
      JTOT = 0
!
!     Determine where the matrix splits and choose QL or QR iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.
!
      L1 = 1
!
   10 CONTINUE
      IF( L1.GT.N ) GO TO 170
      IF( L1.GT.1 ) E( L1-1 ) = ZERO
      DO 20 M = L1, N - 1
         IF( ABS( E( M ) ).LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+1 ) ) ) )*EPS ) THEN
            E( M ) = ZERO
            GO TO 30
         END IF
   20 CONTINUE
      M = N
!
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L ) GO TO 10
!
!     Scale submatrix in rows and columns L to LEND
!
      ANORM = DLANST( 'I', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.GT.SSFMAX ) THEN
         ISCALE = 1
         nik8 = 0; nk8 = 1; CALL DLASCL( 'G', nik8, nik8, ANORM, SSFMAX, LEND-L+1, nk8, D( L ), N, INFO )
         nik8 = 0; nk8 = 1; CALL DLASCL( 'G', nik8, nik8, ANORM, SSFMAX, LEND-L,   nk8, E( L ), N, INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         nik8 = 0; nk8 = 1; CALL DLASCL( 'G', nik8, nik8, ANORM, SSFMIN, LEND-L+1, nk8, D( L ), N, INFO )
         nik8 = 0; nk8 = 1; CALL DLASCL( 'G', nik8, nik8, ANORM, SSFMIN, LEND-L,   nk8, E( L ), N, INFO )
      END IF
!
      DO 40 I = L, LEND - 1
         E( I ) = E( I )**2
   40 CONTINUE
!
!     Choose between QL and QR iteration
!
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
!
      IF( LEND.GE.L ) THEN
!
!        QL Iteration
!
!        Look for small subdiagonal element.
!
   50    CONTINUE
         IF( L.NE.LEND ) THEN
            DO 60 M = L, LEND - 1
               IF( ABS( E( M ) ).LE.EPS2*ABS( D( M )*D( M+1 ) ) ) GO TO 70
   60       CONTINUE
         END IF
         M = LEND
!
   70    CONTINUE
         IF( M.LT.LEND ) E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L ) GO TO 90
!
!        If remaining matrix is 2 by 2, use DLAE2 to compute its
!        eigenvalues.
!
         IF( M.EQ.L+1 ) THEN
            RTE = SQRT( E( L ) )
            CALL DLAE2( D( L ), RTE, D( L+1 ), RT1, RT2 )
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND ) GO TO 50
            GO TO 150
         END IF
!
         IF( JTOT.EQ.NMAXIT ) GO TO 150
         JTOT = JTOT + 1
!
!        Form shift.
!
         RTE = SQRT( E( L ) )
         SIGMA = ( D( L+1 )-P ) / ( TWO*RTE )
         R = DLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
!
         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA
!
!        Inner loop
!
         DO 80 I = M - 1, L, -1
            BB = E( I )
            R = P + BB
            IF( I.NE.M-1 ) E( I+1 ) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I+1 ) = OLDGAM + ( ALPHA-GAMMA )
            IF( C.NE.ZERO ) THEN
               P = ( GAMMA*GAMMA ) / C
            ELSE
               P = OLDC*BB
            END IF
   80    CONTINUE
!
         E( L ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 50
!
!        Eigenvalue found.
!
   90    CONTINUE
         D( L ) = P
!
         L = L + 1
         IF( L.LE.LEND ) GO TO 50
         GO TO 150
!
      ELSE
!
!        QR Iteration
!
!        Look for small superdiagonal element.
!
  100    CONTINUE
         DO 110 M = L, LEND + 1, -1
            IF( ABS( E( M-1 ) ).LE.EPS2*ABS( D( M )*D( M-1 ) ) ) GO TO 120
  110    CONTINUE
         M = LEND
!
  120    CONTINUE
         IF( M.GT.LEND ) E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L ) GO TO 140
!
!        If remaining matrix is 2 by 2, use DLAE2 to compute its
!        eigenvalues.
!
         IF( M.EQ.L-1 ) THEN
            RTE = SQRT( E( L-1 ) )
            CALL DLAE2( D( L ), RTE, D( L-1 ), RT1, RT2 )
            D( L ) = RT1
            D( L-1 ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND ) GO TO 100
            GO TO 150
         END IF
!
         IF( JTOT.EQ.NMAXIT ) GO TO 150
         JTOT = JTOT + 1
!
!        Form shift.
!
         RTE = SQRT( E( L-1 ) )
         SIGMA = ( D( L-1 )-P ) / ( TWO*RTE )
         R = DLAPY2( SIGMA, ONE )
         SIGMA = P - ( RTE / ( SIGMA+SIGN( R, SIGMA ) ) )
!
         C = ONE
         S = ZERO
         GAMMA = D( M ) - SIGMA
         P = GAMMA*GAMMA
!
!        Inner loop
!
         DO 130 I = M, L - 1
            BB = E( I )
            R = P + BB
            IF( I.NE.M ) E( I-1 ) = S*R
            OLDC = C
            C = P / R
            S = BB / R
            OLDGAM = GAMMA
            ALPHA = D( I+1 )
            GAMMA = C*( ALPHA-SIGMA ) - S*OLDGAM
            D( I ) = OLDGAM + ( ALPHA-GAMMA )
            IF( C.NE.ZERO ) THEN
               P = ( GAMMA*GAMMA ) / C
            ELSE
               P = OLDC*BB
            END IF
  130    CONTINUE
!
         E( L-1 ) = S*P
         D( L ) = SIGMA + GAMMA
         GO TO 100
!
!        Eigenvalue found.
!
  140    CONTINUE
         D( L ) = P
!
         L = L - 1
         IF( L.GE.LEND ) GO TO 100
         GO TO 150
!
      END IF
!
!     Undo scaling if necessary
!
150 CONTINUE
        nik8 = 0; nk8 = 1
      IF( ISCALE.EQ.1 ) CALL DLASCL( 'G', nik8, nik8, SSFMAX, ANORM, LENDSV-LSV+1, nk8, D( LSV ), N, INFO )
      IF( ISCALE.EQ.2 ) CALL DLASCL( 'G', nik8, nik8, SSFMIN, ANORM, LENDSV-LSV+1, nk8, D( LSV ), N, INFO )
!
!     Check for no convergence to an eigenvalue after a total
!     of N*MAXIT iterations.
!
      IF( JTOT.LT.NMAXIT ) GO TO 10
      DO 160 I = 1, N - 1
         IF( E( I ).NE.ZERO ) INFO = INFO + 1
  160 CONTINUE
      GO TO 180
!
!     Sort eigenvalues in increasing order.
!
  170 CONTINUE
      CALL DLASRT( 'I', N, D, INFO )
!
  180 CONTINUE
      RETURN
!
!     End of DSTERF
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
END Subroutine DSTERF
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
SUBROUTINE DSYTRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER :: INFO, LDA, LWORK, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DSYTRD reduces a real symmetric matrix A to real symmetric
!  tridiagonal form T by an orthogonal similarity transformation:
!  Q**T * A * Q = T.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          N-by-N upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading N-by-N lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!          On exit, if UPLO = 'U', the diagonal and first superdiagonal
!          of A are overwritten by the corresponding elements of the
!          tridiagonal matrix T, and the elements above the first
!          superdiagonal, with the array TAU, represent the orthogonal
!          matrix Q as a product of elementary reflectors; if UPLO
!          = 'L', the diagonal and first subdiagonal of A are over-
!          written by the corresponding elements of the tridiagonal
!          matrix T, and the elements below the first subdiagonal, with
!          the array TAU, represent the orthogonal matrix Q as a product
!          of elementary reflectors. See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  D       (output) DOUBLE PRECISION array, dimension (N)
!          The diagonal elements of the tridiagonal matrix T:
!          D(i) = A(i,i).
!
!  E       (output) DOUBLE PRECISION array, dimension (N-1)
!          The off-diagonal elements of the tridiagonal matrix T:
!          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
!
!  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  LWORK >= 1.
!          For optimum performance LWORK >= N*NB, where NB is the
!          optimal blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  Further Details
!  ===============
!
!  If UPLO = 'U', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(n-1) . . . H(2) H(1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
!  A(1:i-1,i+1), and tau in TAU(i).
!
!  If UPLO = 'L', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(1) H(2) . . . H(n-1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
!  and tau in TAU(i).
!
!  The contents of A on exit are illustrated by the following examples
!  with n = 5:
!
!  if UPLO = 'U':                       if UPLO = 'L':
!
!    (  d   e   v2  v3  v4 )              (  d                  )
!    (      d   e   v3  v4 )              (  e   d              )
!    (          d   e   v4 )              (  v1  e   d          )
!    (              d   e  )              (  v1  v2  e   d      )
!    (                  d  )              (  v1  v2  v3  e   d  )
!
!  where d and e denote diagonal and off-diagonal elements of T, and vi
!  denotes an element of the vector defining H(i).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER :: I, IINFO, IWS, J, KK, LDWORK, LWKOPT, NB, NBMIN, NX, nik8, nk8
!     ..
!     .. External Subroutines ..
!     EXTERNAL           DLATRD, DSYR2K, DSYTD2, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. External Functions ..
!     LOGICAL            LSAME
!     INTEGER            ILAENV
!     EXTERNAL           LSAME, ILAENV
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -9
      END IF
!
      IF( INFO.EQ.0 ) THEN
!
!        Determine the block size.
!
        nik8 = 1; nk8 = -1
         NB = ILAENV( nik8, 'DSYTRD', UPLO, N, nk8, nk8, nk8 )
         LWKOPT = N*NB
         WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTRD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      NX = N
      IWS = 1
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
!
!        Determine when to cross over from blocked to unblocked code
!        (last block is always handled by unblocked code).
!
          nik8 = 3; nk8 = -1
         NX = MAX( NB, ILAENV( nik8, 'DSYTRD', UPLO, N, nk8, nk8, nk8 ) )
         IF( NX.LT.N ) THEN
!
!           Determine if workspace is large enough for blocked code.
!
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
!
!              Not enough workspace to use optimal NB:  determine the
!              minimum value of NB, and reduce NB or force use of
!              unblocked code by setting NX = N.
!
               NB = MAX( LWORK / LDWORK, 1 )
                nik8 = 2; nk8 = -1
               NBMIN = ILAENV( nik8, 'DSYTRD', UPLO, N, nk8, nk8, nk8 )
               IF( NB.LT.NBMIN ) NX = N
            END IF
         ELSE
            NX = N
         END IF
      ELSE
         NB = 1
      END IF
!
      IF( UPPER ) THEN
!
!        Reduce the upper triangle of A.
!        Columns 1:kk are handled by the unblocked method.
!
         KK = N - ( ( N-NX+NB-1 ) / NB )*NB
         DO 20 I = N - NB + 1, KK + 1, -NB
!
!           Reduce columns i:i+nb-1 to tridiagonal form and form the
!           matrix W which is needed to update the unreduced part of
!           the matrix
!
            CALL DLATRD( UPLO, I+NB-1, NB, A, LDA, E, TAU, WORK, LDWORK )
!
!           Update the unreduced submatrix A(1:i-1,1:i-1), using an
!           update of the form:  A := A - V*W' - W*V'
!
            CALL DSYR2K( UPLO, 'No transpose', I-1, NB, -ONE, A( 1, I ), &
     &                   LDA, WORK, LDWORK, ONE, A, LDA )
!
!           Copy superdiagonal elements back into A, and diagonal
!           elements into D
!
            DO 10 J = I, I + NB - 1
               A( J-1, J ) = E( J-1 )
               D( J ) = A( J, J )
   10       CONTINUE
   20    CONTINUE
!
!        Use unblocked code to reduce the last or only block
!
         CALL DSYTD2( UPLO, KK, A, LDA, D, E, TAU, IINFO )
      ELSE
!
!        Reduce the lower triangle of A
!
         DO 40 I = 1, N - NX, NB
!
!           Reduce columns i:i+nb-1 to tridiagonal form and form the
!           matrix W which is needed to update the unreduced part of
!           the matrix
!
            CALL DLATRD( UPLO, N-I+1, NB, A( I, I ), LDA, E( I ), TAU( I ), WORK, LDWORK )
!
!           Update the unreduced submatrix A(i+ib:n,i+ib:n), using
!           an update of the form:  A := A - V*W' - W*V'
!
            CALL DSYR2K( UPLO, 'No transpose', N-I-NB+1, NB, -ONE,     &
     &                   A( I+NB, I ), LDA, WORK( NB+1 ), LDWORK, ONE, &
     &                   A( I+NB, I+NB ), LDA )
!
!           Copy subdiagonal elements back into A, and diagonal
!           elements into D
!
            DO 30 J = I, I + NB - 1
               A( J+1, J ) = E( J )
               D( J ) = A( J, J )
   30       CONTINUE
   40    CONTINUE
!
!        Use unblocked code to reduce the last or only block
!
         CALL DSYTD2( UPLO, N-I+1, A( I, I ), LDA, D( I ), E( I ), TAU( I ), IINFO )
      END IF
!
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of DSYTRD
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
END Subroutine DSYTRD
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
SUBROUTINE DSYTD2( UPLO, N, A, LDA, D, E, TAU, INFO )
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER :: INFO, LDA, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), D( * ), E( * ), TAU( * )
!     ..
!
!  Purpose
!  =======
!
!  DSYTD2 reduces a real symmetric matrix A to symmetric tridiagonal
!  form T by an orthogonal similarity transformation: Q' * A * Q = T.
!
!  Arguments
!  =========
!
!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored:
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          n-by-n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n-by-n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!          On exit, if UPLO = 'U', the diagonal and first superdiagonal
!          of A are overwritten by the corresponding elements of the
!          tridiagonal matrix T, and the elements above the first
!          superdiagonal, with the array TAU, represent the orthogonal
!          matrix Q as a product of elementary reflectors; if UPLO
!          = 'L', the diagonal and first subdiagonal of A are over-
!          written by the corresponding elements of the tridiagonal
!          matrix T, and the elements below the first subdiagonal, with
!          the array TAU, represent the orthogonal matrix Q as a product
!          of elementary reflectors. See Further Details.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  D       (output) DOUBLE PRECISION array, dimension (N)
!          The diagonal elements of the tridiagonal matrix T:
!          D(i) = A(i,i).
!
!  E       (output) DOUBLE PRECISION array, dimension (N-1)
!          The off-diagonal elements of the tridiagonal matrix T:
!          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
!
!  TAU     (output) DOUBLE PRECISION array, dimension (N-1)
!          The scalar factors of the elementary reflectors (see Further
!          Details).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!
!  Further Details
!  ===============
!
!  If UPLO = 'U', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(n-1) . . . H(2) H(1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
!  A(1:i-1,i+1), and tau in TAU(i).
!
!  If UPLO = 'L', the matrix Q is represented as a product of elementary
!  reflectors
!
!     Q = H(1) H(2) . . . H(n-1).
!
!  Each H(i) has the form
!
!     H(i) = I - tau * v * v'
!
!  where tau is a real scalar, and v is a real vector with
!  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
!  and tau in TAU(i).
!
!  The contents of A on exit are illustrated by the following examples
!  with n = 5:
!
!  if UPLO = 'U':                       if UPLO = 'L':
!
!    (  d   e   v2  v3  v4 )              (  d                  )
!    (      d   e   v3  v4 )              (  e   d              )
!    (          d   e   v4 )              (  v1  e   d          )
!    (              d   e  )              (  v1  v2  e   d      )
!    (                  d  )              (  v1  v2  v3  e   d  )
!
!  where d and e denote diagonal and off-diagonal elements of T, and vi
!  denotes an element of the vector defining H(i).
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO, HALF
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0, HALF = 1.0D0 / 2.0D0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER :: I, nik8
      DOUBLE PRECISION   ALPHA, TAUI
!     ..
!     .. External Subroutines ..
!     EXTERNAL           DAXPY, DLARFG, DSYMV, DSYR2, XERBLA
!     ..
!     .. External Functions ..
!     LOGICAL            LSAME
!     DOUBLE PRECISION   DDOT
!     EXTERNAL           LSAME, DDOT
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYTD2', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( N.LE.0 ) RETURN
!
      IF( UPPER ) THEN
!
!        Reduce the upper triangle of A
!
         DO 10 I = N - 1, 1, -1
!
!           Generate elementary reflector H(i) = I - tau * v * v'
!           to annihilate A(1:i-1,i+1)
!
            nik8 = 1;   CALL DLARFG( I, A( I, I+1 ), A( 1, I+1 ), nik8, TAUI )
            E( I ) = A( I, I+1 )
!
            IF( TAUI.NE.ZERO ) THEN
!
!              Apply H(i) from both sides to A(1:i,1:i)
!
               A( I, I+1 ) = ONE
!
!              Compute  x := tau * A * v  storing x in TAU(1:i)
!
               nik8 = 1;   CALL DSYMV( UPLO, I, TAUI, A, LDA, A( 1, I+1 ), nik8, ZERO, TAU, nik8 )
!
!              Compute  w := x - 1/2 * tau * (x'*v) * v
!
               ALPHA = -HALF*TAUI*DDOT( I, TAU, nik8, A( 1, I+1 ), nik8 )
               CALL DAXPY( I, ALPHA, A( 1, I+1 ), nik8, TAU, nik8 )
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w' - w * v'
!
               CALL DSYR2( UPLO, I, -ONE, A( 1, I+1 ), nik8, TAU, nik8, A, LDA )
!
               A( I, I+1 ) = E( I )
            END IF
            D( I+1 ) = A( I+1, I+1 )
            TAU( I ) = TAUI
   10    CONTINUE
         D( 1 ) = A( 1, 1 )
      ELSE
!
!        Reduce the lower triangle of A
!
         DO 20 I = 1, N - 1
!
!           Generate elementary reflector H(i) = I - tau * v * v'
!           to annihilate A(i+2:n,i)
!
            nik8 = 1;   CALL DLARFG( N-I, A( I+1, I ), A( MIN( I+2, N ), I ), nik8, TAUI )
            E( I ) = A( I+1, I )
!
            IF( TAUI.NE.ZERO ) THEN
!
!              Apply H(i) from both sides to A(i+1:n,i+1:n)
!
               A( I+1, I ) = ONE
!
!              Compute  x := tau * A * v  storing y in TAU(i:n-1)
!
               nik8 = 1;   CALL DSYMV( UPLO, N-I, TAUI, A( I+1, I+1 ), LDA, A( I+1, I ), nik8, ZERO, TAU( I ), nik8 )
!
!              Compute  w := x - 1/2 * tau * (x'*v) * v
!
               nik8 = 1;   ALPHA = -HALF*TAUI*DDOT( N-I, TAU( I ), nik8, A( I+1, I ), nik8 )
               nik8 = 1;   CALL DAXPY( N-I, ALPHA, A( I+1, I ), nik8, TAU( I ), nik8 )
!
!              Apply the transformation as a rank-2 update:
!                 A := A - v * w' - w * v'
!
               nik8 = 1;   CALL DSYR2( UPLO, N-I, -ONE, A( I+1, I ), nik8, TAU( I ), nik8, A( I+1, I+1 ), LDA )
!
               A( I+1, I ) = E( I )
            END IF
            D( I ) = A( I, I )
            TAU( I ) = TAUI
   20    CONTINUE
         D( N ) = A( N, N )
      END IF
!
      RETURN
!
!     End of DSYTD2
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
END Subroutine DSYTD2
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
SUBROUTINE DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, INFO )
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER :: INFO, K, LDA, LDC, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORM2R overwrites the general real m by n matrix C with
!
!        Q * C  if SIDE = 'L' and TRANS = 'N', or
!
!        Q'* C  if SIDE = 'L' and TRANS = 'T', or
!
!        C * Q  if SIDE = 'R' and TRANS = 'N', or
!
!        C * Q' if SIDE = 'R' and TRANS = 'T',
!
!  where Q is a real orthogonal matrix defined as the product of k
!  elementary reflectors
!
!        Q = H(1) H(2) . . . H(k)
!
!  as returned by DGEQRF. Q is of order m if SIDE = 'L' and of order n
!  if SIDE = 'R'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q' from the Left
!          = 'R': apply Q or Q' from the Right
!
!  TRANS   (input) CHARACTER*1
!          = 'N': apply Q  (No transpose)
!          = 'T': apply Q' (Transpose)
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
!          The i-th column must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          DGEQRF in the first k columns of its array argument A.
!          A is modified by the routine but restored on exit.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!          If SIDE = 'L', LDA >= max(1,M);
!          if SIDE = 'R', LDA >= max(1,N).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace) DOUBLE PRECISION array, dimension
!                                   (N) if SIDE = 'L',
!                                   (M) if SIDE = 'R'
!
!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER :: I, I1, I2, I3, IC, JC, MI, NI, NQ, nik8
      DOUBLE PRECISION   AII
!     ..
!     .. External Functions ..
!     LOGICAL            LSAME
!     EXTERNAL           LSAME
!     ..
!     .. External Subroutines ..
!     EXTERNAL           DLARF, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
!
!     NQ is the order of Q
!
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORM2R', -INFO )
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) RETURN
!
      IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) ) THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
!
      IF( LEFT ) THEN
         NI = N
         JC = 1
      ELSE
         MI = M
         IC = 1
      END IF
!
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
!
!           H(i) is applied to C(i:m,1:n)
!
            MI = M - I + 1
            IC = I
         ELSE
!
!           H(i) is applied to C(1:m,i:n)
!
            NI = N - I + 1
            JC = I
         END IF
!
!        Apply H(i)
!
         AII = A( I, I )
         A( I, I ) = ONE
         nik8 = 1; CALL DLARF( SIDE, MI, NI, A( I, I ), nik8, TAU( I ), C( IC, JC ), LDC, WORK )
         A( I, I ) = AII
   10 CONTINUE
      RETURN
!
!     End of DORM2R
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
END Subroutine DORM2R
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
SUBROUTINE DORMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, LWORK, INFO )
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!  -- LAPACK routine (version 3.1) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2006
!
!     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER :: INFO, K, LDA, LDC, LWORK, M, N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!     ..
!
!  Purpose
!  =======
!
!  DORMQR overwrites the general real M-by-N matrix C with
!
!                  SIDE = 'L'     SIDE = 'R'
!  TRANS = 'N':      Q * C          C * Q
!  TRANS = 'T':      Q**T * C       C * Q**T
!
!  where Q is a real orthogonal matrix defined as the product of k
!  elementary reflectors
!
!        Q = H(1) H(2) . . . H(k)
!
!  as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N
!  if SIDE = 'R'.
!
!  Arguments
!  =========
!
!  SIDE    (input) CHARACTER*1
!          = 'L': apply Q or Q**T from the Left;
!          = 'R': apply Q or Q**T from the Right.
!
!  TRANS   (input) CHARACTER*1
!          = 'N':  No transpose, apply Q;
!          = 'T':  Transpose, apply Q**T.
!
!  M       (input) INTEGER
!          The number of rows of the matrix C. M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix C. N >= 0.
!
!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines
!          the matrix Q.
!          If SIDE = 'L', M >= K >= 0;
!          if SIDE = 'R', N >= K >= 0.
!
!  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
!          The i-th column must contain the vector which defines the
!          elementary reflector H(i), for i = 1,2,...,k, as returned by
!          DGEQRF in the first k columns of its array argument A.
!          A is modified by the routine but restored on exit.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.
!          If SIDE = 'L', LDA >= max(1,M);
!          if SIDE = 'R', LDA >= max(1,N).
!
!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.
!
!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the M-by-N matrix C.
!          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
!
!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.
!          If SIDE = 'L', LWORK >= max(1,N);
!          if SIDE = 'R', LWORK >= max(1,M).
!          For optimum performance LWORK >= N*NB if SIDE = 'L', and
!          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
!          blocksize.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER :: NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
!     ..
!     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      INTEGER :: I, I1, I2, I3, IB, IC, IINFO, IWS, JC, LDWORK, LWKOPT, MI, NB, NBMIN, NI, NQ, NW, nik8, nk8
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION   T( LDT, NBMAX )
!     ..
!     .. External Functions ..
!     LOGICAL            LSAME
!     INTEGER            ILAENV
!     EXTERNAL           LSAME, ILAENV
!     ..
!     .. External Subroutines ..
!     EXTERNAL           DLARFB, DLARFT, DORM2R, XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
!
!     NQ is the order of Q and NW is the minimum dimension of WORK
!
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
!
      IF( INFO.EQ.0 ) THEN
!
!        Determine the block size.  NB may be at most NBMAX, where NBMAX
!        is used to define the local array T.
!
         nik8 = 1; NB = MIN( NBMAX, ILAENV( nik8, 'DORMQR', SIDE // TRANS, M, N, K, -nik8 ) )
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      END IF
!
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORMQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
!
!     Quick return if possible
!
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
!
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IWS = NW*NB
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            nik8 = 2; nk8 = -1; NBMIN = MAX( 2, ILAENV( nik8, 'DORMQR', SIDE // TRANS, M, N, K, nk8 ) )
         END IF
      ELSE
         IWS = NW
      END IF
!
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
!
!        Use unblocked code
!
         CALL DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK, IINFO )
      ELSE
!
!        Use blocked code
!
         IF( ( LEFT .AND. .NOT.NOTRAN ) .OR. ( .NOT.LEFT .AND. NOTRAN ) ) THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         END IF
!
         IF( LEFT ) THEN
            NI = N
            JC = 1
         ELSE
            MI = M
            IC = 1
         END IF
!
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
!
!           Form the triangular factor of the block reflector
!           H = H(i) H(i+1) . . . H(i+ib-1)
!
            CALL DLARFT( 'Forward', 'Columnwise', NQ-I+1, IB, A( I, I ), LDA, TAU( I ), T, LDT )
            IF( LEFT ) THEN
!
!              H or H' is applied to C(i:m,1:n)
!
               MI = M - I + 1
               IC = I
            ELSE
!
!              H or H' is applied to C(1:m,i:n)
!
               NI = N - I + 1
               JC = I
            END IF
!
!           Apply H or H'
!
            CALL DLARFB( SIDE, TRANS, 'Forward', 'Columnwise', MI, NI, &
     &                   IB, A( I, I ), LDA, T, LDT, C( IC, JC ), LDC, &
     &                   WORK, LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
!
!     End of DORMQR
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
END Subroutine DORMQR
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!*************************************************************************************************!
!
End Module Eigenlib
!
!*************************************************************************************************!

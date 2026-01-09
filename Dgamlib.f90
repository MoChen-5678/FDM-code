!*****************************************************************************!
!
Module Dgamlib
!
!*****************************************************************************!
!
use DImach
use DSlib
USE ERRORS
implicit none

contains
!=============================================================================!
!
!*DECK DGAMMA
FUNCTION DGAMMA (X)
!***BEGIN PROLOGUE  DGAMMA
!***PURPOSE  Compute the complete Gamma function.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A
!***TYPE      DOUBLE PRECISION (GAMMA-S, DGAMMA-D, CGAMMA-!)
!***KEYWORDS  COMPLETE GAMMA FUNCTION, FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DGAMMA(X) calculates the double precision complete Gamma function
! for double precision argument X.
!
! Series for GAM        on the interval  0.          to  1.00000E+00
!                                        with weighted error   5.79E-32
!                                         log weighted error  31.24
!                               significant figures required  30.00
!                                    decimal places required  32.05
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, D9LGMC, DCSEVL, DGAMLM, INITDS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920618  Removed space from variable name.  (RWC, WRB)
!***END PROLOGUE  DGAMMA
      DOUBLE PRECISION X, GAMCS(42), DXREL, PI, SINPIY, SQ2PIL, XMAX, &
     &  XMIN, Y, DGAMMA !, D9LGMC, DCSEVL, D1MACH
      LOGICAL FIRST
      integer :: NGAM, I, N
!
      SAVE GAMCS, PI, SQ2PIL, NGAM, XMIN, XMAX, DXREL, FIRST
      DATA GAMCS(  1) / +.8571195590989331421920062399942D-2      /
      DATA GAMCS(  2) / +.4415381324841006757191315771652D-2      /
      DATA GAMCS(  3) / +.5685043681599363378632664588789D-1      /
      DATA GAMCS(  4) / -.4219835396418560501012500186624D-2      /
      DATA GAMCS(  5) / +.1326808181212460220584006796352D-2      /
      DATA GAMCS(  6) / -.1893024529798880432523947023886D-3      /
      DATA GAMCS(  7) / +.3606925327441245256578082217225D-4      /
      DATA GAMCS(  8) / -.6056761904460864218485548290365D-5      /
      DATA GAMCS(  9) / +.1055829546302283344731823509093D-5      /
      DATA GAMCS( 10) / -.1811967365542384048291855891166D-6      /
      DATA GAMCS( 11) / +.3117724964715322277790254593169D-7      /
      DATA GAMCS( 12) / -.5354219639019687140874081024347D-8      /
      DATA GAMCS( 13) / +.9193275519859588946887786825940D-9      /
      DATA GAMCS( 14) / -.1577941280288339761767423273953D-9      /
      DATA GAMCS( 15) / +.2707980622934954543266540433089D-10     /
      DATA GAMCS( 16) / -.4646818653825730144081661058933D-11     /
      DATA GAMCS( 17) / +.7973350192007419656460767175359D-12     /
      DATA GAMCS( 18) / -.1368078209830916025799499172309D-12     /
      DATA GAMCS( 19) / +.2347319486563800657233471771688D-13     /
      DATA GAMCS( 20) / -.4027432614949066932766570534699D-14     /
      DATA GAMCS( 21) / +.6910051747372100912138336975257D-15     /
      DATA GAMCS( 22) / -.1185584500221992907052387126192D-15     /
      DATA GAMCS( 23) / +.2034148542496373955201026051932D-16     /
      DATA GAMCS( 24) / -.3490054341717405849274012949108D-17     /
      DATA GAMCS( 25) / +.5987993856485305567135051066026D-18     /
      DATA GAMCS( 26) / -.1027378057872228074490069778431D-18     /
      DATA GAMCS( 27) / +.1762702816060529824942759660748D-19     /
      DATA GAMCS( 28) / -.3024320653735306260958772112042D-20     /
      DATA GAMCS( 29) / +.5188914660218397839717833550506D-21     /
      DATA GAMCS( 30) / -.8902770842456576692449251601066D-22     /
      DATA GAMCS( 31) / +.1527474068493342602274596891306D-22     /
      DATA GAMCS( 32) / -.2620731256187362900257328332799D-23     /
      DATA GAMCS( 33) / +.4496464047830538670331046570666D-24     /
      DATA GAMCS( 34) / -.7714712731336877911703901525333D-25     /
      DATA GAMCS( 35) / +.1323635453126044036486572714666D-25     /
      DATA GAMCS( 36) / -.2270999412942928816702313813333D-26     /
      DATA GAMCS( 37) / +.3896418998003991449320816639999D-27     /
      DATA GAMCS( 38) / -.6685198115125953327792127999999D-28     /
      DATA GAMCS( 39) / +.1146998663140024384347613866666D-28     /
      DATA GAMCS( 40) / -.1967938586345134677295103999999D-29     /
      DATA GAMCS( 41) / +.3376448816585338090334890666666D-30     /
      DATA GAMCS( 42) / -.5793070335782135784625493333333D-31     /
      DATA PI / 3.14159265358979323846264338327950D0 /
      DATA SQ2PIL / 0.91893853320467274178032973640562D0 /
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DGAMMA
      IF (FIRST) THEN
         NGAM = INITDS (GAMCS, 42, 0.1*REAL(D1MACH(3)) )
!
         CALL DGAMLM (XMIN, XMAX)
         DXREL = SQRT(D1MACH(4))
      ENDIF
      FIRST = .FALSE.
!
      Y = ABS(X)
      IF (Y.GT.10.D0) GO TO 50
!
! COMPUTE GAMMA(X) FOR -XBND .LE. X .LE. XBND.  REDUCE INTERVAL AND FIND
! GAMMA(1+Y) FOR 0.0 .LE. Y .LT. 1.0 FIRST OF ALL.
!
      N = X
      IF (X.LT.0.D0) N = N - 1
      Y = X - N
      N = N - 1
      DGAMMA = 0.9375D0 + DCSEVL (2.D0*Y-1.D0, GAMCS, NGAM)
      IF (N.EQ.0) RETURN
!
      IF (N.GT.0) GO TO 30
!
! COMPUTE GAMMA(X) FOR X .LT. 1.0
!
      N = -N
      IF (X .EQ. 0.D0) CALL XERMSG ('SLATEC', 'DGAMMA', 'X IS 0', 4, 2)
      IF (X .LT. 0.0 .AND. X+N-2 .EQ. 0.D0) CALL XERMSG ('SLATEC', &
     &   'DGAMMA', 'X IS A NEGATIVE INTEGER', 4, 2)
      IF (X .LT. (-0.5D0) .AND. ABS((X-AINT(X-0.5D0))/X) .LT. DXREL) &
     &   CALL XERMSG ('SLATEC', 'DGAMMA',  &
     &   'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER',&
     &   1, 1)
!
      DO 20 I=1,N
        DGAMMA = DGAMMA/(X+I-1 )
 20   CONTINUE
      RETURN
!
! GAMMA(X) FOR X .GE. 2.0 AND X .LE. 10.0
!
 30   DO 40 I=1,N
        DGAMMA = (Y+I) * DGAMMA
 40   CONTINUE
      RETURN
!
! GAMMA(X) FOR ABS(X) .GT. 10.0.  RECALL Y = ABS(X).
!
 50   IF (X .GT. XMAX) CALL XERMSG ('SLATEC', 'DGAMMA','X SO BIG GAMMA OVERFLOWS', 3, 2)
!
      DGAMMA = 0.D0
      IF (X .LT. XMIN) CALL XERMSG ('SLATEC', 'DGAMMA', 'X SO SMALL GAMMA UNDERFLOWS', 2, 1)
      IF (X.LT.XMIN) RETURN
!
      DGAMMA = EXP ((Y-0.5D0)*LOG(Y) - Y + SQ2PIL + D9LGMC(Y) )
      IF (X.GT.0.D0) RETURN
!
      IF (ABS((X-AINT(X-0.5D0))/X) .LT. DXREL) CALL XERMSG ('SLATEC', &
     &   'DGAMMA', 'ANSWER LT HALF PRECISION, X TOO NEAR NEGATIVE INTEGER', 1, 1)
!
      SINPIY = SIN (PI*Y)
      IF (SINPIY .EQ. 0.D0) CALL XERMSG ('SLATEC', 'DGAMMA', 'X IS A NEGATIVE INTEGER', 4, 2)
!
      DGAMMA = -PI/(Y*SINPIY*DGAMMA)
!
      RETURN
END function DGAMMA

!=============================================================================!
!
!DECK D9LGMC
FUNCTION D9LGMC (X)
!***BEGIN PROLOGUE  D9LGMC
!***SUBSIDIARY
!***PURPOSE  Compute the log Gamma correction factor so that
!            LOG(DGAMMA(X)) = LOG(SQRT(2*PI)) + (X-5.)*LOG(X) - X
!            + D9LGMC(X).
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7E
!***TYPE      DOUBLE PRECISION (R9LGMC-S, D9LGMC-D, C9LGMC-!)
!***KEYWORDS  COMPLETE GAMMA FUNCTION, CORRECTION TERM, FNLIB,
!             LOG GAMMA, LOGARITHM, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Compute the log gamma correction factor for X .GE. 10. so that
! LOG (DGAMMA(X)) = LOG(SQRT(2*PI)) + (X-.5)*LOG(X) - X + D9lGMC(X)
!
! Series for ALGM       on the interval  0.          to  1.00000E-02
!                                        with weighted error   1.28E-31
!                                         log weighted error  30.89
!                               significant figures required  29.81
!                                    decimal places required  31.48
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900720  Routine changed from user-callable to subsidiary.  (WRB)
!***END PROLOGUE  D9LGMC
      DOUBLE PRECISION X, ALGMCS(15), XBIG, XMAX, D9LGMC !, DCSEVL, D1MACH
      integer :: NALGM
      LOGICAL FIRST
      SAVE ALGMCS, NALGM, XBIG, XMAX, FIRST
      DATA ALGMCS(  1) / +.1666389480451863247205729650822D+0      /
      DATA ALGMCS(  2) / -.1384948176067563840732986059135D-4      /
      DATA ALGMCS(  3) / +.9810825646924729426157171547487D-8      /
      DATA ALGMCS(  4) / -.1809129475572494194263306266719D-10     /
      DATA ALGMCS(  5) / +.6221098041892605227126015543416D-13     /
      DATA ALGMCS(  6) / -.3399615005417721944303330599666D-15     /
      DATA ALGMCS(  7) / +.2683181998482698748957538846666D-17     /
      DATA ALGMCS(  8) / -.2868042435334643284144622399999D-19     /
      DATA ALGMCS(  9) / +.3962837061046434803679306666666D-21     /
      DATA ALGMCS( 10) / -.6831888753985766870111999999999D-23     /
      DATA ALGMCS( 11) / +.1429227355942498147573333333333D-24     /
      DATA ALGMCS( 12) / -.3547598158101070547199999999999D-26     /
      DATA ALGMCS( 13) / +.1025680058010470912000000000000D-27     /
      DATA ALGMCS( 14) / -.3401102254316748799999999999999D-29     /
      DATA ALGMCS( 15) / +.1276642195630062933333333333333D-30     /
      DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  D9LGMC
      IF (FIRST) THEN
         NALGM = INITDS (ALGMCS, 15, REAL(D1MACH(3)) )
         XBIG = 1.0D0/SQRT(D1MACH(3))
         XMAX = EXP (MIN(LOG(D1MACH(2)/12.D0), -LOG(12.D0*D1MACH(1))))
      ENDIF
      FIRST = .FALSE.
!
      IF (X .LT. 10.D0) CALL XERMSG ('SLATEC', 'D9LGMC', 'X MUST BE GE 10', 1, 2)
      IF (X.GE.XMAX) GO TO 20
!
      D9LGMC = 1.D0/(12.D0*X)
      IF (X.LT.XBIG) D9LGMC = DCSEVL (2.0D0*(10.D0/X)**2-1.D0, ALGMCS, NALGM) / X
      RETURN
!
 20   D9LGMC = 0.D0
      CALL XERMSG ('SLATEC', 'D9LGMC', 'X SO BIG D9LGMC UNDERFLOWS', 2, 1)
      RETURN
!
END function D9LGMC


!=============================================================================!
!
!DECK DGAMLM
SUBROUTINE DGAMLM (XMIN, XMAX)
!***BEGIN PROLOGUE  DGAMLM
!***PURPOSE  Compute the minimum and maximum bounds for the argument in
!            the Gamma function.
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A, R2
!***TYPE      DOUBLE PRECISION (GAMLIM-S, DGAMLM-D)
!***KEYWORDS  COMPLETE GAMMA FUNCTION, FNLIB, LIMITS, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Calculate the minimum and maximum legal bounds for X in gamma(X).
! XMIN and XMAX are not the only bounds, but they are the only non-
! trivial ones to calculate.
!
!             Output Arguments --
! XMIN   double precision minimum legal value of X in gamma(X).  Any
!        smaller value of X might result in underflow.
! XMAX   double precision maximum legal value of X in gamma(X).  Any
!        larger value of X might cause overflow.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DGAMLM
      DOUBLE PRECISION XMIN, XMAX, ALNBIG, ALNSML, XLN, XOLD !, D1MACH
      Integer :: I
!***FIRST EXECUTABLE STATEMENT  DGAMLM
      ALNSML = LOG(D1MACH(1))
      XMIN = -ALNSML
      DO 10 I=1,10
        XOLD = XMIN
        XLN = LOG(XMIN)
        XMIN = XMIN - XMIN*((XMIN+0.5D0)*XLN - XMIN - 0.2258D0 + ALNSML)/(XMIN*XLN+0.5D0)
        IF (ABS(XMIN-XOLD).LT.0.005D0) GO TO 20
 10   CONTINUE
      CALL XERMSG ('SLATEC', 'DGAMLM', 'UNABLE TO FIND XMIN', 1, 2)
!
 20   XMIN = -XMIN + 0.01D0
!
      ALNBIG = LOG (D1MACH(2))
      XMAX = ALNBIG
      DO 30 I=1,10
        XOLD = XMAX
        XLN = LOG(XMAX)
        XMAX = XMAX - XMAX*((XMAX-0.5D0)*XLN - XMAX + 0.9189D0 - ALNBIG)/(XMAX*XLN-0.5D0)
        IF (ABS(XMAX-XOLD).LT.0.005D0) GO TO 40
 30   CONTINUE
      CALL XERMSG ('SLATEC', 'DGAMLM', 'UNABLE TO FIND XMAX', 2, 2)
!
 40   XMAX = XMAX - 0.01D0
      XMIN = MAX (XMIN, -XMAX+1.D0)
!
      RETURN
END subroutine DGAMLM

end module Dgamlib
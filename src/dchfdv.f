*DECK DCHFDV
      SUBROUTINE DCHFDV (X1, X2, F1, F2, D1, D2, NE, XE, FE, DE, NEXT,
     +   IERR)
C***BEGIN PROLOGUE  DCHFDV
C***PURPOSE  Evaluate a cubic polynomial given in Hermite form and its
C            first derivative at an array of points.  While designed for
C            use by DPCHFD, it may be useful directly as an evaluator
C            for a piecewise cubic Hermite function in applications,
C            such as graphing, where the interval is known in advance.
C            If only function values are required, use DCHFEV instead.
C***LIBRARY   SLATEC (PCHIP)
C***CATEGORY  E3, H1
C***TYPE      DOUBLE PRECISION (CHFDV-S, DCHFDV-D)
C***KEYWORDS  CUBIC HERMITE DIFFERENTIATION, CUBIC HERMITE EVALUATION,
C             CUBIC POLYNOMIAL EVALUATION, PCHIP
C***AUTHOR  Fritsch, F. N., (LLNL)
C             Lawrence Livermore National Laboratory
C             P.O. Box 808  (L-316)
C             Livermore, CA  94550
C             FTS 532-4275, (510) 422-4275
C***DESCRIPTION
C
C        DCHFDV:  Cubic Hermite Function and Derivative Evaluator
C
C     Evaluates the cubic polynomial determined by function values
C     F1,F2 and derivatives D1,D2 on interval (X1,X2), together with
C     its first derivative, at the points  XE(J), J=1(1)NE.
C
C     If only function values are required, use DCHFEV, instead.
C
C ----------------------------------------------------------------------
C
C  Calling sequence:
C
C        INTEGER  NE, NEXT(2), IERR
C        DOUBLE PRECISION  X1, X2, F1, F2, D1, D2, XE(NE), FE(NE),
C                          DE(NE)
C
C        CALL  DCHFDV (X1,X2, F1,F2, D1,D2, NE, XE, FE, DE, NEXT, IERR)
C
C   Parameters:
C
C     X1,X2 -- (input) endpoints of interval of definition of cubic.
C           (Error return if  X1.EQ.X2 .)
C
C     F1,F2 -- (input) values of function at X1 and X2, respectively.
C
C     D1,D2 -- (input) values of derivative at X1 and X2, respectively.
C
C     NE -- (input) number of evaluation points.  (Error return if
C           NE.LT.1 .)
C
C     XE -- (input) real*8 array of points at which the functions are to
C           be evaluated.  If any of the XE are outside the interval
C           [X1,X2], a warning error is returned in NEXT.
C
C     FE -- (output) real*8 array of values of the cubic function
C           defined by  X1,X2, F1,F2, D1,D2  at the points  XE.
C
C     DE -- (output) real*8 array of values of the first derivative of
C           the same function at the points  XE.
C
C     NEXT -- (output) integer array indicating number of extrapolation
C           points:
C            NEXT(1) = number of evaluation points to left of interval.
C            NEXT(2) = number of evaluation points to right of interval.
C
C     IERR -- (output) error flag.
C           Normal return:
C              IERR = 0  (no errors).
C           "Recoverable" errors:
C              IERR = -1  if NE.LT.1 .
C              IERR = -2  if X1.EQ.X2 .
C                (Output arrays have not been changed in either case.)
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   811019  DATE WRITTEN
C   820803  Minor cosmetic changes for release 1.
C   870707  Corrected XERROR calls for d.p. names(s).
C   870813  Minor cosmetic changes.
C   890411  Added SAVE statements (Vers. 3.2).
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   891006  Cosmetic changes to prologue.  (WRB)
C   891006  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  DCHFDV
C  Programming notes:
C
C     To produce a single precision version, simply:
C        a. Change DCHFDV to CHFDV wherever it occurs,
C        b. Change the double precision declaration to real, and
C        c. Change the constant ZERO to single precision.
C
C  DECLARE ARGUMENTS.
C
      INTEGER  NE, NEXT(2), IERR
      DOUBLE PRECISION  X1, X2, F1, F2, D1, D2, XE(*), FE(*), DE(*)
C
C  DECLARE LOCAL VARIABLES.
C
      INTEGER  I
      DOUBLE PRECISION  C2, C2T2, C3, C3T3, DEL1, DEL2, DELTA, H, X,
     *  XMI, XMA, ZERO
      SAVE ZERO
      DATA  ZERO /0.D0/
C
C  VALIDITY-CHECK ARGUMENTS.
C
C***FIRST EXECUTABLE STATEMENT  DCHFDV
      IF (NE .LT. 1)  GO TO 5001
      H = X2 - X1
      IF (H .EQ. ZERO)  GO TO 5002
C
C  INITIALIZE.
C
      IERR = 0
      NEXT(1) = 0
      NEXT(2) = 0
      XMI = MIN(ZERO, H)
      XMA = MAX(ZERO, H)
C
C  COMPUTE CUBIC COEFFICIENTS (EXPANDED ABOUT X1).
C
      DELTA = (F2 - F1)/H
      DEL1 = (D1 - DELTA)/H
      DEL2 = (D2 - DELTA)/H
C                                           (DELTA IS NO LONGER NEEDED.)
      C2 = -(DEL1+DEL1 + DEL2)
      C2T2 = C2 + C2
      C3 = (DEL1 + DEL2)/H
C                               (H, DEL1 AND DEL2 ARE NO LONGER NEEDED.)
      C3T3 = C3+C3+C3
C
C  EVALUATION LOOP.
C
      DO 500  I = 1, NE
         X = XE(I) - X1
         FE(I) = F1 + X*(D1 + X*(C2 + X*C3))
         DE(I) = D1 + X*(C2T2 + X*C3T3)
C          COUNT EXTRAPOLATION POINTS.
         IF ( X.LT.XMI )  NEXT(1) = NEXT(1) + 1
         IF ( X.GT.XMA )  NEXT(2) = NEXT(2) + 1
C        (NOTE REDUNDANCY--IF EITHER CONDITION IS TRUE, OTHER IS FALSE.)
  500 CONTINUE
C
C  NORMAL RETURN.
C
      RETURN
C
C  ERROR RETURNS.
C
 5001 CONTINUE
C     NE.LT.1 RETURN.
      IERR = -1
      CALL XERMSG ('SLATEC', 'DCHFDV',
     +   'NUMBER OF EVALUATION POINTS LESS THAN ONE', IERR, 1)
      RETURN
C
 5002 CONTINUE
C     X1.EQ.X2 RETURN.
      IERR = -2
      CALL XERMSG ('SLATEC', 'DCHFDV', 'INTERVAL ENDPOINTS EQUAL',
     +   IERR, 1)
      RETURN
C------------- LAST LINE OF DCHFDV FOLLOWS -----------------------------
      END

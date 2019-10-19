!*******************************************************************************************************
!>
!
   module pchip_module

   use, intrinsic :: iso_fortran_env, only: wp => real64

   implicit none

   real(wp),parameter :: zero   = 0.0_wp
   real(wp),parameter :: half   = 0.5_wp
   real(wp),parameter :: one    = 1.0_wp
   real(wp),parameter :: two    = 2.0_wp
   real(wp),parameter :: three  = 3.0_wp
   real(wp),parameter :: four   = 4.0_wp
   real(wp),parameter :: six    = 6.0_wp
   real(wp),parameter :: ten    = 10.0_wp

   real(wp),parameter :: d1mach4 = epsilon(one)   !! d1mach(4) -- the largest relative spacing

   contains
!*******************************************************************************************************

   pure function dpchsd(s1,s2,h1,h2)
   !! inline function for weighted average of slopes.
   real(wp),intent(in) :: s1, s2, h1, h2
   real(wp) :: dpchsd
   dpchsd = (h2/(h1+h2))*s1 + (h1/(h1+h2))*s2
   end function dpchsd



!***PURPOSE  Check a single cubic for monotonicity.
!***LIBRARY   SLATEC (PCHIP)
!***TYPE      real(wp) (CHFCM-S, DCHFCM-D)
!***AUTHOR  Fritsch, F. N., (LLNL)
!***DESCRIPTION
!
! *Usage:
!
!        real(wp)  D1, D2, DELTA
!        INTEGER  ISMON, DCHFCM
!
!        ISMON = DCHFCM (D1, D2, DELTA)
!
! *Arguments:
!
!     D1,D2:IN  are the derivative values at the ends of an interval.
!
!     DELTA:IN  is the data slope over that interval.
!
! *Function Return Values:
!     ISMON : indicates the monotonicity of the cubic segment:
!             ISMON = -3  if function is probably decreasing;
!             ISMON = -1  if function is strictly decreasing;
!             ISMON =  0  if function is constant;
!             ISMON =  1  if function is strictly increasing;
!             ISMON =  2  if function is non-monotonic;
!             ISMON =  3  if function is probably increasing.
!           If ABS(ISMON)=3, the derivative values are too close to the
!           boundary of the monotonicity region to declare monotonicity
!           in the presence of roundoff error.
!
! *Description:
!
!          DCHFCM:  Cubic Hermite Function -- Check Monotonicity.
!
!    Called by  DPCHCM  to determine the monotonicity properties of the
!    cubic with boundary derivative values D1,D2 and chord slope DELTA.
!
! *Cautions:
!     This is essentially the same as old DCHFMC, except that a
!     new output value, -3, was added February 1989.  (Formerly, -3
!     and +3 were lumped together in the single value 3.)  Codes that
!     flag nonmonotonicity by "IF (ISMON==2)" need not be changed.
!     Codes that check via "IF (ISMON>=3)" should change the test to
!     "IF (IABS(ISMON)>=3)".  Codes that declare monotonicity via
!     "IF (ISMON<=1)" should change to "IF (IABS(ISMON)<=1)".
!
!   REFER TO  DPCHCM
!
!***ROUTINES CALLED  D1MACH
!***REVISION HISTORY  (YYMMDD)
!   820518  DATE WRITTEN
!   820805  Converted to SLATEC library version.
!   831201  Changed from  ISIGN  to SIGN  to correct bug that
!           produced wrong sign when -1 < DELTA < 0 .
!   890206  Added SAVE statements.
!   890209  Added sign to returned value ISMON=3 and corrected
!           argument description accordingly.
!   890306  Added caution about changed output.
!   890407  Changed name from DCHFMC to DCHFCM, as requested at the
!           March 1989 SLATEC CML meeting, and made a few other
!           minor modifications necessitated by this change.
!   890407  Converted to new SLATEC format.
!   890407  Modified DESCRIPTION to LDOC format.
!   891214  Moved SAVE statements.  (WRB)

   integer function dchfcm (d1, d2, delta)

   real(wp) d1, d2, delta

   integer ismon, itrue
   real(wp)  a, b, phi

   real(wp),parameter :: eps = ten*d1mach4  !! machine-dependent parameter -- should be about 10*uround.
                                          !! TEN is actually a tuning parameter, which determines the width of
                                          !! the fuzz around the elliptical boundary.

   !  MAKE THE CHECK.

   if (delta == zero)  then
   !        CASE OF CONSTANT DATA.
      if ((d1==zero) .and. (d2==zero))  then
         ismon = 0
      else
         ismon = 2
      endif
   else
   !        DATA IS NOT CONSTANT -- PICK UP SIGN.
      itrue = sign (one, delta)
      a = d1/delta
      b = d2/delta
      if ((a<zero) .or. (b<zero))  then
         ismon = 2
      else if ((a<=three-eps) .and. (b<=three-eps))  then
   !           INSIDE SQUARE (0,3)X(0,3)  IMPLIES   OK.
         ismon = itrue
      else if ((a>four+eps) .and. (b>four+eps))  then
   !           OUTSIDE SQUARE (0,4)X(0,4)  IMPLIES   NONMONOTONIC.
         ismon = 2
      else
   !           MUST CHECK AGAINST BOUNDARY OF ELLIPSE.
         a = a - two
         b = b - two
         phi = ((a*a + b*b) + a*b) - three
         if (phi < -eps)  then
            ismon = itrue
         else if (phi > eps)  then
            ismon = 2
         else
   !              TO CLOSE TO BOUNDARY TO TELL,
   !                  IN THE PRESENCE OF ROUND-OFF ERRORS.
            ismon = 3*itrue
         endif
      endif
   endif

   dchfcm = ismon

   end function dchfcm

subroutine dchfdv (x1, x2, f1, f2, d1, d2, ne, xe, fe, de, next, ierr)

!***PURPOSE  Evaluate a cubic polynomial given in Hermite form and its
!            first derivative at an array of points.  While designed for
!            use by DPCHFD, it may be useful directly as an evaluator
!            for a piecewise cubic Hermite function in applications,
!            such as graphing, where the interval is known in advance.
!            If only function values are required, use DCHFEV instead.
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E3, H1
!***TYPE      real(wp) (CHFDV-S, DCHFDV-D)
!***KEYWORDS  CUBIC HERMITE DIFFERENTIATION, CUBIC HERMITE EVALUATION,
!             CUBIC POLYNOMIAL EVALUATION, PCHIP
!***AUTHOR  Fritsch, F. N., (LLNL)
!             Lawrence Livermore National Laboratory
!             P.O. Box 808  (L-316)
!             Livermore, CA  94550
!             FTS 532-4275, (510) 422-4275
!***DESCRIPTION
!
!        DCHFDV:  Cubic Hermite Function and Derivative Evaluator
!
!     Evaluates the cubic polynomial determined by function values
!     F1,F2 and derivatives D1,D2 on interval (X1,X2), together with
!     its first derivative, at the points  XE(J), J=1(1)NE.
!
!     If only function values are required, use DCHFEV, instead.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        INTEGER  NE, NEXT(2), IERR
!        real(wp)  X1, X2, F1, F2, D1, D2, XE(NE), FE(NE),
!                          DE(NE)
!
!        CALL  DCHFDV (X1,X2, F1,F2, D1,D2, NE, XE, FE, DE, NEXT, IERR)
!
!   Parameters:
!
!     X1,X2 -- (input) endpoints of interval of definition of cubic.
!           (Error return if  X1==X2 .)
!
!     F1,F2 -- (input) values of function at X1 and X2, respectively.
!
!     D1,D2 -- (input) values of derivative at X1 and X2, respectively.
!
!     NE -- (input) number of evaluation points.  (Error return if
!           NE<1 .)
!
!     XE -- (input) real*8 array of points at which the functions are to
!           be evaluated.  If any of the XE are outside the interval
!           [X1,X2], a warning error is returned in NEXT.
!
!     FE -- (output) real*8 array of values of the cubic function
!           defined by  X1,X2, F1,F2, D1,D2  at the points  XE.
!
!     DE -- (output) real*8 array of values of the first derivative of
!           the same function at the points  XE.
!
!     NEXT -- (output) integer array indicating number of extrapolation
!           points:
!            NEXT(1) = number of evaluation points to left of interval.
!            NEXT(2) = number of evaluation points to right of interval.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  if NE<1 .
!              IERR = -2  if X1==X2 .
!                (Output arrays have not been changed in either case.)
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   811019  DATE WRITTEN
!   820803  Minor cosmetic changes for release 1.
!   870707  Corrected XERROR calls for d.p. names(s).
!   870813  Minor cosmetic changes.
!   890411  Added SAVE statements (Vers. 3.2).
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891006  Cosmetic changes to prologue.  (WRB)
!   891006  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DCHFDV

integer  ne, next(2), ierr
real(wp)  x1, x2, f1, f2, d1, d2, xe(*), fe(*), de(*)

integer  i
real(wp)  c2, c2t2, c3, c3t3, del1, del2, delta, h, x, xmi, xma

!  VALIDITY-CHECK ARGUMENTS.

if (ne < 1)  go to 5001
h = x2 - x1
if (h == zero)  go to 5002
!
!  INITIALIZE.
!
ierr = 0
next(1) = 0
next(2) = 0
xmi = min(zero, h)
xma = max(zero, h)
!
!  COMPUTE CUBIC COEFFICIENTS (EXPANDED ABOUT X1).
!
delta = (f2 - f1)/h
del1 = (d1 - delta)/h
del2 = (d2 - delta)/h
!                                           (DELTA IS NO LONGER NEEDED.)
c2 = -(del1+del1 + del2)
c2t2 = c2 + c2
c3 = (del1 + del2)/h
!                               (H, DEL1 AND DEL2 ARE NO LONGER NEEDED.)
c3t3 = c3+c3+c3
!
!  EVALUATION LOOP.
!
do 500  i = 1, ne
   x = xe(i) - x1
   fe(i) = f1 + x*(d1 + x*(c2 + x*c3))
   de(i) = d1 + x*(c2t2 + x*c3t3)
!          COUNT EXTRAPOLATION POINTS.
   if ( x<xmi )  next(1) = next(1) + 1
   if ( x>xma )  next(2) = next(2) + 1
!        (NOTE REDUNDANCY--IF EITHER CONDITION IS TRUE, OTHER IS FALSE.)
500 continue
!
!  NORMAL RETURN.
!
return
!
!  ERROR RETURNS.
!
5001 continue
!     NE<1 RETURN.
ierr = -1
call xermsg ('SLATEC', 'DCHFDV', &
   'NUMBER OF EVALUATION POINTS LESS THAN ONE', ierr, 1)
return
!
5002 continue
!     X1==X2 RETURN.
ierr = -2
call xermsg ('SLATEC', 'DCHFDV', 'INTERVAL ENDPOINTS EQUAL', &
   ierr, 1)

end subroutine dchfdv

subroutine dchfev (x1, x2, f1, f2, d1, d2, ne, xe, fe, next, ierr)

!***PURPOSE  Evaluate a cubic polynomial given in Hermite form at an
!            array of points.  While designed for use by DPCHFE, it may
!            be useful directly as an evaluator for a piecewise cubic
!            Hermite function in applications, such as graphing, where
!            the interval is known in advance.
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E3
!***TYPE      real(wp) (CHFEV-S, DCHFEV-D)
!***KEYWORDS  CUBIC HERMITE EVALUATION, CUBIC POLYNOMIAL EVALUATION,
!             PCHIP
!***AUTHOR  Fritsch, F. N., (LLNL)
!             Lawrence Livermore National Laboratory
!             P.O. Box 808  (L-316)
!             Livermore, CA  94550
!             FTS 532-4275, (510) 422-4275
!***DESCRIPTION
!
!          DCHFEV:  Cubic Hermite Function EValuator
!
!     Evaluates the cubic polynomial determined by function values
!     F1,F2 and derivatives D1,D2 on interval (X1,X2) at the points
!     XE(J), J=1(1)NE.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        INTEGER  NE, NEXT(2), IERR
!        real(wp)  X1, X2, F1, F2, D1, D2, XE(NE), FE(NE)
!
!        CALL  DCHFEV (X1,X2, F1,F2, D1,D2, NE, XE, FE, NEXT, IERR)
!
!   Parameters:
!
!     X1,X2 -- (input) endpoints of interval of definition of cubic.
!           (Error return if  X1==X2 .)
!
!     F1,F2 -- (input) values of function at X1 and X2, respectively.
!
!     D1,D2 -- (input) values of derivative at X1 and X2, respectively.
!
!     NE -- (input) number of evaluation points.  (Error return if
!           NE<1 .)
!
!     XE -- (input) real*8 array of points at which the function is to
!           be evaluated.  If any of the XE are outside the interval
!           [X1,X2], a warning error is returned in NEXT.
!
!     FE -- (output) real*8 array of values of the cubic function
!           defined by  X1,X2, F1,F2, D1,D2  at the points  XE.
!
!     NEXT -- (output) integer array indicating number of extrapolation
!           points:
!            NEXT(1) = number of evaluation points to left of interval.
!            NEXT(2) = number of evaluation points to right of interval.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  if NE<1 .
!              IERR = -2  if X1==X2 .
!                (The FE-array has not been changed in either case.)
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   811019  DATE WRITTEN
!   820803  Minor cosmetic changes for release 1.
!   870813  Corrected XERROR calls for d.p. names(s).
!   890206  Corrected XERROR calls.
!   890411  Added SAVE statements (Vers. 3.2).
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890703  Corrected category record.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891006  Cosmetic changes to prologue.  (WRB)
!   891006  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DCHFEV

integer  ne, next(2), ierr
real(wp)  x1, x2, f1, f2, d1, d2, xe(*), fe(*)

integer  i
real(wp)  c2, c3, del1, del2, delta, h, x, xmi, xma

!  VALIDITY-CHECK ARGUMENTS.

if (ne < 1)  go to 5001
h = x2 - x1
if (h == zero)  go to 5002
!
!  INITIALIZE.
!
ierr = 0
next(1) = 0
next(2) = 0
xmi = min(zero, h)
xma = max(zero, h)
!
!  COMPUTE CUBIC COEFFICIENTS (EXPANDED ABOUT X1).
!
delta = (f2 - f1)/h
del1 = (d1 - delta)/h
del2 = (d2 - delta)/h
!                                           (DELTA IS NO LONGER NEEDED.)
c2 = -(del1+del1 + del2)
c3 = (del1 + del2)/h
!                               (H, DEL1 AND DEL2 ARE NO LONGER NEEDED.)
!
!  EVALUATION LOOP.
!
do 500  i = 1, ne
   x = xe(i) - x1
   fe(i) = f1 + x*(d1 + x*(c2 + x*c3))
!          COUNT EXTRAPOLATION POINTS.
   if ( x<xmi )  next(1) = next(1) + 1
   if ( x>xma )  next(2) = next(2) + 1
!        (NOTE REDUNDANCY--IF EITHER CONDITION IS TRUE, OTHER IS FALSE.)
500 continue
!
!  NORMAL RETURN.
!
return
!
!  ERROR RETURNS.
!
5001 continue
!     NE<1 RETURN.
ierr = -1
call xermsg ('SLATEC', 'DCHFEV', &
   'NUMBER OF EVALUATION POINTS LESS THAN ONE', ierr, 1)
return
!
5002 continue
!     X1==X2 RETURN.
ierr = -2
call xermsg ('SLATEC', 'DCHFEV', 'INTERVAL ENDPOINTS EQUAL', &
   ierr, 1)

end subroutine dchfev

real(wp) function dchfie (x1, x2, f1, f2, d1, d2, a, b)


!***PURPOSE  Evaluates integral of a single cubic for DPCHIA
!***LIBRARY   SLATEC (PCHIP)
!***TYPE      real(wp) (CHFIE-S, DCHFIE-D)
!***AUTHOR  Fritsch, F. N., (LLNL)
!***DESCRIPTION
!
!          DCHFIE:  Cubic Hermite Function Integral Evaluator.
!
!     Called by  DPCHIA  to evaluate the integral of a single cubic (in
!     Hermite form) over an arbitrary interval (A,B).
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        real(wp)  X1, X2, F1, F2, D1, D2, A, B
!        real(wp)  VALUE, DCHFIE
!
!        VALUE = DCHFIE (X1, X2, F1, F2, D1, D2, A, B)
!
!   Parameters:
!
!     VALUE -- (output) value of the requested integral.
!
!     X1,X2 -- (input) endpoints if interval of definition of cubic.
!
!     F1,F2 -- (input) function values at the ends of the interval.
!
!     D1,D2 -- (input) derivative values at the ends of the interval.
!
!     A,B -- (input) endpoints of interval of integration.
!
!***SEE ALSO  DPCHIA
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   820730  DATE WRITTEN
!   820805  Converted to SLATEC library version.
!   870707  Corrected subroutine name from DCHIV to DCHFIV.
!   870813  Minor cosmetic changes.
!   890411  1. Added SAVE statements (Vers. 3.2).
!           2. Added SIX to real(wp) declaration.
!   890411  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900328  Added TYPE section.  (WRB)
!   910408  Updated AUTHOR section in prologue.  (WRB)
!   930503  Corrected to set VALUE=0 when IERR/=0.  (FNF)
!   930504  Eliminated IERR and changed name DCHFIV to DCHFIE.  (FNF)
!***END PROLOGUE  DCHFIE
!
!  Programming notes:
!  1. There is no error return from this routine because zero is
!     indeed the mathematically correct answer when X1==X2 .
!**End


real(wp)  x1, x2, f1, f2, d1, d2, a, b

real(wp)  dterm, fterm, h, phia1, phia2, &
      phib1, phib2, psia1, psia2, psib1, psib2, ta1, ta2, &
      tb1, tb2, ua1, ua2, ub1, ub2

!  VALIDITY CHECK INPUT.

if (x1 == x2)  then
   dchfie = 0
else
   h = x2 - x1
   ta1 = (a - x1) / h
   ta2 = (x2 - a) / h
   tb1 = (b - x1) / h
   tb2 = (x2 - b) / h
!
   ua1 = ta1**3
   phia1 = ua1 * (two - ta1)
   psia1 = ua1 * (three*ta1 - four)
   ua2 = ta2**3
   phia2 =  ua2 * (two - ta2)
   psia2 = -ua2 * (three*ta2 - four)
!
   ub1 = tb1**3
   phib1 = ub1 * (two - tb1)
   psib1 = ub1 * (three*tb1 - four)
   ub2 = tb2**3
   phib2 =  ub2 * (two - tb2)
   psib2 = -ub2 * (three*tb2 - four)
!
   fterm =   f1*(phia2 - phib2) + f2*(phib1 - phia1)
   dterm = ( d1*(psia2 - psib2) + d2*(psib1 - psia1) )*(h/six)
!
   dchfie = (half*h) * (fterm + dterm)
endif
!

end function dchfie

subroutine dpchbs (n, x, f, d, incfd, knotyp, nknots, t, bcoef, &
   ndim, kord, ierr)

!***PURPOSE  Piecewise Cubic Hermite to B-Spline converter.
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E3
!***TYPE      real(wp) (PCHBS-S, DPCHBS-D)
!***KEYWORDS  B-SPLINES, CONVERSION, CUBIC HERMITE INTERPOLATION,
!             PIECEWISE CUBIC INTERPOLATION
!***AUTHOR  Fritsch, F. N., (LLNL)
!             Computing and Mathematics Research Division
!             Lawrence Livermore National Laboratory
!             P.O. Box 808  (L-316)
!             Livermore, CA  94550
!             FTS 532-4275, (510) 422-4275
!***DESCRIPTION
!
! *Usage:
!
!        INTEGER  N, INCFD, KNOTYP, NKNOTS, NDIM, KORD, IERR
!        PARAMETER  (INCFD = ...)
!        real(wp)  X(nmax), F(INCFD,nmax), D(INCFD,nmax),
!       *      T(2*nmax+4), BCOEF(2*nmax)
!
!        CALL DPCHBS (N, X, F, D, INCFD, KNOTYP, NKNOTS, T, BCOEF,
!       *             NDIM, KORD, IERR)
!
! *Arguments:
!
!     N:IN  is the number of data points, N>=2 .  (not checked)
!
!     X:IN  is the real array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) < X(I),  I = 2(1)N.   (not checked)
!           nmax, the dimension of X, must be >=N.
!
!     F:IN  is the real array of dependent variable values.
!           F(1+(I-1)*INCFD) is the value corresponding to X(I).
!           nmax, the second dimension of F, must be >=N.
!
!     D:IN  is the real array of derivative values at the data points.
!           D(1+(I-1)*INCFD) is the value corresponding to X(I).
!           nmax, the second dimension of D, must be >=N.
!
!     INCFD:IN  is the increment between successive values in F and D.
!           This argument is provided primarily for 2-D applications.
!           It may have the value 1 for one-dimensional applications,
!           in which case F and D may be singly-subscripted arrays.
!
!     KNOTYP:IN  is a flag to control the knot sequence.
!           The knot sequence T is normally computed from X by putting
!           a double knot at each X and setting the end knot pairs
!           according to the value of KNOTYP:
!              KNOTYP = 0:  Quadruple knots at X(1) and X(N).  (default)
!              KNOTYP = 1:  Replicate lengths of extreme subintervals:
!                           T( 1 ) = T( 2 ) = X(1) - (X(2)-X(1))  ;
!                           T(M+4) = T(M+3) = X(N) + (X(N)-X(N-1)).
!              KNOTYP = 2:  Periodic placement of boundary knots:
!                           T( 1 ) = T( 2 ) = X(1) - (X(N)-X(N-1));
!                           T(M+4) = T(M+3) = X(N) + (X(2)-X(1))  .
!              Here M=NDIM=2*N.
!           If the input value of KNOTYP is negative, however, it is
!           assumed that NKNOTS and T were set in a previous call.
!           This option is provided for improved efficiency when used
!           in a parametric setting.
!
!     NKNOTS:INOUT  is the number of knots.
!           If KNOTYP>=0, then NKNOTS will be set to NDIM+4.
!           If KNOTYP<0, then NKNOTS is an input variable, and an
!              error return will be taken if it is not equal to NDIM+4.
!
!     T:INOUT  is the array of 2*N+4 knots for the B-representation.
!           If KNOTYP>=0, T will be returned by DPCHBS with the
!              interior double knots equal to the X-values and the
!              boundary knots set as indicated above.
!           If KNOTYP<0, it is assumed that T was set by a
!              previous call to DPCHBS.  (This routine does **not**
!              verify that T forms a legitimate knot sequence.)
!
!     BCOEF:OUT  is the array of 2*N B-spline coefficients.
!
!     NDIM:OUT  is the dimension of the B-spline space.  (Set to 2*N.)
!
!     KORD:OUT  is the order of the B-spline.  (Set to 4.)
!
!     IERR:OUT  is an error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -4  if KNOTYP>2 .
!              IERR = -5  if KNOTYP<0 and NKNOTS/=(2*N+4).
!
! *Description:
!     DPCHBS computes the B-spline representation of the PCH function
!     determined by N,X,F,D.  To be compatible with the rest of PCHIP,
!     DPCHBS includes INCFD, the increment between successive values of
!     the F- and D-arrays.
!
!     The output is the B-representation for the function:  NKNOTS, T,
!     BCOEF, NDIM, KORD.
!
! *Caution:
!     Since it is assumed that the input PCH function has been
!     computed by one of the other routines in the package PCHIP,
!     input arguments N, X, INCFD are **not** checked for validity.
!
! *Restrictions/assumptions:
!     1. N>=2 .  (not checked)
!     2. X(i)<X(i+1), i=1,...,N .  (not checked)
!     3. INCFD>0 .  (not checked)
!     4. KNOTYP<=2 .  (error return if not)
!    *5. NKNOTS = NDIM+4 = 2*N+4 .  (error return if not)
!    *6. T(2*k+1) = T(2*k) = X(k), k=1,...,N .  (not checked)
!
!       * Indicates this applies only if KNOTYP<0 .
!
! *Portability:
!     Argument INCFD is used only to cause the compiler to generate
!     efficient code for the subscript expressions (1+(I-1)*INCFD) .
!     The normal usage, in which DPCHBS is called with one-dimensional
!     arrays F and D, is probably non-Fortran 77, in the strict sense,
!     but it works on all systems on which DPCHBS has been tested.
!
! *See Also:
!     PCHIC, PCHIM, or PCHSP can be used to determine an interpolating
!        PCH function from a set of data.
!     The B-spline routine DBVALU can be used to evaluate the
!        B-representation that is output by DPCHBS.
!        (See BSPDOC for more information.)
!
!***REFERENCES  F. N. Fritsch, "Representations for parametric cubic
!                 splines," Computer Aided Geometric Design 6 (1989),
!                 pp.79-82.
!***ROUTINES CALLED  DPCHKT, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   870701  DATE WRITTEN
!   900405  Converted Fortran to upper case.
!   900405  Removed requirement that X be dimensioned N+1.
!   900406  Modified to make PCHKT a subsidiary routine to simplify
!           usage.  In the process, added argument INCFD to be com-
!           patible with the rest of PCHIP.
!   900410  Converted prologue to SLATEC 4.0 format.
!   900410  Added calls to XERMSG and changed constant 3. to 3 to
!           reduce single/double differences.
!   900411  Added reference.
!   900430  Produced real(wp) version.
!   900501  Corrected declarations.
!   930317  Minor cosmetic changes.  (FNF)
!   930514  Corrected problems with dimensioning of arguments and
!           clarified DESCRIPTION.  (FNF)
!   930604  Removed  NKNOTS from DPCHKT call list.  (FNF)
!***END PROLOGUE  DPCHBS
!
!*Internal Notes:
!
!**End


integer  n, incfd, knotyp, nknots, ndim, kord, ierr
real(wp)  x(*), f(incfd,*), d(incfd,*), t(*), bcoef(*)

integer  k, kk
real(wp)  dov3, hnew, hold
character*8  libnam, subnam

!  Initialize.
ndim = 2*n
kord = 4
ierr = 0
libnam = 'SLATEC'
subnam = 'DPCHBS'
!
!  Check argument validity.  Set up knot sequence if OK.
!
if ( knotyp>2 )  then
   ierr = -1
   call xermsg (libnam, subnam, 'KNOTYP GREATER THAN 2', ierr, 1)
   return
endif
if ( knotyp<0 )  then
   if ( nknots/=ndim+4 )  then
      ierr = -2
      call xermsg (libnam, subnam, &
                    'KNOTYP<0 AND NKNOTS/=(2*N+4)', ierr, 1)
      return
   endif
else
!          Set up knot sequence.
   nknots = ndim + 4
   call dpchkt (n, x, knotyp, t)
endif
!
!  Compute B-spline coefficients.
!
hnew = t(3) - t(1)
do 40  k = 1, n
   kk = 2*k
   hold = hnew
!          The following requires mixed mode arithmetic.
   dov3 = d(1,k)/3
   bcoef(kk-1) = f(1,k) - hold*dov3
!          The following assumes T(2*K+1) = X(K).
   hnew = t(kk+3) - t(kk+1)
   bcoef(kk) = f(1,k) + hnew*dov3
40 continue
!
!  Terminate.
!

end subroutine dpchbs

subroutine dpchce (ic, vc, n, x, h, slope, d, incfd, ierr)


!***PURPOSE  Set boundary conditions for DPCHIC
!***LIBRARY   SLATEC (PCHIP)
!***TYPE      real(wp) (PCHCE-S, DPCHCE-D)
!***AUTHOR  Fritsch, F. N., (LLNL)
!***DESCRIPTION
!
!          DPCHCE:  DPCHIC End Derivative Setter.
!
!    Called by DPCHIC to set end derivatives as requested by the user.
!    It must be called after interior derivative values have been set.
!                      -----
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the D-array.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  IC(2), N, IERR
!        real(wp)  VC(2), X(N), H(N), SLOPE(N), D(INCFD,N)
!
!        CALL  DPCHCE (IC, VC, N, X, H, SLOPE, D, INCFD, IERR)
!
!   Parameters:
!
!     IC -- (input) integer array of length 2 specifying desired
!           boundary conditions:
!           IC(1) = IBEG, desired condition at beginning of data.
!           IC(2) = IEND, desired condition at end of data.
!           ( see prologue to DPCHIC for details. )
!
!     VC -- (input) real*8 array of length 2 specifying desired boundary
!           values.  VC(1) need be set only if IC(1) = 2 or 3 .
!                    VC(2) need be set only if IC(2) = 2 or 3 .
!
!     N -- (input) number of data points.  (assumes N>=2)
!
!     X -- (input) real*8 array of independent variable values.  (the
!           elements of X are assumed to be strictly increasing.)
!
!     H -- (input) real*8 array of interval lengths.
!     SLOPE -- (input) real*8 array of data slopes.
!           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
!                  H(I) =  X(I+1)-X(I),
!              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
!
!     D -- (input) real*8 array of derivative values at the data points.
!           The value corresponding to X(I) must be stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!          (output) the value of D at X(1) and/or X(N) is changed, if
!           necessary, to produce the requested boundary conditions.
!           no other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in D.
!           This argument is provided primarily for 2-D applications.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning errors:
!              IERR = 1  if IBEG<0 and D(1) had to be adjusted for
!                        monotonicity.
!              IERR = 2  if IEND<0 and D(1+(N-1)*INCFD) had to be
!                        adjusted for monotonicity.
!              IERR = 3  if both of the above are true.
!
!    -------
!    WARNING:  This routine does no validity-checking of arguments.
!    -------
!***SEE ALSO  DPCHIC
!***ROUTINES CALLED  DPCHDF, DPCHST, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   820218  DATE WRITTEN
!   820805  Converted to SLATEC library version.
!   870707  Corrected XERROR calls for d.p. name(s).
!   890206  Corrected XERROR calls.
!   890411  Added SAVE statements (Vers. 3.2).
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900328  Added TYPE section.  (WRB)
!   910408  Updated AUTHOR section in prologue.  (WRB)
!   930503  Improved purpose.  (FNF)
!***END PROLOGUE  DPCHCE
!
!  Programming notes:
!     1. The function DPCHST(ARG1,ARG2)  is assumed to return zero if
!        either argument is zero, +1 if they are of the same sign, and
!        -1 if they are of opposite sign.
!     2. One could reduce the number of arguments and amount of local
!        storage, at the expense of reduced code clarity, by passing in
!        the array WK (rather than splitting it into H and SLOPE) and
!        increasing its length enough to incorporate STEMP and XTEMP.
!     3. The two monotonicity checks only use the sufficient conditions.
!        Thus, it is possible (but unlikely) for a boundary condition to
!        be changed, even though the original interpolant was monotonic.
!        (At least the result is a continuous function of the data.)
!**End


integer  ic(2), n, incfd, ierr
real(wp)  vc(2), x(*), h(*), slope(*), d(incfd,*)

integer  ibeg, iend, ierf, index, j, k
real(wp)  stemp(3), xtemp(4)

ibeg = ic(1)
iend = ic(2)
ierr = 0
!
!  SET TO DEFAULT BOUNDARY CONDITIONS IF N IS TOO SMALL.
!
if ( abs(ibeg)>n )  ibeg = 0
if ( abs(iend)>n )  iend = 0
!
!  TREAT BEGINNING BOUNDARY CONDITION.
!
if (ibeg == 0)  go to 2000
k = abs(ibeg)
if (k == 1)  then
!        BOUNDARY VALUE PROVIDED.
   d(1,1) = vc(1)
else if (k == 2)  then
!        BOUNDARY SECOND DERIVATIVE PROVIDED.
   d(1,1) = half*( (three*slope(1) - d(1,2)) - half*vc(1)*h(1) )
else if (k < 5)  then
!        USE K-POINT DERIVATIVE FORMULA.
!        PICK UP FIRST K POINTS, IN REVERSE ORDER.
   do 10  j = 1, k
      index = k-j+1
!           INDEX RUNS FROM K DOWN TO 1.
      xtemp(j) = x(index)
      if (j < k)  stemp(j) = slope(index-1)
10    continue
!                 -----------------------------
   d(1,1) = dpchdf (k, xtemp, stemp, ierf)
!                 -----------------------------
   if (ierf /= 0)  go to 5001
else
!        USE 'NOT A KNOT' CONDITION.
   d(1,1) = ( three*(h(1)*slope(2) + h(2)*slope(1)) &
             - two*(h(1)+h(2))*d(1,2) - h(1)*d(1,3) ) / h(2)
endif
!
if (ibeg > 0)  go to 2000
!
!  CHECK D(1,1) FOR COMPATIBILITY WITH MONOTONICITY.
!
if (slope(1) == zero)  then
   if (d(1,1) /= zero)  then
      d(1,1) = zero
      ierr = ierr + 1
   endif
else if ( dpchst(d(1,1),slope(1)) < zero)  then
   d(1,1) = zero
   ierr = ierr + 1
else if ( abs(d(1,1)) > three*abs(slope(1)) )  then
   d(1,1) = three*slope(1)
   ierr = ierr + 1
endif
!
!  TREAT END BOUNDARY CONDITION.
!
2000 continue
if (iend == 0)  go to 5000
k = abs(iend)
if (k == 1)  then
!        BOUNDARY VALUE PROVIDED.
   d(1,n) = vc(2)
else if (k == 2)  then
!        BOUNDARY SECOND DERIVATIVE PROVIDED.
   d(1,n) = half*( (three*slope(n-1) - d(1,n-1)) + &
                                           half*vc(2)*h(n-1) )
else if (k < 5)  then
!        USE K-POINT DERIVATIVE FORMULA.
!        PICK UP LAST K POINTS.
   do 2010  j = 1, k
      index = n-k+j
!           INDEX RUNS FROM N+1-K UP TO N.
      xtemp(j) = x(index)
      if (j < k)  stemp(j) = slope(index)
2010    continue
!                 -----------------------------
   d(1,n) = dpchdf (k, xtemp, stemp, ierf)
!                 -----------------------------
   if (ierf /= 0)  go to 5001
else
!        USE 'NOT A KNOT' CONDITION.
   d(1,n) = ( three*(h(n-1)*slope(n-2) + h(n-2)*slope(n-1)) &
             - two*(h(n-1)+h(n-2))*d(1,n-1) - h(n-1)*d(1,n-2) ) &
                                                         / h(n-2)
endif
!
if (iend > 0)  go to 5000
!
!  CHECK D(1,N) FOR COMPATIBILITY WITH MONOTONICITY.
!
if (slope(n-1) == zero)  then
   if (d(1,n) /= zero)  then
      d(1,n) = zero
      ierr = ierr + 2
   endif
else if ( dpchst(d(1,n),slope(n-1)) < zero)  then
   d(1,n) = zero
   ierr = ierr + 2
else if ( abs(d(1,n)) > three*abs(slope(n-1)) )  then
   d(1,n) = three*slope(n-1)
   ierr = ierr + 2
endif
!
!  NORMAL RETURN.
!
5000 continue
return
!
!  ERROR RETURN.
!
5001 continue
!     ERROR RETURN FROM DPCHDF.
!   *** THIS CASE SHOULD NEVER OCCUR ***
ierr = -1
call xermsg ('SLATEC', 'DPCHCE', 'ERROR RETURN FROM DPCHDF', &
   ierr, 1)

end subroutine dpchce

subroutine dpchci (n, h, slope, d, incfd)


!***PURPOSE  Set interior derivatives for DPCHIC
!***LIBRARY   SLATEC (PCHIP)
!***TYPE      real(wp) (PCHCI-S, DPCHCI-D)
!***AUTHOR  Fritsch, F. N., (LLNL)
!***DESCRIPTION
!
!          DPCHCI:  DPCHIC Initial Derivative Setter.
!
!    Called by DPCHIC to set derivatives needed to determine a monotone
!    piecewise cubic Hermite interpolant to the data.
!
!    Default boundary conditions are provided which are compatible
!    with monotonicity.  If the data are only piecewise monotonic, the
!    interpolant will have an extremum at each point where monotonicity
!    switches direction.
!
!    To facilitate two-dimensional applications, includes an increment
!    between successive values of the D-array.
!
!    The resulting piecewise cubic Hermite function should be identical
!    (within roundoff error) to that produced by DPCHIM.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  N
!        real(wp)  H(N), SLOPE(N), D(INCFD,N)
!
!        CALL  DPCHCI (N, H, SLOPE, D, INCFD)
!
!   Parameters:
!
!     N -- (input) number of data points.
!           If N=2, simply does linear interpolation.
!
!     H -- (input) real*8 array of interval lengths.
!     SLOPE -- (input) real*8 array of data slopes.
!           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
!                  H(I) =  X(I+1)-X(I),
!              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
!
!     D -- (output) real*8 array of derivative values at data points.
!           If the data are monotonic, these values will determine a
!           a monotone cubic Hermite function.
!           The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in D.
!           This argument is provided primarily for 2-D applications.
!
!    -------
!    WARNING:  This routine does no validity-checking of arguments.
!    -------
!
!***SEE ALSO  DPCHIC
!***ROUTINES CALLED  DPCHST
!***REVISION HISTORY  (YYMMDD)
!   820218  DATE WRITTEN
!   820601  Modified end conditions to be continuous functions of
!           data when monotonicity switches in next interval.
!   820602  1. Modified formulas so end conditions are less prone
!             to over/underflow problems.
!           2. Minor modification to HSUM calculation.
!   820805  Converted to SLATEC library version.
!   870813  Minor cosmetic changes.
!   890411  Added SAVE statements (Vers. 3.2).
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910408  Updated AUTHOR section in prologue.  (WRB)
!   930503  Improved purpose.  (FNF)
!***END PROLOGUE  DPCHCI
!
!  Programming notes:
!     1. The function  DPCHST(ARG1,ARG2)  is assumed to return zero if
!        either argument is zero, +1 if they are of the same sign, and
!        -1 if they are of opposite sign.
!**End


integer  n, incfd
real(wp)  h(*), slope(*), d(incfd,*)

integer  i, nless1
real(wp)  del1, del2, dmax, dmin, drat1, drat2, hsum, &
      hsumt3, w1, w2

nless1 = n - 1
del1 = slope(1)
!
!  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.
!
if (nless1 > 1)  go to 10
d(1,1) = del1
d(1,n) = del1
go to 5000
!
!  NORMAL CASE  (N >= 3).
!
10 continue
del2 = slope(2)
!
!  SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!     SHAPE-PRESERVING.
!
hsum = h(1) + h(2)
w1 = (h(1) + hsum)/hsum
w2 = -h(1)/hsum
d(1,1) = w1*del1 + w2*del2
if ( dpchst(d(1,1),del1) <= zero)  then
   d(1,1) = zero
else if ( dpchst(del1,del2) < zero)  then
!        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
   dmax = three*del1
   if (abs(d(1,1)) > abs(dmax))  d(1,1) = dmax
endif
!
!  LOOP THROUGH INTERIOR POINTS.
!
do 50  i = 2, nless1
   if (i == 2)  go to 40
!
   hsum = h(i-1) + h(i)
   del1 = del2
   del2 = slope(i)
40    continue
!
!        SET D(I)=0 UNLESS DATA ARE STRICTLY MONOTONIC.
!
   d(1,i) = zero
   if ( dpchst(del1,del2) <= zero)  go to 50
!
!        USE BRODLIE MODIFICATION OF BUTLAND FORMULA.
!
   hsumt3 = hsum+hsum+hsum
   w1 = (hsum + h(i-1))/hsumt3
   w2 = (hsum + h(i)  )/hsumt3
   dmax = max( abs(del1), abs(del2) )
   dmin = min( abs(del1), abs(del2) )
   drat1 = del1/dmax
   drat2 = del2/dmax
   d(1,i) = dmin/(w1*drat1 + w2*drat2)
!
50 continue
!
!  SET D(N) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!     SHAPE-PRESERVING.
!
w1 = -h(n-1)/hsum
w2 = (h(n-1) + hsum)/hsum
d(1,n) = w1*del1 + w2*del2
if ( dpchst(d(1,n),del2) <= zero)  then
   d(1,n) = zero
else if ( dpchst(del1,del2) < zero)  then
!        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
   dmax = three*del2
   if (abs(d(1,n)) > abs(dmax))  d(1,n) = dmax
endif
!
!  NORMAL RETURN.
!
5000 continue

end subroutine dpchci

subroutine dpchcm (n, x, f, d, incfd, skip, ismon, ierr)

!***PURPOSE  Check a cubic Hermite function for monotonicity.
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E3
!***TYPE      real(wp) (PCHCM-S, DPCHCM-D)
!***KEYWORDS  CUBIC HERMITE INTERPOLATION, MONOTONE INTERPOLATION,
!             PCHIP, PIECEWISE CUBIC INTERPOLATION, UTILITY ROUTINE
!***AUTHOR  Fritsch, F. N., (LLNL)
!             Computing & Mathematics Research Division
!             Lawrence Livermore National Laboratory
!             P.O. Box 808  (L-316)
!             Livermore, CA  94550
!             FTS 532-4275, (510) 422-4275
!***DESCRIPTION
!
! *Usage:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  N, ISMON(N), IERR
!        real(wp)  X(N), F(INCFD,N), D(INCFD,N)
!        LOGICAL  SKIP
!
!        CALL  DPCHCM (N, X, F, D, INCFD, SKIP, ISMON, IERR)
!
! *Arguments:
!
!     N:IN  is the number of data points.  (Error return if N<2 .)
!
!     X:IN  is a real*8 array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) < X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F:IN  is a real*8 array of function values.  F(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     D:IN  is a real*8 array of derivative values.  D(1+(I-1)*INCFD) is
!           is the value corresponding to X(I).
!
!     INCFD:IN  is the increment between successive values in F and D.
!           (Error return if  INCFD<1 .)
!
!     SKIP:INOUT  is a logical variable which should be set to
!           .TRUE. if the user wishes to skip checks for validity of
!           preceding parameters, or to .FALSE. otherwise.
!           This will save time in case these checks have already
!           been performed.
!           SKIP will be set to .TRUE. on normal return.
!
!     ISMON:OUT  is an integer array indicating on which intervals the
!           PCH function defined by  N, X, F, D  is monotonic.
!           For data interval [X(I),X(I+1)],
!             ISMON(I) = -3  if function is probably decreasing;
!             ISMON(I) = -1  if function is strictly decreasing;
!             ISMON(I) =  0  if function is constant;
!             ISMON(I) =  1  if function is strictly increasing;
!             ISMON(I) =  2  if function is non-monotonic;
!             ISMON(I) =  3  if function is probably increasing.
!                If ABS(ISMON)=3, this means that the D-values are near
!                the boundary of the monotonicity region.  A small
!                increase produces non-monotonicity; decrease, strict
!                monotonicity.
!           The above applies to I=1(1)N-1.  ISMON(N) indicates whether
!              the entire function is monotonic on [X(1),X(N)].
!
!     IERR:OUT  is an error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  if N<2 .
!              IERR = -2  if INCFD<1 .
!              IERR = -3  if the X-array is not strictly increasing.
!          (The ISMON-array has not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
! *Description:
!
!          DPCHCM:  Piecewise Cubic Hermite -- Check Monotonicity.
!
!     Checks the piecewise cubic Hermite function defined by  N,X,F,D
!     for monotonicity.
!
!     To provide compatibility with DPCHIM and DPCHIC, includes an
!     increment between successive values of the F- and D-arrays.
!
! *Cautions:
!     This provides the same capability as old DPCHMC, except that a
!     new output value, -3, was added February 1989.  (Formerly, -3
!     and +3 were lumped together in the single value 3.)  Codes that
!     flag nonmonotonicity by "IF (ISMON==2)" need not be changed.
!     Codes that check via "IF (ISMON>=3)" should change the test to
!     "IF (IABS(ISMON)>=3)".  Codes that declare monotonicity via
!     "IF (ISMON<=1)" should change to "IF (IABS(ISMON)<=1)".
!
!***REFERENCES  F. N. Fritsch and R. E. Carlson, Monotone piecewise
!                 cubic interpolation, SIAM Journal on Numerical Ana-
!                 lysis 17, 2 (April 1980), pp. 238-246.
!***ROUTINES CALLED  DCHFCM, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   820518  DATE WRITTEN
!   820804  Converted to SLATEC library version.
!   831201  Reversed order of subscripts of F and D, so that the
!           routine will work properly when INCFD>1 .  (Bug!)
!   870707  Corrected XERROR calls for d.p. name(s).
!   890206  Corrected XERROR calls.
!   890209  Added possible ISMON value of -3 and modified code so
!           that 1,3,-1 produces ISMON(N)=2, rather than 3.
!   890306  Added caution about changed output.
!   890407  Changed name from DPCHMC to DPCHCM, as requested at the
!           March 1989 SLATEC CML meeting, and made a few other
!           minor modifications necessitated by this change.
!   890407  Converted to new SLATEC format.
!   890407  Modified DESCRIPTION to LDOC format.
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920429  Revised format and order of references.  (WRB,FNF)
!***END PROLOGUE  DPCHCM
!
! ----------------------------------------------------------------------
!
!  Programming notes:
!
!     An alternate organization would have separate loops for computing
!     ISMON(i), i=1,...,NSEG, and for the computation of ISMON(N).  The
!     first loop can be readily parallelized, since the NSEG calls to
!     CHFCM are independent.  The second loop can be cut short if
!     ISMON(N) is ever equal to 2, for it cannot be changed further.


integer n, incfd, ismon(n), ierr
real(wp)  x(n), f(incfd,n), d(incfd,n)
logical  skip

integer i, nseg
real(wp)  delta

!  VALIDITY-CHECK ARGUMENTS.

if (skip)  go to 5
!
if ( n<2 )  go to 5001
if ( incfd<1 )  go to 5002
do 1  i = 2, n
   if ( x(i)<=x(i-1) )  go to 5003
1 continue
skip = .true.
!
!  FUNCTION DEFINITION IS OK -- GO ON.
!
5 continue
nseg = n - 1
do 90  i = 1, nseg
   delta = (f(1,i+1)-f(1,i))/(x(i+1)-x(i))
!                   -------------------------------
   ismon(i) = dchfcm (d(1,i), d(1,i+1), delta)
!                   -------------------------------
   if (i == 1)  then
      ismon(n) = ismon(1)
   else
!           Need to figure out cumulative monotonicity from following
!           "multiplication table":
!
!                    +        I S M O N (I)
!                     +  -3  -1   0   1   3   2
!                      +------------------------+
!               I   -3 I -3  -3  -3   2   2   2 I
!               S   -1 I -3  -1  -1   2   2   2 I
!               M    0 I -3  -1   0   1   3   2 I
!               O    1 I  2   2   1   1   3   2 I
!               N    3 I  2   2   3   3   3   2 I
!              (N)   2 I  2   2   2   2   2   2 I
!                      +------------------------+
!           Note that the 2 row and column are out of order so as not
!           to obscure the symmetry in the rest of the table.
!
!           No change needed if equal or constant on this interval or
!           already declared nonmonotonic.
      if ( (ismon(i)/=ismon(n)) .and. (ismon(i)/=0) &
                                  .and. (ismon(n)/=2) )  then
         if ( (ismon(i)==2) .or. (ismon(n)==0) )  then
            ismon(n) =  ismon(i)
         else if (ismon(i)*ismon(n) < 0)  then
!                 This interval has opposite sense from curve so far.
            ismon(n) = 2
         else
!                 At this point, both are nonzero with same sign, and
!                 we have already eliminated case both +-1.
            ismon(n) = sign (3, ismon(n))
         endif
      endif
   endif
90 continue
!
!  NORMAL RETURN.
!
ierr = 0
return
!
!  ERROR RETURNS.
!
5001 continue
!     N<2 RETURN.
ierr = -1
call xermsg ('SLATEC', 'DPCHCM', &
   'NUMBER OF DATA POINTS LESS THAN TWO', ierr, 1)
return
!
5002 continue
!     INCFD<1 RETURN.
ierr = -2
call xermsg ('SLATEC', 'DPCHCM', 'INCREMENT LESS THAN ONE', ierr, &
   1)
return
!
5003 continue
!     X-ARRAY NOT STRICTLY INCREASING.
ierr = -3
call xermsg ('SLATEC', 'DPCHCM', &
   'X-ARRAY NOT STRICTLY INCREASING', ierr, 1)

end subroutine dpchcm

subroutine dpchcs (switch, n, h, slope, d, incfd, ierr)


!***PURPOSE  Adjusts derivative values for DPCHIC
!***LIBRARY   SLATEC (PCHIP)
!***TYPE      real(wp) (PCHCS-S, DPCHCS-D)
!***AUTHOR  Fritsch, F. N., (LLNL)
!***DESCRIPTION
!
!         DPCHCS:  DPCHIC Monotonicity Switch Derivative Setter.
!
!     Called by  DPCHIC  to adjust the values of D in the vicinity of a
!     switch in direction of monotonicity, to produce a more "visually
!     pleasing" curve than that given by  DPCHIM .
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  N, IERR
!        real(wp)  SWITCH, H(N), SLOPE(N), D(INCFD,N)
!
!        CALL  DPCHCS (SWITCH, N, H, SLOPE, D, INCFD, IERR)
!
!   Parameters:
!
!     SWITCH -- (input) indicates the amount of control desired over
!           local excursions from data.
!
!     N -- (input) number of data points.  (assumes N>2 .)
!
!     H -- (input) real*8 array of interval lengths.
!     SLOPE -- (input) real*8 array of data slopes.
!           If the data are (X(I),Y(I)), I=1(1)N, then these inputs are:
!                  H(I) =  X(I+1)-X(I),
!              SLOPE(I) = (Y(I+1)-Y(I))/H(I),  I=1(1)N-1.
!
!     D -- (input) real*8 array of derivative values at the data points,
!           as determined by DPCHCI.
!          (output) derivatives in the vicinity of switches in direction
!           of monotonicity may be adjusted to produce a more "visually
!           pleasing" curve.
!           The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in D.
!           This argument is provided primarily for 2-D applications.
!
!     IERR -- (output) error flag.  should be zero.
!           If negative, trouble in DPCHSW.  (should never happen.)
!
!    -------
!    WARNING:  This routine does no validity-checking of arguments.
!    -------
!
!***SEE ALSO  DPCHIC
!***ROUTINES CALLED  DPCHST, DPCHSW
!***REVISION HISTORY  (YYMMDD)
!   820218  DATE WRITTEN
!   820617  Redesigned to (1) fix  problem with lack of continuity
!           approaching a flat-topped peak (2) be cleaner and
!           easier to verify.
!           Eliminated subroutines PCHSA and PCHSX in the process.
!   820622  1. Limited fact to not exceed one, so computed D is a
!             convex combination of DPCHCI value and DPCHSD value.
!           2. Changed fudge from 1 to 4 (based on experiments).
!   820623  Moved PCHSD to an inline function (eliminating MSWTYP).
!   820805  Converted to SLATEC library version.
!   870707  Corrected conversion to real(wp).
!   870813  Minor cosmetic changes.
!   890411  Added SAVE statements (Vers. 3.2).
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891006  Modified spacing in computation of DFLOC.  (WRB)
!   891006  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910408  Updated AUTHOR section in prologue.  (WRB)
!   930503  Improved purpose.  (FNF)
!***END PROLOGUE  DPCHCS
!
!  Programming notes:
!     1. The function  DPCHST(ARG1,ARG2)  is assumed to return zero if
!        either argument is zero, +1 if they are of the same sign, and
!        -1 if they are of opposite sign.
!**End


integer  n, incfd, ierr
real(wp)  switch, h(*), slope(*), d(incfd,*)

integer  i, indx, k, nless1
real(wp)  del(3), dext, dfloc, dfmx, fact, slmax, wtave(2)

real(wp),parameter :: fudge = four

ierr = 0
nless1 = n - 1
!
!  LOOP OVER SEGMENTS.
!
do 900  i = 2, nless1
   if ( dpchst(slope(i-1),slope(i)) )  100, 300, 900
!             --------------------------
!
100    continue
!
!....... SLOPE SWITCHES MONOTONICITY AT I-TH POINT .....................
!
!           DO NOT CHANGE D IF 'UP-DOWN-UP'.
      if (i > 2)  then
         if ( dpchst(slope(i-2),slope(i)) > zero)  go to 900
!                   --------------------------
      endif
      if (i < nless1)  then
         if ( dpchst(slope(i+1),slope(i-1)) > zero)  go to 900
!                   ----------------------------
      endif
!
!   ....... COMPUTE PROVISIONAL VALUE FOR D(1,I).
!
      dext = dpchsd (slope(i-1), slope(i), h(i-1), h(i))
!
!   ....... DETERMINE WHICH INTERVAL CONTAINS THE EXTREMUM.
!
      if ( dpchst(dext, slope(i-1)) )  200, 900, 250
!                -----------------------
!
200       continue
!              DEXT AND SLOPE(I-1) HAVE OPPOSITE SIGNS --
!                        EXTREMUM IS IN (X(I-1),X(I)).
         k = i-1
!              SET UP TO COMPUTE NEW VALUES FOR D(1,I-1) AND D(1,I).
         wtave(2) = dext
         if (k > 1) &
            wtave(1) = dpchsd (slope(k-1), slope(k), h(k-1), h(k))
         go to 400
!
250       continue
!              DEXT AND SLOPE(I) HAVE OPPOSITE SIGNS --
!                        EXTREMUM IS IN (X(I),X(I+1)).
         k = i
!              SET UP TO COMPUTE NEW VALUES FOR D(1,I) AND D(1,I+1).
         wtave(1) = dext
         if (k < nless1) &
            wtave(2) = dpchsd (slope(k), slope(k+1), h(k), h(k+1))
         go to 400
!
300    continue
!
!....... AT LEAST ONE OF SLOPE(I-1) AND SLOPE(I) IS ZERO --
!                     CHECK FOR FLAT-TOPPED PEAK .......................
!
      if (i == nless1)  go to 900
      if ( dpchst(slope(i-1), slope(i+1)) >= zero)  go to 900
!                -----------------------------
!
!           WE HAVE FLAT-TOPPED PEAK ON (X(I),X(I+1)).
      k = i
!           SET UP TO COMPUTE NEW VALUES FOR D(1,I) AND D(1,I+1).
      wtave(1) = dpchsd (slope(k-1), slope(k), h(k-1), h(k))
      wtave(2) = dpchsd (slope(k), slope(k+1), h(k), h(k+1))
!
400    continue
!
!....... AT THIS POINT WE HAVE DETERMINED THAT THERE WILL BE AN EXTREMUM
!        ON (X(K),X(K+1)), WHERE K=I OR I-1, AND HAVE SET ARRAY WTAVE--
!           WTAVE(1) IS A WEIGHTED AVERAGE OF SLOPE(K-1) AND SLOPE(K),
!                    IF K>1
!           WTAVE(2) IS A WEIGHTED AVERAGE OF SLOPE(K) AND SLOPE(K+1),
!                    IF K<N-1
!
   slmax = abs(slope(k))
   if (k > 1)    slmax = max( slmax, abs(slope(k-1)) )
   if (k<nless1) slmax = max( slmax, abs(slope(k+1)) )
!
   if (k > 1)  del(1) = slope(k-1) / slmax
   del(2) = slope(k) / slmax
   if (k<nless1)  del(3) = slope(k+1) / slmax
!
   if ((k>1) .and. (k<nless1))  then
!           NORMAL CASE -- EXTREMUM IS NOT IN A BOUNDARY INTERVAL.
      fact = fudge* abs(del(3)*(del(1)-del(2))*(wtave(2)/slmax))
      d(1,k) = d(1,k) + min(fact,one)*(wtave(1) - d(1,k))
      fact = fudge* abs(del(1)*(del(3)-del(2))*(wtave(1)/slmax))
      d(1,k+1) = d(1,k+1) + min(fact,one)*(wtave(2) - d(1,k+1))
   else
!           SPECIAL CASE K=1 (WHICH CAN OCCUR ONLY IF I=2) OR
!                        K=NLESS1 (WHICH CAN OCCUR ONLY IF I=NLESS1).
      fact = fudge* abs(del(2))
      d(1,i) = min(fact,one) * wtave(i-k+1)
!              NOTE THAT I-K+1 = 1 IF K=I  (=NLESS1),
!                        I-K+1 = 2 IF K=I-1(=1).
   endif
!
!
!....... ADJUST IF NECESSARY TO LIMIT EXCURSIONS FROM DATA.
!
   if (switch <= zero)  go to 900
!
   dfloc = h(k)*abs(slope(k))
   if (k > 1)    dfloc = max( dfloc, h(k-1)*abs(slope(k-1)) )
   if (k<nless1) dfloc = max( dfloc, h(k+1)*abs(slope(k+1)) )
   dfmx = switch*dfloc
   indx = i-k+1
!        INDX = 1 IF K=I, 2 IF K=I-1.
!        ---------------------------------------------------------------
   call dpchsw(dfmx, indx, d(1,k), d(1,k+1), h(k), slope(k), ierr)
!        ---------------------------------------------------------------
   if (ierr /= 0)  return
!
!....... END OF SEGMENT LOOP.
!
900 continue
!

end subroutine dpchcs

real(wp) function dpchdf (k, x, s, ierr)


!***PURPOSE  Computes divided differences for DPCHCE and DPCHSP
!***LIBRARY   SLATEC (PCHIP)
!***TYPE      real(wp) (PCHDF-S, DPCHDF-D)
!***AUTHOR  Fritsch, F. N., (LLNL)
!***DESCRIPTION
!
!          DPCHDF:   DPCHIP Finite Difference Formula
!
!     Uses a divided difference formulation to compute a K-point approx-
!     imation to the derivative at X(K) based on the data in X and S.
!
!     Called by  DPCHCE  and  DPCHSP  to compute 3- and 4-point boundary
!     derivative approximations.
!
! ----------------------------------------------------------------------
!
!     On input:
!        K      is the order of the desired derivative approximation.
!               K must be at least 3 (error return if not).
!        X      contains the K values of the independent variable.
!               X need not be ordered, but the values **MUST** be
!               distinct.  (Not checked here.)
!        S      contains the associated slope values:
!                  S(I) = (F(I+1)-F(I))/(X(I+1)-X(I)), I=1(1)K-1.
!               (Note that S need only be of length K-1.)
!
!     On return:
!        S      will be destroyed.
!        IERR   will be set to -1 if K<2 .
!        DPCHDF  will be set to the desired derivative approximation if
!               IERR=0 or to zero if IERR=-1.
!
! ----------------------------------------------------------------------
!
!***SEE ALSO  DPCHCE, DPCHSP
!***REFERENCES  Carl de Boor, A Practical Guide to Splines, Springer-
!                 Verlag, New York, 1978, pp. 10-16.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   820503  DATE WRITTEN
!   820805  Converted to SLATEC library version.
!   870707  Corrected XERROR calls for d.p. name(s).
!   870813  Minor cosmetic changes.
!   890206  Corrected XERROR calls.
!   890411  Added SAVE statements (Vers. 3.2).
!   890411  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900328  Added TYPE section.  (WRB)
!   910408  Updated AUTHOR and DATE WRITTEN sections in prologue.  (WRB)
!   920429  Revised format and order of references.  (WRB,FNF)
!   930503  Improved purpose.  (FNF)
!***END PROLOGUE  DPCHDF
!
!**End


integer  k, ierr
real(wp)  x(k), s(k)

integer  i, j
real(wp)  value

!  CHECK FOR LEGAL VALUE OF K.

if (k < 3)  go to 5001
!
!  COMPUTE COEFFICIENTS OF INTERPOLATING POLYNOMIAL.
!
do 10  j = 2, k-1
   do 9  i = 1, k-j
      s(i) = (s(i+1)-s(i))/(x(i+j)-x(i))
9    continue
10 continue
!
!  EVALUATE DERIVATIVE AT X(K).
!
value = s(1)
do 20  i = 2, k-1
   value = s(i) + value*(x(k)-x(i))
20 continue
!
!  NORMAL RETURN.
!
ierr = 0
dpchdf = value
return
!
!  ERROR RETURN.
!
5001 continue
!     K<3 RETURN.
ierr = -1
call xermsg ('SLATEC', 'DPCHDF', 'K LESS THAN THREE', ierr, 1)
dpchdf = zero

end function dpchdf

subroutine dpchfd (n, x, f, d, incfd, skip, ne, xe, fe, de, ierr)

!***PURPOSE  Evaluate a piecewise cubic Hermite function and its first
!            derivative at an array of points.  May be used by itself
!            for Hermite interpolation, or as an evaluator for DPCHIM
!            or DPCHIC. If only function values are required, use
!            DPCHFE instead.
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E3, H1
!***TYPE      real(wp) (PCHFD-S, DPCHFD-D)
!***KEYWORDS  CUBIC HERMITE DIFFERENTIATION, CUBIC HERMITE EVALUATION,
!             HERMITE INTERPOLATION, PCHIP, PIECEWISE CUBIC EVALUATION
!***AUTHOR  Fritsch, F. N., (LLNL)
!             Lawrence Livermore National Laboratory
!             P.O. Box 808  (L-316)
!             Livermore, CA  94550
!             FTS 532-4275, (510) 422-4275
!***DESCRIPTION
!
!          DPCHFD:  Piecewise Cubic Hermite Function and Derivative
!                  evaluator
!
!     Evaluates the cubic Hermite function defined by  N, X, F, D,  to-
!     gether with its first derivative, at the points  XE(J), J=1(1)NE.
!
!     If only function values are required, use DPCHFE, instead.
!
!     To provide compatibility with DPCHIM and DPCHIC, includes an
!     increment between successive values of the F- and D-arrays.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  N, NE, IERR
!        real(wp)  X(N), F(INCFD,N), D(INCFD,N), XE(NE), FE(NE),
!                          DE(NE)
!        LOGICAL  SKIP
!
!        CALL  DPCHFD (N, X, F, D, INCFD, SKIP, NE, XE, FE, DE, IERR)
!
!   Parameters:
!
!     N -- (input) number of data points.  (Error return if N<2 .)
!
!     X -- (input) real*8 array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) < X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real*8 array of function values.  F(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     D -- (input) real*8 array of derivative values.  D(1+(I-1)*INCFD)
!           is the value corresponding to X(I).
!
!     INCFD -- (input) increment between successive values in F and D.
!           (Error return if  INCFD<1 .)
!
!     SKIP -- (input/output) logical variable which should be set to
!           .TRUE. if the user wishes to skip checks for validity of
!           preceding parameters, or to .FALSE. otherwise.
!           This will save time in case these checks have already
!           been performed (say, in DPCHIM or DPCHIC).
!           SKIP will be set to .TRUE. on normal return.
!
!     NE -- (input) number of evaluation points.  (Error return if
!           NE<1 .)
!
!     XE -- (input) real*8 array of points at which the functions are to
!           be evaluated.
!
!
!          NOTES:
!           1. The evaluation will be most efficient if the elements
!              of XE are increasing relative to X;
!              that is,   XE(J) >= X(I)
!              implies    XE(K) >= X(I),  all K>=J .
!           2. If any of the XE are outside the interval [X(1),X(N)],
!              values are extrapolated from the nearest extreme cubic,
!              and a warning error is returned.
!
!     FE -- (output) real*8 array of values of the cubic Hermite
!           function defined by  N, X, F, D  at the points  XE.
!
!     DE -- (output) real*8 array of values of the first derivative of
!           the same function at the points  XE.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR>0  means that extrapolation was performed at
!                 IERR points.
!           "Recoverable" errors:
!              IERR = -1  if N<2 .
!              IERR = -2  if INCFD<1 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if NE<1 .
!           (Output arrays have not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!              IERR = -5  if an error has occurred in the lower-level
!                         routine DCHFDV.  NB: this should never happen.
!                         Notify the author **IMMEDIATELY** if it does.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DCHFDV, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   811020  DATE WRITTEN
!   820803  Minor cosmetic changes for release 1.
!   870707  Corrected XERROR calls for d.p. name(s).
!   890206  Corrected XERROR calls.
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891006  Cosmetic changes to prologue.  (WRB)
!   891006  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DPCHFD
!  Programming notes:
!
!     2. Most of the coding between the call to DCHFDV and the end of
!        the IR-loop could be eliminated if it were permissible to
!        assume that XE is ordered relative to X.
!
!     3. DCHFDV does not assume that X1 is less than X2.  thus, it would
!        be possible to write a version of DPCHFD that assumes a strict-
!        ly decreasing X-array by simply running the IR-loop backwards
!        (and reversing the order of appropriate tests).
!
!     4. The present code has a minor bug, which I have decided is not
!        worth the effort that would be required to fix it.
!        If XE contains points in [X(N-1),X(N)], followed by points <
!        X(N-1), followed by points >X(N), the extrapolation points
!        will be counted (at least) twice in the total returned in IERR.


integer  n, incfd, ne, ierr
real(wp)  x(*), f(incfd,*), d(incfd,*), xe(*), fe(*), &
 de(*)
logical  skip

integer  i, ierc, ir, j, jfirst, next(2), nj

!  VALIDITY-CHECK ARGUMENTS.

if (skip)  go to 5
!
if ( n<2 )  go to 5001
if ( incfd<1 )  go to 5002
do 1  i = 2, n
   if ( x(i)<=x(i-1) )  go to 5003
1 continue
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
5 continue
if ( ne<1 )  go to 5004
ierr = 0
skip = .true.
!
!  LOOP OVER INTERVALS.        (   INTERVAL INDEX IS  IL = IR-1  . )
!                              ( INTERVAL IS X(IL)<=X<X(IR) . )
jfirst = 1
ir = 2
10 continue
!
!     SKIP OUT OF LOOP IF HAVE PROCESSED ALL EVALUATION POINTS.
!
   if (jfirst > ne)  go to 5000
!
!     LOCATE ALL POINTS IN INTERVAL.
!
   do 20  j = jfirst, ne
      if (xe(j) >= x(ir))  go to 30
20    continue
   j = ne + 1
   go to 40
!
!     HAVE LOCATED FIRST POINT BEYOND INTERVAL.
!
30    continue
   if (ir == n)  j = ne + 1
!
40    continue
   nj = j - jfirst
!
!     SKIP EVALUATION IF NO POINTS IN INTERVAL.
!
   if (nj == 0)  go to 50
!
!     EVALUATE CUBIC AT XE(I),  I = JFIRST (1) J-1 .
!
!       ----------------------------------------------------------------
  call dchfdv (x(ir-1),x(ir), f(1,ir-1),f(1,ir), d(1,ir-1),d(1,ir) &
              ,nj, xe(jfirst), fe(jfirst), de(jfirst), next, ierc)
!       ----------------------------------------------------------------
   if (ierc < 0)  go to 5005
!
   if (next(2) == 0)  go to 42
!        IF (NEXT(2) > 0)  THEN
!           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(2) TO THE
!           RIGHT OF X(IR).
!
      if (ir < n)  go to 41
!           IF (IR == N)  THEN
!              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
         ierr = ierr + next(2)
         go to 42
41       continue
!           ELSE
!              WE SHOULD NEVER HAVE GOTTEN HERE.
         go to 5005
!           ENDIF
!        ENDIF
42    continue
!
   if (next(1) == 0)  go to 49
!        IF (NEXT(1) > 0)  THEN
!           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(1) TO THE
!           LEFT OF X(IR-1).
!
      if (ir > 2)  go to 43
!           IF (IR == 2)  THEN
!              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
         ierr = ierr + next(1)
         go to 49
43       continue
!           ELSE
!              XE IS NOT ORDERED RELATIVE TO X, SO MUST ADJUST
!              EVALUATION INTERVAL.
!
!              FIRST, LOCATE FIRST POINT TO LEFT OF X(IR-1).
         do 44  i = jfirst, j-1
            if (xe(i) < x(ir-1))  go to 45
44          continue
!              NOTE-- CANNOT DROP THROUGH HERE UNLESS THERE IS AN ERROR
!                     IN DCHFDV.
         go to 5005
!
45          continue
!              RESET J.  (THIS WILL BE THE NEW JFIRST.)
         j = i
!
!              NOW FIND OUT HOW FAR TO BACK UP IN THE X-ARRAY.
         do 46  i = 1, ir-1
            if (xe(j) < x(i)) go to 47
46          continue
!              NB-- CAN NEVER DROP THROUGH HERE, SINCE XE(J)<X(IR-1).
!
47          continue
!              AT THIS POINT, EITHER  XE(J) < X(1)
!                 OR      X(I-1) <= XE(J) < X(I) .
!              RESET IR, RECOGNIZING THAT IT WILL BE INCREMENTED BEFORE
!              CYCLING.
         ir = max(1, i-1)
!           ENDIF
!        ENDIF
49    continue
!
   jfirst = j
!
!     END OF IR-LOOP.
!
50 continue
ir = ir + 1
if (ir <= n)  go to 10
!
!  NORMAL RETURN.
!
5000 continue
return
!
!  ERROR RETURNS.
!
5001 continue
!     N<2 RETURN.
ierr = -1
call xermsg ('SLATEC', 'DPCHFD', &
   'NUMBER OF DATA POINTS LESS THAN TWO', ierr, 1)
return
!
5002 continue
!     INCFD<1 RETURN.
ierr = -2
call xermsg ('SLATEC', 'DPCHFD', 'INCREMENT LESS THAN ONE', ierr, &
   1)
return
!
5003 continue
!     X-ARRAY NOT STRICTLY INCREASING.
ierr = -3
call xermsg ('SLATEC', 'DPCHFD', &
   'X-ARRAY NOT STRICTLY INCREASING', ierr, 1)
return
!
5004 continue
!     NE<1 RETURN.
ierr = -4
call xermsg ('SLATEC', 'DPCHFD', &
   'NUMBER OF EVALUATION POINTS LESS THAN ONE', ierr, 1)
return
!
5005 continue
!     ERROR RETURN FROM DCHFDV.
!   *** THIS CASE SHOULD NEVER OCCUR ***
ierr = -5
call xermsg ('SLATEC', 'DPCHFD', &
   'ERROR RETURN FROM DCHFDV -- FATAL', ierr, 2)

end subroutine dpchfd

subroutine dpchfe (n, x, f, d, incfd, skip, ne, xe, fe, ierr)

!***PURPOSE  Evaluate a piecewise cubic Hermite function at an array of
!            points.  May be used by itself for Hermite interpolation,
!            or as an evaluator for DPCHIM or DPCHIC.
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E3
!***TYPE      real(wp) (PCHFE-S, DPCHFE-D)
!***KEYWORDS  CUBIC HERMITE EVALUATION, HERMITE INTERPOLATION, PCHIP,
!             PIECEWISE CUBIC EVALUATION
!***AUTHOR  Fritsch, F. N., (LLNL)
!             Lawrence Livermore National Laboratory
!             P.O. Box 808  (L-316)
!             Livermore, CA  94550
!             FTS 532-4275, (510) 422-4275
!***DESCRIPTION
!
!          DPCHFE:  Piecewise Cubic Hermite Function Evaluator
!
!     Evaluates the cubic Hermite function defined by  N, X, F, D  at
!     the points  XE(J), J=1(1)NE.
!
!     To provide compatibility with DPCHIM and DPCHIC, includes an
!     increment between successive values of the F- and D-arrays.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  N, NE, IERR
!        real(wp)  X(N), F(INCFD,N), D(INCFD,N), XE(NE), FE(NE)
!        LOGICAL  SKIP
!
!        CALL  DPCHFE (N, X, F, D, INCFD, SKIP, NE, XE, FE, IERR)
!
!   Parameters:
!
!     N -- (input) number of data points.  (Error return if N<2 .)
!
!     X -- (input) real*8 array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) < X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real*8 array of function values.  F(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     D -- (input) real*8 array of derivative values.  D(1+(I-1)*INCFD)
!           is the value corresponding to X(I).
!
!     INCFD -- (input) increment between successive values in F and D.
!           (Error return if  INCFD<1 .)
!
!     SKIP -- (input/output) logical variable which should be set to
!           .TRUE. if the user wishes to skip checks for validity of
!           preceding parameters, or to .FALSE. otherwise.
!           This will save time in case these checks have already
!           been performed (say, in DPCHIM or DPCHIC).
!           SKIP will be set to .TRUE. on normal return.
!
!     NE -- (input) number of evaluation points.  (Error return if
!           NE<1 .)
!
!     XE -- (input) real*8 array of points at which the function is to
!           be evaluated.
!
!          NOTES:
!           1. The evaluation will be most efficient if the elements
!              of XE are increasing relative to X;
!              that is,   XE(J) >= X(I)
!              implies    XE(K) >= X(I),  all K>=J .
!           2. If any of the XE are outside the interval [X(1),X(N)],
!              values are extrapolated from the nearest extreme cubic,
!              and a warning error is returned.
!
!     FE -- (output) real*8 array of values of the cubic Hermite
!           function defined by  N, X, F, D  at the points  XE.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR>0  means that extrapolation was performed at
!                 IERR points.
!           "Recoverable" errors:
!              IERR = -1  if N<2 .
!              IERR = -2  if INCFD<1 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if NE<1 .
!             (The FE-array has not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DCHFEV, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   811020  DATE WRITTEN
!   820803  Minor cosmetic changes for release 1.
!   870707  Corrected XERROR calls for d.p. name(s).
!   890206  Corrected XERROR calls.
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891006  Cosmetic changes to prologue.  (WRB)
!   891006  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DPCHFE
!  Programming notes:
!
!     2. Most of the coding between the call to DCHFEV and the end of
!        the IR-loop could be eliminated if it were permissible to
!        assume that XE is ordered relative to X.
!
!     3. DCHFEV does not assume that X1 is less than X2.  thus, it would
!        be possible to write a version of DPCHFE that assumes a
!        decreasing X-array by simply running the IR-loop backwards
!        (and reversing the order of appropriate tests).
!
!     4. The present code has a minor bug, which I have decided is not
!        worth the effort that would be required to fix it.
!        If XE contains points in [X(N-1),X(N)], followed by points <
!        X(N-1), followed by points >X(N), the extrapolation points
!        will be counted (at least) twice in the total returned in IERR.


integer  n, incfd, ne, ierr
real(wp)  x(*), f(incfd,*), d(incfd,*), xe(*), fe(*)
logical  skip

integer  i, ierc, ir, j, jfirst, next(2), nj

!  VALIDITY-CHECK ARGUMENTS.
if (skip)  go to 5
!
if ( n<2 )  go to 5001
if ( incfd<1 )  go to 5002
do 1  i = 2, n
   if ( x(i)<=x(i-1) )  go to 5003
1 continue
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
5 continue
if ( ne<1 )  go to 5004
ierr = 0
skip = .true.
!
!  LOOP OVER INTERVALS.        (   INTERVAL INDEX IS  IL = IR-1  . )
!                              ( INTERVAL IS X(IL)<=X<X(IR) . )
jfirst = 1
ir = 2
10 continue
!
!     SKIP OUT OF LOOP IF HAVE PROCESSED ALL EVALUATION POINTS.
!
   if (jfirst > ne)  go to 5000
!
!     LOCATE ALL POINTS IN INTERVAL.
!
   do 20  j = jfirst, ne
      if (xe(j) >= x(ir))  go to 30
20    continue
   j = ne + 1
   go to 40
!
!     HAVE LOCATED FIRST POINT BEYOND INTERVAL.
!
30    continue
   if (ir == n)  j = ne + 1
!
40    continue
   nj = j - jfirst
!
!     SKIP EVALUATION IF NO POINTS IN INTERVAL.
!
   if (nj == 0)  go to 50
!
!     EVALUATE CUBIC AT XE(I),  I = JFIRST (1) J-1 .
!
!       ----------------------------------------------------------------
  call dchfev (x(ir-1),x(ir), f(1,ir-1),f(1,ir), d(1,ir-1),d(1,ir) &
              ,nj, xe(jfirst), fe(jfirst), next, ierc)
!       ----------------------------------------------------------------
   if (ierc < 0)  go to 5005
!
   if (next(2) == 0)  go to 42
!        IF (NEXT(2) > 0)  THEN
!           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(2) TO THE
!           RIGHT OF X(IR).
!
      if (ir < n)  go to 41
!           IF (IR == N)  THEN
!              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
         ierr = ierr + next(2)
         go to 42
41       continue
!           ELSE
!              WE SHOULD NEVER HAVE GOTTEN HERE.
         go to 5005
!           ENDIF
!        ENDIF
42    continue
!
   if (next(1) == 0)  go to 49
!        IF (NEXT(1) > 0)  THEN
!           IN THE CURRENT SET OF XE-POINTS, THERE ARE NEXT(1) TO THE
!           LEFT OF X(IR-1).
!
      if (ir > 2)  go to 43
!           IF (IR == 2)  THEN
!              THESE ARE ACTUALLY EXTRAPOLATION POINTS.
         ierr = ierr + next(1)
         go to 49
43       continue
!           ELSE
!              XE IS NOT ORDERED RELATIVE TO X, SO MUST ADJUST
!              EVALUATION INTERVAL.
!
!              FIRST, LOCATE FIRST POINT TO LEFT OF X(IR-1).
         do 44  i = jfirst, j-1
            if (xe(i) < x(ir-1))  go to 45
44          continue
!              NOTE-- CANNOT DROP THROUGH HERE UNLESS THERE IS AN ERROR
!                     IN DCHFEV.
         go to 5005
!
45          continue
!              RESET J.  (THIS WILL BE THE NEW JFIRST.)
         j = i
!
!              NOW FIND OUT HOW FAR TO BACK UP IN THE X-ARRAY.
         do 46  i = 1, ir-1
            if (xe(j) < x(i)) go to 47
46          continue
!              NB-- CAN NEVER DROP THROUGH HERE, SINCE XE(J)<X(IR-1).
!
47          continue
!              AT THIS POINT, EITHER  XE(J) < X(1)
!                 OR      X(I-1) <= XE(J) < X(I) .
!              RESET IR, RECOGNIZING THAT IT WILL BE INCREMENTED BEFORE
!              CYCLING.
         ir = max(1, i-1)
!           ENDIF
!        ENDIF
49    continue
!
   jfirst = j
!
!     END OF IR-LOOP.
!
50 continue
ir = ir + 1
if (ir <= n)  go to 10
!
!  NORMAL RETURN.
!
5000 continue
return
!
!  ERROR RETURNS.
!
5001 continue
!     N<2 RETURN.
ierr = -1
call xermsg ('SLATEC', 'DPCHFE', &
   'NUMBER OF DATA POINTS LESS THAN TWO', ierr, 1)
return
!
5002 continue
!     INCFD<1 RETURN.
ierr = -2
call xermsg ('SLATEC', 'DPCHFE', 'INCREMENT LESS THAN ONE', ierr, &
   1)
return
!
5003 continue
!     X-ARRAY NOT STRICTLY INCREASING.
ierr = -3
call xermsg ('SLATEC', 'DPCHFE', &
   'X-ARRAY NOT STRICTLY INCREASING', ierr, 1)
return
!
5004 continue
!     NE<1 RETURN.
ierr = -4
call xermsg ('SLATEC', 'DPCHFE', &
   'NUMBER OF EVALUATION POINTS LESS THAN ONE', ierr, 1)
return
!
5005 continue
!     ERROR RETURN FROM DCHFEV.
!   *** THIS CASE SHOULD NEVER OCCUR ***
ierr = -5
call xermsg ('SLATEC', 'DPCHFE', &
   'ERROR RETURN FROM DCHFEV -- FATAL', ierr, 2)

end subroutine dpchfe

real(wp) function dpchia (n, x, f, d, incfd, skip, a, b, &
   ierr)

!***PURPOSE  Evaluate the definite integral of a piecewise cubic
!            Hermite function over an arbitrary interval.
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E3, H2A1B2
!***TYPE      real(wp) (PCHIA-S, DPCHIA-D)
!***KEYWORDS  CUBIC HERMITE INTERPOLATION, NUMERICAL INTEGRATION, PCHIP,
!             QUADRATURE
!***AUTHOR  Fritsch, F. N., (LLNL)
!             Lawrence Livermore National Laboratory
!             P.O. Box 808  (L-316)
!             Livermore, CA  94550
!             FTS 532-4275, (510) 422-4275
!***DESCRIPTION
!
!          DPCHIA:  Piecewise Cubic Hermite Integrator, Arbitrary limits
!
!     Evaluates the definite integral of the cubic Hermite function
!     defined by  N, X, F, D  over the interval [A, B].
!
!     To provide compatibility with DPCHIM and DPCHIC, includes an
!     increment between successive values of the F- and D-arrays.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  N, IERR
!        real(wp)  X(N), F(INCFD,N), D(INCFD,N), A, B
!        real(wp)  VALUE, DPCHIA
!        LOGICAL  SKIP
!
!        VALUE = DPCHIA (N, X, F, D, INCFD, SKIP, A, B, IERR)
!
!   Parameters:
!
!     VALUE -- (output) value of the requested integral.
!
!     N -- (input) number of data points.  (Error return if N<2 .)
!
!     X -- (input) real*8 array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) < X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real*8 array of function values.  F(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     D -- (input) real*8 array of derivative values.  D(1+(I-1)*INCFD)
!           is the value corresponding to X(I).
!
!     INCFD -- (input) increment between successive values in F and D.
!           (Error return if  INCFD<1 .)
!
!     SKIP -- (input/output) logical variable which should be set to
!           .TRUE. if the user wishes to skip checks for validity of
!           preceding parameters, or to .FALSE. otherwise.
!           This will save time in case these checks have already
!           been performed (say, in DPCHIM or DPCHIC).
!           SKIP will be set to .TRUE. on return with IERR>=0 .
!
!     A,B -- (input) the limits of integration.
!           NOTE:  There is no requirement that [A,B] be contained in
!                  [X(1),X(N)].  However, the resulting integral value
!                  will be highly suspect, if not.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning errors:
!              IERR = 1  if  A  is outside the interval [X(1),X(N)].
!              IERR = 2  if  B  is outside the interval [X(1),X(N)].
!              IERR = 3  if both of the above are true.  (Note that this
!                        means that either [A,B] contains data interval
!                        or the intervals do not intersect at all.)
!           "Recoverable" errors:
!              IERR = -1  if N<2 .
!              IERR = -2  if INCFD<1 .
!              IERR = -3  if the X-array is not strictly increasing.
!                (VALUE will be zero in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!              IERR = -4  in case of an error return from DPCHID (which
!                         should never occur).
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DCHFIE, DPCHID, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   820730  DATE WRITTEN
!   820804  Converted to SLATEC library version.
!   870707  Corrected XERROR calls for d.p. name(s).
!   870707  Corrected conversion to real(wp).
!   870813  Minor cosmetic changes.
!   890206  Corrected XERROR calls.
!   890411  Added SAVE statements (Vers. 3.2).
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890703  Corrected category record.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891006  Cosmetic changes to prologue.  (WRB)
!   891006  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   930503  Corrected to set VALUE=0 when IERR<0.  (FNF)
!   930504  Changed DCHFIV to DCHFIE.  (FNF)
!***END PROLOGUE  DPCHIA
!
!  Programming notes:
!  1. The error flag from DPCHID is tested, because a logic flaw
!     could conceivably result in IERD=-4, which should be reported.
!**End


integer  n, incfd, ierr
real(wp)  x(*), f(incfd,*), d(incfd,*), a, b
logical  skip

integer  i, ia, ib, ierd, il, ir
real(wp)  value, xa, xb

value = zero
!
!  VALIDITY-CHECK ARGUMENTS.
!
if (skip)  go to 5
!
if ( n<2 )  go to 5001
if ( incfd<1 )  go to 5002
do 1  i = 2, n
   if ( x(i)<=x(i-1) )  go to 5003
1 continue
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
5 continue
skip = .true.
ierr = 0
if ( (a<x(1)) .or. (a>x(n)) )  ierr = ierr + 1
if ( (b<x(1)) .or. (b>x(n)) )  ierr = ierr + 2
!
!  COMPUTE INTEGRAL VALUE.
!
if (a /= b)  then
   xa = min (a, b)
   xb = max (a, b)
   if (xb <= x(2))  then
!           INTERVAL IS TO LEFT OF X(2), SO USE FIRST CUBIC.
!                   ---------------------------------------
      value = dchfie (x(1),x(2), f(1,1),f(1,2), &
                                 d(1,1),d(1,2), a, b)
!                   ---------------------------------------
   else if (xa >= x(n-1))  then
!           INTERVAL IS TO RIGHT OF X(N-1), SO USE LAST CUBIC.
!                   ------------------------------------------
      value = dchfie(x(n-1),x(n), f(1,n-1),f(1,n), &
                                  d(1,n-1),d(1,n), a, b)
!                   ------------------------------------------
   else
!           'NORMAL' CASE -- XA<XB, XA<X(N-1), XB>X(2).
!      ......LOCATE IA AND IB SUCH THAT
!               X(IA-1)<XA<=X(IA)<=X(IB)<=XB<=X(IB+1)
      ia = 1
      do 10  i = 1, n-1
         if (xa > x(i))  ia = i + 1
10       continue
!             IA = 1 IMPLIES XA<X(1) .  OTHERWISE,
!             IA IS LARGEST INDEX SUCH THAT X(IA-1)<XA,.
!
      ib = n
      do 20  i = n, ia, -1
         if (xb < x(i))  ib = i - 1
20       continue
!             IB = N IMPLIES XB>X(N) .  OTHERWISE,
!             IB IS SMALLEST INDEX SUCH THAT XB<X(IB+1) .
!
!     ......COMPUTE THE INTEGRAL.
      if (ib < ia)  then
!              THIS MEANS IB = IA-1 AND
!                 (A,B) IS A SUBSET OF (X(IB),X(IA)).
!                      -------------------------------------------
         value = dchfie (x(ib),x(ia), f(1,ib),f(1,ia), &
                                      d(1,ib),d(1,ia), a, b)
!                      -------------------------------------------
      else
!
!              FIRST COMPUTE INTEGRAL OVER (X(IA),X(IB)).
!                (Case (IB == IA) is taken care of by initialization
!                 of VALUE to ZERO.)
         if (ib > ia)  then
!                         ---------------------------------------------
            value = dpchid (n, x, f, d, incfd, skip, ia, ib, ierd)
!                         ---------------------------------------------
            if (ierd < 0)  go to 5004
         endif
!
!              THEN ADD ON INTEGRAL OVER (XA,X(IA)).
         if (xa < x(ia))  then
            il = max(1, ia-1)
            ir = il + 1
!                                 -------------------------------------
            value = value + dchfie (x(il),x(ir), f(1,il),f(1,ir), &
                                      d(1,il),d(1,ir), xa, x(ia))
!                                 -------------------------------------
         endif
!
!              THEN ADD ON INTEGRAL OVER (X(IB),XB).
         if (xb > x(ib))  then
            ir = min (ib+1, n)
            il = ir - 1
!                                 -------------------------------------
            value = value + dchfie (x(il),x(ir), f(1,il),f(1,ir), &
                                      d(1,il),d(1,ir), x(ib), xb)
!                                 -------------------------------------
         endif
!
!              FINALLY, ADJUST SIGN IF NECESSARY.
         if (a > b)  value = -value
      endif
   endif
endif
!
!  NORMAL RETURN.
!
5000 continue
dpchia = value
return
!
!  ERROR RETURNS.
!
5001 continue
!     N<2 RETURN.
ierr = -1
call xermsg ('SLATEC', 'DPCHIA', &
   'NUMBER OF DATA POINTS LESS THAN TWO', ierr, 1)
go to 5000
!
5002 continue
!     INCFD<1 RETURN.
ierr = -2
call xermsg ('SLATEC', 'DPCHIA', 'INCREMENT LESS THAN ONE', ierr, &
   1)
go to 5000
!
5003 continue
!     X-ARRAY NOT STRICTLY INCREASING.
ierr = -3
call xermsg ('SLATEC', 'DPCHIA', &
   'X-ARRAY NOT STRICTLY INCREASING', ierr, 1)
go to 5000
!
5004 continue
!     TROUBLE IN DPCHID.  (SHOULD NEVER OCCUR.)
ierr = -4
call xermsg ('SLATEC', 'DPCHIA', 'TROUBLE IN DPCHID', ierr, 1)
go to 5000
!------------- LAST LINE OF DPCHIA FOLLOWS -----------------------------
end function dpchia

subroutine dpchic (ic, vc, switch, n, x, f, d, incfd, wk, nwk, &
   ierr)

!***PURPOSE  Set derivatives needed to determine a piecewise monotone
!            piecewise cubic Hermite interpolant to given data.
!            User control is available over boundary conditions and/or
!            treatment of points where monotonicity switches direction.
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E1A
!***TYPE      real(wp) (PCHIC-S, DPCHIC-D)
!***KEYWORDS  CUBIC HERMITE INTERPOLATION, MONOTONE INTERPOLATION,
!             PCHIP, PIECEWISE CUBIC INTERPOLATION,
!             SHAPE-PRESERVING INTERPOLATION
!***AUTHOR  Fritsch, F. N., (LLNL)
!             Lawrence Livermore National Laboratory
!             P.O. Box 808  (L-316)
!             Livermore, CA  94550
!             FTS 532-4275, (510) 422-4275
!***DESCRIPTION
!
!         DPCHIC:  Piecewise Cubic Hermite Interpolation Coefficients.
!
!     Sets derivatives needed to determine a piecewise monotone piece-
!     wise cubic interpolant to the data given in X and F satisfying the
!     boundary conditions specified by IC and VC.
!
!     The treatment of points where monotonicity switches direction is
!     controlled by argument SWITCH.
!
!     To facilitate two-dimensional applications, includes an increment
!     between successive values of the F- and D-arrays.
!
!     The resulting piecewise cubic Hermite function may be evaluated
!     by DPCHFE or DPCHFD.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  IC(2), N, NWK, IERR
!        real(wp)  VC(2), SWITCH, X(N), F(INCFD,N), D(INCFD,N),
!                          WK(NWK)
!
!        CALL DPCHIC (IC, VC, SWITCH, N, X, F, D, INCFD, WK, NWK, IERR)
!
!   Parameters:
!
!     IC -- (input) integer array of length 2 specifying desired
!           boundary conditions:
!           IC(1) = IBEG, desired condition at beginning of data.
!           IC(2) = IEND, desired condition at end of data.
!
!           IBEG = 0  for the default boundary condition (the same as
!                     used by DPCHIM).
!           If IBEG/=0, then its sign indicates whether the boundary
!                     derivative is to be adjusted, if necessary, to be
!                     compatible with monotonicity:
!              IBEG>0  if no adjustment is to be performed.
!              IBEG<0  if the derivative is to be adjusted for
!                     monotonicity.
!
!           Allowable values for the magnitude of IBEG are:
!           IBEG = 1  if first derivative at X(1) is given in VC(1).
!           IBEG = 2  if second derivative at X(1) is given in VC(1).
!           IBEG = 3  to use the 3-point difference formula for D(1).
!                     (Reverts to the default b.c. if N<3 .)
!           IBEG = 4  to use the 4-point difference formula for D(1).
!                     (Reverts to the default b.c. if N<4 .)
!           IBEG = 5  to set D(1) so that the second derivative is con-
!              tinuous at X(2). (Reverts to the default b.c. if N<4.)
!              This option is somewhat analogous to the "not a knot"
!              boundary condition provided by DPCHSP.
!
!          NOTES (IBEG):
!           1. An error return is taken if ABS(IBEG)>5 .
!           2. Only in case  IBEG<=0  is it guaranteed that the
!              interpolant will be monotonic in the first interval.
!              If the returned value of D(1) lies between zero and
!              3*SLOPE(1), the interpolant will be monotonic.  This
!              is **NOT** checked if IBEG>0 .
!           3. If IBEG<0 and D(1) had to be changed to achieve mono-
!              tonicity, a warning error is returned.
!
!           IEND may take on the same values as IBEG, but applied to
!           derivative at X(N).  In case IEND = 1 or 2, the value is
!           given in VC(2).
!
!          NOTES (IEND):
!           1. An error return is taken if ABS(IEND)>5 .
!           2. Only in case  IEND<=0  is it guaranteed that the
!              interpolant will be monotonic in the last interval.
!              If the returned value of D(1+(N-1)*INCFD) lies between
!              zero and 3*SLOPE(N-1), the interpolant will be monotonic.
!              This is **NOT** checked if IEND>0 .
!           3. If IEND<0 and D(1+(N-1)*INCFD) had to be changed to
!              achieve monotonicity, a warning error is returned.
!
!     VC -- (input) real*8 array of length 2 specifying desired boundary
!           values, as indicated above.
!           VC(1) need be set only if IC(1) = 1 or 2 .
!           VC(2) need be set only if IC(2) = 1 or 2 .
!
!     SWITCH -- (input) indicates desired treatment of points where
!           direction of monotonicity switches:
!           Set SWITCH to zero if interpolant is required to be mono-
!           tonic in each interval, regardless of monotonicity of data.
!             NOTES:
!              1. This will cause D to be set to zero at all switch
!                 points, thus forcing extrema there.
!              2. The result of using this option with the default boun-
!                 dary conditions will be identical to using DPCHIM, but
!                 will generally cost more compute time.
!                 This option is provided only to facilitate comparison
!                 of different switch and/or boundary conditions.
!           Set SWITCH nonzero to use a formula based on the 3-point
!              difference formula in the vicinity of switch points.
!           If SWITCH is positive, the interpolant on each interval
!              containing an extremum is controlled to not deviate from
!              the data by more than SWITCH*DFLOC, where DFLOC is the
!              maximum of the change of F on this interval and its two
!              immediate neighbors.
!           If SWITCH is negative, no such control is to be imposed.
!
!     N -- (input) number of data points.  (Error return if N<2 .)
!
!     X -- (input) real*8 array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) < X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real*8 array of dependent variable values to be
!           interpolated.  F(1+(I-1)*INCFD) is value corresponding to
!           X(I).
!
!     D -- (output) real*8 array of derivative values at the data
!           points.  These values will determine a monotone cubic
!           Hermite function on each subinterval on which the data
!           are monotonic, except possibly adjacent to switches in
!           monotonicity. The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in F and D.
!           This argument is provided primarily for 2-D applications.
!           (Error return if  INCFD<1 .)
!
!     WK -- (scratch) real*8 array of working storage.  The user may
!           wish to know that the returned values are:
!              WK(I)     = H(I)     = X(I+1) - X(I) ;
!              WK(N-1+I) = SLOPE(I) = (F(1,I+1) - F(1,I)) / H(I)
!           for  I = 1(1)N-1.
!
!     NWK -- (input) length of work array.
!           (Error return if  NWK<2*(N-1) .)
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning errors:
!              IERR = 1  if IBEG<0 and D(1) had to be adjusted for
!                        monotonicity.
!              IERR = 2  if IEND<0 and D(1+(N-1)*INCFD) had to be
!                        adjusted for monotonicity.
!              IERR = 3  if both of the above are true.
!           "Recoverable" errors:
!              IERR = -1  if N<2 .
!              IERR = -2  if INCFD<1 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if ABS(IBEG)>5 .
!              IERR = -5  if ABS(IEND)>5 .
!              IERR = -6  if both of the above are true.
!              IERR = -7  if NWK<2*(N-1) .
!             (The D-array has not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
!***REFERENCES  1. F. N. Fritsch, Piecewise Cubic Hermite Interpolation
!                 Package, Report UCRL-87285, Lawrence Livermore Natio-
!                 nal Laboratory, July 1982.  [Poster presented at the
!                 SIAM 30th Anniversary Meeting, 19-23 July 1982.]
!               2. F. N. Fritsch and J. Butland, A method for construc-
!                 ting local monotone piecewise cubic interpolants, SIAM
!                 Journal on Scientific and Statistical Computing 5, 2
!                 (June 1984), pp. 300-304.
!               3. F. N. Fritsch and R. E. Carlson, Monotone piecewise
!                 cubic interpolation, SIAM Journal on Numerical Ana-
!                 lysis 17, 2 (April 1980), pp. 238-246.
!***ROUTINES CALLED  DPCHCE, DPCHCI, DPCHCS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   820218  DATE WRITTEN
!   820804  Converted to SLATEC library version.
!   870707  Corrected XERROR calls for d.p. name(s).
!   870813  Updated Reference 2.
!   890206  Corrected XERROR calls.
!   890411  Added SAVE statements (Vers. 3.2).
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890703  Corrected category record.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891006  Cosmetic changes to prologue.  (WRB)
!   891006  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920429  Revised format and order of references.  (WRB,FNF)
!***END PROLOGUE  DPCHIC

integer  ic(2), n, incfd, nwk, ierr
real(wp)  vc(2), switch, x(*), f(incfd,*), d(incfd,*), &
 wk(nwk)

integer  i, ibeg, iend, nless1

!  VALIDITY-CHECK ARGUMENTS.
if ( n<2 )  go to 5001
if ( incfd<1 )  go to 5002
do 1  i = 2, n
   if ( x(i)<=x(i-1) )  go to 5003
1 continue
!
ibeg = ic(1)
iend = ic(2)
ierr = 0
if (abs(ibeg) > 5)  ierr = ierr - 1
if (abs(iend) > 5)  ierr = ierr - 2
if (ierr < 0)  go to 5004
!
!  FUNCTION DEFINITION IS OK -- GO ON.
!
nless1 = n - 1
if ( nwk < 2*nless1 )  go to 5007
!
!  SET UP H AND SLOPE ARRAYS.
!
do 20  i = 1, nless1
   wk(i) = x(i+1) - x(i)
   wk(nless1+i) = (f(1,i+1) - f(1,i)) / wk(i)
20 continue
!
!  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.
!
if (nless1 > 1)  go to 1000
d(1,1) = wk(2)
d(1,n) = wk(2)
go to 3000
!
!  NORMAL CASE  (N >= 3) .
!
1000 continue
!
!  SET INTERIOR DERIVATIVES AND DEFAULT END CONDITIONS.
!
!     --------------------------------------
call dpchci (n, wk(1), wk(n), d, incfd)
!     --------------------------------------
!
!  SET DERIVATIVES AT POINTS WHERE MONOTONICITY SWITCHES DIRECTION.
!
if (switch == zero)  go to 3000
!     ----------------------------------------------------
call dpchcs (switch, n, wk(1), wk(n), d, incfd, ierr)
!     ----------------------------------------------------
if (ierr /= 0)  go to 5008
!
!  SET END CONDITIONS.
!
3000 continue
if ( (ibeg==0) .and. (iend==0) )  go to 5000
!     -------------------------------------------------------
call dpchce (ic, vc, n, x, wk(1), wk(n), d, incfd, ierr)
!     -------------------------------------------------------
if (ierr < 0)  go to 5009
!
!  NORMAL RETURN.
!
5000 continue
return
!
!  ERROR RETURNS.
!
5001 continue
!     N<2 RETURN.
ierr = -1
call xermsg ('SLATEC', 'DPCHIC', &
   'NUMBER OF DATA POINTS LESS THAN TWO', ierr, 1)
return
!
5002 continue
!     INCFD<1 RETURN.
ierr = -2
call xermsg ('SLATEC', 'DPCHIC', 'INCREMENT LESS THAN ONE', ierr, &
   1)
return
!
5003 continue
!     X-ARRAY NOT STRICTLY INCREASING.
ierr = -3
call xermsg ('SLATEC', 'DPCHIC', &
   'X-ARRAY NOT STRICTLY INCREASING', ierr, 1)
return
!
5004 continue
!     IC OUT OF RANGE RETURN.
ierr = ierr - 3
call xermsg ('SLATEC', 'DPCHIC', 'IC OUT OF RANGE', ierr, 1)
return
!
5007 continue
!     NWK < 2*(N-1)  RETURN.
ierr = -7
call xermsg ('SLATEC', 'DPCHIC', 'WORK ARRAY TOO SMALL', ierr, 1)
return
!
5008 continue
!     ERROR RETURN FROM DPCHCS.
ierr = -8
call xermsg ('SLATEC', 'DPCHIC', 'ERROR RETURN FROM DPCHCS', &
   ierr, 1)
return
!
5009 continue
!     ERROR RETURN FROM DPCHCE.
!   *** THIS CASE SHOULD NEVER OCCUR ***
ierr = -9
call xermsg ('SLATEC', 'DPCHIC', 'ERROR RETURN FROM DPCHCE', &
   ierr, 1)

end subroutine dpchic

real(wp) function dpchid (n, x, f, d, incfd, skip, ia, ib, &
   ierr)

!***PURPOSE  Evaluate the definite integral of a piecewise cubic
!            Hermite function over an interval whose endpoints are data
!            points.
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E3, H2A1B2
!***TYPE      real(wp) (PCHID-S, DPCHID-D)
!***KEYWORDS  CUBIC HERMITE INTERPOLATION, NUMERICAL INTEGRATION, PCHIP,
!             QUADRATURE
!***AUTHOR  Fritsch, F. N., (LLNL)
!             Lawrence Livermore National Laboratory
!             P.O. Box 808  (L-316)
!             Livermore, CA  94550
!             FTS 532-4275, (510) 422-4275
!***DESCRIPTION
!
!          DPCHID:  Piecewise Cubic Hermite Integrator, Data limits
!
!     Evaluates the definite integral of the cubic Hermite function
!     defined by  N, X, F, D  over the interval [X(IA), X(IB)].
!
!     To provide compatibility with DPCHIM and DPCHIC, includes an
!     increment between successive values of the F- and D-arrays.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  N, IA, IB, IERR
!        real(wp)  X(N), F(INCFD,N), D(INCFD,N)
!        LOGICAL  SKIP
!
!        VALUE = DPCHID (N, X, F, D, INCFD, SKIP, IA, IB, IERR)
!
!   Parameters:
!
!     VALUE -- (output) value of the requested integral.
!
!     N -- (input) number of data points.  (Error return if N<2 .)
!
!     X -- (input) real*8 array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) < X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real*8 array of function values.  F(1+(I-1)*INCFD) is
!           the value corresponding to X(I).
!
!     D -- (input) real*8 array of derivative values.  D(1+(I-1)*INCFD)
!           is the value corresponding to X(I).
!
!     INCFD -- (input) increment between successive values in F and D.
!           (Error return if  INCFD<1 .)
!
!     SKIP -- (input/output) logical variable which should be set to
!           .TRUE. if the user wishes to skip checks for validity of
!           preceding parameters, or to .FALSE. otherwise.
!           This will save time in case these checks have already
!           been performed (say, in DPCHIM or DPCHIC).
!           SKIP will be set to .TRUE. on return with IERR = 0 or -4.
!
!     IA,IB -- (input) indices in X-array for the limits of integration.
!           both must be in the range [1,N].  (Error return if not.)
!           No restrictions on their relative values.
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  if N<2 .
!              IERR = -2  if INCFD<1 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if IA or IB is out of range.
!                (VALUE will be zero in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   820723  DATE WRITTEN
!   820804  Converted to SLATEC library version.
!   870707  Corrected XERROR calls for d.p. name(s).
!   870813  Minor cosmetic changes.
!   890206  Corrected XERROR calls.
!   890411  Added SAVE statements (Vers. 3.2).
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890703  Corrected category record.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891006  Cosmetic changes to prologue.  (WRB)
!   891006  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   930504  Corrected to set VALUE=0 when IERR/=0.  (FNF)
!***END PROLOGUE  DPCHID
!
!  Programming notes:
!  1. This routine uses a special formula that is valid only for
!     integrals whose limits coincide with data values.  This is
!     mathematically equivalent to, but much more efficient than,
!     calls to DCHFIE.
!**End


integer  n, incfd, ia, ib, ierr
real(wp)  x(*), f(incfd,*), d(incfd,*)
logical  skip

integer  i, iup, low
real(wp)  h, sum, value

value = zero
!
!  VALIDITY-CHECK ARGUMENTS.
!
if (skip)  go to 5
!
if ( n<2 )  go to 5001
if ( incfd<1 )  go to 5002
do 1  i = 2, n
   if ( x(i)<=x(i-1) )  go to 5003
1 continue
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
5 continue
skip = .true.
if ((ia<1) .or. (ia>n))  go to 5004
if ((ib<1) .or. (ib>n))  go to 5004
ierr = 0
!
!  COMPUTE INTEGRAL VALUE.
!
if (ia /= ib)  then
   low = min(ia, ib)
   iup = max(ia, ib) - 1
   sum = zero
   do 10  i = low, iup
      h = x(i+1) - x(i)
      sum = sum + h*( (f(1,i) + f(1,i+1)) + &
                      (d(1,i) - d(1,i+1))*(h/six) )
10    continue
   value = half * sum
   if (ia > ib)  value = -value
endif
!
!  NORMAL RETURN.
!
5000 continue
dpchid = value
return
!
!  ERROR RETURNS.
!
5001 continue
!     N<2 RETURN.
ierr = -1
call xermsg ('SLATEC', 'DPCHID', &
   'NUMBER OF DATA POINTS LESS THAN TWO', ierr, 1)
go to 5000
!
5002 continue
!     INCFD<1 RETURN.
ierr = -2
call xermsg ('SLATEC', 'DPCHID', 'INCREMENT LESS THAN ONE', ierr, &
   1)
go to 5000
!
5003 continue
!     X-ARRAY NOT STRICTLY INCREASING.
ierr = -3
call xermsg ('SLATEC', 'DPCHID', &
   'X-ARRAY NOT STRICTLY INCREASING', ierr, 1)
go to 5000
!
5004 continue
!     IA OR IB OUT OF RANGE RETURN.
ierr = -4
call xermsg ('SLATEC', 'DPCHID', 'IA OR IB OUT OF RANGE', ierr, &
   1)
go to 5000
!------------- LAST LINE OF DPCHID FOLLOWS -----------------------------
end function dpchid

subroutine dpchim (n, x, f, d, incfd, ierr)

!***PURPOSE  Set derivatives needed to determine a monotone piecewise
!            cubic Hermite interpolant to given data.  Boundary values
!            are provided which are compatible with monotonicity.  The
!            interpolant will have an extremum at each point where mono-
!            tonicity switches direction.  (See DPCHIC if user control
!            is desired over boundary or switch conditions.)
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E1A
!***TYPE      real(wp) (PCHIM-S, DPCHIM-D)
!***KEYWORDS  CUBIC HERMITE INTERPOLATION, MONOTONE INTERPOLATION,
!             PCHIP, PIECEWISE CUBIC INTERPOLATION
!***AUTHOR  Fritsch, F. N., (LLNL)
!             Lawrence Livermore National Laboratory
!             P.O. Box 808  (L-316)
!             Livermore, CA  94550
!             FTS 532-4275, (510) 422-4275
!***DESCRIPTION
!
!          DPCHIM:  Piecewise Cubic Hermite Interpolation to
!                  Monotone data.
!
!     Sets derivatives needed to determine a monotone piecewise cubic
!     Hermite interpolant to the data given in X and F.
!
!     Default boundary conditions are provided which are compatible
!     with monotonicity.  (See DPCHIC if user control of boundary con-
!     ditions is desired.)
!
!     If the data are only piecewise monotonic, the interpolant will
!     have an extremum at each point where monotonicity switches direc-
!     tion.  (See DPCHIC if user control is desired in such cases.)
!
!     To facilitate two-dimensional applications, includes an increment
!     between successive values of the F- and D-arrays.
!
!     The resulting piecewise cubic Hermite function may be evaluated
!     by DPCHFE or DPCHFD.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  N, IERR
!        real(wp)  X(N), F(INCFD,N), D(INCFD,N)
!
!        CALL  DPCHIM (N, X, F, D, INCFD, IERR)
!
!   Parameters:
!
!     N -- (input) number of data points.  (Error return if N<2 .)
!           If N=2, simply does linear interpolation.
!
!     X -- (input) real*8 array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) < X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real*8 array of dependent variable values to be
!           interpolated.  F(1+(I-1)*INCFD) is value corresponding to
!           X(I).  DPCHIM is designed for monotonic data, but it will
!           work for any F-array.  It will force extrema at points where
!           monotonicity switches direction.  If some other treatment of
!           switch points is desired, DPCHIC should be used instead.
!                                     -----
!     D -- (output) real*8 array of derivative values at the data
!           points.  If the data are monotonic, these values will
!           determine a monotone cubic Hermite function.
!           The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in F and D.
!           This argument is provided primarily for 2-D applications.
!           (Error return if  INCFD<1 .)
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           Warning error:
!              IERR>0  means that IERR switches in the direction
!                 of monotonicity were detected.
!           "Recoverable" errors:
!              IERR = -1  if N<2 .
!              IERR = -2  if INCFD<1 .
!              IERR = -3  if the X-array is not strictly increasing.
!             (The D-array has not been changed in any of these cases.)
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!
!***REFERENCES  1. F. N. Fritsch and J. Butland, A method for construc-
!                 ting local monotone piecewise cubic interpolants, SIAM
!                 Journal on Scientific and Statistical Computing 5, 2
!                 (June 1984), pp. 300-304.
!               2. F. N. Fritsch and R. E. Carlson, Monotone piecewise
!                 cubic interpolation, SIAM Journal on Numerical Ana-
!                 lysis 17, 2 (April 1980), pp. 238-246.
!***ROUTINES CALLED  DPCHST, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   811103  DATE WRITTEN
!   820201  1. Introduced  DPCHST  to reduce possible over/under-
!             flow problems.
!           2. Rearranged derivative formula for same reason.
!   820602  1. Modified end conditions to be continuous functions
!             of data when monotonicity switches in next interval.
!           2. Modified formulas so end conditions are less prone
!             of over/underflow problems.
!   820803  Minor cosmetic changes for release 1.
!   870707  Corrected XERROR calls for d.p. name(s).
!   870813  Updated Reference 1.
!   890206  Corrected XERROR calls.
!   890411  Added SAVE statements (Vers. 3.2).
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890703  Corrected category record.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891006  Cosmetic changes to prologue.  (WRB)
!   891006  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920429  Revised format and order of references.  (WRB,FNF)
!***END PROLOGUE  DPCHIM
!  Programming notes:
!
!     1. The function  DPCHST(ARG1,ARG2)  is assumed to return zero if
!        either argument is zero, +1 if they are of the same sign, and
!        -1 if they are of opposite sign.

integer  n, incfd, ierr
real(wp)  x(*), f(incfd,*), d(incfd,*)

integer  i, nless1
real(wp)  del1, del2, dmax, dmin, drat1, drat2, dsave, &
      h1, h2, hsum, hsumt3, w1, w2

!  VALIDITY-CHECK ARGUMENTS.
if ( n<2 )  go to 5001
if ( incfd<1 )  go to 5002
do 1  i = 2, n
   if ( x(i)<=x(i-1) )  go to 5003
1 continue
!
!  FUNCTION DEFINITION IS OK, GO ON.
!
ierr = 0
nless1 = n - 1
h1 = x(2) - x(1)
del1 = (f(1,2) - f(1,1))/h1
dsave = del1
!
!  SPECIAL CASE N=2 -- USE LINEAR INTERPOLATION.
!
if (nless1 > 1)  go to 10
d(1,1) = del1
d(1,n) = del1
go to 5000
!
!  NORMAL CASE  (N >= 3).
!
10 continue
h2 = x(3) - x(2)
del2 = (f(1,3) - f(1,2))/h2
!
!  SET D(1) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!     SHAPE-PRESERVING.
!
hsum = h1 + h2
w1 = (h1 + hsum)/hsum
w2 = -h1/hsum
d(1,1) = w1*del1 + w2*del2
if ( dpchst(d(1,1),del1) <= zero)  then
   d(1,1) = zero
else if ( dpchst(del1,del2) < zero)  then
!        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
   dmax = three*del1
   if (abs(d(1,1)) > abs(dmax))  d(1,1) = dmax
endif
!
!  LOOP THROUGH INTERIOR POINTS.
!
do 50  i = 2, nless1
   if (i == 2)  go to 40
!
   h1 = h2
   h2 = x(i+1) - x(i)
   hsum = h1 + h2
   del1 = del2
   del2 = (f(1,i+1) - f(1,i))/h2
40    continue
!
!        SET D(I)=0 UNLESS DATA ARE STRICTLY MONOTONIC.
!
   d(1,i) = zero
   if ( dpchst(del1,del2) )  42, 41, 45
!
!        COUNT NUMBER OF CHANGES IN DIRECTION OF MONOTONICITY.
!
41    continue
   if (del2 == zero)  go to 50
   if ( dpchst(dsave,del2) < zero)  ierr = ierr + 1
   dsave = del2
   go to 50
!
42    continue
   ierr = ierr + 1
   dsave = del2
   go to 50
!
!        USE BRODLIE MODIFICATION OF BUTLAND FORMULA.
!
45    continue
   hsumt3 = hsum+hsum+hsum
   w1 = (hsum + h1)/hsumt3
   w2 = (hsum + h2)/hsumt3
   dmax = max( abs(del1), abs(del2) )
   dmin = min( abs(del1), abs(del2) )
   drat1 = del1/dmax
   drat2 = del2/dmax
   d(1,i) = dmin/(w1*drat1 + w2*drat2)
!
50 continue
!
!  SET D(N) VIA NON-CENTERED THREE-POINT FORMULA, ADJUSTED TO BE
!     SHAPE-PRESERVING.
!
w1 = -h2/hsum
w2 = (h2 + hsum)/hsum
d(1,n) = w1*del1 + w2*del2
if ( dpchst(d(1,n),del2) <= zero)  then
   d(1,n) = zero
else if ( dpchst(del1,del2) < zero)  then
!        NEED DO THIS CHECK ONLY IF MONOTONICITY SWITCHES.
   dmax = three*del2
   if (abs(d(1,n)) > abs(dmax))  d(1,n) = dmax
endif
!
!  NORMAL RETURN.
!
5000 continue
return
!
!  ERROR RETURNS.
!
5001 continue
!     N<2 RETURN.
ierr = -1
call xermsg ('SLATEC', 'DPCHIM', &
   'NUMBER OF DATA POINTS LESS THAN TWO', ierr, 1)
return
!
5002 continue
!     INCFD<1 RETURN.
ierr = -2
call xermsg ('SLATEC', 'DPCHIM', 'INCREMENT LESS THAN ONE', ierr, &
   1)
return
!
5003 continue
!     X-ARRAY NOT STRICTLY INCREASING.
ierr = -3
call xermsg ('SLATEC', 'DPCHIM', &
   'X-ARRAY NOT STRICTLY INCREASING', ierr, 1)

end subroutine dpchim

subroutine dpchkt (n, x, knotyp, t)


!***PURPOSE  Compute B-spline knot sequence for DPCHBS.
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E3
!***TYPE      real(wp) (PCHKT-S, DPCHKT-D)
!***AUTHOR  Fritsch, F. N., (LLNL)
!***DESCRIPTION
!
!     Set a knot sequence for the B-spline representation of a PCH
!     function with breakpoints X.  All knots will be at least double.
!     Endknots are set as:
!        (1) quadruple knots at endpoints if KNOTYP=0;
!        (2) extrapolate the length of end interval if KNOTYP=1;
!        (3) periodic if KNOTYP=2.
!
!  Input arguments:  N, X, KNOTYP.
!  Output arguments:  T.
!
!  Restrictions/assumptions:
!     1. N>=2 .  (not checked)
!     2. X(i)<X(i+1), i=1,...,N .  (not checked)
!     3. 0<=KNOTYP<=2 .  (Acts like KNOTYP=0 for any other value.)
!
!***SEE ALSO  DPCHBS
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   870701  DATE WRITTEN
!   900405  Converted Fortran to upper case.
!   900410  Converted prologue to SLATEC 4.0 format.
!   900410  Minor cosmetic changes.
!   900430  Produced real(wp) version.
!   930514  Changed NKNOTS from an output to an input variable.  (FNF)
!   930604  Removed unused variable NKNOTS from argument list.  (FNF)
!***END PROLOGUE  DPCHKT
!
!*Internal Notes:
!
!  Since this is subsidiary to DPCHBS, which validates its input before
!  calling, it is unnecessary for such validation to be done here.
!
!**End


integer  n, knotyp
real(wp)  x(*), t(*)

integer  j, k, ndim
real(wp)  hbeg, hend

!  Initialize.
ndim = 2*n
!
!  Set interior knots.
!
j = 1
do 20  k = 1, n
   j = j + 2
   t(j) = x(k)
   t(j+1) = t(j)
20 continue
!     Assertion:  At this point T(3),...,T(NDIM+2) have been set and
!                 J=NDIM+1.
!
!  Set end knots according to KNOTYP.
!
hbeg = x(2) - x(1)
hend = x(n) - x(n-1)
if (knotyp==1 )  then
!          Extrapolate.
   t(2) = x(1) - hbeg
   t(ndim+3) = x(n) + hend
else if ( knotyp==2 )  then
!          Periodic.
   t(2) = x(1) - hend
   t(ndim+3) = x(n) + hbeg
else
!          Quadruple end knots.
   t(2) = x(1)
   t(ndim+3) = x(n)
endif
t(1) = t(2)
t(ndim+4) = t(ndim+3)
!
!  Terminate.
!

end subroutine dpchkt

subroutine dpchsp (ic, vc, n, x, f, d, incfd, wk, nwk, ierr)

!***PURPOSE  Set derivatives needed to determine the Hermite represen-
!            tation of the cubic spline interpolant to given data, with
!            specified boundary conditions.
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E1A
!***TYPE      real(wp) (PCHSP-S, DPCHSP-D)
!***KEYWORDS  CUBIC HERMITE INTERPOLATION, PCHIP,
!             PIECEWISE CUBIC INTERPOLATION, SPLINE INTERPOLATION
!***AUTHOR  Fritsch, F. N., (LLNL)
!             Lawrence Livermore National Laboratory
!             P.O. Box 808  (L-316)
!             Livermore, CA  94550
!             FTS 532-4275, (510) 422-4275
!***DESCRIPTION
!
!          DPCHSP:   Piecewise Cubic Hermite Spline
!
!     Computes the Hermite representation of the cubic spline inter-
!     polant to the data given in X and F satisfying the boundary
!     conditions specified by IC and VC.
!
!     To facilitate two-dimensional applications, includes an increment
!     between successive values of the F- and D-arrays.
!
!     The resulting piecewise cubic Hermite function may be evaluated
!     by DPCHFE or DPCHFD.
!
!     NOTE:  This is a modified version of C. de Boor's cubic spline
!            routine CUBSPL.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  IC(2), N, NWK, IERR
!        real(wp)  VC(2), X(N), F(INCFD,N), D(INCFD,N), WK(NWK)
!
!        CALL  DPCHSP (IC, VC, N, X, F, D, INCFD, WK, NWK, IERR)
!
!   Parameters:
!
!     IC -- (input) integer array of length 2 specifying desired
!           boundary conditions:
!           IC(1) = IBEG, desired condition at beginning of data.
!           IC(2) = IEND, desired condition at end of data.
!
!           IBEG = 0  to set D(1) so that the third derivative is con-
!              tinuous at X(2).  This is the "not a knot" condition
!              provided by de Boor's cubic spline routine CUBSPL.
!              < This is the default boundary condition. >
!           IBEG = 1  if first derivative at X(1) is given in VC(1).
!           IBEG = 2  if second derivative at X(1) is given in VC(1).
!           IBEG = 3  to use the 3-point difference formula for D(1).
!                     (Reverts to the default b.c. if N<3 .)
!           IBEG = 4  to use the 4-point difference formula for D(1).
!                     (Reverts to the default b.c. if N<4 .)
!          NOTES:
!           1. An error return is taken if IBEG is out of range.
!           2. For the "natural" boundary condition, use IBEG=2 and
!              VC(1)=0.
!
!           IEND may take on the same values as IBEG, but applied to
!           derivative at X(N).  In case IEND = 1 or 2, the value is
!           given in VC(2).
!
!          NOTES:
!           1. An error return is taken if IEND is out of range.
!           2. For the "natural" boundary condition, use IEND=2 and
!              VC(2)=0.
!
!     VC -- (input) real*8 array of length 2 specifying desired boundary
!           values, as indicated above.
!           VC(1) need be set only if IC(1) = 1 or 2 .
!           VC(2) need be set only if IC(2) = 1 or 2 .
!
!     N -- (input) number of data points.  (Error return if N<2 .)
!
!     X -- (input) real*8 array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1) < X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real*8 array of dependent variable values to be
!           interpolated.  F(1+(I-1)*INCFD) is value corresponding to
!           X(I).
!
!     D -- (output) real*8 array of derivative values at the data
!           points.  These values will determine the cubic spline
!           interpolant with the requested boundary conditions.
!           The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in F and D.
!           This argument is provided primarily for 2-D applications.
!           (Error return if  INCFD<1 .)
!
!     WK -- (scratch) real*8 array of working storage.
!
!     NWK -- (input) length of work array.
!           (Error return if NWK<2*N .)
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  if N<2 .
!              IERR = -2  if INCFD<1 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if IBEG<0 or IBEG>4 .
!              IERR = -5  if IEND<0 of IEND>4 .
!              IERR = -6  if both of the above are true.
!              IERR = -7  if NWK is too small.
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!             (The D-array has not been changed in any of these cases.)
!              IERR = -8  in case of trouble solving the linear system
!                         for the interior derivative values.
!             (The D-array may have been changed in this case.)
!             (             Do **NOT** use it!                )
!
!***REFERENCES  Carl de Boor, A Practical Guide to Splines, Springer-
!                 Verlag, New York, 1978, pp. 53-59.
!***ROUTINES CALLED  DPCHDF, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   820503  DATE WRITTEN
!   820804  Converted to SLATEC library version.
!   870707  Corrected XERROR calls for d.p. name(s).
!   890206  Corrected XERROR calls.
!   890411  Added SAVE statements (Vers. 3.2).
!   890703  Corrected category record.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891006  Cosmetic changes to prologue.  (WRB)
!   891006  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920429  Revised format and order of references.  (WRB,FNF)
!***END PROLOGUE  DPCHSP

integer  ic(2), n, incfd, nwk, ierr
real(wp)  vc(2), x(*), f(incfd,*), d(incfd,*), wk(2,*)

integer  ibeg, iend, index, j, nm1
real(wp)  g, stemp(3), xtemp(4)

!  VALIDITY-CHECK ARGUMENTS.
if ( n<2 )  go to 5001
if ( incfd<1 )  go to 5002
do 1  j = 2, n
   if ( x(j)<=x(j-1) )  go to 5003
1 continue
!
ibeg = ic(1)
iend = ic(2)
ierr = 0
if ( (ibeg<0).or.(ibeg>4) )  ierr = ierr - 1
if ( (iend<0).or.(iend>4) )  ierr = ierr - 2
if ( ierr<0 )  go to 5004
!
!  FUNCTION DEFINITION IS OK -- GO ON.
!
if ( nwk < 2*n )  go to 5007
!
!  COMPUTE FIRST DIFFERENCES OF X SEQUENCE AND STORE IN WK(1,.). ALSO,
!  COMPUTE FIRST DIVIDED DIFFERENCE OF DATA AND STORE IN WK(2,.).
do 5  j=2,n
   wk(1,j) = x(j) - x(j-1)
   wk(2,j) = (f(1,j) - f(1,j-1))/wk(1,j)
5 continue
!
!  SET TO DEFAULT BOUNDARY CONDITIONS IF N IS TOO SMALL.
!
if ( ibeg>n )  ibeg = 0
if ( iend>n )  iend = 0
!
!  SET UP FOR BOUNDARY CONDITIONS.
!
if ( (ibeg==1).or.(ibeg==2) )  then
   d(1,1) = vc(1)
else if (ibeg > 2)  then
!        PICK UP FIRST IBEG POINTS, IN REVERSE ORDER.
   do 10  j = 1, ibeg
      index = ibeg-j+1
!           INDEX RUNS FROM IBEG DOWN TO 1.
      xtemp(j) = x(index)
      if (j < ibeg)  stemp(j) = wk(2,index)
10    continue
!                 --------------------------------
   d(1,1) = dpchdf (ibeg, xtemp, stemp, ierr)
!                 --------------------------------
   if (ierr /= 0)  go to 5009
   ibeg = 1
endif
!
if ( (iend==1).or.(iend==2) )  then
   d(1,n) = vc(2)
else if (iend > 2)  then
!        PICK UP LAST IEND POINTS.
   do 15  j = 1, iend
      index = n-iend+j
!           INDEX RUNS FROM N+1-IEND UP TO N.
      xtemp(j) = x(index)
      if (j < iend)  stemp(j) = wk(2,index+1)
15    continue
!                 --------------------------------
   d(1,n) = dpchdf (iend, xtemp, stemp, ierr)
!                 --------------------------------
   if (ierr /= 0)  go to 5009
   iend = 1
endif
!
! --------------------( BEGIN CODING FROM CUBSPL )--------------------
!
!  **** A TRIDIAGONAL LINEAR SYSTEM FOR THE UNKNOWN SLOPES S(J) OF
!  F  AT X(J), J=1,...,N, IS GENERATED AND THEN SOLVED BY GAUSS ELIM-
!  INATION, WITH S(J) ENDING UP IN D(1,J), ALL J.
!     WK(1,.) AND WK(2,.) ARE USED FOR TEMPORARY STORAGE.
!
!  CONSTRUCT FIRST EQUATION FROM FIRST BOUNDARY CONDITION, OF THE FORM
!             WK(2,1)*S(1) + WK(1,1)*S(2) = D(1,1)
!
if (ibeg == 0)  then
   if (n == 2)  then
!           NO CONDITION AT LEFT END AND N = 2.
      wk(2,1) = one
      wk(1,1) = one
      d(1,1) = two*wk(2,2)
   else
!           NOT-A-KNOT CONDITION AT LEFT END AND N > 2.
      wk(2,1) = wk(1,3)
      wk(1,1) = wk(1,2) + wk(1,3)
      d(1,1) =((wk(1,2) + two*wk(1,1))*wk(2,2)*wk(1,3) &
                        + wk(1,2)**2*wk(2,3)) / wk(1,1)
   endif
else if (ibeg == 1)  then
!        SLOPE PRESCRIBED AT LEFT END.
   wk(2,1) = one
   wk(1,1) = zero
else
!        SECOND DERIVATIVE PRESCRIBED AT LEFT END.
   wk(2,1) = two
   wk(1,1) = one
   d(1,1) = three*wk(2,2) - half*wk(1,2)*d(1,1)
endif
!
!  IF THERE ARE INTERIOR KNOTS, GENERATE THE CORRESPONDING EQUATIONS AND
!  CARRY OUT THE FORWARD PASS OF GAUSS ELIMINATION, AFTER WHICH THE J-TH
!  EQUATION READS    WK(2,J)*S(J) + WK(1,J)*S(J+1) = D(1,J).
!
nm1 = n-1
if (nm1 > 1)  then
   do 20 j=2,nm1
      if (wk(2,j-1) == zero)  go to 5008
      g = -wk(1,j+1)/wk(2,j-1)
      d(1,j) = g*d(1,j-1) &
                  + three*(wk(1,j)*wk(2,j+1) + wk(1,j+1)*wk(2,j))
      wk(2,j) = g*wk(1,j-1) + two*(wk(1,j) + wk(1,j+1))
20    continue
endif
!
!  CONSTRUCT LAST EQUATION FROM SECOND BOUNDARY CONDITION, OF THE FORM
!           (-G*WK(2,N-1))*S(N-1) + WK(2,N)*S(N) = D(1,N)
!
!     IF SLOPE IS PRESCRIBED AT RIGHT END, ONE CAN GO DIRECTLY TO BACK-
!     SUBSTITUTION, SINCE ARRAYS HAPPEN TO BE SET UP JUST RIGHT FOR IT
!     AT THIS POINT.
if (iend == 1)  go to 30
!
if (iend == 0)  then
   if (n==2 .and. ibeg==0)  then
!           NOT-A-KNOT AT RIGHT ENDPOINT AND AT LEFT ENDPOINT AND N = 2.
      d(1,2) = wk(2,2)
      go to 30
   else if ((n==2) .or. (n==3 .and. ibeg==0))  then
!           EITHER (N=3 AND NOT-A-KNOT ALSO AT LEFT) OR (N=2 AND *NOT*
!           NOT-A-KNOT AT LEFT END POINT).
      d(1,n) = two*wk(2,n)
      wk(2,n) = one
      if (wk(2,n-1) == zero)  go to 5008
      g = -one/wk(2,n-1)
   else
!           NOT-A-KNOT AND N >= 3, AND EITHER N>3 OR  ALSO NOT-A-
!           KNOT AT LEFT END POINT.
      g = wk(1,n-1) + wk(1,n)
!           DO NOT NEED TO CHECK FOLLOWING DENOMINATORS (X-DIFFERENCES).
      d(1,n) = ((wk(1,n)+two*g)*wk(2,n)*wk(1,n-1) &
                  + wk(1,n)**2*(f(1,n-1)-f(1,n-2))/wk(1,n-1))/g
      if (wk(2,n-1) == zero)  go to 5008
      g = -g/wk(2,n-1)
      wk(2,n) = wk(1,n-1)
   endif
else
!        SECOND DERIVATIVE PRESCRIBED AT RIGHT ENDPOINT.
   d(1,n) = three*wk(2,n) + half*wk(1,n)*d(1,n)
   wk(2,n) = two
   if (wk(2,n-1) == zero)  go to 5008
   g = -one/wk(2,n-1)
endif
!
!  COMPLETE FORWARD PASS OF GAUSS ELIMINATION.
!
wk(2,n) = g*wk(1,n-1) + wk(2,n)
if (wk(2,n) == zero)   go to 5008
d(1,n) = (g*d(1,n-1) + d(1,n))/wk(2,n)
!
!  CARRY OUT BACK SUBSTITUTION
!
30 continue
do 40 j=nm1,1,-1
   if (wk(2,j) == zero)  go to 5008
   d(1,j) = (d(1,j) - wk(1,j)*d(1,j+1))/wk(2,j)
40 continue
! --------------------(  END  CODING FROM CUBSPL )--------------------
!
!  NORMAL RETURN.
!
return
!
!  ERROR RETURNS.
!
5001 continue
!     N<2 RETURN.
ierr = -1
call xermsg ('SLATEC', 'DPCHSP', &
   'NUMBER OF DATA POINTS LESS THAN TWO', ierr, 1)
return
!
5002 continue
!     INCFD<1 RETURN.
ierr = -2
call xermsg ('SLATEC', 'DPCHSP', 'INCREMENT LESS THAN ONE', ierr, &
   1)
return
!
5003 continue
!     X-ARRAY NOT STRICTLY INCREASING.
ierr = -3
call xermsg ('SLATEC', 'DPCHSP', &
   'X-ARRAY NOT STRICTLY INCREASING', ierr, 1)
return
!
5004 continue
!     IC OUT OF RANGE RETURN.
ierr = ierr - 3
call xermsg ('SLATEC', 'DPCHSP', 'IC OUT OF RANGE', ierr, 1)
return
!
5007 continue
!     NWK TOO SMALL RETURN.
ierr = -7
call xermsg ('SLATEC', 'DPCHSP', 'WORK ARRAY TOO SMALL', ierr, 1)
return
!
5008 continue
!     SINGULAR SYSTEM.
!   *** THEORETICALLY, THIS CAN ONLY OCCUR IF SUCCESSIVE X-VALUES   ***
!   *** ARE EQUAL, WHICH SHOULD ALREADY HAVE BEEN CAUGHT (IERR=-3). ***
ierr = -8
call xermsg ('SLATEC', 'DPCHSP', 'SINGULAR LINEAR SYSTEM', ierr, &
   1)
return
!
5009 continue
!     ERROR RETURN FROM DPCHDF.
!   *** THIS CASE SHOULD NEVER OCCUR ***
ierr = -9
call xermsg ('SLATEC', 'DPCHSP', 'ERROR RETURN FROM DPCHDF', &
   ierr, 1)

end subroutine dpchsp

!***PURPOSE  DPCHIP Sign-Testing Routine
!***LIBRARY   SLATEC (PCHIP)
!***TYPE      real(wp) (PCHST-S, DPCHST-D)
!***AUTHOR  Fritsch, F. N., (LLNL)
!***DESCRIPTION
!
!         DPCHST:  DPCHIP Sign-Testing Routine.
!
!
!     Returns:
!        -1. if ARG1 and ARG2 are of opposite sign.
!         0. if either argument is zero.
!        +1. if ARG1 and ARG2 are of the same sign.
!
!     The object is to do this without multiplying ARG1*ARG2, to avoid
!     possible over/underflow problems.
!
!***SEE ALSO  DPCHCE, DPCHCI, DPCHCS, DPCHIM
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   811103  DATE WRITTEN
!   820805  Converted to SLATEC library version.
!   870813  Minor cosmetic changes.
!   890411  Added SAVE statements (Vers. 3.2).
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910408  Updated AUTHOR and DATE WRITTEN sections in prologue.  (WRB)
!   930503  Improved purpose.  (FNF)
!***END PROLOGUE  DPCHST
!
!**End

real(wp) function dpchst (arg1, arg2)

real(wp)  arg1, arg2

!  PERFORM THE TEST.
dpchst = sign(one,arg1) * sign(one,arg2)
if ((arg1==zero) .or. (arg2==zero))  dpchst = zero

end function dpchst

subroutine dpchsw (dfmax, iextrm, d1, d2, h, slope, ierr)

!***PURPOSE  Limits excursion from data for DPCHCS
!***LIBRARY   SLATEC (PCHIP)
!***TYPE      real(wp) (PCHSW-S, DPCHSW-D)
!***AUTHOR  Fritsch, F. N., (LLNL)
!***DESCRIPTION
!
!         DPCHSW:  DPCHCS Switch Excursion Limiter.
!
!     Called by  DPCHCS  to adjust D1 and D2 if necessary to insure that
!     the extremum on this interval is not further than DFMAX from the
!     extreme data value.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        INTEGER  IEXTRM, IERR
!        real(wp)  DFMAX, D1, D2, H, SLOPE
!
!        CALL  DPCHSW (DFMAX, IEXTRM, D1, D2, H, SLOPE, IERR)
!
!   Parameters:
!
!     DFMAX -- (input) maximum allowed difference between F(IEXTRM) and
!           the cubic determined by derivative values D1,D2.  (assumes
!           DFMAX>0.)
!
!     IEXTRM -- (input) index of the extreme data value.  (assumes
!           IEXTRM = 1 or 2 .  Any value /=1 is treated as 2.)
!
!     D1,D2 -- (input) derivative values at the ends of the interval.
!           (Assumes D1*D2 <= 0.)
!          (output) may be modified if necessary to meet the restriction
!           imposed by DFMAX.
!
!     H -- (input) interval length.  (Assumes  H>0.)
!
!     SLOPE -- (input) data slope on the interval.
!
!     IERR -- (output) error flag.  should be zero.
!           If IERR=-1, assumption on D1 and D2 is not satisfied.
!           If IERR=-2, quadratic equation locating extremum has
!                       negative discriminant (should never occur).
!
!    -------
!    WARNING:  This routine does no validity-checking of arguments.
!    -------
!
!***SEE ALSO  DPCHCS
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   820218  DATE WRITTEN
!   820805  Converted to SLATEC library version.
!   870707  Corrected XERROR calls for d.p. name(s).
!   870707  Replaced DATA statement for SMALL with a use of D1MACH.
!   870813  Minor cosmetic changes.
!   890206  Corrected XERROR calls.
!   890411  1. Added SAVE statements (Vers. 3.2).
!           2. Added real(wp) declaration for D1MACH.
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900328  Added TYPE section.  (WRB)
!   910408  Updated AUTHOR and DATE WRITTEN sections in prologue.  (WRB)
!   920526  Eliminated possible divide by zero problem.  (FNF)
!   930503  Improved purpose.  (FNF)
!***END PROLOGUE  DPCHSW
!
!**End

integer  iextrm, ierr
real(wp)  dfmax, d1, d2, h, slope

real(wp)  cp, hphi, lambda, nu, phi, radcal, &
                  rho, sigma, that

real(wp),parameter :: fact = 100.0_wp
real(wp),parameter :: third = one/three - d1mach4  !! third should be slightly less than 1/3 (original code had 0.33333)
real(wp),parameter :: small = fact*d1mach4 !! small should be a few orders of magnitude greater than macheps.

!
!  NOTATION AND GENERAL REMARKS.
!
!     RHO IS THE RATIO OF THE DATA SLOPE TO THE DERIVATIVE BEING TESTED.
!     LAMBDA IS THE RATIO OF D2 TO D1.
!     THAT = T-HAT(RHO) IS THE NORMALIZED LOCATION OF THE EXTREMUM.
!     PHI IS THE NORMALIZED VALUE OF P(X)-F1 AT X = XHAT = X-HAT(RHO),
!           WHERE  THAT = (XHAT - X1)/H .
!        THAT IS, P(XHAT)-F1 = D*H*PHI,  WHERE D=D1 OR D2.
!     SIMILARLY,  P(XHAT)-F2 = D*H*(PHI-RHO) .

!
!  DO MAIN CALCULATION.
!
if (d1 == zero)  then
!
!        SPECIAL CASE -- D1==ZERO .
!
!          IF D2 IS ALSO ZERO, THIS ROUTINE SHOULD NOT HAVE BEEN CALLED.
   if (d2 == zero)  go to 5001
!
   rho = slope/d2
!          EXTREMUM IS OUTSIDE INTERVAL WHEN RHO >= 1/3 .
   if (rho >= third)  go to 5000
   that = (two*(three*rho-one)) / (three*(two*rho-one))
   phi = that**2 * ((three*rho-one)/three)
!
!          CONVERT TO DISTANCE FROM F2 IF IEXTRM/=1 .
   if (iextrm /= 1)  phi = phi - rho
!
!          TEST FOR EXCEEDING LIMIT, AND ADJUST ACCORDINGLY.
   hphi = h * abs(phi)
   if (hphi*abs(d2) > dfmax)  then
!           AT THIS POINT, HPHI>0, SO DIVIDE IS OK.
      d2 = sign (dfmax/hphi, d2)
   endif
else
!
   rho = slope/d1
   lambda = -d2/d1
   if (d2 == zero)  then
!
!           SPECIAL CASE -- D2==ZERO .
!
!             EXTREMUM IS OUTSIDE INTERVAL WHEN RHO >= 1/3 .
      if (rho >= third)  go to 5000
      cp = two - three*rho
      nu = one - two*rho
      that = one / (three*nu)
   else
      if (lambda <= zero)  go to 5001
!
!           NORMAL CASE -- D1 AND D2 BOTH NONZERO, OPPOSITE SIGNS.
!
      nu = one - lambda - two*rho
      sigma = one - rho
      cp = nu + sigma
      if (abs(nu) > small)  then
         radcal = (nu - (two*rho+one))*nu + sigma**2
         if (radcal < zero)  go to 5002
         that = (cp - sqrt(radcal)) / (three*nu)
      else
         that = one/(two*sigma)
      endif
   endif
   phi = that*((nu*that - cp)*that + one)
!
!          CONVERT TO DISTANCE FROM F2 IF IEXTRM/=1 .
   if (iextrm /= 1)  phi = phi - rho
!
!          TEST FOR EXCEEDING LIMIT, AND ADJUST ACCORDINGLY.
   hphi = h * abs(phi)
   if (hphi*abs(d1) > dfmax)  then
!           AT THIS POINT, HPHI>0, SO DIVIDE IS OK.
      d1 = sign (dfmax/hphi, d1)
      d2 = -lambda*d1
   endif
endif
!
!  NORMAL RETURN.
!
5000 continue
ierr = 0
return
!
!  ERROR RETURNS.
!
5001 continue
!     D1 AND D2 BOTH ZERO, OR BOTH NONZERO AND SAME SIGN.
ierr = -1
call xermsg ('SLATEC', 'DPCHSW', 'D1 AND/OR D2 INVALID', ierr, 1)
return
!
5002 continue
!     NEGATIVE VALUE OF RADICAL (SHOULD NEVER OCCUR).
ierr = -2
call xermsg ('SLATEC', 'DPCHSW', 'NEGATIVE RADICAL', ierr, 1)

end subroutine dpchsw

!*******************************************************************************************************
   end module pchip_module
!*******************************************************************************************************
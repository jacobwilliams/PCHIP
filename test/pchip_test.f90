!*******************************************************************************************************
!>
!  Tests for PCHIP.

program pchip_test

   use pchip_module
   use iso_fortran_env, only: Lun => output_unit, wp => real64

   implicit none

   integer, parameter :: Kprint = 5 !! printing flag
   integer :: Ipass !! pass/fail flag

   real(wp), parameter :: d1mach4 = epsilon(1.0_wp)  !! `d1mach(4)` -- the largest relative spacing

   ! run all the tests:

   call dpchq1(Lun, Kprint, Ipass); if (ipass == 0) error stop 'test dpchq1 failed'
   call dpchq2(Lun, Kprint, Ipass); if (ipass == 0) error stop 'test dpchq2 failed'
   call dpchq3(Lun, Kprint, Ipass); if (ipass == 0) error stop 'test dpchq3 failed'
   call dpchq4(Lun, Kprint, Ipass); if (ipass == 0) error stop 'test dpchq4 failed'
   call dpchq5(Lun, Kprint, Ipass); if (ipass == 0) error stop 'test dpchq5 failed'

contains
!*******************************************************************************************************

   !***PURPOSE  Test the PCHIP evaluators DCHFDV, DCHFEV, DPCHFD, DPCHFE.
   !***AUTHOR  Fritsch, F. N., (LLNL)
   !***DESCRIPTION
   !
   !             DPCHIP QUICK CHECK NUMBER 1
   !
   !     TESTS THE EVALUATORS:  DCHFDV, DCHFEV, DPCHFD, DPCHFE.
   ! *Usage:
   !
   !        INTEGER  LUN, KPRINT, IPASS
   !
   !        CALL DPCHQ1 (LUN, KPRINT, IPASS)
   !
   ! *Arguments:
   !
   !     LUN   :IN  is the unit number to which output is to be written.
   !
   !     KPRINT:IN  controls the amount of output, as specified in the
   !                SLATEC Guidelines.
   !
   !     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
   !                IPASS=0 indicates one or more tests failed.
   !
   ! *Description:
   !
   !   This routine carries out three tests of the PCH evaluators:
   !     DEVCHK tests the single-cubic evaluators.
   !     DEVPCK tests the full PCH evaluators.
   !     DEVERK exercises the error returns in all evaluators.
   !
   !***REVISION HISTORY  (YYMMDD)
   !   820601  DATE WRITTEN
   !   890306  Changed IPASS to the more accurate name IFAIL.  (FNF)
   !   890307  Removed conditional on call to DEVERK.
   !   890706  Cosmetic changes to prologue.  (WRB)
   !   891004  Correction in prologue.  (WRB)
   !   891214  Prologue converted to Version 4.0 format.  (BAB)
   !   900309  Added DEVERK to list of routines called.  (FNF)
   !   900314  Improved some output formats.
   !   900315  Revised prologue and improved some output formats.  (FNF)
   !   900316  Additional minor cosmetic changes.  (FNF)
   !   900321  Removed IFAIL from call sequence for SLATEC standards and
   !           made miscellaneous cosmetic changes.  (FNF)
   !   930317  Improved output formats.  (FNF)
   subroutine dpchq1(Lun, Kprint, Ipass)
      implicit none

      !
      !  Declare arguments.
      !
      integer Lun, Kprint, Ipass
      !
      !  DECLARE LOCAL VARIABLES.
      !
      integer i1, i2, i3, i4, i5, i6, i7, i8, i9, ifail, npts
      double precision work(4000)
      logical fail

      if (Kprint >= 2) write (Lun, 99001) Kprint
      !
      !  FORMATS.
      !
99001 format('1'/' ------------ DPCHIP QUICK CHECK OUTPUT',            &
               &' ------------'//20x, '( KPRINT =', i2, ' )')
      !
      !  TEST DCHFDV AND DCHFEV.
      !
      ifail = 0
      npts = 1000
      i1 = 1 + npts
      i2 = i1 + npts
      i3 = i2 + npts
      call devchk(Lun, Kprint, npts, work(1), work(i1), work(i2), work(i3),   &
                & fail)
      if (fail) ifail = ifail + 1
      !
      !  TEST DPCHFD AND DPCHFE.
      !
      i1 = 1 + 10
      i2 = i1 + 10
      i3 = i2 + 100
      i4 = i3 + 100
      i5 = i4 + 100
      i6 = i5 + 51
      i7 = i6 + 51
      i8 = i7 + 51
      i9 = i8 + 51
      call devpck(Lun, Kprint, work(1), work(i1), work(i2), work(i3), work(i4)&
                & , work(i5), work(i6), work(i7), work(i8), work(i9), fail)
      if (fail) ifail = ifail + 2
      !
      !  TEST ERROR RETURNS.
      !
      call deverk(Lun, Kprint, fail)
      if (fail) ifail = ifail + 4
      !
      !  PRINT SUMMARY AND TERMINATE.
      !     At this point, IFAIL has the following value:
      !        IFAIL = 0  IF ALL TESTS PASSED.
      !        IFAIL BETWEEN 1 AND 7 IS THE SUM OF:
      !           IFAIL=1  IF SINGLE CUBIC  TEST FAILED. (SEE DEVCHK OUTPUT.)
      !           IFAIL=2  IF DPCHFD/DPCHFE TEST FAILED. (SEE DEVPCK OUTPUT.)
      !           IFAIL=4  IF ERROR RETURN  TEST FAILED. (SEE DEVERK OUTPUT.)
      !
      if ((Kprint >= 2) .and. (ifail /= 0)) write (Lun, 99002) ifail
99002 format(/' *** TROUBLE ***', i5, ' EVALUATION TESTS FAILED.')
      !
      if (ifail == 0) then
         Ipass = 1
         if (Kprint >= 2) write (Lun, 99003)
99003    format(/' ------------ DPCHIP PASSED  ALL EVALUATION TESTS',  &
                     &' ------------')
      else
         Ipass = 0
         if (Kprint >= 1) write (Lun, 99004)
99004    format(/' ************ DPCHIP FAILED SOME EVALUATION TESTS',  &
                     &' ************')
      end if

   end subroutine dpchq1

   subroutine dpchq2(Lun, Kprint, Ipass)
      implicit none
      !***PURPOSE  Test the PCHIP integrators DPCHIA and DPCHID.
      !***AUTHOR  Fritsch, F. N., (LLNL)
      !***DESCRIPTION
      !
      !             DPCHIP QUICK CHECK NUMBER 2
      !
      !     TESTS THE INTEGRATORS:  DPCHIA, DPCHID.
      ! *Usage:
      !
      !        INTEGER  LUN, KPRINT, IPASS
      !
      !        CALL DPCHQ2 (LUN, KPRINT, IPASS)
      !
      ! *Arguments:
      !
      !     LUN   :IN  is the unit number to which output is to be written.
      !
      !     KPRINT:IN  controls the amount of output, as specified in the
      !                SLATEC Guidelines.
      !
      !     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
      !                IPASS=0 indicates one or more tests failed.
      !
      ! *Description:
      !
      !   This routine constructs data from a cubic, integrates it with DPCHIA
      !   and compares the results with the correct answer.
      !   Since DPCHIA calls DPCHID, this tests both integrators.
      !
      !***REVISION HISTORY  (YYMMDD)
      !   820601  DATE WRITTEN
      !   890306  Changed IPASS to the more accurate name IFAIL.  (FNF)
      !   890316  1. Removed IMPLICIT statement.                  (FNF)
      !           2. Eliminated unnecessary variable N1.          (FNF)
      !           3. Miscellaneous cosmetic changes.              (FNF)
      !   891004  Cosmetic changes to prologue.  (WRB)
      !   891214  Prologue converted to Version 4.0 format.  (BAB)
      !   900314  Improved some output formats.  (FNF)
      !   900315  Revised prologue and improved some output formats.  (FNF)
      !   900316  Additional minor cosmetic changes.  (FNF)
      !   900321  Removed IFAIL from call sequence for SLATEC standards and
      !           made miscellaneous cosmetic changes.  (FNF)
      !   900323  Corrected list of routines called.  (FNF)
      !   901130  Added 1P's to formats; changed to allow KPRINT.gt.3.  (FNF)
      !   910708  Minor modifications in use of KPRINT.  (WRB)
      !   930317  Improved output formats.  (FNF)

      !
      !  Declare arguments.
      !
      integer Lun, Kprint, Ipass
      !
      !  DECLARE VARIABLES.
      !
      integer i, ierexp(17), ierr, ifail, n, npairs
      double precision a(17), b(17), calc, d(7), errmax, error,   &
                     & f(7), machep, one, three, thrqtr, tol,     &
                     & true, two, x(7)
      logical fail, skip

      !
      !  DEFINE TEST FUNCTIONS.
      !
      double precision ax, fcn, deriv, antder
      fcn(ax) = three*ax*ax*(ax - two)
      deriv(ax) = three*ax*(two*(ax - two) + ax)
      antder(ax) = ax**3*(thrqtr*ax - two)
      !
      !  INITIALIZE.
      !
      data thrqtr/0.75d0/, one/1.d0/, two/2.d0/, three/3.d0/
      data n/7/
      data x/-4.d0, -2.d0, -0.9d0, 0.d0, 0.9d0, 2.d0, 4.d0/
      data npairs/17/
      data a/-3.0d0, 3.0d0, -0.5d0, -0.5d0, -0.5d0, -4.0d0,     &
         & -4.0d0, 3.0d0, -5.0d0, -5.0d0, -6.0d0, 6.0d0, -1.5d0, &
         & -1.5d0, -3.0d0, 3.0d0, 0.5d0/
      data b/3.0d0, -3.0d0, 1.0d0, 2.0d0, 5.0d0, -0.5d0, 4.0d0,  &
         & 5.0d0, -3.0d0, 5.0d0, -5.0d0, 5.0d0, -0.5d0, -1.0d0,  &
         & -2.5d0, 3.5d0, 0.5d0/
      data ierexp/0, 0, 0, 0, 2, 0, 0, 2, 1, 3, 3, 3, 0,   &
         & 0, 0, 0, 0/
      !
      !  SET PASS/FAIL TOLERANCE.
      !
      machep = d1mach4
      tol = 100.d0*machep
      !
      !  SET UP PCH FUNCTION DEFINITION.
      !
      do i = 1, n
         f(i) = fcn(x(i))
         d(i) = deriv(x(i))
      end do
      !
      if (Kprint >= 3) write (Lun, 99001)
      !
      !  FORMATS.
      !
99001 format('1'//10x, 'TEST DPCHIP INTEGRATORS')
      if (Kprint >= 2) write (Lun, 99002)
99002 format(//10x, 'DPCHQ2 RESULTS'/10x, '--------------')
      if (Kprint >= 3) write (Lun, 99003) (x(i), f(i), d(i), i=1, n)
99003 format(//5x, 'DATA:'//11x, 'X', 9x, 'F', 9x, 'D'/(5x, 3f10.3))
      !
      !  LOOP OVER (A,B)-PAIRS.
      !
      if (Kprint >= 3) write (Lun, 99004)
99004 format(//5x, 'TEST RESULTS:'//'    A     B    ERR     TRUE', 16x,  &
               &'CALC', 15x, 'ERROR')
      !
      ifail = 0
      !
      skip = .false.
      do i = 1, npairs
         !               ---------------------------------------------
         calc = dpchia(n, x, f, d, 1, skip, a(i), b(i), ierr)
         !               ---------------------------------------------
         if (ierr >= 0) then
            fail = ierr /= ierexp(i)
            true = antder(b(i)) - antder(a(i))
            error = calc - true
            if (Kprint >= 3) then
               if (fail) then
                  write (Lun, 99005) a(i), b(i), ierr, true, calc,  &
                                  & error, ierexp(i)
99005             format(2f6.1, i5, 1p, 2d20.10, d15.5, '  (', i1, ') *****')
               else
                  write (Lun, 99010) a(i), b(i), ierr, true, calc,  &
                                  & error
               end if
            end if
            !
            error = abs(error)/max(one, abs(true))
            if (fail .or. (error > tol)) ifail = ifail + 1
            if (i == 1) then
               errmax = error
            else
               errmax = max(errmax, error)
            end if
         else
            if (Kprint >= 3) write (Lun, 99010) a(i), b(i), ierr
            ifail = ifail + 1
         end if
      end do
      !
      !  PRINT SUMMARY.
      !
      if (Kprint >= 2) then
         write (Lun, 99006) errmax, tol
99006    format(/'  MAXIMUM RELATIVE ERROR IS:', 1p, d15.5,              &
                     &',   TOLERANCE:', 1p, d15.5)
         if (ifail /= 0) write (Lun, 99007) ifail
99007    format(/' *** TROUBLE ***', i5, ' INTEGRATION TESTS FAILED.')
      end if
      !
      !  TERMINATE.
      !
      if (ifail == 0) then
         Ipass = 1
         if (Kprint >= 2) write (Lun, 99008)
99008    format(/' ------------ DPCHIP PASSED  ALL INTEGRATION TESTS', &
                     &' ------------')
      else
         Ipass = 0
         if (Kprint >= 1) write (Lun, 99009)
99009    format(/' ************ DPCHIP FAILED SOME INTEGRATION TESTS', &
                     &' ************')
      end if
      !
      return
99010 format(2f6.1, i5, 1p, 2d20.10, d15.5)

   end

   subroutine dpchq3(Lun, Kprint, Ipass)
      implicit none
      !***PURPOSE  Test the PCHIP interpolators DPCHIC, DPCHIM, DPCHSP.
      !***AUTHOR  Fritsch, F. N., (LLNL)
      !***DESCRIPTION
      !
      !             DPCHIP QUICK CHECK NUMBER 3
      !
      !     TESTS THE INTERPOLATORS:  DPCHIC, DPCHIM, DPCHSP.
      ! *Usage:
      !
      !        INTEGER  LUN, KPRINT, IPASS
      !
      !        CALL DPCHQ3 (LUN, KPRINT, IPASS)
      !
      ! *Arguments:
      !
      !     LUN   :IN  is the unit number to which output is to be written.
      !
      !     KPRINT:IN  controls the amount of output, as specified in the
      !                SLATEC Guidelines.
      !
      !     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
      !                IPASS=0 indicates one or more tests failed.
      !
      ! *Description:
      !
      !   This routine interpolates a constructed data set with all three
      !   DPCHIP interpolators and compares the results with those obtained
      !   on a Cray X/MP.  Two different values of the DPCHIC parameter SWITCH
      !   are used.
      !
      ! *Remarks:
      !     1. The Cray results are given only to nine significant figures,
      !        so don't expect them to match to more.
      !     2. The results will depend to some extent on the accuracy of
      !        the EXP function.
      !
      !***REVISION HISTORY  (YYMMDD)
      !   900309  DATE WRITTEN
      !   900314  Converted to a subroutine and added a SLATEC 4.0 prologue.
      !   900315  Revised prologue and improved some output formats.  (FNF)
      !   900316  Made TOLD machine-dependent and added extra output when
      !           KPRINT=3.  (FNF)
      !   900320  Added E0's to DATA statement for X to reduce single/double
      !           differences, and other minor cosmetic changes.
      !   900320  Converted to double precision.
      !   900321  Removed IFAIL from call sequence for SLATEC standards and
      !           made miscellaneous cosmetic changes.  (FNF)
      !   900322  Minor changes to reduce single/double differences.  (FNF)
      !   900530  Tolerance (TOLD) and argument to DPCHIC changed.  (WRB)
      !   900802  Modified TOLD formula and constants in DPCHIC calls to
      !           correct DPCHQ3 failures.  (FNF)
      !   901130  Several significant changes:  (FNF)
      !           1. Changed comparison between DPCHIM and DPCHIC to only
      !              require agreement to machine precision.
      !           2. Revised to print more output when KPRINT=3.
      !           3. Added 1P's to formats.
      !   910708  Minor modifications in use of KPRINT.  (WRB)
      !   930317  Improved output formats.  (FNF)

      !
      !*Internal Notes:
      !
      !     TOLD is used to compare with stored Cray results.  Its value
      !          should be consistent with significance of stored values.
      !     TOLZ is used for cases in which exact equality is expected.
      !     TOL  is used for cases in which agreement to machine precision
      !          is expected.
      !**End
      !
      !  Declare arguments.
      !
      integer Lun, Kprint, Ipass
      !
      !  Declare variables.
      !
      integer i, ic(2), ierr, ifail, n, nbad, nbadz, nwk
      parameter(n=9, nwk=2*n)
      double precision d(n), dc(n), dc5, dc6, dm(n), ds(n), err, &
                     & f(n), mone, tol, told, tolz, vc(2), x(n), &
                     & wk(nwk), zero
      parameter(zero=0.0d0, mone=-1.0d0)
      character*6 result
      !
      !  Initialize.
      !
      !       Data.
      data ic/0, 0/
      data x/-2.2d0, -1.2d0, -1.0d0, -0.5d0, -0.01d0, 0.5d0,    &
         & 1.0d0, 2.0d0, 2.2d0/
      !
      !       Results generated on Cray X/MP (9 sign. figs.)
      data dm/0., 3.80027352d-01, 7.17253009d-01, 5.82014161d-01,   &
         & 0., -5.68208031d-01, -5.13501618d-01, -7.77910977d-02,   &
         & -2.45611117d-03/
      data dc5, dc6/1.76950158d-02, -5.69579814d-01/
      data ds/-5.16830792d-02, 5.71455855d-01, 7.40530225d-01,     &
         & 7.63864934d-01, 1.92614386d-02, -7.65324380d-01,          &
         & -7.28209035d-01, -7.98445427d-02, -2.85983446d-02/
      !

      ifail = 0
      !
      !        Set tolerances.
      tol = 10*d1mach4
      told = max(1.0d-7, 10*tol)
      tolz = zero
      !
      if (Kprint >= 3) write (Lun, 99001)
      !
      !  FORMATS.
      !
99001 format('1'//10x, 'TEST DPCHIP INTERPOLATORS')
      if (Kprint >= 2) write (Lun, 99002)
99002 format(//10x, 'DPCHQ3 RESULTS'/10x, '--------------')
      !
      !  Set up data.
      !
      do i = 1, n
         f(i) = exp(-x(i)**2)
      end do
      !
      if (Kprint >= 3) then
         write (Lun, 99003)
99003    format(//5x, 'DATA:'/39x,                                      &
                     &'---------- EXPECTED D-VALUES ----------'/12x, 'X', 9x,  &
                     &'F', 18x, 'DM', 13x, 'DC', 13x, 'DS')
         do i = 1, 4
            write (Lun, 99009) x(i), f(i), dm(i), ds(i)
         end do
         write (Lun, 99010) x(5), f(5), dm(5), dc5, ds(5)
         write (Lun, 99010) x(6), f(6), dm(6), dc6, ds(6)
         do i = 7, n
            write (Lun, 99009) x(i), f(i), dm(i), ds(i)
         end do
      end if
      !
      !  Test DPCHIM.
      !
      if (Kprint >= 3) write (Lun, 99011) 'IM'
      !     --------------------------------
      call dpchim(n, x, f, d, 1, ierr)
      !     --------------------------------
      !        Expect IERR=1 (one monotonicity switch).
      if (Kprint >= 3) write (Lun, 99012) 1
      if (.not. comp(ierr, 1, Lun, Kprint)) then
         ifail = ifail + 1
      else
         if (Kprint >= 3) write (Lun, 99013)
         nbad = 0
         nbadz = 0
         do i = 1, n
            result = '  OK'
            !             D-values should agree with stored values.
            !               (Zero values should agree exactly.)
            if (dm(i) == zero) then
               err = abs(d(i))
               if (err > tolz) then
                  nbadz = nbadz + 1
                  result = '**BADZ'
               end if
            else
               err = abs((d(i) - dm(i))/dm(i))
               if (err > told) then
                  nbad = nbad + 1
                  result = '**BAD'
               end if
            end if
            if (Kprint >= 3) write (Lun, 99014) i, x(i), d(i), err,  &
                                  & result
         end do
         if ((nbadz /= 0) .or. (nbad /= 0)) then
            ifail = ifail + 1
            if ((nbadz /= 0) .and. (Kprint >= 2)) write (Lun, 99004) nbad
99004       format(/'    **', i5,                                       &
                           &' DPCHIM RESULTS FAILED TO BE EXACTLY ZERO.')
            if ((nbad /= 0) .and. (Kprint >= 2)) write (Lun, 99015) nbad, &
                &'IM', told
         else
            if (Kprint >= 2) write (Lun, 99016) 'IM'
         end if
      end if
      !
      !  Test DPCHIC -- options set to reproduce DPCHIM.
      !
      if (Kprint >= 3) write (Lun, 99011) 'IC'
      !     --------------------------------------------------------
      call dpchic(ic, vc, zero, n, x, f, dc, 1, wk, nwk, ierr)
      !     --------------------------------------------------------
      !        Expect IERR=0 .
      if (Kprint >= 3) write (Lun, 99012) 0
      if (.not. comp(ierr, 0, Lun, Kprint)) then
         ifail = ifail + 1
      else
         if (Kprint >= 3) write (Lun, 99013)
         nbad = 0
         do i = 1, n
            result = '  OK'
            !           D-values should agree exactly with those computed by DPCHIM.
            !            (To be generous, will only test to machine precision.)
            err = abs(d(i) - dc(i))
            if (err > tol) then
               nbad = nbad + 1
               result = '**BAD'
            end if
            if (Kprint >= 3) write (Lun, 99014) i, x(i), dc(i), err, &
                                  & result
         end do
         if (nbad /= 0) then
            ifail = ifail + 1
            if (Kprint >= 2) write (Lun, 99015) nbad, 'IC', tol
         else
            if (Kprint >= 2) write (Lun, 99016) 'IC'
         end if
      end if
      !
      !  Test DPCHIC -- default nonzero switch derivatives.
      !
      if (Kprint >= 3) write (Lun, 99011) 'IC'
      !     -------------------------------------------------------
      call dpchic(ic, vc, mone, n, x, f, d, 1, wk, nwk, ierr)
      !     -------------------------------------------------------
      !        Expect IERR=0 .
      if (Kprint >= 3) write (Lun, 99012) 0
      if (.not. comp(ierr, 0, Lun, Kprint)) then
         ifail = ifail + 1
      else
         if (Kprint >= 3) write (Lun, 99013)
         nbad = 0
         nbadz = 0
         do i = 1, n
            result = '  OK'
            !            D-values should agree exactly with those computed in
            !            previous call, except at points 5 and 6.
            if ((i < 5) .or. (i > 6)) then
               err = abs(d(i) - dc(i))
               if (err > tolz) then
                  nbadz = nbadz + 1
                  result = '**BADA'
               end if
            else
               if (i == 5) then
                  err = abs((d(i) - dc5)/dc5)
               else
                  err = abs((d(i) - dc6)/dc6)
               end if
               if (err > told) then
                  nbad = nbad + 1
                  result = '**BAD'
               end if
            end if
            if (Kprint >= 3) write (Lun, 99014) i, x(i), d(i), err,  &
                                  & result
         end do
         if ((nbadz /= 0) .or. (nbad /= 0)) then
            ifail = ifail + 1
            if ((nbadz /= 0) .and. (Kprint >= 2)) write (Lun, 99005) nbad
99005       format(/'    **', i5, ' DPCHIC RESULTS FAILED TO AGREE WITH',&
                           &' PREVIOUS CALL.')
            if ((nbad /= 0) .and. (Kprint >= 2)) write (Lun, 99015) nbad, &
                &'IC', told
         else
            if (Kprint >= 2) write (Lun, 99016) 'IC'
         end if
      end if
      !
      !  Test DPCHSP.
      !
      if (Kprint >= 3) write (Lun, 99011) 'SP'
      !     -------------------------------------------------
      call dpchsp(ic, vc, n, x, f, d, 1, wk, nwk, ierr)
      !     -------------------------------------------------
      !        Expect IERR=0 .
      if (Kprint >= 3) write (Lun, 99012) 0
      if (.not. comp(ierr, 0, Lun, Kprint)) then
         ifail = ifail + 1
      else
         if (Kprint >= 3) write (Lun, 99013)
         nbad = 0
         do i = 1, n
            result = '  OK'
            !             D-values should agree with stored values.
            err = abs((d(i) - ds(i))/ds(i))
            if (err > told) then
               nbad = nbad + 1
               result = '**BAD'
            end if
            if (Kprint >= 3) write (Lun, 99014) i, x(i), d(i), err,  &
                                  & result
         end do
         if (nbad /= 0) then
            ifail = ifail + 1
            if (Kprint >= 2) write (Lun, 99015) nbad, 'SP', told
         else
            if (Kprint >= 2) write (Lun, 99016) 'SP'
         end if
      end if
      !
      !  PRINT SUMMARY AND TERMINATE.
      !
      if ((Kprint >= 2) .and. (ifail /= 0)) write (Lun, 99006) ifail
99006 format(/' *** TROUBLE ***', i5, ' INTERPOLATION TESTS FAILED.')
      !
      if (ifail == 0) then
         Ipass = 1
         if (Kprint >= 2) write (Lun, 99007)
99007    format(/' ------------ DPCHIP PASSED  ALL INTERPOLATION TESTS'&
                    & , ' ------------')
      else
         Ipass = 0
         if (Kprint >= 1) write (Lun, 99008)
99008    format(/' ************ DPCHIP FAILED SOME INTERPOLATION TESTS'&
                    & , ' ************')
      end if
      !
      return
99009 format(5x, f10.2, 1p, d15.5, 4x, d15.5, 15x, d15.5)
99010 format(5x, f10.2, 1p, d15.5, 4x, 3d15.5)
99011 format(/5x, 'DPCH', a2, ' TEST:')
99012 format(15x, 'EXPECT  IERR =', i5)
99013 format(/9x, 'I', 7x, 'X', 9x, 'D', 13x, 'ERR')
99014 format(5x, i5, f10.2, 1p, 2d15.5, 2x, a)
99015 format(/'    **', i5, ' DPCH', a2, ' RESULTS FAILED TOLERANCE TEST.',&
               &'  TOL =', 1p, d10.3)
99016 format(/5x, '  ALL DPCH', a2, ' RESULTS OK.')

   end

   subroutine dpchq4(Lun, Kprint, Ipass)
      implicit none
      !***PURPOSE  Test the PCHIP monotonicity checker DPCHCM.
      !***AUTHOR  Fritsch, F. N., (LLNL)
      !***DESCRIPTION
      !
      !             DPCHIP QUICK CHECK NUMBER 4
      !
      !     TESTS THE MONOTONICITY CHECKER:  DPCHCM.
      ! *Usage:
      !
      !        INTEGER  LUN, KPRINT, IPASS
      !
      !        CALL DPCHQ4 (LUN, KPRINT, IPASS)
      !
      ! *Arguments:
      !
      !     LUN   :IN  is the unit number to which output is to be written.
      !
      !     KPRINT:IN  controls the amount of output, as specified in the
      !                SLATEC Guidelines.
      !
      !     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
      !                IPASS=0 indicates one or more tests failed.
      !
      ! *Description:
      !
      !   This routine tests a constructed data set with three different
      !   INCFD settings and compares with the expected results.  It then
      !   runs a special test to check for bug in overall monotonicity found
      !   in DPCHMC.  Finally, it reverses the data and repeats all tests.
      !
      !***REVISION HISTORY  (YYMMDD)
      !   890208  DATE WRITTEN
      !   890306  Changed LOUT to LUN and added it to call list.  (FNF)
      !   890316  Removed DATA statements to suit new quick check standards.
      !   890410  Changed PCHMC to PCHCM.
      !   890410  Added a SLATEC 4.0 format prologue.
      !   900314  Changed name from PCHQK3 to PCHQK4 and improved some output
      !           formats.
      !   900315  Revised prologue and improved some output formats.  (FNF)
      !   900320  Converted to double precision.
      !   900321  Removed IFAIL from call sequence for SLATEC standards and
      !           made miscellaneous cosmetic changes.  (FNF)
      !   900322  Added declarations so all variables are declared.  (FNF)
      !   910708  Minor modifications in use of KPRINT.  (WRB)
      !   930317  Improved output formats.  (FNF)

      !
      !*Internal Notes:
      !
      !     Data set-up is done via assignment statements to avoid modifying
      !     DATA-loaded arrays, as required by the 1989 SLATEC Guidelines.
      !     Run with KPRINT=3 to display the data.
      !**End
      !
      !  Declare arguments.
      !
      integer Lun, Kprint, Ipass
      !
      !  DECLARE VARIABLES.
      !
      integer maxn, maxn2, maxn3, nb
      parameter(maxn=16, maxn2=8, maxn3=6, nb=7)
      integer i, ierr, ifail, incfd, ismex1(maxn), ismex2(maxn2), &
            & ismex3(maxn3), ismexb(nb), ismon(maxn), k, n, ns(3)
      double precision d(maxn), db(nb), f(maxn), fb(nb), x(maxn)
      logical skip
      !
      !  DEFINE EXPECTED RESULTS.
      !
      data ismex1/1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, -1,   &
         & 1, 1, -1, 2/
      data ismex2/1, 2, 2, 1, 2, 2, 1, 2/
      data ismex3/1, 1, 1, 1, 1, 1/
      data ismexb/1, 3, 1, -1, -3, -1, 2/
      !
      !  DEFINE TEST DATA.
      !
      data ns/16, 8, 6/
      !

      if (Kprint >= 3) write (Lun, 99001)
      !
      !  FORMATS.
      !
99001 format('1'//10x, 'TEST DPCHIP MONOTONICITY CHECKER')
      if (Kprint >= 2) write (Lun, 99002)
99002 format(//10x, 'DPCHQ4 RESULTS'/10x, '--------------')
      !
      !       Define X, F, D.
      do i = 1, maxn
         x(i) = i
         d(i) = 0.d0
      end do
      do i = 2, maxn, 3
         d(i) = 2.d0
      end do
      do i = 1, 3
         f(i) = x(i)
         f(i + 3) = f(i) + 1.d0
         f(i + 6) = f(i + 3) + 1.d0
         f(i + 9) = f(i + 6) + 1.d0
         f(i + 12) = f(i + 9) + 1.d0
      end do
      f(16) = 6.d0
      !       Define FB, DB.
      fb(1) = 0.d0
      fb(2) = 2.d0
      fb(3) = 3.d0
      fb(4) = 5.d0
      db(1) = 1.d0
      db(2) = 3.d0
      db(3) = 3.d0
      db(4) = 0.d0
      do i = 1, 3
         fb(nb - i + 1) = fb(i)
         db(nb - i + 1) = -db(i)
      end do
      !
      !  INITIALIZE.
      !
      ifail = 0
      !
      if (Kprint >= 3) then
         write (Lun, 99003)
99003    format(//5x, 'DATA:'//9x, 'I', 4x, 'X', 5x, 'F', 5x, 'D', 5x, 'FB', 4x,  &
                     &'DB')
         do i = 1, nb
            write (Lun, 99010) i, x(i), f(i), d(i), fb(i), db(i)
         end do
         do i = nb + 1, maxn
            write (Lun, 99010) i, x(i), f(i), d(i)
         end do
      end if
      !
      !  TRANSFER POINT FOR SECOND SET OF TESTS.
      !
      !
      !  Loop over a series of values of INCFD.
      !
100   do incfd = 1, 3
         n = ns(incfd)
         skip = .false.
         !        -------------------------------------------------
         call dpchcm(n, x, f, d, incfd, skip, ismon, ierr)
         !        -------------------------------------------------
         if (Kprint >= 3) write (Lun, 99004) incfd, ierr,              &
                               & (ismon(i), i=1, n)
99004    format(/4x, 'INCFD =', i2, ':  IERR =', i3/15x, 'ISMON =', 16i3)
         if (ierr /= 0) then
            ifail = ifail + 1
            if (Kprint >= 3) write (Lun, 99011)
         else
            do i = 1, n
               if (incfd == 1) then
                  if (ismon(i) /= ismex1(i)) then
                     ifail = ifail + 1
                     if (Kprint >= 3) write (Lun, 99012)                 &
                        & (ismex1(k), k=1, n)
                     goto 200
                  end if
               elseif (incfd == 2) then
                  if (ismon(i) /= ismex2(i)) then
                     ifail = ifail + 1
                     if (Kprint >= 3) write (Lun, 99012)                 &
                        & (ismex2(k), k=1, n)
                     goto 200
                  end if
               elseif (ismon(i) /= ismex3(i)) then
                  ifail = ifail + 1
                  if (Kprint >= 3) write (Lun, 99012) (ismex3(k), k=1, n)
                  goto 200
               end if
            end do
         end if
200   end do
      !
      !  Test for -1,3,1 bug.
      !
      skip = .false.
      !     ------------------------------------------------
      call dpchcm(nb, x, fb, db, 1, skip, ismon, ierr)
      !     ------------------------------------------------
      if (Kprint >= 3) write (Lun, 99005) ierr, (ismon(i), i=1, nb)
99005 format(/4x, ' Bug test:  IERR =', i3/15x, 'ISMON =', 7i3)
      if (ierr /= 0) then
         ifail = ifail + 1
         if (Kprint >= 3) write (Lun, 99011)
      else
         do i = 1, nb
            if (ismon(i) /= ismexb(i)) then
               ifail = ifail + 1
               if (Kprint >= 3) write (Lun, 99012) (ismexb(k), k=1, nb)
               goto 300
            end if
         end do
      end if
      !
300   if (f(1) < 0.) then
         !
         !  PRINT SUMMARY AND TERMINATE.
         !
         if ((Kprint >= 2) .and. (ifail /= 0)) write (Lun, 99006) ifail
99006    format(/' *** TROUBLE ***', i5, ' MONOTONICITY TESTS FAILED.')
         !
         if (ifail == 0) then
            Ipass = 1
            if (Kprint >= 2) write (Lun, 99007)
99007       format(/                                                   &
                          &' ------------ DPCHIP PASSED  ALL MONOTONICITY TESTS'&
                         & , ' ------------')
         else
            Ipass = 0
            if (Kprint >= 1) write (Lun, 99008)
99008       format(/                                                   &
                          &' ************ DPCHIP FAILED SOME MONOTONICITY TESTS'&
                         & , ' ************')
         end if
         !
         return
      else
         !
         !  Change sign and do again.
         !
         if (Kprint >= 3) write (Lun, 99009)
99009    format(/4x, 'Changing sign of data.....')
         do i = 1, maxn
            f(i) = -f(i)
            d(i) = -d(i)
            if (ismex1(i) /= 2) ismex1(i) = -ismex1(i)
         end do
         do i = 1, maxn2
            if (ismex2(i) /= 2) ismex2(i) = -ismex2(i)
         end do
         do i = 1, maxn3
            if (ismex3(i) /= 2) ismex3(i) = -ismex3(i)
         end do
         do i = 1, nb
            fb(i) = -fb(i)
            db(i) = -db(i)
            if (ismexb(i) /= 2) ismexb(i) = -ismexb(i)
         end do
         goto 100
      end if
99010 format(5x, i5, 5f6.1)
99011 format(' *** Failed -- bad IERR value.')
99012 format(' *** Failed -- expect:', 16i3)

   end

   SUBROUTINE DPCHQ5(Lun, Kprint, Ipass)
      IMPLICIT NONE
!***PURPOSE  Test the PCH to B-spline conversion routine DPCHBS.
!***AUTHOR  Fritsch, F. N., (LLNL)
!***DESCRIPTION
!
!             DPCHIP QUICK CHECK NUMBER 5
!
!     TESTS THE CONVERSION ROUTINE:  DPCHBS.
! *Usage:
!
!        INTEGER  LUN, KPRINT, IPASS
!
!        CALL DPCHQ5 (LUN, KPRINT, IPASS)
!
! *Arguments:
!
!     LUN   :IN  is the unit number to which output is to be written.
!
!     KPRINT:IN  controls the amount of output, as specified in the
!                SLATEC Guidelines.
!
!     IPASS:OUT  will contain a pass/fail flag.  IPASS=1 is good.
!                IPASS=0 indicates one or more tests failed.
!
! *Description:
!
!   This routine tests a constructed data set with four different
!   KNOTYP settings.  It computes the function and derivatives of the
!   resulting B-representation via DBVALU and compares with PCH data.
!
! *Caution:
!   This routine assumes DBVALU has already been successfully tested.
!
!***ROUTINES CALLED  DBVALU, DPCHBS, D1MACH
!***REVISION HISTORY  (YYMMDD)
!   900411  DATE WRITTEN
!   900412  Corrected minor errors in initial implementation.
!   900430  Produced double precision version.
!   900501  Corrected declarations.
!   930317  Improved output formats.  (FNF)
!***END PROLOGUE  DPCHQ5
!
!*Internal Notes:
!  TOL  is the tolerance to use for quantities that should only
!       theoretically be equal.
!  TOLZ is the tolerance to use for quantities that should be exactly
!       equal.
!
!**End
!
!  Declare arguments.
!
      INTEGER Lun, Kprint, Ipass
!
!  Declare externals.
!
      !  DOUBLE PRECISION DBVALU , D1MACH
      !  EXTERNAL DBVALU , DPCHBS , D1MACH
!
!  Declare variables.
!
      INTEGER i, ierr, ifail, inbv, j, knotyp, k, N, ndim,     &
            & nknots
      PARAMETER(N=9)
      DOUBLE PRECISION bcoef(2*N), d(N), dcalc, derr, dermax,      &
                     & f(N), fcalc, ferr, fermax, t(2*N + 4), terr, &
                     & termax, tol, tolz, tsave(2*N + 4), work(16*N),&
                     & x(N), ZERO
      PARAMETER(ZERO=0.0D0)
      LOGICAL fail
!
!  Define relative error function.
!
      DOUBLE PRECISION ans, err, RELERR
      RELERR(err, ans) = ABS(err)/MAX(1.0D-5, ABS(ans))
!
!  Define test data.
!
      DATA x/-2.2D0, -1.2D0, -1.0D0, -0.5D0, -0.01D0, 0.5D0,    &
         & 1.0D0, 2.0D0, 2.2D0/
      DATA f/0.0079D0, 0.2369D0, 0.3679D0, 0.7788D0, 0.9999D0,     &
         & 0.7788D0, 0.3679D0, 0.1083D0, 0.0079D0/
      DATA d/0.0000D0, 0.3800D0, 0.7173D0, 0.5820D0, 0.0177D0,     &
         & -0.5696D0, -0.5135D0, -0.0778D0, -0.0025D0/
!
!  Initialize.
!
!***FIRST EXECUTABLE STATEMENT  DPCHQ5
      ifail = 0
      tol = 100*D1MACH4
      tolz = ZERO
!
      IF (Kprint >= 3) WRITE (Lun, 99001)
!
!  FORMATS.
!
99001 FORMAT('1'//10X, 'TEST PCH TO B-SPLINE CONVERTER')
      IF (Kprint >= 2) WRITE (Lun, 99002)
99002 FORMAT(//10X, 'DPCHQ5 RESULTS'/10X, '--------------')
!
!  Loop over a series of values of KNOTYP.
!
      IF (Kprint >= 3) WRITE (Lun, 99003)
99003 FORMAT(/4X, '(Results should be the same for all KNOTYP values.)')
      DO knotyp = 2, -1, -1
!        ------------
         CALL DPCHBS(N, x, f, d, 1, knotyp, nknots, t, bcoef, ndim, k, ierr)
!        ------------
         IF (Kprint >= 3) WRITE (Lun, 99004) knotyp, nknots, ndim, k,&
                               & ierr
99004    FORMAT(/4X, 'KNOTYP =', I2, ':  NKNOTS =', I3, ',  NDIM =', I3,     &
                 &',  K =', I2, ',  IERR =', I3)
         IF (ierr /= 0) THEN
            ifail = ifail + 1
            IF (Kprint >= 3) WRITE (Lun, 99005)
99005       FORMAT(' *** Failed -- bad IERR value.')
         ELSE
!             Compare evaluated results with inputs to DPCHBS.
            inbv = 1
            fermax = ZERO
            dermax = ZERO
            IF (Kprint >= 3) THEN
               WRITE (Lun, 99006)
99006          FORMAT(/15X, 'X', 9X, 'KNOTS', 10X, 'F', 7X, 'FERR', 8X, 'D', 7X, &
                             &'DERR')
               WRITE (Lun, 99013) t(1), t(2)
               j = 1
            END IF
            DO i = 1, N
               fcalc = DBVALU(t, bcoef, ndim, k, 0, x(i), inbv, work)
               ferr = f(i) - fcalc
               fermax = MAX(fermax, RELERR(ferr, f(i)))
               dcalc = DBVALU(t, bcoef, ndim, k, 1, x(i), inbv, work)
               derr = d(i) - dcalc
               dermax = MAX(dermax, RELERR(derr, d(i)))
               IF (Kprint >= 3) THEN
                  j = j + 2
                  WRITE (Lun, 99007) x(i), t(j), t(j + 1), f(i), ferr,&
                                  & d(i), derr
99007             FORMAT(10X, 3F8.2, F10.4, 1P, D10.2, 0P, F10.4, 1P, D10.2)
               END IF
            END DO
            IF (Kprint >= 3) THEN
               j = j + 2
               WRITE (Lun, 99013) t(j), t(j + 1)
            END IF
            fail = (fermax > tol) .OR. (dermax > tol)
            IF (fail) ifail = ifail + 1
            IF ((Kprint >= 3) .OR. (Kprint >= 2) .AND. fail)              &
               & WRITE (Lun, 99008) fermax, dermax, tol
99008       FORMAT(/5X, 'Maximum relative errors:'/15X, 'F-error =', 1P,  &
                      & D13.5, 5X, 'D-error =', D13.5/5X,                      &
                       &'Both should be less than  TOL =', D13.5)
         END IF
!
!          Special check for KNOTYP=-1.
         IF (knotyp == 0) THEN
!             Save knot vector for next test.
            DO i = 1, nknots
               tsave(i) = t(i)
            END DO
         ELSEIF (knotyp == -1) THEN
!             Check that knot vector is unchanged.
            termax = ZERO
            DO i = 1, nknots
               terr = ABS(t(i) - tsave(i))
               termax = MAX(termax, terr)
            END DO
            IF (termax > tolz) THEN
               ifail = ifail + 1
               IF (Kprint >= 2) WRITE (Lun, 99009) termax, tolz
99009          FORMAT(/' *** T-ARRAY MAXIMUM CHANGE =', 1P, D13.5,       &
                             &';  SHOULD NOT EXCEED TOLZ =', D13.5)
            END IF
         END IF
      END DO
!
!  PRINT SUMMARY AND TERMINATE.
!
      IF ((Kprint >= 2) .AND. (ifail /= 0)) WRITE (Lun, 99010) ifail
99010 FORMAT(/' *** TROUBLE ***', I5, ' CONVERSION TESTS FAILED.')
!
      IF (ifail == 0) THEN
         Ipass = 1
         IF (Kprint >= 2) WRITE (Lun, 99011)
99011    FORMAT(/' ------------ DPCHIP PASSED  ALL CONVERSION TESTS',  &
                 &' ------------')
      ELSE
         Ipass = 0
         IF (Kprint >= 1) WRITE (Lun, 99012)
99012    FORMAT(/' ************ DPCHIP FAILED SOME CONVERSION TESTS',  &
                 &' ************')
      END IF
!
      RETURN
99013 FORMAT(18X, 2F8.2)
!------------- LAST LINE OF DPCHQ5 FOLLOWS -----------------------------
   END

   subroutine devchk(Lout, Kprint, Npts, Xev, Fev, Dev, Fev2, Fail)
      implicit none
      !***PURPOSE  Test evaluation accuracy of DCHFDV and DCHFEV for DPCHQ1.
      !***AUTHOR  Fritsch, F. N., (LLNL)
      !***DESCRIPTION
      !
      ! -------- CODE TO TEST EVALUATION ACCURACY OF DCHFDV AND DCHFEV -------
      !
      !     USING FUNCTION AND DERIVATIVE VALUES FROM A CUBIC (COMPUTED IN
      !     DOUBLE PRECISION) AT NINT DIFFERENT (X1,X2) PAIRS:
      !     1. CHECKS THAT DCHFDV AND DCHFEV BOTH REPRODUCE ENDPOINT VALUES.
      !     2. EVALUATES AT NPTS POINTS, 10 OF WHICH ARE OUTSIDE THE INTERVAL
      !        AND:
      !        A. CHECKS ACCURACY OF DCHFDV FUNCTION AND DERIVATIVE VALUES
      !           AGAINST EXACT VALUES.
      !        B. CHECKS THAT RETURNED VALUES OF NEXT SUM TO 10.
      !        C. CHECKS THAT FUNCTION VALUES FROM DCHFEV AGREE WITH THOSE
      !           FROM DCHFDV.
      !
      !
      !     FORTRAN INTRINSICS USED:  ABS, MAX, MIN.
      !     FORTRAN LIBRARY ROUTINES USED:  SQRT, (READ), (WRITE).
      !     SLATEC LIBRARY ROUTINES USED:  DCHFDV, DCHFEV, D1MACH, RAND.
      !     OTHER ROUTINES USED:  DFDTRU.
      !
      !***REVISION HISTORY  (YYMMDD)
      !   820601  DATE WRITTEN
      !   820624  CONVERTED TO QUICK CHECK FOR SLATEC LIBRARY.
      !   820630  1. MODIFIED DEFINITIONS OF RELATIVE ERROR AND TEST
      !             TOLERANCES.
      !           2. VARIOUS IMPROVEMENTS TO OUTPUT FORMATS.
      !   820716  1. SET MACHEP VIA A CALL TO D1MACH.
      !           2. CHANGED FROM FORTLIB'S RANF TO SLATEC'S RAND.
      !   890628  1. Removed unnecessary IMPLICIT declaration.
      !           2. Removed unnecessary variable NEV.
      !           3. Other changes to reduce S.P./D.P. differences.
      !   890629  Added RERR to DOUBLE PRECISION declaration.
      !   890706  Cosmetic changes to prologue.  (WRB)
      !   890831  Modified array declarations.  (WRB)
      !   890911  Removed unnecessary intrinsics.  (WRB)
      !   891214  Prologue converted to Version 4.0 format.  (BAB)
      !   900315  Revised prologue and improved some output formats.  (FNF)
      !           Also moved formats to end to be consistent with other PCHIP
      !           quick checks.
      !   900316  Additional minor cosmetic changes.  (FNF)
      !   900321  Changed name of DFTRUE to DFDTRU and made additional minor
      !           cosmetic changes.  (FNF)
      !   901130  Added 1P's to formats and revised some to reduce maximum
      !           line length.  (FNF)
      !   910708  Minor modifications in use of KPRINT.  (WRB)
      !   910801  Added EXTERNAL statement for RAND due to problem on IBM
      !           RS 6000.  (WRB)
      !   910819  Changed argument to RAND function from a D.P. zero to a
      !           S.P. zero.  (WRB)

      !
      !  Declare arguments.
      !
      integer Lout, Kprint, Npts
      double precision Xev(*), Fev(*), Dev(*), Fev2(*)
      logical Fail
      !
      !  DECLARATIONS.
      !
      integer i, ierr, iint, next(2), next2(2), nint
      double precision aed, aed2, aedmax, aedmin, aef, aef2,      &
                     & aefmax, aefmin, check(2), checkf(2),         &
                     & checkd(2), d1, d2, dermax, dtrue, dx,      &
                     & eps1, eps2, f1, f2, fact, fermax, floord, &
                     & floorf, four, ftrue, left(3), machep, one, &
                     & red, red2, redmax, redmin, ref, ref2,      &
                     & refmax, refmin, right(3), small, ten, tol1,&
                     & tol2, x1, x2, xadmax, xadmin, xafmax,      &
                     & xafmin, xrdmax, xrdmin, xrfmax, xrfmin, zero
      logical failoc, failnx
      !
      !       The following should stay REAL (no D.P. equivalent).
      !real rand
      !external rand
      !
      !  DEFINE RELATIVE ERROR WITH FLOOR.
      !
      double precision rerr, err, value, floor
      rerr(err, value, floor) = err/max(abs(value), floor)
      !
      !  INITIALIZE.
      !
      data zero/0.d0/, one/1.d0/, four/4.d0/, ten/10.d0/
      data small/1.0d-10/
      data nint/3/
      data left/-1.5d0, 2.0d-10, 1.0d0/
      data right/2.5d0, 3.0d-10, 1.0d+8/
      !

      machep = d1mach4
      eps1 = four*machep
      eps2 = ten*machep
      !
      Fail = .false.
      !
      if (Kprint >= 2) write (Lout, 99001)
99001 format(//10x, 'DEVCHK RESULTS'/10x, '--------------')
      !
      !  CYCLE OVER INTERVALS.
      !
      do iint = 1, nint
         x1 = left(iint)
         x2 = right(iint)
         !
         fact = max(sqrt(x2 - x1), one)
         tol1 = eps1*fact
         tol2 = eps2*fact
         !
         !  COMPUTE AND PRINT ENDPOINT VALUES.
         !
         call dfdtru(x1, f1, d1)
         call dfdtru(x2, f2, d2)
         !
         if (Kprint >= 3) then
            if (iint == 1) write (Lout, 99002)
            !
            !  FORMATS.
            !
99002       format(/10x, 'DCHFDV ACCURACY TEST')
            write (Lout, '(/)')
            write (Lout, 99017) 'X1', x1, 'X2', x2
            write (Lout, 99017) 'F1', f1, 'F2', f2
            write (Lout, 99017) 'D1', d1, 'D2', d2
         end if
         !
         if (Kprint >= 2) write (Lout, 99003) x1, x2
99003    format(/10x, 'INTERVAL = (', 1p, d12.5, ',', d12.5, ' ):')
         !
         !  COMPUTE FLOORS FOR RELATIVE ERRORS.
         !
         floorf = max(min(abs(f1), abs(f2)), small)
         floord = max(min(abs(d1), abs(d2)), small)
         !
         !  CHECK REPRODUCTION OF ENDPOINT VALUES.
         !
         Xev(1) = x1
         Xev(2) = x2
         !     -----------------------------------------------------------
         call dchfdv(x1, x2, f1, f2, d1, d2, 2, Xev, checkf, checkd, next, ierr)
         !     -----------------------------------------------------------
         aef = checkf(1) - f1
         ref = rerr(aef, f1, floorf)
         aef2 = checkf(2) - f2
         ref2 = rerr(aef2, f2, floorf)
         aed = checkd(1) - d1
         red = rerr(aed, d1, floord)
         aed2 = checkd(2) - d2
         red2 = rerr(aed2, d2, floord)
         !
         failoc = max(abs(ref), abs(ref2), abs(red), abs(red2)) > tol1
         Fail = Fail .or. failoc
         !
         if (Kprint >= 3) then
            write (Lout, 99004) next, aef, aef2, aed, aed2
99004       format(/' ERRORS AT ENDPOINTS:', 40x, '(NEXT =', 2i3, ')'//1p, &
                           & 4x, 'F1:', d13.5, 4x, 'F2:', d13.5, 4x, 'D1:', d13.5, 4x,    &
                            &'D2:', d13.5)
            write (Lout, 99005) ref, ref2, red, red2
99005       format(1p, 4(7x, d13.5))
         end if
         !
         if (failoc .and. (Kprint >= 2)) write (Lout, 99006)
99006    format(/' ***** DCHFDV FAILED TO REPRODUCE ENDPOINT VALUES.')
         !
         !  DCHFEV SHOULD AGREE EXACTLY WITH DCHFDV.
         !                     -------
         !     --------------------------------------------------------------
         call dchfev(x1, x2, f1, f2, d1, d2, 2, Xev, check, next, ierr)
         !     --------------------------------------------------------------
         failoc = (check(1) /= checkf(1)) .or. (check(2) /= checkf(2))
         Fail = Fail .or. failoc
         !
         if (failoc .and. (Kprint >= 2)) write (Lout, 99007)
99007    format(/                                                      &
                     &' ***** DCHFEV DOES NOT AGREE WITH DCHFDV AT ENDPOINTS.'&
                    & )
         !
         !  EVALUATE AT NPTS 'UNIFORMLY RANDOM' POINTS IN (X1,X2).
         !     THIS VERSION EXTENDS EVALUATION DOMAIN BY ADDING 4 SUBINTERVALS
         !     TO LEFT AND 6 TO RIGHT OF [X1,X2].
         !
         dx = (x2 - x1)/(Npts - 10)
         do i = 1, Npts
            Xev(i) = (x1 + (i - 5)*dx) + dx*rand(1) !! JW mod - replace with intrinsic random number generator -TODO
         end do
         !     --------------------------------------------------------
         call dchfdv(x1, x2, f1, f2, d1, d2, Npts, Xev, Fev, Dev, next, ierr)
         !     --------------------------------------------------------
         if (ierr /= 0) then
            failoc = .true.
            if (Kprint >= 2) write (Lout, 99008) ierr
99008       format(/' ***** ERROR ***** DCHFDV RETURNED IERR =', i5)
         else
            !
            !     CUMULATE LARGEST AND SMALLEST ERRORS FOR SUMMARY.
            !
            do i = 1, Npts
               call dfdtru(Xev(i), ftrue, dtrue)
               aef = Fev(i) - ftrue
               ref = rerr(aef, ftrue, floorf)
               aed = Dev(i) - dtrue
               red = rerr(aed, dtrue, floord)
               !
               if (i == 1) then
                  !            INITIALIZE.
                  aefmin = aef
                  aefmax = aef
                  aedmin = aed
                  aedmax = aed
                  refmin = ref
                  refmax = ref
                  redmin = red
                  redmax = red
                  xafmin = Xev(1)
                  xafmax = Xev(1)
                  xadmin = Xev(1)
                  xadmax = Xev(1)
                  xrfmin = Xev(1)
                  xrfmax = Xev(1)
                  xrdmin = Xev(1)
                  xrdmax = Xev(1)
               else
                  !            SELECT.
                  if (aef < aefmin) then
                     aefmin = aef
                     xafmin = Xev(i)
                  elseif (aef > aefmax) then
                     aefmax = aef
                     xafmax = Xev(i)
                  end if
                  if (aed < aedmin) then
                     aedmin = aed
                     xadmin = Xev(i)
                  elseif (aed > aedmax) then
                     aedmax = aed
                     xadmax = Xev(i)
                  end if
                  if (ref < refmin) then
                     refmin = ref
                     xrfmin = Xev(i)
                  elseif (ref > refmax) then
                     refmax = ref
                     xrfmax = Xev(i)
                  end if
                  if (red < redmin) then
                     redmin = red
                     xrdmin = Xev(i)
                  elseif (red > redmax) then
                     redmax = red
                     xrdmax = Xev(i)
                  end if
               end if
            end do
            !
            fermax = max(abs(refmax), abs(refmin))
            dermax = max(abs(redmax), abs(redmin))
            !
            failnx = (next(1) + next(2)) /= 10
            failoc = failnx .or. (max(fermax, dermax) > tol2)
         end if
         Fail = Fail .or. failoc
         !
         !  PRINT SUMMARY.
         !
         if (Kprint >= 3) then
            write (Lout, 99009) Npts - 10, next
99009       format(/' ERRORS AT ', i5, ' INTERIOR POINTS + 10 OUTSIDE:', &
                           & 15x, '(NEXT =', 2i3, ')'//30x, 'FUNCTION', 17x,          &
                            &'DERIVATIVE'/15x, 2(11x, 'ABS', 9x, 'REL'))
            !
            write (Lout, 99018) 'MIN', aefmin, refmin, aedmin, redmin
            write (Lout, 99019) xafmin, xrfmin, xadmin, xrdmin
            write (Lout, 99018) 'MAX', aefmax, refmax, aedmax, redmax
            write (Lout, 99019) xafmax, xrfmax, xadmax, xrdmax
         end if
         !
         if (Kprint >= 2) then
            if (failoc) then
               if (fermax > tol2) write (Lout, 99020) 'F', fermax, tol2
               if (dermax > tol2) write (Lout, 99020) 'D', dermax, tol2
               if (failnx) write (Lout, 99010) next
99010          format(/' ***** REPORTED NEXT =', 2i5,                   &
                                  &'   RATHER THAN    4    6')
            else
               write (Lout, 99011)
99011          format(/' DCHFDV RESULTS OK.')
            end if
         end if
         !
         !  CHECK THAT DCHFEV AGREES WITH DCHFDV.
         !
         !     -----------------------------------------------------------------
         call dchfev(x1, x2, f1, f2, d1, d2, Npts, Xev, Fev2, next2, ierr)
         !     -----------------------------------------------------------------
         if (ierr /= 0) then
            failoc = .true.
            if (Kprint >= 2) write (Lout, 99012) ierr
99012       format(/' ***** ERROR ***** DCHFEV RETURNED IERR =', i5)
         else
            aefmax = abs(Fev2(1) - Fev(1))
            xafmax = Xev(1)
            do i = 2, Npts
               aef = abs(Fev2(i) - Fev(i))
               if (aef > aefmax) then
                  aefmax = aef
                  xafmax = Xev(i)
               end if
            end do
            failnx = (next2(1) /= next(1)) .or. (next2(2) /= next(2))
            failoc = failnx .or. (aefmax /= zero)
            if (Kprint >= 2) then
               if (failoc) then
                  write (Lout, 99013)
99013             format(/' ***** DCHFEV DID NOT AGREE WITH DCHFDV:')
                  if (aefmax /= zero) write (Lout, 99014) aefmax, xafmax
99014             format(7x, 'MAXIMUM DIFFERENCE ', 1p, d12.5,            &
                                        &'; OCCURRED AT X =', d12.5)
                  if (failnx) write (Lout, 99015) next2, next
99015             format(7x, 'REPORTED NEXT =', 2i3, '   RATHER THAN ',   &
                                       & 2i3)
               else
                  write (Lout, 99016)
99016             format(/' DCHFEV AGREES WITH DCHFDV.')
               end if
            end if
         end if
         !
         Fail = Fail .or. failoc
         !
         !  GO BACK FOR ANOTHER INTERVAL.
         !
      end do
      !
      return
99017 format(10x, a2, ' =', 1p, d18.10, 5x, a2, ' =', d18.10)
99018 format(/5x, a3, 'IMUM ERROR:  ', 1p, 2d12.4, 2x, 2d12.4)
99019 format(5x, 'LOCATED AT X =  ', 1p, 2d12.4, 2x, 2d12.4)
99020 format(/' ***** MAXIMUM RELATIVE ERROR IN ', a1, ' =', 1p, d12.5,    &
               & ','/17x, 'EXCEEDS TOLERANCE =', d12.5)

   end

   subroutine deverk(Lout, Kprint, Fail)
      implicit none
      !***PURPOSE  Test error returns from DPCHIP evaluators for DPCHQ1.
      !***AUTHOR  Fritsch, F. N., (LLNL)
      !***DESCRIPTION
      !
      ! --------- CODE TO TEST ERROR RETURNS FROM DPCHIP EVALUATORS. ---------
      !
      !
      !     FORTRAN LIBRARY ROUTINES USED:  (WRITE).
      !     SLATEC LIBRARY ROUTINES USED:  DCHFDV, DCHFEV, DPCHFD, DPCHFE,
      !                                    XERDMP, XGETF, XSETF.
      !     OTHER ROUTINES USED:  COMP.
      !
      !***REVISION HISTORY  (YYMMDD)
      !   820601  DATE WRITTEN
      !   820715  CONVERTED TO QUICK CHECK FOR SLATEC LIBRARY.
      !   890207  ADDED CALLS TO ERROR HANDLER.
      !   890316  Added call to XERDMP if KPRINT.GT.2 (FNF).
      !   890706  Cosmetic changes to prologue.  (WRB)
      !   890911  Removed unnecessary intrinsics.  (WRB)
      !   891009  Removed unreferenced statement label.  (WRB)
      !   891214  Prologue converted to Version 4.0 format.  (BAB)
      !   900309  Added COMP to list of routines called.  (FNF)
      !   900315  Revised prologue and improved some output formats.  (FNF)
      !   900316  Deleted INCFD tests because some compilers object to them,
      !           and made additional minor cosmetic changes.  (FNF)
      !   900322  Made miscellaneous cosmetic changes.  (FNF)
      !   910708  Minor modifications in use of KPRINT.  (WRB)
      !   930504  Removed parens from constants in WRITE statements.  (FNF)

      !
      !  Declare arguments.
      !
      integer Lout, Kprint
      logical Fail
      !
      !  DECLARATIONS.
      !
      integer i, ierr, kontrl, n, nerr, next(2)
      double precision d(10), dum(2), f(10), temp, x(10)
      logical skip
      !
      !  INITIALIZE.
      !
      parameter(n=10)

      nerr = 0
      !
      ! call xgetf(kontrl)
      ! if ( Kprint<=2 ) then
      !    call xsetf(0)
      ! else
      !    call xsetf(1)
      ! endif
      !
      if (Kprint >= 3) write (Lout, 99001)
      !
      !  FORMATS.
      !
99001 format('1'//10x, 'TEST ERROR RETURNS')
      if (Kprint >= 2) write (Lout, 99002)
99002 format(//10x, 'DEVERK RESULTS'/10x, '--------------')
      !
      !  FIRST, TEST DCHFEV AND DCHFDV.
      !
      if (Kprint >= 3) write (Lout, 99005) - 1
      call dchfev(0.d0, 1.d0, 3.d0, 7.d0, 3.d0, 6.d0, 0, dum, dum, next, ierr)
      if (.not. comp(ierr, -1, Lout, Kprint)) nerr = nerr + 1
      !
      if (Kprint >= 3) write (Lout, 99005) - 2
      call dchfev(1.d0, 1.d0, 3.d0, 7.d0, 3.d0, 6.d0, 1, dum, dum, next, ierr)
      if (.not. comp(ierr, -2, Lout, Kprint)) nerr = nerr + 1
      !
      if (Kprint >= 3) write (Lout, 99005) - 1
      call dchfdv(0.d0, 1.d0, 3.d0, 7.d0, 3.d0, 6.d0, 0, dum, dum, dum, next, ierr)
      if (.not. comp(ierr, -1, Lout, Kprint)) nerr = nerr + 1
      !
      if (Kprint >= 3) write (Lout, 99005) - 2
      call dchfdv(1.d0, 1.d0, 3.d0, 7.d0, 3.d0, 6.d0, 1, dum, dum, dum, next, ierr)
      if (.not. comp(ierr, -2, Lout, Kprint)) nerr = nerr + 1
      !
      !  SET UP PCH DEFINITION.
      !
      do i = 1, n
         x(i) = i
         f(i) = i + 2
         d(i) = 1.d0
      end do
      !
      !  SWAP POINTS 4 AND 7, SO X-ARRAY IS OUT OF ORDER.
      !
      temp = x(4)
      x(4) = x(7)
      x(7) = temp
      !
      !  NOW, TEST DPCHFE AND DPCHFD.
      !
      if (Kprint >= 3) write (Lout, 99005) - 1
      skip = .false.
      call dpchfe(1, x, f, d, 1, skip, 0, dum, dum, ierr)
      if (.not. comp(ierr, -1, Lout, Kprint)) nerr = nerr + 1
      !
      if (Kprint >= 3) write (Lout, 99005) - 3
      skip = .false.
      call dpchfe(n, x, f, d, 1, skip, 0, dum, dum, ierr)
      if (.not. comp(ierr, -3, Lout, Kprint)) nerr = nerr + 1
      !
      if (Kprint >= 3) write (Lout, 99005) - 4
      skip = .true.
      call dpchfe(n, x, f, d, 1, skip, 0, dum, dum, ierr)
      if (.not. comp(ierr, -4, Lout, Kprint)) nerr = nerr + 1
      !
      if (Kprint >= 3) write (Lout, 99005) - 1
      skip = .false.
      call dpchfd(1, x, f, d, 1, skip, 0, dum, dum, dum, ierr)
      if (.not. comp(ierr, -1, Lout, Kprint)) nerr = nerr + 1
      !
      if (Kprint >= 3) write (Lout, 99005) - 3
      skip = .false.
      call dpchfd(n, x, f, d, 1, skip, 0, dum, dum, dum, ierr)
      if (.not. comp(ierr, -3, Lout, Kprint)) nerr = nerr + 1
      !
      if (Kprint >= 3) write (Lout, 99005) - 4
      skip = .true.
      call dpchfd(n, x, f, d, 1, skip, 0, dum, dum, dum, ierr)
      if (.not. comp(ierr, -4, Lout, Kprint)) nerr = nerr + 1
      !
      !  SUMMARIZE RESULTS.
      !
      !  if ( Kprint>2 ) call xerdmp
      if (nerr == 0) then
         Fail = .false.
         if (Kprint >= 2) write (Lout, 99003)
99003    format(/' ALL ERROR RETURNS OK.')
      else
         Fail = .true.
         if (Kprint >= 2) write (Lout, 99004) nerr
99004    format(//' ***** TROUBLE IN DEVERK *****'//5x, i5,             &
                      &' TESTS FAILED TO GIVE EXPECTED RESULTS.')
      end if
      !
      !  TERMINATE.
      !
      ! call xsetf(kontrl)
      return
99005 format(/' THIS CALL SHOULD RETURN IERR =', i3)

   end

   subroutine devpck(Lout, Kprint, x, y, f, Fx, Fy, Xe, Ye, Fe, De, Fe2, Fail)
      implicit none
      !***PURPOSE  Test usage of increment argument in DPCHFD and DPCHFE for
      !            DPCHQ1.
      !***AUTHOR  Fritsch, F. N., (LLNL)
      !***DESCRIPTION
      !
      ! ---- CODE TO TEST USAGE OF INCREMENT ARGUMENT IN DPCHFD AND DPCHFE ---
      !
      !     EVALUATES A BICUBIC FUNCTION AND ITS FIRST PARTIAL DERIVATIVES
      !     ON A 4X6 MESH CONTAINED IN A 10X10 ARRAY.
      !
      !     INTERPOLATION OF THESE DATA ALONG MESH LINES IN EITHER DIMENSION
      !     SHOULD AGREE WITH CORRECT FUNCTION WITHIN ROUNDOFF ERROR.
      !
      !     ARRAYS ARE ARGUMENTS ONLY TO ALLOW SHARING STORAGE WITH OTHER
      !     TEST ROUTINES.
      !
      !     NOTE:  RUN WITH KPRINT=4 FOR FULL GORY DETAILS (10 PAGES WORTH).
      !
      !
      !     FORTRAN INTRINSICS USED:  ABS.
      !     FORTRAN LIBRARY ROUTINES USED:  (WRITE).
      !     SLATEC LIBRARY ROUTINES USED:  DPCHFD, DPCHFE, D1MACH.
      !
      !***REVISION HISTORY  (YYMMDD)
      !   820601  DATE WRITTEN
      !   820714  CONVERTED TO QUICK CHECK FOR SLATEC LIBRARY.
      !   820715  1. CORRECTED SOME FORMATS.
      !           2. ADDED CALL TO D1MACH TO SET MACHEP.
      !   890406  1. Modified to make sure final elements of X and XE
      !             agree, to avoid possible failure due to roundoff
      !             error.
      !           2. Added printout of TOL in case of failure.
      !           3. Removed unnecessary IMPLICIT declaration.
      !           4. Corrected a few S.P. constants to D.P.
      !           5. Minor cosmetic changes.
      !   890706  Cosmetic changes to prologue.  (WRB)
      !   890911  Removed unnecessary intrinsics.  (WRB)
      !   891004  Cosmetic changes to prologue.  (WRB)
      !   891214  Prologue converted to Version 4.0 format.  (BAB)
      !   900315  Revised prologue and improved some output formats.  (FNF)
      !   900316  Additional minor cosmetic changes.  (FNF)
      !   900321  Made miscellaneous cosmetic changes.  (FNF)
      !   901130  Made many changes to output:  (FNF)
      !           1. Reduced amount of output for KPRINT=3.  (Now need to
      !              use KPRINT=4 for full output.)
      !           2. Added 1P's to formats and revised some to reduce maximum
      !              line length.
      !   910708  Minor modifications in use of KPRINT.  (WRB)
      !   930317  Improved output formats.  (FNF)

      !
      !  Declare arguments.
      !
      integer Lout, Kprint
      logical Fail
      double precision x(10), y(10), f(10, 10), Fx(10, 10), Fy(10, 10) &
                     & , Xe(51), Ye(51), Fe(51), De(51), Fe2(51)
      !
      !  DECLARATIONS.
      !
      integer i, ier2, ierr, inc, j, k, ne, nerr, nmax, nx, ny
      logical faild, faile, failoc, skip
      double precision dermax, derr, dtrue, dx, fdiff, fdifmx,    &
                     & fermax, ferr, ftrue, machep, tol, pdermx,  &
                     & pdifmx, pfermx, zero
      !
      !  DEFINE TEST FUNCTION AND DERIVATIVES.
      !
      double precision ax, ay, fcn, dfdx, dfdy
      fcn(ax, ay) = ax*(ay*ay)*(ax*ax + 1.d0)
      dfdx(ax, ay) = (ay*ay)*(3.d0*ax*ax + 1.d0)
      dfdy(ax, ay) = 2.d0*ax*ay*(ax*ax + 1.d0)
      !
      data nmax/10/, nx/4/, ny/6/
      data ne/51/
      data zero/0.d0/
      !
      !  INITIALIZE.
      !

      machep = d1mach4
      !       Following tolerance is looser than S.P. version to avoid
      !       spurious failures on some systems.
      tol = 25.d0*machep
      !
      Fail = .false.
      !
      !  SET UP 4-BY-6 MESH IN A 10-BY-10 ARRAY:
      !     X =  0.25(0.25)1.   ;
      !     Y = -0.75(0.5 )1.75 .
      !
      do i = 1, nx - 1
         x(i) = 0.25d0*i
      end do
      x(nx) = 1.d0
      do j = 1, ny
         y(j) = 0.5d0*j - 1.25d0
         do i = 1, nx
            f(i, j) = fcn(x(i), y(j))
            Fx(i, j) = dfdx(x(i), y(j))
            Fy(i, j) = dfdy(x(i), y(j))
         end do
      end do
      !
      !  SET UP EVALUATION POINTS:
      !     XE =  0.(0.02)1. ;
      !     YE = -2.(0.08)2. .
      !
      dx = 1.d0/(ne - 1)
      do k = 1, ne - 1
         Xe(k) = dx*(k - 1)
         Ye(k) = 4.d0*Xe(k) - 2.d0
      end do
      Xe(ne) = 1.d0
      Ye(ne) = 2.d0
      !
      if (Kprint >= 3) write (Lout, 99001)
      !
      !  FORMATS.
      !
99001 format('1'//10x, 'TEST DPCHFE AND DPCHFD')
      if (Kprint >= 2) write (Lout, 99002)
99002 format(//10x, 'DEVPCK RESULTS'/10x, '--------------')
      !
      !  EVALUATE ON HORIZONTAL MESH LINES (Y FIXED, X RUNNING) ..............
      !
      nerr = 0
      inc = 1
      skip = .false.
      do j = 1, ny
         !        --------------------------------------------------------------
         call dpchfd(nx, x, f(1, j), Fx(1, j), inc, skip, ne, Xe, Fe, De, ierr)
         !        --------------------------------------------------------------
         if (Kprint >= 3) write (Lout, 99003) inc, 'J', j, 'Y',      &
                               & y(j), ierr
         if (ierr < 0) then
            !
            failoc = .true.
            if (Kprint >= 2) write (Lout, 99011) ierr
         else
            if (Kprint > 3) write (Lout, 99004) 'X'
            !
            !        DPCHFE SHOULD AGREE EXACTLY WITH DPCHFD.
            !
            !        -----------------------------------------------------------
            call dpchfe(nx, x, f(1, j), Fx(1, j), inc, skip, ne, Xe, Fe2, ier2)
            !        -----------------------------------------------------------
            !
            do k = 1, ne
               ftrue = fcn(Xe(k), y(j))
               ferr = Fe(k) - ftrue
               dtrue = dfdx(Xe(k), y(j))
               derr = De(k) - dtrue
               if (Kprint > 3) write (Lout, 99005) Xe(k), ftrue,       &
                                    & Fe(k), ferr, dtrue, De(k),    &
                                    & derr
               if (k == 1) then
                  !              INITIALIZE.
                  fermax = abs(ferr)
                  pfermx = Xe(1)
                  dermax = abs(derr)
                  pdermx = Xe(1)
                  fdifmx = abs(Fe2(1) - Fe(1))
                  pdifmx = Xe(1)
               else
                  !              SELECT.
                  ferr = abs(ferr)
                  if (ferr > fermax) then
                     fermax = ferr
                     pfermx = Xe(k)
                  end if
                  derr = abs(derr)
                  if (derr > dermax) then
                     dermax = derr
                     pdermx = Xe(k)
                  end if
                  fdiff = abs(Fe2(k) - Fe(k))
                  if (fdiff > fdifmx) then
                     fdifmx = fdiff
                     pdifmx = Xe(k)
                  end if
               end if
            end do
            !
            faild = (fermax > tol) .or. (dermax > tol)
            faile = fdifmx /= zero
            failoc = faild .or. faile .or. (ierr /= 13) .or. (ier2 /= ierr)
            !
            if (failoc .and. (Kprint >= 2)) write (Lout, 99006) 'J', j,&
                &'Y', y(j)
            !
            if ((Kprint >= 3) .or. (faild .and. (Kprint == 2)))           &
               & write (Lout, 99007) fermax, pfermx, dermax, pdermx
            if (faild .and. (Kprint >= 2)) write (Lout, 99010) tol
            !
            if ((Kprint >= 3) .or. (faile .and. (Kprint == 2)))           &
               & write (Lout, 99008) fdifmx, pdifmx
            !
            if ((ierr /= 13) .and. (Kprint >= 2)) write (Lout, 99009) 'D',&
               & ierr, 13
            !
            if ((ier2 /= ierr) .and. (Kprint >= 2)) write (Lout, 99009)    &
               & 'E', ier2, ierr
         end if
         !
         if (failoc) nerr = nerr + 1
         Fail = Fail .or. failoc
      end do
      !
      if (Kprint >= 2) then
         if (nerr > 0) then
            write (Lout, 99012) nerr, 'J'
         else
            write (Lout, 99013) 'J'
         end if
      end if
      !
      !  EVALUATE ON VERTICAL MESH LINES (X FIXED, Y RUNNING) ................
      !
      nerr = 0
      inc = nmax
      skip = .false.
      do i = 1, nx
         !        --------------------------------------------------------------
         call dpchfd(ny, y, f(i, 1), Fy(i, 1), inc, skip, ne, Ye, Fe, De, ierr)
         !        --------------------------------------------------------------
         if (Kprint >= 3) write (Lout, 99003) inc, 'I', i, 'X',      &
                               & x(i), ierr
         if (ierr < 0) then
            !
            failoc = .true.
            if (Kprint >= 2) write (Lout, 99011) ierr
         else
            if (Kprint > 3) write (Lout, 99004) 'Y'
            !
            !        DPCHFE SHOULD AGREE EXACTLY WITH DPCHFD.
            !
            !        -----------------------------------------------------------
            call dpchfe(ny, y, f(i, 1), Fy(i, 1), inc, skip, ne, Ye, Fe2, ier2)
            !        -----------------------------------------------------------
            !
            do k = 1, ne
               ftrue = fcn(x(i), Ye(k))
               ferr = Fe(k) - ftrue
               dtrue = dfdy(x(i), Ye(k))
               derr = De(k) - dtrue
               if (Kprint > 3) write (Lout, 99005) Ye(k), ftrue,       &
                                    & Fe(k), ferr, dtrue, De(k),    &
                                    & derr
               if (k == 1) then
                  !              INITIALIZE.
                  fermax = abs(ferr)
                  pfermx = Ye(1)
                  dermax = abs(derr)
                  pdermx = Ye(1)
                  fdifmx = abs(Fe2(1) - Fe(1))
                  pdifmx = Ye(1)
               else
                  !              SELECT.
                  ferr = abs(ferr)
                  if (ferr > fermax) then
                     fermax = ferr
                     pfermx = Ye(k)
                  end if
                  derr = abs(derr)
                  if (derr > dermax) then
                     dermax = derr
                     pdermx = Ye(k)
                  end if
                  fdiff = abs(Fe2(k) - Fe(k))
                  if (fdiff > fdifmx) then
                     fdifmx = fdiff
                     pdifmx = Ye(k)
                  end if
               end if
            end do
            !
            faild = (fermax > tol) .or. (dermax > tol)
            faile = fdifmx /= zero
            failoc = faild .or. faile .or. (ierr /= 20) .or. (ier2 /= ierr)
            !
            if (failoc .and. (Kprint >= 2)) write (Lout, 99006) 'I', i,&
                &'X', x(i)
            !
            if ((Kprint >= 3) .or. (faild .and. (Kprint == 2)))           &
               & write (Lout, 99007) fermax, pfermx, dermax, pdermx
            if (faild .and. (Kprint >= 2)) write (Lout, 99010) tol
            !
            if ((Kprint >= 3) .or. (faile .and. (Kprint == 2)))           &
               & write (Lout, 99008) fdifmx, pdifmx
            !
            if ((ierr /= 20) .and. (Kprint >= 2)) write (Lout, 99009) 'D',&
               & ierr, 20
            !
            if ((ier2 /= ierr) .and. (Kprint >= 2)) write (Lout, 99009)    &
               & 'E', ier2, ierr
         end if
         !
         if (failoc) nerr = nerr + 1
         Fail = Fail .or. failoc
      end do
      !
      if (Kprint >= 2) then
         if (nerr > 0) then
            write (Lout, 99012) nerr, 'I'
         else
            write (Lout, 99013) 'I'
         end if
      end if
      !
      !  TERMINATE.
      !
      return
99003 format(//20x, 'DPCHFD INCREMENT TEST -- INCFD = ', i2/15x, 'ON ', a1,&
                &'-LINE ', i2, ',  ', a1, ' =', f8.4, '  --  IERR =', i3)
99004 format(/3x, a1, 'E', 10x, 'F', 8x, 'FE', 9x, 'DIFF', 13x, 'D', 8x, 'DE', 9x,  &
                &'DIFF')
99005 format(f7.2, 2(2x, 2f10.5, 1p, d15.5, 0p))
99006 format(/' ***** DPCHFD AND/OR DPCHFE FAILED ON ', a1, '-LINE ', i1, &
                &',  ', a1, ' =', f8.4)
99007 format(/19x, '  MAXIMUM ERROR IN FUNCTION =', 1p, 1p, d13.5, 0p,      &
               & ' (AT', f6.2, '),'/33x, 'IN DERIVATIVE =', 1p, d13.5, 0p, ' (AT',&
               & f6.2, ').')
99008 format('  MAXIMUM DIFFERENCE BETWEEN DPCHFE AND DPCHFD =', 1p,    &
               & d13.5, 0p, ' (AT', f6.2, ').')
99009 format(/'  DPCHF', a1, ' RETURNED IERR = ', i2, ' INSTEAD OF ', i2)
99010 format('  *** BOTH SHOULD BE .LE. TOL =', 1p, d12.5, ' ***')
99011 format(//' ***** ERROR ***** DPCHFD RETURNED IERR =', i5//)
99012 format(//' ***** ERROR ***** DPCHFD AND/OR DPCHFE FAILED ON', i2, &
               & 1x, a1, '-LINES.'//)
99013 format(/' DPCHFD AND DPCHFE OK ON ', a1, '-LINES.')

   end

   subroutine dfdtru(x, f, d)
      implicit none
      !***PURPOSE  Compute exact function values for DEVCHK.
      !***AUTHOR  Fritsch, F. N., (LLNL)
      !***DESCRIPTION
      !
      !        COMPUTE EXACT FUNCTION VALUES IN DOUBLE PRECISION.
      !
      !                   F(X) = X*(X+1)*(X-2)
      !
      !***REVISION HISTORY  (YYMMDD)
      !   820601  DATE WRITTEN
      !   890618  REVISION DATE from Version 3.2
      !   890706  Cosmetic changes to prologue.  (WRB)
      !   891214  Prologue converted to Version 4.0 format.  (BAB)
      !   900315  Revised prologue.  (FNF)
      !   900316  Deleted variables ONE and TWO.  (FNF)
      !   900321  Changed name of d.p. version from DFTRUE to DFDTRU.

      double precision x, f, d
      double precision fact1, fact2, xx
      !

      xx = x
      fact1 = xx + 1
      fact2 = xx - 2
      f = xx*fact1*fact2
      d = fact1*fact2 + xx*(fact1 + fact2)
      !

   end

   LOGICAL FUNCTION COMP(Ieract, Ierexp, Lout, Kprint)
      IMPLICIT NONE
      !***PURPOSE  Compare actual and expected values of error flag.
      !***AUTHOR  Fritsch, F. N., (LLNL)
      !***DESCRIPTION
      !
      !     COMPARE ACTUAL VALUE OF IERR WITH EXPECTED VALUE.
      !        PRINT ERROR MESSAGE IF THEY DON'T AGREE.
      !
      !***REVISION HISTORY  (YYMMDD)
      !   820601  DATE WRITTEN
      !   890618  REVISION DATE from Version 3.2
      !   890706  Cosmetic changes to prologue.  (WRB)
      !   891214  Prologue converted to Version 4.0 format.  (BAB)
      !   900315  Revised prologue.  (FNF)
      !   900316  Minor modification to format 5010.  (FNF)
      !   910708  Minor modifications in use of KPRINT.  (WRB)

      INTEGER Ieract, Ierexp, Lout, Kprint

      IF (Ieract == Ierexp) THEN
         COMP = .TRUE.
         IF (Kprint >= 3) WRITE (Lout, 99001)
99001    FORMAT('     OK.')
      ELSE
         COMP = .FALSE.
         IF (Kprint >= 3) WRITE (Lout, 99002) Ieract
99002    FORMAT(' *** COMPARE FAILED -- IERR =', I5)
      END IF

   END

!*****************************************************************************************
!>
!***PURPOSE  Evaluate the B-representation of a B-spline at X for the
!            function value or any of its derivatives.
!***KEYWORDS  DIFFERENTIATION OF B-SPLINE, EVALUATION OF B-SPLINE
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Written by Carl de Boor and modified by D. E. Amos
!
!     Abstract   **** a double precision routine ****
!         DBVALU is the BVALUE function of the reference.
!
!         DBVALU evaluates the B-representation (T,A,N,K) of a B-spline
!         at X for the function value on IDERIV=0 or any of its
!         derivatives on IDERIV=1,2,...,K-1.  Right limiting values
!         (right derivatives) are returned except at the right end
!         point X=T(N+1) where left limiting values are computed.  The
!         spline is defined on T(K) .LE. X .LE. T(N+1).  DBVALU returns
!         a fatal error message when X is outside of this interval.
!
!         To compute left derivatives or left limiting values at a
!         knot T(I), replace N by I-1 and set X=T(I), I=K+1,N+1.
!
!         DBVALU calls DINTRV
!
!     Description of Arguments
!
!         Input      T,A,X are double precision
!          T       - knot vector of length N+K
!          A       - B-spline coefficient vector of length N
!          N       - number of B-spline coefficients
!                    N = sum of knot multiplicities-K
!          K       - order of the B-spline, K .GE. 1
!          IDERIV  - order of the derivative, 0 .LE. IDERIV .LE. K-1
!                    IDERIV = 0 returns the B-spline value
!          X       - argument, T(K) .LE. X .LE. T(N+1)
!          INBV    - an initialization parameter which must be set
!                    to 1 the first time DBVALU is called.
!
!         Output     WORK,DBVALU are double precision
!          INBV    - INBV contains information for efficient process-
!                    ing after the initial call and INBV must not
!                    be changed by the user.  Distinct splines require
!                    distinct INBV parameters.
!          WORK    - work vector of length 3*K.
!          DBVALU  - value of the IDERIV-th derivative at X
!
!     Error Conditions
!         An improper input is a fatal error
!
!***REFERENCES  Carl de Boor, Package for calculating with B-splines,
!                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
!                 pp. 441-472.
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920501  Reformatted the REFERENCES section.  (WRB)

   double precision function dbvalu(t, a, n, k, Ideriv, x, Inbv, Work)
      implicit none

      integer i, Ideriv, iderp1, ihi, ihmkmj, ilo, imk, imkpj, &
         Inbv, ipj, ip1, ip1mj, j, jj, j1, j2, k, kmider, &
         kmj, km1, kpk, mflag, n
      double precision a, fkmj, t, Work, x
      dimension t(*), a(*), Work(*)
      dbvalu = 0.0d0
      if (k < 1) then
         error stop 'K DOES NOT SATISFY K.GE.1'
         return
      elseif (n < k) then
         !
         !
         error stop 'N DOES NOT SATISFY N.GE.K'
         return
      elseif (Ideriv < 0 .or. Ideriv >= k) then
         error stop 'IDERIV DOES NOT SATISFY 0.LE.IDERIV.LT.K'
         return
      else
         kmider = k - Ideriv
         !
         ! *** FIND *I* IN (K,N) SUCH THAT T(I) .LE. X .LT. T(I+1)
         !     (OR, .LE. T(I+1) IF T(I) .LT. T(I+1) = T(N+1)).
         km1 = k - 1
         call dintrv(t, n + 1, x, Inbv, i, mflag)
         if (x < t(k)) then
            error stop 'X IS N0T GREATER THAN OR EQUAL TO T(K)'
            return
         else
            if (mflag /= 0) then
               if (x > t(i)) then
                  error stop 'X IS NOT LESS THAN OR EQUAL TO T(N+1)'
                  return
               else
                  do
5                    if (i == k) then
                        error stop 'A LEFT LIMITING VALUE CANNOT BE OBTAINED AT T(K)'
                        return
                     else
                        i = i - 1
                        if (x == t(i)) cycle
                     end if
                     exit
                  end do
               end if
            end if
            !
            ! *** DIFFERENCE THE COEFFICIENTS *IDERIV* TIMES
            !     WORK(I) = AJ(I), WORK(K+I) = DP(I), WORK(K+K+I) = DM(I), I=1.K
            !
            imk = i - k
            do j = 1, k
               imkpj = imk + j
               Work(j) = a(imkpj)
            end do
            if (Ideriv /= 0) then
               do j = 1, Ideriv
                  kmj = k - j
                  fkmj = kmj
                  do jj = 1, kmj
                     ihi = i + jj
                     ihmkmj = ihi - kmj
                     Work(jj) = (Work(jj + 1) - Work(jj))/(t(ihi) - t(ihmkmj))*fkmj
                  end do
               end do
            end if
            !
            ! *** COMPUTE VALUE AT *X* IN (T(I),(T(I+1)) OF IDERIV-TH DERIVATIVE,
            !     GIVEN ITS RELEVANT B-SPLINE COEFF. IN AJ(1),...,AJ(K-IDERIV).
            if (Ideriv /= km1) then
               ip1 = i + 1
               kpk = k + k
               j1 = k + 1
               j2 = kpk + 1
               do j = 1, kmider
                  ipj = i + j
                  Work(j1) = t(ipj) - x
                  ip1mj = ip1 - j
                  Work(j2) = x - t(ip1mj)
                  j1 = j1 + 1
                  j2 = j2 + 1
               end do
               iderp1 = Ideriv + 1
               do j = iderp1, km1
                  kmj = k - j
                  ilo = kmj
                  do jj = 1, kmj
                     Work(jj) = (Work(jj + 1)*Work(kpk + ilo) + Work(jj)*Work(k + jj)) &
                                /(Work(kpk + ilo) + Work(k + jj))
                     ilo = ilo - 1
                  end do
               end do
            end if
            dbvalu = Work(1)
            return
         end if
      end if
   end function dbvalu

!*****************************************************************************************
!>
!***PURPOSE  Compute the largest integer ILEFT in 1 .LE. ILEFT .LE. LXT
!            such that XT(ILEFT) .LE. X where XT(*) is a subdivision of
!            the X interval.
!***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, SPLINES
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Written by Carl de Boor and modified by D. E. Amos
!
!     Abstract    **** a double precision routine ****
!         DINTRV is the INTERV routine of the reference.
!
!         DINTRV computes the largest integer ILEFT in 1 .LE. ILEFT .LE.
!         LXT such that XT(ILEFT) .LE. X where XT(*) is a subdivision of
!         the X interval.  Precisely,
!
!                      X .LT. XT(1)                1         -1
!         if  XT(I) .LE. X .LT. XT(I+1)  then  ILEFT=I  , MFLAG=0
!           XT(LXT) .LE. X                         LXT        1,
!
!         That is, when multiplicities are present in the break point
!         to the left of X, the largest index is taken for ILEFT.
!
!     Description of Arguments
!
!         Input      XT,X are double precision
!          XT      - XT is a knot or break point vector of length LXT
!          LXT     - length of the XT vector
!          X       - argument
!          ILO     - an initialization parameter which must be set
!                    to 1 the first time the spline array XT is
!                    processed by DINTRV.
!
!         Output
!          ILO     - ILO contains information for efficient process-
!                    ing after the initial call and ILO must not be
!                    changed by the user.  Distinct splines require
!                    distinct ILO parameters.
!          ILEFT   - largest integer satisfying XT(ILEFT) .LE. X
!          MFLAG   - signals when X lies out of bounds
!
!     Error Conditions
!         None
!
!***REFERENCES  Carl de Boor, Package for calculating with B-splines,
!                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
!                 pp. 441-472.
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)

   subroutine dintrv(Xt, Lxt, x, Ilo, Ileft, Mflag)
      implicit none

      integer ihi, Ileft, Ilo, istep, Lxt, Mflag, middle, gt(4)
      double precision x, Xt
      dimension Xt(*)
      gt = 0
      ihi = Ilo + 1
      do
         if (ihi >= Lxt) then
            if (x >= Xt(Lxt)) then
               gt(4) = 1
               exit
            end if
            if (Lxt <= 1) then
               gt(3) = 1
               exit
            end if
            Ilo = Lxt - 1
            ihi = Lxt
         end if
         !
         if (x >= Xt(ihi)) then
            ! *** NOW X .GE. XT(ILO) . FIND UPPER BOUND
            istep = 1
            do while (.true.)
               Ilo = ihi
               ihi = Ilo + istep
               if (ihi >= Lxt) exit
               if (x < Xt(ihi)) then
                  gt(1) = 1
                  exit
               end if
               istep = istep*2
            end do
            if (any(gt == 1)) exit
            if (x >= Xt(Lxt)) then
               gt(4) = 1
               exit
            end if
            ihi = Lxt
         elseif (x >= Xt(Ilo)) then
            gt(2) = 1
            exit
         else
            !
            ! *** NOW X .LT. XT(IHI) . FIND LOWER BOUND
            istep = 1
            do while (.true.)
               ihi = Ilo
               Ilo = ihi - istep
               if (Ilo <= 1) exit
               if (x >= Xt(Ilo)) then
                  gt(1) = 1
                  exit
               end if
               istep = istep*2
            end do
            if (any(gt == 1)) exit
            Ilo = 1
            if (x < Xt(1)) then
               gt(3) = 1
               exit
            end if
         end if
         exit
      end do

      ! *** NOW XT(ILO) .LE. X .LT. XT(IHI) . NARROW THE INTERVAL
      if (gt(4) == 0) then
         if (gt(3) == 0) then
            if (gt(2) == 0) then
100            do while (.true.)
                  middle = (Ilo + ihi)/2
                  if (middle == Ilo) exit
                  !     NOTE. IT IS ASSUMED THAT MIDDLE = ILO IN CASE IHI = ILO+1
                  if (x < Xt(middle)) then
                     ihi = middle
                  else
                     Ilo = middle
                  end if
               end do
            end if
            gt(2) = 0
200         Mflag = 0
            Ileft = Ilo
            return
            ! *** SET OUTPUT AND RETURN
         end if
         gt(3) = 0
300      Mflag = -1
         Ileft = 1
         return
      end if
      gt(4) = 0
400   Mflag = 1
      Ileft = Lxt
   end subroutine dintrv

end program pchip_test
!*******************************************************************************************************

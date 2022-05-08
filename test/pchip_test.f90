!*******************************************************************************************************
!>
!  Tests for PCHIP.

    program pchip_test

    use pchip_module 
    use iso_fortran_env, only: Lun => output_unit, wp => real64

    implicit none

    integer,parameter :: Kprint = 5 !! printing flag
    integer :: Ipass !! pass/fail flag

    real(wp),parameter :: d1mach4 = epsilon(1.0_wp)  !! `d1mach(4)` -- the largest relative spacing

    ! run all the tests:

    call dpchq1(Lun,Kprint,Ipass); if (ipass==0) error stop 'test dpchq1 failed'
    call dpchq2(Lun,Kprint,Ipass); if (ipass==0) error stop 'test dpchq2 failed'
    call dpchq3(Lun,Kprint,Ipass); if (ipass==0) error stop 'test dpchq3 failed'
    call dpchq4(Lun,Kprint,Ipass); if (ipass==0) error stop 'test dpchq4 failed'

    contains 
!*******************************************************************************************************

    subroutine dpchq1(Lun,Kprint,Ipass)
        implicit none
  !*--DPCHQ15
  !***BEGIN PROLOGUE  DPCHQ1
  !***PURPOSE  Test the PCHIP evaluators DCHFDV, DCHFEV, DPCHFD, DPCHFE.
  !***LIBRARY   SLATEC (PCHIP)
  !***TYPE      DOUBLE PRECISION (PCHQK1-S, DPCHQ1-D)
  !***KEYWORDS  PCHIP EVALUATOR QUICK CHECK
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
  !***ROUTINES CALLED  DEVCHK, DEVERK, DEVPCK
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
  !***END PROLOGUE  DPCHQ1
  !
  !  Declare arguments.
  !
        integer Lun , Kprint , Ipass
  !
  !  DECLARE LOCAL VARIABLES.
  !
        integer i1 , i2 , i3 , i4 , i5 , i6 , i7 , i8 , i9 , ifail , npts
        double precision work(4000)
        logical fail
  !
  !***FIRST EXECUTABLE STATEMENT  DPCHQ1
        if ( Kprint>=2 ) write (Lun,99001) Kprint
  !
  !  FORMATS.
  !
  99001 format ('1'/' ------------ DPCHIP QUICK CHECK OUTPUT',            &
               &' ------------'//20x,'( KPRINT =',i2,' )')
  !
  !  TEST DCHFDV AND DCHFEV.
  !
        ifail = 0
        npts = 1000
        i1 = 1 + npts
        i2 = i1 + npts
        i3 = i2 + npts
        call devchk(Lun,Kprint,npts,work(1),work(i1),work(i2),work(i3),   &
                  & fail)
        if ( fail ) ifail = ifail + 1
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
        call devpck(Lun,Kprint,work(1),work(i1),work(i2),work(i3),work(i4)&
                  & ,work(i5),work(i6),work(i7),work(i8),work(i9),fail)
        if ( fail ) ifail = ifail + 2
  !
  !  TEST ERROR RETURNS.
  !
        call deverk(Lun,Kprint,fail)
        if ( fail ) ifail = ifail + 4
  !
  !  PRINT SUMMARY AND TERMINATE.
  !     At this point, IFAIL has the following value:
  !        IFAIL = 0  IF ALL TESTS PASSED.
  !        IFAIL BETWEEN 1 AND 7 IS THE SUM OF:
  !           IFAIL=1  IF SINGLE CUBIC  TEST FAILED. (SEE DEVCHK OUTPUT.)
  !           IFAIL=2  IF DPCHFD/DPCHFE TEST FAILED. (SEE DEVPCK OUTPUT.)
  !           IFAIL=4  IF ERROR RETURN  TEST FAILED. (SEE DEVERK OUTPUT.)
  !
        if ( (Kprint>=2) .and. (ifail/=0) ) write (Lun,99002) ifail
  99002 format (/' *** TROUBLE ***',i5,' EVALUATION TESTS FAILED.')
  !
        if ( ifail==0 ) then
           Ipass = 1
           if ( Kprint>=2 ) write (Lun,99003)
  99003    format (/' ------------ DPCHIP PASSED  ALL EVALUATION TESTS',  &
                  &' ------------')
        else
           Ipass = 0
           if ( Kprint>=1 ) write (Lun,99004)
  99004    format (/' ************ DPCHIP FAILED SOME EVALUATION TESTS',  &
                  &' ************')
        endif
  !
        return
  !------------- LAST LINE OF DPCHQ1 FOLLOWS -----------------------------
        end

        subroutine dpchq2(Lun,Kprint,Ipass)
        implicit none
  !*--DPCHQ2136
  !***BEGIN PROLOGUE  DPCHQ2
  !***PURPOSE  Test the PCHIP integrators DPCHIA and DPCHID.
  !***LIBRARY   SLATEC (PCHIP)
  !***TYPE      DOUBLE PRECISION (PCHQK2-S, DPCHQ2-D)
  !***KEYWORDS  PCHIP INTEGRATOR QUICK CHECK
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
  !***ROUTINES CALLED  D1MACH, DPCHIA
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
  !***END PROLOGUE  DPCHQ2
  !
  !  Declare arguments.
  !
        integer Lun , Kprint , Ipass
  !
  !  DECLARE VARIABLES.
  !
        integer i , ierexp(17) , ierr , ifail , n , npairs
        double precision a(17) , b(17) , calc , d(7) , errmax , error ,   &
                       & f(7) , machep , one , three , thrqtr , tol ,     &
                       & true , two , x(7)
        logical fail , skip
  !
  !  DECLARE EXTERNALS.
  !
  !      double precision dpchia , d1mach
  !
  !  DEFINE TEST FUNCTIONS.
  !
        double precision ax , fcn , deriv , antder
        fcn(ax) = three*ax*ax*(ax-two)
        deriv(ax) = three*ax*(two*(ax-two)+ax)
        antder(ax) = ax**3*(thrqtr*ax-two)
  !
  !  INITIALIZE.
  !
        data thrqtr/0.75d0/ , one/1.d0/ , two/2.d0/ , three/3.d0/
        data n/7/
        data x/ - 4.d0 , -2.d0 , -0.9d0 , 0.d0 , 0.9d0 , 2.d0 , 4.d0/
        data npairs/17/
        data a/ - 3.0d0 , 3.0d0 , -0.5d0 , -0.5d0 , -0.5d0 , -4.0d0 ,     &
           & -4.0d0 , 3.0d0 , -5.0d0 , -5.0d0 , -6.0d0 , 6.0d0 , -1.5d0 , &
           & -1.5d0 , -3.0d0 , 3.0d0 , 0.5d0/
        data b/3.0d0 , -3.0d0 , 1.0d0 , 2.0d0 , 5.0d0 , -0.5d0 , 4.0d0 ,  &
           & 5.0d0 , -3.0d0 , 5.0d0 , -5.0d0 , 5.0d0 , -0.5d0 , -1.0d0 ,  &
           & -2.5d0 , 3.5d0 , 0.5d0/
        data ierexp/0 , 0 , 0 , 0 , 2 , 0 , 0 , 2 , 1 , 3 , 3 , 3 , 0 ,   &
           & 0 , 0 , 0 , 0/
  !
  !  SET PASS/FAIL TOLERANCE.
  !
  !***FIRST EXECUTABLE STATEMENT  DPCHQ2
        machep = d1mach4
        tol = 100.d0*machep
  !
  !  SET UP PCH FUNCTION DEFINITION.
  !
        do i = 1 , n
           f(i) = fcn(x(i))
           d(i) = deriv(x(i))
        enddo
  !
        if ( Kprint>=3 ) write (Lun,99001)
  !
  !  FORMATS.
  !
  99001 format ('1'//10x,'TEST DPCHIP INTEGRATORS')
        if ( Kprint>=2 ) write (Lun,99002)
  99002 format (//10x,'DPCHQ2 RESULTS'/10x,'--------------')
        if ( Kprint>=3 ) write (Lun,99003) (x(i),f(i),d(i),i=1,n)
  99003 format (//5x,'DATA:'//11x,'X',9x,'F',9x,'D'/(5x,3f10.3))
  !
  !  LOOP OVER (A,B)-PAIRS.
  !
        if ( Kprint>=3 ) write (Lun,99004)
  99004 format (//5x,'TEST RESULTS:'//'    A     B    ERR     TRUE',16x,  &
               &'CALC',15x,'ERROR')
  !
        ifail = 0
  !
        skip = .false.
        do i = 1 , npairs
  !               ---------------------------------------------
           calc = dpchia(n,x,f,d,1,skip,a(i),b(i),ierr)
  !               ---------------------------------------------
           if ( ierr>=0 ) then
              fail = ierr/=ierexp(i)
              true = antder(b(i)) - antder(a(i))
              error = calc - true
              if ( Kprint>=3 ) then
                 if ( fail ) then
                    write (Lun,99005) a(i) , b(i) , ierr , true , calc ,  &
                                    & error , ierexp(i)
  99005             format (2f6.1,i5,1p,2d20.10,d15.5,'  (',i1,') *****')
                 else
                    write (Lun,99010) a(i) , b(i) , ierr , true , calc ,  &
                                    & error
                 endif
              endif
  !
              error = abs(error)/max(one,abs(true))
              if ( fail .or. (error>tol) ) ifail = ifail + 1
              if ( i==1 ) then
                 errmax = error
              else
                 errmax = max(errmax,error)
              endif
           else
              if ( Kprint>=3 ) write (Lun,99010) a(i) , b(i) , ierr
              ifail = ifail + 1
           endif
        enddo
  !
  !  PRINT SUMMARY.
  !
        if ( Kprint>=2 ) then
           write (Lun,99006) errmax , tol
  99006    format (/'  MAXIMUM RELATIVE ERROR IS:',1p,d15.5,              &
                  &',   TOLERANCE:',1p,d15.5)
           if ( ifail/=0 ) write (Lun,99007) ifail
  99007    format (/' *** TROUBLE ***',i5,' INTEGRATION TESTS FAILED.')
        endif
  !
  !  TERMINATE.
  !
        if ( ifail==0 ) then
           Ipass = 1
           if ( Kprint>=2 ) write (Lun,99008)
  99008    format (/' ------------ DPCHIP PASSED  ALL INTEGRATION TESTS', &
                  &' ------------')
        else
           Ipass = 0
           if ( Kprint>=1 ) write (Lun,99009)
  99009    format (/' ************ DPCHIP FAILED SOME INTEGRATION TESTS', &
                  &' ************')
        endif
  !
        return
  99010 format (2f6.1,i5,1p,2d20.10,d15.5)
  !------------- LAST LINE OF DPCHQ2 FOLLOWS -----------------------------
        end

        subroutine dpchq3(Lun,Kprint,Ipass)
        implicit none
  !*--DPCHQ3324
  !***BEGIN PROLOGUE  DPCHQ3
  !***PURPOSE  Test the PCHIP interpolators DPCHIC, DPCHIM, DPCHSP.
  !***LIBRARY   SLATEC (PCHIP)
  !***TYPE      DOUBLE PRECISION (PCHQK3-S, DPCHQ3-D)
  !***KEYWORDS  PCHIP INTERPOLATOR QUICK CHECK
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
  !***ROUTINES CALLED  COMP, D1MACH, DPCHIC, DPCHIM, DPCHSP
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
  !***END PROLOGUE  DPCHQ3
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
        integer Lun , Kprint , Ipass
        double precision d1mach
  !
  !  Declare variables.
  !
        integer i , ic(2) , ierr , ifail , n , nbad , nbadz , nwk
        parameter (n=9,nwk=2*n)
        double precision d(n) , dc(n) , dc5 , dc6 , dm(n) , ds(n) , err , &
                       & f(n) , mone , tol , told , tolz , vc(2) , x(n) , &
                       & wk(nwk) , zero
        parameter (zero=0.0d0,mone=-1.0d0)
        character*6 result
  !
  !  Initialize.
  !
  !       Data.
        data ic/0 , 0/
        data x/ - 2.2d0 , -1.2d0 , -1.0d0 , -0.5d0 , -0.01d0 , 0.5d0 ,    &
           & 1.0d0 , 2.0d0 , 2.2d0/
  !
  !       Results generated on Cray X/MP (9 sign. figs.)
        data dm/0. , 3.80027352d-01 , 7.17253009d-01 , 5.82014161d-01 ,   &
           & 0. , -5.68208031d-01 , -5.13501618d-01 , -7.77910977d-02 ,   &
           & -2.45611117d-03/
        data dc5 , dc6/1.76950158d-02 , -5.69579814d-01/
        data ds/ - 5.16830792d-02 , 5.71455855d-01 , 7.40530225d-01 ,     &
           & 7.63864934d-01 , 1.92614386d-02 , -7.65324380d-01 ,          &
           & -7.28209035d-01 , -7.98445427d-02 , -2.85983446d-02/
  !
  !***FIRST EXECUTABLE STATEMENT  DPCHQ3
        ifail = 0
  !
  !        Set tolerances.
        tol = 10*d1mach4
        told = max(1.0d-7,10*tol)
        tolz = zero
  !
        if ( Kprint>=3 ) write (Lun,99001)
  !
  !  FORMATS.
  !
  99001 format ('1'//10x,'TEST DPCHIP INTERPOLATORS')
        if ( Kprint>=2 ) write (Lun,99002)
  99002 format (//10x,'DPCHQ3 RESULTS'/10x,'--------------')
  !
  !  Set up data.
  !
        do i = 1 , n
           f(i) = exp(-x(i)**2)
        enddo
  !
        if ( Kprint>=3 ) then
           write (Lun,99003)
  99003    format (//5x,'DATA:'/39x,                                      &
                  &'---------- EXPECTED D-VALUES ----------'/12x,'X',9x,  &
                  &'F',18x,'DM',13x,'DC',13x,'DS')
           do i = 1 , 4
              write (Lun,99009) x(i) , f(i) , dm(i) , ds(i)
           enddo
           write (Lun,99010) x(5) , f(5) , dm(5) , dc5 , ds(5)
           write (Lun,99010) x(6) , f(6) , dm(6) , dc6 , ds(6)
           do i = 7 , n
              write (Lun,99009) x(i) , f(i) , dm(i) , ds(i)
           enddo
        endif
  !
  !  Test DPCHIM.
  !
        if ( Kprint>=3 ) write (Lun,99011) 'IM'
  !     --------------------------------
        call dpchim(n,x,f,d,1,ierr)
  !     --------------------------------
  !        Expect IERR=1 (one monotonicity switch).
        if ( Kprint>=3 ) write (Lun,99012) 1
        if ( .not.comp(ierr,1,Lun,Kprint) ) then
           ifail = ifail + 1
        else
           if ( Kprint>=3 ) write (Lun,99013)
           nbad = 0
           nbadz = 0
           do i = 1 , n
              result = '  OK'
  !             D-values should agree with stored values.
  !               (Zero values should agree exactly.)
              if ( dm(i)==zero ) then
                 err = abs(d(i))
                 if ( err>tolz ) then
                    nbadz = nbadz + 1
                    result = '**BADZ'
                 endif
              else
                 err = abs((d(i)-dm(i))/dm(i))
                 if ( err>told ) then
                    nbad = nbad + 1
                    result = '**BAD'
                 endif
              endif
              if ( Kprint>=3 ) write (Lun,99014) i , x(i) , d(i) , err ,  &
                                    & result
           enddo
           if ( (nbadz/=0) .or. (nbad/=0) ) then
              ifail = ifail + 1
              if ( (nbadz/=0) .and. (Kprint>=2) ) write (Lun,99004) nbad
  99004       format (/'    **',i5,                                       &
                     &' DPCHIM RESULTS FAILED TO BE EXACTLY ZERO.')
              if ( (nbad/=0) .and. (Kprint>=2) ) write (Lun,99015) nbad , &
                  &'IM' , told
           else
              if ( Kprint>=2 ) write (Lun,99016) 'IM'
           endif
        endif
  !
  !  Test DPCHIC -- options set to reproduce DPCHIM.
  !
        if ( Kprint>=3 ) write (Lun,99011) 'IC'
  !     --------------------------------------------------------
        call dpchic(ic,vc,zero,n,x,f,dc,1,wk,nwk,ierr)
  !     --------------------------------------------------------
  !        Expect IERR=0 .
        if ( Kprint>=3 ) write (Lun,99012) 0
        if ( .not.comp(ierr,0,Lun,Kprint) ) then
           ifail = ifail + 1
        else
           if ( Kprint>=3 ) write (Lun,99013)
           nbad = 0
           do i = 1 , n
              result = '  OK'
  !           D-values should agree exactly with those computed by DPCHIM.
  !            (To be generous, will only test to machine precision.)
              err = abs(d(i)-dc(i))
              if ( err>tol ) then
                 nbad = nbad + 1
                 result = '**BAD'
              endif
              if ( Kprint>=3 ) write (Lun,99014) i , x(i) , dc(i) , err , &
                                    & result
           enddo
           if ( nbad/=0 ) then
              ifail = ifail + 1
              if ( Kprint>=2 ) write (Lun,99015) nbad , 'IC' , tol
           else
              if ( Kprint>=2 ) write (Lun,99016) 'IC'
           endif
        endif
  !
  !  Test DPCHIC -- default nonzero switch derivatives.
  !
        if ( Kprint>=3 ) write (Lun,99011) 'IC'
  !     -------------------------------------------------------
        call dpchic(ic,vc,mone,n,x,f,d,1,wk,nwk,ierr)
  !     -------------------------------------------------------
  !        Expect IERR=0 .
        if ( Kprint>=3 ) write (Lun,99012) 0
        if ( .not.comp(ierr,0,Lun,Kprint) ) then
           ifail = ifail + 1
        else
           if ( Kprint>=3 ) write (Lun,99013)
           nbad = 0
           nbadz = 0
           do i = 1 , n
              result = '  OK'
  !            D-values should agree exactly with those computed in
  !            previous call, except at points 5 and 6.
              if ( (i<5) .or. (i>6) ) then
                 err = abs(d(i)-dc(i))
                 if ( err>tolz ) then
                    nbadz = nbadz + 1
                    result = '**BADA'
                 endif
              else
                 if ( i==5 ) then
                    err = abs((d(i)-dc5)/dc5)
                 else
                    err = abs((d(i)-dc6)/dc6)
                 endif
                 if ( err>told ) then
                    nbad = nbad + 1
                    result = '**BAD'
                 endif
              endif
              if ( Kprint>=3 ) write (Lun,99014) i , x(i) , d(i) , err ,  &
                                    & result
           enddo
           if ( (nbadz/=0) .or. (nbad/=0) ) then
              ifail = ifail + 1
              if ( (nbadz/=0) .and. (Kprint>=2) ) write (Lun,99005) nbad
  99005       format (/'    **',i5,' DPCHIC RESULTS FAILED TO AGREE WITH',&
                     &' PREVIOUS CALL.')
              if ( (nbad/=0) .and. (Kprint>=2) ) write (Lun,99015) nbad , &
                  &'IC' , told
           else
              if ( Kprint>=2 ) write (Lun,99016) 'IC'
           endif
        endif
  !
  !  Test DPCHSP.
  !
        if ( Kprint>=3 ) write (Lun,99011) 'SP'
  !     -------------------------------------------------
        call dpchsp(ic,vc,n,x,f,d,1,wk,nwk,ierr)
  !     -------------------------------------------------
  !        Expect IERR=0 .
        if ( Kprint>=3 ) write (Lun,99012) 0
        if ( .not.comp(ierr,0,Lun,Kprint) ) then
           ifail = ifail + 1
        else
           if ( Kprint>=3 ) write (Lun,99013)
           nbad = 0
           do i = 1 , n
              result = '  OK'
  !             D-values should agree with stored values.
              err = abs((d(i)-ds(i))/ds(i))
              if ( err>told ) then
                 nbad = nbad + 1
                 result = '**BAD'
              endif
              if ( Kprint>=3 ) write (Lun,99014) i , x(i) , d(i) , err ,  &
                                    & result
           enddo
           if ( nbad/=0 ) then
              ifail = ifail + 1
              if ( Kprint>=2 ) write (Lun,99015) nbad , 'SP' , told
           else
              if ( Kprint>=2 ) write (Lun,99016) 'SP'
           endif
        endif
  !
  !  PRINT SUMMARY AND TERMINATE.
  !
        if ( (Kprint>=2) .and. (ifail/=0) ) write (Lun,99006) ifail
  99006 format (/' *** TROUBLE ***',i5,' INTERPOLATION TESTS FAILED.')
  !
        if ( ifail==0 ) then
           Ipass = 1
           if ( Kprint>=2 ) write (Lun,99007)
  99007    format (/' ------------ DPCHIP PASSED  ALL INTERPOLATION TESTS'&
                 & ,' ------------')
        else
           Ipass = 0
           if ( Kprint>=1 ) write (Lun,99008)
  99008    format (/' ************ DPCHIP FAILED SOME INTERPOLATION TESTS'&
                 & ,' ************')
        endif
  !
        return
  99009 format (5x,f10.2,1p,d15.5,4x,d15.5,15x,d15.5)
  99010 format (5x,f10.2,1p,d15.5,4x,3d15.5)
  99011 format (/5x,'DPCH',a2,' TEST:')
  99012 format (15x,'EXPECT  IERR =',i5)
  99013 format (/9x,'I',7x,'X',9x,'D',13x,'ERR')
  99014 format (5x,i5,f10.2,1p,2d15.5,2x,a)
  99015 format (/'    **',i5,' DPCH',a2,' RESULTS FAILED TOLERANCE TEST.',&
               &'  TOL =',1p,d10.3)
  99016 format (/5x,'  ALL DPCH',a2,' RESULTS OK.')
  !------------- LAST LINE OF DPCHQ3 FOLLOWS -----------------------------
        end

        subroutine dpchq4(Lun,Kprint,Ipass)
        implicit none
  !*--DPCHQ4662
  !***BEGIN PROLOGUE  DPCHQ4
  !***PURPOSE  Test the PCHIP monotonicity checker DPCHCM.
  !***LIBRARY   SLATEC (PCHIP)
  !***TYPE      DOUBLE PRECISION (PCHQK4-S, DPCHQ4-D)
  !***KEYWORDS  PCHIP MONOTONICITY CHECKER QUICK CHECK
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
  !***ROUTINES CALLED  DPCHCM
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
  !***END PROLOGUE  DPCHQ4
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
        integer Lun , Kprint , Ipass
  !
  !  DECLARE VARIABLES.
  !
        integer maxn , maxn2 , maxn3 , nb
        parameter (maxn=16,maxn2=8,maxn3=6,nb=7)
        integer i , ierr , ifail , incfd , ismex1(maxn) , ismex2(maxn2) , &
              & ismex3(maxn3) , ismexb(nb) , ismon(maxn) , k , n , ns(3)
        double precision d(maxn) , db(nb) , f(maxn) , fb(nb) , x(maxn)
        logical skip
  !
  !  DEFINE EXPECTED RESULTS.
  !
        data ismex1/1 , 1 , -1 , 1 , 1 , -1 , 1 , 1 , -1 , 1 , 1 , -1 ,   &
           & 1 , 1 , -1 , 2/
        data ismex2/1 , 2 , 2 , 1 , 2 , 2 , 1 , 2/
        data ismex3/1 , 1 , 1 , 1 , 1 , 1/
        data ismexb/1 , 3 , 1 , -1 , -3 , -1 , 2/
  !
  !  DEFINE TEST DATA.
  !
        data ns/16 , 8 , 6/
  !
  !***FIRST EXECUTABLE STATEMENT  DPCHQ4
        if ( Kprint>=3 ) write (Lun,99001)
  !
  !  FORMATS.
  !
  99001 format ('1'//10x,'TEST DPCHIP MONOTONICITY CHECKER')
        if ( Kprint>=2 ) write (Lun,99002)
  99002 format (//10x,'DPCHQ4 RESULTS'/10x,'--------------')
  !
  !       Define X, F, D.
        do i = 1 , maxn
           x(i) = i
           d(i) = 0.d0
        enddo
        do i = 2 , maxn , 3
           d(i) = 2.d0
        enddo
        do i = 1 , 3
           f(i) = x(i)
           f(i+3) = f(i) + 1.d0
           f(i+6) = f(i+3) + 1.d0
           f(i+9) = f(i+6) + 1.d0
           f(i+12) = f(i+9) + 1.d0
        enddo
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
        do i = 1 , 3
           fb(nb-i+1) = fb(i)
           db(nb-i+1) = -db(i)
        enddo
  !
  !  INITIALIZE.
  !
        ifail = 0
  !
        if ( Kprint>=3 ) then
           write (Lun,99003)
  99003    format (//5x,'DATA:'//9x,'I',4x,'X',5x,'F',5x,'D',5x,'FB',4x,  &
                  &'DB')
           do i = 1 , nb
              write (Lun,99010) i , x(i) , f(i) , d(i) , fb(i) , db(i)
           enddo
           do i = nb + 1 , maxn
              write (Lun,99010) i , x(i) , f(i) , d(i)
           enddo
        endif
  !
  !  TRANSFER POINT FOR SECOND SET OF TESTS.
  !
  !
  !  Loop over a series of values of INCFD.
  !
   100  do incfd = 1 , 3
           n = ns(incfd)
           skip = .false.
  !        -------------------------------------------------
           call dpchcm(n,x,f,d,incfd,skip,ismon,ierr)
  !        -------------------------------------------------
           if ( Kprint>=3 ) write (Lun,99004) incfd , ierr ,              &
                                 & (ismon(i),i=1,n)
  99004    format (/4x,'INCFD =',i2,':  IERR =',i3/15x,'ISMON =',16i3)
           if ( ierr/=0 ) then
              ifail = ifail + 1
              if ( Kprint>=3 ) write (Lun,99011)
           else
              do i = 1 , n
                 if ( incfd==1 ) then
                    if ( ismon(i)/=ismex1(i) ) then
                       ifail = ifail + 1
                       if ( Kprint>=3 ) write (Lun,99012)                 &
                          & (ismex1(k),k=1,n)
                       goto 200
                    endif
                 elseif ( incfd==2 ) then
                    if ( ismon(i)/=ismex2(i) ) then
                       ifail = ifail + 1
                       if ( Kprint>=3 ) write (Lun,99012)                 &
                          & (ismex2(k),k=1,n)
                       goto 200
                    endif
                 elseif ( ismon(i)/=ismex3(i) ) then
                    ifail = ifail + 1
                    if ( Kprint>=3 ) write (Lun,99012) (ismex3(k),k=1,n)
                    goto 200
                 endif
              enddo
           endif
   200  enddo
  !
  !  Test for -1,3,1 bug.
  !
        skip = .false.
  !     ------------------------------------------------
        call dpchcm(nb,x,fb,db,1,skip,ismon,ierr)
  !     ------------------------------------------------
        if ( Kprint>=3 ) write (Lun,99005) ierr , (ismon(i),i=1,nb)
  99005 format (/4x,' Bug test:  IERR =',i3/15x,'ISMON =',7i3)
        if ( ierr/=0 ) then
           ifail = ifail + 1
           if ( Kprint>=3 ) write (Lun,99011)
        else
           do i = 1 , nb
              if ( ismon(i)/=ismexb(i) ) then
                 ifail = ifail + 1
                 if ( Kprint>=3 ) write (Lun,99012) (ismexb(k),k=1,nb)
                 goto 300
              endif
           enddo
        endif
  !
   300  if ( f(1)<0. ) then
  !
  !  PRINT SUMMARY AND TERMINATE.
  !
           if ( (Kprint>=2) .and. (ifail/=0) ) write (Lun,99006) ifail
  99006    format (/' *** TROUBLE ***',i5,' MONOTONICITY TESTS FAILED.')
  !
           if ( ifail==0 ) then
              Ipass = 1
              if ( Kprint>=2 ) write (Lun,99007)
  99007       format (/                                                   &
                    &' ------------ DPCHIP PASSED  ALL MONOTONICITY TESTS'&
                   & ,' ------------')
           else
              Ipass = 0
              if ( Kprint>=1 ) write (Lun,99008)
  99008       format (/                                                   &
                    &' ************ DPCHIP FAILED SOME MONOTONICITY TESTS'&
                   & ,' ************')
           endif
  !
           return
        else
  !
  !  Change sign and do again.
  !
           if ( Kprint>=3 ) write (Lun,99009)
  99009    format (/4x,'Changing sign of data.....')
           do i = 1 , maxn
              f(i) = -f(i)
              d(i) = -d(i)
              if ( ismex1(i)/=2 ) ismex1(i) = -ismex1(i)
           enddo
           do i = 1 , maxn2
              if ( ismex2(i)/=2 ) ismex2(i) = -ismex2(i)
           enddo
           do i = 1 , maxn3
              if ( ismex3(i)/=2 ) ismex3(i) = -ismex3(i)
           enddo
           do i = 1 , nb
              fb(i) = -fb(i)
              db(i) = -db(i)
              if ( ismexb(i)/=2 ) ismexb(i) = -ismexb(i)
           enddo
           goto 100
        endif
  99010 format (5x,i5,5f6.1)
  99011 format (' *** Failed -- bad IERR value.')
  99012 format (' *** Failed -- expect:',16i3)
  !------------- LAST LINE OF DPCHQ4 FOLLOWS -----------------------------
        end
  

        subroutine devchk(Lout,Kprint,Npts,Xev,Fev,Dev,Fev2,Fail)
         implicit none
   !*--DEVCHK5
   !***BEGIN PROLOGUE  DEVCHK
   !***SUBSIDIARY
   !***PURPOSE  Test evaluation accuracy of DCHFDV and DCHFEV for DPCHQ1.
   !***LIBRARY   SLATEC (PCHIP)
   !***TYPE      DOUBLE PRECISION (EVCHCK-S, DEVCHK-D)
   !***KEYWORDS  PCHIP EVALUATOR QUICK CHECK
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
   !***ROUTINES CALLED  D1MACH, DCHFDV, DCHFEV, DFDTRU, RAND
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
   !***END PROLOGUE  DEVCHK
   !
   !  Declare arguments.
   !
         integer Lout , Kprint , Npts
         double precision Xev(*) , Fev(*) , Dev(*) , Fev2(*)
         logical Fail
   !
   !  DECLARATIONS.
   !
         integer i , ierr , iint , next(2) , next2(2) , nint
         double precision aed , aed2 , aedmax , aedmin , aef , aef2 ,      &
                        & aefmax , aefmin , check(2) , checkf(2) ,         &
                        & checkd(2) , d1 , d2 , dermax , dtrue , dx ,      &
                        & eps1 , eps2 , f1 , f2 , fact , fermax , floord , &
                        & floorf , four , ftrue , left(3) , machep , one , &
                        & red , red2 , redmax , redmin , ref , ref2 ,      &
                        & refmax , refmin , right(3) , small , ten , tol1 ,&
                        & tol2 , x1 , x2 , xadmax , xadmin , xafmax ,      &
                        & xafmin , xrdmax , xrdmin , xrfmax , xrfmin , zero
         logical failoc , failnx
   !
         double precision d1mach
   !       The following should stay REAL (no D.P. equivalent).
         !real rand
         !external rand
   !
   !  DEFINE RELATIVE ERROR WITH FLOOR.
   !
         double precision rerr , err , value , floor
         rerr(err,value,floor) = err/max(abs(value),floor)
   !
   !  INITIALIZE.
   !
         data zero/0.d0/ , one/1.d0/ , four/4.d0/ , ten/10.d0/
         data small/1.0d-10/
         data nint/3/
         data left/ - 1.5d0 , 2.0d-10 , 1.0d0/
         data right/2.5d0 , 3.0d-10 , 1.0d+8/
   !
   !***FIRST EXECUTABLE STATEMENT  DEVCHK
         machep = d1mach4
         eps1 = four*machep
         eps2 = ten*machep
   !
         Fail = .false.
   !
         if ( Kprint>=2 ) write (Lout,99001)
   99001 format (//10x,'DEVCHK RESULTS'/10x,'--------------')
   !
   !  CYCLE OVER INTERVALS.
   !
         do iint = 1 , nint
            x1 = left(iint)
            x2 = right(iint)
   !
            fact = max(sqrt(x2-x1),one)
            tol1 = eps1*fact
            tol2 = eps2*fact
   !
   !  COMPUTE AND PRINT ENDPOINT VALUES.
   !
            call dfdtru(x1,f1,d1)
            call dfdtru(x2,f2,d2)
   !
            if ( Kprint>=3 ) then
               if ( iint==1 ) write (Lout,99002)
   !
   !  FORMATS.
   !
   99002       format (/10x,'DCHFDV ACCURACY TEST')
               write (Lout,'(/)')
               write (Lout,99017) 'X1' , x1 , 'X2' , x2
               write (Lout,99017) 'F1' , f1 , 'F2' , f2
               write (Lout,99017) 'D1' , d1 , 'D2' , d2
            endif
   !
            if ( Kprint>=2 ) write (Lout,99003) x1 , x2
   99003    format (/10x,'INTERVAL = (',1p,d12.5,',',d12.5,' ):')
   !
   !  COMPUTE FLOORS FOR RELATIVE ERRORS.
   !
            floorf = max(min(abs(f1),abs(f2)),small)
            floord = max(min(abs(d1),abs(d2)),small)
   !
   !  CHECK REPRODUCTION OF ENDPOINT VALUES.
   !
            Xev(1) = x1
            Xev(2) = x2
   !     -----------------------------------------------------------
            call dchfdv(x1,x2,f1,f2,d1,d2,2,Xev,checkf,checkd,next,ierr)
   !     -----------------------------------------------------------
            aef = checkf(1) - f1
            ref = rerr(aef,f1,floorf)
            aef2 = checkf(2) - f2
            ref2 = rerr(aef2,f2,floorf)
            aed = checkd(1) - d1
            red = rerr(aed,d1,floord)
            aed2 = checkd(2) - d2
            red2 = rerr(aed2,d2,floord)
   !
            failoc = max(abs(ref),abs(ref2),abs(red),abs(red2))>tol1
            Fail = Fail .or. failoc
   !
            if ( Kprint>=3 ) then
               write (Lout,99004) next , aef , aef2 , aed , aed2
   99004       format (/' ERRORS AT ENDPOINTS:',40x,'(NEXT =',2i3,')'//1p, &
                     & 4x,'F1:',d13.5,4x,'F2:',d13.5,4x,'D1:',d13.5,4x,    &
                      &'D2:',d13.5)
               write (Lout,99005) ref , ref2 , red , red2
   99005       format (1p,4(7x,d13.5))
            endif
   !
            if ( failoc .and. (Kprint>=2) ) write (Lout,99006)
   99006    format (/' ***** DCHFDV FAILED TO REPRODUCE ENDPOINT VALUES.')
   !
   !  DCHFEV SHOULD AGREE EXACTLY WITH DCHFDV.
   !                     -------
   !     --------------------------------------------------------------
            call dchfev(x1,x2,f1,f2,d1,d2,2,Xev,check,next,ierr)
   !     --------------------------------------------------------------
            failoc = (check(1)/=checkf(1)) .or. (check(2)/=checkf(2))
            Fail = Fail .or. failoc
   !
            if ( failoc .and. (Kprint>=2) ) write (Lout,99007)
   99007    format (/                                                      &
                  &' ***** DCHFEV DOES NOT AGREE WITH DCHFDV AT ENDPOINTS.'&
                 & )
   !
   !  EVALUATE AT NPTS 'UNIFORMLY RANDOM' POINTS IN (X1,X2).
   !     THIS VERSION EXTENDS EVALUATION DOMAIN BY ADDING 4 SUBINTERVALS
   !     TO LEFT AND 6 TO RIGHT OF [X1,X2].
   !
            dx = (x2-x1)/(Npts-10)
            do i = 1 , Npts
               Xev(i) = (x1+(i-5)*dx) + dx*rand(1) !! JW mod - replace with intrinsic random number generator -TODO
            enddo
   !     --------------------------------------------------------
            call dchfdv(x1,x2,f1,f2,d1,d2,Npts,Xev,Fev,Dev,next,ierr)
   !     --------------------------------------------------------
            if ( ierr/=0 ) then
               failoc = .true.
               if ( Kprint>=2 ) write (Lout,99008) ierr
   99008       format (/' ***** ERROR ***** DCHFDV RETURNED IERR =',i5)
            else
   !
   !     CUMULATE LARGEST AND SMALLEST ERRORS FOR SUMMARY.
   !
               do i = 1 , Npts
                  call dfdtru(Xev(i),ftrue,dtrue)
                  aef = Fev(i) - ftrue
                  ref = rerr(aef,ftrue,floorf)
                  aed = Dev(i) - dtrue
                  red = rerr(aed,dtrue,floord)
   !
                  if ( i==1 ) then
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
                     if ( aef<aefmin ) then
                        aefmin = aef
                        xafmin = Xev(i)
                     elseif ( aef>aefmax ) then
                        aefmax = aef
                        xafmax = Xev(i)
                     endif
                     if ( aed<aedmin ) then
                        aedmin = aed
                        xadmin = Xev(i)
                     elseif ( aed>aedmax ) then
                        aedmax = aed
                        xadmax = Xev(i)
                     endif
                     if ( ref<refmin ) then
                        refmin = ref
                        xrfmin = Xev(i)
                     elseif ( ref>refmax ) then
                        refmax = ref
                        xrfmax = Xev(i)
                     endif
                     if ( red<redmin ) then
                        redmin = red
                        xrdmin = Xev(i)
                     elseif ( red>redmax ) then
                        redmax = red
                        xrdmax = Xev(i)
                     endif
                  endif
               enddo
   !
               fermax = max(abs(refmax),abs(refmin))
               dermax = max(abs(redmax),abs(redmin))
   !
               failnx = (next(1)+next(2))/=10
               failoc = failnx .or. (max(fermax,dermax)>tol2)
            endif
            Fail = Fail .or. failoc
   !
   !  PRINT SUMMARY.
   !
            if ( Kprint>=3 ) then
               write (Lout,99009) Npts - 10 , next
   99009       format (/' ERRORS AT ',i5,' INTERIOR POINTS + 10 OUTSIDE:', &
                     & 15x,'(NEXT =',2i3,')'//30x,'FUNCTION',17x,          &
                      &'DERIVATIVE'/15x,2(11x,'ABS',9x,'REL'))
   !
               write (Lout,99018) 'MIN' , aefmin , refmin , aedmin , redmin
               write (Lout,99019) xafmin , xrfmin , xadmin , xrdmin
               write (Lout,99018) 'MAX' , aefmax , refmax , aedmax , redmax
               write (Lout,99019) xafmax , xrfmax , xadmax , xrdmax
            endif
   !
            if ( Kprint>=2 ) then
               if ( failoc ) then
                  if ( fermax>tol2 ) write (Lout,99020) 'F' , fermax , tol2
                  if ( dermax>tol2 ) write (Lout,99020) 'D' , dermax , tol2
                  if ( failnx ) write (Lout,99010) next
   99010          format (/' ***** REPORTED NEXT =',2i5,                   &
                         &'   RATHER THAN    4    6')
               else
                  write (Lout,99011)
   99011          format (/' DCHFDV RESULTS OK.')
               endif
            endif
   !
   !  CHECK THAT DCHFEV AGREES WITH DCHFDV.
   !
   !     -----------------------------------------------------------------
            call dchfev(x1,x2,f1,f2,d1,d2,Npts,Xev,Fev2,next2,ierr)
   !     -----------------------------------------------------------------
            if ( ierr/=0 ) then
               failoc = .true.
               if ( Kprint>=2 ) write (Lout,99012) ierr
   99012       format (/' ***** ERROR ***** DCHFEV RETURNED IERR =',i5)
            else
               aefmax = abs(Fev2(1)-Fev(1))
               xafmax = Xev(1)
               do i = 2 , Npts
                  aef = abs(Fev2(i)-Fev(i))
                  if ( aef>aefmax ) then
                     aefmax = aef
                     xafmax = Xev(i)
                  endif
               enddo
               failnx = (next2(1)/=next(1)) .or. (next2(2)/=next(2))
               failoc = failnx .or. (aefmax/=zero)
               if ( Kprint>=2 ) then
                  if ( failoc ) then
                     write (Lout,99013)
   99013             format (/' ***** DCHFEV DID NOT AGREE WITH DCHFDV:')
                     if ( aefmax/=zero ) write (Lout,99014) aefmax , xafmax
   99014             format (7x,'MAXIMUM DIFFERENCE ',1p,d12.5,            &
                            &'; OCCURRED AT X =',d12.5)
                     if ( failnx ) write (Lout,99015) next2 , next
   99015             format (7x,'REPORTED NEXT =',2i3,'   RATHER THAN ',   &
                           & 2i3)
                  else
                     write (Lout,99016)
   99016             format (/' DCHFEV AGREES WITH DCHFDV.')
                  endif
               endif
            endif
   !
            Fail = Fail .or. failoc
   !
   !  GO BACK FOR ANOTHER INTERVAL.
   !
         enddo
   !
         return
   99017 format (10x,a2,' =',1p,d18.10,5x,a2,' =',d18.10)
   99018 format (/5x,a3,'IMUM ERROR:  ',1p,2d12.4,2x,2d12.4)
   99019 format (5x,'LOCATED AT X =  ',1p,2d12.4,2x,2d12.4)
   99020 format (/' ***** MAXIMUM RELATIVE ERROR IN ',a1,' =',1p,d12.5,    &
               & ','/17x,'EXCEEDS TOLERANCE =',d12.5)
   !------------- LAST LINE OF DEVCHK FOLLOWS -----------------------------
         end

        subroutine deverk(Lout,Kprint,Fail)
         implicit none
   !*--DEVERK5
   !***BEGIN PROLOGUE  DEVERK
   !***SUBSIDIARY
   !***PURPOSE  Test error returns from DPCHIP evaluators for DPCHQ1.
   !***LIBRARY   SLATEC (PCHIP)
   !***TYPE      DOUBLE PRECISION (EVERCK-S, DEVERK-D)
   !***KEYWORDS  PCHIP EVALUATOR QUICK CHECK
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
   !***ROUTINES CALLED  COMP, DCHFDV, DCHFEV, DPCHFD, DPCHFE, XERDMP,
   !                    XGETF, XSETF
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
   !***END PROLOGUE  DEVERK
   !
   !  Declare arguments.
   !
         integer Lout , Kprint
         logical Fail
   !
   !  DECLARATIONS.
   !
         integer i , ierr , kontrl , n , nerr , next(2)
         double precision d(10) , dum(2) , f(10) , temp , x(10)
         logical skip
   !
   !  INITIALIZE.
   !
         parameter (n=10)
   !***FIRST EXECUTABLE STATEMENT  DEVERK
         nerr = 0
   !
         ! call xgetf(kontrl)
         ! if ( Kprint<=2 ) then
         !    call xsetf(0)
         ! else
         !    call xsetf(1)
         ! endif
   !
         if ( Kprint>=3 ) write (Lout,99001)
   !
   !  FORMATS.
   !
   99001 format ('1'//10x,'TEST ERROR RETURNS')
         if ( Kprint>=2 ) write (Lout,99002)
   99002 format (//10x,'DEVERK RESULTS'/10x,'--------------')
   !
   !  FIRST, TEST DCHFEV AND DCHFDV.
   !
         if ( Kprint>=3 ) write (Lout,99005) -1
         call dchfev(0.d0,1.d0,3.d0,7.d0,3.d0,6.d0,0,dum,dum,next,ierr)
         if ( .not.comp(ierr,-1,Lout,Kprint) ) nerr = nerr + 1
   !
         if ( Kprint>=3 ) write (Lout,99005) -2
         call dchfev(1.d0,1.d0,3.d0,7.d0,3.d0,6.d0,1,dum,dum,next,ierr)
         if ( .not.comp(ierr,-2,Lout,Kprint) ) nerr = nerr + 1
   !
         if ( Kprint>=3 ) write (Lout,99005) -1
         call dchfdv(0.d0,1.d0,3.d0,7.d0,3.d0,6.d0,0,dum,dum,dum,next,ierr)
         if ( .not.comp(ierr,-1,Lout,Kprint) ) nerr = nerr + 1
   !
         if ( Kprint>=3 ) write (Lout,99005) -2
         call dchfdv(1.d0,1.d0,3.d0,7.d0,3.d0,6.d0,1,dum,dum,dum,next,ierr)
         if ( .not.comp(ierr,-2,Lout,Kprint) ) nerr = nerr + 1
   !
   !  SET UP PCH DEFINITION.
   !
         do i = 1 , n
            x(i) = i
            f(i) = i + 2
            d(i) = 1.d0
         enddo
   !
   !  SWAP POINTS 4 AND 7, SO X-ARRAY IS OUT OF ORDER.
   !
         temp = x(4)
         x(4) = x(7)
         x(7) = temp
   !
   !  NOW, TEST DPCHFE AND DPCHFD.
   !
         if ( Kprint>=3 ) write (Lout,99005) -1
         skip = .false.
         call dpchfe(1,x,f,d,1,skip,0,dum,dum,ierr)
         if ( .not.comp(ierr,-1,Lout,Kprint) ) nerr = nerr + 1
   !
         if ( Kprint>=3 ) write (Lout,99005) -3
         skip = .false.
         call dpchfe(n,x,f,d,1,skip,0,dum,dum,ierr)
         if ( .not.comp(ierr,-3,Lout,Kprint) ) nerr = nerr + 1
   !
         if ( Kprint>=3 ) write (Lout,99005) -4
         skip = .true.
         call dpchfe(n,x,f,d,1,skip,0,dum,dum,ierr)
         if ( .not.comp(ierr,-4,Lout,Kprint) ) nerr = nerr + 1
   !
         if ( Kprint>=3 ) write (Lout,99005) -1
         skip = .false.
         call dpchfd(1,x,f,d,1,skip,0,dum,dum,dum,ierr)
         if ( .not.comp(ierr,-1,Lout,Kprint) ) nerr = nerr + 1
   !
         if ( Kprint>=3 ) write (Lout,99005) -3
         skip = .false.
         call dpchfd(n,x,f,d,1,skip,0,dum,dum,dum,ierr)
         if ( .not.comp(ierr,-3,Lout,Kprint) ) nerr = nerr + 1
   !
         if ( Kprint>=3 ) write (Lout,99005) -4
         skip = .true.
         call dpchfd(n,x,f,d,1,skip,0,dum,dum,dum,ierr)
         if ( .not.comp(ierr,-4,Lout,Kprint) ) nerr = nerr + 1
   !
   !  SUMMARIZE RESULTS.
   !
       !  if ( Kprint>2 ) call xerdmp
         if ( nerr==0 ) then
            Fail = .false.
            if ( Kprint>=2 ) write (Lout,99003)
   99003    format (/' ALL ERROR RETURNS OK.')
         else
            Fail = .true.
            if ( Kprint>=2 ) write (Lout,99004) nerr
   99004    format (//' ***** TROUBLE IN DEVERK *****'//5x,i5,             &
                   &' TESTS FAILED TO GIVE EXPECTED RESULTS.')
         endif
   !
   !  TERMINATE.
   !
         ! call xsetf(kontrl)
         return
   99005 format (/' THIS CALL SHOULD RETURN IERR =',i3)
   !------------- LAST LINE OF DEVERK FOLLOWS -----------------------------
         end

         subroutine devpck(Lout,Kprint,x,y,f,Fx,Fy,Xe,Ye,Fe,De,Fe2,Fail)
         implicit none
   !*--DEVPCK164
   !***BEGIN PROLOGUE  DEVPCK
   !***SUBSIDIARY
   !***PURPOSE  Test usage of increment argument in DPCHFD and DPCHFE for
   !            DPCHQ1.
   !***LIBRARY   SLATEC (PCHIP)
   !***TYPE      DOUBLE PRECISION (EVPCCK-S, DEVPCK-D)
   !***KEYWORDS  PCHIP EVALUATOR QUICK CHECK
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
   !***ROUTINES CALLED  D1MACH, DPCHFD, DPCHFE
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
   !***END PROLOGUE  DEVPCK
   !
   !  Declare arguments.
   !
         integer Lout , Kprint
         logical Fail
         double precision x(10) , y(10) , f(10,10) , Fx(10,10) , Fy(10,10) &
                        & , Xe(51) , Ye(51) , Fe(51) , De(51) , Fe2(51)
   !
   !  DECLARATIONS.
   !
         integer i , ier2 , ierr , inc , j , k , ne , nerr , nmax , nx , ny
         logical faild , faile , failoc , skip
         double precision dermax , derr , dtrue , dx , fdiff , fdifmx ,    &
                        & fermax , ferr , ftrue , machep , tol , pdermx ,  &
                        & pdifmx , pfermx , zero
         double precision d1mach
   !
   !  DEFINE TEST FUNCTION AND DERIVATIVES.
   !
         double precision ax , ay , fcn , dfdx , dfdy
         fcn(ax,ay) = ax*(ay*ay)*(ax*ax+1.d0)
         dfdx(ax,ay) = (ay*ay)*(3.d0*ax*ax+1.d0)
         dfdy(ax,ay) = 2.d0*ax*ay*(ax*ax+1.d0)
   !
         data nmax/10/ , nx/4/ , ny/6/
         data ne/51/
         data zero/0.d0/
   !
   !  INITIALIZE.
   !
   !***FIRST EXECUTABLE STATEMENT  DEVPCK
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
         do i = 1 , nx - 1
            x(i) = 0.25d0*i
         enddo
         x(nx) = 1.d0
         do j = 1 , ny
            y(j) = 0.5d0*j - 1.25d0
            do i = 1 , nx
               f(i,j) = fcn(x(i),y(j))
               Fx(i,j) = dfdx(x(i),y(j))
               Fy(i,j) = dfdy(x(i),y(j))
            enddo
         enddo
   !
   !  SET UP EVALUATION POINTS:
   !     XE =  0.(0.02)1. ;
   !     YE = -2.(0.08)2. .
   !
         dx = 1.d0/(ne-1)
         do k = 1 , ne - 1
            Xe(k) = dx*(k-1)
            Ye(k) = 4.d0*Xe(k) - 2.d0
         enddo
         Xe(ne) = 1.d0
         Ye(ne) = 2.d0
   !
         if ( Kprint>=3 ) write (Lout,99001)
   !
   !  FORMATS.
   !
   99001 format ('1'//10x,'TEST DPCHFE AND DPCHFD')
         if ( Kprint>=2 ) write (Lout,99002)
   99002 format (//10x,'DEVPCK RESULTS'/10x,'--------------')
   !
   !  EVALUATE ON HORIZONTAL MESH LINES (Y FIXED, X RUNNING) ..............
   !
         nerr = 0
         inc = 1
         skip = .false.
         do j = 1 , ny
   !        --------------------------------------------------------------
            call dpchfd(nx,x,f(1,j),Fx(1,j),inc,skip,ne,Xe,Fe,De,ierr)
   !        --------------------------------------------------------------
            if ( Kprint>=3 ) write (Lout,99003) inc , 'J' , j , 'Y' ,      &
                                  & y(j) , ierr
            if ( ierr<0 ) then
   !
               failoc = .true.
               if ( Kprint>=2 ) write (Lout,99011) ierr
            else
               if ( Kprint>3 ) write (Lout,99004) 'X'
   !
   !        DPCHFE SHOULD AGREE EXACTLY WITH DPCHFD.
   !
   !        -----------------------------------------------------------
               call dpchfe(nx,x,f(1,j),Fx(1,j),inc,skip,ne,Xe,Fe2,ier2)
   !        -----------------------------------------------------------
   !
               do k = 1 , ne
                  ftrue = fcn(Xe(k),y(j))
                  ferr = Fe(k) - ftrue
                  dtrue = dfdx(Xe(k),y(j))
                  derr = De(k) - dtrue
                  if ( Kprint>3 ) write (Lout,99005) Xe(k) , ftrue ,       &
                                       & Fe(k) , ferr , dtrue , De(k) ,    &
                                       & derr
                  if ( k==1 ) then
   !              INITIALIZE.
                     fermax = abs(ferr)
                     pfermx = Xe(1)
                     dermax = abs(derr)
                     pdermx = Xe(1)
                     fdifmx = abs(Fe2(1)-Fe(1))
                     pdifmx = Xe(1)
                  else
   !              SELECT.
                     ferr = abs(ferr)
                     if ( ferr>fermax ) then
                        fermax = ferr
                        pfermx = Xe(k)
                     endif
                     derr = abs(derr)
                     if ( derr>dermax ) then
                        dermax = derr
                        pdermx = Xe(k)
                     endif
                     fdiff = abs(Fe2(k)-Fe(k))
                     if ( fdiff>fdifmx ) then
                        fdifmx = fdiff
                        pdifmx = Xe(k)
                     endif
                  endif
               enddo
   !
               faild = (fermax>tol) .or. (dermax>tol)
               faile = fdifmx/=zero
               failoc = faild .or. faile .or. (ierr/=13) .or. (ier2/=ierr)
   !
               if ( failoc .and. (Kprint>=2) ) write (Lout,99006) 'J' , j ,&
                   &'Y' , y(j)
   !
               if ( (Kprint>=3) .or. (faild .and. (Kprint==2)) )           &
                  & write (Lout,99007) fermax , pfermx , dermax , pdermx
               if ( faild .and. (Kprint>=2) ) write (Lout,99010) tol
   !
               if ( (Kprint>=3) .or. (faile .and. (Kprint==2)) )           &
                  & write (Lout,99008) fdifmx , pdifmx
   !
               if ( (ierr/=13) .and. (Kprint>=2) ) write (Lout,99009) 'D' ,&
                  & ierr , 13
   !
               if ( (ier2/=ierr) .and. (Kprint>=2) ) write (Lout,99009)    &
                  & 'E' , ier2 , ierr
            endif
   !
            if ( failoc ) nerr = nerr + 1
            Fail = Fail .or. failoc
         enddo
   !
         if ( Kprint>=2 ) then
            if ( nerr>0 ) then
               write (Lout,99012) nerr , 'J'
            else
               write (Lout,99013) 'J'
            endif
         endif
   !
   !  EVALUATE ON VERTICAL MESH LINES (X FIXED, Y RUNNING) ................
   !
         nerr = 0
         inc = nmax
         skip = .false.
         do i = 1 , nx
   !        --------------------------------------------------------------
            call dpchfd(ny,y,f(i,1),Fy(i,1),inc,skip,ne,Ye,Fe,De,ierr)
   !        --------------------------------------------------------------
            if ( Kprint>=3 ) write (Lout,99003) inc , 'I' , i , 'X' ,      &
                                  & x(i) , ierr
            if ( ierr<0 ) then
   !
               failoc = .true.
               if ( Kprint>=2 ) write (Lout,99011) ierr
            else
               if ( Kprint>3 ) write (Lout,99004) 'Y'
   !
   !        DPCHFE SHOULD AGREE EXACTLY WITH DPCHFD.
   !
   !        -----------------------------------------------------------
               call dpchfe(ny,y,f(i,1),Fy(i,1),inc,skip,ne,Ye,Fe2,ier2)
   !        -----------------------------------------------------------
   !
               do k = 1 , ne
                  ftrue = fcn(x(i),Ye(k))
                  ferr = Fe(k) - ftrue
                  dtrue = dfdy(x(i),Ye(k))
                  derr = De(k) - dtrue
                  if ( Kprint>3 ) write (Lout,99005) Ye(k) , ftrue ,       &
                                       & Fe(k) , ferr , dtrue , De(k) ,    &
                                       & derr
                  if ( k==1 ) then
   !              INITIALIZE.
                     fermax = abs(ferr)
                     pfermx = Ye(1)
                     dermax = abs(derr)
                     pdermx = Ye(1)
                     fdifmx = abs(Fe2(1)-Fe(1))
                     pdifmx = Ye(1)
                  else
   !              SELECT.
                     ferr = abs(ferr)
                     if ( ferr>fermax ) then
                        fermax = ferr
                        pfermx = Ye(k)
                     endif
                     derr = abs(derr)
                     if ( derr>dermax ) then
                        dermax = derr
                        pdermx = Ye(k)
                     endif
                     fdiff = abs(Fe2(k)-Fe(k))
                     if ( fdiff>fdifmx ) then
                        fdifmx = fdiff
                        pdifmx = Ye(k)
                     endif
                  endif
               enddo
   !
               faild = (fermax>tol) .or. (dermax>tol)
               faile = fdifmx/=zero
               failoc = faild .or. faile .or. (ierr/=20) .or. (ier2/=ierr)
   !
               if ( failoc .and. (Kprint>=2) ) write (Lout,99006) 'I' , i ,&
                   &'X' , x(i)
   !
               if ( (Kprint>=3) .or. (faild .and. (Kprint==2)) )           &
                  & write (Lout,99007) fermax , pfermx , dermax , pdermx
               if ( faild .and. (Kprint>=2) ) write (Lout,99010) tol
   !
               if ( (Kprint>=3) .or. (faile .and. (Kprint==2)) )           &
                  & write (Lout,99008) fdifmx , pdifmx
   !
               if ( (ierr/=20) .and. (Kprint>=2) ) write (Lout,99009) 'D' ,&
                  & ierr , 20
   !
               if ( (ier2/=ierr) .and. (Kprint>=2) ) write (Lout,99009)    &
                  & 'E' , ier2 , ierr
            endif
   !
            if ( failoc ) nerr = nerr + 1
            Fail = Fail .or. failoc
         enddo
   !
         if ( Kprint>=2 ) then
            if ( nerr>0 ) then
               write (Lout,99012) nerr , 'I'
            else
               write (Lout,99013) 'I'
            endif
         endif
   !
   !  TERMINATE.
   !
         return
   99003 format (//20x,'DPCHFD INCREMENT TEST -- INCFD = ',i2/15x,'ON ',a1,&
                &'-LINE ',i2,',  ',a1,' =',f8.4,'  --  IERR =',i3)
   99004 format (/3x,a1,'E',10x,'F',8x,'FE',9x,'DIFF',13x,'D',8x,'DE',9x,  &
                &'DIFF')
   99005 format (f7.2,2(2x,2f10.5,1p,d15.5,0p))
   99006 format (/' ***** DPCHFD AND/OR DPCHFE FAILED ON ',a1,'-LINE ',i1, &
                &',  ',a1,' =',f8.4)
   99007 format (/19x,'  MAXIMUM ERROR IN FUNCTION =',1p,1p,d13.5,0p,      &
               & ' (AT',f6.2,'),'/33x,'IN DERIVATIVE =',1p,d13.5,0p,' (AT',&
               & f6.2,').')
   99008 format ('  MAXIMUM DIFFERENCE BETWEEN DPCHFE AND DPCHFD =',1p,    &
               & d13.5,0p,' (AT',f6.2,').')
   99009 format (/'  DPCHF',a1,' RETURNED IERR = ',i2,' INSTEAD OF ',i2)
   99010 format ('  *** BOTH SHOULD BE .LE. TOL =',1p,d12.5,' ***')
   99011 format (//' ***** ERROR ***** DPCHFD RETURNED IERR =',i5//)
   99012 format (//' ***** ERROR ***** DPCHFD AND/OR DPCHFE FAILED ON',i2, &
               & 1x,a1,'-LINES.'//)
   99013 format (/' DPCHFD AND DPCHFE OK ON ',a1,'-LINES.')
   !------------- LAST LINE OF DEVPCK FOLLOWS -----------------------------
         end

         subroutine dfdtru(x,f,d)
         implicit none
   !*--DFDTRU509
   !***BEGIN PROLOGUE  DFDTRU
   !***SUBSIDIARY
   !***PURPOSE  Compute exact function values for DEVCHK.
   !***LIBRARY   SLATEC (PCHIP)
   !***TYPE      DOUBLE PRECISION (FDTRUE-S, DFDTRU-D)
   !***KEYWORDS  PCHIP EVALUATOR QUICK CHECK
   !***AUTHOR  Fritsch, F. N., (LLNL)
   !***DESCRIPTION
   !
   !        COMPUTE EXACT FUNCTION VALUES IN DOUBLE PRECISION.
   !
   !                   F(X) = X*(X+1)*(X-2)
   !
   !***ROUTINES CALLED  (NONE)
   !***REVISION HISTORY  (YYMMDD)
   !   820601  DATE WRITTEN
   !   890618  REVISION DATE from Version 3.2
   !   890706  Cosmetic changes to prologue.  (WRB)
   !   891214  Prologue converted to Version 4.0 format.  (BAB)
   !   900315  Revised prologue.  (FNF)
   !   900316  Deleted variables ONE and TWO.  (FNF)
   !   900321  Changed name of d.p. version from DFTRUE to DFDTRU.
   !***END PROLOGUE  DFDTRU
         double precision x , f , d
         double precision fact1 , fact2 , xx
   !
   !***FIRST EXECUTABLE STATEMENT  DFDTRU
         xx = x
         fact1 = xx + 1
         fact2 = xx - 2
         f = xx*fact1*fact2
         d = fact1*fact2 + xx*(fact1+fact2)
   !
   !------------- LAST LINE OF DFDTRU FOLLOWS -----------------------------
         end
   
         LOGICAL FUNCTION COMP(Ieract,Ierexp,Lout,Kprint)
         IMPLICIT NONE
   !*--COMP5
   !***BEGIN PROLOGUE  COMP
   !***SUBSIDIARY
   !***PURPOSE  Compare actual and expected values of error flag.
   !***LIBRARY   SLATEC
   !***KEYWORDS  QUICK CHECK SERVICE ROUTINE
   !***AUTHOR  Fritsch, F. N., (LLNL)
   !***DESCRIPTION
   !
   !     COMPARE ACTUAL VALUE OF IERR WITH EXPECTED VALUE.
   !        PRINT ERROR MESSAGE IF THEY DON'T AGREE.
   !
   !***ROUTINES CALLED  (NONE)
   !***REVISION HISTORY  (YYMMDD)
   !   820601  DATE WRITTEN
   !   890618  REVISION DATE from Version 3.2
   !   890706  Cosmetic changes to prologue.  (WRB)
   !   891214  Prologue converted to Version 4.0 format.  (BAB)
   !   900315  Revised prologue.  (FNF)
   !   900316  Minor modification to format 5010.  (FNF)
   !   910708  Minor modifications in use of KPRINT.  (WRB)
   !***END PROLOGUE  COMP
         INTEGER Ieract , Ierexp , Lout , Kprint
   !***FIRST EXECUTABLE STATEMENT  COMP
         IF ( Ieract==Ierexp ) THEN
            COMP = .TRUE.
            IF ( Kprint>=3 ) WRITE (Lout,99001)
   99001    FORMAT ('     OK.')
         ELSE
            COMP = .FALSE.
            IF ( Kprint>=3 ) WRITE (Lout,99002) Ieract
   99002    FORMAT (' *** COMPARE FAILED -- IERR =',I5)
         ENDIF
   !
   !------------- LAST LINE OF COMP FOLLOWS -----------------------------
         END
               
    end program pchip_test
!*******************************************************************************************************

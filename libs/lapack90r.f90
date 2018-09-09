!LAPACK90R -reduced >the set of modules from the followinig files:
!l_auxmod.f90, dblas1.f90, dblas2.f90, dblas3.f90, dla.f90, 
!dor.f90, dge.f90, lapack90.f90, la_dgesv.f90, dsy.f90, la_dsyev.f90


MODULE LA_PRECISION

!  -- LAPACK90 interface driver routine (version 1.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     May 31, 1997

!  DEFINES SINGLE AND DOUBLE PRECISION PARAMETERS, SP AND DP.
!  THESE VALUES ARE COMPILER DEPENDENT.

!  N.B. With some compilers, it may be necessary to split this file into two
!       separate files so that the first module is available for use by the
!       second.

IMPLICIT NONE
INTEGER, PARAMETER :: sp = KIND(1.0), dp = KIND(1.0D0)
!
END MODULE LA_PRECISION


MODULE LA_AUXMOD

!  -- LAPACK90 interface driver routine (version 1.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     May 31, 1997

! ILAENV, DLAMCH and DLASRT have been added to this module.
! For a single-precision version, change wp below to point to sp.

! ELF90 translation by Alan Miller    03-Sep-1997
! Latest revision - 3 November 1997

USE la_precision, ONLY: wp => dp
IMPLICIT NONE

CONTAINS


FUNCTION LSAME( CA, CB ) RESULT(fn_val)
!
!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994
!
!  Purpose
!  =======
!
!  LSAME  tests if CA is the same letter as CB regardless of case.
!
!  Parameters
!  ==========
!
!  CA      (input) CHARACTER*1
!  CB      (input) CHARACTER*1
!           Characters to be compared.
!
!  .. Scalar Arguments ..
CHARACTER (LEN=1), INTENT(IN) :: CA, CB
LOGICAL                       :: fn_val

!  .. Parameters ..
!  INTEGER, PARAMETER :: IOFF=32       ! Not used!
!  .. Local Scalars ..
   INTEGER :: INTA, INTB, ZCODE
!  .. Intrinsic Functions ..
!  INTRINSIC             ICHAR
!
!  .. Executable Statements ..
!
!  Test if the characters are equal
!
   fn_val = CA == CB
!
!  Now test for equivalence
!
   IF ( .NOT. fn_val ) THEN
!
!     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
!     machines, on which ICHAR returns a value with bit 8 set.
!     ICHAR('A') on Prime machines returns 193 which is the same as
!     ICHAR('A') on an EBCDIC machine.
!
      ZCODE = ICHAR( 'Z' )
!
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
!
      IF( ZCODE == 90 .OR. ZCODE == 122 ) THEN
!
!        ASCII is assumed - ZCODE is the ASCII code of either lower or
!        upper case 'Z'.
!
         IF( INTA >= 97 .AND. INTA <= 122 ) INTA = INTA - 32
         IF( INTB >= 97 .AND. INTB <= 122 ) INTB = INTB - 32
!
      ELSE IF( ZCODE == 233 .OR. ZCODE == 169 ) THEN
!
!        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
!        upper case 'Z'.
!
         IF( INTA >= 129 .AND. INTA <= 137 .OR. &
!            INTA >= 145 .AND. INTA <= 153 .OR. &
             INTA >= 162 .AND. INTA <= 169 ) INTA = INTA + 64
         IF( INTB >= 129 .AND. INTB <= 137 .OR. &
             INTB >= 145 .AND. INTB <= 153 .OR. &
             INTB >= 162 .AND. INTB <= 169 ) INTB = INTB + 64
!
      ELSE IF( ZCODE == 218 .OR. ZCODE == 250 ) THEN
!
!        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
!        plus 128 of either lower or upper case 'Z'.
!
         IF( INTA >= 225 .AND. INTA <= 250 ) INTA = INTA - 32
         IF( INTB >= 225 .AND. INTB <= 250 ) INTB = INTB - 32
      END IF
      fn_val = INTA == INTB
   END IF

RETURN
END FUNCTION LSAME


SUBROUTINE ERINFO(LINFO, SRNAME, INFO, ISTAT)

!  -- LAPACK90 interface driver routine (version 1.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     May 31, 1997

!  .. Scalar Arguments ..
CHARACTER( LEN = * ), INTENT(IN)               :: SRNAME
INTEGER             , INTENT(IN)               :: LINFO
INTEGER             , INTENT(IN OUT), OPTIONAL :: INFO
INTEGER             , INTENT(IN), OPTIONAL     :: ISTAT
!
!  .. Executable Statements ..
!
IF ( ( ( LINFO < 0 .AND. LINFO > -200 ) .OR. LINFO > 0 )           &
     &            .AND. .NOT.PRESENT(INFO) ) THEN
   WRITE (*,*) 'Program terminated in LAPACK_90 subroutine ',SRNAME
   WRITE (*,*) 'Error indicator, INFO = ',LINFO
   IF ( PRESENT(ISTAT) )THEN
     IF ( ISTAT /= 0 ) THEN
       IF ( LINFO == -100 ) THEN
         WRITE (*,*) 'The statement ALLOCATE causes STATUS = ', ISTAT
       ELSE
         WRITE (*,*) 'LINFO = ', LINFO, ' not expected'
       END IF
     END IF
   END IF
   STOP
    ELSE IF ( LINFO <= -200 ) THEN
      WRITE(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE(*,*) '*** WARNING, INFO = ', LINFO, ' WARNING ***'
      IF ( LINFO == -200 ) THEN
        WRITE(*,*) 'Could not allocate sufficient workspace for the optimum'
        WRITE(*,*) 'blocksize, hence the routine may not have performed as'
        WRITE(*,*) 'efficiently as possible'
    ELSE
      WRITE(*,*) 'Unexpected warning'
    END IF
      WRITE(*,*) '++++++++++++++++++++++++++++++++++++++++++++++++'
   END IF
   IF ( PRESENT(INFO) ) THEN
     INFO = LINFO
   END IF

RETURN
END SUBROUTINE ERINFO


FUNCTION ilaenv( ispec, NAME, opts, n1, n2, n3, n4 ) RESULT(fn_val)

!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994

! This version, which is compatible with ELF90, is by Alan Miller
! N.B. Arguments opts and n3 were not used in the F77 version.
! Latest revision - 2 September 1997

!     .. Scalar Arguments ..
CHARACTER (LEN = *), INTENT(IN) :: NAME, opts
INTEGER, INTENT(IN)             :: ispec, n1, n2, n3, n4
INTEGER                         :: fn_val
!     ..

!  Purpose
!  =======

!  ILAENV is called from the LAPACK routines to choose problem-dependent
!  parameters for the local environment.  See ISPEC for a description of
!  the parameters.

!  This version provides a set of parameters which should give good,
!  but not optimal, performance on many of the currently available
!  computers.  Users are encouraged to modify this subroutine to set
!  the tuning parameters for their particular machine using the option
!  and problem size information in the arguments.

!  This routine will not function correctly if it is converted to all
!  lower case.  Converting it to all upper case is allowed.

!  Arguments
!  =========

!  ISPEC   (input) INTEGER
!          Specifies the parameter to be returned as the value of ILAENV.
!          = 1: the optimal blocksize; if this value is 1, an unblocked
!               algorithm will give the best performance.
!          = 2: the minimum block size for which the block routine
!               should be used; if the usable block size is less than
!               this value, an unblocked routine should be used.
!          = 3: the crossover point (in a block routine, for N less
!               than this value, an unblocked routine should be used)
!          = 4: the number of shifts, used in the nonsymmetric
!               eigenvalue routines
!          = 5: the minimum column dimension for blocking to be used;
!               rectangular blocks must have dimension at least k by m,
!               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
!          = 6: the crossover point for the SVD (when reducing an m by n
!               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
!               this value, a QR factorization is used first to reduce
!               the matrix to a triangular form.)
!          = 7: the number of processors
!          = 8: the crossover point for the multishift QR and QZ methods
!               for nonsymmetric eigenvalue problems.

!  NAME    (input) CHARACTER*(*)
!          The name of the calling subroutine, in either upper or lower case.

!  OPTS    (input) CHARACTER*(*)
!          The character options to the subroutine NAME, concatenated
!          into a single character string.  For example, UPLO = 'U',
!          TRANS = 'T', and DIAG = 'N' for a triangular routine would
!          be specified as OPTS = 'UTN'.

!  N1      (input) INTEGER
!  N2      (input) INTEGER
!  N3      (input) INTEGER
!  N4      (input) INTEGER
!          Problem dimensions for the subroutine NAME; these may not all
!          be required.

! (ILAENV) (output) INTEGER
!          >= 0: the value of the parameter specified by ISPEC
!          < 0:  if ILAENV = -k, the k-th argument had an illegal value.

!  Further Details
!  ===============

!  The following conventions have been used when calling ILAENV from the
!  LAPACK routines:
!  1)  OPTS is a concatenation of all of the character options to
!      subroutine NAME, in the same order that they appear in the
!      argument list for NAME, even if they are not used in determining
!      the value of the parameter specified by ISPEC.
!  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
!      that they appear in the argument list for NAME.  N1 is used
!      first, N2 second, and so on, and unused problem dimensions are
!      passed a value of -1.
!  3)  The parameter value returned by ILAENV is checked for validity in
!      the calling subroutine.  For example, ILAENV is used to retrieve
!      the optimal blocksize for STRTRI as follows:

!      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
!      IF( NB.LE.1 ) NB = MAX( 1, N )

!  =====================================================================

!     .. Local Scalars ..
LOGICAL           :: cname, sname
CHARACTER (LEN=1) :: c1
CHARACTER (LEN=2) :: c2, c4
CHARACTER (LEN=3) :: c3
CHARACTER (LEN=6) :: subnam
INTEGER           :: i, ic, iz, nb, nbmin, nx
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
!     ..
!     .. Executable Statements ..

SELECT CASE (ispec)
  CASE (1:3)
    GO TO 100
  CASE (4)
    GO TO 400
  CASE (5)
    GO TO 500
  CASE (6)
    GO TO 600
  CASE (7)
    GO TO 700
  CASE (8)
    GO TO 800
  CASE DEFAULT
    fn_val = -1              !     Invalid value for ISPEC
    RETURN
END SELECT

!     Convert NAME to upper case if the first character is lower case.

100 fn_val = 1
subnam = NAME
ic = ICHAR( subnam( 1:1 ) )
iz = ICHAR( 'Z' )
IF( iz == 90 .OR. iz == 122 ) THEN
  
!        ASCII character set
  
  IF( ic >= 97 .AND. ic <= 122 ) THEN
    subnam( 1:1 ) = CHAR( ic-32 )
    DO i = 2, 6
      ic = ICHAR( subnam( i:i ) )
      IF( ic >= 97 .AND. ic <= 122 ) subnam( i:i ) = CHAR( ic-32 )
    END DO
  END IF
  
ELSE IF( iz == 233 .OR. iz == 169 ) THEN
  
!        EBCDIC character set
  
  IF( ( ic >= 129 .AND. ic <= 137 ) .OR.  &
      ( ic >= 145 .AND. ic <= 153 ) .OR.  &
      ( ic >= 162 .AND. ic <= 169 ) ) THEN
    subnam( 1:1 ) = CHAR( ic+64 )
    DO i = 2, 6
      ic = ICHAR( subnam( i:i ) )
      IF( ( ic >= 129 .AND. ic <= 137 ) .OR.  &
          ( ic >= 145 .AND. ic <= 153 ) .OR.  &
          ( ic >= 162 .AND. ic <= 169 ) ) subnam( i:i ) = CHAR( ic+64 )
    END DO
  END IF
  
ELSE IF( iz == 218 .OR. iz == 250 ) THEN
  
!        Prime machines:  ASCII+128
  
  IF( ic >= 225 .AND. ic <= 250 ) THEN
    subnam( 1:1 ) = CHAR( ic-32 )
    DO i = 2, 6
      ic = ICHAR( subnam( i:i ) )
      IF( ic >= 225 .AND. ic <= 250 ) subnam( i:i ) = CHAR( ic-32 )
    END DO
  END IF
END IF

c1 = subnam( 1:1 )
sname = c1 == 'S' .OR. c1 == 'D'
cname = c1 == 'C' .OR. c1 == 'Z'
IF( .NOT.( cname .OR. sname ) ) RETURN
c2 = subnam( 2:3 )
c3 = subnam( 4:6 )
c4 = c3( 2:3 )

IF (ispec == 2) THEN
  GO TO 200
ELSE IF (ispec == 3) THEN
  GO TO 300
END IF

!     ISPEC = 1:  block size

!     In these examples, separate code is provided for setting NB for
!     real and complex.  We assume that NB will take the same value in
!     single or double precision.

nb = 1

IF( c2 == 'GE' ) THEN
  IF( c3 == 'TRF' ) THEN
    IF( sname ) THEN
      nb = 64
    ELSE
      nb = 64
    END IF
  ELSE IF( c3 == 'QRF' .OR. c3 == 'RQF' .OR. c3 == 'LQF' .OR.  &
    c3 == 'QLF' ) THEN
    IF( sname ) THEN
      nb = 32
    ELSE
      nb = 32
    END IF
  ELSE IF( c3 == 'HRD' ) THEN
    IF( sname ) THEN
      nb = 32
    ELSE
      nb = 32
    END IF
  ELSE IF( c3 == 'BRD' ) THEN
    IF( sname ) THEN
      nb = 32
    ELSE
      nb = 32
    END IF
  ELSE IF( c3 == 'TRI' ) THEN
    IF( sname ) THEN
      nb = 64
    ELSE
      nb = 64
    END IF
  END IF
ELSE IF( c2 == 'PO' ) THEN
  IF( c3 == 'TRF' ) THEN
    IF( sname ) THEN
      nb = 64
    ELSE
      nb = 64
    END IF
  END IF
ELSE IF( c2 == 'SY' ) THEN
  IF( c3 == 'TRF' ) THEN
    IF( sname ) THEN
      nb = 64
    ELSE
      nb = 64
    END IF
  ELSE IF( sname .AND. c3 == 'TRD' ) THEN
    nb = 1
  ELSE IF( sname .AND. c3 == 'GST' ) THEN
    nb = 64
  END IF
ELSE IF( cname .AND. c2 == 'HE' ) THEN
  IF( c3 == 'TRF' ) THEN
    nb = 64
  ELSE IF( c3 == 'TRD' ) THEN
    nb = 1
  ELSE IF( c3 == 'GST' ) THEN
    nb = 64
  END IF
ELSE IF( sname .AND. c2 == 'OR' ) THEN
  IF( c3( 1:1 ) == 'G' ) THEN
    IF( c4 == 'QR' .OR. c4 == 'RQ' .OR. c4 == 'LQ' .OR.  &
      c4 == 'QL' .OR. c4 == 'HR' .OR. c4 == 'TR' .OR.c4 == 'BR' ) THEN
      nb = 32
    END IF
  ELSE IF( c3( 1:1 ) == 'M' ) THEN
    IF( c4 == 'QR' .OR. c4 == 'RQ' .OR. c4 == 'LQ' .OR.  &
      c4 == 'QL' .OR. c4 == 'HR' .OR. c4 == 'TR' .OR.c4 == 'BR' ) THEN
      nb = 32
    END IF
  END IF
ELSE IF( cname .AND. c2 == 'UN' ) THEN
  IF( c3( 1:1 ) == 'G' ) THEN
    IF( c4 == 'QR' .OR. c4 == 'RQ' .OR. c4 == 'LQ' .OR.  &
      c4 == 'QL' .OR. c4 == 'HR' .OR. c4 == 'TR' .OR.c4 == 'BR' ) THEN
      nb = 32
    END IF
  ELSE IF( c3( 1:1 ) == 'M' ) THEN
    IF( c4 == 'QR' .OR. c4 == 'RQ' .OR. c4 == 'LQ' .OR.  &
      c4 == 'QL' .OR. c4 == 'HR' .OR. c4 == 'TR' .OR.c4 == 'BR' ) THEN
      nb = 32
    END IF
  END IF
ELSE IF( c2 == 'GB' ) THEN
  IF( c3 == 'TRF' ) THEN
    IF( sname ) THEN
      IF( n4 <= 64 ) THEN
        nb = 1
      ELSE
        nb = 32
      END IF
    ELSE
      IF( n4 <= 64 ) THEN
        nb = 1
      ELSE
        nb = 32
      END IF
    END IF
  END IF
ELSE IF( c2 == 'PB' ) THEN
  IF( c3 == 'TRF' ) THEN
    IF( sname ) THEN
      IF( n2 <= 64 ) THEN
        nb = 1
      ELSE
        nb = 32
      END IF
    ELSE
      IF( n2 <= 64 ) THEN
        nb = 1
      ELSE
        nb = 32
      END IF
    END IF
  END IF
ELSE IF( c2 == 'TR' ) THEN
  IF( c3 == 'TRI' ) THEN
    IF( sname ) THEN
      nb = 64
    ELSE
      nb = 64
    END IF
  END IF
ELSE IF( c2 == 'LA' ) THEN
  IF( c3 == 'UUM' ) THEN
    IF( sname ) THEN
      nb = 64
    ELSE
      nb = 64
    END IF
  END IF
ELSE IF( sname .AND. c2 == 'ST' ) THEN
  IF( c3 == 'EBZ' ) THEN
    nb = 1
  END IF
END IF
fn_val = nb
RETURN

!     ISPEC = 2:  minimum block size

200 nbmin = 2
IF( c2 == 'GE' ) THEN
  IF( c3 == 'QRF' .OR. c3 == 'RQF' .OR. c3 == 'LQF' .OR. c3 == 'QLF' ) THEN
    IF( sname ) THEN
      nbmin = 2
    ELSE
      nbmin = 2
    END IF
  ELSE IF( c3 == 'HRD' ) THEN
    IF( sname ) THEN
      nbmin = 2
    ELSE
      nbmin = 2
    END IF
  ELSE IF( c3 == 'BRD' ) THEN
    IF( sname ) THEN
      nbmin = 2
    ELSE
      nbmin = 2
    END IF
  ELSE IF( c3 == 'TRI' ) THEN
    IF( sname ) THEN
      nbmin = 2
    ELSE
      nbmin = 2
    END IF
  END IF
ELSE IF( c2 == 'SY' ) THEN
  IF( c3 == 'TRF' ) THEN
    IF( sname ) THEN
      nbmin = 8
    ELSE
      nbmin = 8
    END IF
  ELSE IF( sname .AND. c3 == 'TRD' ) THEN
    nbmin = 2
  END IF
ELSE IF( cname .AND. c2 == 'HE' ) THEN
  IF( c3 == 'TRD' ) THEN
    nbmin = 2
  END IF
ELSE IF( sname .AND. c2 == 'OR' ) THEN
  IF( c3( 1:1 ) == 'G' ) THEN
    IF( c4 == 'QR' .OR. c4 == 'RQ' .OR. c4 == 'LQ' .OR.  &
      c4 == 'QL' .OR. c4 == 'HR' .OR. c4 == 'TR' .OR.c4 == 'BR' ) THEN
      nbmin = 2
    END IF
  ELSE IF( c3( 1:1 ) == 'M' ) THEN
    IF( c4 == 'QR' .OR. c4 == 'RQ' .OR. c4 == 'LQ' .OR.  &
      c4 == 'QL' .OR. c4 == 'HR' .OR. c4 == 'TR' .OR.c4 == 'BR' ) THEN
      nbmin = 2
    END IF
  END IF
ELSE IF( cname .AND. c2 == 'UN' ) THEN
  IF( c3( 1:1 ) == 'G' ) THEN
    IF( c4 == 'QR' .OR. c4 == 'RQ' .OR. c4 == 'LQ' .OR.  &
      c4 == 'QL' .OR. c4 == 'HR' .OR. c4 == 'TR' .OR.c4 == 'BR' ) THEN
      nbmin = 2
    END IF
  ELSE IF( c3( 1:1 ) == 'M' ) THEN
    IF( c4 == 'QR' .OR. c4 == 'RQ' .OR. c4 == 'LQ' .OR.  &
      c4 == 'QL' .OR. c4 == 'HR' .OR. c4 == 'TR' .OR.c4 == 'BR' ) THEN
      nbmin = 2
    END IF
  END IF
END IF
fn_val = nbmin
RETURN

!     ISPEC = 3:  crossover point

300 nx = 0
IF( c2 == 'GE' ) THEN
  IF( c3 == 'QRF' .OR. c3 == 'RQF' .OR. c3 == 'LQF' .OR. c3 == 'QLF' ) THEN
    IF( sname ) THEN
      nx = 128
    ELSE
      nx = 128
    END IF
  ELSE IF( c3 == 'HRD' ) THEN
    IF( sname ) THEN
      nx = 128
    ELSE
      nx = 128
    END IF
  ELSE IF( c3 == 'BRD' ) THEN
    IF( sname ) THEN
      nx = 128
    ELSE
      nx = 128
    END IF
  END IF
ELSE IF( c2 == 'SY' ) THEN
  IF( sname .AND. c3 == 'TRD' ) THEN
    nx = 1
  END IF
ELSE IF( cname .AND. c2 == 'HE' ) THEN
  IF( c3 == 'TRD' ) THEN
    nx = 1
  END IF
ELSE IF( sname .AND. c2 == 'OR' ) THEN
  IF( c3( 1:1 ) == 'G' ) THEN
    IF( c4 == 'QR' .OR. c4 == 'RQ' .OR. c4 == 'LQ' .OR.  &
      c4 == 'QL' .OR. c4 == 'HR' .OR. c4 == 'TR' .OR.c4 == 'BR' ) THEN
      nx = 128
    END IF
  END IF
ELSE IF( cname .AND. c2 == 'UN' ) THEN
  IF( c3( 1:1 ) == 'G' ) THEN
    IF( c4 == 'QR' .OR. c4 == 'RQ' .OR. c4 == 'LQ' .OR.  &
      c4 == 'QL' .OR. c4 == 'HR' .OR. c4 == 'TR' .OR.c4 == 'BR' ) THEN
      nx = 128
    END IF
  END IF
END IF
fn_val = nx
RETURN

!     ISPEC = 4:  number of shifts (used by xHSEQR)

400 fn_val = 6
RETURN

!     ISPEC = 5:  minimum column dimension (not used)

500 fn_val = 2

IF (opts == '$#&?') fn_val = 0         ! These two lines are here to make
IF (n3 > 999) fn_val = n3              ! sure that opts and n3 are referenced.
RETURN

!     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)

600 fn_val = INT( REAL( MIN( n1, n2 ) )*1.6E0 )
RETURN

!     ISPEC = 7:  number of processors (not used)

700 fn_val = 1
RETURN

!     ISPEC = 8:  crossover point for multishift (used by xHSEQR)

800 fn_val = 50
RETURN

!     End of ILAENV

END FUNCTION ilaenv


FUNCTION dlamch( cmach ) RESULT(rmach)

!  -- LAPACK auxiliary routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller    02-Sep-1997
! The version of DLAMCH distributed with LAPACK90 causes Lahey's LF90
! compiler (version 3.5) to hang.   This drastically cut down version
! does not, and gives the same results as ELF90 on the original.

!     .. Scalar Arguments ..
CHARACTER (LEN=1), INTENT(IN) :: cmach
REAL (wp)                     :: rmach
!     ..

!  Purpose
!  =======

!  DLAMCH determines double precision machine parameters.

!  Arguments
!  =========

!  CMACH   (input) CHARACTER*1
!          Specifies the value to be returned by DLAMCH:
!          = 'E' or 'e',   DLAMCH := eps
!          = 'S' or 's ,   DLAMCH := sfmin
!          = 'B' or 'b',   DLAMCH := base
!          = 'P' or 'p',   DLAMCH := eps*base
!          = 'N' or 'n',   DLAMCH := t
!          = 'R' or 'r',   DLAMCH := rnd
!          = 'M' or 'm',   DLAMCH := emin
!          = 'U' or 'u',   DLAMCH := rmin
!          = 'L' or 'l',   DLAMCH := emax
!          = 'O' or 'o',   DLAMCH := rmax

!          where

!          eps   = relative machine precision
!          sfmin = safe minimum, such that 1/sfmin does not overflow
!          base  = base of the machine
!          prec  = eps*base
!          t     = number of (base) digits in the mantissa
!          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!          emin  = minimum exponent before (gradual) underflow
!          rmin  = underflow threshold - base**(emin-1)
!          emax  = largest exponent before overflow
!          rmax  = overflow threshold  - (base**emax)*(1-eps)

!     .. Local Scalars ..
LOGICAL, SAVE   :: first = .TRUE.
REAL (wp), SAVE :: base, emax, emin, eps, prec, rmax, rmin, rnd, sfmin, t
REAL (wp)       :: small
!     ..
!     .. External Functions ..
! LOGICAL ::            lsame
! EXTERNAL           lsame
!     ..
!     .. External Subroutines ..
! EXTERNAL           dlamc2
!     ..
!     .. Save statement ..
! SAVE   first, eps, sfmin, base, t, rnd, emin, rmin, emax, rmax, prec
!     ..
!     .. Data statements ..
! DATA               first / .true. /
!     ..
!     .. Executable Statements ..

IF( first ) THEN
  first = .false.
!  CALL dlamc2( beta, it, lrnd, eps, imin, rmin, imax, rmax )
  base = RADIX(1.0_wp)
  t = DIGITS(1.0_wp)
  rnd = 1.0_wp
  eps = EPSILON(1.0_wp)
  prec = eps*base
  emin = MINEXPONENT(1.0_wp)
  emax = MAXEXPONENT(1.0_wp)
  rmin = TINY(1.0_wp)
  rmax = HUGE(1.0_wp)
  sfmin = rmin
  small =1.0_wp / rmax
  IF( small >= sfmin ) THEN

!           Use SMALL plus a bit, to avoid the possibility of rounding
!           causing overflow when computing  1/sfmin.

    sfmin = small*(1.0_wp + eps )
  END IF
END IF

IF( lsame( cmach, 'E' ) ) THEN
  rmach = eps
ELSE IF( lsame( cmach, 'S' ) ) THEN
  rmach = sfmin
ELSE IF( lsame( cmach, 'B' ) ) THEN
  rmach = base
ELSE IF( lsame( cmach, 'P' ) ) THEN
  rmach = prec
ELSE IF( lsame( cmach, 'N' ) ) THEN
  rmach = t
ELSE IF( lsame( cmach, 'R' ) ) THEN
  rmach = rnd
ELSE IF( lsame( cmach, 'M' ) ) THEN
  rmach = emin
ELSE IF( lsame( cmach, 'U' ) ) THEN
  rmach = rmin
ELSE IF( lsame( cmach, 'L' ) ) THEN
  rmach = emax
ELSE IF( lsame( cmach, 'O' ) ) THEN
  rmach = rmax
END IF

RETURN

!     End of DLAMCH

END FUNCTION dlamch



SUBROUTINE dlasrt( id, n, d, info )

!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994

! ELF90 translation by Alan Miller    11-Sep-1997

!     .. Scalar Arguments ..
CHARACTER (LEN=1), INTENT(IN) :: id
INTEGER, INTENT(IN)           :: n
INTEGER, INTENT(OUT)          :: info
!     ..
!     .. Array Arguments ..
REAL (wp), INTENT(IN OUT)     :: d(:)
!     ..

!  Purpose
!  =======

!  Sort the numbers in D in increasing order (if ID = 'I') or
!  in decreasing order (if ID = 'D' ).

!  Use Quick Sort, reverting to Insertion sort on arrays of
!  size <= 20. Dimension of STACK limits N to about 2**32.

!  Arguments
!  =========

!  ID      (input) CHARACTER*1
!          = 'I': sort D in increasing order;
!          = 'D': sort D in decreasing order.

!  N       (input) INTEGER
!          The length of the array D.

!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the array to be sorted.
!          On exit, D has been sorted into increasing order
!          (D(1) <= ... <= D(N) ) or into decreasing order
!          (D(1) >= ... >= D(N) ), depending on ID.

!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value

!  =====================================================================

!     .. Parameters ..
INTEGER, PARAMETER :: select = 20
!     ..
!     .. Local Scalars ..
INTEGER   :: dir, endd, i, j, start, stkpnt
REAL (wp) :: d1, d2, d3, dmnmx, tmp
!     ..
!     .. Local Arrays ..
INTEGER   :: stack( 2, 32 )
!     ..
!     .. External Functions ..
! LOGICAL ::            lsame
! EXTERNAL           lsame
!     ..
!     .. External Subroutines ..
! EXTERNAL           xerbla
!     ..
!     .. Executable Statements ..

!     Test the input paramters.

info = 0
dir = -1
IF( lsame( id, 'D' ) ) THEN
  dir = 0
ELSE IF( lsame( id, 'I' ) ) THEN
  dir = 1
END IF
IF( dir == -1 ) THEN
  info = -1
ELSE IF( n < 0 ) THEN
  info = -2
END IF
IF( info /= 0 ) THEN
  CALL erinfo(-info, 'DLASRT')
  RETURN
END IF

!     Quick return if possible

IF( n <= 1 ) RETURN

stkpnt = 1
stack( 1, 1 ) = 1
stack( 2, 1 ) = n
10 start = stack( 1, stkpnt )
endd = stack( 2, stkpnt )
stkpnt = stkpnt - 1
IF( endd-start <= select .AND. endd-start > 0 ) THEN

!        Do Insertion sort on D( START:ENDD )

  IF( dir == 0 ) THEN

!           Sort into decreasing order

    DO i = start + 1, endd
      DO j = i, start + 1, -1
        IF( d( j ) > d( j-1 ) ) THEN
          dmnmx = d( j )
          d( j ) = d( j-1 )
          d( j-1 ) = dmnmx
        ELSE
          CYCLE
        END IF
      END DO
    END DO

  ELSE

!           Sort into increasing order

    DO i = start + 1, endd
      DO j = i, start + 1, -1
        IF( d( j ) < d( j-1 ) ) THEN
          dmnmx = d( j )
          d( j ) = d( j-1 )
          d( j-1 ) = dmnmx
        ELSE
          CYCLE
        END IF
      END DO
    END DO

  END IF

ELSE IF( endd-start > select ) THEN

!        Partition D( START:ENDD ) and stack parts, largest one first

!        Choose partition entry as median of 3

  d1 = d( start )
  d2 = d( endd )
  i = ( start + endd ) / 2
  d3 = d( i )
  IF( d1 < d2 ) THEN
    IF( d3 < d1 ) THEN
      dmnmx = d1
    ELSE IF( d3 < d2 ) THEN
      dmnmx = d3
    ELSE
      dmnmx = d2
    END IF
  ELSE
    IF( d3 < d2 ) THEN
      dmnmx = d2
    ELSE IF( d3 < d1 ) THEN
      dmnmx = d3
    ELSE
      dmnmx = d1
    END IF
  END IF

  IF( dir == 0 ) THEN

!           Sort into decreasing order

    i = start - 1
    j = endd + 1

    60 j = j - 1
    IF( d( j ) < dmnmx ) GO TO 60

    80 i = i + 1
    IF( d( i ) > dmnmx ) GO TO 80
    IF( i < j ) THEN
      tmp = d( i )
      d( i ) = d( j )
      d( j ) = tmp
      GO TO 60
    END IF
    IF( j-start > endd-j-1 ) THEN
      stkpnt = stkpnt + 1
      stack( 1, stkpnt ) = start
      stack( 2, stkpnt ) = j
      stkpnt = stkpnt + 1
      stack( 1, stkpnt ) = j + 1
      stack( 2, stkpnt ) = endd
    ELSE
      stkpnt = stkpnt + 1
      stack( 1, stkpnt ) = j + 1
      stack( 2, stkpnt ) = endd
      stkpnt = stkpnt + 1
      stack( 1, stkpnt ) = start
      stack( 2, stkpnt ) = j
    END IF
  ELSE

!           Sort into increasing order

    i = start - 1
    j = endd + 1

    90 j = j - 1
    IF( d( j ) > dmnmx ) GO TO 90

    110 i = i + 1
    IF( d( i ) < dmnmx ) GO TO 110
    IF( i < j ) THEN
      tmp = d( i )
      d( i ) = d( j )
      d( j ) = tmp
      GO TO 90
    END IF
    IF( j-start > endd-j-1 ) THEN
      stkpnt = stkpnt + 1
      stack( 1, stkpnt ) = start
      stack( 2, stkpnt ) = j
      stkpnt = stkpnt + 1
      stack( 1, stkpnt ) = j + 1
      stack( 2, stkpnt ) = endd
    ELSE
      stkpnt = stkpnt + 1
      stack( 1, stkpnt ) = j + 1
      stack( 2, stkpnt ) = endd
      stkpnt = stkpnt + 1
      stack( 1, stkpnt ) = start
      stack( 2, stkpnt ) = j
    END IF
  END IF
END IF
IF( stkpnt > 0 ) GO TO 10

RETURN

!     End of DLASRT

END SUBROUTINE dlasrt


END MODULE LA_AUXMOD

MODULE dblas
! Contains only the INTEGER, REAL (wp) & COMPLEX (wp)
! parts of BLAS level 1.
! For single precision, change wp to point to sp instead of dp
! This very much simplified BLAS module is by Alan Miller
! alan @ vic.cmis.csiro.au    URL: www.ozemail.com.au/~milleraj

! Latest revision - 25 August 1998

USE la_precision, ONLY: wp => dp             ! From file l_auxmod.f90
IMPLICIT NONE

REAL (wp), PARAMETER :: zero = 0._wp, one = 1._wp

CONTAINS

! ============= dscal.f ==============
SUBROUTINE dscal(n, a, x, incx)

!     scales a vector by a constant.

INTEGER, INTENT(IN)       :: n, incx
REAL (wp), INTENT(IN)     :: a
REAL (wp), INTENT(IN OUT) :: x(:)

IF( n <= 0 .OR. incx <= 0 ) RETURN
IF(incx == 1) THEN
  x(:n) = a * x(:n)
  RETURN
END IF
x(:n*incx:incx) = a * x(:n*incx:incx)

RETURN
END SUBROUTINE dscal

! ============= ddot.f ==============
FUNCTION ddot(n, x, incx, y, incy) RESULT(fn_val)

!     forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.

INTEGER, INTENT(IN)   :: n, incx, incy
REAL (wp), INTENT(IN) :: x(:), y(:)
REAL (wp)             :: fn_val

fn_val = zero
IF(n <= 0) RETURN
IF(incx == 1 .AND. incy == 1) THEN
  fn_val = DOT_PRODUCT( x(:n), y(:n) )
  RETURN
END IF
fn_val = DOT_PRODUCT( x(:n*incx:incx), y(:n*incy:incy) )

RETURN
END FUNCTION ddot
! ============= daxpy.f ==============
SUBROUTINE daxpy(n, a, x, incx, y, incy)

!     constant times a vector plus a vector.
!     new y() = old y() + a.x()

INTEGER, INTENT(IN)       :: n, incx, incy
REAL (wp), INTENT(IN)     :: x(:), a
REAL (wp), INTENT(IN OUT) :: y(:)

IF (n <= 0) RETURN
IF (a == zero) RETURN
IF (incx == 1 .AND. incy == 1) THEN
  y(:n) = y(:n) + a * x(:n)
  RETURN
END IF
y(:n*incy:incy) = y(:n*incy:incy) + a * x(:n*incx:incx)

RETURN
END SUBROUTINE daxpy

! ============= idamax.f ==============
FUNCTION idamax(n, x, incx) RESULT(fn_val)

!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)

INTEGER, INTENT(IN)   :: n, incx
REAL (wp), INTENT(IN) :: x(:)
INTEGER               :: fn_val

INTEGER :: imax(1)

fn_val = 0
IF( n < 1 .OR. incx <= 0 ) RETURN
fn_val = 1
IF(n == 1) RETURN
IF(incx == 1) THEN
  imax = MAXLOC( ABS(x(:n)) )
ELSE
  imax = MAXLOC( ABS(x(:n*incx:incx)) )
END IF
fn_val = imax(1)

RETURN
END FUNCTION idamax
! ============= dswap.f ==============
SUBROUTINE dswap (n, x, incx, y, incy)

!     interchanges two vectors.

INTEGER, INTENT(IN)       :: n, incx, incy
REAL (wp), INTENT(IN OUT) :: x(:), y(:)

! Local variables
REAL (wp) :: temp(n)

IF(n <= 0) RETURN
IF(incx == 1 .AND. incy == 1) THEN
  temp = x(:n)
  x(:n) = y(:n)
  y(:n) = temp
  RETURN
END IF

temp = x(:n*incx:incx)
x(:n*incx:incx) = y(:n*incy:incy)
y(:n*incy:incy) = temp

RETURN
END SUBROUTINE dswap

! ============= dnrm2.f ==============
FUNCTION dnrm2 ( n, x, incx) RESULT( norm )
!     .. Scalar Arguments ..

INTEGER, INTENT(IN)   :: n
REAL (wp), INTENT(IN) :: x(:)
INTEGER, INTENT(IN)   :: incx
REAL (wp)             :: norm

!  DNRM2 returns the Euclidean norm of a vector via the function
!  name, so that

!     DNRM2 := sqrt( x'*x )


!  -- This version written on 25-October-1982.
!     Modified on 14-October-1993 to inline the call to DLASSQ.
!     Sven Hammarling, Nag Ltd.

!     .. Local Scalars ..
INTEGER   :: ix
REAL (wp) :: absxi, scale, ssq
!     ..
!     .. Executable Statements ..
IF( n < 1 .OR. incx < 1 ) THEN
  norm  = zero
ELSE IF( n == 1 ) THEN
  norm  = ABS( x( 1 ) )
ELSE
  scale = zero
  ssq   = one
!        The following loop is equivalent to this call to the LAPACK
!        auxiliary routine:
!        CALL DLASSQ( N, X, INCX, SCALE, SSQ )

  DO  ix = 1, 1 + ( n - 1 )*incx, incx
    IF( x( ix ) /= zero ) THEN
      absxi = ABS( x( ix ) )
      IF( scale < absxi ) THEN
        ssq   = one   + ssq*( scale/absxi )**2
        scale = absxi
      ELSE
        ssq   = ssq   +     ( absxi/scale )**2
      END IF
    END IF
  END DO
  norm  = scale * SQRT( ssq )
END IF

RETURN

!     End of DNRM2.

END FUNCTION dnrm2



END MODULE dblas

MODULE dblas2
! For single precision, change wp to point to sp instead of dp
! This very much simplified BLAS2 module is by Alan Miller
! alan @ vic.cmis.csiro.au    URL: www.ozemail.com.au/~milleraj

! Latest revision - 23 October 1997

USE la_precision, ONLY: wp => dp             ! From file l_auxmod.f90
USE la_auxmod
USE dblas
IMPLICIT NONE

CONTAINS


SUBROUTINE dger  ( m, n, alpha, x, incx, y, incy, a, lda )
! ELF90 translation by Alan Miller   31-Aug-1997

!     .. Scalar Arguments ..
REAL (wp), INTENT(IN)     :: alpha
INTEGER, INTENT(IN)       :: incx, incy, lda, m, n
!     .. Array Arguments ..
REAL (wp), INTENT(IN)     :: x(:), y(:)
REAL (wp), INTENT(IN OUT) :: a(:,:)
!     ..

!  Purpose
!  =======

!  DGER   performs the rank 1 operation

!     A := alpha*x*y' + A,

!  where alpha is a scalar, x is an m element vector, y is an n element
!  vector and A is an m by n matrix.

!  Parameters
!  ==========

!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.

!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.

!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.

!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( m - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the m
!           element vector x.
!           Unchanged on exit.

!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.

!  Y      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.

!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.

!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients. On exit, A is
!           overwritten by the updated matrix.

!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, m ).
!           Unchanged on exit.


!  Level 2 Blas routine.

!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.


!     .. Local Scalars ..
REAL (wp) :: temp
INTEGER   ::  i, info, ix, j, jy, kx
!     .. External Subroutines ..
! EXTERNAL           erinfo
!     .. Intrinsic Functions ..
! INTRINSIC          MAX
!     ..
!     .. Executable Statements ..

!     Test the input parameters.

info = 0
IF     ( m < 0 ) THEN
  info = 1
ELSE IF( n < 0 ) THEN
  info = 2
ELSE IF( incx == 0 ) THEN
  info = 5
ELSE IF( incy == 0 ) THEN
  info = 7
ELSE IF( lda < MAX( 1, m ) ) THEN
  info = 9
END IF
IF( info /= 0 ) THEN
  CALL erinfo(info,  'DGER  ')
  RETURN
END IF

!     Quick return if possible.

IF( ( m == 0 ).OR.( n == 0 ).OR.( alpha == zero ) ) RETURN

!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.

IF( incy > 0 ) THEN
  jy = 1
ELSE
  jy = 1 - ( n - 1 )*incy
END IF
IF( incx == 1 ) THEN
  DO j = 1, n
    IF( y( jy ) /= zero ) THEN
      temp = alpha*y( jy )
      a( 1:m, j ) = a( 1:m, j ) + x( 1:m )*temp
    END IF
    jy = jy + incy
  END DO
ELSE
  IF( incx > 0 ) THEN
    kx = 1
  ELSE
    kx = 1 - ( m - 1 )*incx
  END IF
  DO j = 1, n
    IF( y( jy ) /= zero ) THEN
      temp = alpha*y( jy )
      ix   = kx
      DO i = 1, m
        a( i, j ) = a( i, j ) + x( ix )*temp
        ix        = ix        + incx
      END DO
    END IF
    jy = jy + incy
  END DO
END IF

RETURN

!     End of DGER  .

END SUBROUTINE dger


SUBROUTINE dsymv ( uplo, n, alpha, a, lda, x, incx, beta, y, incy )
! ELF90 translation by Alan Miller   31-Aug-1997

!     .. Scalar Arguments ..
REAL (wp), INTENT(IN)         :: alpha, beta
INTEGER, INTENT(IN)           :: incx, incy, lda, n
CHARACTER (LEN=1), INTENT(IN) :: uplo
!     .. Array Arguments ..
REAL (wp), INTENT(IN)         :: a(:,:), x(:)
REAL (wp), INTENT(IN OUT)     :: y(:)
!     ..

!  Purpose
!  =======

!  DSYMV  performs the matrix-vector  operation

!     y := alpha*A*x + beta*y,

!  where alpha and beta are scalars, x and y are n element vectors and
!  A is an n by n symmetric matrix.

!  Parameters
!  ==========

!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the array A is to be referenced as
!           follows:

!              UPLO = 'U' or 'u'   Only the upper triangular part of A
!                                  is to be referenced.

!              UPLO = 'L' or 'l'   Only the lower triangular part of A
!                                  is to be referenced.

!           Unchanged on exit.

!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.

!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.

!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular part of the symmetric matrix and the strictly
!           lower triangular part of A is not referenced.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular part of the symmetric matrix and the strictly
!           upper triangular part of A is not referenced.
!           Unchanged on exit.

!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.

!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.

!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.

!  BETA   - DOUBLE PRECISION.
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.

!  Y      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y. On exit, Y is overwritten by the updated
!           vector y.

!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.


!  Level 2 Blas routine.

!  -- Written on 20-July-1986.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.

!     .. Local Scalars ..
REAL (wp) :: temp1, temp2
INTEGER   ::  i, info, ix, iy, j, jx, jy, kx, ky
!     .. External Functions ..
! LOGICAL ::  lsame
! EXTERNAL           lsame
!     .. External Subroutines ..
! EXTERNAL           erinfo
!     .. Intrinsic Functions ..
! INTRINSIC          MAX
!     ..
!     .. Executable Statements ..

!     Test the input parameters.

info = 0
IF      ( .NOT.lsame( uplo, 'U' ).AND. .NOT.lsame( uplo, 'L' )      ) THEN
  info = 1
ELSE IF ( n < 0 ) THEN
  info = 2
ELSE IF ( lda < MAX(1,n) ) THEN
  info = 5
ELSE IF ( incx == 0 ) THEN
  info = 7
ELSE IF ( incy == 0 ) THEN
  info = 10
END IF
IF( info /= 0 ) THEN
  CALL erinfo(info, 'DSYMV ')
  RETURN
END IF

!     Quick return if possible.

IF( ( n == 0 ).OR.( ( alpha == zero ).AND.( beta == one ) ) ) RETURN

!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.

!     First form  y := beta*y  and set up the start points in X and Y if
!     the increments are not both unity.

IF( ( incx == 1 ).AND.( incy == 1 ) ) THEN
  IF( beta /= one ) THEN
    IF( beta == zero ) THEN
      DO i = 1, n
        y( i ) = zero
      END DO
    ELSE
      DO i = 1, n
        y( i ) = beta*y( i )
      END DO
    END IF
  END IF
ELSE
  IF( incx > 0 ) THEN
    kx = 1
  ELSE
    kx = 1 - ( n - 1 )*incx
  END IF
  IF( incy > 0 ) THEN
    ky = 1
  ELSE
    ky = 1 - ( n - 1 )*incy
  END IF
  IF( beta /= one ) THEN
    iy = ky
    IF( beta == zero ) THEN
      DO i = 1, n
        y( iy ) = zero
        iy      = iy   + incy
      END DO
    ELSE
      DO i = 1, n
        y( iy ) = beta*y( iy )
        iy      = iy           + incy
      END DO
    END IF
  END IF
END IF
IF( alpha == zero ) RETURN
IF( lsame( uplo, 'U' ) ) THEN

!        Form  y  when A is stored in upper triangle.

  IF( ( incx == 1 ).AND.( incy == 1 ) ) THEN
    DO j = 1, n
      temp1 = alpha*x( j )
      temp2 = zero
      DO i = 1, j - 1
        y( i ) = y( i ) + temp1*a( i, j )
        temp2  = temp2  + a( i, j )*x( i )
      END DO
      y( j ) = y( j ) + temp1*a( j, j ) + alpha*temp2
    END DO
  ELSE
    ix = kx - incx
    DO j = 1, n
      temp1 = alpha*x( ix + incx )
      temp2 = zero
      ix    = kx
      iy    = ky
      DO i = 1, j - 1
        y( iy ) = y( iy ) + temp1*a( i, j )
        temp2   = temp2   + a( i, j )*x( ix )
        ix      = ix      + incx
        iy      = iy      + incy
      END DO
      y( iy ) = y( iy ) + temp1*a( j, j ) + alpha*temp2
    END DO
  END IF
ELSE

!        Form  y  when A is stored in lower triangle.

  IF( ( incx == 1 ).AND.( incy == 1 ) ) THEN
    DO j = 1, n
      temp1  = alpha*x( j )
      temp2  = zero
      y( j ) = y( j )       + temp1*a( j, j )
      DO i = j + 1, n
        y( i ) = y( i ) + temp1*a( i, j )
        temp2  = temp2  + a( i, j )*x( i )
      END DO
      y( j ) = y( j ) + alpha*temp2
    END DO
  ELSE
    jx = kx
    jy = ky
    DO j = 1, n
      temp1   = alpha*x( jx )
      temp2   = zero
      y( jy ) = y( jy )       + temp1*a( j, j )
      ix      = jx
      iy      = jy
      DO i = j + 1, n
        ix      = ix      + incx
        iy      = iy      + incy
        y( iy ) = y( iy ) + temp1*a( i, j )
        temp2   = temp2   + a( i, j )*x( ix )
      END DO
      y( jy ) = y( jy ) + alpha*temp2
      jx      = jx      + incx
      jy      = jy      + incy
    END DO
  END IF
END IF

RETURN

!     End of DSYMV .

END SUBROUTINE dsymv


SUBROUTINE dsyr2 ( uplo, n, alpha, x, incx, y, incy, a, lda )
! ELF90 translation by Alan Miller   31-Aug-1997

!     .. Scalar Arguments ..
REAL (wp), INTENT(IN)         :: alpha
INTEGER, INTENT(IN)           :: incx, incy, lda, n
CHARACTER (LEN=1), INTENT(IN) :: uplo
!     .. Array Arguments ..
REAL (wp), INTENT(IN)         :: x(:), y(:)
REAL (wp), INTENT(IN OUT)     :: a(:,:)
!     ..

!  Purpose
!  =======

!  DSYR2  performs the symmetric rank 2 operation

!     A := alpha*x*y' + alpha*y*x' + A,

!  where alpha is a scalar, x and y are n element vectors and A is an n
!  by n symmetric matrix.

!  Parameters
!  ==========

!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the upper or lower
!           triangular part of the array A is to be referenced as
!           follows:

!              UPLO = 'U' or 'u'   Only the upper triangular part of A
!                                  is to be referenced.

!              UPLO = 'L' or 'l'   Only the lower triangular part of A
!                                  is to be referenced.

!           Unchanged on exit.

!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.

!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.

!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x.
!           Unchanged on exit.

!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.

!  Y      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCY ) ).
!           Before entry, the incremented array Y must contain the n
!           element vector y.
!           Unchanged on exit.

!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.

!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular part of the symmetric matrix and the strictly
!           lower triangular part of A is not referenced. On exit, the
!           upper triangular part of the array A is overwritten by the
!           upper triangular part of the updated matrix.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular part of the symmetric matrix and the strictly
!           upper triangular part of A is not referenced. On exit, the
!           lower triangular part of the array A is overwritten by the
!           lower triangular part of the updated matrix.

!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.


!  Level 2 Blas routine.

!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.


!     .. Local Scalars ..
REAL (wp) :: temp1, temp2
INTEGER ::  i, info, ix, iy, j, jx, jy, kx, ky
!     .. External Functions ..
! LOGICAL ::  lsame
! EXTERNAL           lsame
!     .. External Subroutines ..
! EXTERNAL           erinfo
!     .. Intrinsic Functions ..
! INTRINSIC          MAX
!     ..
!     .. Executable Statements ..

!     Test the input parameters.

info = 0
IF     ( .NOT.lsame( uplo, 'U' ).AND. .NOT.lsame( uplo, 'L' )      ) THEN
  info = 1
ELSE IF( n < 0 ) THEN
  info = 2
ELSE IF( incx == 0 ) THEN
  info = 5
ELSE IF( incy == 0 ) THEN
  info = 7
ELSE IF( lda < MAX( 1, n ) ) THEN
  info = 9
END IF
IF( info /= 0 ) THEN
  CALL erinfo(info, 'DSYR2 ')
  RETURN
END IF

!     Quick return if possible.

IF( ( n == 0 ).OR.( alpha == zero ) ) RETURN

!     Set up the start points in X and Y if the increments are not both
!     unity.

IF( ( incx /= 1 ).OR.( incy /= 1 ) ) THEN
  IF( incx > 0 ) THEN
    kx = 1
  ELSE
    kx = 1 - ( n - 1 )*incx
  END IF
  IF( incy > 0 ) THEN
    ky = 1
  ELSE
    ky = 1 - ( n - 1 )*incy
  END IF
  jx = kx
  jy = ky
END IF

!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.

IF( lsame( uplo, 'U' ) ) THEN

!        Form  A  when A is stored in the upper triangle.

  IF( ( incx == 1 ).AND.( incy == 1 ) ) THEN
    DO j = 1, n
      IF( ( x( j ) /= zero ).OR.( y( j ) /= zero ) ) THEN
        temp1 = alpha*y( j )
        temp2 = alpha*x( j )
        a( 1:j, j ) = a( 1:j, j ) + x( 1:j )*temp1 + y( 1:j )*temp2
      END IF
    END DO
  ELSE
    DO j = 1, n
      IF( ( x( jx ) /= zero ).OR.( y( jy ) /= zero ) ) THEN
        temp1 = alpha*y( jy )
        temp2 = alpha*x( jx )
        ix    = kx
        iy    = ky
        DO i = 1, j
          a( i, j ) = a( i, j ) + x( ix )*temp1+ y( iy )*temp2
          ix        = ix        + incx
          iy        = iy        + incy
        END DO
      END IF
      jx = jx + incx
      jy = jy + incy
    END DO
  END IF
ELSE

!        Form  A  when A is stored in the lower triangle.

  IF( ( incx == 1 ).AND.( incy == 1 ) ) THEN
    DO j = 1, n
      IF( ( x( j ) /= zero ).OR.( y( j ) /= zero ) ) THEN
        temp1 = alpha*y( j )
        temp2 = alpha*x( j )
        a( j:n, j ) = a( j:n, j ) + x( j:n )*temp1 + y( j:n )*temp2
      END IF
    END DO
  ELSE
    DO j = 1, n
      IF( ( x( jx ) /= zero ).OR.( y( jy ) /= zero ) ) THEN
        temp1 = alpha*y( jy )
        temp2 = alpha*x( jx )
        ix    = jx
        iy    = jy
        DO i = j, n
          a( i, j ) = a( i, j ) + x( ix )*temp1+ y( iy )*temp2
          ix        = ix        + incx
          iy        = iy        + incy
        END DO
      END IF
      jx = jx + incx
      jy = jy + incy
    END DO
  END IF
END IF

RETURN

!     End of DSYR2 .

END SUBROUTINE dsyr2


SUBROUTINE dgemv ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )
! ELF90 translation by Alan Miller   31-Aug-1997

!     .. Scalar Arguments ..
REAL (wp), INTENT(IN)         :: alpha, beta
INTEGER, INTENT(IN)           :: incx, incy, lda, m, n
CHARACTER (LEN=1), INTENT(IN) :: trans
!     .. Array Arguments ..
REAL (wp), INTENT(IN)         :: a(:,:), x(:)
REAL (wp), INTENT(IN OUT)     :: y(:)
!     ..

!  Purpose
!  =======

!  DGEMV  performs one of the matrix-vector operations

!     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,

!  where alpha and beta are scalars, x and y are vectors and A is an
!  m by n matrix.

!  Parameters
!  ==========

!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:

!              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.

!              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.

!              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.

!           Unchanged on exit.

!  M      - INTEGER.
!           On entry, M specifies the number of rows of the matrix A.
!           M must be at least zero.
!           Unchanged on exit.

!  N      - INTEGER.
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.

!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.

!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry, the leading m by n part of the array A must
!           contain the matrix of coefficients.
!           Unchanged on exit.

!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least max( 1, m ).
!           Unchanged on exit.

!  X      - DOUBLE PRECISION array of DIMENSION at least
!           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
!           Before entry, the incremented array X must contain the
!           vector x.
!           Unchanged on exit.

!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.

!  BETA   - DOUBLE PRECISION.
!           On entry, BETA specifies the scalar beta. When BETA is
!           supplied as zero then Y need not be set on input.
!           Unchanged on exit.

!  Y      - DOUBLE PRECISION array of DIMENSION at least
!           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!           and at least
!           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
!           Before entry with BETA non-zero, the incremented array Y
!           must contain the vector y. On exit, Y is overwritten by the
!           updated vector y.

!  INCY   - INTEGER.
!           On entry, INCY specifies the increment for the elements of
!           Y. INCY must not be zero.
!           Unchanged on exit.


!  Level 2 Blas routine.

!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.


!     .. Local Scalars ..
REAL (wp) :: temp
INTEGER   ::  i, info, ix, iy, j, jx, jy, kx, ky, lenx, leny
!     .. External Functions ..
! LOGICAL ::  lsame
! EXTERNAL           lsame
!     .. External Subroutines ..
! EXTERNAL           erinfo
!     .. Intrinsic Functions ..
! INTRINSIC          MAX
!     ..
!     .. Executable Statements ..

!     Test the input parameters.

info = 0
IF     ( .NOT.lsame( trans, 'N' ).AND. .NOT.lsame( trans, 'T' ).AND.  &
  .NOT.lsame( trans, 'C' )      ) THEN
  info = 1
ELSE IF( m < 0 ) THEN
  info = 2
ELSE IF( n < 0 ) THEN
  info = 3
ELSE IF( lda < MAX( 1, m ) ) THEN
  info = 6
ELSE IF( incx == 0 ) THEN
  info = 8
ELSE IF( incy == 0 ) THEN
  info = 11
END IF
IF( info /= 0 ) THEN
  CALL erinfo(info,  'DGEMV ')
  RETURN
END IF

!     Quick return if possible.

IF( ( m == 0 ).OR.( n == 0 ).OR.  &
( ( alpha == zero ).AND.( beta == one ) ) )RETURN

!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.

IF( lsame( trans, 'N' ) ) THEN
  lenx = n
  leny = m
ELSE
  lenx = m
  leny = n
END IF
IF( incx > 0 ) THEN
  kx = 1
ELSE
  kx = 1 - ( lenx - 1 )*incx
END IF
IF( incy > 0 ) THEN
  ky = 1
ELSE
  ky = 1 - ( leny - 1 )*incy
END IF

!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.

!     First form  y := beta*y.

IF( beta /= one ) THEN
  IF( incy == 1 ) THEN
    IF( beta == zero ) THEN
      y( :leny ) = zero
    ELSE
      y( :leny ) = beta*y( :leny )
    END IF
  ELSE
    iy = ky
    IF( beta == zero ) THEN
      DO i = 1, leny
        y( iy ) = zero
        iy      = iy   + incy
      END DO
    ELSE
      DO i = 1, leny
        y( iy ) = beta*y( iy )
        iy      = iy           + incy
      END DO
    END IF
  END IF
END IF
IF( alpha == zero ) RETURN
IF( lsame( trans, 'N' ) ) THEN

!        Form  y := alpha*A*x + y.

  jx = kx
  IF( incy == 1 ) THEN
    DO j = 1, n
      IF( x( jx ) /= zero ) THEN
        temp = alpha*x( jx )
        y( 1:m ) = y( 1:m ) + temp*a( 1:m, j )
      END IF
      jx = jx + incx
    END DO
  ELSE
    DO j = 1, n
      IF( x( jx ) /= zero ) THEN
        temp = alpha*x( jx )
        iy   = ky
        DO i = 1, m
          y( iy ) = y( iy ) + temp*a( i, j )
          iy      = iy      + incy
        END DO
      END IF
      jx = jx + incx
    END DO
  END IF
ELSE

!        Form  y := alpha*A'*x + y.

  jy = ky
  IF( incx == 1 ) THEN
    DO j = 1, n
      temp = DOT_PRODUCT( a(1:m,j), x(1:m) )
      y( jy ) = y( jy ) + alpha*temp
      jy      = jy      + incy
    END DO
  ELSE
    DO j = 1, n
      temp = zero
      ix   = kx
      DO i = 1, m
        temp = temp + a( i, j )*x( ix )
        ix   = ix   + incx
      END DO
      y( jy ) = y( jy ) + alpha*temp
      jy      = jy      + incy
    END DO
  END IF
END IF

RETURN

!     End of DGEMV .

END SUBROUTINE dgemv


SUBROUTINE dtrmv ( uplo, trans, diag, n, a, lda, x, incx )
! ELF90 translation by Alan Miller   31-Aug-1997

!     .. Scalar Arguments ..
INTEGER, INTENT(IN)           :: incx, lda, n
CHARACTER (LEN=1), INTENT(IN) :: diag, trans, uplo
!     .. Array Arguments ..
REAL (wp), INTENT(IN)         :: a(:,:)
REAL (wp), INTENT(IN OUT)     :: x(:)
!     ..

!  Purpose
!  =======

!  DTRMV  performs one of the matrix-vector operations

!     x := A*x,   or   x := A'*x,

!  where x is n element vector and A is an n by n unit, or non-unit,
!  upper or lower triangular matrix.

!  Parameters
!  ==========

!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:

!              UPLO = 'U' or 'u'   A is an upper triangular matrix.

!              UPLO = 'L' or 'l'   A is a lower triangular matrix.

!           Unchanged on exit.

!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the operation to be performed as
!           follows:

!              TRANS = 'N' or 'n'   x := A*x.

!              TRANS = 'T' or 't'   x := A'*x.

!              TRANS = 'C' or 'c'   x := A'*x.

!           Unchanged on exit.

!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:

!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.

!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.

!           Unchanged on exit.

!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.

!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular matrix and the strictly lower triangular part of
!           A is not referenced.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular matrix and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!           A are not referenced either, but are assumed to be unity.
!           Unchanged on exit.

!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.

!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element vector x. On exit, X is overwritten with the
!           tranformed vector x.

!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.


!  Level 2 Blas routine.

!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.


!     .. Local Scalars ..
REAL (wp) :: temp
INTEGER ::  i, info, ix, j, jx, kx
LOGICAL ::  nounit
!     .. External Functions ..
! LOGICAL ::  lsame
! EXTERNAL           lsame
!     .. External Subroutines ..
! EXTERNAL           erinfo
!     .. Intrinsic Functions ..
! INTRINSIC          MAX
!     ..
!     .. Executable Statements ..

!     Test the input parameters.

info = 0
IF     ( .NOT.lsame( uplo , 'U' ).AND. .NOT.lsame( uplo , 'L' )      ) THEN
  info = 1
ELSE IF( .NOT.lsame( trans, 'N' ).AND.  &
  .NOT.lsame( trans, 'T' ).AND..NOT.lsame( trans, 'C' )      ) THEN
  info = 2
ELSE IF( .NOT.lsame( diag , 'U' ).AND.  &
  .NOT.lsame( diag , 'N' )      ) THEN
  info = 3
ELSE IF( n < 0 ) THEN
  info = 4
ELSE IF( lda < MAX( 1, n ) ) THEN
  info = 6
ELSE IF( incx == 0 ) THEN
  info = 8
END IF
IF( info /= 0 ) THEN
  CALL erinfo(info, 'DTRMV ')
  RETURN
END IF

!     Quick return if possible.

IF( n == 0 ) RETURN

nounit = lsame( diag, 'N' )

!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.

IF( incx <= 0 ) THEN
  kx = 1 - ( n - 1 )*incx
ELSE IF( incx /= 1 ) THEN
  kx = 1
END IF

!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.

IF( lsame( trans, 'N' ) ) THEN

!        Form  x := A*x.

  IF( lsame( uplo, 'U' ) ) THEN
    IF( incx == 1 ) THEN
      DO j = 1, n
        IF( x( j ) /= zero ) THEN
          temp = x( j )
          x( 1:j-1 ) = x( 1:j-1 ) + temp*a( 1:j-1, j )
          IF( nounit ) x( j ) = x( j )*a( j, j )
        END IF
      END DO
    ELSE
      jx = kx
      DO j = 1, n
        IF( x( jx ) /= zero ) THEN
          temp = x( jx )
          ix   = kx
          DO i = 1, j - 1
            x( ix ) = x( ix ) + temp*a( i, j )
            ix      = ix      + incx
          END DO
          IF( nounit ) x( jx ) = x( jx )*a( j, j )
        END IF
        jx = jx + incx
      END DO
    END IF
  ELSE
    IF( incx == 1 ) THEN
      DO j = n, 1, -1
        IF( x( j ) /= zero ) THEN
          temp = x( j )
          x( n:j+1:-1 ) = x( n:j+1:-1 ) + temp*a( n:j+1:-1, j )
          IF( nounit ) x( j ) = x( j )*a( j, j )
        END IF
      END DO
    ELSE
      kx = kx + ( n - 1 )*incx
      jx = kx
      DO j = n, 1, -1
        IF( x( jx ) /= zero ) THEN
          temp = x( jx )
          ix   = kx
          DO i = n, j + 1, -1
            x( ix ) = x( ix ) + temp*a( i, j )
            ix      = ix      - incx
          END DO
          IF( nounit ) x( jx ) = x( jx )*a( j, j )
        END IF
        jx = jx - incx
      END DO
    END IF
  END IF
ELSE

!        Form  x := A'*x.

  IF( lsame( uplo, 'U' ) ) THEN
    IF( incx == 1 ) THEN
      DO j = n, 1, -1
        temp = x( j )
        IF( nounit ) temp = temp*a( j, j )
        temp = temp + DOT_PRODUCT( a(j-1:1:-1,j), x(j-1:1:-1) )
        x( j ) = temp
      END DO
    ELSE
      jx = kx + ( n - 1 )*incx
      DO j = n, 1, -1
        temp = x( jx )
        ix   = jx
        IF( nounit ) temp = temp*a( j, j )
        DO i = j - 1, 1, -1
          ix   = ix   - incx
          temp = temp + a( i, j )*x( ix )
        END DO
        x( jx ) = temp
        jx      = jx   - incx
      END DO
    END IF
  ELSE
    IF( incx == 1 ) THEN
      DO j = 1, n
        temp = x( j )
        IF( nounit ) temp = temp*a( j, j )
        DO i = j + 1, n
          temp = temp + a( i, j )*x( i )
        END DO
        x( j ) = temp
      END DO
    ELSE
      jx = kx
      DO j = 1, n
        temp = x( jx )
        ix   = jx
        IF( nounit ) temp = temp*a( j, j )
        DO i = j + 1, n
          ix   = ix   + incx
          temp = temp + a( i, j )*x( ix )
        END DO
        x( jx ) = temp
        jx      = jx   + incx
      END DO
    END IF
  END IF
END IF

RETURN

!     End of DTRMV .

END SUBROUTINE dtrmv

SUBROUTINE dtrsv ( uplo, trans, diag, n, a, lda, x, incx )
! ELF90 translation by Alan Miller   31-Aug-1997

!     .. Scalar Arguments ..
INTEGER, INTENT(IN)           :: incx, lda, n
CHARACTER (LEN=1), INTENT(IN) :: diag, trans, uplo
!     .. Array Arguments ..
REAL (wp), INTENT(IN)         :: a(:,:)
REAL (wp), INTENT(IN OUT)     :: x(:)
!     ..

!  Purpose
!  =======

!  DTRSV  solves one of the systems of equations

!     A*x = b,   or   A'*x = b,

!  where b and x are n element vectors and A is an n by n unit, or
!  non-unit, upper or lower triangular matrix.

!  No test for singularity or near-singularity is included in this
!  routine. Such tests must be performed before calling this routine.

!  Parameters
!  ==========

!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix is an upper or
!           lower triangular matrix as follows:

!              UPLO = 'U' or 'u'   A is an upper triangular matrix.

!              UPLO = 'L' or 'l'   A is a lower triangular matrix.

!           Unchanged on exit.

!  TRANS  - CHARACTER*1.
!           On entry, TRANS specifies the equations to be solved as
!           follows:

!              TRANS = 'N' or 'n'   A*x = b.

!              TRANS = 'T' or 't'   A'*x = b.

!              TRANS = 'C' or 'c'   A'*x = b.

!           Unchanged on exit.

!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit
!           triangular as follows:

!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.

!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.

!           Unchanged on exit.

!  N      - INTEGER.
!           On entry, N specifies the order of the matrix A.
!           N must be at least zero.
!           Unchanged on exit.

!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
!           Before entry with  UPLO = 'U' or 'u', the leading n by n
!           upper triangular part of the array A must contain the upper
!           triangular matrix and the strictly lower triangular part of
!           A is not referenced.
!           Before entry with UPLO = 'L' or 'l', the leading n by n
!           lower triangular part of the array A must contain the lower
!           triangular matrix and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!           A are not referenced either, but are assumed to be unity.
!           Unchanged on exit.

!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least
!           max( 1, n ).
!           Unchanged on exit.

!  X      - DOUBLE PRECISION array of dimension at least
!           ( 1 + ( n - 1 )*abs( INCX ) ).
!           Before entry, the incremented array X must contain the n
!           element right-hand side vector b. On exit, X is overwritten
!           with the solution vector x.

!  INCX   - INTEGER.
!           On entry, INCX specifies the increment for the elements of
!           X. INCX must not be zero.
!           Unchanged on exit.


!  Level 2 Blas routine.

!  -- Written on 22-October-1986.
!     Jack Dongarra, Argonne National Lab.
!     Jeremy Du Croz, Nag Central Office.
!     Sven Hammarling, Nag Central Office.
!     Richard Hanson, Sandia National Labs.


!     .. Local Scalars ..
REAL (wp) :: temp
INTEGER   :: i, info, ix, j, jx, kx
LOGICAL   :: nounit
!     .. External Functions ..
! LOGICAL :: lsame
! EXTERNAL           lsame
!     .. External Subroutines ..
! EXTERNAL           erinfo
!     .. Intrinsic Functions ..
! INTRINSIC          MAX
!     ..
!     .. Executable Statements ..

!     Test the input parameters.

info = 0
IF     ( .NOT.lsame( uplo , 'U' ).AND. .NOT.lsame( uplo , 'L' )      ) THEN
  info = 1
ELSE IF( .NOT.lsame( trans, 'N' ).AND.  &
  .NOT.lsame( trans, 'T' ).AND..NOT.lsame( trans, 'C' )      ) THEN
  info = 2
ELSE IF( .NOT.lsame( diag , 'U' ).AND.  &
  .NOT.lsame( diag , 'N' )      ) THEN
  info = 3
ELSE IF( n < 0 ) THEN
  info = 4
ELSE IF( lda < MAX( 1, n ) ) THEN
  info = 6
ELSE IF( incx == 0 ) THEN
  info = 8
END IF
IF( info /= 0 ) THEN
  CALL erinfo(info, 'DTRSV ')
  RETURN
END IF

!     Quick return if possible.

IF( n == 0 ) RETURN

nounit = lsame( diag, 'N' )

!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.

IF( incx <= 0 ) THEN
  kx = 1 - ( n - 1 )*incx
ELSE IF( incx /= 1 ) THEN
  kx = 1
END IF

!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.

IF( lsame( trans, 'N' ) ) THEN

!        Form  x := inv( A )*x.

  IF( lsame( uplo, 'U' ) ) THEN
    IF( incx == 1 ) THEN
      DO j = n, 1, -1
        IF( x( j ) /= zero ) THEN
          IF( nounit ) x( j ) = x( j )/a( j, j )
          temp = x( j )
          x( j-1:1:-1 ) = x( j-1:1:-1 ) - temp*a( j-1:1:-1, j )
        END IF
      END DO
    ELSE
      jx = kx + ( n - 1 )*incx
      DO j = n, 1, -1
        IF( x( jx ) /= zero ) THEN
          IF( nounit ) x( jx ) = x( jx )/a( j, j )
          temp = x( jx )
          ix   = jx
          DO i = j - 1, 1, -1
            ix      = ix      - incx
            x( ix ) = x( ix ) - temp*a( i, j )
          END DO
        END IF
        jx = jx - incx
      END DO
    END IF
  ELSE
    IF( incx == 1 ) THEN
      DO j = 1, n
        IF( x( j ) /= zero ) THEN
          IF( nounit ) x( j ) = x( j )/a( j, j )
          temp = x( j )
          x( j+1:n ) = x( j+1:n ) - temp*a( j+1:n, j )
        END IF
      END DO
    ELSE
      jx = kx
      DO j = 1, n
        IF( x( jx ) /= zero ) THEN
          IF( nounit ) x( jx ) = x( jx )/a( j, j )
          temp = x( jx )
          ix   = jx
          DO i = j + 1, n
            ix      = ix      + incx
            x( ix ) = x( ix ) - temp*a( i, j )
          END DO
        END IF
        jx = jx + incx
      END DO
    END IF
  END IF
ELSE

!        Form  x := inv( A' )*x.

  IF( lsame( uplo, 'U' ) ) THEN
    IF( incx == 1 ) THEN
      DO j = 1, n
        temp = x( j ) - DOT_PRODUCT( a( 1:j-1, j ), x( 1:j-1 ) )
        IF( nounit ) temp = temp/a( j, j )
        x( j ) = temp
      END DO
    ELSE
      jx = kx
      DO j = 1, n
        temp = x( jx )
        ix   = kx
        DO i = 1, j - 1
          temp = temp - a( i, j )*x( ix )
          ix   = ix   + incx
        END DO
        IF( nounit ) temp = temp/a( j, j )
        x( jx ) = temp
        jx      = jx   + incx
      END DO
    END IF
  ELSE
    IF( incx == 1 ) THEN
      DO j = n, 1, -1
        temp = x( j ) - DOT_PRODUCT( a( n:j+1:-1, j), x( n:j+1:-1 ) )
        IF( nounit ) temp = temp/a( j, j )
        x( j ) = temp
      END DO
    ELSE
      kx = kx + ( n - 1 )*incx
      jx = kx
      DO j = n, 1, -1
        temp = x( jx )
        ix   = kx
        DO i = n, j + 1, -1
          temp = temp - a( i, j )*x( ix )
          ix   = ix   - incx
        END DO
        IF( nounit ) temp = temp/a( j, j )
        x( jx ) = temp
        jx      = jx   - incx
      END DO
    END IF
  END IF

END IF

RETURN

!     End of DTRSV .

END SUBROUTINE dtrsv

END MODULE dblas2


MODULE dblas3
! For single precision, change wp to point to sp instead of dp
! This very much simplified BLAS2 module is by Alan Miller
! alan @ vic.cmis.csiro.au    URL: www.ozemail.com.au/~milleraj

! Latest revision - 3 September 1997

USE la_precision, ONLY: wp => dp             ! From file l_auxmod.f90
USE la_auxmod
USE dblas
IMPLICIT NONE

CONTAINS




SUBROUTINE dsyr2k( uplo, trans, n, k, alpha, a, lda, b, ldb,beta, c, ldc )
! ELF90 translation by Alan Miller   31-Aug-1997

!     .. Scalar Arguments ..
CHARACTER (LEN=1), INTENT(IN) :: uplo, trans
INTEGER, INTENT(IN)           :: n, k, lda, ldb, ldc
REAL (wp), INTENT(IN)         :: alpha, beta
!     .. Array Arguments ..
REAL (wp), INTENT(IN)         :: a(:,:), b(:,:)
REAL (wp), INTENT(IN OUT)     :: c(:,:)
!     ..

!  Purpose
!  =======

!  DSYR2K  performs one of the symmetric rank 2k operations

!     C := alpha*A*B' + alpha*B*A' + beta*C,

!  or

!     C := alpha*A'*B + alpha*B'*A + beta*C,

!  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
!  and  A and B  are  n by k  matrices  in the  first  case  and  k by n
!  matrices in the second case.

!  Parameters
!  ==========

!  UPLO   - CHARACTER*1.
!           On  entry,   UPLO  specifies  whether  the  upper  or  lower
!           triangular  part  of the  array  C  is to be  referenced  as
!           follows:

!              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
!                                  is to be referenced.

!              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
!                                  is to be referenced.

!           Unchanged on exit.

!  TRANS  - CHARACTER*1.
!           On entry,  TRANS  specifies the operation to be performed as
!           follows:

!              TRANS = 'N' or 'n'   C := alpha*A*B' + alpha*B*A' +
!                                        beta*C.

!              TRANS = 'T' or 't'   C := alpha*A'*B + alpha*B'*A +
!                                        beta*C.

!              TRANS = 'C' or 'c'   C := alpha*A'*B + alpha*B'*A +
!                                        beta*C.

!           Unchanged on exit.

!  N      - INTEGER.
!           On entry,  N specifies the order of the matrix C.  N must be
!           at least zero.
!           Unchanged on exit.

!  K      - INTEGER.
!           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
!           of  columns  of the  matrices  A and B,  and on  entry  with
!           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
!           of rows of the matrices  A and B.  K must be at least  zero.
!           Unchanged on exit.

!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.

!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
!           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by n  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.

!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
!           then  LDA must be at least  max( 1, n ), otherwise  LDA must
!           be at least  max( 1, k ).
!           Unchanged on exit.

!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
!           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
!           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  k by n  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.

!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
!           then  LDB must be at least  max( 1, n ), otherwise  LDB must
!           be at least  max( 1, k ).
!           Unchanged on exit.

!  BETA   - DOUBLE PRECISION.
!           On entry, BETA specifies the scalar beta.
!           Unchanged on exit.

!  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
!           upper triangular part of the array C must contain the upper
!           triangular part  of the  symmetric matrix  and the strictly
!           lower triangular part of C is not referenced.  On exit, the
!           upper triangular part of the array  C is overwritten by the
!           upper triangular part of the updated matrix.
!           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
!           lower triangular part of the array C must contain the lower
!           triangular part  of the  symmetric matrix  and the strictly
!           upper triangular part of C is not referenced.  On exit, the
!           lower triangular part of the array  C is overwritten by the
!           lower triangular part of the updated matrix.

!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, n ).
!           Unchanged on exit.


!  Level 3 Blas routine.


!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.


!     .. External Functions ..
! LOGICAL ::            lsame
! EXTERNAL           lsame
!     .. External Subroutines ..
! EXTERNAL           xerbla
!     .. Intrinsic Functions ..
! INTRINSIC          MAX
!     .. Local Scalars ..
LOGICAL   :: upper
INTEGER   :: i, info, j, l, nrowa
REAL (wp) :: temp1, temp2
!     ..
!     .. Executable Statements ..

!     Test the input parameters.

IF( lsame( trans, 'N' ) )THEN
  nrowa = n
ELSE
  nrowa = k
END IF
upper = lsame( uplo, 'U' )

info = 0
IF(      ( .NOT.upper               ).AND.  &
  ( .NOT.lsame( uplo , 'L' ) )      )THEN
  info = 1
ELSE IF( ( .NOT.lsame( trans, 'N' ) ).AND.  &
  ( .NOT.lsame( trans, 'T' ) ).AND.( .NOT.lsame( trans, 'C' ) )      )THEN
  info = 2
ELSE IF( n  < 0               )THEN
  info = 3
ELSE IF( k  < 0               )THEN
  info = 4
ELSE IF( lda < MAX( 1, nrowa ) )THEN
  info = 7
ELSE IF( ldb < MAX( 1, nrowa ) )THEN
  info = 9
ELSE IF( ldc < MAX( 1, n     ) )THEN
  info = 12
END IF
IF( info /= 0 )THEN
  CALL erinfo(info, 'DSYR2K')
  RETURN
END IF

!     Quick return if possible.

IF( ( n == 0 ).OR.  &
( ( ( alpha == zero ).OR.( k == 0 ) ).AND.( beta == one ) ) )RETURN

!     And when  alpha.eq.zero.

IF( alpha == zero )THEN
  IF( upper )THEN
    IF( beta == zero )THEN
      DO j = 1, n
        c( 1:j, j ) = zero
      END DO
    ELSE
      DO j = 1, n
        c( 1:j, j ) = beta*c( 1:j, j )
      END DO
    END IF
  ELSE
    IF( beta == zero )THEN
      DO j = 1, n
        c( j:n, j ) = zero
      END DO
    ELSE
      DO j = 1, n
        c( j:n, j ) = beta*c( j:n, j )
      END DO
    END IF
  END IF
  RETURN
END IF

!     Start the operations.

IF( lsame( trans, 'N' ) )THEN

!        Form  C := alpha*A*B' + alpha*B*A' + C.

  IF( upper )THEN
    DO j = 1, n
      IF( beta == zero )THEN
        c( 1:j, j ) = zero
      ELSE IF( beta /= one )THEN
        c( 1:j, j ) = beta*c( 1:j, j )
      END IF
      DO l = 1, k
        IF( ( a( j, l ) /= zero ).OR. ( b( j, l ) /= zero )     )THEN
          temp1 = alpha*b( j, l )
          temp2 = alpha*a( j, l )
          c( 1:j, j ) = c( 1:j, j ) + a( 1:j, l )*temp1 + b( 1:j, l )*temp2
        END IF
      END DO
    END DO
  ELSE
    DO j = 1, n
      IF( beta == zero )THEN
        c( j:n, j ) = zero
      ELSE IF( beta /= one )THEN
        c( j:n, j ) = beta*c( j:n, j )
      END IF
      DO l = 1, k
        IF( ( a( j, l ) /= zero ).OR. ( b( j, l ) /= zero )     )THEN
          temp1 = alpha*b( j, l )
          temp2 = alpha*a( j, l )
          c( j:n, j ) = c( j:n, j ) + a( j:n, l )*temp1 + b( j:n, l )*temp2
        END IF
      END DO
    END DO
  END IF
ELSE

!        Form  C := alpha*A'*B + alpha*B'*A + C.

  IF( upper )THEN
    DO j = 1, n
      DO i = 1, j
        temp1 = zero
        temp2 = zero
        DO l = 1, k
          temp1 = temp1 + a( l, i )*b( l, j )
          temp2 = temp2 + b( l, i )*a( l, j )
        END DO
        IF( beta == zero )THEN
          c( i, j ) = alpha*temp1 + alpha*temp2
        ELSE
          c( i, j ) = beta *c( i, j ) +alpha*temp1 + alpha*temp2
        END IF
      END DO
    END DO
  ELSE
    DO j = 1, n
      DO i = j, n
        temp1 = zero
        temp2 = zero
        DO l = 1, k
          temp1 = temp1 + a( l, i )*b( l, j )
          temp2 = temp2 + b( l, i )*a( l, j )
        END DO
        IF( beta == zero )THEN
          c( i, j ) = alpha*temp1 + alpha*temp2
        ELSE
          c( i, j ) = beta *c( i, j ) +alpha*temp1 + alpha*temp2
        END IF
      END DO
    END DO
  END IF
END IF

RETURN

!     End of DSYR2K.

END SUBROUTINE dsyr2k


SUBROUTINE dtrsm ( side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb )
! ELF90 translation by Alan Miller   31-Aug-1997

!     .. Scalar Arguments ..
CHARACTER (LEN=1), INTENT(IN) :: side, uplo, transa, diag
INTEGER, INTENT(IN)           :: m, n, lda, ldb
REAL (wp), INTENT(IN)         :: alpha
!     .. Array Arguments ..
REAL (wp), INTENT(IN)         :: a(:,:)
REAL (wp), INTENT(IN OUT)     :: b(:,:)
!     ..

!  Purpose
!  =======

!  DTRSM  solves one of the matrix equations

!     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,

!  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of

!     op( A ) = A   or   op( A ) = A'.

!  The matrix X is overwritten on B.

!  Parameters
!  ==========

!  SIDE   - CHARACTER*1.
!           On entry, SIDE specifies whether op( A ) appears on the left
!           or right of X as follows:

!              SIDE = 'L' or 'l'   op( A )*X = alpha*B.

!              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.

!           Unchanged on exit.

!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:

!              UPLO = 'U' or 'u'   A is an upper triangular matrix.

!              UPLO = 'L' or 'l'   A is a lower triangular matrix.

!           Unchanged on exit.

!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:

!              TRANSA = 'N' or 'n'   op( A ) = A.

!              TRANSA = 'T' or 't'   op( A ) = A'.

!              TRANSA = 'C' or 'c'   op( A ) = A'.

!           Unchanged on exit.

!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:

!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.

!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.

!           Unchanged on exit.

!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.

!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.

!  ALPHA  - DOUBLE PRECISION.
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.

!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.

!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.

!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain  the  right-hand  side  matrix  B,  and  on exit  is
!           overwritten by the solution matrix  X.

!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.


!  Level 3 Blas routine.


!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.


!     .. External Functions ..
! LOGICAL ::            lsame
! EXTERNAL           lsame
!     .. External Subroutines ..
! EXTERNAL           xerbla
!     .. Intrinsic Functions ..
! INTRINSIC          MAX
!     .. Local Scalars ..
LOGICAL   :: lside, nounit, upper
INTEGER   :: i, info, j, k, nrowa
REAL (wp) :: temp
!     ..
!     .. Executable Statements ..

!     Test the input parameters.

lside  = lsame( side  , 'L' )
IF( lside )THEN
  nrowa = m
ELSE
  nrowa = n
END IF
nounit = lsame( diag  , 'N' )
upper  = lsame( uplo  , 'U' )

info   = 0
IF(      ( .NOT.lside                ).AND.  &
  ( .NOT.lsame( side  , 'R' ) )      )THEN
  info = 1
ELSE IF( ( .NOT.upper                ).AND.  &
  ( .NOT.lsame( uplo  , 'L' ) )      )THEN
  info = 2
ELSE IF( ( .NOT.lsame( transa, 'N' ) ).AND.  &
  ( .NOT.lsame( transa, 'T' ) ).AND.( .NOT.lsame( transa, 'C' ) )      )THEN
  info = 3
ELSE IF( ( .NOT.lsame( diag  , 'U' ) ).AND.  &
  ( .NOT.lsame( diag  , 'N' ) )      )THEN
  info = 4
ELSE IF( m  < 0               )THEN
  info = 5
ELSE IF( n  < 0               )THEN
  info = 6
ELSE IF( lda < MAX( 1, nrowa ) )THEN
  info = 9
ELSE IF( ldb < MAX( 1, m     ) )THEN
  info = 11
END IF
IF( info /= 0 )THEN
  CALL erinfo(info, 'DTRSM ')
  RETURN
END IF

!     Quick return if possible.

IF( n == 0 ) RETURN

!     And when  alpha.eq.zero.

IF( alpha == zero )THEN
  b( 1:m, 1:n ) = zero
  RETURN
END IF

!     Start the operations.

IF( lside )THEN
  IF( lsame( transa, 'N' ) )THEN

!           Form  B := alpha*inv( A )*B.

    IF( upper )THEN
      DO j = 1, n
        IF( alpha /= one )THEN
          b( 1:m, j ) = alpha*b( 1:m, j )
        END IF
        DO k = m, 1, -1
          IF( b( k, j ) /= zero )THEN
            IF( nounit ) b( k, j ) = b( k, j )/a( k, k )
            b( 1:k-1, j ) = b( 1:k-1, j ) - b( k, j )*a( 1:k-1, k )
          END IF
        END DO
      END DO
    ELSE
      DO j = 1, n
        IF( alpha /= one )THEN
          b( 1:m, j ) = alpha*b( 1:m, j )
        END IF
        DO k = 1, m
          IF( b( k, j ) /= zero )THEN
            IF( nounit ) b( k, j ) = b( k, j )/a( k, k )
            b( k+1:m, j ) = b( k+1:m, j ) - b( k, j )*a( k+1:m, k )
          END IF
        END DO
      END DO
    END IF
  ELSE

!           Form  B := alpha*inv( A' )*B.

    IF( upper )THEN
      DO j = 1, n
        DO i = 1, m
          temp = alpha*b( i, j ) - DOT_PRODUCT( a(:i-1,i), b(:i-1,j) )
          IF( nounit ) temp = temp/a( i, i )
          b( i, j ) = temp
        END DO
      END DO
    ELSE
      DO j = 1, n
        DO i = m, 1, -1
          temp = alpha*b( i, j ) - DOT_PRODUCT( a(i+1:m,i), b(i+1:m,j) )
          IF( nounit ) temp = temp/a( i, i )
          b( i, j ) = temp
        END DO
      END DO
    END IF
  END IF
ELSE
  IF( lsame( transa, 'N' ) )THEN

!           Form  B := alpha*B*inv( A ).

    IF( upper )THEN
      DO j = 1, n
        IF( alpha /= one )THEN
          b( 1:m, j ) = alpha*b( 1:m, j )
        END IF
        DO k = 1, j - 1
          IF( a( k, j ) /= zero )THEN
            b( 1:m, j ) = b( 1:m, j ) - a( k, j )*b( 1:m, k )
          END IF
        END DO
        IF( nounit )THEN
          temp = one/a( j, j )
          b( 1:m, j ) = temp*b( 1:m, j )
        END IF
      END DO
    ELSE
      DO j = n, 1, -1
        IF( alpha /= one )THEN
          b( 1:m, j ) = alpha*b( 1:m, j )
        END IF
        DO k = j + 1, n
          IF( a( k, j ) /= zero )THEN
            b( 1:m, j ) = b( 1:m, j ) - a( k, j )*b( 1:m, k )
          END IF
        END DO
        IF( nounit )THEN
          temp = one/a( j, j )
          b( 1:m, j ) = temp*b( 1:m, j )
        END IF
      END DO
    END IF
  ELSE

!           Form  B := alpha*B*inv( A' ).

    IF( upper )THEN
      DO k = n, 1, -1
        IF( nounit )THEN
          temp = one/a( k, k )
          b( 1:m, k ) = temp*b( 1:m, k )
        END IF
        DO j = 1, k - 1
          IF( a( j, k ) /= zero )THEN
            temp = a( j, k )
            b( 1:m, j ) = b( 1:m, j ) - temp*b( 1:m, k )
          END IF
        END DO
        IF( alpha /= one )THEN
          b( 1:m, k ) = alpha*b( 1:m, k )
        END IF
      END DO
    ELSE
      DO k = 1, n
        IF( nounit )THEN
          temp = one/a( k, k )
          b( 1:m, k ) = temp*b( 1:m, k )
        END IF
        DO j = k + 1, n
          IF( a( j, k ) /= zero )THEN
            temp = a( j, k )
            b( 1:m, j ) = b( 1:m, j ) - temp*b( 1:m, k )
          END IF
        END DO
        IF( alpha /= one )THEN
          b( 1:m, k ) = alpha*b( 1:m, k )
        END IF
      END DO
    END IF
  END IF
END IF

RETURN

!     End of DTRSM .

END SUBROUTINE dtrsm


SUBROUTINE dtrmm ( side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb )
! ELF90 translation by Alan Miller   31-Aug-1997

!     .. Scalar Arguments ..
CHARACTER (LEN=1), INTENT(IN) :: side, uplo, transa, diag
INTEGER, INTENT(IN)           :: m, n, lda, ldb
REAL (wp), INTENT(IN)         :: alpha
!     .. Array Arguments ..
REAL (wp), INTENT(IN)         :: a(:,:)
REAL (wp), INTENT(IN OUT)     :: b(:,:)
!     ..

!  Purpose
!  =======

!  DTRMM  performs one of the matrix-matrix operations

!     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),

!  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
!  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of

!     op( A ) = A   or   op( A ) = A'.

!  Parameters
!  ==========

!  SIDE   - CHARACTER*1.
!           On entry,  SIDE specifies whether  op( A ) multiplies B from
!           the left or right as follows:

!              SIDE = 'L' or 'l'   B := alpha*op( A )*B.

!              SIDE = 'R' or 'r'   B := alpha*B*op( A ).

!           Unchanged on exit.

!  UPLO   - CHARACTER*1.
!           On entry, UPLO specifies whether the matrix A is an upper or
!           lower triangular matrix as follows:

!              UPLO = 'U' or 'u'   A is an upper triangular matrix.

!              UPLO = 'L' or 'l'   A is a lower triangular matrix.

!           Unchanged on exit.

!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:

!              TRANSA = 'N' or 'n'   op( A ) = A.

!              TRANSA = 'T' or 't'   op( A ) = A'.

!              TRANSA = 'C' or 'c'   op( A ) = A'.

!           Unchanged on exit.

!  DIAG   - CHARACTER*1.
!           On entry, DIAG specifies whether or not A is unit triangular
!           as follows:

!              DIAG = 'U' or 'u'   A is assumed to be unit triangular.

!              DIAG = 'N' or 'n'   A is not assumed to be unit
!                                  triangular.

!           Unchanged on exit.

!  M      - INTEGER.
!           On entry, M specifies the number of rows of B. M must be at
!           least zero.
!           Unchanged on exit.

!  N      - INTEGER.
!           On entry, N specifies the number of columns of B.  N must be
!           at least zero.
!           Unchanged on exit.

!  ALPHA  - DOUBLE PRECISION.
!           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!           zero then  A is not referenced and  B need not be set before
!           entry.
!           Unchanged on exit.

!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
!           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
!           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!           upper triangular part of the array  A must contain the upper
!           triangular matrix  and the strictly lower triangular part of
!           A is not referenced.
!           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!           lower triangular part of the array  A must contain the lower
!           triangular matrix  and the strictly upper triangular part of
!           A is not referenced.
!           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!           A  are not referenced either,  but are assumed to be  unity.
!           Unchanged on exit.

!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!           then LDA must be at least max( 1, n ).
!           Unchanged on exit.

!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
!           Before entry,  the leading  m by n part of the array  B must
!           contain the matrix  B,  and  on exit  is overwritten  by the
!           transformed matrix.

!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in  the  calling  (sub)  program.   LDB  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.


!  Level 3 Blas routine.

!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.


!     .. External Functions ..
! LOGICAL ::            lsame
! EXTERNAL           lsame
!     .. External Subroutines ..
! EXTERNAL           xerbla
!     .. Intrinsic Functions ..
! INTRINSIC          MAX
!     .. Local Scalars ..
LOGICAL ::            lside, nounit, upper
INTEGER ::            i, info, j, k, nrowa
REAL (wp) ::   temp
!     ..
!     .. Executable Statements ..

!     Test the input parameters.

lside  = lsame( side  , 'L' )
IF( lside )THEN
  nrowa = m
ELSE
  nrowa = n
END IF
nounit = lsame( diag  , 'N' )
upper  = lsame( uplo  , 'U' )

info   = 0
IF(      ( .NOT.lside                ).AND.  &
  ( .NOT.lsame( side  , 'R' ) )      )THEN
  info = 1
ELSE IF( ( .NOT.upper                ).AND.  &
  ( .NOT.lsame( uplo  , 'L' ) )      )THEN
  info = 2
ELSE IF( ( .NOT.lsame( transa, 'N' ) ).AND.  &
  ( .NOT.lsame( transa, 'T' ) ).AND.( .NOT.lsame( transa, 'C' ) )      )THEN
  info = 3
ELSE IF( ( .NOT.lsame( diag  , 'U' ) ).AND.  &
  ( .NOT.lsame( diag  , 'N' ) )      )THEN
  info = 4
ELSE IF( m  < 0               )THEN
  info = 5
ELSE IF( n  < 0               )THEN
  info = 6
ELSE IF( lda < MAX( 1, nrowa ) )THEN
  info = 9
ELSE IF( ldb < MAX( 1, m     ) )THEN
  info = 11
END IF
IF( info /= 0 )THEN
  CALL erinfo(info, 'DTRMM ')
  RETURN
END IF

!     Quick return if possible.

IF( n == 0 ) RETURN

!     And when  alpha = zero.

IF( alpha == zero )THEN
  b( 1:m, 1:n ) = zero
  RETURN
END IF

!     Start the operations.

IF( lside )THEN
  IF( lsame( transa, 'N' ) )THEN

!           Form  B := alpha*A*B.

    IF( upper )THEN
      DO j = 1, n
        DO k = 1, m
          IF( b( k, j ) /= zero )THEN
            temp = alpha*b( k, j )
            b( 1:k-1, j ) = b( 1:k-1, j ) + temp*a( 1:k-1, k )
            IF( nounit ) temp = temp*a( k, k )
            b( k, j ) = temp
          END IF
        END DO
      END DO
    ELSE
      DO j = 1, n
        DO k = m, 1, -1
          IF( b( k, j ) /= zero )THEN
            temp      = alpha*b( k, j )
            b( k, j ) = temp
            IF( nounit ) b( k, j ) = b( k, j )*a( k, k )
            b( k+1:m, j ) = b( k+1:m, j ) + temp*a( k+1:m, k )
          END IF
        END DO
      END DO
    END IF
  ELSE

!           Form  B := alpha*B*A'.

    IF( upper )THEN
      DO j = 1, n
        DO i = m, 1, -1
          temp = b( i, j )
          IF( nounit ) temp = temp*a( i, i )
          temp = temp + DOT_PRODUCT( a( 1:i-1, i ), b( 1:i-1, j ) )
          b( i, j ) = alpha*temp
        END DO
      END DO
    ELSE
      DO j = 1, n
        DO i = 1, m
          temp = b( i, j )
          IF( nounit ) temp = temp*a( i, i )
          temp = temp + DOT_PRODUCT( a( i+1:m, i ), b( i+1:m, j ) )
          b( i, j ) = alpha*temp
        END DO
      END DO
    END IF
  END IF
ELSE
  IF( lsame( transa, 'N' ) )THEN

!           Form  B := alpha*B*A.

    IF( upper )THEN
      DO j = n, 1, -1
        temp = alpha
        IF( nounit ) temp = temp*a( j, j )
        b( 1:m, j ) = temp*b( 1:m, j )
        DO k = 1, j - 1
          IF( a( k, j ) /= zero )THEN
            temp = alpha*a( k, j )
            b( 1:m, j ) = b( 1:m, j ) + temp*b( 1:m, k )
          END IF
        END DO
      END DO
    ELSE
      DO j = 1, n
        temp = alpha
        IF( nounit ) temp = temp*a( j, j )
        b( 1:m, j ) = temp*b( 1:m, j )
        DO k = j + 1, n
          IF( a( k, j ) /= zero )THEN
            temp = alpha*a( k, j )
            b( 1:m, j ) = b( 1:m, j ) + temp*b( 1:m, k )
          END IF
        END DO
      END DO
    END IF
  ELSE

!           Form  B := alpha*B*A'.

    IF( upper )THEN
      DO k = 1, n
        DO j = 1, k - 1
          IF( a( j, k ) /= zero )THEN
            temp = alpha*a( j, k )
            b( 1:m, j ) = b( 1:m, j ) + temp*b( 1:m, k )
          END IF
        END DO
        temp = alpha
        IF( nounit ) temp = temp*a( k, k )
        IF( temp /= one )THEN
          b( 1:m, k ) = temp*b( 1:m, k )
        END IF
      END DO
    ELSE
      DO k = n, 1, -1
        DO j = k + 1, n
          IF( a( j, k ) /= zero )THEN
            temp = alpha*a( j, k )
            b( 1:m, j ) = b( 1:m, j ) + temp*b( 1:m, k )
          END IF
        END DO
        temp = alpha
        IF( nounit ) temp = temp*a( k, k )
        IF( temp /= one )THEN
          b( 1:m, k ) = temp*b( 1:m, k )
        END IF
      END DO
    END IF
  END IF
END IF

RETURN

!     End of DTRMM .

END SUBROUTINE dtrmm


SUBROUTINE dgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb,  &
                   beta, c, ldc )
! ELF90 translation by Alan Miller   31-Aug-1997

!     .. Scalar Arguments ..
CHARACTER (LEN=1), INTENT(IN) :: transa, transb
INTEGER, INTENT(IN)           :: m, n, k, lda, ldb, ldc
REAL (wp), INTENT(IN)         :: alpha, beta
!     .. Array Arguments ..
REAL (wp), INTENT(IN)         :: a(:,:), b(:,:)
REAL (wp), INTENT(IN OUT)     :: c(:,:)
!     ..

!  Purpose
!  =======

!  DGEMM  performs one of the matrix-matrix operations

!     C := alpha*op( A )*op( B ) + beta*C,

!  where  op( X ) is one of

!     op( X ) = X   or   op( X ) = X',

!  alpha and beta are scalars, and A, B and C are matrices, with op( A )
!  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.

!  Parameters
!  ==========

!  TRANSA - CHARACTER*1.
!           On entry, TRANSA specifies the form of op( A ) to be used in
!           the matrix multiplication as follows:

!              TRANSA = 'N' or 'n',  op( A ) = A.

!              TRANSA = 'T' or 't',  op( A ) = A'.

!              TRANSA = 'C' or 'c',  op( A ) = A'.

!           Unchanged on exit.

!  TRANSB - CHARACTER*1.
!           On entry, TRANSB specifies the form of op( B ) to be used in
!           the matrix multiplication as follows:

!              TRANSB = 'N' or 'n',  op( B ) = B.

!              TRANSB = 'T' or 't',  op( B ) = B'.

!              TRANSB = 'C' or 'c',  op( B ) = B'.

!           Unchanged on exit.

!  M      - INTEGER.
!           On entry,  M  specifies  the number  of rows  of the  matrix
!           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!           Unchanged on exit.

!  N      - INTEGER.
!           On entry,  N  specifies the number  of columns of the matrix
!           op( B ) and the number of columns of the matrix C. N must be
!           at least zero.
!           Unchanged on exit.

!  K      - INTEGER.
!           On entry,  K  specifies  the number of columns of the matrix
!           op( A ) and the number of rows of the matrix op( B ). K must
!           be at least  zero.
!           Unchanged on exit.

!  ALPHA  - DOUBLE PRECISION.
!           On entry, ALPHA specifies the scalar alpha.
!           Unchanged on exit.

!  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
!           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!           part of the array  A  must contain the matrix  A,  otherwise
!           the leading  k by m  part of the array  A  must contain  the
!           matrix A.
!           Unchanged on exit.

!  LDA    - INTEGER.
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!           least  max( 1, k ).
!           Unchanged on exit.

!  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
!           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!           part of the array  B  must contain the matrix  B,  otherwise
!           the leading  n by k  part of the array  B  must contain  the
!           matrix B.
!           Unchanged on exit.

!  LDB    - INTEGER.
!           On entry, LDB specifies the first dimension of B as declared
!           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!           least  max( 1, n ).
!           Unchanged on exit.

!  BETA   - DOUBLE PRECISION.
!           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!           supplied as zero then C need not be set on input.
!           Unchanged on exit.

!  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
!           Before entry, the leading  m by n  part of the array  C must
!           contain the matrix  C,  except when  beta  is zero, in which
!           case C need not be set on entry.
!           On exit, the array  C  is overwritten by the  m by n  matrix
!           ( alpha*op( A )*op( B ) + beta*C ).

!  LDC    - INTEGER.
!           On entry, LDC specifies the first dimension of C as declared
!           in  the  calling  (sub)  program.   LDC  must  be  at  least
!           max( 1, m ).
!           Unchanged on exit.


!  Level 3 Blas routine.

!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.


!     .. External Functions ..
! LOGICAL ::            lsame
! EXTERNAL           lsame
!     .. External Subroutines ..
! EXTERNAL           xerbla
!     .. Intrinsic Functions ..
! INTRINSIC          MAX
!     .. Local Scalars ..
LOGICAL ::            nota, notb
INTEGER ::            i, info, j, l, nrowa, nrowb
REAL (wp) ::   temp
!     ..
!     .. Executable Statements ..

!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
!     and  columns of  A  and the  number of  rows  of  B  respectively.

nota  = lsame( transa, 'N' )
notb  = lsame( transb, 'N' )
IF( nota )THEN
  nrowa = m
ELSE
  nrowa = k
END IF
IF( notb )THEN
  nrowb = k
ELSE
  nrowb = n
END IF

!     Test the input parameters.

info = 0
IF(      ( .NOT.nota                 ).AND.  &
  ( .NOT.lsame( transa, 'C' ) ).AND.( .NOT.lsame( transa, 'T' ) )      )THEN
  info = 1
ELSE IF( ( .NOT.notb                 ).AND.  &
  ( .NOT.lsame( transb, 'C' ) ).AND.( .NOT.lsame( transb, 'T' ) )      )THEN
  info = 2
ELSE IF( m  < 0               )THEN
  info = 3
ELSE IF( n  < 0               )THEN
  info = 4
ELSE IF( k  < 0               )THEN
  info = 5
ELSE IF( lda < MAX( 1, nrowa ) )THEN
  info = 8
ELSE IF( ldb < MAX( 1, nrowb ) )THEN
  info = 10
ELSE IF( ldc < MAX( 1, m     ) )THEN
  info = 13
END IF
IF( info /= 0 )THEN
  CALL erinfo(info, 'DGEMM ')
  RETURN
END IF

!     Quick return if possible.

IF( ( m == 0 ).OR.( n == 0 ).OR.  &
( ( ( alpha == zero ).OR.( k == 0 ) ).AND.( beta == one ) ) )RETURN

!     And if  alpha.eq.zero.

IF( alpha == zero ) THEN
  IF( beta == zero ) THEN
    c(1:m,1:n) = zero
  ELSE
    c(1:m,1:n) = beta * c(1:m,1:n)
  END IF
  RETURN
END IF

!     Start the operations.

IF( notb )THEN
  IF( nota )THEN

!           Form  C := alpha*A*B + beta*C.

    DO j = 1, n
      IF( beta == zero )THEN
        DO i = 1, m
          c( i, j ) = zero
        END DO
      ELSE IF( beta /= one )THEN
        c( 1:m, j ) = beta*c( 1:m, j )
      END IF
      DO l = 1, k
        IF( b( l, j ) /= zero )THEN
          temp = alpha*b( l, j )
          c( 1:m, j ) = c( 1:m, j ) + temp*a( 1:m, l )
        END IF
      END DO
    END DO
  ELSE

!           Form  C := alpha*A'*B + beta*C

    DO j = 1, n
      DO i = 1, m
        temp = DOT_PRODUCT( a(1:k,i), b(1:k,j) )
        IF( beta == zero )THEN
          c( i, j ) = alpha*temp
        ELSE
          c( i, j ) = alpha*temp + beta*c( i, j )
        END IF
      END DO
    END DO
  END IF
ELSE
  IF( nota )THEN

!           Form  C := alpha*A*B' + beta*C

    DO j = 1, n
      IF( beta == zero )THEN
        c( 1:m, j ) = zero
      ELSE IF( beta /= one )THEN
        c( 1:m, j ) = beta*c( 1:m, j )
      END IF
      DO l = 1, k
        IF( b( j, l ) /= zero )THEN
          temp = alpha*b( j, l )
          c( 1:m, j ) = c( 1:m, j ) + temp*a( 1:m, l )
        END IF
      END DO
    END DO
  ELSE

!           Form  C := alpha*A'*B' + beta*C

    DO j = 1, n
      DO i = 1, m
        temp = DOT_PRODUCT( a(1:k,i), b(j,1:k) )
        IF( beta == zero )THEN
          c( i, j ) = alpha*temp
        ELSE
          c( i, j ) = alpha*temp + beta*c( i, j )
        END IF
      END DO
    END DO
  END IF
END IF

RETURN

!     End of DGEMM .

END SUBROUTINE dgemm



END MODULE dblas3


MODULE dla
! LAPACK auxiliary routines

! Translated by Alan Miller
! Alan.Miller @ vic.cmis.csiro.au    URL: www.ozemail.com.au/~milleraj

! Latest revision - 17 August 1998

USE la_precision, ONLY: wp => dp
USE la_auxmod
USE dblas
USE dblas2, ONLY: dgemv, dsymv, dtrmv, dtrsv, dger
USE dblas3, ONLY: dgemm, dtrmm
IMPLICIT NONE

CONTAINS



SUBROUTINE dlarfg( n, alpha, x, incx, tau )

!  -- LAPACK auxiliary routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997

!     .. Scalar Arguments ..
INTEGER, INTENT(IN)       :: incx, n
REAL (wp), INTENT(IN OUT) :: alpha
REAL (wp), INTENT(OUT)    :: tau
!     ..
!     .. Array Arguments ..
REAL (wp), INTENT(IN OUT) :: x(:)
!     ..

!  Purpose
!  =======

!  DLARFG generates a real elementary reflector H of order n, such
!  that

!        H * ( alpha ) = ( beta ),   H' * H = I.
!            (   x   )   (   0  )

!  where alpha and beta are scalars, and x is an (n-1)-element real
!  vector. H is represented in the form

!        H = I - tau * ( 1 ) * ( 1 v' ) ,
!                      ( v )

!  where tau is a real scalar and v is a real (n-1)-element
!  vector.

!  If the elements of x are all zero, then tau = 0 and H is taken to be
!  the unit matrix.

!  Otherwise  1 <= tau <= 2.

!  Arguments
!  =========

!  N       (input) INTEGER
!          The order of the elementary reflector.

!  ALPHA   (input/output) DOUBLE PRECISION
!          On entry, the value alpha.
!          On exit, it is overwritten with the value beta.

!  X       (input/output) DOUBLE PRECISION array, dimension
!                         (1+(N-2)*abs(INCX))
!          On entry, the vector x.
!          On exit, it is overwritten with the vector v.

!  INCX    (input) INTEGER
!          The increment between elements of X. INCX <> 0.

!  TAU     (output) DOUBLE PRECISION
!          The value tau.

!  =====================================================================

!     .. Local Scalars ..
INTEGER   :: knt
REAL (wp) :: beta, rsafmn, safmin, xnorm
!     ..
!     .. External Functions ..
! DOUBLE PRECISION ::   dlamch, dlapy2, dnrm2
! EXTERNAL           dlamch, dlapy2, dnrm2
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC          ABS, SIGN
!     ..
!     .. External Subroutines ..
! EXTERNAL           dscal
!     ..
!     .. Executable Statements ..

IF( n <= 1 ) THEN
  tau = zero
  RETURN
END IF

xnorm = dnrm2( n-1, x, incx )

IF( xnorm == zero ) THEN

!        H  =  I

  tau = zero
ELSE

!        general case

  beta = -SIGN( dlapy2( alpha, xnorm ), alpha )
  safmin = dlamch( 'S' )
  IF( ABS( beta ) < safmin ) THEN

!           XNORM, BETA may be inaccurate; scale X and recompute them

    rsafmn = one / safmin
    knt = 0

    10 knt = knt + 1
    CALL dscal( n-1, rsafmn, x, incx )
    beta = beta*rsafmn
    alpha = alpha*rsafmn
    IF( ABS( beta ) < safmin ) GO TO 10

!           New BETA is at most 1, at least SAFMIN

    xnorm = dnrm2( n-1, x, incx )
    beta = -SIGN( dlapy2( alpha, xnorm ), alpha )
    tau = ( beta-alpha ) / beta
    CALL dscal( n-1, one / ( alpha-beta ), x, incx )

!           If ALPHA is subnormal, it may lose relative accuracy

    alpha = beta * safmin**knt
  ELSE
    tau = ( beta-alpha ) / beta
    CALL dscal( n-1, one / ( alpha-beta ), x, incx )
    alpha = beta
  END IF
END IF

RETURN

!     End of DLARFG

END SUBROUTINE dlarfg



SUBROUTINE dlartg( f, g, cs, sn, r )

!  -- LAPACK auxiliary routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994

! ELF90 translation by Alan Miller    30-Nov-1997

!     .. Scalar Arguments ..
REAL (wp), INTENT(IN)  :: f, g
REAL (wp), INTENT(OUT) :: cs, r, sn
!     ..

!  Purpose
!  =======

!  DLARTG generate a plane rotation so that

!     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
!     [ -SN  CS  ]     [ G ]     [ 0 ]

!  This is a slower, more accurate version of the BLAS1 routine DROTG,
!  with the following other differences:
!     F and G are unchanged on return.
!     If G=0, then CS=1 and SN=0.
!     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
!        floating point operations (saves work in DBDSQR when
!        there are zeros on the diagonal).

!  If F exceeds G in magnitude, CS will be positive.

!  Arguments
!  =========

!  F       (input) DOUBLE PRECISION
!          The first component of vector to be rotated.

!  G       (input) DOUBLE PRECISION
!          The second component of vector to be rotated.

!  CS      (output) DOUBLE PRECISION
!          The cosine of the rotation.

!  SN      (output) DOUBLE PRECISION
!          The sine of the rotation.

!  R       (output) DOUBLE PRECISION
!          The nonzero component of the rotated vector.

!  =====================================================================

!     .. Local Scalars ..
LOGICAL, SAVE        :: first = .TRUE.
INTEGER              :: count, i
REAL (wp)            :: eps, f1, g1, scale
REAL (wp), SAVE      :: safmin, safmn2, safmx2
REAL (wp), PARAMETER :: two = 2.0_wp
!     ..
!     .. External Functions ..
! DOUBLE PRECISION ::   dlamch
! EXTERNAL           dlamch
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC          ABS, INT, LOG, MAX, SQRT
!     ..
!     .. Save statement ..
! SAVE               first, safmx2, safmin, safmn2
!     ..
!     .. Data statements ..
! DATA               first / .true. /
!     ..
!     .. Executable Statements ..

IF( first ) THEN
  first = .false.
  safmin = dlamch( 'S' )
  eps = dlamch( 'E' )
  safmn2 = dlamch( 'B' )**INT( LOG( safmin / eps ) /  &
  LOG( dlamch( 'B' ) ) / two )
  safmx2 = one / safmn2
END IF
IF( g == zero ) THEN
  cs = one
  sn = zero
  r = f
ELSE IF( f == zero ) THEN
  cs = zero
  sn = one
  r = g
ELSE
  f1 = f
  g1 = g
  scale = MAX( ABS( f1 ), ABS( g1 ) )
  IF( scale >= safmx2 ) THEN
    count = 0

    10 count = count + 1
    f1 = f1*safmn2
    g1 = g1*safmn2
    scale = MAX( ABS( f1 ), ABS( g1 ) )
    IF( scale >= safmx2 ) GO TO 10
    r = SQRT( f1**2 + g1**2 )
    cs = f1 / r
    sn = g1 / r
    DO i = 1, count
      r = r*safmx2
    END DO
  ELSE IF( scale <= safmn2 ) THEN
    count = 0

    30 count = count + 1
    f1 = f1*safmx2
    g1 = g1*safmx2
    scale = MAX( ABS( f1 ), ABS( g1 ) )
    IF( scale <= safmn2 ) GO TO 30
    r = SQRT( f1**2 + g1**2 )
    cs = f1 / r
    sn = g1 / r
    DO i = 1, count
      r = r*safmn2
    END DO
  ELSE
    r = SQRT( f1**2 + g1**2 )
    cs = f1 / r
    sn = g1 / r
  END IF
  IF( ABS( f ) > ABS( g ) .AND. cs < zero ) THEN
    cs = -cs
    sn = -sn
    r = -r
  END IF
END IF
RETURN

!     End of DLARTG

END SUBROUTINE dlartg



SUBROUTINE dlae2( a, b, c, rt1, rt2 )

!  -- LAPACK auxiliary routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997

!     .. Scalar Arguments ..
REAL (wp), INTENT(IN)  :: a, b, c
REAL (wp), INTENT(OUT) :: rt1, rt2
!     ..

!  Purpose
!  =======

!  DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix
!     [  A   B  ]
!     [  B   C  ].
!  On return, RT1 is the eigenvalue of larger absolute value, and RT2
!  is the eigenvalue of smaller absolute value.

!  Arguments
!  =========

!  A       (input) DOUBLE PRECISION
!          The (1,1) entry of the 2-by-2 matrix.

!  B       (input) DOUBLE PRECISION
!          The (1,2) and (2,1) entries of the 2-by-2 matrix.

!  C       (input) DOUBLE PRECISION
!          The (2,2) entry of the 2-by-2 matrix.

!  RT1     (output) DOUBLE PRECISION
!          The eigenvalue of larger absolute value.

!  RT2     (output) DOUBLE PRECISION
!          The eigenvalue of smaller absolute value.

!  Further Details
!  ===============

!  RT1 is accurate to a few ulps barring over/underflow.

!  RT2 may be inaccurate if there is massive cancellation in the
!  determinant A*C-B*B; higher precision or correctly rounded or
!  correctly truncated arithmetic would be needed to compute RT2
!  accurately in all cases.

!  Overflow is possible only if RT1 is within a factor of 5 of overflow.
!  Underflow is harmless if the input data is 0 or exceeds
!     underflow_threshold / macheps.


!     .. Parameters ..
REAL (wp), PARAMETER :: two = 2._wp, half = 0.5_wp
!     ..
!     .. Local Scalars ..
REAL (wp) :: ab, acmn, acmx, adf, df, rt, sm, tb
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC          ABS, SQRT
!     ..
!     .. Executable Statements ..

!     Compute the eigenvalues

sm = a + c
df = a - c
adf = ABS( df )
tb = b + b
ab = ABS( tb )
IF( ABS( a ) > ABS( c ) ) THEN
  acmx = a
  acmn = c
ELSE
  acmx = c
  acmn = a
END IF
IF( adf > ab ) THEN
  rt = adf*SQRT( one+( ab / adf )**2 )
ELSE IF( adf < ab ) THEN
  rt = ab*SQRT( one+( adf / ab )**2 )
ELSE

!        Includes case AB=ADF=0

  rt = ab*SQRT( two )
END IF
IF( sm < zero ) THEN
  rt1 = half*( sm-rt )

!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.

  rt2 = ( acmx / rt1 )*acmn - ( b / rt1 )*b
ELSE IF( sm > zero ) THEN
  rt1 = half*( sm+rt )

!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.

  rt2 = ( acmx / rt1 )*acmn - ( b / rt1 )*b
ELSE

!        Includes case RT1 = RT2 = 0

  rt1 = half*rt
  rt2 = -half*rt
END IF
RETURN

!     End of DLAE2

END SUBROUTINE dlae2




SUBROUTINE dlaev2( a, b, c, rt1, rt2, cs1, sn1 )

!  -- LAPACK auxiliary routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997

!     .. Scalar Arguments ..
REAL (wp), INTENT(IN)  :: a, b, c
REAL (wp), INTENT(OUT) :: cs1, rt1, rt2, sn1
!     ..

!  Purpose
!  =======

!  DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix
!     [  A   B  ]
!     [  B   C  ].
!  On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
!  eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
!  eigenvector for RT1, giving the decomposition

!     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
!     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].

!  Arguments
!  =========

!  A       (input) DOUBLE PRECISION
!          The (1,1) entry of the 2-by-2 matrix.

!  B       (input) DOUBLE PRECISION
!          The (1,2) entry and the conjugate of the (2,1) entry of the
!          2-by-2 matrix.

!  C       (input) DOUBLE PRECISION
!          The (2,2) entry of the 2-by-2 matrix.

!  RT1     (output) DOUBLE PRECISION
!          The eigenvalue of larger absolute value.

!  RT2     (output) DOUBLE PRECISION
!          The eigenvalue of smaller absolute value.

!  CS1     (output) DOUBLE PRECISION
!  SN1     (output) DOUBLE PRECISION
!          The vector (CS1, SN1) is a unit right eigenvector for RT1.

!  Further Details
!  ===============

!  RT1 is accurate to a few ulps barring over/underflow.

!  RT2 may be inaccurate if there is massive cancellation in the
!  determinant A*C-B*B; higher precision or correctly rounded or
!  correctly truncated arithmetic would be needed to compute RT2
!  accurately in all cases.

!  CS1 and SN1 are accurate to a few ulps barring over/underflow.

!  Overflow is possible only if RT1 is within a factor of 5 of overflow.
!  Underflow is harmless if the input data is 0 or exceeds
!     underflow_threshold / macheps.


!     .. Parameters ..
REAL (wp), PARAMETER :: two = 2.0_wp, half = 0.5_wp
!     ..
!     .. Local Scalars ..
INTEGER   :: sgn1, sgn2
REAL (wp) :: ab, acmn, acmx, acs, adf, cs, ct, df, rt, sm,tb, tn
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC          ABS, SQRT
!     ..
!     .. Executable Statements ..

!     Compute the eigenvalues

sm = a + c
df = a - c
adf = ABS( df )
tb = b + b
ab = ABS( tb )
IF( ABS( a ) > ABS( c ) ) THEN
  acmx = a
  acmn = c
ELSE
  acmx = c
  acmn = a
END IF
IF( adf > ab ) THEN
  rt = adf*SQRT( one+( ab / adf )**2 )
ELSE IF( adf < ab ) THEN
  rt = ab*SQRT( one+( adf / ab )**2 )
ELSE

!        Includes case AB=ADF=0

  rt = ab*SQRT( two )
END IF
IF( sm < zero ) THEN
  rt1 = half*( sm-rt )
  sgn1 = -1

!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.

  rt2 = ( acmx / rt1 )*acmn - ( b / rt1 )*b
ELSE IF( sm > zero ) THEN
  rt1 = half*( sm+rt )
  sgn1 = 1

!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.

  rt2 = ( acmx / rt1 )*acmn - ( b / rt1 )*b
ELSE

!        Includes case RT1 = RT2 = 0

  rt1 = half*rt
  rt2 = -half*rt
  sgn1 = 1
END IF

!     Compute the eigenvector

IF( df >= zero ) THEN
  cs = df + rt
  sgn2 = 1
ELSE
  cs = df - rt
  sgn2 = -1
END IF
acs = ABS( cs )
IF( acs > ab ) THEN
  ct = -tb / cs
  sn1 = one / SQRT( one+ct*ct )
  cs1 = ct*sn1
ELSE
  IF( ab == zero ) THEN
    cs1 = one
    sn1 = zero
  ELSE
    tn = -cs / tb
    cs1 = one / SQRT( one+tn*tn )
    sn1 = tn*cs1
  END IF
END IF
IF( sgn1 == sgn2 ) THEN
  tn = cs1
  cs1 = -sn1
  sn1 = tn
END IF
RETURN

!     End of DLAEV2

END SUBROUTINE dlaev2




SUBROUTINE dlarf( side, m, n, v, incv, tau, c, ldc )

!  -- LAPACK auxiliary routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997
! N.B. Argument WORK has been removed.

!     .. Scalar Arguments ..
CHARACTER (LEN=1), INTENT(IN) :: side
INTEGER, INTENT(IN)           :: incv, ldc, m, n
REAL (wp), INTENT(IN)         :: tau
!     ..
!     .. Array Arguments ..
REAL (wp), INTENT(IN)         :: v(:)
REAL (wp), INTENT(IN OUT)     :: c(:,:)
!     ..

!  Purpose
!  =======

!  DLARF applies a real elementary reflector H to a real m by n matrix
!  C, from either the left or the right. H is represented in the form

!        H = I - tau * v * v'

!  where tau is a real scalar and v is a real vector.

!  If tau = 0, then H is taken to be the unit matrix.

!  Arguments
!  =========

!  SIDE    (input) CHARACTER*1
!          = 'L': form  H * C
!          = 'R': form  C * H

!  M       (input) INTEGER
!          The number of rows of the matrix C.

!  N       (input) INTEGER
!          The number of columns of the matrix C.

!  V       (input) DOUBLE PRECISION array, dimension
!                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
!                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
!          The vector v in the representation of H. V is not used if
!          TAU = 0.

!  INCV    (input) INTEGER
!          The increment between elements of v. INCV <> 0.

!  TAU     (input) DOUBLE PRECISION
!          The value tau in the representation of H.

!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
!          or C * H if SIDE = 'R'.

!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDC >= max(1,M).

!  =====================================================================

REAL (wp), ALLOCATABLE,save :: work(:)
INTEGER, save :: workd1=-1

!     .. External Subroutines ..
! EXTERNAL           dgemv, dger
!     ..
!     .. External Functions ..
! LOGICAL ::            lsame
! EXTERNAL           lsame
!     ..
!     .. Executable Statements ..

IF( lsame( side, 'L' ) ) THEN
  if (n.ne.workd1) then
    !   only perform expensive allocation if necessary
    if (allocated(work))deallocate(work)
    ALLOCATE( work(n) )
    workd1=n
  endif

!        Form  H * C

  IF( tau /= zero ) THEN

!           w := C' * v

    CALL dgemv( 'Transpose', m, n, one, c, ldc, v, incv, zero, work, 1 )

!           C := C - v * w'

    CALL dger( m, n, -tau, v, incv, work, 1, c, ldc )
  END IF
ELSE
  ALLOCATE( work(m) )

!        Form  C * H

  IF( tau /= zero ) THEN

!           w := C * v

    CALL dgemv( 'No transpose', m, n, one, c, ldc, v, incv, zero, work, 1 )

!           C := C - w * v'

    CALL dger( m, n, -tau, work, 1, v, incv, c, ldc )
  END IF
END IF

!DEALLOCATE( work )
RETURN

!     End of DLARF

END SUBROUTINE dlarf



SUBROUTINE dlarfb( side, trans, DIRECT, storev, m, n, k, v, ldv,  &
                   t, ldt, c, ldc )

!  -- LAPACK auxiliary routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997
! N.B. Arguments WORK and LDWORK have been removed.

!     .. Scalar Arguments ..
CHARACTER (LEN=1), INTENT(IN) :: side, trans, direct, storev
INTEGER, INTENT(IN)           :: k, ldc, ldt, ldv, m, n
!     ..
!     .. Array Arguments ..
REAL (wp), INTENT(IN)         :: t(:,:), v(:,:)
REAL (wp), INTENT(IN OUT)     :: c(:,:)
!     ..

!  Purpose
!  =======

!  DLARFB applies a real block reflector H or its transpose H' to a
!  real m by n matrix C, from either the left or the right.

!  Arguments
!  =========

!  SIDE    (input) CHARACTER*1
!          = 'L': apply H or H' from the Left
!          = 'R': apply H or H' from the Right

!  TRANS   (input) CHARACTER*1
!          = 'N': apply H (No transpose)
!          = 'T': apply H' (Transpose)

!  DIRECT  (input) CHARACTER*1
!          Indicates how H is formed from a product of elementary
!          reflectors
!          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!          = 'B': H = H(k) . . . H(2) H(1) (Backward)

!  STOREV  (input) CHARACTER*1
!          Indicates how the vectors which define the elementary
!          reflectors are stored:
!          = 'C': Columnwise
!          = 'R': Rowwise

!  M       (input) INTEGER
!          The number of rows of the matrix C.

!  N       (input) INTEGER
!          The number of columns of the matrix C.

!  K       (input) INTEGER
!          The order of the matrix T (= the number of elementary
!          reflectors whose product defines the block reflector).

!  V       (input) DOUBLE PRECISION array, dimension
!                                (LDV,K) if STOREV = 'C'
!                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
!                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
!          The matrix V. See further details.

!  LDV     (input) INTEGER
!          The leading dimension of the array V.
!          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
!          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
!          if STOREV = 'R', LDV >= K.

!  T       (input) DOUBLE PRECISION array, dimension (LDT,K)
!          The triangular k by k matrix T in the representation of the
!          block reflector.

!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= K.

!  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
!          On entry, the m by n matrix C.
!          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.

!  LDC     (input) INTEGER
!          The leading dimension of the array C. LDA >= max(1,M).

!  =====================================================================

!     .. Local Scalars ..
INTEGER                :: j, ldwork
CHARACTER (LEN=1)      :: transt
REAL (wp), ALLOCATABLE :: work(:,:)
!     ..
!     .. External Functions ..
! LOGICAL ::            lsame
! EXTERNAL           lsame
!     ..
!     .. External Subroutines ..
! EXTERNAL           dgemm, dtrmm
!     ..
!     .. Executable Statements ..

!     Quick return if possible

IF( m <= 0 .OR. n <= 0 ) RETURN

!     Allocate workspace
IF( lsame( side, 'L') ) THEN
  ldwork = MAX(1, n)
ELSE
  ldwork = MAX(1, m)
END IF
ALLOCATE( work(ldwork, k) )

IF( lsame( trans, 'N' ) ) THEN
  transt = 'T'
ELSE
  transt = 'N'
END IF

IF( lsame( storev, 'C' ) ) THEN

  IF( lsame( DIRECT, 'F' ) ) THEN

!           Let  V =  ( V1 )    (first K rows)
!                     ( V2 )
!           where  V1  is unit lower triangular.

    IF( lsame( side, 'L' ) ) THEN

!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )

!              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)

!              W := C1'

      DO j = 1, k
        work(1:n, j) = c(j, 1:n)
      END DO

!              W := W * V1

      CALL dtrmm( 'Right', 'Lower', 'No transpose', 'Unit', n,  &
                  k, one, v, ldv, work, ldwork )
      IF( m > k ) THEN

!                 W := W + C2'*V2

        CALL dgemm( 'Transpose', 'No transpose', n, k, m-k, one,  &
                    c( k+1:, : ), ldc, v( k+1:, : ), ldv, one, work, ldwork )
      END IF

!              W := W * T'  or  W * T

      CALL dtrmm( 'Right', 'Upper', transt, 'Non-unit', n, k,  &
                  one, t, ldt, work, ldwork )

!              C := C - V * W'

      IF( m > k ) THEN

!                 C2 := C2 - V2 * W'

        CALL dgemm( 'No transpose', 'Transpose', m-k, n, k, -one,  &
                    v( k+1:, : ), ldv, work, ldwork, one, c( k+1:, : ), ldc )
      END IF

!              W := W * V1'

      CALL dtrmm( 'Right', 'Lower', 'Transpose', 'Unit', n, k,  &
                  one, v, ldv, work, ldwork )

!              C1 := C1 - W'

      DO j = 1, k
        c( j, 1:n ) = c( j, 1:n ) - work( 1:n, j )
      END DO

    ELSE IF( lsame( side, 'R' ) ) THEN

!              Form  C * H  or  C * H'  where  C = ( C1  C2 )

!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)

!              W := C1

      work(1:m, 1:k) = c(1:m, 1:k)

!              W := W * V1

      CALL dtrmm( 'Right', 'Lower', 'No transpose', 'Unit', m,  &
                  k, one, v, ldv, work, ldwork )
      IF( n > k ) THEN

!                 W := W + C2 * V2

        CALL dgemm( 'No transpose', 'No transpose', m, k, n-k, one,  &
                    c( :, k+1: ), ldc, v( k+1:, : ), ldv, one, work, ldwork )
      END IF

!              W := W * T  or  W * T'

      CALL dtrmm( 'Right', 'Upper', trans, 'Non-unit', m, k,  &
                  one, t, ldt, work, ldwork )

!              C := C - W * V'

      IF( n > k ) THEN

!                 C2 := C2 - W * V2'

        CALL dgemm( 'No transpose', 'Transpose', m, n-k, k, -one,  &
                    work, ldwork, v( k+1:, : ), ldv, one, c( :, k+1: ), ldc )
      END IF

!              W := W * V1'

      CALL dtrmm( 'Right', 'Lower', 'Transpose', 'Unit', m, k,  &
                  one, v, ldv, work, ldwork )

!              C1 := C1 - W

      DO j = 1, k
        c( 1:m, j ) = c( 1:m, j ) - work( 1:m, j )
      END DO
    END IF

  ELSE

!           Let  V =  ( V1 )
!                     ( V2 )    (last K rows)
!           where  V2  is unit upper triangular.

    IF( lsame( side, 'L' ) ) THEN

!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )

!              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)

!              W := C2'

      DO j = 1, k
        work(1:n, j) = c(m-k+j, 1:n)
      END DO

!              W := W * V2

      CALL dtrmm( 'Right', 'Upper', 'No transpose', 'Unit', n,  &
                  k, one, v( m-k+1:, : ), ldv, work, ldwork )
      IF( m > k ) THEN

!                 W := W + C1'*V1

        CALL dgemm( 'Transpose', 'No transpose', n, k, m-k,  &
                    one, c, ldc, v, ldv, one, work, ldwork )
      END IF

!              W := W * T'  or  W * T

      CALL dtrmm( 'Right', 'Lower', transt, 'Non-unit', n, k,  &
                  one, t, ldt, work, ldwork )

!              C := C - V * W'

      IF( m > k ) THEN

!                 C1 := C1 - V1 * W'

        CALL dgemm( 'No transpose', 'Transpose', m-k, n, k, -one, v,  &
                    ldv, work, ldwork, one, c, ldc )
      END IF

!              W := W * V2'

      CALL dtrmm( 'Right', 'Upper', 'Transpose', 'Unit', n, k,  &
                  one, v( m-k+1:, : ), ldv, work, ldwork )

!              C2 := C2 - W'

      DO j = 1, k
        c( m-k+j, 1:n ) = c( m-k+j, 1:n ) - work( 1:n, j )
      END DO

    ELSE IF( lsame( side, 'R' ) ) THEN

!              Form  C * H  or  C * H'  where  C = ( C1  C2 )

!              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)

!              W := C2

      work(1:m, 1:k) = c(1:m, n-k+1:n)

!              W := W * V2

      CALL dtrmm( 'Right', 'Upper', 'No transpose', 'Unit', m,  &
                  k, one, v( n-k+1:, : ), ldv, work, ldwork )
      IF( n > k ) THEN

!                 W := W + C1 * V1

        CALL dgemm( 'No transpose', 'No transpose', m, k, n-k, one, c,  &
                    ldc, v, ldv, one, work, ldwork )
      END IF

!              W := W * T  or  W * T'

      CALL dtrmm( 'Right', 'Lower', trans, 'Non-unit', m, k,  &
                  one, t, ldt, work, ldwork )

!              C := C - W * V'

      IF( n > k ) THEN

!                 C1 := C1 - W * V1'

        CALL dgemm( 'No transpose', 'Transpose', m, n-k, k,  &
                    -one, work, ldwork, v, ldv, one, c, ldc )
      END IF

!              W := W * V2'

      CALL dtrmm( 'Right', 'Upper', 'Transpose', 'Unit', m, k,  &
                  one, v( n-k+1:, : ), ldv, work, ldwork )

!              C2 := C2 - W

      DO j = 1, k
        c( 1:m, n-k+j ) = c( 1:m, n-k+j ) - work( 1:m, j )
      END DO
    END IF
  END IF

ELSE IF( lsame( storev, 'R' ) ) THEN

  IF( lsame( DIRECT, 'F' ) ) THEN

!           Let  V =  ( V1  V2 )    (V1: first K columns)
!           where  V1  is unit upper triangular.

    IF( lsame( side, 'L' ) ) THEN

!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )

!              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)

!              W := C1'

      DO j = 1, k
        work(1:n, j) = c(j, 1:n)
      END DO

!              W := W * V1'

      CALL dtrmm( 'Right', 'Upper', 'Transpose', 'Unit', n, k,  &
                  one, v, ldv, work, ldwork )
      IF( m > k ) THEN

!                 W := W + C2'*V2'

        CALL dgemm( 'Transpose', 'Transpose', n, k, m-k, one,  &
                    c( k+1:, : ), ldc, v( :, k+1: ), ldv, one, work, ldwork )
      END IF

!              W := W * T'  or  W * T

      CALL dtrmm( 'Right', 'Upper', transt, 'Non-unit', n, k,  &
                  one, t, ldt, work, ldwork )

!              C := C - V' * W'

      IF( m > k ) THEN

!                 C2 := C2 - V2' * W'

        CALL dgemm( 'Transpose', 'Transpose', m-k, n, k, -one,  &
                    v( :, k+1: ), ldv, work, ldwork, one, c( k+1:, : ), ldc )
      END IF

!              W := W * V1

      CALL dtrmm( 'Right', 'Upper', 'No transpose', 'Unit', n,  &
                  k, one, v, ldv, work, ldwork )

!              C1 := C1 - W'

      DO j = 1, k
        c( j, 1:n ) = c( j, 1:n ) - work( 1:n, j )
      END DO

    ELSE IF( lsame( side, 'R' ) ) THEN

!              Form  C * H  or  C * H'  where  C = ( C1  C2 )

!              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)

!              W := C1

      work(1:m, 1:k) = c(1:m, 1:k)

!              W := W * V1'

      CALL dtrmm( 'Right', 'Upper', 'Transpose', 'Unit', m, k,  &
                  one, v, ldv, work, ldwork )
      IF( n > k ) THEN

!                 W := W + C2 * V2'

        CALL dgemm( 'No transpose', 'Transpose', m, k, n-k, one,  &
                    c( :, k+1: ), ldc, v( :, k+1: ), ldv, one, work, ldwork )
      END IF

!              W := W * T  or  W * T'

      CALL dtrmm( 'Right', 'Upper', trans, 'Non-unit', m, k,  &
                  one, t, ldt, work, ldwork )

!              C := C - W * V

      IF( n > k ) THEN

!                 C2 := C2 - W * V2

        CALL dgemm( 'No transpose', 'No transpose', m, n-k, k, -one,  &
                    work, ldwork, v( :, k+1: ), ldv, one, c( :, k+1: ), ldc )
      END IF

!              W := W * V1

      CALL dtrmm( 'Right', 'Upper', 'No transpose', 'Unit', m,  &
                  k, one, v, ldv, work, ldwork )

!              C1 := C1 - W

      c( 1:m, 1:k ) = c( 1:m, 1:k ) - work( 1:m, 1:k )

    END IF

  ELSE

!           Let  V =  ( V1  V2 )    (V2: last K columns)
!           where  V2  is unit lower triangular.

    IF( lsame( side, 'L' ) ) THEN

!              Form  H * C  or  H' * C  where  C = ( C1 )
!                                                  ( C2 )

!              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)

!              W := C2'

      DO j = 1, k
        work(1:n, j) = c(m-k+j, 1:n)
      END DO

!              W := W * V2'

      CALL dtrmm( 'Right', 'Lower', 'Transpose', 'Unit', n, k,  &
                  one, v( :, m-k+1: ), ldv, work, ldwork )
      IF( m > k ) THEN

!                 W := W + C1'*V1'

        CALL dgemm( 'Transpose', 'Transpose', n, k, m-k, one,  &
                    c, ldc, v, ldv, one, work, ldwork )
      END IF

!              W := W * T'  or  W * T

      CALL dtrmm( 'Right', 'Lower', transt, 'Non-unit', n, k,  &
                  one, t, ldt, work, ldwork )

!              C := C - V' * W'

      IF( m > k ) THEN

!                 C1 := C1 - V1' * W'

        CALL dgemm( 'Transpose', 'Transpose', m-k, n, k, -one,  &
                    v, ldv, work, ldwork, one, c, ldc )
      END IF

!              W := W * V2

      CALL dtrmm( 'Right', 'Lower', 'No transpose', 'Unit', n,  &
                  k, one, v( :, m-k+1: ), ldv, work, ldwork )

!              C2 := C2 - W'

      DO j = 1, k
        c( m-k+j, 1:n ) = c( m-k+j, 1:n ) - work( 1:n, j )
      END DO

    ELSE IF( lsame( side, 'R' ) ) THEN

!              Form  C * H  or  C * H'  where  C = ( C1  C2 )

!              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)

!              W := C2

      work(1:m, 1:k) = c(1:m, n-k+1:n)

!              W := W * V2'

      CALL dtrmm( 'Right', 'Lower', 'Transpose', 'Unit', m, k,  &
                  one, v( :, n-k+1: ), ldv, work, ldwork )
      IF( n > k ) THEN

!                 W := W + C1 * V1'

        CALL dgemm( 'No transpose', 'Transpose', m, k, n-k,  &
                    one, c, ldc, v, ldv, one, work, ldwork )
      END IF

!              W := W * T  or  W * T'

      CALL dtrmm( 'Right', 'Lower', trans, 'Non-unit', m, k,  &
                  one, t, ldt, work, ldwork )

!              C := C - W * V

      IF( n > k ) THEN

!                 C1 := C1 - W * V1

        CALL dgemm( 'No transpose', 'No transpose', m, n-k, k,  &
                    -one, work, ldwork, v, ldv, one, c, ldc )
      END IF

!              W := W * V2

      CALL dtrmm( 'Right', 'Lower', 'No transpose', 'Unit', m,  &
                  k, one, v( :, n-k+1: ), ldv, work, ldwork )

!              C1 := C1 - W

      DO j = 1, k
        c( 1:m, n-k+j ) = c( 1:m, n-k+j ) - work( 1:m, j )
      END DO

    END IF

  END IF
END IF

DEALLOCATE( work )
RETURN

!     End of DLARFB

END SUBROUTINE dlarfb



SUBROUTINE dlarft( DIRECT, storev, n, k, v, ldv, tau, t, ldt )

!  -- LAPACK auxiliary routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997

!     .. Scalar Arguments ..
CHARACTER (LEN=1), INTENT(IN) :: direct, storev
INTEGER, INTENT(IN)           :: k, ldt, ldv, n
!     ..
!     .. Array Arguments ..
REAL (wp), INTENT(IN)     :: tau(:)
REAL (wp), INTENT(IN OUT) :: v(:,:)
REAL (wp), INTENT(OUT)    :: t(:,:)
!     ..

!  Purpose
!  =======

!  DLARFT forms the triangular factor T of a real block reflector H
!  of order n, which is defined as a product of k elementary reflectors.

!  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;

!  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.

!  If STOREV = 'C', the vector which defines the elementary reflector
!  H(i) is stored in the i-th column of the array V, and

!     H  =  I - V * T * V'

!  If STOREV = 'R', the vector which defines the elementary reflector
!  H(i) is stored in the i-th row of the array V, and

!     H  =  I - V' * T * V

!  Arguments
!  =========

!  DIRECT  (input) CHARACTER*1
!          Specifies the order in which the elementary reflectors are
!          multiplied to form the block reflector:
!          = 'F': H = H(1) H(2) . . . H(k) (Forward)
!          = 'B': H = H(k) . . . H(2) H(1) (Backward)

!  STOREV  (input) CHARACTER*1
!          Specifies how the vectors which define the elementary
!          reflectors are stored (see also Further Details):
!          = 'C': columnwise
!          = 'R': rowwise

!  N       (input) INTEGER
!          The order of the block reflector H. N >= 0.

!  K       (input) INTEGER
!          The order of the triangular factor T (= the number of
!          elementary reflectors). K >= 1.

!  V       (input/output) DOUBLE PRECISION array, dimension
!                               (LDV,K) if STOREV = 'C'
!                               (LDV,N) if STOREV = 'R'
!          The matrix V. See further details.

!  LDV     (input) INTEGER
!          The leading dimension of the array V.
!          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.

!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i).

!  T       (output) DOUBLE PRECISION array, dimension (LDT,K)
!          The k by k triangular factor T of the block reflector.
!          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
!          lower triangular. The rest of the array is not used.

!  LDT     (input) INTEGER
!          The leading dimension of the array T. LDT >= K.

!  Further Details
!  ===============

!  The shape of the matrix V and the storage of the vectors which define
!  the H(i) is best illustrated by the following example with n = 5 and
!  k = 3. The elements equal to 1 are not stored; the corresponding
!  array elements are modified but restored on exit. The rest of the
!  array is not used.

!  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':

!               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
!                   ( v1  1    )                     (     1 v2 v2 v2 )
!                   ( v1 v2  1 )                     (        1 v3 v3 )
!                   ( v1 v2 v3 )
!                   ( v1 v2 v3 )

!  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':

!               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
!                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
!                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
!                   (     1 v3 )
!                   (        1 )

!  =====================================================================

!     .. Local Scalars ..
INTEGER   :: i
REAL (wp) :: vii
!     ..
!     .. External Subroutines ..
! EXTERNAL           dgemv, dtrmv
!     ..
!     .. External Functions ..
! LOGICAL ::            lsame
! EXTERNAL           lsame
!     ..
!     .. Executable Statements ..

!     Quick return if possible

IF( n == 0 ) RETURN

IF( lsame( DIRECT, 'F' ) ) THEN
  DO i = 1, k
    IF( tau( i ) == zero ) THEN

!              H(i)  =  I

      t( 1:i, i ) = zero
    ELSE

!              general case

      vii = v( i, i )
      v( i, i ) = one
      IF( lsame( storev, 'C' ) ) THEN

!                 T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' * V(i:n,i)

        CALL dgemv( 'Transpose', n-i+1, i-1, -tau( i ), v( i:, : ), ldv,  &
                    v( i:, i ), 1, zero, t( :, i ), 1 )
      ELSE

!                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) * V(i,i:n)'

        CALL dgemv( 'No transpose', i-1, n-i+1, -tau( i ), v( :, i: ), ldv, &
                    v( i, i: ), 1, zero, t( :, i ), 1 )
      END IF
      v( i, i ) = vii

!              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)

      CALL dtrmv( 'Upper', 'No transpose', 'Non-unit', i-1, t,  &
                  ldt, t( :, i ), 1 )
      t( i, i ) = tau( i )
    END IF
  END DO
ELSE
  DO i = k, 1, -1
    IF( tau( i ) == zero ) THEN

!              H(i)  =  I

      t( i:k, i ) = zero
    ELSE

!              general case

      IF( i < k ) THEN
        IF( lsame( storev, 'C' ) ) THEN
          vii = v( n-k+i, i )
          v( n-k+i, i ) = one

!                    T(i+1:k,i) :=
!                            - tau(i) * V(1:n-k+i,i+1:k)' * V(1:n-k+i,i)

          CALL dgemv( 'Transpose', n-k+i, k-i, -tau( i ), v( :, i+1: ),  &
                      ldv, v( :, i ), 1, zero, t( i+1:, i ), 1 )
          v( n-k+i, i ) = vii
        ELSE
          vii = v( i, n-k+i )
          v( i, n-k+i ) = one

!                    T(i+1:k,i) :=
!                            - tau(i) * V(i+1:k,1:n-k+i) * V(i,1:n-k+i)'

          CALL dgemv( 'No transpose', k-i, n-k+i, -tau( i ), v( i+1:, : ),  &
                      ldv, v( i, : ), 1, zero, t( i+1:, i ), 1 )
          v( i, n-k+i ) = vii
        END IF

!                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)

        CALL dtrmv( 'Lower', 'No transpose', 'Non-unit', k-i,  &
                    t( i+1:, i+1: ), ldt, t( i+1:, i ), 1 )
      END IF
      t( i, i ) = tau( i )
    END IF
  END DO
END IF
RETURN

!     End of DLARFT

END SUBROUTINE dlarft



SUBROUTINE dlasr( side, pivot, DIRECT, m, n, c, s, a, lda )

!  -- LAPACK auxiliary routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997

!     .. Scalar Arguments ..
CHARACTER (LEN=1), INTENT(IN) :: side, pivot, direct
INTEGER, INTENT(IN)           :: lda, m, n
!     ..
!     .. Array Arguments ..
REAL (wp), INTENT(IN)         :: c(:), s(:)
REAL (wp), INTENT(IN OUT)     :: a(:,:)
!     ..

!  Purpose
!  =======

!  DLASR   performs the transformation

!     A := P*A,   when SIDE = 'L' or 'l'  (  Left-hand side )

!     A := A*P',  when SIDE = 'R' or 'r'  ( Right-hand side )

!  where A is an m by n real matrix and P is an orthogonal matrix,
!  consisting of a sequence of plane rotations determined by the
!  parameters PIVOT and DIRECT as follows ( z = m when SIDE = 'L' or 'l'
!  and z = n when SIDE = 'R' or 'r' ):

!  When  DIRECT = 'F' or 'f'  ( Forward sequence ) then

!     P = P( z - 1 )*...*P( 2 )*P( 1 ),

!  and when DIRECT = 'B' or 'b'  ( Backward sequence ) then

!     P = P( 1 )*P( 2 )*...*P( z - 1 ),

!  where  P( k ) is a plane rotation matrix for the following planes:

!     when  PIVOT = 'V' or 'v'  ( Variable pivot ),
!        the plane ( k, k + 1 )

!     when  PIVOT = 'T' or 't'  ( Top pivot ),
!        the plane ( 1, k + 1 )

!     when  PIVOT = 'B' or 'b'  ( Bottom pivot ),
!        the plane ( k, z )

!  c( k ) and s( k )  must contain the  cosine and sine that define the
!  matrix  P( k ).  The two by two plane rotation part of the matrix
!  P( k ), R( k ), is assumed to be of the form

!     R( k ) = (  c( k )  s( k ) ).
!              ( -s( k )  c( k ) )

!  This version vectorises across rows of the array A when SIDE = 'L'.

!  Arguments
!  =========

!  SIDE    (input) CHARACTER*1
!          Specifies whether the plane rotation matrix P is applied to
!          A on the left or the right.
!          = 'L':  Left, compute A := P*A
!          = 'R':  Right, compute A:= A*P'

!  DIRECT  (input) CHARACTER*1
!          Specifies whether P is a forward or backward sequence of
!          plane rotations.
!          = 'F':  Forward, P = P( z - 1 )*...*P( 2 )*P( 1 )
!          = 'B':  Backward, P = P( 1 )*P( 2 )*...*P( z - 1 )

!  PIVOT   (input) CHARACTER*1
!          Specifies the plane for which P(k) is a plane rotation
!          matrix.
!          = 'V':  Variable pivot, the plane (k,k+1)
!          = 'T':  Top pivot, the plane (1,k+1)
!          = 'B':  Bottom pivot, the plane (k,z)

!  M       (input) INTEGER
!          The number of rows of the matrix A.  If m <= 1, an immediate
!          return is effected.

!  N       (input) INTEGER
!          The number of columns of the matrix A.  If n <= 1, an
!          immediate return is effected.

!  C, S    (input) DOUBLE PRECISION arrays, dimension
!                  (M-1) if SIDE = 'L'
!                  (N-1) if SIDE = 'R'
!          c(k) and s(k) contain the cosine and sine that define the
!          matrix P(k).  The two by two plane rotation part of the
!          matrix P(k), R(k), is assumed to be of the form
!          R( k ) = (  c( k )  s( k ) ).
!                   ( -s( k )  c( k ) )

!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          The m by n matrix A.  On exit, A is overwritten by P*A if
!          SIDE = 'R' or by A*P' if SIDE = 'L'.

!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).

!     .. Local Scalars ..
INTEGER ::            i, info, j
REAL (wp) ::   ctemp, stemp, temp
!     ..
!     .. External Functions ..
! LOGICAL ::            lsame
! EXTERNAL           lsame
!     ..
!     .. External Subroutines ..
! EXTERNAL           xerbla
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC          MAX
!     ..
!     .. Executable Statements ..

!     Test the input parameters

info = 0
IF( .NOT.( lsame( side, 'L' ) .OR. lsame( side, 'R' ) ) ) THEN
  info = 1
ELSE IF( .NOT.( lsame( pivot, 'V' ) .OR. lsame( pivot,  &
  'T' ) .OR. lsame( pivot, 'B' ) ) ) THEN
  info = 2
ELSE IF( .NOT.( lsame( DIRECT, 'F' ) .OR. lsame( DIRECT, 'B' ) ) )  &
  THEN
  info = 3
ELSE IF( m < 0 ) THEN
  info = 4
ELSE IF( n < 0 ) THEN
  info = 5
ELSE IF( lda < MAX( 1, m ) ) THEN
  info = 9
END IF
IF( info /= 0 ) THEN
  CALL erinfo(info, 'DLASR ')
  RETURN
END IF

!     Quick return if possible

IF( ( m == 0 ) .OR. ( n == 0 ) ) RETURN
IF( lsame( side, 'L' ) ) THEN

!        Form  P * A

  IF( lsame( pivot, 'V' ) ) THEN
    IF( lsame( DIRECT, 'F' ) ) THEN
      DO j = 1, m - 1
        ctemp = c( j )
        stemp = s( j )
        IF( ( ctemp /= one ) .OR. ( stemp /= zero ) ) THEN
          DO i = 1, n
            temp = a( j+1, i )
            a( j+1, i ) = ctemp*temp - stemp*a( j, i )
            a( j, i ) = stemp*temp + ctemp*a( j, i )
          END DO
        END IF
      END DO
    ELSE IF( lsame( DIRECT, 'B' ) ) THEN
      DO j = m - 1, 1, -1
        ctemp = c( j )
        stemp = s( j )
        IF( ( ctemp /= one ) .OR. ( stemp /= zero ) ) THEN
          DO i = 1, n
            temp = a( j+1, i )
            a( j+1, i ) = ctemp*temp - stemp*a( j, i )
            a( j, i ) = stemp*temp + ctemp*a( j, i )
          END DO
        END IF
      END DO
    END IF
  ELSE IF( lsame( pivot, 'T' ) ) THEN
    IF( lsame( DIRECT, 'F' ) ) THEN
      DO j = 2, m
        ctemp = c( j-1 )
        stemp = s( j-1 )
        IF( ( ctemp /= one ) .OR. ( stemp /= zero ) ) THEN
          DO i = 1, n
            temp = a( j, i )
            a( j, i ) = ctemp*temp - stemp*a( 1, i )
            a( 1, i ) = stemp*temp + ctemp*a( 1, i )
          END DO
        END IF
      END DO
    ELSE IF( lsame( DIRECT, 'B' ) ) THEN
      DO j = m, 2, -1
        ctemp = c( j-1 )
        stemp = s( j-1 )
        IF( ( ctemp /= one ) .OR. ( stemp /= zero ) ) THEN
          DO i = 1, n
            temp = a( j, i )
            a( j, i ) = ctemp*temp - stemp*a( 1, i )
            a( 1, i ) = stemp*temp + ctemp*a( 1, i )
          END DO
        END IF
      END DO
    END IF
  ELSE IF( lsame( pivot, 'B' ) ) THEN
    IF( lsame( DIRECT, 'F' ) ) THEN
      DO j = 1, m - 1
        ctemp = c( j )
        stemp = s( j )
        IF( ( ctemp /= one ) .OR. ( stemp /= zero ) ) THEN
          DO i = 1, n
            temp = a( j, i )
            a( j, i ) = stemp*a( m, i ) + ctemp*temp
            a( m, i ) = ctemp*a( m, i ) - stemp*temp
          END DO
        END IF
      END DO
    ELSE IF( lsame( DIRECT, 'B' ) ) THEN
      DO j = m - 1, 1, -1
        ctemp = c( j )
        stemp = s( j )
        IF( ( ctemp /= one ) .OR. ( stemp /= zero ) ) THEN
          DO i = 1, n
            temp = a( j, i )
            a( j, i ) = stemp*a( m, i ) + ctemp*temp
            a( m, i ) = ctemp*a( m, i ) - stemp*temp
          END DO
        END IF
      END DO
    END IF
  END IF
ELSE IF( lsame( side, 'R' ) ) THEN

!        Form A * P'

  IF( lsame( pivot, 'V' ) ) THEN
    IF( lsame( DIRECT, 'F' ) ) THEN
      DO j = 1, n - 1
        ctemp = c( j )
        stemp = s( j )
        IF( ( ctemp /= one ) .OR. ( stemp /= zero ) ) THEN
          DO i = 1, m
            temp = a( i, j+1 )
            a( i, j+1 ) = ctemp*temp - stemp*a( i, j )
            a( i, j ) = stemp*temp + ctemp*a( i, j )
          END DO
        END IF
      END DO
    ELSE IF( lsame( DIRECT, 'B' ) ) THEN
      DO j = n - 1, 1, -1
        ctemp = c( j )
        stemp = s( j )
        IF( ( ctemp /= one ) .OR. ( stemp /= zero ) ) THEN
          DO i = 1, m
            temp = a( i, j+1 )
            a( i, j+1 ) = ctemp*temp - stemp*a( i, j )
            a( i, j ) = stemp*temp + ctemp*a( i, j )
          END DO
        END IF
      END DO
    END IF
  ELSE IF( lsame( pivot, 'T' ) ) THEN
    IF( lsame( DIRECT, 'F' ) ) THEN
      DO j = 2, n
        ctemp = c( j-1 )
        stemp = s( j-1 )
        IF( ( ctemp /= one ) .OR. ( stemp /= zero ) ) THEN
          DO i = 1, m
            temp = a( i, j )
            a( i, j ) = ctemp*temp - stemp*a( i, 1 )
            a( i, 1 ) = stemp*temp + ctemp*a( i, 1 )
          END DO
        END IF
      END DO
    ELSE IF( lsame( DIRECT, 'B' ) ) THEN
      DO j = n, 2, -1
        ctemp = c( j-1 )
        stemp = s( j-1 )
        IF( ( ctemp /= one ) .OR. ( stemp /= zero ) ) THEN
          DO i = 1, m
            temp = a( i, j )
            a( i, j ) = ctemp*temp - stemp*a( i, 1 )
            a( i, 1 ) = stemp*temp + ctemp*a( i, 1 )
          END DO
        END IF
      END DO
    END IF
  ELSE IF( lsame( pivot, 'B' ) ) THEN
    IF( lsame( DIRECT, 'F' ) ) THEN
      DO j = 1, n - 1
        ctemp = c( j )
        stemp = s( j )
        IF( ( ctemp /= one ) .OR. ( stemp /= zero ) ) THEN
          DO i = 1, m
            temp = a( i, j )
            a( i, j ) = stemp*a( i, n ) + ctemp*temp
            a( i, n ) = ctemp*a( i, n ) - stemp*temp
          END DO
        END IF
      END DO
    ELSE IF( lsame( DIRECT, 'B' ) ) THEN
      DO j = n - 1, 1, -1
        ctemp = c( j )
        stemp = s( j )
        IF( ( ctemp /= one ) .OR. ( stemp /= zero ) ) THEN
          DO i = 1, m
            temp = a( i, j )
            a( i, j ) = stemp*a( i, n ) + ctemp*temp
            a( i, n ) = ctemp*a( i, n ) - stemp*temp
          END DO
        END IF
      END DO
    END IF
  END IF
END IF

RETURN

!     End of DLASR

END SUBROUTINE dlasr




SUBROUTINE dlaswp( n, a, k1, k2, ipiv, incx )

!  -- LAPACK auxiliary routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997
! Argument LDA has been removed.

!     .. Scalar Arguments ..
INTEGER, INTENT(IN)       :: incx, k1, k2, n
!     ..
!     .. Array Arguments ..
INTEGER, INTENT(IN)       :: ipiv(:)
REAL (wp), INTENT(IN OUT) :: a(:,:)
!     ..

!  Purpose
!  =======

!  DLASWP performs a series of row interchanges on the matrix A.
!  One row interchange is initiated for each of rows K1 through K2 of A.

!  Arguments
!  =========

!  N       (input) INTEGER
!          The number of columns of the matrix A.

!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the matrix of column dimension N to which the row
!          interchanges will be applied.
!          On exit, the permuted matrix.

!  K1      (input) INTEGER
!          The first element of IPIV for which a row interchange will
!          be done.

!  K2      (input) INTEGER
!          The last element of IPIV for which a row interchange will
!          be done.

!  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
!          The vector of pivot indices.  Only the elements in positions
!          K1 through K2 of IPIV are accessed.
!          IPIV(K) = L implies rows K and L are to be interchanged.

!  INCX    (input) INTEGER
!          The increment between successive values of IPIV.  If IPIV
!          is negative, the pivots are applied in reverse order.


!     .. Local Scalars ..
INTEGER ::            i, ip, ix
!     ..
!     .. External Subroutines ..
! EXTERNAL           dswap
!     ..
!     .. Executable Statements ..

!     Interchange row I with row IPIV(I) for each of rows K1 through K2.

IF( incx == 0 ) RETURN
IF( incx > 0 ) THEN
  ix = k1
ELSE
  ix = 1 + ( 1-k2 )*incx
END IF
IF( incx == 1 ) THEN
  DO i = k1, k2
    ip = ipiv( i )
    IF( ip /= i ) CALL dswap( n, a( i, : ), 1, a( ip, : ), 1 )
  END DO
ELSE IF( incx > 1 ) THEN
  DO i = k1, k2
    ip = ipiv( ix )
    IF( ip /= i ) CALL dswap( n, a( i, : ), 1, a( ip, : ), 1 )
    ix = ix + incx
  END DO
ELSE IF( incx < 0 ) THEN
  DO i = k2, k1, -1
    ip = ipiv( ix )
    IF( ip /= i ) CALL dswap( n, a( i, : ), 1, a( ip, : ), 1 )
    ix = ix + incx
  END DO
END IF

RETURN

!     End of DLASWP

END SUBROUTINE dlaswp


SUBROUTINE dlatrd( uplo, n, nb, a, lda, e, tau, w, ldw )

!  -- LAPACK auxiliary routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997

!     .. Scalar Arguments ..
CHARACTER (LEN=1), INTENT(IN) :: uplo
INTEGER, INTENT(IN)           :: lda, ldw, n, nb
!     ..
!     .. Array Arguments ..
REAL (wp), INTENT(IN OUT) :: a(:,:)
REAL (wp), INTENT(OUT)    :: e(:), tau(:), w(:,:)
!     ..

!  Purpose
!  =======

!  DLATRD reduces NB rows and columns of a real symmetric matrix A to
!  symmetric tridiagonal form by an orthogonal similarity
!  transformation Q' * A * Q, and returns the matrices V and W which are
!  needed to apply the transformation to the unreduced part of A.

!  If UPLO = 'U', DLATRD reduces the last NB rows and columns of a
!  matrix, of which the upper triangle is supplied;
!  if UPLO = 'L', DLATRD reduces the first NB rows and columns of a
!  matrix, of which the lower triangle is supplied.

!  This is an auxiliary routine called by DSYTRD.

!  Arguments
!  =========

!  UPLO    (input) CHARACTER
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored:
!          = 'U': Upper triangular
!          = 'L': Lower triangular

!  N       (input) INTEGER
!          The order of the matrix A.

!  NB      (input) INTEGER
!          The number of rows and columns to be reduced.

!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!          n-by-n upper triangular part of A contains the upper
!          triangular part of the matrix A, and the strictly lower
!          triangular part of A is not referenced.  If UPLO = 'L', the
!          leading n-by-n lower triangular part of A contains the lower
!          triangular part of the matrix A, and the strictly upper
!          triangular part of A is not referenced.
!          On exit:
!          if UPLO = 'U', the last NB columns have been reduced to
!            tridiagonal form, with the diagonal elements overwriting
!            the diagonal elements of A; the elements above the diagonal
!            with the array TAU, represent the orthogonal matrix Q as a
!            product of elementary reflectors;
!          if UPLO = 'L', the first NB columns have been reduced to
!            tridiagonal form, with the diagonal elements overwriting
!            the diagonal elements of A; the elements below the diagonal
!            with the array TAU, represent the  orthogonal matrix Q as a
!            product of elementary reflectors.
!          See Further Details.

!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= (1,N).

!  E       (output) DOUBLE PRECISION array, dimension (N-1)
!          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal
!          elements of the last NB columns of the reduced matrix;
!          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of
!          the first NB columns of the reduced matrix.

!  TAU     (output) DOUBLE PRECISION array, dimension (N)
!          The scalar factors of the elementary reflectors, stored in
!          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.
!          See Further Details.

!  W       (output) DOUBLE PRECISION array, dimension (LDW,NB)
!          The n-by-nb matrix W required to update the unreduced part
!          of A.

!  LDW     (input) INTEGER
!          The leading dimension of the array W. LDW >= max(1,N).

!  Further Details
!  ===============

!  If UPLO = 'U', the matrix Q is represented as a product of elementary
!  reflectors

!     Q = H(n) H(n-1) . . . H(n-nb+1).

!  Each H(i) has the form

!     H(i) = I - tau * v * v'

!  where tau is a real scalar, and v is a real vector with
!  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),
!  and tau in TAU(i-1).

!  If UPLO = 'L', the matrix Q is represented as a product of elementary
!  reflectors

!     Q = H(1) H(2) . . . H(nb).

!  Each H(i) has the form

!     H(i) = I - tau * v * v'

!  where tau is a real scalar, and v is a real vector with
!  v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),
!  and tau in TAU(i).

!  The elements of the vectors v together form the n-by-nb matrix V
!  which is needed, with W, to apply the transformation to the unreduced
!  part of the matrix, using a symmetric rank-2k update of the form:
!  A := A - V*W' - W*V'.

!  The contents of A on exit are illustrated by the following examples
!  with n = 5 and nb = 2:

!  if UPLO = 'U':                       if UPLO = 'L':

!    (  a   a   a   v4  v5 )              (  d                  )
!    (      a   a   v4  v5 )              (  1   d              )
!    (          a   1   v5 )              (  v1  1   a          )
!    (              d   1  )              (  v1  v2  a   a      )
!    (                  d  )              (  v1  v2  a   a   a  )

!  where d denotes a diagonal element of the reduced matrix, a denotes
!  an element of the original matrix that is unchanged, and vi denotes
!  an element of the vector defining H(i).

!  =====================================================================

!     .. Parameters ..
REAL (wp), PARAMETER :: half = 0.5_wp
!     ..
!     .. Local Scalars ..
INTEGER   :: i, iw
REAL (wp) :: alpha
!     ..
!     .. External Subroutines ..
! EXTERNAL           daxpy, dgemv, dlarfg, dscal, dsymv
!     ..
!     .. External Functions ..
! LOGICAL ::            lsame
! DOUBLE PRECISION ::   ddot
! EXTERNAL           lsame, ddot
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC          MIN
!     ..
!     .. Executable Statements ..

!     Quick return if possible

IF( n <= 0 ) RETURN

IF( lsame( uplo, 'U' ) ) THEN

!        Reduce last NB columns of upper triangle

  DO i = n, n - nb + 1, -1
    iw = i - n + nb
    IF( i < n ) THEN

!              Update A(1:i,i)

      CALL dgemv( 'No transpose', i, n-i, -one, a( :, i+1: ),  &
                  lda, w( i, iw+1: ), 1, one, a( :, i ), 1 )
      CALL dgemv( 'No transpose', i, n-i, -one, w( :, iw+1: ),  &
                  ldw, a( i, i+1: ), 1, one, a( :, i ), 1 )
    END IF
    IF( i > 1 ) THEN

!              Generate elementary reflector H(i) to annihilate
!              A(1:i-2,i)

      CALL dlarfg( i-1, a( i-1, i ), a( :, i ), 1, tau( i-1 ) )
      e( i-1 ) = a( i-1, i )
      a( i-1, i ) = one

!              Compute W(1:i-1,i)

      CALL dsymv( 'Upper', i-1, one, a, lda, a( :, i ), 1, zero, w( :, iw ), &
                  1 )
      IF( i < n ) THEN
        CALL dgemv( 'Transpose', i-1, n-i, one, w( :, iw+1: ),  &
                    ldw, a( :, i ), 1, zero, w( i+1:, iw ), 1 )
        CALL dgemv( 'No transpose', i-1, n-i, -one, a( :, i+1: ),  &
                    lda, w( i+1:, iw ), 1, one, w( :, iw ), 1 )
        CALL dgemv( 'Transpose', i-1, n-i, one, a( :, i+1: ),  &
                    lda, a( :, i ), 1, zero, w( i+1:, iw ), 1 )
        CALL dgemv( 'No transpose', i-1, n-i, -one, w( :, iw+1: ),  &
                    ldw, w( i+1:, iw ), 1, one, w( :, iw ), 1 )
      END IF
      w(1:i-1, iw) = tau(i-1) * w(1:i-1, iw)
      alpha = -half*tau( i-1 ) * ddot( i-1, w( :, iw ), 1, a( :, i ), 1 )
      CALL daxpy( i-1, alpha, a( :, i ), 1, w( :, iw ), 1 )
    END IF

  END DO
ELSE

!        Reduce first NB columns of lower triangle

  DO i = 1, nb

!           Update A(i:n,i)

    CALL dgemv( 'No transpose', n-i+1, i-1, -one, a( i:, : ),  &
                lda, w( i, : ), 1, one, a( i:, i ), 1 )
    CALL dgemv( 'No transpose', n-i+1, i-1, -one, w( i:, : ),  &
                ldw, a( i, : ), 1, one, a( i:, i ), 1 )
    IF( i < n ) THEN

!              Generate elementary reflector H(i) to annihilate
!              A(i+2:n,i)

      CALL dlarfg( n-i, a( i+1, i ), a( MIN( i+2, n ):, i ), 1, tau( i ) )
      e( i ) = a( i+1, i )
      a( i+1, i ) = one

!              Compute W(i+1:n,i)

      CALL dsymv( 'Lower', n-i, one, a( i+1:, i+1: ), lda,  &
                  a( i+1:, i ), 1, zero, w( i+1:, i ), 1 )
      CALL dgemv( 'Transpose', n-i, i-1, one, w( i+1:, : ), ldw,  &
                  a( i+1:, i ), 1, zero, w( :, i ), 1 )
      CALL dgemv( 'No transpose', n-i, i-1, -one, a( i+1:, : ),  &
                  lda, w( :, i ), 1, one, w( i+1:, i ), 1 )
      CALL dgemv( 'Transpose', n-i, i-1, one, a( i+1:, : ), lda,  &
                  a( i+1:, i ), 1, zero, w( :, i ), 1 )
      CALL dgemv( 'No transpose', n-i, i-1, -one, w( i+1:, : ),  &
                  ldw, w( :, i ), 1, one, w( i+1:, i ), 1 )
      w(i+1:n, i) = tau(i) * w(i+1:n, i)
      alpha = -half*tau( i ) * ddot( n-i, w( i+1:, i ), 1, a( i+1:, i ), 1 )
      CALL daxpy( n-i, alpha, a( i+1:, i ), 1, w( i+1:, i ), 1 )
    END IF

  END DO
END IF

RETURN

!     End of DLATRD

END SUBROUTINE dlatrd



SUBROUTINE dlaset( uplo, m, n, alpha, beta, a)

!  -- LAPACK auxiliary routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997
! Last argument, LDA, has been deleted.

!     .. Scalar Arguments ..
CHARACTER (LEN=1), INTENT(IN) :: uplo
INTEGER, INTENT(IN)           :: m, n
REAL (wp), INTENT(IN)         :: alpha, beta
!     ..
!     .. Array Arguments ..
REAL (wp), INTENT(IN OUT)     :: a(:,:)
!     ..

!  Purpose
!  =======

!  DLASET initializes an m-by-n matrix A to BETA on the diagonal and
!  ALPHA on the offdiagonals.

!  Arguments
!  =========

!  UPLO    (input) CHARACTER*1
!          Specifies the part of the matrix A to be set.
!          = 'U':      Upper triangular part is set; the strictly lower
!                      triangular part of A is not changed.
!          = 'L':      Lower triangular part is set; the strictly upper
!                      triangular part of A is not changed.
!          Otherwise:  All of the matrix A is set.

!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.

!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.

!  ALPHA   (input) DOUBLE PRECISION
!          The constant to which the offdiagonal elements are to be set.

!  BETA    (input) DOUBLE PRECISION
!          The constant to which the diagonal elements are to be set.

!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On exit, the leading m-by-n submatrix of A is set as follows:

!          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
!          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
!          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,

!          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).

!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).


!     .. Local Scalars ..
INTEGER ::            i, j
!     ..
!     .. External Functions ..
! LOGICAL ::            lsame
! EXTERNAL           lsame
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC          MIN
!     ..
!     .. Executable Statements ..

IF( lsame( uplo, 'U' ) ) THEN

!        Set the strictly upper triangular or trapezoidal part of the
!        array to ALPHA.

  DO j = 2, n
    DO i = 1, MIN( j-1, m )
      a( i, j ) = alpha
    END DO
  END DO

ELSE IF( lsame( uplo, 'L' ) ) THEN

!        Set the strictly lower triangular or trapezoidal part of the
!        array to ALPHA.

  DO j = 1, MIN( m, n )
    a( j+1:m, j ) = alpha
  END DO

ELSE

!        Set the leading m-by-n submatrix to ALPHA.

  a( 1:m, 1:n ) = alpha
END IF

!     Set the first min(M,N) diagonal elements to BETA.

DO i = 1, MIN( m, n )
  a( i, i ) = beta
END DO

RETURN

!     End of DLASET

END SUBROUTINE dlaset


SUBROUTINE dlassq( n, x, incx, scale, sumsq )


!  -- LAPACK auxiliary routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997

!     .. Scalar Arguments ..
INTEGER, INTENT(IN)       :: incx, n
REAL (wp), INTENT(IN OUT) :: scale, sumsq
!     ..
!     .. Array Arguments ..
REAL (wp), INTENT(IN)     ::   x(:)
!     ..

!  Purpose
!  =======

!  DLASSQ  returns the values  scl  and  smsq  such that

!     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,

!  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
!  assumed to be non-negative and  scl  returns the value

!     scl = max( scale, abs( x( i ) ) ).

!  scale and sumsq must be supplied in SCALE and SUMSQ and
!  scl and smsq are overwritten on SCALE and SUMSQ respectively.

!  The routine makes only one pass through the vector x.

!  Arguments
!  =========

!  N       (input) INTEGER
!          The number of elements to be used from the vector X.

!  X       (input) DOUBLE PRECISION
!          The vector for which a scaled sum of squares is computed.
!             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.

!  INCX    (input) INTEGER
!          The increment between successive values of the vector X.
!          INCX > 0.

!  SCALE   (input/output) DOUBLE PRECISION
!          On entry, the value  scale  in the equation above.
!          On exit, SCALE is overwritten with  scl , the scaling factor
!          for the sum of squares.

!  SUMSQ   (input/output) DOUBLE PRECISION
!          On entry, the value  sumsq  in the equation above.
!          On exit, SUMSQ is overwritten with  smsq , the basic sum of
!          squares from which  scl  has been factored out.


!     .. Local Scalars ..
INTEGER   :: ix
REAL (wp) :: absxi
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC          ABS
!     ..
!     .. Executable Statements ..

IF( n > 0 ) THEN
  DO ix = 1, 1 + ( n-1 )*incx, incx
    IF( x( ix ) /= zero ) THEN
      absxi = ABS( x( ix ) )
      IF( scale < absxi ) THEN
        sumsq = 1 + sumsq*( scale / absxi )**2
        scale = absxi
      ELSE
        sumsq = sumsq + ( absxi / scale )**2
      END IF
    END IF
  END DO
END IF
RETURN

!     End of DLASSQ

END SUBROUTINE dlassq



FUNCTION dlapy2( x, y ) RESULT(fn_val)

!  -- LAPACK auxiliary routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997

!     .. Scalar Arguments ..
REAL (wp)             :: fn_val
REAL (wp), INTENT(IN) :: x, y
!     ..

!  Purpose
!  =======

!  DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
!  overflow.

!  Arguments
!  =========

!  X       (input) DOUBLE PRECISION
!  Y       (input) DOUBLE PRECISION
!          X and Y specify the values x and y.

!     .. Local Scalars ..
REAL (wp) :: w, xabs, yabs, z
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC          ABS, MAX, MIN, SQRT
!     ..
!     .. Executable Statements ..

xabs = ABS( x )
yabs = ABS( y )
w = MAX( xabs, yabs )
z = MIN( xabs, yabs )
IF( z == zero ) THEN
  fn_val = w
ELSE
  fn_val = w*SQRT( one+( z / w )**2 )
END IF
RETURN

!     End of DLAPY2

END FUNCTION dlapy2


FUNCTION dlanst( norm, n, d, e ) RESULT(fn_val)

!  -- LAPACK auxiliary routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997

!     .. Scalar Arguments ..
CHARACTER (LEN=1), INTENT(IN) :: norm
REAL (wp)                     :: fn_val
INTEGER, INTENT(IN)           :: n
!     ..
!     .. Array Arguments ..
REAL (wp), INTENT(IN)         :: d(:), e(:)
!     ..

!  Purpose
!  =======

!  DLANST  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  real symmetric tridiagonal matrix A.

!  Description
!  ===========

!  DLANST returns the value

!     DLANST = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'

!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.

!  Arguments
!  =========

!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in DLANST as described
!          above.

!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.  When N = 0, DLANST is
!          set to zero.

!  D       (input) DOUBLE PRECISION array, dimension (N)
!          The diagonal elements of A.

!  E       (input) DOUBLE PRECISION array, dimension (N-1)
!          The (n-1) sub-diagonal or super-diagonal elements of A.

!  =====================================================================

!     .. Local Scalars ..
INTEGER   :: i
REAL (wp) :: anorm, scale, sum
!     ..
!     .. External Functions ..
! LOGICAL ::            lsame
! EXTERNAL           lsame
!     ..
!     .. External Subroutines ..
! EXTERNAL           dlassq
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..

IF( n <= 0 ) THEN
  anorm = zero
ELSE IF( lsame( norm, 'M' ) ) THEN

!        Find max(abs(A(i,j))).

  anorm = ABS( d( n ) )
  DO i = 1, n - 1
    anorm = MAX( anorm, ABS( d( i ) ) )
    anorm = MAX( anorm, ABS( e( i ) ) )
  END DO
ELSE IF( lsame( norm, 'O' ) .OR. norm == '1' .OR.  &
  lsame( norm, 'I' ) ) THEN

!        Find norm1(A).

  IF( n == 1 ) THEN
    anorm = ABS( d( 1 ) )
  ELSE
    anorm = MAX( ABS( d( 1 ) )+ABS( e( 1 ) ),ABS( e( n-1 ) )+ABS( d( n ) ) )
    DO i = 2, n - 1
      anorm = MAX( anorm, ABS( d( i ) )+ABS( e( i ) )+ABS( e( i-1 ) ) )
    END DO
  END IF
ELSE IF( ( lsame( norm, 'F' ) ) .OR. ( lsame( norm, 'E' ) ) ) THEN

!        Find normF(A).

  scale = zero
  sum = one
  IF( n > 1 ) THEN
    CALL dlassq( n-1, e, 1, scale, sum )
    sum = 2*sum
  END IF
  CALL dlassq( n, d, 1, scale, sum )
  anorm = scale*SQRT( sum )
END IF

fn_val = anorm
RETURN

!     End of DLANST

END FUNCTION dlanst


FUNCTION dlansy( norm, uplo, n, a, lda ) RESULT(value)

!  -- LAPACK auxiliary routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997
! N.B. The last argument, WORK, has been removed.

!     .. Scalar Arguments ..
CHARACTER (LEN=1), INTENT(IN) :: norm, uplo
REAL (wp)                     :: value
INTEGER, INTENT(IN)           :: lda, n
!     ..
!     .. Array Arguments ..
REAL (wp), INTENT(IN)         :: a(:,:)
!     ..

!  Purpose
!  =======

!  DLANSY  returns the value of the one norm,  or the Frobenius norm, or
!  the  infinity norm,  or the  element of  largest absolute value  of a
!  real symmetric matrix A.

!  Description
!  ===========

!  DLANSY returns the value

!     DLANSY = ( max(abs(A(i,j))), NORM = 'M' or 'm'
!              (
!              ( norm1(A),         NORM = '1', 'O' or 'o'
!              (
!              ( normI(A),         NORM = 'I' or 'i'
!              (
!              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'

!  where  norm1  denotes the  one norm of a matrix (maximum column sum),
!  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
!  normF  denotes the  Frobenius norm of a matrix (square root of sum of
!  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.

!  Arguments
!  =========

!  NORM    (input) CHARACTER*1
!          Specifies the value to be returned in DLANSY as described
!          above.

!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is to be referenced.
!          = 'U':  Upper triangular part of A is referenced
!          = 'L':  Lower triangular part of A is referenced

!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.  When N = 0, DLANSY is
!          set to zero.

!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The symmetric matrix A.  If UPLO = 'U', the leading n by n
!          upper triangular part of A contains the upper triangular part
!          of the matrix A, and the strictly lower triangular part of A
!          is not referenced.  If UPLO = 'L', the leading n by n lower
!          triangular part of A contains the lower triangular part of
!          the matrix A, and the strictly upper triangular part of A is
!          not referenced.

!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(N,1).

!     .. Local Scalars ..
INTEGER   :: i, j
REAL (wp) :: absa, scale, sum
REAL (wp) :: work(n)
!     ..
!     .. External Subroutines ..
! EXTERNAL           dlassq
!     ..
!     .. External Functions ..
! LOGICAL ::            lsame
! EXTERNAL           lsame
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC          ABS, MAX, SQRT
!     ..
!     .. Executable Statements ..

IF( n == 0 ) THEN
  value = zero
ELSE IF( lsame( norm, 'M' ) ) THEN

!        Find max(abs(A(i,j))).

  value = zero
  IF( lsame( uplo, 'U' ) ) THEN
    DO j = 1, n
      DO i = 1, j
        value = MAX( value, ABS( a( i, j ) ) )
      END DO
    END DO
  ELSE
    DO j = 1, n
      DO i = j, n
        value = MAX( value, ABS( a( i, j ) ) )
      END DO
    END DO
  END IF
ELSE IF( ( lsame( norm, 'I' ) ) .OR. ( lsame( norm, 'O' ) ) .OR.  &
  ( norm == '1' ) ) THEN

!        Find normI(A) ( = norm1(A), since A is symmetric).

  value = zero
  IF( lsame( uplo, 'U' ) ) THEN
    DO j = 1, n
      sum = zero
      DO i = 1, j - 1
        absa = ABS( a( i, j ) )
        sum = sum + absa
        work( i ) = work( i ) + absa
      END DO
      work( j ) = sum + ABS( a( j, j ) )
    END DO
    value = MAXVAL( work(1:n) )
  ELSE
    work( 1:n ) = zero
    DO j = 1, n
      sum = work( j ) + ABS( a( j, j ) )
      DO i = j + 1, n
        absa = ABS( a( i, j ) )
        sum = sum + absa
        work( i ) = work( i ) + absa
      END DO
      value = MAX( value, sum )
    END DO
  END IF
ELSE IF( ( lsame( norm, 'F' ) ) .OR. ( lsame( norm, 'E' ) ) ) THEN

!        Find normF(A).

  scale = zero
  sum = one
  IF( lsame( uplo, 'U' ) ) THEN
    DO j = 2, n
      CALL dlassq( j-1, a( :, j ), 1, scale, sum )
    END DO
  ELSE
    DO j = 1, n - 1
      CALL dlassq( n-j, a( j+1:, j ), 1, scale, sum )
    END DO
  END IF
  sum = 2*sum
  CALL dlassq( n, a(:,1), lda+1, scale, sum )
  value = scale*SQRT( sum )
END IF

RETURN

!     End of DLANSY

END FUNCTION dlansy



END MODULE dla

MODULE dor
! This module contains translations to ELF90 compatability of the
! LAPACK Orthogonal Real routines

! Translated by Alan Miller
! alan @ vic.cmis.csiro.au    URL: www.ozemail.com.au/~milleraj

! Latest revision - 30 October 1997

USE la_precision, ONLY: wp => dp
USE la_auxmod
USE dblas
USE dla, ONLY: dlarf, dlarft, dlarfb
IMPLICIT NONE

CONTAINS



SUBROUTINE dorg2l( m, n, k, a, lda, tau, info )

!  -- LAPACK routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997
! N.B. The workspace has been removed from the list of arguments.

!     .. Scalar Arguments ..
INTEGER, INTENT(IN)       :: k, lda, m, n
INTEGER, INTENT(OUT)      :: info
!     ..
!     .. Array Arguments ..
REAL (wp), INTENT(IN OUT) :: a(:,:)
REAL (wp), INTENT(IN)     :: tau(:)
!     ..

!  Purpose
!  =======

!  DORG2L generates an m by n real matrix Q with orthonormal columns,
!  which is defined as the last n columns of a product of k elementary
!  reflectors of order m

!        Q  =  H(k) . . . H(2) H(1)

!  as returned by DGEQLF.

!  Arguments
!  =========

!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.

!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.

!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.

!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the (n-k+i)-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by DGEQLF in the last k columns of its array
!          argument A.
!          On exit, the m by n matrix Q.

!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).

!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQLF.

!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument has an illegal value

!  =====================================================================

!     .. Local Scalars ..
INTEGER   :: i, ii, j, l
!     ..
!     .. External Subroutines ..
! EXTERNAL           dlarf, dscal, xerbla
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC          MAX
!     ..
!     .. Executable Statements ..

!     Test the input arguments

info = 0
IF( m < 0 ) THEN
  info = -1
ELSE IF( n < 0 .OR. n > m ) THEN
  info = -2
ELSE IF( k < 0 .OR. k > n ) THEN
  info = -3
ELSE IF( lda < MAX( 1, m ) ) THEN
  info = -5
END IF
IF( info /= 0 ) THEN
  CALL erinfo(-info, 'DORG2L')
  RETURN
END IF

!     Quick return if possible

IF( n <= 0 ) RETURN

!     Initialise columns 1:n-k to columns of the unit matrix

DO j = 1, n - k
  DO l = 1, m
    a( l, j ) = zero
  END DO
  a( m-n+j, j ) = one
END DO

DO i = 1, k
  ii = n - k + i

!        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left

  a( m-n+ii, ii ) = one
  CALL dlarf( 'Left', m-n+ii, ii-1, a( :, ii ), 1, tau( i ), a, lda )
  CALL dscal( m-n+ii-1, -tau( i ), a( :, ii ), 1 )
  a( m-n+ii, ii ) = one - tau( i )

!        Set A(m-k+i+1:m,n-k+i) to zero

  DO l = m - n + ii + 1, m
    a( l, ii ) = zero
  END DO
END DO
RETURN

!     End of DORG2L

END SUBROUTINE dorg2l


SUBROUTINE dorg2r( m, n, k, a, lda, tau, info )

!  -- LAPACK routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997
! N.B. The workspace has been removed from the list of arguments.

!     .. Scalar Arguments ..
INTEGER, INTENT(IN)       :: k, lda, m, n
INTEGER, INTENT(OUT)      :: info
!     ..
!     .. Array Arguments ..
REAL (wp), INTENT(IN)     :: tau(:)
REAL (wp), INTENT(IN OUT) :: a(:,:)
!     ..

!  Purpose
!  =======

!  DORG2R generates an m by n real matrix Q with orthonormal columns,
!  which is defined as the first n columns of a product of k elementary
!  reflectors of order m

!        Q  =  H(1) H(2) . . . H(k)

!  as returned by DGEQRF.

!  Arguments
!  =========

!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.

!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.

!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.

!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the i-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by DGEQRF in the first k columns of its array
!          argument A.
!          On exit, the m-by-n matrix Q.

!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).

!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.

!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument has an illegal value

!  =====================================================================

!     .. Local Scalars ..
INTEGER   :: i, j, l
!     ..
!     .. External Subroutines ..
! EXTERNAL           dlarf, dscal, xerbla
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC          MAX
!     ..
!     .. Executable Statements ..

!     Test the input arguments

info = 0
IF( m < 0 ) THEN
  info = -1
ELSE IF( n < 0 .OR. n > m ) THEN
  info = -2
ELSE IF( k < 0 .OR. k > n ) THEN
  info = -3
ELSE IF( lda < MAX( 1, m ) ) THEN
  info = -5
END IF
IF( info /= 0 ) THEN
  CALL erinfo(-info, 'DORG2R')
  RETURN
END IF

!     Quick return if possible

IF( n <= 0 ) RETURN

!     Initialise columns k+1:n to columns of the unit matrix

DO j = k + 1, n
  DO l = 1, m
    a( l, j ) = zero
  END DO
  a( j, j ) = one
END DO

DO i = k, 1, -1

!        Apply H(i) to A(i:m,i:n) from the left

  IF( i < n ) THEN
    a( i, i ) = one
    CALL dlarf( 'Left', m-i+1, n-i, a( i:, i ), 1, tau( i ),  &
                a( i:, i+1: ), lda )
  END IF
  IF( i < m ) CALL dscal( m-i, -tau( i ), a( i+1:, i ), 1 )
  a( i, i ) = one - tau( i )

!        Set A(1:i-1,i) to zero

  DO l = 1, i - 1
    a( l, i ) = zero
  END DO
END DO
RETURN

!     End of DORG2R

END SUBROUTINE dorg2r



SUBROUTINE dorgql( m, n, k, a, lda, tau, info )

!  -- LAPACK routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997
! N.B. Arguments WORK and LWORK have been removed.

!     .. Scalar Arguments ..
INTEGER, INTENT(IN)       :: k, lda, m, n
INTEGER, INTENT(OUT)      :: info
!     ..
!     .. Array Arguments ..
REAL (wp), INTENT(IN)     :: tau(:)
REAL (wp), INTENT(IN OUT) :: a(:,:)
!     ..

!  Purpose
!  =======

!  DORGQL generates an m by n real matrix Q with orthonormal columns,
!  which is defined as the last n columns of a product of k elementary
!  reflectors of order m

!        Q  =  H(k) . . . H(2) H(1)

!  as returned by DGEQLF.

!  Arguments
!  =========

!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.

!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.

!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.

!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the (n-k+i)-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by DGEQLF in the last k columns of its array
!          argument A.
!          On exit, the m by n matrix Q.

!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).

!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQLF.

!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument has an illegal value

!  =====================================================================

!     .. Local Scalars ..
INTEGER                :: i, ib, j, kk, l, nb, nbmin, nx
REAL (wp)              :: t(k,k)
!     ..
!     .. External Subroutines ..
! EXTERNAL           dlarfb, dlarft, dorg2l, xerbla
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
! INTEGER ::            ilaenv
! EXTERNAL           ilaenv
!     ..
!     .. Executable Statements ..

!     Test the input arguments

info = 0
IF( m < 0 ) THEN
  info = -1
ELSE IF( n < 0 .OR. n > m ) THEN
  info = -2
ELSE IF( k < 0 .OR. k > n ) THEN
  info = -3
ELSE IF( lda < MAX( 1, m ) ) THEN
  info = -5
END IF
IF( info /= 0 ) THEN
  CALL erinfo(-info, 'DORGQL')
  RETURN
END IF

!     Quick return if possible

IF( n <= 0 ) THEN
  RETURN
END IF

!     Determine the block size.

nb = ilaenv( 1, 'DORGQL', ' ', m, n, k, -1 )
nbmin = 2
nx = 0
IF( nb > 1 .AND. nb < k ) THEN

!        Determine when to cross over from blocked to unblocked code.

  nx = MAX( 0, ilaenv( 3, 'DORGQL', ' ', m, n, k, -1 ) )
END IF

IF( nb >= nbmin .AND. nb < k .AND. nx < k ) THEN

!        Use blocked code after the first block.
!        The last kk columns are handled by the block method.

  kk = MIN( k, ( ( k-nx+nb-1 ) / nb )*nb )

!        Set A(m-kk+1:m,1:n-kk) to zero.

  DO j = 1, n - kk
    DO i = m - kk + 1, m
      a( i, j ) = zero
    END DO
  END DO
ELSE
  kk = 0
END IF

!     Use unblocked code for the first or only block.

CALL dorg2l( m-kk, n-kk, k-kk, a, lda, tau, info )

IF( kk > 0 ) THEN

!        Use blocked code

  DO i = k - kk + 1, k, nb
    ib = MIN( nb, k-i+1 )
    IF( n-k+i > 1 ) THEN

!              Form the triangular factor of the block reflector
!              H = H(i+ib-1) . . . H(i+1) H(i)

      CALL dlarft( 'Backward', 'Columnwise', m-k+i+ib-1, ib,  &
                   a( :, n-k+i: ), lda, tau( i: ), t, k )

!              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left

      CALL dlarfb( 'Left', 'No transpose', 'Backward', 'Columnwise',  &
                  m-k+i+ib-1, n-k+i-1, ib, a( :, n-k+i: ), lda, t, k, a, lda )
    END IF

!           Apply H to rows 1:m-k+i+ib-1 of current block

    CALL dorg2l( m-k+i+ib-1, ib, ib, a( :, n-k+i: ), lda, tau( i: ), info )

!           Set rows m-k+i+ib:m of current block to zero

    DO j = n - k + i, n - k + i + ib - 1
      DO l = m - k + i + ib, m
        a( l, j ) = zero
      END DO
    END DO
  END DO
END IF

RETURN

!     End of DORGQL

END SUBROUTINE dorgql


SUBROUTINE dorgqr( m, n, k, a, lda, tau, info )

!  -- LAPACK routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997
! N.B. Arguments WORK and LWORK have been removed.

!     .. Scalar Arguments ..
INTEGER, INTENT(IN)       :: k, lda, m, n
INTEGER, INTENT(OUT)      :: info
!     ..
!     .. Array Arguments ..
REAL (wp), INTENT(IN)     :: tau(:)
REAL (wp), INTENT(IN OUT) :: a(:,:)
!     ..

!  Purpose
!  =======

!  DORGQR generates an m by n real matrix Q with orthonormal columns,
!  which is defined as the first n columns of a product of k elementary
!  reflectors of order m

!        Q  =  H(1) H(2) . . . H(k)

!  as returned by DGEQRF.

!  Arguments
!  =========

!  M       (input) INTEGER
!          The number of rows of the matrix Q. M >= 0.

!  N       (input) INTEGER
!          The number of columns of the matrix Q. M >= N >= 0.

!  K       (input) INTEGER
!          The number of elementary reflectors whose product defines the
!          matrix Q. N >= K >= 0.

!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the i-th column must contain the vector which
!          defines the elementary reflector H(i), for i = 1,2,...,k, as
!          returned by DGEQRF in the first k columns of its array
!          argument A.
!          On exit, the m by n matrix Q.

!  LDA     (input) INTEGER
!          The first dimension of the array A. LDA >= max(1,M).

!  TAU     (input) DOUBLE PRECISION array, dimension (K)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DGEQRF.

!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument has an illegal value

!  =====================================================================

!     .. Local Scalars ..
INTEGER                :: i, ib, j, ki, kk, l, nb, nbmin, nx
REAL (wp)              :: t(k,k)
!     ..
!     .. External Subroutines ..
! EXTERNAL           dlarfb, dlarft, dorg2r, xerbla
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC          MAX, MIN
!     ..
!     .. External Functions ..
! INTEGER ::            ilaenv
! EXTERNAL           ilaenv
!     ..
!     .. Executable Statements ..

!     Test the input arguments

info = 0
IF( m < 0 ) THEN
  info = -1
ELSE IF( n < 0 .OR. n > m ) THEN
  info = -2
ELSE IF( k < 0 .OR. k > n ) THEN
  info = -3
ELSE IF( lda < MAX( 1, m ) ) THEN
  info = -5
END IF
IF( info /= 0 ) THEN
  CALL erinfo(-info, 'DORGQR')
  RETURN
END IF

!     Quick return if possible

IF( n <= 0 ) THEN
  RETURN
END IF

!     Determine the block size.

nb = ilaenv( 1, 'DORGQR', ' ', m, n, k, -1 )
nbmin = 2
nx = 0
IF( nb > 1 .AND. nb < k ) THEN

!        Determine when to cross over from blocked to unblocked code.

  nx = MAX( 0, ilaenv( 3, 'DORGQR', ' ', m, n, k, -1 ) )
END IF

IF( nb >= nbmin .AND. nb < k .AND. nx < k ) THEN

!        Use blocked code after the last block.
!        The first kk columns are handled by the block method.

  ki = ( ( k-nx-1 ) / nb )*nb
  kk = MIN( k, ki+nb )

!        Set A(1:kk,kk+1:n) to zero.

  DO j = kk + 1, n
    DO i = 1, kk
      a( i, j ) = zero
    END DO
  END DO
ELSE
  kk = 0
END IF

!     Use unblocked code for the last or only block.

IF( kk < n ) CALL dorg2r( m-kk, n-kk, k-kk, a( kk+1:, kk+1: ), lda,  &
                          tau( kk+1: ), info )

IF( kk > 0 ) THEN

!        Use blocked code

  DO i = ki + 1, 1, -nb
    ib = MIN( nb, k-i+1 )
    IF( i+ib <= n ) THEN

!              Form the triangular factor of the block reflector
!              H = H(i) H(i+1) . . . H(i+ib-1)

      CALL dlarft( 'Forward', 'Columnwise', m-i+1, ib,  &
                   a( i:, i: ), lda, tau( i: ), t, k )

!              Apply H to A(i:m,i+ib:n) from the left

      CALL dlarfb( 'Left', 'No transpose', 'Forward', 'Columnwise', m-i+1,  &
                   n-i-ib+1, ib, a( i:, i: ), lda, t, k, a( i:, i+ib: ), lda )
    END IF

!           Apply H to rows i:m of current block

    CALL dorg2r( m-i+1, ib, ib, a( i:, i: ), lda, tau( i: ), info )

!           Set rows 1:i-1 of current block to zero

    DO j = i, i + ib - 1
      DO l = 1, i - 1
        a( l, j ) = zero
      END DO
    END DO
  END DO
END IF

RETURN

!     End of DORGQR

END SUBROUTINE dorgqr



SUBROUTINE dorgtr( uplo, n, a, lda, tau, info )

!  -- LAPACK routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997
! N.B. Arguments WORK and LWORK have been removed.

!     .. Scalar Arguments ..
CHARACTER (LEN=1), INTENT(IN) :: uplo
INTEGER, INTENT(IN)           :: lda, n
INTEGER, INTENT(OUT)          :: info
!     ..
!     .. Array Arguments ..
REAL (wp), INTENT(IN)         :: tau(:)
REAL (wp), INTENT(IN OUT)     :: a(:,:)
!     ..

!  Purpose
!  =======

!  DORGTR generates a real orthogonal matrix Q which is defined as the
!  product of n-1 elementary reflectors of order n, as returned by
!  DSYTRD:

!  if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),

!  if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).

!  Arguments
!  =========

!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangle of the array A
!          holds details of the elementary reflectors, as returned by
!          DSYTRD:
!          = 'U': Upper triangle;
!          = 'L': Lower triangle.

!  N       (input) INTEGER
!          The order of the matrix Q. N >= 0.

!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the vectors which define the elementary reflectors,
!          as returned by DSYTRD.
!          On exit, the n by n orthogonal matrix Q.

!  LDA     (input) INTEGER
!          The leading dimension of the array A. LDA >= max(1,N).

!  TAU     (input) DOUBLE PRECISION array, dimension (N-1)
!          TAU(i) must contain the scalar factor of the elementary
!          reflector H(i), as returned by DSYTRD.

!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -i, the i-th argument had an illegal value

!  =====================================================================

!     .. Local Scalars ..
LOGICAL :: upper
INTEGER :: i, j
!     ..
!     .. External Functions ..
! LOGICAL ::            lsame
! EXTERNAL           lsame
!     ..
!     .. External Subroutines ..
! EXTERNAL           dorgql, dorgqr, xerbla
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC          MAX
!     ..
!     .. Executable Statements ..

!     Test the input arguments

info = 0
upper = lsame( uplo, 'U' )
IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
  info = -1
ELSE IF( n < 0 ) THEN
  info = -2
ELSE IF( lda < MAX( 1, n ) ) THEN
  info = -4
END IF
IF( info /= 0 ) THEN
  CALL erinfo(-info, 'DORGTR')
  RETURN
END IF

!     Quick return if possible

IF( n == 0 ) THEN
  RETURN
END IF

IF( upper ) THEN

!        Q was determined by a call to DSYTRD with UPLO = 'U'

!        Shift the vectors which define the elementary reflectors one
!        column to the left, and set the last row and column of Q to
!        those of the unit matrix

  DO j = 1, n - 1
    DO i = 1, j - 1
      a( i, j ) = a( i, j+1 )
    END DO
    a( n, j ) = zero
  END DO
  DO i = 1, n - 1
    a( i, n ) = zero
  END DO
  a( n, n ) = one

!        Generate Q(1:n-1,1:n-1)

  CALL dorgql( n-1, n-1, n-1, a, lda, tau, info )

ELSE

!        Q was determined by a call to DSYTRD with UPLO = 'L'.

!        Shift the vectors which define the elementary reflectors one
!        column to the right, and set the first row and column of Q to
!        those of the unit matrix

  DO j = n, 2, -1
    a( 1, j ) = zero
    DO i = j + 1, n
      a( i, j ) = a( i, j-1 )
    END DO
  END DO
  a( 1, 1 ) = one
  DO i = 2, n
    a( i, 1 ) = zero
  END DO
  IF( n > 1 ) THEN

!           Generate Q(2:n,2:n)

    CALL dorgqr( n-1, n-1, n-1, a( 2:, 2: ), lda, tau, info )
  END IF
END IF

RETURN

!     End of DORGTR

END SUBROUTINE dorgtr


END MODULE dor

MODULE dge
! This module contains translations to ELF90 compatability of the
! LAPACK GEneral routines + DBDSQR

! Translated by Alan Miller
! alan @ vic.cmis.csiro.au    URL: www.ozemail.com.au/~milleraj

! Latest revision - 10 November 1997

USE la_precision, ONLY: wp => dp
USE la_auxmod
USE dblas

CONTAINS



SUBROUTINE dgetrf( m, n, a, lda, ipiv, info )

!  -- LAPACK routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997

USE dblas3, ONLY: dgemm, dtrsm
USE dla,    ONLY: dlaswp
IMPLICIT NONE

!     .. Scalar Arguments ..
INTEGER, INTENT(IN)       :: lda, m, n
INTEGER, INTENT(OUT)      :: info
!     ..
!     .. Array Arguments ..
INTEGER, INTENT(OUT)      :: ipiv(:)
REAL (wp), INTENT(IN OUT) :: a(:,:)
!     ..

!  Purpose
!  =======

!  DGETRF computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.

!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).

!  This is the right-looking Level 3 BLAS version of the algorithm.

!  Arguments
!  =========

!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.

!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.

!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.

!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).

!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).

!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.

!  =====================================================================

!     .. Local Scalars ..
INTEGER :: i, iinfo, j, jb, nb
!     ..
!     .. External Subroutines ..
! EXTERNAL           dgemm, dgetf2, dlaswp, dtrsm, xerbla
!     ..
!     .. External Functions ..
! INTEGER ::            ilaenv
! EXTERNAL           ilaenv
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..

!     Test the input parameters.

info = 0
IF( m < 0 ) THEN
  info = -1
ELSE IF( n < 0 ) THEN
  info = -2
ELSE IF( lda < MAX( 1, m ) ) THEN
  info = -4
END IF
IF( info /= 0 ) THEN
  CALL erinfo(-info, 'DGETRF')
  RETURN
END IF

!     Quick return if possible

IF( m == 0 .OR. n == 0 ) RETURN

!     Determine the block size for this environment.

nb = ilaenv( 1, 'DGETRF', ' ', m, n, -1, -1 )
IF( nb <= 1 .OR. nb >= MIN( m, n ) ) THEN

!        Use unblocked code.

  CALL dgetf2( m, n, a, lda, ipiv, info )
ELSE

!        Use blocked code.

  DO j = 1, MIN( m, n ), nb
    jb = MIN( MIN( m, n )-j+1, nb )

!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.

    CALL dgetf2( m-j+1, jb, a( j:, j: ), lda, ipiv( j: ), iinfo )

!           Adjust INFO and the pivot indices.

    IF( info == 0 .AND. iinfo > 0 ) info = iinfo + j - 1
    DO i = j, MIN( m, j+jb-1 )
      ipiv( i ) = j - 1 + ipiv( i )
    END DO

!           Apply interchanges to columns 1:J-1.

    CALL dlaswp( j-1, a, j, j+jb-1, ipiv, 1 )

    IF( j+jb <= n ) THEN

!              Apply interchanges to columns J+JB:N.

      CALL dlaswp( n-j-jb+1, a( :, j+jb: ), j, j+jb-1, ipiv, 1 )

!              Compute block row of U.

      CALL dtrsm( 'Left', 'Lower', 'No transpose', 'Unit', jb,  &
                  n-j-jb+1, one, a( j:, j: ), lda, a( j:, j+jb: ), lda )
      IF( j+jb <= m ) THEN

!                 Update trailing submatrix.

        CALL dgemm( 'No transpose', 'No transpose', m-j-jb+1,  &
                    n-j-jb+1, jb, -one, a( j+jb:, j: ), lda,  &
                    a( j:, j+jb: ), lda, one, a( j+jb:, j+jb: ), lda )
      END IF
    END IF
  END DO
END IF
RETURN

!     End of DGETRF

END SUBROUTINE dgetrf



SUBROUTINE dgetrs( trans, n, nrhs, a, lda, ipiv, b, ldb, info )

!  -- LAPACK routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997

USE dblas3, ONLY: dtrsm
USE dla,    ONLY: dlaswp
IMPLICIT NONE

!     .. Scalar Arguments ..
CHARACTER (LEN=1), INTENT(IN) :: trans
INTEGER, INTENT(IN)           :: lda, ldb, n, nrhs
INTEGER, INTENT(OUT)          :: info
!     ..
!     .. Array Arguments ..
INTEGER, INTENT(IN)           :: ipiv(:)
REAL (wp), INTENT(IN)         :: a(:,:)
REAL (wp), INTENT(IN OUT)     :: b(:,:)
!     ..

!  Purpose
!  =======

!  DGETRS solves a system of linear equations
!     A * X = B  or  A' * X = B
!  with a general n by n matrix A using the LU factorization computed
!  by DGETRF.

!  Arguments
!  =========

!  TRANS   (input) CHARACTER*1
!          Specifies the form of the system of equations.
!          = 'N':  A * X = B  (No transpose)
!          = 'T':  A'* X = B  (Transpose)
!          = 'C':  A'* X = B  (Conjugate transpose = Transpose)

!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.

!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.

!  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
!          The factors L and U from the factorization A = P*L*U
!          as computed by DGETRF.

!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).

!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices from DGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).

!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the right hand side vectors B for the system of
!          linear equations.
!          On exit, the solution vectors, X.

!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).

!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value

!  =====================================================================

!     ..
!     .. Local Scalars ..
LOGICAL ::            notran
!     ..
!     .. External Functions ..
! LOGICAL ::            lsame
! EXTERNAL           lsame
!     ..
!     .. External Subroutines ..
! EXTERNAL           dlaswp, dtrsm, xerbla
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC          MAX
!     ..
!     .. Executable Statements ..

!     Test the input parameters.

info = 0
notran = lsame( trans, 'N' )
IF( .NOT.notran .AND. .NOT.lsame( trans, 'T' ) .AND. .NOT.  &
  lsame( trans, 'C' ) ) THEN
  info = -1
ELSE IF( n < 0 ) THEN
  info = -2
ELSE IF( nrhs < 0 ) THEN
  info = -3
ELSE IF( lda < MAX( 1, n ) ) THEN
  info = -5
ELSE IF( ldb < MAX( 1, n ) ) THEN
  info = -8
END IF
IF( info /= 0 ) THEN
  CALL erinfo(-info, 'DGETRS')
  RETURN
END IF

!     Quick return if possible

IF( n == 0 .OR. nrhs == 0 ) RETURN

IF( notran ) THEN

!        Solve A * X = B.

!        Apply row interchanges to the right hand sides.

  CALL dlaswp( nrhs, b, 1, n, ipiv, 1 )

!        Solve L*X = B, overwriting B with X.

  CALL dtrsm( 'Left', 'Lower', 'No transpose', 'Unit', n, nrhs,  &
              one, a, lda, b, ldb )

!        Solve U*X = B, overwriting B with X.

  CALL dtrsm( 'Left', 'Upper', 'No transpose', 'Non-unit', n,  &
              nrhs, one, a, lda, b, ldb )
ELSE

!        Solve A' * X = B.

!        Solve U'*X = B, overwriting B with X.

  CALL dtrsm( 'Left', 'Upper', 'Transpose', 'Non-unit', n, nrhs,  &
              one, a, lda, b, ldb )

!        Solve L'*X = B, overwriting B with X.

  CALL dtrsm( 'Left', 'Lower', 'Transpose', 'Unit', n, nrhs, one,  &
              a, lda, b, ldb )

!        Apply row interchanges to the solution vectors.

  CALL dlaswp( nrhs, b, 1, n, ipiv, -1 )
END IF

RETURN

!     End of DGETRS

END SUBROUTINE dgetrs




SUBROUTINE dgetf2( m, n, a, lda, ipiv, info )

!  -- LAPACK routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997

USE dblas2, ONLY: dger
IMPLICIT NONE

!     .. Scalar Arguments ..
INTEGER, INTENT(IN)       :: lda, m, n
INTEGER, INTENT(OUT)      :: info
!     ..
!     .. Array Arguments ..
INTEGER, INTENT(OUT)      :: ipiv(:)
REAL (wp), INTENT(IN OUT) :: a(:,:)
!     ..

!  Purpose
!  =======

!  DGETF2 computes an LU factorization of a general m-by-n matrix A
!  using partial pivoting with row interchanges.

!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).

!  This is the right-looking Level 2 BLAS version of the algorithm.

!  Arguments
!  =========

!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.

!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.

!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the m by n matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.

!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).

!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).

!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!               has been completed, but the factor U is exactly
!               singular, and division by zero will occur if it is used
!               to solve a system of equations.

!  =====================================================================

!     .. Local Scalars ..
INTEGER :: j, jp
!     ..
!     .. External Functions ..
! INTEGER ::            idamax
! EXTERNAL           idamax
!     ..
!     .. External Subroutines ..
! EXTERNAL           dger, dscal, dswap, xerbla
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..

!     Test the input parameters.

info = 0
IF( m < 0 ) THEN
  info = -1
ELSE IF( n < 0 ) THEN
  info = -2
ELSE IF( lda < MAX( 1, m ) ) THEN
  info = -4
END IF
IF( info /= 0 ) THEN
  CALL erinfo(-info, 'DGETF2')
  RETURN
END IF

!     Quick return if possible

IF( m == 0 .OR. n == 0 ) RETURN

DO j = 1, MIN( m, n )

!        Find pivot and test for singularity.

  jp = j - 1 + idamax( m-j+1, a( j:, j ), 1 )
  ipiv( j ) = jp
  
  
  IF( a( jp, j ) /= zero ) THEN

!           Apply the interchange to columns 1:N.

    IF( jp /= j ) CALL dswap( n, a( j, : ), 1, a( jp, : ), 1 )

!           Compute elements J+1:M of J-th column.

    IF( j < m ) CALL dscal( m-j, one / a( j, j ), a( j+1:, j ), 1 )

  ELSE IF( info == 0 ) THEN

    info = j
  END IF

  IF( j+1 <= n ) THEN

!           Update trailing submatrix.

    CALL dger( m-j, n-j, -one, a( j+1:, j ), 1, a( j, j+1: ), 1,  &
               a( j+1:, j+1: ), lda )
  END IF
END DO

RETURN

!     End of DGETF2

END SUBROUTINE dgetf2



END MODULE dge

! Module contains interface to all necessary subroutines from LAPACK

module lapack90


INTERFACE
  SUBROUTINE DSYEV_F90( A, W, JOBZ, UPLO, INFO )
    USE LA_PRECISION, ONLY: WP => DP
    USE LA_AUXMOD
!    USE dsy, ONLY: dsyev
    IMPLICIT NONE
    CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: JOBZ, UPLO
    INTEGER, INTENT(OUT), OPTIONAL         :: INFO
    REAL(WP), INTENT(IN OUT)               :: A(:,:)
    REAL(WP), INTENT(OUT)                  :: W(:)
  END SUBROUTINE DSYEV_F90
END INTERFACE

INTERFACE
  SUBROUTINE DGESV_F90( A, B, IPIV, INFO )
    USE la_precision, ONLY: wp => dp
    USE la_auxmod
    USE dblas
    USE dblas2
    USE dblas3
    USE dla
    USE dge, ONLY: dgetrf, dgetrs, dgetf2
    IMPLICIT NONE
    INTEGER, INTENT(OUT), OPTIONAL         :: INFO
    INTEGER, INTENT(OUT), OPTIONAL, TARGET :: IPIV(:)
    REAL(WP), INTENT(IN OUT)               :: A(:,:), B(:,:)
  END SUBROUTINE DGESV_F90
END INTERFACE


end module

SUBROUTINE DGESV_F90( A, B, IPIV, INFO )
!
!  -- LAPACK90 interface driver routine (version 1.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     May 31, 1997
!
! ELF90 version by Alan Miller  23 October 1997
! Latest revision - 14 November 1997
!
!     .. "Use Statements" ..
!      USE LA_PRECISION, ONLY: WP => DP
!      USE LA_AUXMOD, ONLY: ERINFO
!      USE F77_LAPACK, ONLY: GESV_F77 => LA_GESV
USE la_precision, ONLY: wp => dp
USE la_auxmod
!
!     .. "Implicit Statement" ..
      IMPLICIT NONE
!     .. "Scalar Arguments" ..
      INTEGER, INTENT(OUT), OPTIONAL         :: INFO
!     .. "Array Arguments" ..
      INTEGER, INTENT(OUT), OPTIONAL, TARGET :: IPIV(:)
      REAL(WP), INTENT(IN OUT)               :: A(:,:), B(:,:)
!-----------------------------------------------------------------
!
! Purpose
! =======
!
! LA_GESV computes the solution to either a real or complex system of
! linear equations  AX = B,
! where A is a square matrix and X and B are either rectangular
! matrices or vectors.
!
! The LU decomposition with partial pivoting and row interchanges is
! used to factor A as  A = PLU,
! where P is a permutation matrix, L is unit lower triangular, and U is
! upper triangular. The factored form of A is then used to solve the
! system of equations AX = B.
!
! Arguments
! =========
!
! SUBROUTINE LA_GESV ( A, B, IPIV, INFO )
!    <type>(<wp>), INTENT(INOUT) :: A(:,:), <rhs>
!    INTEGER, INTENT(OUT), OPTIONAL :: IPIV(:)
!    INTEGER, INTENT(OUT), OPTIONAL :: INFO
!    where
!    <type> ::= REAL | COMPLEX
!    <wp>   ::= KIND(1.0) | KIND(1.0D0)
!    <rhs>  ::= B(:,:) | B(:)
!
! =====================
!
! A    (input/output) either REAL or COMPLEX square array,
!      shape (:,:), size(A,1) == size(A,2).
!      On entry, the matrix A.
!      On exit, the factors L and U from the factorization A = PLU;
!         the unit diagonal elements of L are not stored.
!
! B    (input/output) either REAL or COMPLEX rectangular array,
!      shape either (:,:) or (:), size(B,1) or size(B) == size(A,1).
!      On entry, the right hand side vector(s) of matrix B for the
!         system of equations AX = B.
!      On exit, if there is no error, the matrix of solution
!         vector(s) X.
!
! IPIV Optional (output) INTEGER array, shape (:),
!      size(IPIV) == size(A,1). If IPIV is present it indice define
!      the permutation matrix P; row i of the matrix was interchanged
!      with row IPIV(i).
!
! INFO Optional (output) INTEGER.
!      If INFO is present
!         = 0: successful exit
!         < 0: if INFO = -k, the k-th argument had an illegal value
!         > 0: if INFO = k, U(k,k) is exactly zero.  The factorization
!              has been completed, but the factor U is exactly
!              singular, so the solution could not be computed.
!      If INFO is not present and an error occurs, then the program is
!         terminated with an error message.
!-------------------------------------
!     .. "Parameters" ..
      CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_GESV'
!     .. LOCAL SCALARS ..
      INTEGER                     :: LINFO, ISTAT, ISTAT1, SIPIV
!     .. "Local Pointers" ..
      INTEGER, POINTER            :: LPIV(:)
!     .. "Intrinsic Functions" ..
!      INTRINSIC SIZE, PRESENT, MAX
!     ..
!     .. "Executable Statements" ..
!     ..
      LINFO = 0
      ISTAT = 0
      IF ( PRESENT(IPIV) )THEN
         SIPIV = SIZE(IPIV)
      ELSE
         SIPIV = SIZE(A,1)
      END IF
!     ..
!     .. "Test the arguments" ..
!     ..
    !  print*,A
    !  print*,size(A,2),size(A,1)
      IF( SIZE( A, 2 ) /= SIZE(A,1) .OR. SIZE(A,1) < 0 ) THEN
         LINFO = -1
      ELSE IF( SIZE( B, 1 ) /= SIZE(A,1) .OR. SIZE(B,2) < 0 ) THEN
         LINFO = -2
      ELSE IF( SIPIV /= SIZE(A,1) ) THEN
            LINFO = -3
      ELSE IF( SIZE(A,1) > 0 ) THEN
         IF( PRESENT(IPIV) ) THEN
            LPIV => IPIV
         ELSE
            ALLOCATE( LPIV(SIZE(A,1)), STAT = ISTAT )
         END IF
         IF( ISTAT == 0 ) THEN
!     ..
!        .. "Call LAPACK90 routine" ..
!     ..
            CALL DGESV( SIZE(A,1), SIZE(B,2), A, MAX(1,SIZE(A,1)), LPIV, &
                           B, MAX(1,SIZE(A,1)), LINFO )
         ELSE
            LINFO = -100
         END IF
         IF( .NOT. PRESENT(IPIV) ) THEN
            DEALLOCATE( LPIV, STAT = ISTAT1 )
            IF ( istat1 /= 0 )  &
                 WRITE(*, *) ' Error in deallocating LPIV in LA_DGESV'
         END IF
      END IF

      CALL ERINFO( LINFO, SRNAME, INFO, ISTAT )
!
RETURN

CONTAINS


SUBROUTINE dgesv( n, nrhs, a, lda, ipiv, b, ldb, info )

!  -- LAPACK driver routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997

USE dge, ONLY: dgetrf, dgetrs
IMPLICIT NONE

!     .. Scalar Arguments ..
INTEGER, INTENT(IN)       :: lda, ldb, n, nrhs
INTEGER, INTENT(OUT)      :: info
!     ..
!     .. Array Arguments ..
INTEGER, INTENT(OUT)      :: ipiv(:)
REAL (wp), INTENT(IN OUT) :: a(:,:), b(:,:)
!     ..

!  Purpose
!  =======

!  DGESV computes the solution to a real system of linear equations
!     A * X = B,
!  where A is an N by N matrix and X and B are N by NRHS matrices.

!  The LU decomposition with partial pivoting and row interchanges is
!  used to factor A as
!     A = P * L * U,
!  where P is a permutation matrix, L is unit lower triangular, and U is
!  upper triangular.  The factored form of A is then used to solve the
!  system of equations A * X = B.

!  Arguments
!  =========

!  N       (input) INTEGER
!          The number of linear equations, i.e., the order of the
!          matrix A.  N >= 0.

!  NRHS    (input) INTEGER
!          The number of right hand sides, i.e., the number of columns
!          of the matrix B.  NRHS >= 0.

!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the N by N matrix of coefficients A.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.

!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).

!  IPIV    (output) INTEGER array, dimension (N)
!          The pivot indices that define the permutation matrix P;
!          row i of the matrix was interchanged with row IPIV(i).

!  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
!          On entry, the N by NRHS matrix of right hand side vectors B
!          for the system of equations A*X = B.
!          On exit, if INFO = 0, the N by NRHS matrix of solution vectors X.

!  LDB     (input) INTEGER
!          The leading dimension of the array B.  LDB >= max(1,N).

!  INFO    (output) INTEGER
!          = 0: successful exit
!          < 0: if INFO = -k, the k-th argument had an illegal value
!          > 0: if INFO = k, U(k,k) is exactly zero.  The factorization
!               has been completed, but the factor U is exactly
!               singular, so the solution could not be computed.

!  =====================================================================

!     .. External Subroutines ..
! EXTERNAL           dgetrf, dgetrs, xerbla
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC          MAX
!     ..
!     .. Executable Statements ..

!     Test the input parameters.

info = 0
IF( n < 0 ) THEN
  info = -1
ELSE IF( nrhs < 0 ) THEN
  info = -2
ELSE IF( lda < MAX( 1, n ) ) THEN
  info = -4
ELSE IF( ldb < MAX( 1, n ) ) THEN
  info = -7
END IF
IF( info /= 0 ) THEN
  CALL erinfo(-info, 'DGESV ')
  RETURN
END IF

!     Compute the LU factorization of A.

CALL dgetrf( n, n, a, lda, ipiv, info )
IF( info == 0 ) THEN

!        Solve the system A*X = B, overwriting B with X.

  CALL dgetrs( 'No transpose', n, nrhs, a, lda, ipiv, b, ldb,info )
END IF
RETURN

!     End of DGESV

END SUBROUTINE dgesv



END SUBROUTINE DGESV_F90


MODULE dsy
! This module contains translations to ELF90 compatability of the
! LAPACK SYmmetric routines, including DSY (symmetric), DST (symmetric
! tri-diagonal) and DPO (symmetric positive definite) routines.

! Translated by Alan Miller
! alan @ vic.cmis.csiro.au    URL: www.ozemail.com.au/~milleraj

! Latest revision - 3 April 1998

USE la_precision, ONLY: wp => dp
USE la_auxmod
USE dblas
USE dblas2, ONLY: dsyr2, dsymv
USE dblas3, ONLY: dsyr2k
USE dla,    ONLY: dlarfg, dlaev2, dlatrd, dlae2, dlartg, dlasr, dlanst, dlapy2, &
		  dlaset
IMPLICIT NONE


CONTAINS

SUBROUTINE dsytd2( uplo, n, a, lda, d, e, tau, info )

!  -- LAPACK routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997

!     .. Scalar Arguments ..
CHARACTER (LEN=1), INTENT(IN) :: uplo
INTEGER, INTENT(IN)           :: lda, n
INTEGER, INTENT(OUT)          :: info
!     ..
!     .. Array Arguments ..
REAL (wp), INTENT(IN OUT)     :: a( :, : )
REAL (wp), INTENT(OUT)        :: d( : ), e( : ), tau( : )
!     ..

!  Purpose
!  =======

!  DSYTD2 reduces a real symmetric matrix A to symmetric tridiagonal
!  form T by an orthogonal similarity transformation: Q' * A * Q = T.

!  Arguments
!  =========

!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored:
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular

!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.

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

!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).

!  D       (output) DOUBLE PRECISION array, dimension (N)
!          The diagonal elements of the tridiagonal matrix T:
!          D(i) = A(i,i).

!  E       (output) DOUBLE PRECISION array, dimension (N-1)
!          The off-diagonal elements of the tridiagonal matrix T:
!          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.

!  TAU     (output) DOUBLE PRECISION array, dimension (N)
!          The scalar factors of the elementary reflectors (see Further
!          Details).

!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value.

!  Further Details
!  ===============

!  If UPLO = 'U', the matrix Q is represented as a product of elementary
!  reflectors

!     Q = H(n-1) . . . H(2) H(1).

!  Each H(i) has the form

!     H(i) = I - tau * v * v'

!  where tau is a real scalar, and v is a real vector with
!  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
!  A(1:i-1,i+1), and tau in TAU(i).

!  If UPLO = 'L', the matrix Q is represented as a product of elementary
!  reflectors

!     Q = H(1) H(2) . . . H(n-1).

!  Each H(i) has the form

!     H(i) = I - tau * v * v'

!  where tau is a real scalar, and v is a real vector with
!  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
!  and tau in TAU(i).

!  The contents of A on exit are illustrated by the following examples
!  with n = 5:

!  if UPLO = 'U':                       if UPLO = 'L':

!    (  d   e   v2  v3  v4 )              (  d                  )
!    (      d   e   v3  v4 )              (  e   d              )
!    (          d   e   v4 )              (  v1  e   d          )
!    (              d   e  )              (  v1  v2  e   d      )
!    (                  d  )              (  v1  v2  v3  e   d  )

!  where d and e denote diagonal and off-diagonal elements of T, and vi
!  denotes an element of the vector defining H(i).

!  =====================================================================

!     .. Parameters ..
REAL (wp), PARAMETER :: half = 0.5D0
!     ..
!     .. Local Scalars ..
LOGICAL   :: upper
INTEGER   :: i
REAL (wp) :: alpha, taui
!     ..
!     .. External Subroutines ..
! EXTERNAL           daxpy, dlarfg, dsymv, dsyr2, xerbla
!     ..
!     .. External Functions ..
! LOGICAL ::            lsame
! DOUBLE PRECISION ::   ddot
! EXTERNAL           lsame, ddot
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC          MAX, MIN
!     ..
!     .. Executable Statements ..

!     Test the input parameters

info = 0
upper = lsame( uplo, 'U' )
IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
  info = -1
ELSE IF( n < 0 ) THEN
  info = -2
ELSE IF( lda < MAX( 1, n ) ) THEN
  info = -4
END IF
IF( info /= 0 ) THEN
  CALL erinfo(-info, 'DSYTD2')
  RETURN
END IF

!     Quick return if possible

IF( n <= 0 ) RETURN

IF( upper ) THEN

!        Reduce the upper triangle of A

  DO i = n - 1, 1, -1

!           Generate elementary reflector H(i) = I - tau * v * v'
!           to annihilate A(1:i-1,i+1)

    CALL dlarfg( i, a( i, i+1 ), a( :, i+1 ), 1, taui )
    e( i ) = a( i, i+1 )

    IF( taui /= zero ) THEN

!              Apply H(i) from both sides to A(1:i,1:i)

      a( i, i+1 ) = one

!              Compute  x := tau * A * v  storing x in TAU(1:i)

      CALL dsymv( uplo, i, taui, a, lda, a( :, i+1 ), 1, zero, tau, 1 )

!              Compute  w := x - 1/2 * tau * (x'*v) * v

      alpha = -half*taui*ddot( i, tau, 1, a( :, i+1 ), 1 )
      CALL daxpy( i, alpha, a( :, i+1 ), 1, tau, 1 )

!              Apply the transformation as a rank-2 update:
!                 A := A - v * w' - w * v'

      CALL dsyr2( uplo, i, -one, a( :, i+1 ), 1, tau, 1, a, lda )

      a( i, i+1 ) = e( i )
    END IF
    d( i+1 ) = a( i+1, i+1 )
    tau( i ) = taui
  END DO
  d( 1 ) = a( 1, 1 )
ELSE

!        Reduce the lower triangle of A

  DO i = 1, n - 1

!           Generate elementary reflector H(i) = I - tau * v * v'
!           to annihilate A(i+2:n,i)

    CALL dlarfg( n-i, a( i+1, i ), a( MIN( i+2, n ):, i ), 1, taui )
    e( i ) = a( i+1, i )

    IF( taui /= zero ) THEN

!              Apply H(i) from both sides to A(i+1:n,i+1:n)

      a( i+1, i ) = one

!              Compute  x := tau * A * v  storing y in TAU(i+1:n)

      CALL dsymv( uplo, n-i, taui, a( i+1:, i+1: ), lda,  &
                  a( i+1:, i ), 1, zero, tau( i+1: ), 1 )

!              Compute  w := x - 1/2 * tau * (x'*v) * v

      alpha = -half*taui*ddot( n-i, tau( i+1: ), 1, a( i+1:, i ),1 )
      CALL daxpy( n-i, alpha, a( i+1:, i ), 1, tau( i+1: ), 1 )

!              Apply the transformation as a rank-2 update:
!                 A := A - v * w' - w * v'

      CALL dsyr2( uplo, n-i, -one, a( i+1:, i ), 1, tau( i+1: ),  &
                  1, a( i+1:, i+1: ), lda )

      a( i+1, i ) = e( i )
    END IF
    d( i ) = a( i, i )
    tau( i ) = taui
  END DO
  d( n ) = a( n, n )
END IF

RETURN

!     End of DSYTD2

END SUBROUTINE dsytd2


SUBROUTINE dsytrd( uplo, n, a, lda, d, e, tau, info )

!  -- LAPACK routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997
! N.B. Arguments WORK & LWORK have been removed.

!     .. Scalar Arguments ..
CHARACTER (LEN=1), INTENT(IN) :: uplo
INTEGER, INTENT(IN)           :: lda, n
INTEGER, INTENT(OUT)          :: info
!     ..
!     .. Array Arguments ..
REAL (wp), INTENT(IN OUT)     :: a( :, : )
REAL (wp), INTENT(OUT)        :: d( : ), e( : ), tau( : )
!     ..

!  Purpose
!  =======

!  DSYTRD reduces a real symmetric matrix A to symmetric tridiagonal
!  form T by an orthogonal similarity transformation: Q' * A * Q = T.

!  Arguments
!  =========

!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored:
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular

!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.

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

!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).

!  D       (output) DOUBLE PRECISION array, dimension (N)
!          The diagonal elements of the tridiagonal matrix T:
!          D(i) = A(i,i).

!  E       (output) DOUBLE PRECISION array, dimension (N-1)
!          The off-diagonal elements of the tridiagonal matrix T:
!          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.

!  TAU     (output) DOUBLE PRECISION array, dimension (N)
!          The scalar factors of the elementary reflectors (see Further
!          Details).

!  INFO    (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.

!  Further Details
!  ===============

!  If UPLO = 'U', the matrix Q is represented as a product of elementary
!  reflectors

!     Q = H(n-1) . . . H(2) H(1).

!  Each H(i) has the form

!     H(i) = I - tau * v * v'

!  where tau is a real scalar, and v is a real vector with
!  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
!  A(1:i-1,i+1), and tau in TAU(i).

!  If UPLO = 'L', the matrix Q is represented as a product of elementary
!  reflectors

!     Q = H(1) H(2) . . . H(n-1).

!  Each H(i) has the form

!     H(i) = I - tau * v * v'

!  where tau is a real scalar, and v is a real vector with
!  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
!  and tau in TAU(i).

!  The contents of A on exit are illustrated by the following examples
!  with n = 5:

!  if UPLO = 'U':                       if UPLO = 'L':

!    (  d   e   v2  v3  v4 )              (  d                  )
!    (      d   e   v3  v4 )              (  e   d              )
!    (          d   e   v4 )              (  v1  e   d          )
!    (              d   e  )              (  v1  v2  e   d      )
!    (                  d  )              (  v1  v2  v3  e   d  )

!  where d and e denote diagonal and off-diagonal elements of T, and vi
!  denotes an element of the vector defining H(i).

!  =====================================================================

!     .. Local Scalars ..
LOGICAL :: upper
INTEGER :: i, j, kk, ldwork, nb, nx
!     ..
!     .. Local allocatable arrays ..
REAL (wp), ALLOCATABLE, save :: work(:,:)
INTEGER,save :: workd1=-1,workd2=-1
!     ..
!     .. External Subroutines ..
! EXTERNAL           dlatrd, dsyr2k, dsytd2, xerbla
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC          MAX
!     ..
!     .. External Functions ..
! LOGICAL ::            lsame
! INTEGER ::            ilaenv
! EXTERNAL           lsame, ilaenv
!     ..
!     .. Executable Statements ..

!     Test the input parameters

info = 0
upper = lsame( uplo, 'U' )
IF( .NOT.upper .AND. .NOT.lsame( uplo, 'L' ) ) THEN
  info = -1
ELSE IF( n < 0 ) THEN
  info = -2
ELSE IF( lda < MAX( 1, n ) ) THEN
  info = -4
END IF
IF( info /= 0 ) THEN
  CALL erinfo(-info, 'DSYTRD')
  RETURN
END IF

!     Quick return if possible

IF( n == 0 ) THEN
  RETURN
END IF

!     Determine the block size.

nb = ilaenv( 1, 'DSYTRD', uplo, n, -1, -1, -1 )
nx = n
IF( nb > 1 .AND. nb < n ) THEN

!        Determine when to cross over from blocked to unblocked code
!        (last block is always handled by unblocked code).

  nx = MAX( nb, ilaenv( 3, 'DSYTRD', uplo, n, -1, -1, -1 ) )
ELSE
  nb = 1
END IF

ldwork = n
if (ldwork.ne.workd1.or.nb.ne.workd2) then
! repeated allocates are expensive so only do if essential
  if (allocated(work)) deallocate(work)
  ALLOCATE( work(ldwork,nb) )
  workd1=ldwork
  workd2=nb
endif

IF( upper ) THEN

!        Reduce the upper triangle of A.
!        Columns 1:kk are handled by the unblocked method.

  kk = n - ( ( n-nx+nb-1 ) / nb )*nb
  DO i = n - nb + 1, kk + 1, -nb

!           Reduce columns i:i+nb-1 to tridiagonal form and form the
!           matrix W which is needed to update the unreduced part of
!           the matrix

    CALL dlatrd( uplo, i+nb-1, nb, a, lda, e, tau, work, ldwork )

!           Update the unreduced submatrix A(1:i-1,1:i-1), using an
!           update of the form:  A := A - V*W' - W*V'

    CALL dsyr2k( uplo, 'No transpose', i-1, nb, -one, a( :, i: ),  &
                 lda, work, ldwork, one, a, lda )

!           Copy superdiagonal elements back into A, and diagonal
!           elements into D

    DO j = i, i + nb - 1
      a( j-1, j ) = e( j-1 )
      d( j ) = a( j, j )
    END DO
  END DO

!        Use unblocked code to reduce the last or only block

  CALL dsytd2( uplo, kk, a, lda, d, e, tau, info )
ELSE

!        Reduce the lower triangle of A

  DO i = 1, n - nx, nb

!           Reduce columns i:i+nb-1 to tridiagonal form and form the
!           matrix W which is needed to update the unreduced part of
!           the matrix

    CALL dlatrd( uplo, n-i+1, nb, a( i:, i: ), lda, e( i: ),  &
                 tau( i: ), work, ldwork )

!           Update the unreduced submatrix A(i+ib:n,i+ib:n), using
!           an update of the form:  A := A - V*W' - W*V'

    CALL dsyr2k( uplo, 'No transpose', n-i-nb+1, nb, -one, a( i+nb:, i: ),  &
                 lda, work( nb+1:, : ), ldwork, one, a( i+nb:, i+nb: ), lda )

!           Copy subdiagonal elements back into A, and diagonal
!           elements into D

    DO j = i, i + nb - 1
      a( j+1, j ) = e( j )
      d( j ) = a( j, j )
    END DO
  END DO

!        Use unblocked code to reduce the last or only block

  CALL dsytd2( uplo, n-i+1, a( i:, i: ), lda, d( i: ), e( i: ), tau( i: ), &
               info )
END IF

! dont deallocate to avoid expensive allocate/deallocate
!DEALLOCATE( work )

RETURN

!     End of DSYTRD

END SUBROUTINE dsytrd




SUBROUTINE dsteqr( compz, n, d, e, z, ldz, info )

!  -- LAPACK routine (version 2.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     September 30, 1994

! ELF90 translation by Alan Miller   31-Aug-1997
! N.B. Argument WORK has been removed.

!     .. Scalar Arguments ..
CHARACTER (LEN=1), INTENT(IN) :: compz
INTEGER, INTENT(IN)           :: ldz, n
INTEGER, INTENT(OUT)          :: info
!     ..
!     .. Array Arguments ..
REAL (wp), INTENT(IN OUT)     :: d( : ), e( : ), z( :, : )
!     ..

!  Purpose
!  =======

!  DSTEQR computes all eigenvalues and, optionally, eigenvectors of a
!  symmetric tridiagonal matrix using the implicit QL or QR method.
!  The eigenvectors of a full or band symmetric matrix can also be found
!  if DSYTRD or DSPTRD or DSBTRD has been used to reduce this matrix to
!  tridiagonal form.

!  Arguments
!  =========

!  COMPZ   (input) CHARACTER*1
!          = 'N':  Compute eigenvalues only.
!          = 'V':  Compute eigenvalues and eigenvectors of the original
!                  symmetric matrix.  On entry, Z must contain the orthogonal
!                  matrix used to reduce the original matrix to tridiagonal
!                  form.
!          = 'I':  Compute eigenvalues and eigenvectors of the tridiagonal
!                  matrix.  Z is initialized to the identity matrix.

!  N       (input) INTEGER
!          The order of the matrix.  N >= 0.

!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, the diagonal elements of the tridiagonal matrix.
!          On exit, if INFO = 0, the eigenvalues in ascending order.

!  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
!          On entry, the (n-1) subdiagonal elements of the tridiagonal matrix.
!          On exit, E has been destroyed.

!  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)
!          On entry, if  COMPZ = 'V', then Z contains the orthogonal
!          matrix used in the reduction to tridiagonal form.
!          On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the
!          orthonormal eigenvectors of the original symmetric matrix,
!          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
!          of the symmetric tridiagonal matrix.
!          If COMPZ = 'N', then Z is not referenced.

!  LDZ     (input) INTEGER
!          The leading dimension of the array Z.  LDZ >= 1, and if
!          eigenvectors are desired, then  LDZ >= max(1,N).

!  WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,2*N-2))
!          If COMPZ = 'N', then WORK is not referenced.

!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  the algorithm has failed to find all the eigenvalues in a
!                total of 30*N iterations; if INFO = i, then i elements of E
!                have not converged to zero; on exit, D and E contain the
!                elements of a symmetric tridiagonal matrix which is
!                orthogonally similar to the original matrix.

!  =====================================================================

!     .. Parameters ..
REAL (wp), PARAMETER :: two = 2.0D0, three = 3.0D0
INTEGER, PARAMETER   :: maxit = 30
!     ..
!     .. Local Scalars ..
INTEGER :: i, icompz, ii, iscale, j, jtot, k, l, l1, lend,  &
           lendm1, lendp1, lendsv, lm1, lsv, m, mm, mm1,nm1, nmaxit
REAL (wp) :: anorm, b, c, eps, eps2, f, g, p, r, rt1, rt2,  &
             s, safmax, safmin, ssfmax, ssfmin, tst
REAL (wp), ALLOCATABLE, save :: work(:)
INTEGER, save :: workd1=-1
!     ..
!     .. External Functions ..
! LOGICAL ::            lsame
! DOUBLE PRECISION ::   dlamch, dlanst, dlapy2
! EXTERNAL           lsame, dlamch, dlanst, dlapy2
!     ..
!     .. External Subroutines ..
! EXTERNAL           dlae2, dlaev2, dlartg, dlascl, dlaset, dlasr,  &
!                    dlasrt, dswap, xerbla
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC          ABS, MAX, SIGN, SQRT
!     ..
!     .. Executable Statements ..

!     Test the input parameters.

info = 0

IF( lsame( compz, 'N' ) ) THEN
  icompz = 0
ELSE IF( lsame( compz, 'V' ) ) THEN
  icompz = 1
ELSE IF( lsame( compz, 'I' ) ) THEN
  icompz = 2
ELSE
  icompz = -1
END IF

IF( icompz < 0 ) THEN
  info = -1
ELSE IF( n < 0 ) THEN
  info = -2
ELSE IF( ( ldz < 1 ) .OR. ( icompz > 0 .AND. ldz < MAX( 1, n ) ) ) THEN
  info = -6
END IF
IF( info /= 0 ) THEN
  CALL erinfo(-info, 'DSTEQR')
  RETURN
END IF

!     Quick return if possible

IF( n == 0 ) RETURN

IF( n == 1 ) THEN
  IF( icompz == 2 ) z( 1, 1 ) = one
  RETURN
END IF

IF( icompz > 0) THEN
  if (n.ne.workd1) then
    if (allocated(work))deallocate(work)
    ALLOCATE( work( MAX(1, 2*n-2) ) )
    workd1=n
  endif
END IF  

!     Determine the unit roundoff and over/underflow thresholds.

eps = dlamch( 'E' )
eps2 = eps**2
safmin = dlamch( 'S' )
safmax = one / safmin
ssfmax = SQRT( safmax ) / three
ssfmin = SQRT( safmin ) / eps2

!     Compute the eigenvalues and eigenvectors of the tridiagonal
!     matrix.

IF( icompz == 2 ) CALL dlaset( 'Full', n, n, zero, one, z )

nmaxit = n*maxit
jtot = 0

!     Determine where the matrix splits and choose QL or QR iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.

l1 = 1
nm1 = n - 1

10 IF( l1 > n ) GO TO 160
IF( l1 > 1 ) e( l1-1 ) = zero
IF( l1 <= nm1 ) THEN
  DO m = l1, nm1
    tst = ABS( e( m ) )
    IF( tst == zero ) GO TO 30
    IF( tst <= ( SQRT( ABS( d( m ) ) )*SQRT( ABS( d( m+ 1 ) ) ) )*eps ) THEN
      e( m ) = zero
      GO TO 30
    END IF
  END DO
END IF
m = n

30 l = l1
lsv = l
lend = m
lendsv = lend
l1 = m + 1
IF( lend == l ) GO TO 10

!     Scale submatrix in rows and columns L to LEND

anorm = dlanst( 'I', lend-l+1, d( l: ), e( l: ) )
iscale = 0
IF( anorm == zero ) GO TO 10
IF( anorm > ssfmax ) THEN
  iscale = 1
!  CALL dlascl( 'G', 0, 0, anorm, ssfmax, lend-l+1, 1, d( l: ), n, info )
  d( l:lend ) = (ssfmax / anorm) * d( l:lend )
!  CALL dlascl( 'G', 0, 0, anorm, ssfmax, lend-l, 1, e( l: ), n, info )
  e( l:lend-1 ) = (ssfmax / anorm) * e( l:lend-1 )

ELSE IF( anorm < ssfmin ) THEN
  iscale = 2
!  CALL dlascl( 'G', 0, 0, anorm, ssfmin, lend-l+1, 1, d( l: ), n, info )
  d( l:lend ) = (ssfmin / anorm) * d( l:lend )
!  CALL dlascl( 'G', 0, 0, anorm, ssfmin, lend-l, 1, e( l: ), n, info )
  e( l:lend-1 ) = (ssfmin / anorm) * e( l:lend-1 )
END IF

!     Choose between QL and QR iteration

IF( ABS( d( lend ) ) < ABS( d( l ) ) ) THEN
  lend = lsv
  l = lendsv
END IF

IF( lend > l ) THEN

!        QL Iteration

!        Look for small subdiagonal element.

  40 IF( l /= lend ) THEN
    lendm1 = lend - 1
    DO m = l, lendm1
      tst = ABS( e( m ) )**2
      IF( tst <= ( eps2*ABS( d( m ) ) )*ABS( d( m+1 ) ) + safmin ) GO TO 60
    END DO
  END IF

  m = lend

  60 IF( m < lend ) e( m ) = zero
  p = d( l )
  IF( m == l ) GO TO 80

!        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
!        to compute its eigensystem.

  IF( m == l+1 ) THEN
    IF( icompz > 0 ) THEN
      CALL dlaev2( d( l ), e( l ), d( l+1 ), rt1, rt2, c, s )
      work( l ) = c
      work( n-1+l ) = s
      CALL dlasr( 'R', 'V', 'B', n, 2, work( l: ),  &
                  work( n-1+l: ), z( :, l: ), ldz )
    ELSE
      CALL dlae2( d( l ), e( l ), d( l+1 ), rt1, rt2 )
    END IF
    d( l ) = rt1
    d( l+1 ) = rt2
    e( l ) = zero
    l = l + 2
    IF( l <= lend ) GO TO 40
    GO TO 140
  END IF

  IF( jtot == nmaxit ) GO TO 140
  jtot = jtot + 1

!        Form shift.

  g = ( d( l+1 )-p ) / ( two*e( l ) )
  r = dlapy2( g, one )
  g = d( m ) - p + ( e( l ) / ( g+SIGN( r, g ) ) )

  s = one
  c = one
  p = zero

!        Inner loop

  mm1 = m - 1
  DO i = mm1, l, -1
    f = s*e( i )
    b = c*e( i )
    CALL dlartg( g, f, c, s, r )
    IF( i /= m-1 ) e( i+1 ) = r
    g = d( i+1 ) - p
    r = ( d( i )-g )*s + two*c*b
    p = s*r
    d( i+1 ) = g + p
    g = c*r - b

!           If eigenvectors are desired, then save rotations.

    IF( icompz > 0 ) THEN
      work( i ) = c
      work( n-1+i ) = -s
    END IF

  END DO

!        If eigenvectors are desired, then apply saved rotations.

  IF( icompz > 0 ) THEN
    mm = m - l + 1
    CALL dlasr( 'R', 'V', 'B', n, mm, work( l: ), work( n-1+l: ), &
                z( :, l: ), ldz )
  END IF

  d( l ) = d( l ) - p
  e( l ) = g
  GO TO 40

!        Eigenvalue found.

  80 d( l ) = p

  l = l + 1
  IF( l <= lend ) GO TO 40
  GO TO 140

ELSE

!        QR Iteration

!        Look for small superdiagonal element.

  90 IF( l /= lend ) THEN
    lendp1 = lend + 1
    DO m = l, lendp1, -1
      tst = ABS( e( m-1 ) )**2
      IF( tst <= ( eps2*ABS( d( m ) ) )*ABS( d( m-1 ) )+ safmin )GO TO 110
    END DO
  END IF

  m = lend

  110 IF( m > lend ) e( m-1 ) = zero
  p = d( l )
  IF( m == l ) GO TO 130

!        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
!        to compute its eigensystem.

  IF( m == l-1 ) THEN
    IF( icompz > 0 ) THEN
      CALL dlaev2( d( l-1 ), e( l-1 ), d( l ), rt1, rt2, c, s )
      work( m ) = c
      work( n-1+m ) = s
      CALL dlasr( 'R', 'V', 'F', n, 2, work( m: ),  &
                  work( n-1+m: ), z( :, l-1: ), ldz )
    ELSE
      CALL dlae2( d( l-1 ), e( l-1 ), d( l ), rt1, rt2 )
    END IF
    d( l-1 ) = rt1
    d( l ) = rt2
    e( l-1 ) = zero
    l = l - 2
    IF( l >= lend ) GO TO 90
    GO TO 140
  END IF

  IF( jtot == nmaxit ) GO TO 140
  jtot = jtot + 1

!        Form shift.

  g = ( d( l-1 )-p ) / ( two*e( l-1 ) )
  r = dlapy2( g, one )
  g = d( m ) - p + ( e( l-1 ) / ( g+SIGN( r, g ) ) )

  s = one
  c = one
  p = zero

!        Inner loop

  lm1 = l - 1
  DO i = m, lm1
    f = s*e( i )
    b = c*e( i )
    CALL dlartg( g, f, c, s, r )
    IF( i /= m ) e( i-1 ) = r
    g = d( i ) - p
    r = ( d( i+1 )-g )*s + two*c*b
    p = s*r
    d( i ) = g + p
    g = c*r - b

!           If eigenvectors are desired, then save rotations.

    IF( icompz > 0 ) THEN
      work( i ) = c
      work( n-1+i ) = s
    END IF

  END DO

!        If eigenvectors are desired, then apply saved rotations.

  IF( icompz > 0 ) THEN
    mm = l - m + 1
    CALL dlasr( 'R', 'V', 'F', n, mm, work( m: ), work( n-1+m: ),  &
                z( :, m: ), ldz )
  END IF

  d( l ) = d( l ) - p
  e( lm1 ) = g
  GO TO 90

!        Eigenvalue found.

  130 d( l ) = p

  l = l - 1
  IF( l >= lend ) GO TO 90
  GO TO 140

END IF

!     Undo scaling if necessary

140 IF( iscale == 1 ) THEN
!  CALL dlascl( 'G', 0, 0, ssfmax, anorm, lendsv-lsv+1, 1, d( lsv: ), n, info )
   d( lsv:lendsv ) = (anorm / ssfmax) * d( lsv:lendsv )
!  CALL dlascl( 'G', 0, 0, ssfmax, anorm, lendsv-lsv, 1, e( lsv: ), n, info )
   e( lsv:lendsv-1 ) = (anorm / ssfmax) * e( lsv:lendsv-1 )

ELSE IF( iscale == 2 ) THEN
!  CALL dlascl( 'G', 0, 0, ssfmin, anorm, lendsv-lsv+1, 1, d( lsv: ), n, info )
   d( lsv:lendsv ) = (anorm / ssfmin) * d( lsv:lendsv )
!  CALL dlascl( 'G', 0, 0, ssfmin, anorm, lendsv-lsv, 1, e( lsv: ), n, info )
   e( lsv:lendsv-1 ) = (anorm / ssfmin) * e( lsv:lendsv-1 )
END IF

!     Check for no convergence to an eigenvalue after a total
!     of N*MAXIT iterations.

IF( jtot < nmaxit ) GO TO 10
DO i = 1, n - 1
  IF( e( i ) /= zero ) info = info + 1
END DO
GO TO 190

!     Order eigenvalues and eigenvectors.

160 IF( icompz == 0 ) THEN

!        Use Quick Sort

  CALL dlasrt( 'I', n, d, info )

ELSE

!        Use Selection Sort to minimize swaps of eigenvectors

  DO ii = 2, n
    i = ii - 1
    k = i
    p = d( i )
    DO j = ii, n
      IF( d( j ) < p ) THEN
        k = j
        p = d( j )
      END IF
    END DO
    IF( k /= i ) THEN
      d( k ) = d( i )
      d( i ) = p
      CALL dswap( n, z( :, i ), 1, z( :, k ), 1 )
    END IF
  END DO

!  DEALLOCATE( work )
END IF

190 RETURN

!     End of DSTEQR

END SUBROUTINE dsteqr




SUBROUTINE dsterf( n, d, e, info )

!  -- LAPACK routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997

!     .. Scalar Arguments ..
INTEGER, INTENT(IN)       :: n
INTEGER, INTENT(OUT)      :: info
!     ..
!     .. Array Arguments ..
REAL (wp), INTENT(IN OUT) :: d( : ), e( : )
!     ..

!  Purpose
!  =======

!  DSTERF computes all eigenvalues of a symmetric tridiagonal matrix
!  using the Pal-Walker-Kahan variant of the QL or QR algorithm.

!  Arguments
!  =========

!  N       (input) INTEGER
!          The number of rows and columns in the matrix.  N >= 0.

!  D       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, D contains the diagonal elements of the
!          tridiagonal matrix.
!          On exit, D contains the eigenvalues in ascending order.
!          If an error exit is made, the eigenvalues are correct
!          but unordered for indices 1,2,...,INFO-1.

!  E       (input/output) DOUBLE PRECISION array, dimension (N)
!          On entry, E contains the subdiagonal elements of the
!          tridiagonal matrix in positions 1 through N-1.
!          E(N) is arbitrary.
!          On exit, E has been destroyed.

!  INFO    (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if INFO = +i, the i-th eigenvalue has not converged
!                after a total of  30*N  iterations.

!     .. Parameters ..
REAL (wp), PARAMETER :: two = 2.0D0
INTEGER, PARAMETER   :: maxit = 30
!     ..
!     .. Local Scalars ..
INTEGER :: i, ii, j, jtot, k, l, l1, lend, lendm1, lendp1,  &
           lm1, m, mm1, nconv, nm1, nmaxit
REAL (wp) :: alpha, bb, c, eps, gamma, oldc, oldgam, p, r,  &
             rt1, rt2, rte, s, sigma, tst
!     ..
!     .. External Functions ..
! DOUBLE PRECISION ::   dlamch, dlapy2
! EXTERNAL           dlamch, dlapy2
!     ..
!     .. External Subroutines ..
! EXTERNAL           dlae2, xerbla
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC          ABS, SIGN, SQRT
!     ..
!     .. Executable Statements ..

!     Test the input parameters.

info = 0

!     Quick return if possible

IF( n < 0 ) THEN
  info = -1
  CALL erinfo(-info, 'DSTERF')
  RETURN
END IF
IF( n <= 1 ) RETURN

!     Determine the unit roundoff for this environment.

eps = dlamch( 'E' )

!     Compute the eigenvalues of the tridiagonal matrix.

DO i = 1, n - 1
  e( i ) = e( i )**2
END DO
e( n ) = zero

nmaxit = n*maxit
sigma = zero
jtot = 0
nconv = 0

!     Determine where the matrix splits and choose QL or QR iteration
!     for each block, according to whether top or bottom diagonal
!     element is smaller.

l1 = 1
nm1 = n - 1

20 IF( l1 > n ) GO TO 160
IF( l1 <= nm1 ) THEN
  DO m = l1, nm1
    tst = SQRT( ABS( e( m ) ) )
    IF( tst <= eps*( ABS( d( m ) )+ABS( d( m+1 ) ) ) ) GO TO 40
  END DO
END IF
m = n

40 l = l1
lend = m
IF( ABS( d( lend ) ) < ABS( d( l ) ) ) THEN
  l = lend
  lend = l1
END IF
l1 = m + 1

IF( lend >= l ) THEN

!        QL Iteration

!        Look for small subdiagonal element.

  50 IF( l /= lend ) THEN
    lendm1 = lend - 1
    DO m = l, lendm1
      tst = SQRT( ABS( e( m ) ) )
      IF( tst <= eps*( ABS( d( m ) )+ABS( d( m+1 ) ) ) ) GO TO 70
    END DO
  END IF

  m = lend

  70 p = d( l )
  IF( m == l ) GO TO 90

!        If remaining matrix is 2 by 2, use DLAE2 to compute its
!        eigenvalues.

  IF( m == l+1 ) THEN
    rte = SQRT( e( l ) )
    CALL dlae2( d( l ), rte, d( l+1 ), rt1, rt2 )
    d( l ) = rt1
    d( l+1 ) = rt2
    nconv = nconv + 2
    l = l + 2
    IF( l <= lend ) GO TO 50
    GO TO 20
  END IF

  IF( jtot == nmaxit ) GO TO 150
  jtot = jtot + 1

!        Form shift.

  rte = SQRT( e( l ) )
  sigma = ( d( l+1 )-p ) / ( two*rte )
  r = dlapy2( sigma, one )
  sigma = p - ( rte / ( sigma+SIGN( r, sigma ) ) )

  c = one
  s = zero
  gamma = d( m ) - sigma
  p = gamma*gamma

!        Inner loop

  mm1 = m - 1
  DO i = mm1, l, -1
    bb = e( i )
    r = p + bb
    e( i+1 ) = s*r
    oldc = c
    c = p / r
    s = bb / r
    oldgam = gamma
    alpha = d( i )
    gamma = c*( alpha-sigma ) - s*oldgam
    d( i+1 ) = oldgam + ( alpha-gamma )
    IF( c /= zero ) THEN
      p = ( gamma*gamma ) / c
    ELSE
      p = oldc*bb
    END IF
  END DO

  e( l ) = s*p
  d( l ) = sigma + gamma
  e( m ) = zero
  GO TO 50

!        Eigenvalue found.

  90 d( l ) = p
  nconv = nconv + 1

  l = l + 1
  IF( l <= lend ) GO TO 50
  GO TO 20

ELSE

!        QR Iteration

!        Look for small superdiagonal element.

  100 IF( l /= lend ) THEN
    lendp1 = lend + 1
    DO m = l, lendp1, -1
      tst = SQRT( ABS( e( m-1 ) ) )
      IF( tst <= eps*( ABS( d( m ) )+ABS( d( m-1 ) ) ) ) GO TO 120
    END DO
  END IF

  m = lend

  120 p = d( l )
  IF( m == l ) GO TO 140

!        If remaining matrix is 2 by 2, use DLAE2 to compute its eigenvalues.

  IF( m == l-1 ) THEN
    rte = SQRT( e( l-1 ) )
    CALL dlae2( d( l ), rte, d( l-1 ), rt1, rt2 )
    d( l ) = rt1
    d( l-1 ) = rt2
    nconv = nconv + 2
    l = l - 2
    IF( l >= lend ) GO TO 100
    GO TO 20
  END IF

  IF( jtot == nmaxit ) GO TO 150
  jtot = jtot + 1

!        Form shift.

  rte = SQRT( e( l-1 ) )
  sigma = ( d( l-1 )-p ) / ( two*rte )
  r = dlapy2( sigma, one )
  sigma = p - ( rte / ( sigma+SIGN( r, sigma ) ) )

  c = one
  s = zero
  gamma = d( m ) - sigma
  p = gamma*gamma

!        Inner loop

  lm1 = l - 1
  DO i = m, lm1
    bb = e( i )
    r = p + bb
    IF( i /= 1 ) e( i-1 ) = s*r
    oldc = c
    c = p / r
    s = bb / r
    oldgam = gamma
    alpha = d( i+1 )
    gamma = c*( alpha-sigma ) - s*oldgam
    d( i ) = oldgam + ( alpha-gamma )
    IF( c /= zero ) THEN
      p = ( gamma*gamma ) / c
    ELSE
      p = oldc*bb
    END IF
  END DO

  e( lm1 ) = s*p
  d( l ) = sigma + gamma
  IF( m /= 1 ) e( m-1 ) = zero
  GO TO 100

!        Eigenvalue found.

  140 d( l ) = p
  nconv = nconv + 1

  l = l - 1
  IF( l >= lend ) GO TO 100
  GO TO 20

END IF

!     Set error -- no convergence to an eigenvalue after a total
!     of N*MAXIT iterations.

150 info = nconv
RETURN

!     Sort eigenvalues in increasing order.

160 DO ii = 2, n
  i = ii - 1
  k = i
  p = d( i )
  DO j = ii, n
    IF( d( j ) < p ) THEN
      k = j
      p = d( j )
    END IF
  END DO
  IF( k /= i ) THEN
    d( k ) = d( i )
    d( i ) = p
  END IF
END DO

RETURN

!     End of DSTERF

END SUBROUTINE dsterf

END MODULE dsy
SUBROUTINE DSYEV_F90( A, W, JOBZ, UPLO, INFO )
!
!  -- LAPACK90 interface driver routine (version 1.0) --
!     UNI-C, Denmark; Univ. of Tennessee, USA; NAG Ltd., UK
!     May 31, 1997
!
! ELF90 translation by Alan Miller    10-Nov-1997
! Latest revision - 23 November 1997

!  .. USE STATEMENTS ..
   USE LA_PRECISION, ONLY: WP => DP
   USE LA_AUXMOD
!   USE F77_LAPACK, ONLY: SYEV_F77 => LA_SYEV, ILAENV_F77 => ILAENV

!  .. IMPLICIT STATEMENT ..
   IMPLICIT NONE
!  .. CHARACTER ARGUMENTS ..
   CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: JOBZ, UPLO
!  .. SCALAR ARGUMENTS ..
   INTEGER, INTENT(OUT), OPTIONAL         :: INFO
!  .. ARRAY ARGUMENTS ..
   REAL(WP), INTENT(IN OUT)               :: A(:,:)
   REAL(WP), INTENT(OUT)                  :: W(:)
!-----------------------------------------------------------------
!
! Purpose
! =======
!
! LA_SYEV / LA_HEEV computes all eigenvalues and, optionally,
! eigenvectors of a real symmetric or Hermitian matrix A.
!
! =======
!
!    SUBROUTINE LA_SYEV /  LA_HEEV( A, W, JOBZ, UPLO, INFO )
!       <type>(<wp>), INTENT(INOUT) :: A(:,:)
!       REAL(<wp>), INTENT(OUT) :: W(:)
!       CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: JOBZ, UPLO
!       INTEGER, INTENT(OUT), OPTIONAL :: INFO
!    where
!       <type> ::= REAL | COMPLEX
!       <wp>   ::= KIND(1.0) | KIND(1.0D0)
!
! Defaults
! ========
!
! 1. If JOBZ is not present then JOBZ = 'N' is assumed.
!
! 2. If UPLO is not present then UPLO = 'U' is assumed.
!
! Arguments
! =========
!
! A       (input/output) either REAL/COMPLEX square array,
!         shape (:,:), size(A,1) == size(A,2).
!         On entry, the symmetric (Hermitian) matrix A.
!            If UPLO = 'U', the upper triangular part of A contains
!               the upper triangular part of the matrix A.
!            If UPLO = 'L', the lower triangular part of A contains
!               the lower triangular part of the matrix A.
!         On exit:
!            If JOBZ = 'V', then if INFO = 0, A contains the
!               orthonormal eigenvectors of the matrix A.
!            If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
!               or the upper triangle (if UPLO='U') of A, including the
!               diagonal, is destroyed.
!
! W       (output) REAL array,  shape (:), size(W) == size(A,1) >= 0.
!         If INFO = 0, the eigenvalues in ascending order.
!
! JOBZ    Optional, (input) CHARACTER*1
!         If JOBZ is present then:
!            = 'N':  Compute eigenvalues only;
!            = 'V':  Compute eigenvalues and eigenvectors.
!         otherwise JOBZ = 'N' is assumed.
!
! UPLO    Optional, (input) CHARACTER*1
!         If UPLO is present then:
!            = 'U':  Upper triangle of A is stored;
!            = 'L':  Lower triangle of A is stored.
!         otherwise UPLO = 'U' is assumed.
!
! INFO    Optional, (output) INTEGER
!         If INFO is present:
!            = 0: successful exit
!            < 0: if INFO = -i, the i-th argument had an illegal value
!            > 0: if INFO = i, the algorithm failed to converge; i
!                 off-diagonal elements of an intermediate tridiagonal
!                 form did not converge to zero.
!         If INFO is not present and an error occurs, then the program
!            is terminated with an error message.
!
!-------------------------------------------
!  .. LOCAL PARAMETERS ..
   CHARACTER(LEN=7), PARAMETER :: SRNAME = 'LA_SYEV'
!  .. LOCAL SCALARS ..
   CHARACTER(LEN=1)  :: LJOBZ, LUPLO
   INTEGER           :: N, LINFO, LD, ISTAT
!  .. INTRINSIC FUNCTIONS ..
!   INTRINSIC MAX, PRESENT

!  .. EXECUTABLE STATEMENTS ..
   N = SIZE( A, 1 )
   LINFO = 0
   ISTAT = 0
   LD = MAX(1,N)

   IF( PRESENT(JOBZ) ) THEN
      LJOBZ = JOBZ
   ELSE
      LJOBZ = 'N'
   END IF

   IF( PRESENT(UPLO) ) THEN
      LUPLO = UPLO
   ELSE
      LUPLO = 'U'
   END IF

!  .. TEST THE ARGUMENTS
   IF( SIZE( A, 2 ) /= N .OR. N < 0 )THEN
      LINFO = -1
   ELSE IF( SIZE( W ) /= N )THEN
      LINFO = -2
   ELSE IF( .NOT.LSAME(LJOBZ,'N') .AND. .NOT.LSAME(LJOBZ,'V') )THEN
      LINFO = -3
   ELSE IF( .NOT.LSAME(LUPLO,'U') .AND. .NOT.LSAME(LUPLO,'L') )THEN
      LINFO = -4
   ELSE IF( N > 0 )THEN
!
      IF( LINFO == 0 )THEN
         CALL DSYEV( LJOBZ, LUPLO, N, A, LD, W, LINFO )
      END IF
   END IF
   IF( PRESENT( info ) ) CALL ERINFO(LINFO, SRNAME, INFO, ISTAT)
!
RETURN

CONTAINS


SUBROUTINE dsyev( jobz, uplo, n, a, lda, w, info )

!  -- LAPACK driver routine (version 1.0) --
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
!     Courant Institute, Argonne National Lab, and Rice University
!     February 29, 1992

! ELF90 translation by Alan Miller   31-Aug-1997
! N.B. Arguments WORK & LWORK have been removed

USE dblas
USE dla, ONLY: dlansy
USE dsy, ONLY: dsteqr, dsterf, dsytrd
USE dor, ONLY: dorgtr

!     .. Scalar Arguments ..
CHARACTER (LEN=1), INTENT(IN) :: jobz, uplo
INTEGER, INTENT(IN)           :: lda, n
INTEGER, INTENT(OUT)          :: info
!     ..
!     .. Array Arguments ..
REAL (wp), INTENT(IN OUT)     :: a( :, : )
REAL (wp), INTENT(OUT)        :: w( : )
!     ..

!  Purpose
!  =======

!  DSYEV  computes all eigenvalues and, optionally, eigenvectors of a
!  real symmetric matrix A by calling the recommended sequence of LAPACK
!  routines.

!  Arguments
!  =========

!  JOBZ    (input) CHARACTER*1
!          Specifies whether or not to compute the eigenvectors:
!          = 'N':  Compute eigenvalues only.
!          = 'V':  Compute eigenvalues and eigenvectors.

!  UPLO    (input) CHARACTER*1
!          Specifies whether the upper or lower triangular part of the
!          symmetric matrix A is stored:
!          = 'U':  Upper triangular
!          = 'L':  Lower triangular

!  N       (input) INTEGER
!          The number of rows and columns of the matrix A.  N >= 0.

!  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', only the
!          upper triangular part of A is used to define the elements of
!          the symmetric matrix.  If UPLO = 'L', only the lower
!          triangular part of A is used to define the elements of the
!          symmetric matrix.

!          If JOBZ = 'V', then if INFO = 0 on exit, A contains the
!          orthonormal eigenvectors of the matrix A.  If INFO > 0, A
!          contains the eigenvectors associated with only the stored
!          eigenvalues.
!          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
!          or the upper triangle (if UPLO='U') of A, including the
!          diagonal, is destroyed.

!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).

!  W       (output) DOUBLE PRECISION array, dimension (N)
!          On exit, if INFO = 0, W contains the eigenvalues in ascending
!          order.  If INFO > 0, the eigenvalues are correct for indices
!          1, 2, ..., INFO-1, but they are unordered and may not be the
!          smallest eigenvalues of the matrix.

!  INFO    (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if INFO = +i, the algorithm terminated before finding
!                the i-th eigenvalue.

!     .. Local Scalars ..
LOGICAL   :: lower, wantz
INTEGER   :: imax, inde, indtau, iscale, j
REAL (wp) :: anrm, bignum, eps, rmax, rmin, safmin, sigma, smlnum
!     ..
!     .. Local array
REAL (wp) :: work( 4*n )
!     ..
!     .. External Functions ..
! LOGICAL ::            lsame
! DOUBLE PRECISION ::   dlamch, dlansy
! EXTERNAL           lsame, dlamch, dlansy
!     ..
!     .. External Subroutines ..
! EXTERNAL           dorgtr, dscal, dsteqr, dsterf, dsytrd, xerbla
!     ..
!     .. Intrinsic Functions ..
! INTRINSIC          MAX, SQRT
!     ..
!     .. Executable Statements ..

!     Test the input parameters.

wantz = lsame( jobz, 'V' )
lower = lsame( uplo, 'L' )

info = 0
IF( .NOT.( wantz .OR. lsame( jobz, 'N' ) ) ) THEN
  info = -1
ELSE IF( .NOT.( lower .OR. lsame( uplo, 'U' ) ) ) THEN
  info = -2
ELSE IF( n < 0 ) THEN
  info = -3
ELSE IF( lda < MAX( 1, n ) ) THEN
  info = -5
END IF

IF( info /= 0 ) THEN
  CALL erinfo(-info, 'DSYEV ')
  RETURN
END IF

!     Quick return if possible

IF( n == 0 ) THEN
  RETURN
END IF

IF( n == 1 ) THEN
  w( 1 ) = a( 1, 1 )
  IF( wantz ) a( 1, 1 ) = one
  RETURN
END IF

!     Get machine constants.

safmin = dlamch( 'Safe minimum' )
eps = dlamch( 'Precision' )
smlnum = safmin / eps
bignum = one / smlnum
rmin = SQRT( smlnum )
rmax = SQRT( bignum )

!     Scale matrix to allowable range, if necessary.

anrm = dlansy( 'M', uplo, n, a, lda )
iscale = 0
IF( anrm > zero .AND. anrm < rmin ) THEN
  iscale = 1
  sigma = rmin / anrm
ELSE IF( anrm > rmax ) THEN
  iscale = 1
  sigma = rmax / anrm
END IF
IF( iscale == 1 ) THEN
  IF( lower ) THEN
    DO j = 1, n
      CALL dscal( n-j+1, sigma, a( j:, j ), 1 )
    END DO
  ELSE
    DO j = 1, n
      CALL dscal( j, sigma, a( :, j ), 1 )
    END DO
  END IF
END IF

!     Call DSYTRD to reduce symmetric matrix to tridiagonal form.

inde = 1
indtau = inde + n
CALL dsytrd( uplo, n, a, lda, w, work( inde: ), work( indtau: ), info )

!     For eigenvalues only, call DSTERF.  For eigenvectors, first call
!     DORGTR to generate the orthogonal matrix, then call DSTEQR.

IF( .NOT.wantz ) THEN
  CALL dsterf( n, w, work( inde: ), info )
ELSE
  CALL dorgtr( uplo, n, a, lda, work( indtau: ), info )
  CALL dsteqr( jobz, n, w, work( inde: ), a, lda, info )
END IF

!     If matrix was scaled, then rescale eigenvalues appropriately.

IF( iscale == 1 ) THEN
  IF( info == 0 ) THEN
    imax = n
  ELSE
    imax = info - 1
  END IF
  CALL dscal( imax, one / sigma, w, 1 )
END IF

RETURN

!     End of DSYEV

END SUBROUTINE dsyev


END SUBROUTINE DSYEV_F90


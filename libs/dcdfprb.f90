!  dcdfprb.f90  20 June 2000
!
program dcdfprb
!
!*******************************************************************************
!
!! DCDFPRB calls the DCDFLIB tests.
!
  character ( len = 8 ) date
  character ( len = 10 ) time
!
  call date_and_time ( date, time )

  write ( *, * ) ' '
  write ( *, * ) 'DCDFPRB'
  write ( *, * ) '  Tests for the DCDFLIB package.'
  write ( *, * ) ' '
  write ( *, * ) '  Today''s date: ', date
  write ( *, * ) '  Today''s time: ', time
 
  call test01
  call test02
  call test03
  call test04
  call test05
  call test06
  call test07
  call test08
  call test09
  call test10

  call test11
  call test12
  call test13
  call test14

  write ( *, * ) ' '
  write ( *, * ) 'DCDFPRB'
  write ( *, * ) '  Normal end of DCDFLIB tests.'

  stop
end
subroutine test01
!
!*******************************************************************************
!
!! TEST01 tests CUMBET.
!
  integer, parameter :: test_num = 16
!
  double precision a
  double precision b
  double precision ccum
  double precision cum
  integer i
  double precision x
  double precision y
!
  a = 3.0
  b = 2.0

  write ( *, * ) ' '
  write ( *, * ) 'TEST01'
  write ( *, * ) '  CUMBET computes the cumulative'
  write ( *, * ) '  density function for the beta distribution.'
  write ( *, * ) ' '
  write ( *, * ) '  Parameters:'
  write ( *, * ) ' '
  write ( *, * ) '    A = ', a
  write ( *, * ) '    B = ', b
  write ( *, * ) ' '
  write ( *, * ) '  X  CDF  1-CDF'
  write ( *, * ) ' '

  do i = 0, test_num

    x = real ( i ) / real ( test_num )
    y = 1.0 - x

    call cumbet ( x, y, a, b, cum, ccum )

    write ( *,     '(3g14.6)' ) x, cum, ccum

  end do

  return
end
subroutine test02
!
!*******************************************************************************
!
!! TEST02 tests CUMBIN.
!
  double precision ccum
  double precision cum
  integer i
  double precision ompr
  double precision pr
  double precision s
  double precision xn
!
  pr = 0.75
  ompr = 1.0 - pr
  xn = 10

  write ( *, * ) ' '
  write ( *, * ) 'TEST02'
  write ( *, * ) '  CUMBIN computes the cumulative density'
  write ( *, * ) '  function for the binomial distribution.'
  write ( *, * ) ' '
  write ( *, * ) '  Parameters:'
  write ( *, * ) ' '
  write ( *, * ) '    PR = ', pr
  write ( *, * ) '    XN = ', xn
  write ( *, * ) ' '
  write ( *, * ) '  X  CDF  1-CDF'
  write ( *, * ) ' '

  do i = 0, int ( xn )

    s = real ( i ) 

    call cumbin ( s, xn, pr, ompr, cum, ccum )

    write ( *,     '(3g14.6)' ) s, cum, ccum

  end do

  return
end
subroutine test03
!
!*******************************************************************************
!
!! TEST03 tests CUMCHI.
!
  double precision ccum
  double precision cum
  double precision df
  integer i
  double precision x
!
  df = 4.0

  write ( *, * ) ' '
  write ( *, * ) 'TEST03'
  write ( *, * ) '  CUMCHI computes the cumulative density'
  write ( *, * ) '  function for the chi-squared distribution.'
  write ( *, * ) ' '
  write ( *, * ) '  Parameters:'
  write ( *, * ) ' '
  write ( *, * ) '    DF = ', df
  write ( *, * ) ' '
  write ( *, * ) '  X  CDF  1-CDF'
  write ( *, * ) ' '

  do i = 0, 10

    x = real ( i ) / 5.0

    call cumchi ( x, df, cum, ccum )

    write ( *,     '(3g14.6)' ) x, cum, ccum

  end do

  return
end
subroutine test04
!
!*******************************************************************************
!
!! TEST04 tests CUMCHN.
!
  double precision ccum
  double precision cum
  double precision df
  integer i
  double precision pnonc
  double precision x
!
  df = 4.0
  pnonc = 0.5

  write ( *, * ) ' '
  write ( *, * ) 'TEST04'
  write ( *, * ) '  CUMCHN computes the cumulative density'
  write ( *, * ) '  function for the noncentral chi-squared '
  write ( *, * ) '  distribution.'
  write ( *, * ) ' '
  write ( *, * ) '  Parameters:'
  write ( *, * ) ' '
  write ( *, * ) '    DF = ', df
  write ( *, * ) '    PNONC = ', pnonc
  write ( *, * ) ' '
  write ( *, * ) '  X  CDF  1-CDF'
  write ( *, * ) ' '

  do i = 0, 10

    x = real ( i ) / 5.0

    call cumchn ( x, df, pnonc, cum, ccum )

    write ( *,     '(3g14.6)' ) x, cum, ccum

  end do

  return
end
subroutine test05
!
!*******************************************************************************
!
!! TEST05 tests CUMF.
!
  double precision ccum
  double precision cum
  double precision dfd
  double precision dfn
  integer i
  double precision x
!
  dfd = 4.0
  dfn = 0.5

  write ( *, * ) ' '
  write ( *, * ) 'TEST05'
  write ( *, * ) '  CUMF computes the cumulative density'
  write ( *, * ) '  function for the F distribution.'
  write ( *, * ) ' '
  write ( *, * ) '  Parameters:'
  write ( *, * ) ' '
  write ( *, * ) '    DFN = ', dfn
  write ( *, * ) '    DFD = ', dfd
  write ( *, * ) ' '
  write ( *, * ) '  X  CDF  1-CDF'
  write ( *, * ) ' '

  do i = 0, 10

    x = real ( i ) / 5.0

    call cumf ( x, dfn, dfd, cum, ccum )

    write ( *,     '(3g14.6)' ) x, cum, ccum

  end do

  return
end
subroutine test06
!
!*******************************************************************************
!
!! TEST06 tests CUMFNC.
!
  double precision ccum
  double precision cum
  double precision dfd
  double precision dfn
  integer i
  double precision pnonc
  double precision x
!
  dfd = 4.0
  dfn = 0.5
  pnonc = 0.5

  write ( *, * ) ' '
  write ( *, * ) 'TEST06'
  write ( *, * ) '  CUMFNC computes the cumulative density'
  write ( *, * ) '  function for the noncentral F distribution.'
  write ( *, * ) ' '
  write ( *, * ) '  Parameters:'
  write ( *, * ) ' '
  write ( *, * ) '    DFN = ', dfn
  write ( *, * ) '    DFD = ', dfd
  write ( *, * ) '    PNONC = ', pnonc
  write ( *, * ) ' '
  write ( *, * ) '  X  CDF  1-CDF'
  write ( *, * ) ' '

  do i = 0, 10

    x = real ( i ) / 5.0

    call cumfnc ( x, dfn, dfd, pnonc, cum, ccum )

    write ( *,     '(3g14.6)' ) x, cum, ccum

  end do

  return
end
subroutine test07
!
!*******************************************************************************
!
!! TEST07 tests CUMGAM.
!
  double precision a
  double precision ccum
  double precision cum
  integer i
  double precision x
!
  a = 2.0

  write ( *, * ) ' '
  write ( *, * ) 'TEST07'
  write ( *, * ) '  CUMGAM computes the cumulative density'
  write ( *, * ) '  function for the gamma distribution.'
  write ( *, * ) ' '
  write ( *, * ) '  Parameters:'
  write ( *, * ) ' '
  write ( *, * ) '    A = ', a
  write ( *, * ) ' '
  write ( *, * ) '  X  CDF  1-CDF'
  write ( *, * ) ' '

  do i = 0, 10

    x = real ( i ) / 5.0

    call cumgam ( x, a, cum, ccum )

    write ( *,     '(3g14.6)' ) x, cum, ccum

  end do

  return
end
subroutine test08
!
!*******************************************************************************
!
!! TEST08 tests CUMNBN.
!
  double precision ccum
  double precision cum
  integer i
  double precision ompr
  double precision pr
  double precision s
  double precision xn
!
  pr = 0.75
  ompr = 1.0 - pr
  xn = 10

  write ( *, * ) ' '
  write ( *, * ) 'TEST08'
  write ( *, * ) '  CUMNBN computes the cumulative density'
  write ( *, * ) '  function for the negative binomial '
  write ( *, * ) '  distribution.'
  write ( *, * ) ' '
  write ( *, * ) '  Parameters:'
  write ( *, * ) ' '
  write ( *, * ) '    PR = ', pr
  write ( *, * ) '    XN = ', xn
  write ( *, * ) ' '
  write ( *, * ) '  X  CDF  1-CDF'
  write ( *, * ) ' '

  do i = 0, int ( xn )

    s = real ( i ) 

    call cumnbn ( s, xn, pr, ompr, cum, ccum )

    write ( *,     '(3g14.6)' ) s, cum, ccum

  end do

  return
end
subroutine test09
!
!*******************************************************************************
!
!! TEST09 tests CUMNOR.
!
  integer, parameter :: test_num = 16
!
  double precision ccum
  double precision cum
  integer i
  double precision x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST09'
  write ( *, * ) '  CUMNOR computes the cumulative'
  write ( *, * ) '  density function for the normal distribution.'
  write ( *, * ) ' '
  write ( *, * ) '  X  CDF  1-CDF'
  write ( *, * ) ' '

  do i = 1, test_num

    x = 3.0 * real ( i - 1 ) / real ( test_num - 1 )

    call cumnor ( x, cum, ccum )

    write ( *,     '(3g14.6)' ) x, cum, ccum

  end do

  return
end
subroutine test10
!
!*******************************************************************************
!
!! TEST10 tests CUMPOI.
!
  double precision ccum
  double precision cum
  integer i
  double precision x
  double precision xlam
!
  xlam = 3.0

  write ( *, * ) ' '
  write ( *, * ) 'TEST10'
  write ( *, * ) '  CUMPOI computes the cumulative'
  write ( *, * ) '  density function for the Poisson distribution.'
  write ( *, * ) ' '
  write ( *, * ) '  Parameters:'
  write ( *, * ) ' '
  write ( *, * ) '    XLAM = ', xlam
  write ( *, * ) ' '
  write ( *, * ) '  X  CDF  1-CDF'
  write ( *, * ) ' '

  do i = 0, 10

    x = real ( i )

    call cumpoi ( x, xlam, cum, ccum )

    write ( *,     '(3g14.6)' ) x, cum, ccum

  end do

  return
end
subroutine test11
!
!*******************************************************************************
!
!! TEST11 tests CUMT.
!
  double precision ccum
  double precision cum
  double precision df
  integer i
  double precision x
!
  df = 4.0

  write ( *, * ) ' '
  write ( *, * ) 'TEST11'
  write ( *, * ) '  CUMT computes the cumulative density'
  write ( *, * ) '  function for the T distribution.'
  write ( *, * ) ' '
  write ( *, * ) '  Parameters:'
  write ( *, * ) ' '
  write ( *, * ) '    DF = ', df
  write ( *, * ) ' '
  write ( *, * ) '  X  CDF  1-CDF'
  write ( *, * ) ' '

  do i = 0, 10

    x = real ( i ) / 5.0

    call cumt ( x, df, cum, ccum )

    write ( *,     '(3g14.6)' ) x, cum, ccum

  end do

  return
end
subroutine test12
!
!*******************************************************************************
!
!! TEST12 tests BETA;
!! TEST12 tests GAMMA.
!
  double precision a
  double precision b
  double precision beta
  double precision beta1
  double precision beta2
  double precision gamma
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST12'
  write ( *, * ) '  BETA evaluates the Beta function;'
  write ( *, * ) '  GAMMA evaluates the Gamma function.'

  a = 2.2
  b = 3.7

  beta1 = beta ( a, b )
  beta2 = gamma ( a ) * gamma ( b ) / gamma ( a + b )

  write ( *, * ) ' '
  write ( *, * ) '  Argument A =                   ', a
  write ( *, * ) '  Argument B =                   ', b
  write ( *, * ) '  Beta(A,B) =                    ', beta1
  write ( *, * ) '  (Expected value = 0.0454 )'
  write ( *, * ) ' '
  write ( *, * ) '  Gamma(A)*Gamma(B)/Gamma(A+B) = ', beta2

  return
end
subroutine test13
!
!*******************************************************************************
!
!! TEST13 tests GAMMA_INC;
!! TEST13 tests GAMMA_INC_INV.
!
  integer, parameter :: test_num = 10
!
  double precision a
  integer i
  integer ierror
  integer ind
  double precision p
  double precision q
  double precision x
  double precision x0
  double precision x2
!
  a = 3.0
  ind = 1
  x0 = 0

  write ( *, * ) ' '
  write ( *, * ) 'TEST13'
  write ( *, * ) '  GAMMA_INC evaluates the incomplete gamma ratio;'
  write ( *, * ) '  GAMMA_INC_INV inverts it.'
  write ( *, * ) ' '
  write ( *, * ) '  Parameters:'
  write ( *, * ) ' '
  write ( *, * ) '    A = ', a
  write ( *, * ) ' '
  write ( *, * ) '  X  P  Q  Inverse'
  write ( *, * ) ' '

  do i = 0, test_num

    x = real ( i ) / real ( test_num )
    call gamma_inc ( a, x, p, q, ind )

    call gamma_inc_inv ( a, x2, x0, p, q, ierror )

    write ( *,     '(4g14.6)' ) x, p, q, x2

  end do

  return
end
subroutine test14
!
!*******************************************************************************
!
!! TEST14 tests ERF.
!! TEST14 tests ERFC1.
!
  integer, parameter :: test_num = 20
!
  double precision e1
  double precision e2
  double precision erf
  double precision erfc1
  integer i
  integer ind
  double precision x
!
  write ( *, * ) ' '
  write ( *, * ) 'TEST14'
  write ( *, * ) '  ERF computes the error function;'
  write ( *, * ) '  ERFC1 the complementary error function.'
  write ( *, * ) ' '
  write ( *, * ) ' '
  write ( *, * ) '  X  ERF  ERFC1'
  write ( *, * ) ' '

  ind = 0

  do i = 0, test_num

    x = real ( i ) / real ( 10 )

    e1 = erf ( x )
    e2 = erfc1 ( ind, x )
    write ( *,     '(3g14.6)' ) x, e1, e2

  end do

  return
end

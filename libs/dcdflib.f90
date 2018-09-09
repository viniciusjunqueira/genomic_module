function algdiv ( a, b )
!
!*******************************************************************************
!
!! ALGDIV computes LN(GAMMA(B)/GAMMA(A+B)) when B >= 8.
!
!
!  Discussion:
!
!    In this algorithm, DEL(X) is the function defined by
!    LN(GAMMA(X)) = (X - 0.5)*LN(X) - X + 0.5*LN(2*PI) + DEL(X).
!
!  Parameters:
!
!    Input, double precision A, B, the arguments.
!
!    Output, double precision ALGDIV, the value of LN(GAMMA(B)/GAMMA(A+B)).
!
  double precision, parameter :: c0 =   0.833333333333333d-01
  double precision, parameter :: c1 = - 0.277777777760991d-02
  double precision, parameter :: c2 =   0.793650666825390d-03
  double precision, parameter :: c3 = - 0.595202931351870d-03
  double precision, parameter :: c4 =   0.837308034031215d-03
  double precision, parameter :: c5 = - 0.165322962780713d-02
!
  double precision a
  double precision algdiv
  double precision alnrel
  double precision b
  double precision c
  double precision d
  double precision h
  double precision s11
  double precision s3
  double precision s5
  double precision s7
  double precision s9
  double precision t
  double precision u
  double precision v
  double precision w
  double precision x
  double precision x2
!
  if ( a > b ) then
    h = b / a
    c = 1.0d0 / ( 1.0d0 + h )
    x = h / ( 1.0d0 + h )
    d = a + ( b - 0.5d0 )
  else
    h = a / b
    c = h / ( 1.0d0 + h )
    x = 1.0d0 / ( 1.0d0 + h )
    d = b + ( a - 0.5d0 )
  end if
!
!  Set SN = (1 - X**N)/(1 - X).
!
  x2 = x * x
  s3 = 1.0d0 + ( x + x2 )
  s5 = 1.0d0 + ( x + x2 * s3 )
  s7 = 1.0d0 + ( x + x2 * s5 )
  s9 = 1.0d0 + ( x + x2 * s7 )
  s11 = 1.0d0 + ( x + x2 * s9 )
!
!  Set W = DEL(B) - DEL(A + B).
!
  t = ( 1.0d0 / b )**2
  w = (((( &
          c5 * s11  * t &
        + c4 * s9 ) * t &
        + c3 * s7 ) * t &
        + c2 * s5 ) * t &
        + c1 * s3 ) * t &
        + c0

  w = w * ( c / b )
!
!  Combine the results.
!
  u = d * alnrel ( a / b )
  v = a * ( log ( b ) - 1.0d0 )

  if ( u > v ) then
    algdiv = ( w - v ) - u
  else
    algdiv = ( w - u ) - v
  end if

  return
end
function alnrel ( a )
!
!*******************************************************************************
!
!! ALNREL evaluates the function LN(1 + A).
!
!
!  Reference:
!
!    A R DiDinato and A H Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
  double precision, parameter :: p1 = - 0.129418923021993d+01
  double precision, parameter :: p2 =   0.405303492862024d+00
  double precision, parameter :: p3 = - 0.178874546012214d-01
  double precision, parameter :: q1 = - 0.162752256355323d+01
  double precision, parameter :: q2 =   0.747811014037616d+00
  double precision, parameter :: q3 = - 0.845104217945565d-01
!
  double precision a
  double precision alnrel
  double precision t
  double precision t2
  double precision w
  double precision x
!
  if ( abs ( a ) <= 0.375d0 ) then

    t = a / ( a + 2.0d0 )
    t2 = t * t

    w = ((( p3 * t2 + p2 ) * t2 + p1 ) * t2 + 1.0d0 ) &
      / ((( q3 * t2 + q2 ) * t2 + q1 ) * t2 + 1.0d0 )

    alnrel = 2.0d0 * t * w

  else

    x = 1.0d0 + dble ( a )
    alnrel = log ( x )

  end if

  return
end
function apser ( a, b, x, eps )
!
!*******************************************************************************
!
!! APSER computes the incomplete beta ratio I(SUB(1-X))(B,A).
!
!
!  APSER is used for
!     A .LE. MIN(EPS,EPS*B), B*X .LE. 1, AND X .LE. 0.5. Used when
!     A is very small. Use only if above inequalities are satisfied.
!
  double precision, parameter :: g = 0.577215664901533d0
!
  double precision a
  double precision aj
  double precision apser
  double precision b
  double precision bx
  double precision c
  double precision eps
  double precision j
  double precision psi
  double precision s
  double precision t
  double precision tol
  double precision x
!
  bx = b * x
  t = x - bx

  if ( b * eps <= 2.0d-2 ) then
    c = log ( x ) + psi ( b ) + g + t
  else
    c = log ( bx ) + g + t
  end if

  tol = 5.0d0 * eps * abs ( c )
  j = 1.0d0
  s = 0.0d0

  do

    j = j + 1.0d0
    t = t * ( x - bx / j )
    aj = t / j
    s = s + aj

    if ( abs ( aj ) <= tol ) then
      exit
    end if

  end do

  apser = - a * ( c + s )

  return
end
function beta_asym ( a, b, lambda, eps )
!
!*******************************************************************************
!
!! BETA_ASYM computes an asymptotic expansion for IX(A,B), for large A and B.
!
!
!  Discussion:
!
!    LAMBDA = (A + B)*Y - B  and EPS is the tolerance used.
!    It is assumed that LAMBDA is nonnegative and that
!    A and B are greater than or equal to 15.
!
  double precision, parameter :: e0 = 1.12837916709551d0
  double precision, parameter :: e1 = 0.353553390593274d0
  integer, parameter :: num = 20
!
  double precision a
  double precision a0(num+1)
  double precision b
  double precision b0(num+1)
  double precision bcorr
  double precision beta_asym
  double precision bsum
  double precision c(num+1)
  double precision d(num+1)
  double precision dsum
  double precision eps
  double precision erfc1
  double precision f
  double precision h
  double precision h2
  double precision hn
  integer i
  integer j
  double precision j0
  double precision j1
  double precision lambda
  integer m
  integer mm1
  integer mmj
  integer n
  integer np1
  double precision r
  double precision r0
  double precision r1
  double precision rlog1
  double precision s
  double precision sum
  double precision t
  double precision t0
  double precision t1
  double precision u
  double precision w
  double precision w0
  double precision z
  double precision z0
  double precision z2
  double precision zn
  double precision znm1
!
  beta_asym = 0.0d0

  if ( a < b ) then
    h = a / b
    r0 = 1.0d0 / ( 1.0d0 + h )
    r1 = ( b - a ) / b
    w0 = 1.0d0 / sqrt ( a * ( 1.0d0 + h ))
  else
    h = b / a
    r0 = 1.0d0 / ( 1.0d0 + h )
    r1 = ( b - a ) / a
    w0 = 1.0d0 / sqrt ( b * ( 1.0d0 + h ))
  end if

  f = a * rlog1 ( - lambda / a ) + b * rlog1 ( lambda / b )
  t = exp ( - f )
  if ( t == 0.0d0 ) then
    return
  end if

  z0 = sqrt ( f )
  z = 0.5d0 * ( z0 / e1 )
  z2 = f + f

  a0(1) = ( 2.0d0 / 3.0d0 ) * r1
  c(1) = - 0.5d0 * a0(1)
  d(1) = - c(1)
  j0 = ( 0.5d0 / e0 ) * erfc1 ( 1, z0 )
  j1 = e1
  sum = j0 + d(1) * w0 * j1

  s = 1.0d0
  h2 = h * h
  hn = 1.0d0
  w = w0
  znm1 = z
  zn = z2

  do n = 2, num, 2

    hn = h2 * hn
    a0(n) = 2.0d0 * r0 * ( 1.0d0 + h * hn ) / ( n + 2.0d0 )
    np1 = n + 1
    s = s + hn
    a0(np1) = 2.0d0 * r1 * s / ( n + 3.0d0 )

    do i = n, np1

      r = - 0.5d0 * ( i + 1.0d0 )
      b0(1) = r * a0(1)
      do m = 2, i
        bsum = 0.0d0
        mm1 = m - 1
        do j = 1, mm1
          mmj = m - j
          bsum = bsum + ( j * r - mmj ) * a0(j) * b0(mmj)
        end do
        b0(m) = r * a0(m) + bsum / m
      end do

      c(i) = b0(i) / ( i + 1.0d0 )

      dsum = 0.0d0
      do j = 1, i-1
        dsum = dsum + d(i-j) * c(j)
      end do
      d(i) = - ( dsum + c(i) )

    end do

    j0 = e1 * znm1 + ( n - 1.0d0 ) * j0
    j1 = e1 * zn + n * j1
    znm1 = z2 * znm1
    zn = z2 * zn
    w = w0 * w
    t0 = d(n) * w * j0
    w = w0 * w
    t1 = d(np1) * w * j1
    sum = sum + ( t0 + t1 )

    if (( abs ( t0 ) + abs ( t1 )) <= eps * sum ) then
      u = exp ( - bcorr ( a, b ) )
      beta_asym = e0 * t * u * sum
      return
    end if

  end do

  u = exp ( - bcorr ( a, b ) )
  beta_asym = e0 * t * u * sum

  return
end
function bcorr ( a0, b0 )
!
!*******************************************************************************
!
!! BCORR evaluates DEL(A0) + DEL(B0) - DEL(A0 + B0)
!
!
!     LN(GAMMA(A)) = (A - 0.5)*LN(A) - A + 0.5*LN(2*PI) + DEL(A).
!     It is assumed that A0 .GE. 8 AND B0 .GE. 8.
!
  double precision a
  double precision a0
  double precision b
  double precision b0
  double precision bcorr
  double precision c
  double precision, parameter :: c0 =   0.833333333333333d-01
  double precision, parameter :: c1 = - 0.277777777760991d-02
  double precision, parameter :: c2 =   0.793650666825390d-03
  double precision, parameter :: c3 = - 0.595202931351870d-03
  double precision, parameter :: c4 =   0.837308034031215d-03
  double precision, parameter :: c5 = - 0.165322962780713d-02
  double precision h
  double precision s11
  double precision s3
  double precision s5
  double precision s7
  double precision s9
  double precision t
  double precision w
  double precision x
  double precision x2
!
  a = min ( a0, b0 )
  b = max ( a0, b0 )

  h = a / b
  c = h / ( 1.0d0 + h )
  x = 1.0d0 / ( 1.0d0 + h )
  x2 = x * x
!
!  Set SN = (1 - X**N)/(1 - X)
!
  s3 = 1.0d0 + ( x + x2 )
  s5 = 1.0d0 + ( x + x2 * s3 )
  s7 = 1.0d0 + ( x + x2 * s5 )
  s9 = 1.0d0 + ( x + x2 * s7 )
  s11 = 1.0d0 + ( x + x2 * s9 )
!
!  Set W = DEL(B) - DEL(A + B)
!
  t = ( 1.0d0 / b )**2

  w = (((( &
       c5 * s11  * t &
     + c4 * s9 ) * t &
     + c3 * s7 ) * t &
     + c2 * s5 ) * t &
     + c1 * s3 ) * t &
     + c0

  w = w * ( c / b )
!
!  Compute  DEL(A) + W.
!
  t = ( 1.0d0 / a )**2

  bcorr = ((((( &
         c5   * t &
       + c4 ) * t &
       + c3 ) * t &
       + c2 ) * t &
       + c1 ) * t &
       + c0 ) / a + w

  return
end
function beta ( a, b )
!
!*******************************************************************************
!
!! BETA evaluates the beta function.
!
!
!  Modified:
!
!    03 December 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, double precision A, B, the arguments of the beta function.
!
!    Output, double precision BETA, the value of the beta function.
!
  double precision a
  double precision b
  double precision beta
  double precision beta_log
!
  beta = exp ( beta_log ( a, b ) )

  return
end
function beta_log ( a0, b0 )
!
!*******************************************************************************
!
!! BETA_LOG evaluates the logarithm of the beta function.
!
!
!  Reference:
!
!    A R DiDinato and A H Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Trans. Math.  Software, 
!    Volume 18, 1993, pages 360-373.
!
  double precision, parameter :: e = 0.918938533204673d0
!
  double precision a
  double precision a0
  double precision algdiv
  double precision alnrel
  double precision b
  double precision b0
  double precision bcorr
  double precision beta_log
  double precision c
  double precision gamma_log
  double precision gsumln
  double precision h
  integer i
  integer n
  double precision u
  double precision v
  double precision w
  double precision z
!
  a = min ( a0, b0 )
  b = max ( a0, b0 )
  if (a>=8.0d0) go to 100
  if (a>=1.0d0) go to 20
!
!  Procedure when A .LT. 1
!
  if ( b < 8.0d0 ) then
    beta_log = gamma_log ( a ) + ( gamma_log ( b ) - gamma_log ( a + b ) )
  else
    beta_log = gamma_log ( a ) + algdiv ( a, b )
  end if

  return
!
!  Procedure when 1 .LE. A .LT. 8
!
   20 continue

  if (a>2.0d0) go to 40

  if ( b <= 2.0d0 ) then
    beta_log = gamma_log ( a ) + gamma_log ( b ) - gsumln ( a, b )
    return
  end if

   30 continue

  w = 0.0d0
  if ( b < 8.0d0 ) go to 60
  beta_log = gamma_log ( a ) + algdiv ( a, b )
  return
!
!  Reduction of A when B .LE. 1000
!
   40 continue

  if ( b > 1000.0d0 ) go to 80
  n = a - 1.0d0
  w = 1.0d0
  do i = 1, n
    a = a - 1.0d0
    h = a / b
    w = w * ( h / ( 1.0d0 + h ))
  end do
  w = log ( w )

  if ( b >= 8.0d0 ) then
    beta_log = w + gamma_log ( a ) + algdiv ( a, b )
    return
  end if
!
!  Reduction of B when B < 8.
!
60 continue

  n = b - 1.0d0
  z = 1.0d0
  do i = 1, n
    b = b - 1.0d0
    z = z * ( b / ( a + b ))
  end do

  beta_log = w + log ( z ) + ( gamma_log ( a ) + ( gamma_log ( b ) &
    - gsumln ( a, b ) ) )
  return
!
!  Reduction of A when B > 1000
!
80 continue

  n = a - 1.0d0
  w = 1.0d0
  do i = 1, n
    a = a - 1.0d0
    w = w * ( a / ( 1.0d0 + a / b ))
  end do

  beta_log = ( log ( w ) - n * log ( b )) + ( gamma_log ( a ) + algdiv ( a, b ))

  return
!
!  Procedure when A .GE. 8
!
100 continue

  w = bcorr ( a, b )
  h = a / b
  c = h / ( 1.0d0 + h )
  u = - ( a - 0.5d0 ) * log ( c )
  v = b * alnrel ( h )

  if ( u > v ) then
    beta_log = ((( -0.5d0 * log ( b ) + e ) + w ) - v ) - u
  else
    beta_log = ((( -0.5d0 * log ( b ) + e ) + w ) - u ) - v
  end if

  return
end
subroutine beta_ratio ( a, b, x, y, w, w1, ierr )
!
!*******************************************************************************
!
!! BETA_RATIO evaluates the incomplete beta function IX(A,B).
!
!
!  Author:
!
!    Alfred H Morris, Jr,
!    Naval Surface Weapons Center,
!    Dahlgren, Virginia.
!
!     It is assumed that A and B are nonnegative, and that X .LE. 1
!     and Y = 1 - X.  BETA_RATIO assigns W and W1 the values
!
!                      W  = IX(A,B)
!                      W1 = 1 - IX(A,B)
!
!     IERR is a variable that reports the status of the results.
!     If no input errors are detected then IERR is set to 0 and
!     W and W1 are computed. Otherwise, if an error is detected,
!     then W and W1 are assigned the value 0 and IERR is set to
!     one of the following values:
!
!        IERR = 1  IF A OR B IS NEGATIVE
!        IERR = 2  IF A = B = 0
!        IERR = 3  IF X .LT. 0 OR X .GT. 1
!        IERR = 4  IF Y .LT. 0 OR Y .GT. 1
!        IERR = 5  IF X + Y .NE. 1
!        IERR = 6  IF X = A = 0
!        IERR = 7  IF Y = B = 0
!
  double precision a
  double precision a0
  double precision apser
  double precision b
  double precision b0
  double precision beta_asym
  double precision beta_frac
  double precision beta_pser
  double precision beta_up
  double precision dpmpar
  double precision eps
  double precision fpser
  integer ierr
  integer ierr1
  integer ind
  double precision lambda
  integer n
  double precision t
  double precision w
  double precision w1
  double precision x
  double precision x0
  double precision y
  double precision y0
  double precision z
!
!     EPS IS A MACHINE DEPendENT CONSTANT. EPS IS THE SMALLEST
!     FLOATING POINT NUMBER FOR WHICH 1.0 + EPS .GT. 1.0
!
  eps = epsilon ( 1.0d0 )
  w = 0.0d0
  w1 = 0.0d0
  if (a<0.0d0 .or. b<0.0d0) go to 270
  if (a==0.0d0 .and. b==0.0d0) go to 280
  if (x<0.0d0 .or. x>1.0d0) go to 290
  if (y<0.0d0 .or. y>1.0d0) go to 300
  z = ((x+y)-0.5d0) - 0.5d0
  if ( abs ( z ) > 3.0d0 * eps ) go to 310

  ierr = 0
  if (x==0.0d0) go to 210
  if (y==0.0d0) go to 230
  if (a==0.0d0) go to 240
  if (b==0.0d0) go to 220

  eps = max ( eps, 1.0d-15 )
  if ( max ( a, b ) < 1.0d-3*eps ) go to 260

  ind = 0
  a0 = a
  b0 = b
  x0 = x
  y0 = y

  if ( min ( a0, b0 ) > 1.0d0 ) then
    go to 40
  end if
!
!  Procedure for A0 .LE. 1 OR B0 .LE. 1
!
  if ( x >0.5d0) then
    ind = 1
    a0 = b
    b0 = a
    x0 = y
    y0 = x
  end if

  if ( b0 < min ( eps, eps*a0 ) ) go to 90
  if ( a0 < min (eps,eps*b0) .and. b0*x0<=1.0d0) go to 100
  if ( max (a0,b0)>1.0d0) go to 20
  if ( a0 >= min (0.2d0,b0)) go to 110
  if ( x0**a0 <= 0.9d0 ) go to 110
  if ( x0 >= 0.3d0 ) go to 120
  n = 20
  go to 140

20 continue

  if ( b0 <= 1.0d0 ) go to 110
  if ( x0 >= 0.3d0 ) go to 120
  if ( x0 >= 0.1d0 ) go to 30
  if ( (x0*b0)**a0 <= 0.7d0 ) go to 110

30 continue

  if (b0>15.0d0) go to 150
  n = 20
  go to 140
!
!  PROCEDURE FOR A0 .GT. 1 AND B0 .GT. 1
!
40 continue

  if ( a <= b ) then
    lambda = a - ( a + b ) * x
  else
    lambda = ( a + b ) * y - b
  end if

  if ( lambda < 0.0d0 ) then
    ind = 1
    a0 = b
    b0 = a
    x0 = y
    y0 = x
    lambda = abs ( lambda )
  end if

70 continue

  if (b0<40.0d0 .and. b0*x0<=0.7d0) go to 110
  if (b0<40.0d0) go to 160
  if (a0>b0) go to 80
  if (a0<=100.0d0) go to 130
  if (lambda>0.03d0*a0) go to 130
  go to 200

80 continue

  if (b0<=100.0d0) go to 130
  if (lambda>0.03d0*b0) go to 130
  go to 200
!
!  Evaluation of the appropriate algorithm
!
90 continue

  w = fpser ( a0, b0, x0, eps )
  w1 = 0.5d0 + ( 0.5d0 - w )
  go to 250

100 continue

  w1 = apser ( a0, b0, x0, eps )
  w = 0.5d0 + ( 0.5d0 - w1 )
  go to 250

110 continue

  w = beta_pser ( a0, b0, x0, eps )
  w1 = 0.5d0 + ( 0.5d0 - w )
  go to 250

120 continue

  w1 = beta_pser ( b0, a0, y0, eps )
  w = 0.5d0 + ( 0.5d0 - w1 )
  go to 250

130 continue

  w = beta_frac ( a0, b0, x0, y0, lambda, 15.0d0*eps )
  w1 = 0.5d0 + ( 0.5d0 - w )
  go to 250

140 continue

  w1 = beta_up ( b0, a0, y0, x0, n, eps )
  b0 = b0 + n

150 continue

  call beta_grat ( b0, a0, y0, x0, w1, 15.0d0 * eps, ierr1 )
  w = 0.5d0 + ( 0.5d0 - w1 )
  go to 250

160 continue

  n = b0
  b0 = b0 - n

  if ( b0 == 0.0d0 ) then
    n = n - 1
    b0 = 1.0d0
  end if

170 continue

  w = beta_up ( b0, a0, y0, x0, n, eps )
  if ( x0 > 0.7d0 ) go to 180
  w = w + beta_pser ( a0, b0, x0, eps )
  w1 = 0.5d0 + ( 0.5d0 - w )
  go to 250

180 continue

  if ( a0 <= 15.0d0 ) then
    n = 20
    w = w + beta_up ( a0, b0, x0, y0, n, eps )
    a0 = a0 + n
  end if

190 continue

  call beta_grat ( a0, b0, x0, y0, w, 15.0d0 * eps, ierr1 )
  w1 = 0.5d0 + ( 0.5d0 - w )
  go to 250

200 continue

  w = beta_asym ( a0, b0, lambda, 100.0d0 * eps )
  w1 = 0.5d0 + ( 0.5d0 - w )
  go to 250
!
!  Termination of the procedure.
!
210 continue

  if (a==0.0d0) go to 320
220 continue
  w = 0.0d0
  w1 = 1.0d0
  return

230 continue

  if (b==0.0d0) go to 330
240 continue
  w = 1.0d0
  w1 = 0.0d0
  return

250 continue

  if ( ind /= 0 ) then
    t = w
    w = w1
    w1 = t
  end if

  return
!
!  Procedure for A AND B .LT. 1.E-3*EPS
!
260 continue
  w = b / ( a + b )
  w1 = a / ( a + b )
  return
!
!  Error return
!
270 continue
  ierr = 1
  return

  280 ierr = 2
  return

  290 ierr = 3
  return

  300 ierr = 4
  return

  310 ierr = 5
  return

  320 ierr = 6
  return

  330 ierr = 7

  return
end
function beta_frac ( a, b, x, y, lambda, eps )
!
!*******************************************************************************
!
!! BETA_FRAC evaluates a continued fraction expansion for IX(A,B) when A,B > 1.
!
!
!  IT IS ASSUMED THAT  LAMBDA = (A + B)*Y - B.
!
  double precision a
  double precision alpha
  double precision an
  double precision anp1
  double precision b
  double precision beta
  double precision beta_frac
  double precision beta_rcomp
  double precision bn
  double precision bnp1
  double precision c
  double precision c0
  double precision c1
  double precision e
  double precision eps
  double precision lambda
  double precision n
  double precision p
  double precision r
  double precision r0
  double precision s
  double precision t
  double precision w
  double precision x
  double precision y
  double precision yp1
!
  beta_frac = beta_rcomp ( a, b, x, y )

  if ( beta_frac == 0.0d0 ) then
    return
  end if

  c = 1.0d0 + lambda
  c0 = b / a
  c1 = 1.0d0 + 1.0d0 / a
  yp1 = y + 1.0d0

  n = 0.0d0
  p = 1.0d0
  s = a + 1.0d0
  an = 0.0d0
  bn = 1.0d0
  anp1 = 1.0d0
  bnp1 = c / c1
  r = c1 / c
!
!  Continued fraction calculation.
!
   10 continue

  n = n + 1.0d0
  t = n / a
  w = n * ( b - n ) * x
  e = a / s
  alpha = ( p * ( p + c0 ) * e * e ) * ( w * x )
  e = ( 1.0d0 + t ) / ( c1 + t + t )
  beta = n + w / s + e * ( c + n * yp1 )
  p = 1.0d0 + t
  s = s + 2.0d0
!
!  Update AN, BN, ANP1, and BNP1.
!
  t = alpha * an + beta * anp1
  an = anp1
  anp1 = t
  t = alpha * bn + beta * bnp1
  bn = bnp1
  bnp1 = t

  r0 = r
  r = anp1 / bnp1

  if ( abs ( r - r0 ) <= eps * r ) then
    beta_frac = beta_frac * r
    return
  end if
!
!  Rescale AN, BN, ANP1, and BNP1.
!
  an = an / bnp1
  bn = bn / bnp1
  anp1 = r
  bnp1 = 1.0d0
  go to 10

end
subroutine beta_grat ( a, b, x, y, w, eps, ierr )
!
!*******************************************************************************
!
!! BETA_GRAT evaluates an asymptotic expansion for IX(A,B).
!
!
!  This routine is used when A is larger than B.
!
!
!     The result of the expansion is added to W. IT IS ASSUMED
!     THAT A .GE. 15 AND B .LE. 1.  EPS IS THE TOLERANCE USED.
!     IERR IS A VARIABLE THAT REPORTS THE STATUS OF THE RESULTS.
!
  double precision a
  double precision algdiv
  double precision alnrel
  double precision b
  double precision bm1
  double precision bp2n
  double precision c(30)
  double precision cn
  double precision coef
  double precision d(30)
  double precision dj
  double precision eps
  double precision gam1
  integer i
  integer ierr
  double precision j
  double precision l
  double precision lnx
  integer n
  double precision n2
  double precision nu
  double precision p
  double precision q
  double precision r
  double precision s
  double precision sum
  double precision t
  double precision t2
  double precision u
  double precision v
  double precision w
  double precision x
  double precision y
  double precision z
!
  bm1 = ( b - 0.5d0 ) - 0.5d0
  nu = a + 0.5d0 * bm1

  if ( y <= 0.375d0 ) then
    lnx = alnrel ( - y )
  else
    lnx = log ( x )
  end if

  z = - nu * lnx

  if ( b * z == 0.0d0 ) then
    ierr = 1
    return
  end if
!
!  Computation of the expansion
!  SET R = EXP(-Z)*Z**B/GAMMA(B)
!
  r = b * ( 1.0d0 + gam1 ( b ) ) * exp ( b * log ( z ))
  r = r * exp ( a * lnx ) * exp ( 0.5d0 * bm1 * lnx )
  u = algdiv ( b, a ) + b * log ( nu )
  u = r * exp ( - u )

  if ( u == 0.0d0 ) then
    ierr = 1
    return
  end if

  call gamma_rat1 ( b, z, r, p, q, eps )

  v = 0.25d0 * ( 1.0d0 / nu )**2
  t2 = 0.25d0 * lnx * lnx
  l = w / u
  j = q / r
  sum = j
  t = 1.0d0
  cn = 1.0d0
  n2 = 0.0d0

  do n = 1, 30

    bp2n = b + n2
    j = ( bp2n * ( bp2n + 1.0d0 ) * j + ( z + bp2n + 1.0d0 ) * t ) * v
    n2 = n2 + 2.0d0
    t = t * t2
    cn = cn / ( n2 * ( n2 + 1.0d0 ))
    c(n) = cn
    s = 0.0d0

    coef = b - n
    do i = 1, n-1
      s = s + coef * c(i) * d(n-i)
      coef = coef + b
    end do

    d(n) = bm1 * cn + s / n
    dj = d(n) * j
    sum = sum + dj

    if ( sum <= 0.0d0 ) then
      ierr = 1
      return
    end if

    if ( abs ( dj ) <= eps * ( sum + l )) then
      ierr = 0
      w = w + u * sum
      return
    end if

  end do

  ierr = 0
  w = w + u * sum

  return
end
function beta_pser ( a, b, x, eps )
!
!*******************************************************************************
!
!! BETA_PSER uses a power series expansion to evaluate IX(A,B).
!
!
!  BETA_PSER is used WHEN B .LE. 1
!     OR B*X .LE. 0.7.  EPS IS THE TOLERANCE USED.
!
  double precision a
  double precision a0
  double precision algdiv
  double precision apb
  double precision b
  double precision b0
  double precision beta_log
  double precision beta_pser
  double precision c
  double precision eps
  double precision gam1
  double precision gamma_ln1
  integer i
  integer m
  double precision n
  double precision sum
  double precision t
  double precision tol
  double precision u
  double precision w
  double precision x
  double precision z
!
  beta_pser = 0.0d0
  if ( x == 0.0d0 ) then
    return
  end if
!
!  Compute the factor X**A/(A*BETA(A,B))
!
  a0 = min ( a, b )

  if ( a0 >= 1.0d0 ) then

    z = a * log ( x ) - beta_log ( a, b )
    beta_pser = exp ( z ) / a
  
  else

    b0 = max ( a, b )

    if ( b0 <= 1.0 ) then

      beta_pser = x**a
      if ( beta_pser == 0.0d0 ) then
        return
      end if

      apb = a + b

      if ( apb <= 1.0d0 ) then
        z = 1.0d0 + gam1 ( apb )
      else
        u = a + b - 1.d0
        z = ( 1.0d0 + gam1 ( u ) ) / apb
      end if

      c = ( 1.0d0 + gam1 ( a ) ) * ( 1.0d0 + gam1 ( b ) ) / z
      beta_pser = beta_pser * c * ( b / apb )

    else if ( b0 < 8.0 ) then

      u = gamma_ln1 ( a0 )
      m = b0 - 1.0d0

      c = 1.0d0
      do i = 1, m
        b0 = b0 - 1.0d0
        c = c * ( b0 / ( a0 + b0 ))
      end do

      u = log ( c ) + u
      z = a * log ( x ) - u
      b0 = b0 - 1.0d0
      apb = a0 + b0

      if ( apb <= 1.0d0 ) then
        t = 1.0d0 + gam1 ( apb )
      else
        u = a0 + b0 - 1.d0
        t = ( 1.0d0 + gam1 ( u ) ) / apb
      end if

      beta_pser = exp ( z ) * ( a0 / a )* ( 1.0d0 + gam1 ( b0 )) / t

    else if ( b0 >= 8.0 ) then

      u = gamma_ln1 ( a0 ) + algdiv ( a0, b0 )
      z = a * log ( x ) - u
      beta_pser = ( a0 / a ) * exp ( z )

    end if

  end if

  if ( beta_pser == 0.0d0 .or. a <= 0.1d0 * eps ) then
    return
  end if
!
!  Compute the series
!
  sum = 0.0d0
  n = 0.0d0
  c = 1.0d0
  tol = eps / a

10 continue

  n = n + 1.0d0
  c = c * ( 0.5d0 + ( 0.5d0 - b / n )) * x
  w = c / ( a + n )
  sum = sum + w
  if ( abs ( w ) > tol ) then
    go to 10
  end if

  beta_pser = beta_pser * ( 1.0d0 + a * sum )

  return
end
function beta_rcomp1 ( mu, a, b, x, y )
!
!*******************************************************************************
!
!! BETA_RCOMP1 evaluates EXP(MU) * (X**A*Y**B/BETA(A,B)).
!
  double precision, parameter :: const = 0.398942280401433d0
!
  double precision a
  double precision a0
  double precision apb
  double precision b
  double precision b0
  double precision beta_rcomp1
  double precision c
  double precision e
  double precision esum
  double precision gam1
  double precision gamma_ln1
  double precision h
  integer i
  double precision lambda
  double precision lnx
  double precision lny
  integer mu
  integer n
  double precision rlog1
  double precision t
  double precision u
  double precision v
  double precision x
  double precision x0
  double precision y
  double precision y0
  double precision z

  double precision algdiv,alnrel,bcorr,beta_log
!
  a0 = min ( a, b )

  if ( a0 >= 8.0d0 ) go to 130

  if ( x <= 0.375d0 ) then
    lnx = log ( x )
    lny = alnrel ( - x )
  else if ( y <= 0.375d0 ) then
    lnx = alnrel ( - y )
    lny = log ( y )
  else
    lnx = log ( x )
    lny = log ( y )
  end if

  z = a * lnx + b * lny

  if ( a0 >= 1.0d0 ) then
    z = z - beta_log(a,b)
    beta_rcomp1 = esum ( mu, z )
    return
  end if
!
!  Procedure for A .LT. 1 OR B .LT. 1
!
   40 continue

  b0 = max ( a, b )
  if (b0>=8.0d0) go to 120
  if (b0>1.0d0) go to 70
!
!  Algorithm for B0 .LE. 1
!
  beta_rcomp1 = esum ( mu, z )
  if ( beta_rcomp1 == 0.0d0 ) then
    return
  end if

  apb = a + b

  if ( apb <= 1.0d0 ) then
    z = 1.0d0 + gam1(apb)
  else
    u = dble(a) + dble(b) - 1.d0
    z = ( 1.0d0 + gam1 ( u )) / apb
  end if

  c = ( 1.0d0 + gam1 ( a ) )* ( 1.0d0 + gam1 ( b ) ) / z
  beta_rcomp1 = beta_rcomp1* ( a0 * c )/ ( 1.0d0 + a0 / b0 )
  return
!
!  Algorithm for 1 .LT. B0 .LT. 8
!
   70 continue

  u = gamma_ln1 ( a0 )
  n = b0 - 1.0d0

  c = 1.0d0
  do i = 1, n
    b0 = b0 - 1.0d0
    c = c * ( b0 / ( a0 + b0 ))
  end do
  u = log ( c ) + u

  z = z - u
  b0 = b0 - 1.0d0
  apb = a0 + b0

  if ( apb <= 1.0d0 ) then
    t = 1.0d0 + gam1 ( apb )
  else
    u = a0 + b0 - 1.d0
    t = ( 1.0d0 + gam1 ( u ) ) / apb
  end if

  beta_rcomp1 = a0 * esum ( mu, z ) * (1.0d0+gam1(b0))/t
  return
!
!  Algorithm for B0 .GE. 8
!
  120 continue

  u = gamma_ln1 ( a0 ) + algdiv ( a0, b0 )
  beta_rcomp1 = a0 * esum ( mu, z-u )
  return
!
!  Procedure for A .GE. 8 and B .GE. 8
!
  130 continue

  if ( a <= b ) then
    h = a / b
    x0 = h / ( 1.0d0 + h )
    y0 = 1.0d0 / ( 1.0d0 + h )
    lambda = a - ( a + b ) * x
  else
    h = b / a
    x0 = 1.0d0 / ( 1.0d0 + h )
    y0 = h / ( 1.0d0 + h )
    lambda = ( a + b ) * y - b
  end if

  e = - lambda / a

  if ( abs ( e ) <= 0.6d0 ) then
    u = rlog1 ( e )
  else
    u = e - log ( x / x0 )
  end if

  e = lambda / b

  if ( abs ( e ) <= 0.6d0 ) then
    v = rlog1 ( e )
  else
    v = e - log ( y / y0 )
  end if

  z = esum ( mu, - ( a * u + b * v ))
  beta_rcomp1 = const * sqrt ( b * x0 ) * z * exp ( - bcorr ( a, b ))

  return
end
function beta_rcomp ( a, b, x, y )
!
!*******************************************************************************
!
!! BETA_RCOMP evaluates X**A * Y**B / BETA(A,B).
!
  double precision, parameter :: const = 0.398942280401433d0
!
  double precision a
  double precision a0
  double precision algdiv
  double precision alnrel
  double precision apb
  double precision b
  double precision b0
  double precision bcorr
  double precision beta_log
  double precision beta_rcomp
  double precision c
  double precision e
  double precision gam1
  double precision gamma_ln1
  double precision h
  integer i
  double precision lambda
  double precision lnx
  double precision lny
  integer n
  double precision rlog1
  double precision t
  double precision u
  double precision v
  double precision x
  double precision x0
  double precision y
  double precision y0
  double precision z
!
  beta_rcomp = 0.0d0
  if ( x == 0.0d0 .or. y == 0.0d0 ) then
    return
  end if

  a0 = min ( a, b )

  if ( a0 < 8.0d0 ) then

    if ( x <= 0.375d0 ) then
      lnx = log ( x )
      lny = alnrel ( - x )      
    else if ( y <= 0.375d0 ) then
      lnx = alnrel ( - y )
      lny = log ( y )
    else
      lnx = log ( x )
      lny = log ( y )
    end if

    z = a * lnx + b * lny

    if ( a0 >= 1.0d0 ) then
      z = z - beta_log ( a, b )
      beta_rcomp = exp ( z )
      return
    end if
!
!  Procedure for A .LT. 1 or B .LT. 1
!
    b0 = max ( a, b )

    if ( b0 <= 1.0 ) then

      beta_rcomp = exp ( z )
      if ( beta_rcomp == 0.0d0 ) then
        return
      end if

      apb = a + b

      if ( apb <= 1.0d0 ) then
        z = 1.0d0 + gam1 ( apb )
      else
        u = a + b - 1.d0
        z = ( 1.0d0 + gam1 ( u ) ) / apb
      end if

      c = ( 1.0d0 + gam1 ( a ) )* ( 1.0d0 + gam1 ( b ) ) / z
      beta_rcomp = beta_rcomp * ( a0 * c ) / ( 1.0d0 + a0 / b0 )

    else if ( b0 < 8.0 ) then

      u = gamma_ln1 ( a0 )
      n = b0 - 1.0d0

      c = 1.0d0
      do i = 1, n
        b0 = b0 - 1.0d0
        c = c * ( b0 / ( a0 + b0 ))
      end do
      u = log ( c ) + u

      z = z - u
      b0 = b0 - 1.0d0
      apb = a0 + b0

      if ( apb <= 1.0d0 ) then
        t = 1.0d0 + gam1 ( apb )
      else
        u = a0 + b0 - 1.d0
        t = ( 1.0d0 + gam1 ( u ) ) / apb
      end if

      beta_rcomp = a0 * exp ( z ) * ( 1.0d0 + gam1 ( b0 ) ) / t

    else if ( b0 >= 8.0 ) then

      u = gamma_ln1 ( a0 ) + algdiv ( a0, b0 )
      beta_rcomp = a0 * exp ( z - u ) 

    end if

  else

    if ( a <= b ) then
      h = a / b
      x0 = h / ( 1.0d0 + h )
      y0 = 1.0d0 / ( 1.0d0 + h )
      lambda = a - ( a + b ) * x
    else
      h = b / a
      x0 = 1.0d0 / ( 1.0d0 + h )
      y0 = h / ( 1.0d0 + h )
      lambda = ( a + b ) * y - b
    end if

    e = - lambda / a

    if ( abs ( e ) <= 0.6d0 ) then
      u = rlog1 ( e )
    else
      u = e - log ( x / x0 )
    end if

    e = lambda / b
 
    if ( abs ( e ) <= 0.6d0 ) then
      v = rlog1 ( e )
    else
      v = e - log ( y / y0 )
    end if

    z = exp ( - ( a * u + b * v ))
    beta_rcomp = const * sqrt ( b * x0 ) * z * exp ( - bcorr ( a, b ))

  end if

  return
end
function beta_up ( a, b, x, y, n, eps )
!
!*******************************************************************************
!
!! BETA_UP evaluates IX(A,B) - IX(A+N,B) where N is a positive integer.
!
!
!  EPS is the tolerance used.
!
  double precision a
  double precision ap1
  double precision apb
  double precision b
  double precision beta_rcomp1
  double precision beta_up
  double precision d
  double precision eps
  double precision exparg
  integer i
  integer k
  double precision l
  integer mu
  integer n
  double precision r
  double precision t
  double precision w
  double precision x
  double precision y
!
!  Obtain the scaling factor EXP(-MU) AND
!  EXP(MU)*(X**A*Y**B/BETA(A,B))/A
!
  apb = a + b
  ap1 = a + 1.0d0
  mu = 0
  d = 1.0d0

  if ( n /= 1 ) then

    if ( a >= 1.0 ) then

      if ( apb >= 1.1d0 * ap1 ) then
        mu = abs ( exparg(1) )
        k = exparg(0)
        if ( k < mu ) then
          mu = k
        end if
        t = mu
        d = exp ( - t )
      end if

    end if

  end if

  beta_up = beta_rcomp1 ( mu, a, b, x, y ) / a

  if ( n == 1 .or. beta_up == 0.0d0 ) then
    return
  end if

  w = d
!
!  Let K be the index of the maximum term.
!
  k = 0
  if ( b <= 1.0d0 ) go to 50

  if ( y <= 1.0d-4 ) then
    k = n - 1
  else
    r = ( b - 1.0d0 ) * x / y - a
    if ( r < 1.0d0 ) go to 50
    k = n - 1
    t = n - 1
    if ( r < t ) then
      k = r
    end if
  end if
!
!  Add the increasing terms of the series.
!
  do i = 1, k
    l = i - 1
    d = (( apb + l ) / ( ap1 + l ) ) * x * d
    w = w + d
  end do
!
!  Add the remaining terms of the series.
!
   50 continue

  do i = k+1, n-1
    l = i - 1
    d = (( apb + l ) / ( ap1 + l )) * x * d
    w = w + d
    if ( d <= eps * w ) then
      beta_up = beta_up * w
      return
    end if
  end do
!
!  Terminate the procedure.
!
  beta_up = beta_up * w

  return
end
subroutine cdfbet ( which, p, q, x, y, a, b, status, bound )
!
!*******************************************************************************
!
!! CDFBET evaluates the Cumulative Distribution Function of the Beta Distribution.
!
!
!  Reference:
!
!    A R DiDinato and A H Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Trans. Math.  Software, 
!    Volume 18, 1993, pages 360-373.
!
!  Discussion:
!
!    Calculates any one parameter of the beta distribution given the others.
!
!  Method:
!
!    Cumulative distribution function  (P)  is calculated directly by
!    code associated with the following reference.
!
!    Computation of other parameters involve a seach for a value that
!    produces  the desired  value  of P.   The search relies  on  the
!    monotonicity of P with the other parameter.
!
!   The beta density is proportional to t^(A-1) * (1-t)^(B-1)
!
!  Parameters:
!
!     WHICH --> Integer indicating which of the next four argument
!               values is to be calculated from the others.
!               Legal range: 1..4
!               1 : Calculate P and Q from X,Y,A and B
!               2 : Calculate X and Y from P,Q,A and B
!               3 : Calculate A from P,Q,X,Y and B
!               4 : Calculate B from P,Q,X,Y and A
!
!                    integer WHICH
!
!     P <--> The integral from 0 to X of the chi-square
!            distribution.
!            Input range: [0, 1].
!                    double precision P
!
!     Q <--> 1-P.
!            Input range: [0, 1].
!            P + Q = 1.0.
!                    double precision Q
!
!     X <--> Upper limit of integration of beta density.
!            Input range: [0,1].
!            Search range: [0,1]
!                    double precision X
!
!     Y <--> 1-X.
!            Input range: [0,1].
!            Search range: [0,1]
!            X + Y = 1.0.
!                    double precision Y
!
!     A <--> The first parameter of the beta density.
!            Input range: (0, +infinity).
!            Search range: [1D-300,1D300]
!                    double precision A
!
!     B <--> The second parameter of the beta density.
!            Input range: (0, +infinity).
!            Search range: [1D-300,1D300]
!                    double precision B
!
!     STATUS <-- 0 if calculation completed correctly
!               -I if input parameter number I is out of range
!                1 if answer appears to be lower than lowest
!                  search bound
!                2 if answer appears to be higher than greatest
!                  search bound
!                3 if P + Q /= 1
!                4 if X + Y /= 1
!                    integer STATUS
!
!     BOUND <-- Undefined if STATUS is 0
!
!               Bound exceeded by parameter number I if STATUS is negative.
!
!               Lower search bound if STATUS is 1.
!
!               Upper search bound if STATUS is 2.
!
  double precision, parameter :: tol = 1.0d-8
  double precision, parameter :: atol = 1.0d-50
  double precision, parameter :: inf = 1.0d300
!
  double precision a
  double precision b
  double precision bound
  double precision fx
  double precision p
  double precision q
  integer status
  integer which
  double precision x
  double precision xhi
  double precision xlo
  double precision y

  double precision cum,ccum,xy,pq
  logical qhi,qleft,qporq
  double precision dpmpar
!
  if ( which < 1 ) then
    bound = 1.0d0
    status = -1
    return
  end if

  if ( which > 4 ) then
    bound = 4.0
    status = -1
    return
  end if

  if ( which /= 1 ) then

    if ( p < 0.0d0 ) then
      bound = 0.0d0
      status = - 2
      return
    else if ( p > 1.0 ) then
      bound = 1.0d0
      status = -2
      return
    end if

    if ( q < 0.0d0 ) then
      bound = 0.0d0
      status = - 3
      return
    else
      bound = 1.0d0
      status = -3
      return
    end if

  end if
!
!     X
!
  if ( which /= 2 ) then

    if (.not. ((x<0.0d0).or. (x>1.0d0))) go to 140
    if (.not. (x<0.0d0)) go to 120
    bound = 0.0d0
    go to 130

  120   bound = 1.0d0
  130   status = -4
    return

  140   continue

  end if
!
!     Y
!
  if ( which /= 2 ) then

    if ( y<0.0d0 ) then
      bound = 0.0d0
      status = -5
      return
    else if ( y > 1.0 ) then
      bound = 1.0d0
      status = -5
      return
    end if

  end if

  if ( which /= 3 ) then

    if ( a <= 0.0d0 ) then
      bound = 0.0d0
      status = -6
      return
    end if

  end if

  if ( which /= 4 ) then

    if ( b<=0.0d0 ) then
      bound = 0.0d0
      status = -7
      return
    end if

  end if
!
!     P + Q
!
  if ( which /= 1 ) then

    pq = p + q

    if (.not. ( abs ( ((pq)-0.5d0)-0.5d0 ) > (3.0d0* epsilon ( 1.0d0 ) ))) then
      go to 260
    end if

    if (.not. (pq<0.0d0)) go to 240
    bound = 0.0d0
    go to 250

  240   bound = 1.0d0
  250   status = 3
    return

  260   continue

  end if
!
!     X + Y
!
  if ( which /= 2 ) then
    xy = x + y
    if (.not. ( abs ( ((xy)-0.5d0)-0.5d0)> (3.0d0*epsilon ( 1.0d0 ) ))) then
      go to 300
    end if

    if (.not. (xy<0.0d0)) go to 280
    bound = 0.0d0
    go to 290

  280   bound = 1.0d0
  290   status = 4
    return

  300   continue

  end if

  if ( which /= 1 ) then
    qporq = p <= q
  end if

  if ( which == 1 ) then

    call cumbet ( x, y, a, b, p, q )
    status = 0

  else if ( which == 2 ) then

    call dstzr(0.0d0,1.0d0,atol,tol)

    if ( qporq ) then

      status = 0
      fx = 0.0
      call dzror(status,x,fx,xlo,xhi,qleft,qhi)
      y = 1.0 - x

  320     continue

      if ( status == 1 ) then
        call cumbet(x,y,a,b,cum,ccum)
        fx = cum - p
        call dzror(status,x,fx,xlo,xhi,qleft,qhi)
        y = 1.0 - x
        go to 320
      end if

    else

      status = 0
      fx = 0.0
      call dzror(status,y,fx,xlo,xhi,qleft,qhi)
      x = 1.0 - y

  350     continue

      if ( status == 1 ) then

        call cumbet(x,y,a,b,cum,ccum)
        fx = ccum - q
        call dzror(status,y,fx,xlo,xhi,qleft,qhi)
        x = 1.0 - y
        go to 350

      end if

    end if

    if ( status == -1 ) then

      if ( qleft ) then
        status = 1
        bound = 0.0d0
      else
        status = 2
        bound = 1.0d0
      end if

  end  if

  else if ( which == 3 ) then

    a = 5.0d0
    call dstinv(0.0d0,inf,0.5d0,0.5d0,5.0d0,atol,tol)
    status = 0
    call dinvr ( status, a, fx, qleft, qhi )

  410   continue

    if ( status == 1 ) then

      call cumbet(x,y,a,b,cum,ccum)

      if ( qporq ) then
        fx = cum - p
      else
        fx = ccum - q
      end if

      call dinvr ( status, a, fx, qleft, qhi )
      go to 410

    else if ( status == -1 ) then

      if ( qleft ) then
        status = 1
        bound = 0.0
      else
        status = 2
        bound = inf
      end if

    end if

  else if ( which == 4 ) then

    b = 5.0d0
    call dstinv(0.0d0,inf,0.5d0,0.5d0,5.0d0,atol,tol)
    status = 0
    call dinvr ( status, b, fx, qleft, qhi )

  480   continue

    if ( status == 1 ) then

      call cumbet(x,y,a,b,cum,ccum)

      if ( qporq ) then
        fx = cum - p
      else
        fx = ccum - q
      end if

      call dinvr ( status, b, fx, qleft, qhi )
      go to 480

    else if ( status == -1 ) then

      if ( qleft ) then
        status = 1
        bound = 0.0
      else
        status = 2
        bound = inf
      end if

    end if

  end if

  return
end
subroutine cdfbin ( which, p, q, s, xn, pr, ompr, status, bound )
!
!*******************************************************************************
!
!! CDFBIN evaluates the Cumulative Distribution Function of the Binomial distribution.
!
!
!  Method:
!
!    Computation of other parameters involve a seach for a value that
!    produces the desired value of P.   The search relies on the
!    monotonicity of P with the other parameter.
!
!  Reference:
!
!    Abramowitz and Stegun,  
!    Handbook of Mathematical Functions 
!    1966, Formula 26.5.24.
!
!  Parameters:
!
!     WHICH --> Integer indicating which of the next four argument
!               values is to be calculated from the others.
!               Legal range: 1..4
!               iwhich = 1 : Calculate P and Q from S,XN,PR and OMPR
!               iwhich = 2 : Calculate S from P,Q,XN,PR and OMPR
!               iwhich = 3 : Calculate XN from P,Q,S,PR and OMPR
!               iwhich = 4 : Calculate PR and OMPR from P,Q,S and XN
!                    integer WHICH
!
!     P <--> The cumulation from 0 to S of the binomial distribution.
!            (Probablility of S or fewer successes in XN trials each
!            with probability of success PR.)
!            Input range: [0,1].
!                    double precision P
!
!     Q <--> 1-P.
!            Input range: [0, 1].
!            P + Q = 1.0.
!                    double precision Q
!
!     S <--> The number of successes observed.
!            Input range: [0, XN]
!            Search range: [0, XN]
!                    double precision S
!
!     XN  <--> The number of binomial trials.
!              Input range: (0, +infinity).
!              Search range: [1E-300, 1E300]
!                    double precision XN
!
!     PR  <--> The probability of success in each binomial trial.
!              Input range: [0,1].
!              Search range: [0,1]
!                    double precision PR
!
!     OMPR  <--> 1-PR
!              Input range: [0,1].
!              Search range: [0,1]
!              PR + OMPR = 1.0
!                    double precision OMPR
!
!     STATUS <-- 0 if calculation completed correctly
!               -I if input parameter number I is out of range
!                1 if answer appears to be lower than lowest
!                  search bound
!                2 if answer appears to be higher than greatest
!                  search bound
!                3 if P + Q /= 1
!                4 if PR + OMPR /= 1
!                    integer STATUS
!
!     BOUND <-- Undefined if STATUS is 0
!
!               Bound exceeded by parameter number I if STATUS
!               is negative.
!
!               Lower search bound if STATUS is 1.
!
!               Upper search bound if STATUS is 2.
!
!
  double precision, parameter :: atol = 1.0d-50
  double precision, parameter :: tol = 1.0d-8
  double precision, parameter :: inf = 1.0d300
!
  double precision bound
  double precision ccum
  double precision cum
  double precision dpmpar
  double precision fx
  double precision ompr
  double precision p
  double precision pq
  double precision pr
  double precision prompr
  double precision q
  logical qhi
  logical qleft
  logical qporq
  double precision s
  integer status
  integer which
  double precision xhi
  double precision xlo
  double precision xn
!
!     Check arguments
!
  if ( which < 1 ) then
    bound = 1.0d0
    status = - 1
    return
  end if

  if ( which > 4 ) then
    bound = 4.0d0
    status = -1
    return
  end if
!
!     P
!
  if ( which /= 1 ) then

    if (.not. ((p<0.0d0).or. (p>1.0d0))) go to 60
    if (.not. (p<0.0d0)) go to 40
    bound = 0.0d0
    go to 50

   40   bound = 1.0d0
   50   status = -2
    return

   60   continue
  end if
!
!     Q
!
  if ( which /= 1 ) then

  if (.not. ((q<0.0d0).or. (q>1.0d0))) go to 100
  if (.not. (q<0.0d0)) go to 80
  bound = 0.0d0
  go to 90

   80 bound = 1.0d0
   90 status = -3
  return

  100 continue
  end if
!
!     XN
!
  if ( which /= 3 ) then
    if ( xn > 0.0 ) then
      bound = 0.0d0
      status = -5
      return
    end if
  end if
!
!     S
!
  if ( which /= 2 ) then
    if (.not. ((s<0.0d0).or. ((which/=3).and. (s>xn)))) then
      go to 160
    end if
    if (.not. (s<0.0d0)) go to 140
    bound = 0.0d0
    go to 150

  140   bound = xn
  150   status = -4
    return

  160   continue
  end if
!
!     PR
!
  if ( which /= 4 ) then
    if (.not. ((pr<0.0d0).or. (pr>1.0d0))) go to 200
    if (.not. (pr<0.0d0)) go to 180
    bound = 0.0d0
    go to 190

  180   bound = 1.0d0
  190   status = -6
    return

  200   continue
  end if
!
!     OMPR
!
  if ( which /= 4 ) then
    if (.not. ((ompr<0.0d0).or. (ompr>1.0d0))) go to 240
    if (.not. (ompr<0.0d0)) go to 220
    bound = 0.0d0
    go to 230

  220   bound = 1.0d0
  230   status = -7
    return

  240   continue
  end if
!
!     P + Q
!
  if ( which /= 1 ) then
    pq = p + q
    if (.not. (abs(((pq)-0.5d0)-0.5d0)> (3.0d0* epsilon ( 1.0d0 ) ))) then
      go to 280
    end if

    if (.not. (pq<0.0d0)) go to 260
    bound = 0.0d0
    go to 270

  260   bound = 1.0d0
  270   status = 3
    return

  280   continue
  end if

  if (which==4) go to 330
!
!     PR + OMPR
!
  prompr = pr + ompr

  if (.not. (abs(((prompr)-0.5d0)-0.5d0)> (3.0d0* epsilon ( 1.0d0 ) ))) then
    go to 320
  end if

  if (.not. (prompr<0.0d0)) go to 300
  bound = 0.0d0
  go to 310

  300 bound = 1.0d0
  310 status = 4
  return

  320 continue
  330 continue

  if ( which /= 1 ) then
    qporq = p <= q
  end if
!
!     Calculate ANSWERS
!
  if ( which == 1 ) then

    call cumbin(s,xn,pr,ompr,p,q)
    status = 0

  else if ( which == 2 ) then

      s = 5.0d0
      call dstinv(0.0d0,xn,0.5d0,0.5d0,5.0d0,atol,tol)
      status = 0
      fx = 0.0
      call dinvr ( status, s, fx, qleft, qhi )

  340     if (.not. (status==1)) go to 370
      call cumbin(s,xn,pr,ompr,cum,ccum)

      if ( qporq ) then
        fx = cum - p
      else
        fx = ccum - q   
      end if

      call dinvr ( status, s, fx, qleft, qhi )
      go to 340

  370     if (.not. (status==-1)) go to 400
      if (.not. (qleft)) go to 380
      status = 1
      bound = 0.0d0
      go to 390

  380     status = 2
      bound = xn
  390     continue
  400     continue

  else if ( which == 3 ) then
!
!     Calculating XN
!
      xn = 5.0d0
      call dstinv(0.0d0,inf,0.5d0,0.5d0,5.0d0,atol,tol)
      status = 0
      call dinvr ( status, xn, fx, qleft, qhi )

  410     if (.not. (status==1)) go to 440
      call cumbin(s,xn,pr,ompr,cum,ccum)

      if ( qporq ) then
        fx = cum - p
      else
        fx = ccum - q
      end if

      call dinvr ( status, xn, fx, qleft, qhi )
      go to 410

  440     if (.not. (status==-1)) go to 470
      if (.not. (qleft)) go to 450
      status = 1
      bound = 0.0
      go to 460

  450     status = 2
      bound = inf
  460     continue
  470     continue

  else if ( which == 4 ) then
!
!     Calculating PR and OMPR
!
      call dstzr(0.0d0,1.0d0,atol,tol)
      if (.not. (qporq)) go to 500
      status = 0
      call dzror(status,pr,fx,xlo,xhi,qleft,qhi)
      ompr = 1.0 - pr
  480     continue

      if ( status == 1 ) then
        call cumbin ( s, xn, pr, ompr, cum, ccum )
        fx = cum - p
        call dzror ( status, pr, fx, xlo, xhi, qleft, qhi )
        ompr = 1.0 - pr
        go to 480
      end if

      go to 530

  500     status = 0
      call dzror(status,ompr,fx,xlo,xhi,qleft,qhi)
      pr = 1.0 - ompr

  510     continue

      if ( status == 1 ) then
        call cumbin(s,xn,pr,ompr,cum,ccum)
        fx = ccum - q
        call dzror(status,ompr,fx,xlo,xhi,qleft,qhi)
        pr = 1.0 - ompr
        go to 510
      end if

  530     if (.not. (status==-1)) go to 560

      if ( qleft ) then
        status = 1
        bound = 0.0d0
      else
        status = 2
        bound = 1.0d0
      end if

  560   end if

  return
end
subroutine cdfchi ( which, p, q, x, df, status, bound )
!
!*******************************************************************************
!
!! CDFCHI evaluates the Cumulative Distribution Function of the CHI-Square distribution.
!
!
!  Method:
!
!     Computation of other parameters involve a seach for a value that
!     produces  the desired  value  of P.   The search relies  on  the
!     monotonicity of P with the other parameter.
!
!  Reference:
!
!    Abramowitz and Stegun,  
!    Handbook of Mathematical Functions 
!    1966, Formula 26.4.19.
!
!  Parameters:
!
!     WHICH --> Integer indicating which of the next three argument
!               values is to be calculated from the others.
!               Legal range: 1..3
!               iwhich = 1 : Calculate P and Q from X and DF
!               iwhich = 2 : Calculate X from P,Q and DF
!               iwhich = 3 : Calculate DF from P,Q and X
!                    integer WHICH
!
!     P <--> The integral from 0 to X of the chi-square
!            distribution.
!            Input range: [0, 1].
!                    double precision P
!
!     Q <--> 1-P.
!            Input range: (0, 1].
!            P + Q = 1.0.
!                    double precision Q
!
!     X <--> Upper limit of integration of the non-central
!            chi-square distribution.
!            Input range: [0, +infinity).
!            Search range: [0,1E300]
!                    double precision X
!
!     DF <--> Degrees of freedom of the
!             chi-square distribution.
!             Input range: (0, +infinity).
!             Search range: [ 1E-300, 1E300]
!                    double precision DF
!
!     STATUS <-- 0 if calculation completed correctly
!               -I if input parameter number I is out of range
!                1 if answer appears to be lower than lowest
!                  search bound
!                2 if answer appears to be higher than greatest
!                  search bound
!                3 if P + Q /= 1
!               10 indicates error returned from cumgam.  See
!                  references in cdfgam
!                    integer STATUS
!
!     BOUND <-- Undefined if STATUS is 0
!
!               Bound exceeded by parameter number I if STATUS
!               is negative.
!
!               Lower search bound if STATUS is 1.
!
!               Upper search bound if STATUS is 2.
!
  double precision, parameter :: tol = 1.0d-8
  double precision, parameter :: atol = 1.0d-50
  double precision, parameter :: inf = 1.0d300
!
  double precision bound
  double precision ccum
  double precision cum
  double precision df
  double precision dpmpar
  double precision fx
  double precision p
  double precision porq
  double precision pq
  double precision q
  logical qhi
  logical qleft
  logical qporq
  integer status
  integer which
  double precision x
!
!     Check arguments
!
  if ( which < 1 ) then
    bound = 1.0d0
    status = - 1
    return
  end if

  if ( which > 3 ) then
    bound = 3.0d0
    status = -1
    return
  end if

  if (which==1) go to 70
!
!     P
!
  if (.not. ((p<0.0d0).or. (p>1.0d0))) go to 60
  if (.not. (p<0.0d0)) go to 40
  bound = 0.0d0
  go to 50

   40 bound = 1.0d0
   50 status = -2
  return

   60 continue
   70 if (which==1) go to 110
!
!     Q
!
  if (.not. ((q<=0.0d0).or. (q>1.0d0))) go to 100
  if ( q > 0.0d0 ) go to 80
  bound = 0.0d0
  go to 90

   80 bound = 1.0d0
   90 status = -3
  return

  100 continue
  110 if (which==2) go to 130
!
!     X
!
  if ( x < 0.0d0 ) then
    bound = 0.0d0
    status = -4
    return
  end if

  130 if (which==3) go to 150
!
!     DF
!
  if ( df<=0.0d0 ) then
    bound = 0.0d0
    status = -5
    return
  end if

  150 if (which==1) go to 190
!
!     P + Q
!
  pq = p + q
  if (.not. (abs(((pq)-0.5d0)-0.5d0)> (3.0d0* epsilon ( 1.0d0 ) ))) then
    go to 180
  end if

  if (.not. (pq<0.0d0)) go to 160
  bound = 0.0d0
  go to 170

  160 bound = 1.0d0
  170 status = 3
  return

  180 continue
  190 continue
!
!     Select the minimum of P or Q
!
  if ( which /= 1 ) then

    qporq = p <= q
    porq = min ( p, q )
   
  end if
!
!     Calculate ANSWERS
!
  if ( which == 1 ) then

    status = 0
    call cumchi ( x, df, p, q )

    if ( porq > 1.5d0 ) then
      status = 10
      return
    end if

  else if ( which == 2 ) then

    x = 5.0d0
    call dstinv(0.0d0,inf,0.5d0,0.5d0,5.0d0,atol,tol)
    status = 0
    fx = 0.0
    call dinvr ( status, x, fx, qleft, qhi )

  230   continue

    if ( status == 1 ) then

      call cumchi(x,df,cum,ccum)

      if ( qporq ) then
        fx = cum - p
      else
        fx = ccum - q
      end if

      if ( fx + porq > 1.5d0 ) then
        status = 10
        return
      end if

      call dinvr ( status, x, fx, qleft, qhi )
      go to 230

    else if ( status == - 1 ) then

      if ( qleft ) then
        status = 1
        bound = 0.0d0
      else
        status = 2
        bound = inf
      end if

    end if

  else if ( which == 3 ) then

    df = 5.0d0
    call dstinv(0.0d0,inf,0.5d0,0.5d0,5.0d0,atol,tol)
    status = 0
    call dinvr ( status, df, fx, qleft, qhi )

  310   continue

    if ( status == 1 ) then
 
      call cumchi(x,df,cum,ccum)

      if ( qporq ) then
        fx = cum - p
      else
        fx = ccum - q
      end if

      if ( fx + porq > 1.5d0 ) then
        status = 10
        return
      end if

      call dinvr ( status, df, fx, qleft, qhi )
      go to 310

    else if ( status == - 1 ) then

      if ( qleft ) then
        status = 1
        bound = 0.0
      else
        status = 2
        bound = inf
      end if

    end if

  end if

  return
end
subroutine cdfchn ( which, p, q, x, df, pnonc, status, bound )
!
!*******************************************************************************
!
!! CDFCHN evaluates the Cumulative Distribution Function of the Non-central Chi-Square.
!
!
!  Function:
!
!    Calculates any one parameter of the non-central chi-square
!    distribution given values for the others.
!
!  Method:
!
!     Computation of other parameters involve a seach for a value that
!     produces  the desired  value  of P.   The search relies  on  the
!     monotonicity of P with the other parameter.
!
!  WARNING:
!
!     The computation time  required for this  routine is proportional
!     to the noncentrality  parameter  (PNONC).  Very large  values of
!     this parameter can consume immense  computer resources.  This is
!     why the search range is bounded by 10,000.
!
!  Reference:
!
!    Abramowitz and Stegun,  
!    Handbook of Mathematical Functions 
!    1966, Formula 26.5.25.
!
!  Parameters:
!
!
!     WHICH --> Integer indicating which of the next three argument
!               values is to be calculated from the others.
!               Input range: 1..4
!               iwhich = 1 : Calculate P and Q from X and DF
!               iwhich = 2 : Calculate X from P,DF and PNONC
!               iwhich = 3 : Calculate DF from P,X and PNONC
!               iwhich = 4 : Calculate PNONC from P,X and DF
!                    integer WHICH
!
!     P <--> The integral from 0 to X of the non-central chi-square
!            distribution.
!            Input range: [0, 1-1E-16).
!                    double precision P
!
!     Q <--> 1-P.
!            Q is not used by this subroutine and is only included
!            for similarity with other cdf* routines.
!                    double precision Q
!
!     X <--> Upper limit of integration of the non-central
!            chi-square distribution.
!            Input range: [0, +infinity).
!            Search range: [0,1E300]
!                    double precision X
!
!     DF <--> Degrees of freedom of the non-central
!             chi-square distribution.
!             Input range: (0, +infinity).
!             Search range: [ 1E-300, 1E300]
!                    double precision DF
!
!     PNONC <--> Non-centrality parameter of the non-central
!                chi-square distribution.
!                Input range: [0, +infinity).
!                Search range: [0,1E4]
!                    double precision PNONC
!
!     STATUS <-- 0 if calculation completed correctly
!               -I if input parameter number I is out of range
!                1 if answer appears to be lower than lowest
!                  search bound
!                2 if answer appears to be higher than greatest
!                  search bound
!                    integer STATUS
!
!     BOUND <-- Undefined if STATUS is 0
!
!               Bound exceeded by parameter number I if STATUS
!               is negative.
!
!               Lower search bound if STATUS is 1.
!
!               Upper search bound if STATUS is 2.
!
  double precision, parameter :: tent4=1.0d4
  double precision, parameter :: tol=1.0d-8
  double precision, parameter :: atol=1.0d-50
  double precision, parameter :: inf=1.0d300
!
  double precision bound
  double precision df
  double precision p,q,pnonc,x
  integer status,which
  double precision fx,cum,ccum
  logical qhi,qleft
!
!     Check arguments
!
  if ( which < 1 ) then
    bound = 1.0d0
    status = - 1
    return
  end if

  if ( which > 4 ) then
    bound = 4.0d0
    status = -1
    return
  end if

  if (which==1) go to 70
!
!     P
!
  if (.not. ((p<0.0d0).or. (p>1.0))) go to 60
  if (.not. (p<0.0d0)) go to 40
  bound = 0.0d0
  go to 50

   40 bound = 1.0
   50 status = -2
  return

   60 continue
   70 if (which==2) go to 90
!
!     X
!
  if (.not. (x<0.0d0)) go to 80
  bound = 0.0d0
  status = -4
  return

   80 continue
   90 if (which==3) go to 110
!
!     DF
!
  if ( df <= 0.0d0 ) then
    bound = 0.0d0
    status = -5
    return
  end if

  110 if (which==4) go to 130
!
!     PNONC
!
  if ( pnonc < 0.0d0 ) then
    bound = 0.0d0
    status = -6
    return
  end if
!
!     Calculate ANSWERS
!
  130 if ( which == 1 ) then

    call cumchn ( x, df, pnonc, p, q )
    status = 0

  else if ( which == 2 ) then

      x = 5.0d0
      call dstinv(0.0d0,inf,0.5d0,0.5d0,5.0d0,atol,tol)
      status = 0
      fx = 0.0
      call dinvr ( status, x, fx, qleft, qhi )

  140     if ( status /= 1 ) go to 150
      call cumchn(x,df,pnonc,cum,ccum)
      fx = cum - p
      call dinvr ( status, x, fx, qleft, qhi )
      go to 140

  150     if ( status == -1 ) then

        if ( qleft ) then
          status = 1
          bound = 0.0d0
        else
          status = 2
          bound = inf
        end if

      end if

  else if ( which == 3 ) then

      df = 5.0d0
      call dstinv(0.0d0,inf,0.5d0,0.5d0,5.0d0,atol,tol)
      status = 0
      call dinvr ( status, df, fx, qleft, qhi )

  190     if ( status /= 1 ) go to 200
      call cumchn(x,df,pnonc,cum,ccum)
      fx = cum - p
      call dinvr ( status, df, fx, qleft, qhi )
      go to 190

  200     if (.not. (status==-1)) go to 230
      if (.not. (qleft)) go to 210
      status = 1
      bound = 0.0
      go to 220

  210     status = 2
      bound = inf
  220     continue
  230     continue

  else if ( which == 4 ) then

      pnonc = 5.0d0
      call dstinv(0.0d0,tent4,0.5d0,0.5d0,5.0d0,atol,tol)
      status = 0
      call dinvr ( status, pnonc, fx, qleft, qhi )

  240     if (.not. (status==1)) go to 250
      call cumchn(x,df,pnonc,cum,ccum)
      fx = cum - p
      call dinvr ( status, pnonc, fx, qleft, qhi )
      go to 240

  250     if (.not. (status==-1)) go to 280

      if ( qleft ) then
        status = 1
        bound = 0.0
      else
        status = 2
        bound = tent4
      end if

  280   end if

  return
end
subroutine cdff ( which, p, q, f, dfn, dfd, status, bound )
!
!*******************************************************************************
!
!! CDFF evaluates the Cumulative Distribution Function of the F distribution.
!
!
!  Method:
!
!     Computation of other parameters involve a seach for a value that
!     produces  the desired  value  of P.   The search relies  on  the
!     monotonicity of P with the other parameter.
!
!  WARNING
!
!     The value of the  cumulative  F distribution is  not necessarily
!     monotone in  either degrees of freedom.  There  thus may  be two
!     values  that  provide a given CDF  value.   This routine assumes
!     monotonicity and will find an arbitrary one of the two values.
!
!  Reference:
!
!    Abramowitz and Stegun,  
!    Handbook of Mathematical Functions 
!    1966, Formula 26.6.2.
!
!  Parameters:
!
!     WHICH --> Integer indicating which of the next four argument
!               values is to be calculated from the others.
!               Legal range: 1..4
!               iwhich = 1 : Calculate P and Q from F,DFN and DFD
!               iwhich = 2 : Calculate F from P,Q,DFN and DFD
!               iwhich = 3 : Calculate DFN from P,Q,F and DFD
!               iwhich = 4 : Calculate DFD from P,Q,F and DFN
!                    integer WHICH
!
!       P <--> The integral from 0 to F of the f-density.
!              Input range: [0,1].
!                    double precision P
!
!       Q <--> 1-P.
!              Input range: (0, 1].
!              P + Q = 1.0.
!                    double precision Q
!
!       F <--> Upper limit of integration of the f-density.
!              Input range: [0, +infinity).
!              Search range: [0,1E300]
!                    double precision F
!
!     DFN < --> Degrees of freedom of the numerator sum of squares.
!               Input range: (0, +infinity).
!               Search range: [ 1E-300, 1E300]
!                    double precision DFN
!
!     DFD < --> Degrees of freedom of the denominator sum of squares.
!               Input range: (0, +infinity).
!               Search range: [ 1E-300, 1E300]
!                    double precision DFD
!
!     STATUS <-- 0 if calculation completed correctly
!               -I if input parameter number I is out of range
!                1 if answer appears to be lower than lowest
!                  search bound
!                2 if answer appears to be higher than greatest
!                  search bound
!                3 if P + Q /= 1
!                    integer STATUS
!
!     BOUND <-- Undefined if STATUS is 0
!
!               Bound exceeded by parameter number I if STATUS
!               is negative.
!
!               Lower search bound if STATUS is 1.
!
!               Upper search bound if STATUS is 2.
!
  double precision, parameter :: tol=1.0d-8
  double precision, parameter :: atol=1.0d-50
  double precision, parameter :: inf=1.0d300
!
  double precision bound
  double precision dfd,dfn,f,p,q
  integer status,which
  double precision pq,fx,cum,ccum
  logical qhi,qleft,qporq
  double precision dpmpar
!
  if ( which < 1 ) then
    bound = 1.0d0
    status = - 1
    return
  end if

  if ( which > 4 ) then
    bound = 4.0d0
    status = -1
    return
  end if

  if (which==1) go to 70
!
!     P
!
  if (.not. ((p<0.0d0).or. (p>1.0d0))) go to 60
  if (.not. (p<0.0d0)) go to 40
  bound = 0.0d0
  go to 50

   40 bound = 1.0d0
   50 status = -2
  return

   60 continue
   70 if (which==1) go to 110
!
!     Q
!
  if (.not. ((q<=0.0d0).or. (q>1.0d0))) go to 100
  if (.not. (q<=0.0d0)) go to 80
  bound = 0.0d0
  go to 90

   80 bound = 1.0d0
   90 status = -3
  return

  100 continue
  110 if (which==2) go to 130
!
!     F
!
  if ( f<0.0d0 ) then
    bound = 0.0d0
    status = -4
    return
  end if

  130 continue
!
!     DFN
!
  if ( which /= 3 ) then

    if ( dfn <= 0.0d0 ) then
      bound = 0.0d0
      status = -5
      return
    end if

  end if

  if (which==4) go to 170
!
!     DFD
!
  if (.not. (dfd<=0.0d0)) go to 160
  bound = 0.0d0
  status = -6
  return

  160 continue
  170 if (which==1) go to 210
!
!     P + Q
!
  pq = p + q
  if (.not. (abs(((pq)-0.5d0)-0.5d0)> (3.0d0* epsilon ( 1.0d0 ) ))) then
    go to 200
  end if

  if (.not. (pq<0.0d0)) go to 180
  bound = 0.0d0
  go to 190

  180 bound = 1.0d0
  190 status = 3
  return

  200 continue
  210 continue

  if ( which /= 1 ) then
    qporq = p <= q
  end if
!
!     Calculate ANSWERS
!
  if ( which == 1 ) then

    call cumf ( f, dfn, dfd, p, q )
    status = 0

  else if ( which == 2 ) then

    f = 5.0d0
    call dstinv(0.0d0,inf,0.5d0,0.5d0,5.0d0,atol,tol)
    status = 0
    fx = 0.0
    call dinvr ( status, f, fx, qleft, qhi )

  220   continue

    if ( status == 1 ) then
      call cumf ( f, dfn, dfd, cum, ccum )
      if ( qporq ) then
        fx = cum - p
      else
        fx = ccum - q
      end if
      call dinvr ( status, f, fx, qleft, qhi )
      go to 220
    end if

    if (.not. (status==-1)) go to 280
      if (.not. (qleft)) go to 260
      status = 1
      bound = 0.0d0
      go to 270

  260     status = 2
      bound = inf
  270     continue
  280     continue

  else if ( which == 3 ) then

      dfn = 5.0d0
      call dstinv(0.0d0,inf,0.5d0,0.5d0,5.0d0,atol,tol)
      status = 0
      call dinvr ( status, dfn, fx, qleft, qhi )

  290     continue

      if (.not. (status==1)) go to 320
      call cumf ( f, dfn, dfd, cum, ccum )

      if ( qporq ) then
        fx = cum - p
      else
        fx = ccum - q
      end if

      call dinvr ( status, dfn, fx, qleft, qhi )
      go to 290

  320     if (.not. (status==-1)) go to 350
      if (.not. (qleft)) go to 330
      status = 1
      bound = 0.0
      go to 340

  330     status = 2
      bound = inf
  340     continue
  350     continue

  else if ( which == 4 ) then

      dfd = 5.0d0
      call dstinv(0.0d0,inf,0.5d0,0.5d0,5.0d0,atol,tol)
      status = 0
      call dinvr ( status, dfd, fx, qleft, qhi )

  360     continue

      if (.not. (status==1)) go to 390
      call cumf ( f, dfn, dfd, cum, ccum )

      if ( qporq ) then
        fx = cum - p
      else
        fx = ccum - q
      end if

      call dinvr ( status, dfd, fx, qleft, qhi )
      go to 360

  390     if (.not. (status==-1)) go to 420

      if ( qleft ) then
        status = 1
        bound = 0.0
      else
        status = 2
        bound = inf
      end if

  420   end if

  return
end
subroutine cdffnc ( which, p, q, f, dfn, dfd, phonc, status, bound )
!
!*******************************************************************************
!
!! CDFFNC evaluates the Cumulative Distribution Function of the Non-central F distribution.
!
!
!
!  Method:
!
!     Computation of other parameters involve a seach for a value that
!     produces  the desired  value  of P.   The search relies  on  the
!     monotonicity of P with the other parameter.
!
!  WARNING
!
!     The computation time  required for this  routine is proportional
!     to the noncentrality  parameter  (PNONC).  Very large  values of
!     this parameter can consume immense  computer resources.  This is
!     why the search range is bounded by 10,000.
!
!     The  value  of the  cumulative  noncentral F distribution is not
!     necessarily monotone in either degrees  of freedom.  There  thus
!     may be two values that provide a given  CDF value.  This routine
!     assumes monotonicity  and will find  an arbitrary one of the two
!     values.
!
!  Reference:
!
!    Abramowitz and Stegun,  
!    Handbook of Mathematical Functions 
!    1966, Formula 26.6.20.
!
!  Parameters:
!
!     WHICH --> Integer indicating which of the next five argument
!               values is to be calculated from the others.
!               Legal range: 1..5
!               iwhich = 1 : Calculate P and Q from F,DFN,DFD and PNONC
!               iwhich = 2 : Calculate F from P,Q,DFN,DFD and PNONC
!               iwhich = 3 : Calculate DFN from P,Q,F,DFD and PNONC
!               iwhich = 4 : Calculate DFD from P,Q,F,DFN and PNONC
!               iwhich = 5 : Calculate PNONC from P,Q,F,DFN and DFD
!                    integer WHICH
!
!       P <--> The integral from 0 to F of the non-central f-density.
!              Input range: [0,1-1E-16).
!                    double precision P
!
!       Q <--> 1-P.
!            Q is not used by this subroutine and is only included
!            for similarity with other cdf* routines.
!                    double precision Q
!
!       F <--> Upper limit of integration of the non-central f-density.
!              Input range: [0, +infinity).
!              Search range: [0,1E300]
!                    double precision F
!
!     DFN < --> Degrees of freedom of the numerator sum of squares.
!               Input range: (0, +infinity).
!               Search range: [ 1E-300, 1E300]
!                    double precision DFN
!
!     DFD < --> Degrees of freedom of the denominator sum of squares.
!               Must be in range: (0, +infinity).
!               Input range: (0, +infinity).
!               Search range: [ 1E-300, 1E300]
!                    double precision DFD
!
!     PNONC <-> The non-centrality parameter
!               Input range: [0,infinity)
!               Search range: [0,1E4]
!                    double precision PHONC
!
!     STATUS <-- 0 if calculation completed correctly
!               -I if input parameter number I is out of range
!                1 if answer appears to be lower than lowest
!                  search bound
!                2 if answer appears to be higher than greatest
!                  search bound
!                3 if P + Q /= 1
!                    integer STATUS
!
!     BOUND <-- Undefined if STATUS is 0
!
!               Bound exceeded by parameter number I if STATUS
!               is negative.
!
!               Lower search bound if STATUS is 1.
!
!               Upper search bound if STATUS is 2.
!
  double precision, parameter :: tent4=1.0d4
  double precision, parameter :: tol=1.0d-8
  double precision, parameter :: atol=1.0d-50
  double precision, parameter :: inf=1.0d300
!
  double precision bound,dfd,dfn,f,p,q,phonc
  integer status,which
  double precision fx,cum,ccum
  logical qhi,qleft
!
!     Check arguments
!
  if ( which < 1 ) then
    bound = 1.0d0
    status = - 1
    return
  end if

  if ( which > 5 ) then
    bound = 5.0d0
    status = -1
    return
  end if

  if (which==1) go to 70
!
!     P
!
  if (.not. ((p<0.0d0).or. (p>1.0))) go to 60
  if (.not. (p<0.0d0)) go to 40
  bound = 0.0d0
  go to 50

   40 bound = 1.0
   50 status = -2
  return

   60 continue
   70 if (which==2) go to 90
!
!     F
!
  if ( f < 0.0d0 ) then
    bound = 0.0d0
    status = -4
    return
  end if

   90 if (which==3) go to 110
!
!     DFN
!
  if (.not. (dfn<=0.0d0)) go to 100
  bound = 0.0d0
  status = -5
  return

  100 continue
  110 if (which==4) go to 130
!
!     DFD
!
  if ( dfd <= 0.0d0 ) then
    bound = 0.0d0
    status = -6
    return
  end if

  130 if (which==5) go to 150
!
!     PHONC
!
  if ( phonc<0.0d0 ) then
    bound = 0.0d0
    status = -7
    return
  end if
!
!     Calculate ANSWERS
!
  150 if ( which == 1 ) then

    call cumfnc ( f, dfn, dfd, phonc, p, q )
    status = 0

  else if ( which == 2 ) then

      f = 5.0d0
      call dstinv(0.0d0,inf,0.5d0,0.5d0,5.0d0,atol,tol)
      status = 0
      fx = 0.0
      call dinvr ( status, f, fx, qleft, qhi )

  160     continue

      if ( status /= 1 ) go to 170
      call cumfnc ( f, dfn, dfd, phonc, cum, ccum )
      fx = cum - p
      call dinvr ( status, f, fx, qleft, qhi )
      go to 160

  170     if ( status /= -1 ) go to 200
      if ( .not. qleft ) go to 180
      status = 1
      bound = 0.0d0
      go to 190

  180     status = 2
      bound = inf
  190     continue
  200     continue

  else if ( which == 3 ) then

      dfn = 5.0d0
      call dstinv(0.0d0,inf,0.5d0,0.5d0,5.0d0,atol,tol)
      status = 0
      call dinvr ( status, dfn, fx, qleft, qhi )

  210     continue

      if ( status /= 1 ) go to 220
      call cumfnc ( f, dfn, dfd, phonc, cum, ccum )
      fx = cum - p
      call dinvr ( status, dfn, fx, qleft, qhi )
      go to 210

  220     if ( status /= -1 ) go to 250
      if (.not. qleft ) go to 230
      status = 1
      bound = 0.0
      go to 240

  230     status = 2
      bound = inf
  240     continue
  250     continue

  else if ( which == 4 ) then

      dfd = 5.0d0
      call dstinv(0.0d0,inf,0.5d0,0.5d0,5.0d0,atol,tol)
      status = 0
      call dinvr ( status, dfd, fx, qleft, qhi )

  260     continue

      if ( status /= 1 ) go to 270
      call cumfnc ( f, dfn, dfd, phonc, cum, ccum )
      fx = cum - p
      call dinvr ( status, dfd, fx, qleft, qhi )
      go to 260

  270     if ( status /= -1 ) go to 300
      if (.not. qleft ) go to 280
      status = 1
      bound = 0.0
      go to 290

  280     status = 2
      bound = inf
  290     continue
  300     continue

  else if ( which == 5 ) then

      phonc = 5.0d0
      call dstinv(0.0d0,tent4,0.5d0,0.5d0,5.0d0,atol,tol)
      status = 0
      call dinvr ( status, phonc, fx, qleft, qhi )

  310     if (.not. (status==1)) go to 320
      call cumfnc ( f, dfn, dfd, phonc, cum, ccum )
      fx = cum - p
      call dinvr ( status, phonc, fx, qleft, qhi )
      go to 310

  320     if ( status /= -1 ) go to 350
      if (.not. qleft ) go to 330
      status = 1
      bound = 0.0d0
      go to 340

  330     status = 2
      bound = tent4
  340     continue
  350   end if

  return
end
subroutine cdfgam ( which, p, q, x, shape, scale, status, bound )
!
!*******************************************************************************
!
!! CDFGAM evaluates the Cumulative Distribution Function of the GAMma Distribution.
!
!
!  Parameters:
!
!
!     WHICH --> Integer indicating which of the next four argument
!               values is to be calculated from the others.
!               Legal range: 1..4
!               iwhich = 1 : Calculate P and Q from X,SHAPE and SCALE
!               iwhich = 2 : Calculate X from P,Q,SHAPE and SCALE
!               iwhich = 3 : Calculate SHAPE from P,Q,X and SCALE
!               iwhich = 4 : Calculate SCALE from P,Q,X and SHAPE
!                    integer WHICH
!
!     P <--> The integral from 0 to X of the gamma density.
!            Input range: [0,1].
!                    double precision P
!
!     Q <--> 1-P.
!            Input range: (0, 1].
!            P + Q = 1.0.
!                    double precision Q
!
!
!     X <--> The upper limit of integration of the gamma density.
!            Input range: [0, +infinity).
!            Search range: [0,1E300]
!                    double precision X
!
!     SHAPE <--> The shape parameter of the gamma density.
!                Input range: (0, +infinity).
!                Search range: [1E-300,1E300]
!                  double precision SHAPE
!
!
!     SCALE <--> The scale parameter of the gamma density.
!                Input range: (0, +infinity).
!                Search range: (1E-300,1E300]
!                   double precision SCALE
!
!     STATUS <-- 0 if calculation completed correctly
!               -I if input parameter number I is out of range
!                1 if answer appears to be lower than lowest
!                  search bound
!                2 if answer appears to be higher than greatest
!                  search bound
!                3 if P + Q /= 1
!                10 if the gamma or inverse gamma routine cannot
!                   compute the answer.  Usually happens only for
!                   X and SHAPE very large (gt 1E10 or more)
!                    integer STATUS
!
!     BOUND <-- Undefined if STATUS is 0
!
!               Bound exceeded by parameter number I if STATUS
!               is negative.
!
!               Lower search bound if STATUS is 1.
!
!               Upper search bound if STATUS is 2.
!
!
!                              Method
!
!     Cumulative distribution function (P) is calculated directly by
!     the code associated with:
!
!     DiDinato, A. R. and Morris, A. H. Computation of the  incomplete
!     gamma function  ratios  and their  inverse.   ACM  Trans.  Math.
!     Softw. 12 (1986), 377-393.
!
!     Computation of other parameters involve a seach for a value that
!     produces  the desired  value  of P.   The search relies  on  the
!     monotonicity of P with the other parameter.
!
!
!                              Note
!
!     The gamma density is proportional to
!       T**(SHAPE - 1) * EXP(- SCALE * T)
!
  double precision, parameter :: tol=1.0d-8
  double precision, parameter :: atol=1.0d-50
  double precision, parameter :: inf=1.0d300
!
  double precision bound
  double precision p,q,scale,shape,x
  integer status,which
  double precision xx
  double precision fx,xscale,cum,ccum,pq,porq
  integer ierr
  logical qhi,qleft,qporq
  double precision dpmpar
!
!     Check arguments
!
  if ( which < 1 ) then
    bound = 1.0d0
    status = - 1
    return
  end if

  if ( which > 4 ) then
    bound = 4.0d0
    status = -1
    return
  end if
!
!     P
!
  if ( which /= 1 ) then

    if (.not. ((p<0.0d0).or. (p>1.0d0))) go to 60
    if (.not. (p<0.0d0)) go to 40
    bound = 0.0d0
    go to 50

   40   bound = 1.0d0
   50   status = -2
    return

   60   continue
  end if

  if (which==1) go to 110
!
!     Q
!
  if (.not. ((q<=0.0d0).or. (q>1.0d0))) go to 100
  if (.not. (q<=0.0d0)) go to 80
  bound = 0.0d0
  go to 90

   80 bound = 1.0d0
   90 status = -3
  return

  100 continue
  110 if (which==2) go to 130
!
!     X
!
  if (.not. (x<0.0d0)) go to 120
  bound = 0.0d0
  status = -4
  return

  120 continue
  130 if (which==3) go to 150
!
!     SHAPE
!
  if ( shape <= 0.0d0 ) then
    bound = 0.0d0
    status = -5
    return
  end if

  140 continue
  150 if (which==4) go to 170
!
!     SCALE
!
  if ( scale <= 0.0d0 ) then
    bound = 0.0d0
    status = -6
    return
  end if

  170 if (which==1) go to 210
!
!     P + Q
!
  pq = p + q
  if (.not. (abs(((pq)-0.5d0)-0.5d0)> (3.0d0* epsilon ( 1.0d0 ) ))) then
    go to 200
  end if

  if (.not. (pq<0.0d0)) go to 180
  bound = 0.0d0
  go to 190

  180 bound = 1.0d0
  190 status = 3
  return

  200 continue
  210 if (which==1) go to 240
!
!  Select the minimum of P or Q
!
  qporq = p <= q
  if ( .not. qporq ) go to 220
  porq = p
  go to 230

  220 porq = q
  230 continue
!
!  Calculate ANSWERS
!
  240 if ( which == 1 ) then

    status = 0
    xscale = x * scale
    call cumgam(xscale,shape,p,q)
    if (porq>1.5d0) status = 10

  else if ( which == 2 ) then

    call gamma_inc_inv ( shape, xx, -1.0d0, p, q, ierr )

    if (ierr<0.0d0) then
      status = 10
      return
    else
      x = xx/scale
      status = 0
    end if

  else if ( which == 3 ) then

    shape = 5.0d0
    xscale = x * scale
    call dstinv(0.0d0,inf,0.5d0,0.5d0,5.0d0,atol,tol)
    status = 0
    fx = 0.0
    call dinvr ( status, shape, fx, qleft, qhi )

  250   continue

    if ( status /= 1 ) go to 290
    call cumgam(xscale,shape,cum,ccum)

    if ( qporq ) then
      fx = cum - p
    else
      fx = ccum - q
    end if

    if (.not. ((qporq.and. (cum>1.5d0)).or. &
      ((.not.qporq).and. (ccum>1.5d0)))) then
      go to 280
    end if

    status = 10
    return

  280   call dinvr ( status, shape, fx, qleft, qhi )
    go to 250

  290   continue

    if ( status == -1 ) then

      if ( qleft ) then
        status = 1
        bound = 0.0
      else
        status = 2
        bound = inf
      end if

    end if

  else if ( which == 4 ) then

    call gamma_inc_inv ( shape, xx, -1.0d0, p, q, ierr )

    if ( ierr < 0.0d0 ) then
      status = 10
    else
      scale = xx / x
      status = 0
    end if

  end if

  return
end
subroutine cdfnbn ( which, p, q, s, xn, pr, ompr, status, bound )
!
!*******************************************************************************
!
!! CDFNBN evaluates the Cumulative Distribution Function of the Negative BiNomial distribution
!
!
!  Function:
!
!     Calculates any one parameter of the negative binomial
!     distribution given values for the others.
!
!     The  cumulative  negative   binomial  distribution  returns  the
!     probability that there  will be  F or fewer failures before  the
!     XNth success in binomial trials each of which has probability of
!     success PR.
!
!     The individual term of the negative binomial is the probability of
!     S failures before XN successes and is
!          Choose( S, XN+S-1 ) * PR^(XN) * (1-PR)^S
!
!  Reference:
!
!    Abramowitz and Stegun,  
!    Handbook of Mathematical Functions 
!    1966, Formula 26.5.26.
!
!  Parameters:
!
!     WHICH --> Integer indicating which of the next four argument
!               values is to be calculated from the others.
!               Legal range: 1..4
!               iwhich = 1 : Calculate P and Q from S,XN,PR and OMPR
!               iwhich = 2 : Calculate S from P,Q,XN,PR and OMPR
!               iwhich = 3 : Calculate XN from P,Q,S,PR and OMPR
!               iwhich = 4 : Calculate PR and OMPR from P,Q,S and XN
!                    integer WHICH
!
!     P <--> The cumulation from 0 to S of the  negative
!            binomial distribution.
!            Input range: [0,1].
!                    double precision P
!
!     Q <--> 1-P.
!            Input range: (0, 1].
!            P + Q = 1.0.
!                    double precision Q
!
!     S <--> The upper limit of cumulation of the binomial distribution.
!            There are F or fewer failures before the XNth success.
!            Input range: [0, +infinity).
!            Search range: [0, 1E300]
!                    double precision S
!
!     XN  <--> The number of successes.
!              Input range: [0, +infinity).
!              Search range: [0, 1E300]
!                    double precision XN
!
!     PR  <--> The probability of success in each binomial trial.
!              Input range: [0,1].
!              Search range: [0,1].
!                    double precision PR
!
!     OMPR  <--> 1-PR
!              Input range: [0,1].
!              Search range: [0,1]
!              PR + OMPR = 1.0
!                    double precision OMPR
!
!     STATUS <-- 0 if calculation completed correctly
!               -I if input parameter number I is out of range
!                1 if answer appears to be lower than lowest
!                  search bound
!                2 if answer appears to be higher than greatest
!                  search bound
!                3 if P + Q /= 1
!                4 if PR + OMPR /= 1
!                    integer STATUS
!
!     BOUND <-- Undefined if STATUS is 0
!
!               Bound exceeded by parameter number I if STATUS
!               is negative.
!
!               Lower search bound if STATUS is 1.
!
!               Upper search bound if STATUS is 2.
!
!
!                              Method
!
!     Computation of other parameters involve a seach for a value that
!     produces  the desired  value  of P.   The search relies  on  the
!     monotonicity of P with the other parameter.
!
  double precision, parameter :: tol=1.0d-8
  double precision, parameter :: atol=1.0d-50
  double precision, parameter :: inf=1.0d300
!
  double precision bound
  double precision p,q,pr,ompr,s,xn
  integer status,which
  double precision fx,xhi,xlo,pq,prompr,cum,ccum
  logical qhi,qleft,qporq
!
  double precision dpmpar
!
!     Check arguments
!
  if ( which < 1 ) then
    bound = 1.0d0
    status = - 1
    return
  end if

  if ( which > 4 ) then
    bound = 4.0d0
    status = -1
    return
  end if

  if (which==1) go to 70
!
!     P
!
  if (.not. ((p<0.0d0).or. (p>1.0d0))) go to 60
  if (.not. (p<0.0d0)) go to 40
  bound = 0.0d0
  go to 50

   40 bound = 1.0d0
   50 status = -2
  return

   60 continue
   70 if (which==1) go to 110
!
!     Q
!
  if (.not. ((q<=0.0d0).or. (q>1.0d0))) go to 100
  if (.not. (q<=0.0d0)) go to 80
  bound = 0.0d0
  go to 90

   80 bound = 1.0d0
   90 status = -3
  return

  100 continue
  110 if (which==2) go to 130
!
!     S
!
  if ( s < 0.0d0 ) then
    bound = 0.0d0
    status = -4
    return
  end if

  120 continue
  130 if (which==3) go to 150
!
!     XN
!
  if (.not. (xn<0.0d0)) go to 140
  bound = 0.0d0
  status = -5
  return

  140 continue
  150 if (which==4) go to 190
!
!     PR
!
  if (.not. ((pr<0.0d0).or. (pr>1.0d0))) go to 180
  if (.not. (pr<0.0d0)) go to 160
  bound = 0.0d0
  go to 170

  160 bound = 1.0d0
  170 status = -6
  return

  180 continue
  190 if (which==4) go to 230
!
!     OMPR
!
  if (.not. ((ompr<0.0d0).or. (ompr>1.0d0))) go to 220
  if (.not. (ompr<0.0d0)) go to 200
  bound = 0.0d0
  go to 210

  200 bound = 1.0d0
  210 status = -7
  return

  220 continue
  230 if (which==1) go to 270
!
!     P + Q
!
  pq = p + q
  if (.not. (abs(((pq)-0.5d0)-0.5d0)> (3.0d0* epsilon ( 1.0d0 ) ))) then
    go to 260
  end if

  if (.not. (pq<0.0d0)) go to 240
  bound = 0.0d0
  go to 250

  240 bound = 1.0d0
  250 status = 3
  return

  260 continue
  270 if (which==4) go to 310
!
!     PR + OMPR
!
  prompr = pr + ompr
  if (.not. (abs(((prompr)-0.5d0)-0.5d0)> (3.0d0* epsilon ( 1.0d0 ) ))) then
    go to 300
  end if

  if (.not. (prompr<0.0d0)) go to 280
  bound = 0.0d0
  go to 290

  280 bound = 1.0d0
  290 status = 4
  return

  300 continue
  310 continue

  if ( which /= 1 ) then
    qporq = p <= q
  end if
!
!  Calculate ANSWERS
!
  if ( which == 1 ) then

    call cumnbn(s,xn,pr,ompr,p,q)
    status = 0

  else if ( which == 2 ) then

    s = 5.0d0
    call dstinv(0.0d0,inf,0.5d0,0.5d0,5.0d0,atol,tol)
    status = 0
    fx = 0.0
    call dinvr ( status, s, fx, qleft, qhi )

  320   continue

    if ( status == 1 ) then

      call cumnbn(s,xn,pr,ompr,cum,ccum)

      if ( p <= q ) then
        fx = cum - p
      else
        fx = ccum - q
      end if

      call dinvr ( status, s, fx, qleft, qhi )
      go to 320

    else if ( status == -1 ) then

      if ( qleft ) then
        status = 1
        bound = 0.0d0
      else
        status = 2
        bound = inf
      end if

    end if

  else if ( which == 3 ) then
!
!     Calculating XN
!
      xn = 5.0d0
      call dstinv(0.0d0,inf,0.5d0,0.5d0,5.0d0,atol,tol)
      status = 0
      call dinvr ( status, xn, fx, qleft, qhi )

  390     continue

      if (.not. (status==1)) go to 420
      call cumnbn(s,xn,pr,ompr,cum,ccum)

      if ( qporq ) then
        fx = cum - p
      else
        fx = ccum - q
      end if

      call dinvr ( status, xn, fx, qleft, qhi )
      go to 390

  420     if (.not. (status==-1)) go to 450

      if ( qleft ) then
        status = 1
        bound = 0.0d0
      else
        status = 2
        bound = inf
      end if

  450     continue

  else if ( which == 4 ) then
!
!     Calculating PR and OMPR
!
    call dstzr(0.0d0,1.0d0,atol,tol)
    if (.not. (qporq)) go to 480
    status = 0
    call dzror(status,pr,fx,xlo,xhi,qleft,qhi)
    ompr = 1.0 - pr
  460   if (.not. (status==1)) go to 470
    call cumnbn(s,xn,pr,ompr,cum,ccum)
    fx = cum - p
    call dzror(status,pr,fx,xlo,xhi,qleft,qhi)
    ompr = 1.0 - pr
    go to 460

  470   go to 510

  480   status = 0
    call dzror(status,ompr,fx,xlo,xhi,qleft,qhi)
    pr = 1.0 - ompr

  490   continue

    if ( status == 1 ) then
      call cumnbn(s,xn,pr,ompr,cum,ccum)
      fx = ccum - q
      call dzror(status,ompr,fx,xlo,xhi,qleft,qhi)
      pr = 1.0 - ompr
      go to 490
    end if

  510   continue

    if ( status == - 1 ) then

      if ( qleft ) then
        status = 1
        bound = 0.0d0
      else
        status = 2
        bound = 1.0d0
      end if

    end if

  540   end if

  return
end
subroutine cdfnor ( which, p, q, x, mean, sd, status, bound )
!
!*******************************************************************************
!
!! CDFNOR evaluates the Cumulative Distribution Function of the Normal distribution
!
!
!  Reference:
!
!    Abramowitz and Stegun,  
!    Handbook of Mathematical Functions 
!    1966, Formula 26.2.12.
!
!  Parameters:
!
!     WHICH  --> Integer indicating  which of the  next  parameter
!     values is to be calculated using values  of the others.
!     Legal range: 1..4
!               iwhich = 1 : Calculate P and Q from X,MEAN and SD
!               iwhich = 2 : Calculate X from P,Q,MEAN and SD
!               iwhich = 3 : Calculate MEAN from P,Q,X and SD
!               iwhich = 4 : Calculate SD from P,Q,X and MEAN
!                    integer WHICH
!
!     P <--> The integral from -infinity to X of the normal density.
!            Input range: (0,1].
!                    double precision P
!
!     Q <--> 1-P.
!            Input range: (0, 1].
!            P + Q = 1.0.
!                    double precision Q
!
!     X < --> Upper limit of integration of the normal-density.
!             Input range: ( -infinity, +infinity)
!                    double precision X
!
!     MEAN <--> The mean of the normal density.
!               Input range: (-infinity, +infinity)
!                    double precision MEAN
!
!     SD <--> Standard Deviation of the normal density.
!             Input range: (0, +infinity).
!                    double precision SD
!
!     STATUS <-- 0 if calculation completed correctly
!               -I if input parameter number I is out of range
!                1 if answer appears to be lower than lowest
!                  search bound
!                2 if answer appears to be higher than greatest
!                  search bound
!                3 if P + Q /= 1
!                    integer STATUS
!
!     BOUND <-- Undefined if STATUS is 0
!               Bound exceeded by parameter number I if STATUS is negative.
!               Lower search bound if STATUS is 1.
!               Upper search bound if STATUS is 2.
!
!  Method
!
!     A slightly modified version of ANORM from
!
!     Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portable FORTRAN
!     Package of Special Function Routines and Test Drivers"
!     acm Transactions on Mathematical Software. 19, 22-32.
!
!     is used to calculate the  cumulative standard normal distribution.
!
!     The rational functions from pages  90-95  of Kennedy and Gentle,
!     Statistical  Computing,  Marcel  Dekker, NY,  1980 are  used  as
!     starting values to Newton's Iterations which compute the inverse
!     standard normal.  Therefore no  searches  are necessary for  any
!     parameter.
!
!     For X < -15, the asymptotic expansion for the normal is used  as
!     the starting value in finding the inverse standard normal.
!
!                              Note
!
!      The normal density is proportional to
!      exp( - 0.5 * (( X - MEAN)/SD)**2)
!
  double precision bound
  double precision mean,p,q,sd,x
  integer status,which
  double precision z,pq

  double precision dinvnr,dpmpar
!
!  Check.
!
  status = 0

  if ( which < 1 ) then
    status = - 1
    bound = 1.0d0
    return
  else if ( which > 4 ) then
    status = - 1
    bound = 4.0d0
    return
  end if
!
!     P
!
  if ( which /= 1 ) then

    if ( p <= 0.0 ) then
      status = - 2
      bound = 0.0
      return
    else if ( p > 1.0 ) then
      status = - 2
      bound = 1.0
      return
    end if

  end if
!
!     Q
!
  if ( which /= 1 ) then

    if ( q <= 0.0d0 ) then
      bound = 0.0d0
      status = - 3
      return
    else if ( q > 1.0 ) then
      bound = 1.0d0
      status = -3
      return
    end if

  end if
!
!     P + Q
!
  if ( which /= 1 ) then

    pq = p + q

    if ( abs ( pq - 1.0 ) > 3.0d0 * epsilon ( 1.0d0 ) ) then

      status = 3

      if ( pq < 0.0d0 ) then
        bound = 0.0d0
      else
        bound = 1.0d0
      end if

      return

    end if

  end if

  if ( which /= 4 ) then

    if ( sd <= 0.0d0 ) then
      bound = 0.0d0
      status = -6
      return
    end if

  end if
!
!  Calculate.
!
  if ( which == 1 ) then

    z = ( x - mean ) / sd
    call cumnor ( z, p, q )

  else if ( which == 2 ) then

    z = dinvnr ( p, q )
    x = sd * z + mean

  else if ( which == 3 ) then

    z = dinvnr ( p, q )
    mean = x - sd * z

  else if ( which == 4 ) then

    z = dinvnr ( p, q )
    sd = ( x - mean ) / z

  end if

  return
end
subroutine cdfpoi ( which, p, q, s, xlam, status, bound )
!
!*******************************************************************************
!
!! CDFPOI evaluates the Cumulative Distribution Function of the POIsson distribution
!
!
!  Method:
!
!     Cumulative  distribution function  (P) is  calculated  directly.
!     Computation of other parameters involve a seach for a value that
!     produces  the desired value of  P.   The  search relies  on  the
!     monotonicity of P with the other parameter.
!
!  Reference:
!
!    Abramowitz and Stegun,  
!    Handbook of Mathematical Functions 
!    1966, Formula 26.4.21.
!
!  Parameters:
!
!     WHICH --> Integer indicating which  argument
!               value is to be calculated from the others.
!               Legal range: 1..3
!               iwhich = 1 : Calculate P and Q from S and XLAM
!               iwhich = 2 : Calculate A from P,Q and XLAM
!               iwhich = 3 : Calculate XLAM from P,Q and S
!                    integer WHICH
!
!        P <--> The cumulation from 0 to S of the poisson density.
!               Input range: [0,1].
!                    double precision P
!
!        Q <--> 1-P.
!               Input range: (0, 1].
!               P + Q = 1.0.
!                    double precision Q
!
!        S <--> Upper limit of cumulation of the Poisson.
!               Input range: [0, +infinity).
!               Search range: [0,1E300]
!                    double precision S
!
!     XLAM <--> Mean of the Poisson distribution.
!               Input range: [0, +infinity).
!               Search range: [0,1E300]
!                    double precision XLAM
!
!     STATUS <-- 0 if calculation completed correctly
!               -I if input parameter number I is out of range
!                1 if answer appears to be lower than lowest
!                  search bound
!                2 if answer appears to be higher than greatest
!                  search bound
!                3 if P + Q /= 1
!                    integer STATUS
!
!     BOUND <-- Undefined if STATUS is 0
!
!               Bound exceeded by parameter number I if STATUS
!               is negative.
!
!               Lower search bound if STATUS is 1.
!
!               Upper search bound if STATUS is 2.
!
  double precision, parameter :: tol=1.0d-8
  double precision, parameter :: atol=1.0d-50
  double precision, parameter :: inf=1.0d300
!
  double precision bound,p,q,s,xlam
  integer status,which
  double precision fx,cum,ccum,pq
  logical qhi,qleft,qporq

  double precision dpmpar
!
!     Check arguments
!
  if ( which < 1 ) then
    bound = 1.0d0
    status = - 1
    return
  end if

  if ( which > 3 ) then
    bound = 3.0d0
    status = -1
    return
  end if

  if (which==1) go to 70
!
!     P
!
  if (.not. ((p<0.0d0).or. (p>1.0d0))) go to 60
  if (.not. (p<0.0d0)) go to 40
  bound = 0.0d0
  go to 50

   40 bound = 1.0d0
   50 status = -2
  return

   60 continue
   70 if ( which == 1 ) go to 110
!
!     Q
!
  if (.not. ((q<=0.0d0).or. (q>1.0d0))) go to 100
  if (.not. (q<=0.0d0)) go to 80
  bound = 0.0d0
  go to 90

   80 bound = 1.0d0
   90 status = -3
  return

  100 continue
  110 if (which==2) go to 130
!
!     S
!
  if (.not. (s<0.0d0)) go to 120
  bound = 0.0d0
  status = -4
  return

  120 continue
  130 if (which==3) go to 150
!
!     XLAM
!
  if (.not. (xlam<0.0d0)) go to 140
  bound = 0.0d0
  status = -5
  return

  140 continue
  150 if (which==1) go to 190
!
!     P + Q
!
  pq = p + q
  if (.not. (abs(((pq)-0.5d0)-0.5d0)> (3.0d0* epsilon ( 1.0d0 ) ))) then
    go to 180
  end if

  if (.not. (pq<0.0d0)) go to 160
  bound = 0.0d0
  go to 170

  160 bound = 1.0d0
  170 status = 3
  return

  180 continue
  190 continue

  if ( which /= 1 ) then
    qporq = p <= q
  end if
!
!     Calculate ANSWERS
!
  if ( which == 1 ) then

    call cumpoi(s,xlam,p,q)
    status = 0

  else if ( which == 2 ) then

    s = 5.0d0
    call dstinv(0.0d0,inf,0.5d0,0.5d0,5.0d0,atol,tol)
    status = 0
    fx = 0.0
    call dinvr ( status, s, fx, qleft, qhi )

  200   continue

    if ( status==1 ) then

      call cumpoi(s,xlam,cum,ccum)

      if ( qporq ) then
        fx = cum - p
      else
        fx = ccum - q
      end if

      call dinvr ( status, s, fx, qleft, qhi )
      go to 200

    else if ( status == -1 ) then

      if ( qleft ) then
        status = 1
        bound = 0.0d0
      else
        status = 2
        bound = inf
      end if

    end if

  else if ( which == 3 ) then

      xlam = 5.0d0
      call dstinv(0.0d0,inf,0.5d0,0.5d0,5.0d0,atol,tol)
      status = 0
      call dinvr ( status, xlam, fx, qleft, qhi )

  270     continue

      if ( status == 1 ) then

        call cumpoi(s,xlam,cum,ccum)

        if ( qporq ) then
          fx = cum - p
        else
          fx = ccum - q
        end if

        call dinvr ( status, xlam, fx, qleft, qhi )
        go to 270

      else if ( status == -1 ) then

        if ( qleft ) then
          status = 1
          bound = 0.0d0
        else
          status = 2
          bound = inf
        end if

      end if

  end if

  return
end
subroutine cdft ( which, p, q, t, df, status, bound )
!
!*******************************************************************************
!
!! CDFT evaluates the cumulative distribution functionof the T distribution.
!
!
!  Method:
!
!     Computation of other parameters involve a seach for a value that
!     produces  the desired  value  of P.   The search relies  on  the
!     monotonicity of P with the other parameter.
!
!  Reference:
!
!    Abramowitz and Stegun,  
!    Handbook of Mathematical Functions 
!    1966, Formula 26.5.27.
!
!  Parameters:
!
!     WHICH --> Integer indicating which  argument
!               values is to be calculated from the others.
!               Legal range: 1..3
!               iwhich = 1 : Calculate P and Q from T and DF
!               iwhich = 2 : Calculate T from P,Q and DF
!               iwhich = 3 : Calculate DF from P,Q and T
!                    integer WHICH
!
!        P <--> The integral from -infinity to t of the t-density.
!              Input range: (0,1].
!                    double precision P
!
!     Q <--> 1-P.
!            Input range: (0, 1].
!            P + Q = 1.0.
!                    double precision Q
!
!        T <--> Upper limit of integration of the t-density.
!               Input range: ( -infinity, +infinity).
!               Search range: [ -1E300, 1E300 ]
!                    double precision T
!
!        DF <--> Degrees of freedom of the t-distribution.
!                Input range: (0 , +infinity).
!                Search range: [1e-300, 1E10]
!                    double precision DF
!
!     STATUS <-- 0 if calculation completed correctly
!               -I if input parameter number I is out of range
!                1 if answer appears to be lower than lowest
!                  search bound
!                2 if answer appears to be higher than greatest
!                  search bound
!                3 if P + Q /= 1
!                    integer STATUS
!
!     BOUND <-- Undefined if STATUS is 0
!
!               Bound exceeded by parameter number I if STATUS
!               is negative.
!
!               Lower search bound if STATUS is 1.
!
!               Upper search bound if STATUS is 2.
!
  double precision, parameter :: tol=1.0d-8
  double precision, parameter :: atol=1.0d-50
  double precision, parameter :: inf=1.0d300
  double precision, parameter :: maxdf=1.0d10
!
  double precision bound
  double precision df
  double precision p
  double precision q
  integer status
  double precision t
  integer which

  double precision fx,cum,ccum,pq
  logical qhi,qleft,qporq
  double precision dpmpar,dt1
!
!     Check arguments
!
  if ( which < 1 ) then
    bound = 1.0d0
    status = - 1
    return
  end if

  if ( which > 3 ) then
    bound = 3.0d0
    status = -1
    return
  end if

  if ( which /= 1 ) then

    if ( p <= 0.0d0 ) then
      bound = 0.0d0
      status = -2
      return
    else if ( p > 1.0 ) then
      bound = 1.0d0
      status = -2
      return
    end if

  end if
!
!  Q
!
  if ( which /= 1 ) then

    if (.not. ((q<=0.0d0).or. (q>1.0d0))) go to 100
    if (.not. (q<=0.0d0)) go to 80
    bound = 0.0d0
    go to 90

   80   bound = 1.0d0
   90   status = -3
    return

  100   continue

  end if
!
!  DF
!
  if ( which /= 3 ) then

    if ( df <= 0.0d0 ) then
      bound = 0.0d0
      status = -5
      return
    end if

  end if

  130 continue
!
!     P + Q
!
  if ( which /= 1 ) then

    pq = p + q

    if ( abs((pq-0.5d0)-0.5d0) > 3.0d0 * epsilon ( 1.0d0 ) ) then

      if ( pq < 0.0d0 ) then
        bound = 0.0d0
      else
        bound = 1.0d0
      end if

      status = 3
      return

    end if

  end if

  if ( which /= 1 ) then
    qporq = p <= q
  end if
!
!     Calculate ANSWERS
!
  if ( which == 1 ) then

    call cumt ( t, df, p, q )
    status = 0

  else if ( which == 2 ) then

    t = dt1 ( p, q, df )
    call dstinv(-inf,inf,0.5d0,0.5d0,5.0d0,atol,tol)
    status = 0
    fx = 0.0
    call dinvr ( status, t, fx, qleft, qhi )

  180   continue

    if ( status == 1 ) then

      call cumt(t,df,cum,ccum)

      if ( qporq ) then
        fx = cum - p
      else
        fx = ccum - q
      end if

      call dinvr ( status, t, fx, qleft, qhi )
      go to 180

    else if ( status == -1 ) then

      if ( qleft )then
        status = 1
        bound = - inf
      else
        status = 2
        bound = inf
      end if

    end if

  else if ( which == 3 ) then

    df = 5.0d0
    call dstinv(0.0d0,maxdf,0.5d0,0.5d0,5.0d0,atol,tol)
    status = 0
    call dinvr ( status, df, fx, qleft, qhi )

  250   continue

    if ( status == 1 ) then

      call cumt ( t, df, cum, ccum )

      if ( qporq ) then
        fx = cum - p
      else
        fx = ccum - q
      end if

      call dinvr ( status, df, fx, qleft, qhi )

      go to 250

    else if ( status == -1 ) then

      if ( qleft ) then
        status = 1
        bound = 0.0
      else
        status = 2
        bound = maxdf
      end if

    end if

  end if

  return
end
subroutine cumbet ( x, y, a, b, cum, ccum )
!
!*******************************************************************************
!
!! CUMBET evaluates the cumulative incomplete beta distribution.
!
!
!  Discussion:
!
!    Calculates the cdf to X of the incomplete beta distribution
!    with parameters a and b.  This is the integral from 0 to x
!    of (1/B(a,b))*f(t)) where f(t) = t**(a-1) * (1-t)**(b-1)
!
!  Reference:
!
!    Didonato, Armido R. and Morris, Alfred H. Jr. 
!    (1992) 
!    Algorithm 708 
!    Significant Digit Computation of the Incomplete Beta Function Ratios. 
!    ACM ToMS, Vol.18, No. 3, Sept. 1992, 360-373.
!
!  Parameters:
!
!
!     X --> Upper limit of integration.
!                                        X is double precision
!
!     Y --> 1 - X.
!                                        Y is double precision
!
!     A --> First parameter of the beta distribution.
!                                        A is double precision
!
!     B --> Second parameter of the beta distribution.
!                                        B is double precision
!
!     CUM <-- Cumulative incomplete beta distribution.
!                                        CUM is double precision
!
!     CCUM <-- Complement of Cumulative incomplete beta distribution.
!                                        CCUM is double precision
!
  double precision a
  double precision b
  double precision ccum
  double precision cum
  integer ierr
  double precision x
  double precision y
!
  if ( x <= 0.0d0 ) then

    cum = 0.0d0
    ccum = 1.0d0

  else if ( y <= 0.0d0 ) then

    cum = 1.0d0
    ccum = 0.0d0

  else

    call beta_ratio ( a, b, x, y, cum, ccum, ierr )

  end if

  return
end
subroutine cumbin ( s, xn, pr, ompr, cum, ccum )
!
!*******************************************************************************
!
!! CUMBIN evaluates the cumulative binomial distribution.
!
!
!  Discussion:
!
!    Returns the probability of 0 to S successes in XN binomial
!    trials, each of which has a probability of success, PR.
!
!  Reference:
!
!    Abramowitz and Stegun,  
!    Handbook of Mathematical Functions 
!    1966, Formula 26.5.24.
!
!  Parameters:
!
!     S --> The upper limit of cumulation of the binomial distribution.
!                                                  S is double precision
!
!     XN --> The number of binomial trials.
!                                                  XN is DOUBLE PRECISIO
!
!     PR --> The probability of success in each binomial trial.
!                                                  PR is DOUBLE PRECIS
!
!     OMPR --> 1 - PR
!                                                  OMPR is DOUBLE PRECIS
!
!     CUM <-- Cumulative binomial distribution.
!                                                  CUM is DOUBLE PRECISI
!
!     CCUM <-- Complement of Cumulative binomial distribution.
!                                                  CCUM is DOUBLE PRECIS
  double precision ccum
  double precision cum
  double precision ompr
  double precision pr
  double precision s
  double precision xn
!
  if ( s < xn ) then

    call cumbet ( pr, ompr, s + 1.0d0, xn - s, ccum, cum )

  else

    cum = 1.0d0
    ccum = 0.0d0

  end if

  return
end
subroutine cumchi ( x, df, cum, ccum )
!
!*******************************************************************************
!
!! CUMCHI evaluates the cumulative chi-square distribution.
!
!
!  Parameters:
!
!     X       --> Upper limit of integration of the
!                 chi-square distribution.
!                                                 X is double precision
!
!     DF      --> Degrees of freedom of the
!                 chi-square distribution.
!                                                 DF is double precision
!
!     CUM <-- Cumulative chi-square distribution.
!                                                 CUM is DOUBLE PRECISIO
!
!     CCUM <-- Complement of Cumulative chi-square distribution.
!                                                 CCUM is DOUBLE PRECISI
!
  double precision a
  double precision ccum
  double precision cum
  double precision df
  double precision x
  double precision xx
!
  a = df * 0.5d0
  xx = x * 0.5d0

  call cumgam ( xx, a, cum, ccum )

  return
end
subroutine cumchn ( x, df, pnonc, cum, ccum )
!
!*******************************************************************************
!
!! CUMCHN evaluates the cumulative non-central chi-square distribution.
!
!
!  Function:
!
!    Calculates the cumulative non-central chi-square
!    distribution, i.e., the probability that a random variable
!    which follows the non-central chi-square distribution, with
!    non-centrality parameter PNONC and continuous  degrees of
!    freedom DF, is less than or equal to X.
!
!  Reference:
!
!    Abramowitz and Stegun,  
!    Handbook of Mathematical Functions 
!    1966, Formula 26.4.25.
!
!  Parameters:
!
!     X       --> Upper limit of integration of the non-central
!                 chi-square distribution.
!                                                 X is double precision
!
!     DF      --> Degrees of freedom of the non-central
!                 chi-square distribution.
!                                                 DF is double precision
!
!     PNONC   --> Non-centrality parameter of the non-central
!                 chi-square distribution.
!                                                 PNONC is DOUBLE PRECIS
!
!     CUM <-- Cumulative non-central chi-square distribution.
!                                                 CUM is DOUBLE PRECISIO
!
!     CCUM <-- Complement of Cumulative non-central chi-square distribut
!                                                 CCUM is DOUBLE PRECISI
!
!                              Variables
!
!
!     EPS     --- Convergence criterion.  The sum stops when a
!                 term is less than EPS*SUM.
!                                                 EPS is DOUBLE PRECISIO
!
!     NTIRED  --- Maximum number of terms to be evaluated
!                 in each sum.
!                                                 NTIRED is integer
!
!     QCONV   --- .TRUE. if convergence achieved -
!                 i.e., program did not stop on NTIRED criterion.
!                                                 QCONV is LOGICAL
!
  double precision adj
  double precision ccum
  double precision centaj
  double precision centwt
  double precision chid2
  double precision cum
  double precision df
  double precision dfd2
  double precision dg
  double precision, parameter :: eps = 0.00001
  double precision gamma_log
  integer i
  integer icent
  integer iterb
  integer iterf
  double precision lcntaj
  double precision lcntwt      
  double precision lfact
  integer, parameter :: ntired = 1000
  double precision pcent
  double precision pnonc
  double precision pterm
  logical qsmall
  logical qtired
  double precision sum
  double precision sumadj
  double precision term
  double precision wt
  double precision x
  double precision xnonc
  double precision xx
!
  qtired(i) = i > ntired
  qsmall(xx) = sum < 1.0d-20 .or. xx < eps*sum
  dg(i) = df + 2.0d0 * dble(i)

  if ( x <= 0.0d0 ) then
    cum = 0.0d0
    ccum = 1.0d0
    return
  end if
!
!  When non-centrality parameter is (essentially) zero,
!  use cumulative chi-square distribution
!
  if ( pnonc <= 1.0d-10 ) then
    call cumchi ( x, df, cum, ccum )
    return
  end if

  xnonc = pnonc / 2.0d0
!
!  The following code calculates the weight, chi-square, and
!  adjustment term for the central term in the infinite series.
!  The central term is the one in which the poisson weight is
!  greatest.  The adjustment term is the amount that must
!  be subtracted from the chi-square to move up two degrees
!  of freedom.
!
  icent = int ( xnonc )
  if ( icent == 0 ) then
    icent = 1
  end if

  chid2 = x / 2.0d0
!
!  Calculate central weight term
!
  lfact = gamma_log ( dble ( icent + 1 ))
  lcntwt = - xnonc + icent * log ( xnonc ) - lfact
  centwt = exp ( lcntwt )
!
!  Calculate central chi-square.
!
  call cumchi ( x, dg(icent), pcent, ccum )
!
!  Calculate central adjustment term
!
  dfd2 = dg(icent) / 2.0d0
  lfact = gamma_log ( 1.0d0 + dfd2 )
  lcntaj = dfd2 * log ( chid2 ) - chid2 - lfact
  centaj = exp ( lcntaj )
  sum = centwt * pcent
!
!  Sum backwards from the central term towards zero.
!  Quit whenever either
!  (1) the zero term is reached, or
!  (2) the term gets small relative to the sum, or
!  (3) More than NTIRED terms are totaled.
!
  iterb = 0
  sumadj = 0.0d0
  adj = centaj
  wt = centwt
  i = icent
  term = 0.0
  go to 40

   30 if (qtired(iterb) .or. qsmall(term) .or. i==0) go to 50
   40 dfd2 = dg(i)/2.0d0
!
!  Adjust chi-square for two fewer degrees of freedom.
!  The adjusted value ends up in PTERM.
!
  adj = adj * dfd2 / chid2
  sumadj = sumadj + adj
  pterm = pcent + sumadj
!
!  Adjust Poisson weight for J decreased by one
!
  wt = wt * ( i / xnonc )
  term = wt * pterm
  sum = sum + term
  i = i - 1
  iterb = iterb + 1
  go to 30

   50 iterf = 0
!
!  Now sum forward from the central term towards infinity.
!  Quit when either
!    (1) the term gets small relative to the sum, or
!    (2) More than NTIRED terms are totaled.
!
  sumadj = centaj
  adj = centaj
  wt = centwt
  i = icent

  go to 70

   60 if (qtired(iterf) .or. qsmall(term)) go to 80
!
!  Update weights for next higher J
!
   70 continue

  wt = wt * ( xnonc / ( i + 1 ) )
!
!  Calculate PTERM and add term to sum
!
  pterm = pcent - sumadj
  term = wt * pterm
  sum = sum + term
!
!  Update adjustment term for DF for next iteration
!
  i = i + 1
  dfd2 = dg(i) / 2.0d0
  adj = adj * chid2 / dfd2
  sumadj = sumadj + adj
  iterf = iterf + 1
  go to 60

   80 continue

  cum = sum
  ccum = 0.5d0 + ( 0.5d0 - cum )

  return
end
subroutine cumf ( f, dfn, dfd, cum, ccum )
!
!*******************************************************************************
!
!! CUMF evaluates the cumulative F distribution.
!
!
!  Discussion:
!
!    CUMF computes the integral from 0 to F of the F density with DFN
!    numerator and DFD denominator degrees of freedom.
!
!  Reference:
!
!    Abramowitz and Stegun,  
!    Handbook of Mathematical Functions 
!    1966, Formula 26.5.28.
!
!  Parameters:
!
!     F --> Upper limit of integration of the f-density.
!                                                  F is double precision
!
!     DFN --> Degrees of freedom of the numerator sum of squares.
!                                                  DFN is DOUBLE PRECISI
!
!     DFD --> Degrees of freedom of the denominator sum of squares.
!                                                  DFD is DOUBLE PRECISI
!
!     CUM <-- Cumulative f distribution.
!                                                  CUM is DOUBLE PRECISI
!
!     CCUM <-- Complement of Cumulative f distribution.
!                                                  CCUM is DOUBLE PRECIS
!
  double precision ccum
  double precision cum
  double precision dfd
  double precision dfn
  double precision dsum
  double precision f
  integer ierr
  double precision prod
  double precision xx
  double precision yy
!
  if ( f <= 0.0d0 ) then
    cum = 0.0d0
    ccum = 1.0d0
    return
  end if

  prod = dfn * f
!
!  XX is such that the incomplete beta with parameters
!  DFD/2 and DFN/2 evaluated at XX is 1 - CUM or CCUM
!
!  YY is 1 - XX
!
!  Calculate the smaller of XX and YY accurately
!
  dsum = dfd + prod
  xx = dfd / dsum

  if ( xx > 0.5 ) then
    yy = prod / dsum
    xx = 1.0 - yy
  else
    yy = 1.0 - xx
  end if

  call beta_ratio ( 0.5d0*dfd, 0.5d0*dfn, xx, yy, ccum, cum, ierr )

  return
end
subroutine cumfnc ( f, dfn, dfd, pnonc, cum, ccum )
!
!*******************************************************************************
!
!! CUMFNC evaluates the cumulative noncentral F distribution.
!
!
!  Discussion:
!
!    Computes noncentral F distribution with DFN and DFD
!    degrees of freedom and noncentrality parameter PNONC.
!
!  Method:
!
!     The series is calculated backward and forward from J = LAMBDA/2
!     (this is the term with the largest Poisson weight) until
!     the convergence criterion is met.
!
!     The sum continues until a succeeding term is less than EPS
!     times the sum (or the sum is less than 1.0e-20).  EPS is
!     set to 1.0e-4 in a data statement which can be changed.
!
!  Reference:
!
!    Abramowitz and Stegun,  
!    Handbook of Mathematical Functions 
!    1966, Formula 26.5.16, 26.6.17, 26.6.18, 26.6.20.
!
!  Parameters:
!
!     X --> Upper limit of integration of noncentral F in equation.
!
!     DFN --> Degrees of freedom of numerator
!
!     DFD -->  Degrees of freedom of denominator
!
!     PNONC --> Noncentrality parameter.
!
!     CUM <-- Cumulative noncentral F distribution
!
!     CCUM <-- Complement of cummulative
!
  double precision adn
  double precision aup
  double precision b
  double precision betdn
  double precision betup
  double precision centwt
  double precision ccum
  double precision cum
  double precision dfd
  double precision dfn
  double precision dnterm
  double precision dsum
  double precision dummy
  double precision, parameter :: eps = 0.0001
  double precision f
  double precision gamma_log
  integer i
  integer icent
  integer ierr
  double precision pnonc
  double precision prod
  logical qsmall
  double precision sum
  double precision upterm
  double precision x
  double precision xmult
  double precision xnonc
  double precision xx
  double precision yy
!
  qsmall(x) = sum < 1.0d-20 .or. x < eps*sum
!
  if ( f <= 0.0d0 ) then
    cum = 0.0d0
    ccum = 1.0d0
    return
  end if
!
!  Handle case in which the non-centrality parameter is essentially zero.
!
  if ( pnonc < 1.0d-10 ) then
    call cumf ( f, dfn, dfd, cum, ccum )
    return
  end if

  xnonc = pnonc / 2.0d0
!
!  Calculate the central term of the poisson weighting factor.
!
  icent = xnonc

  if ( icent == 0 ) then
    icent = 1
  end if
!
!  Compute central weight term
!
  centwt = exp(-xnonc+icent*log(xnonc)-gamma_log(dble(icent+1)))
!
!  Compute central incomplete beta term.
!  Ensure that minimum of arg to beta and 1 - arg is computed accurately.
!
  prod = dfn * f
  dsum = dfd + prod
  yy = dfd / dsum

  if ( yy > 0.5 ) then
    xx = prod / dsum
    yy = 1.0 - xx
  else
    xx = 1.0 - yy
  end if

  call beta_ratio ( 0.5d0*dfn+dble(icent), 0.5d0*dfd, xx, yy, betdn, &
    dummy, ierr )

  adn = dfn / 2.0d0 + dble ( icent )
  aup = adn
  b = dfd / 2.0d0
  betup = betdn
  sum = centwt * betdn
!
!  Now sum terms backward from icent until convergence or all done
!
  xmult = centwt
  i = icent
  dnterm = exp ( gamma_log ( adn + b )- gamma_log ( adn + 1.0d0 ) &
    - gamma_log ( b ) + adn * log ( xx ) + b * log ( yy ) )

   30 continue

  if ( qsmall(xmult*betdn) .or. i <= 0 ) go to 40
  xmult = xmult * (i/xnonc)
  i = i - 1
  adn = adn - 1
  dnterm = ( adn + 1 ) / (( adn + b ) * xx ) * dnterm
  betdn = betdn + dnterm
  sum = sum + xmult * betdn
  go to 30

40 continue

  i = icent + 1
!
!  Now sum forwards until convergence
!
  xmult = centwt
  if (( aup - 1 + b ) == 0 ) then

    upterm = exp ( - gamma_log ( aup ) - gamma_log ( b ) &
      + ( aup - 1 ) * log ( xx ) + b * log ( yy ) )

  else

    upterm = exp ( gamma_log ( aup-1+b ) - gamma_log ( aup ) &
      - gamma_log ( b ) + ( aup - 1 ) * log ( xx ) + b * log ( yy ) )

  end if

  go to 60

   50 continue

  if ( qsmall(xmult*betup)) go to 70

   60 continue

  xmult = xmult * ( xnonc / i )
  i = i + 1
  aup = aup + 1
  upterm = ( aup + b - 2.0d0 ) * xx / ( aup - 1 ) * upterm
  betup = betup - upterm
  sum = sum + xmult * betup
  go to 50

   70 continue

  cum = sum

  ccum = 0.5d0 + ( 0.5d0 - cum )

  return
end
subroutine cumgam ( x, a, cum, ccum )
!
!*******************************************************************************
!
!! CUMGAM evaluates the cumulative incomplete gamma distribution.
!
!
!  Discussion:
!
!    Computes the cumulative of the incomplete gamma distribution, i.e., 
!    the integral from 0 to X of
!          (1/GAM(A))*EXP(-T)*T**(A-1) DT
!    where GAM(A) is the complete gamma function of A, i.e.,
!          GAM(A) = integral from 0 to infinity of
!                    EXP(-T)*T**(A-1) DT
!
!
!  Parameters:
!
!     X --> The upper limit of integration of the incomplete gamma.
!                                                X is double precision
!
!     A --> The shape parameter of the incomplete gamma.
!                                                A is double precision
!
!     CUM <-- Cumulative incomplete gamma distribution.
!                                        CUM is double precision
!
!     CCUM <-- Complement of Cumulative incomplete gamma distribution.
!                                                CCUM is DOUBLE PRECISIO
!
  double precision a
  double precision ccum
  double precision cum
  double precision x
!
  if ( x <= 0.0 ) then

    cum = 0.0d0
    ccum = 1.0d0

  else

    call gamma_inc ( a, x, cum, ccum, 0 )

  end if

  return
end
subroutine cumnbn ( s, xn, pr, ompr, cum, ccum )
!
!*******************************************************************************
!
!! CUMNBN evaluates the cumulative negative binomial distribution.
!
!
!  Function:
!
!     Returns the probability that it there will be S or fewer failures
!     before there are XN successes, with each binomial trial having
!     a probability of success PR.
!
!     Prob(# failures = S | XN successes, PR)  =
!                        ( XN + S - 1 )
!                        (            ) * PR^XN * (1-PR)^S
!                        (      S     )
!
!
!  Reference:
!
!    Abramowitz and Stegun,  
!    Handbook of Mathematical Functions 
!    1966, Formula 26.5.26.
!
!  Parameters:
!
!
!     S --> The number of failures
!                                                  S is double precision
!
!     XN --> The number of successes
!                                                  XN is DOUBLE PRECISIO
!
!     PR --> The probability of success in each binomial trial.
!                                                  PR is DOUBLE PRECISIO
!
!     OMPR --> 1 - PR
!                                                  OMPR is DOUBLE PRECIS
!
!     CUM <-- Cumulative negative binomial distribution.
!                                                  CUM is DOUBLE PRECISI
!
!     CCUM <-- Complement of Cumulative negative binomial distribution.
!                                                  CCUM is DOUBLE PRECIS
!
  double precision ccum
  double precision cum
  double precision ompr
  double precision pr
  double precision s
  double precision xn
!
  call cumbet ( pr, ompr, xn, s+1.d0, cum, ccum )

  return
end
subroutine cumnor ( arg, result, ccum )
!
!*******************************************************************************
!
!! CUMNOR computes the cumulative normal distribution.
!
!
!     the integral from -infinity to x of
!          (1/sqrt(2*pi)) exp(-u*u/2) du
!
!  Reference:
!
!    W D Cody, 
!    "ALGORITHM 715: SPECFUN - A Portable FORTRAN Package of Special 
!    Function Routines and Test Drivers"
!    ACM Transactions on Mathematical Software,
!    Volume 19, 1993, pages 22-32.
!
!  Parameters:
!
!     ARG --> Upper limit of integration.
!                                        X is double precision
!
!     RESULT <-- Cumulative normal distribution.
!                                        RESULT is double precision
!
!     CCUM <-- Complement of Cumulative normal distribution.
!                                        CCUM is double precision
!
!
! Original Comments:
!
!
! This function evaluates the normal distribution function:
!
!                              / x
!                     1       |       -t*t/2
!          P(x) = ----------- |      e       dt
!                 sqrt(2 pi)  |
!                             /-oo
!
!   The main computation evaluates near-minimax approximations
!   derived from those in "Rational Chebyshev approximations for
!   the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
!   This transportable program uses rational functions that
!   theoretically approximate the normal distribution function to
!   at least 18 significant decimal digits.  The accuracy achieved
!   depends on the arithmetic system, the compiler, the intrinsic
!   functions, and proper selection of the machine-dependent
!   constants.
!
!  Explanation of machine-dependent constants.
!
!   MIN   = smallest machine representable number.
!
!   EPS   = argument below which anorm(x) may be represented by
!           0.5  and above which  x*x  will not underflow.
!           A conservative value is the largest machine number X
!           such that   1.0 + X = 1.0   to machine precision.
!
!  Error returns
!
!  The program returns  ANORM = 0     for  ARG .LE. XLOW.
!
!  Author: 
!
!    W. J. Cody
!    Mathematics and Computer Science Division
!    Argonne National Laboratory
!    Argonne, IL 60439
!
!  Latest modification: March 15, 1992
!
  double precision, parameter, dimension ( 5 ) :: a = (/ &
    2.2352520354606839287d00, &
    1.6102823106855587881d02, &
    1.0676894854603709582d03, &
    1.8154981253343561249d04, &
    6.5682337918207449113d-2 /)
  double precision arg
  double precision, parameter, dimension ( 4 ) :: b = (/ &
    4.7202581904688241870d01, &
    9.7609855173777669322d02, &
    1.0260932208618978205d04, &
    4.5507789335026729956d04 /)
  double precision, parameter, dimension ( 9 ) :: c = (/ &
    3.9894151208813466764d-1, &
    8.8831497943883759412d00, &
    9.3506656132177855979d01, &
    5.9727027639480026226d02, &
    2.4945375852903726711d03, &
    6.8481904505362823326d03, &
    1.1602651437647350124d04, &
    9.8427148383839780218d03, &
    1.0765576773720192317d-8 /)
  double precision ccum
  double precision, parameter, dimension ( 8 ) :: d = (/ &
    2.2266688044328115691d01, &
    2.3538790178262499861d02, &
    1.5193775994075548050d03, &
    6.4855582982667607550d03, &
    1.8615571640885098091d04, &
    3.4900952721145977266d04, &
    3.8912003286093271411d04, &
    1.9685429676859990727d04 /)
  double precision del
  double precision dpmpar
  double precision eps
  integer i
  double precision min
  double precision, parameter, dimension ( 6 ) :: p = (/ &
    2.1589853405795699d-1, &
    1.274011611602473639d-1, &
    2.2235277870649807d-2, &
    1.421619193227893466d-3, &
    2.9112874951168792d-5, &
    2.307344176494017303d-2 /)
  double precision, parameter, dimension ( 5 ) :: q = (/ &
    1.28426009614491121d00, &
    4.68238212480865118d-1, &
    6.59881378689285515d-2, &
    3.78239633202758244d-3, &
    7.29751555083966205d-5 /)
  double precision result
  double precision, parameter :: root32 = 5.656854248d0
  double precision, parameter :: sixten = 16.0
  double precision temp
  double precision, parameter :: sqrpi = 3.9894228040143267794d-1
  double precision, parameter :: thrsh = 0.66291d0
  double precision x
  double precision xden
  double precision xnum
  double precision y
  double precision xsq
!
!  Machine dependent constants
!
  eps = epsilon ( 1.0d0 ) * 0.5d0
  min = dpmpar(2)

  x = arg
  y = abs ( x )

  if ( y <= thrsh ) then
!
!  Evaluate  anorm  for  |X| <= 0.66291
!
    if ( y > eps ) then
      xsq = x * x
    else
      xsq = 0.0
    end if

    xnum = a(5) * xsq
    xden = xsq
    do i = 1, 3
      xnum = ( xnum + a(i) ) * xsq
      xden = ( xden + b(i) ) * xsq
    end do
    result = x * ( xnum + a(4) ) / ( xden + b(4) )
    temp = result
    result = 0.5 + temp
    ccum = 0.5 - temp
!
!  Evaluate ANORM for 0.66291 <= |X| <= sqrt(32)
!
  else if ( y <= root32 ) then

    xnum = c(9) * y
    xden = y
    do i = 1, 7
      xnum = ( xnum + c(i) ) * y
      xden = ( xden + d(i) ) * y
    end do
    result = ( xnum + c(8) ) / ( xden + d(8) )
    xsq = aint ( y * sixten ) / sixten
    del = ( y - xsq ) * ( y + xsq )
    result = exp(-xsq*xsq*0.5) * exp(-del*0.5) * result
    ccum = 1.0 - result

    if ( x > 0.0 ) then
      call d_swap ( result, ccum )
    end if
!
!  Evaluate  anorm  for |X| > sqrt(32).
!
  else

    result = 0.0
    xsq = 1.0 / ( x * x )
    xnum = p(6) * xsq
    xden = xsq
    do i = 1, 4
      xnum = ( xnum + p(i) ) * xsq
      xden = ( xden + q(i) ) * xsq
    end do

    result = xsq * ( xnum + p(5) ) / ( xden + q(5) )
    result = ( sqrpi - result ) / y
    xsq = aint ( x * sixten ) / sixten
    del = ( x - xsq ) * ( x + xsq )
    result = exp ( - xsq * xsq * 0.5 ) * exp ( - del * 0.5 ) * result
    ccum = 1.0 - result

    if ( x > 0.0 ) then
      call d_swap ( result, ccum )
    end if

  end if

  if ( result < min ) then
    result = 0.0d0
  end if

  if ( ccum < min ) then
    ccum = 0.0d0
  end if

  return
end
subroutine cumpoi ( s, xlam, cum, ccum )
!
!*******************************************************************************
!
!! CUMPOI evaluates the cumulative Poisson distribution.
!
!
!  Discussion:
!
!    CUMPOI returns the  probability of S or fewer events in a Poisson
!    distribution with mean XLAM.
!
!  Reference:
!
!    Abramowitz and Stegun, 
!    Handbook of Mathematical Functions,
!    Formula 26.4.21.
!
!  Parameters:
!
!     S --> Upper limit of cumulation of the Poisson.
!                                                  S is double precision
!
!     XLAM --> Mean of the Poisson distribution.
!                                                  XLAM is DOUBLE PRECIS
!
!     CUM <-- Cumulative poisson distribution.
!                                        CUM is double precision
!
!     CCUM <-- Complement of Cumulative poisson distribution.
!                                                  CCUM is DOUBLE PRECIS
!
  double precision ccum
  double precision chi
  double precision cum
  double precision df
  double precision s
  double precision xlam
!
  df = 2.0d0 * ( s + 1.0d0 )
  chi = 2.0d0 * xlam

  call cumchi ( chi, df, ccum, cum )

  return
end
subroutine cumt ( t, df, cum, ccum )
!
!*******************************************************************************
!
!! CUMT evaluates the cumulative T distribution.
!
!
!  Reference:
!
!    Abramowitz and Stegun, 
!    Handbook of Mathematical Functions,
!    Formula 26.5.27.
!
!  Parameters:
!
!    Input, double precision T, the upper limit of integration of the 
!    t-density.
!
!    Input, double precision DF, the number of degrees of freedom of 
!    the t-distribution.
!
!    Output, double precision CUM, the cumulative t-distribution.
!
!    Output, double precision CCUM, the complement of the cumulative 
!    t-distribution.
!
  double precision a
  double precision ccum
  double precision cum
  double precision df
  double precision oma
  double precision t
  double precision xx
  double precision yy
!
  xx = df / ( df + t**2 )
  yy = t**2 / ( df + t**2 )
 
  call cumbet ( xx, yy, 0.5d0*df, 0.5d0, a, oma )

  if ( t <= 0.0d0 ) then
    cum = 0.5d0 * a
    ccum = oma + cum
  else
    ccum = 0.5d0 * a
    cum = oma + ccum
  end if

  return
end
function dbetrm ( a, b )
!
!*******************************************************************************
!
!! DBETRM computes the Sterling remainder for the complete beta function.
!
!
!  Function:
!
!     Log(Beta(A,B)) = Lgamma(A) + Lgamma(B) - Lgamma(A+B)
!     where Lgamma is the log of the (complete) gamma function
!
!     Let ZZ be approximation obtained if each log gamma is approximated
!     by Sterling's formula, i.e.,
!     Sterling(Z) = LOG( SQRT( 2*PI ) ) + ( Z-0.5 ) * LOG( Z ) - Z
!
!     Returns Log(Beta(A,B)) - ZZ
!
!  Parameters:
!
!     A --> One argument of the Beta
!                    double precision A
!
!     B --> The other argument of the Beta
!                    double precision B
!
  double precision a
  double precision b
  double precision dbetrm
  double precision dstrem
!
!  Try to sum from smallest to largest
!
  dbetrm = - dstrem ( a + b )
  dbetrm = dbetrm + dstrem ( max ( a, b ) )
  dbetrm = dbetrm + dstrem ( min ( a, b ) )

  return
end
function eval_pol ( a, n, x )
!
!*******************************************************************************
!
!! EVAL_POL evaluates a polynomial at X.
!
!
!  Discussion:
!
!    EVAL_POL = A(0) + A(1)*X + ... + A(N)*X**N
!
!  Modified:
!
!    15 December 1999
!
!  Parameters:
!
!    Input, double precision A(0:N), coefficients of the polynomial.
!
!    Input, integer N, length of A.
!
!    Input, double precision X, the point at which the polynomial 
!    is to be evaluated.
!
!    Output, double precision EVAL_POL, the value of the polynomial at X.
!
  integer n
!
  double precision a(0:n)
  double precision eval_pol
  integer i
  double precision term
  double precision x
!
  term = a(n)
  do i = n - 1, 0, -1
    term = term * x + a(i)
  end do

  eval_pol = term

  return
end
function dexpm1 ( x )
!
!*******************************************************************************
!
!! DEXPM1 evaluates the function EXP(X) - 1.
!
!
!  Reference:
!
!    A R DiDinato and A H Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Trans. Math.  Software, 
!    Volume 18, 1993, pages 360-373.
!
!  Parameters:
!
!     X --> Argument at which exp(x)-1 desired
!                    double precision X
!
  double precision, parameter :: p1 =   0.914041914819518d-09
  double precision, parameter :: p2 =   0.238082361044469d-01
  double precision, parameter :: q1 = - 0.499999999085958d+00
  double precision, parameter :: q2 =   0.107141568980644d+00
  double precision, parameter :: q3 = - 0.119041179760821d-01
  double precision, parameter :: q4 =   0.595130811860248d-03
!
  double precision bot
  double precision dexpm1
  double precision top
  double precision w
  double precision x
!
  if ( abs ( x ) <= 0.15d0 ) then

    top = ( p2 * x + p1 ) * x + 1.0d0
    bot = ((( q4 * x + q3 ) * x + q2 ) * x + q1 ) * x + 1.0d0
    dexpm1 = x * ( top / bot )

  else

    w = exp ( x )

    if ( x <= 0.0d0 ) then
      dexpm1 = ( w - 0.5d0 ) - 0.5d0
    else
      dexpm1 = w * ( 0.5d0 + ( 0.5d0 - 1.0d0 / w ))
    end if

  end if

  return
end
function dinvnr ( p, q )
!
!*******************************************************************************
!
!! DINVNR computes the inverse of the normal distribution.
!
!
!  Function:
!
!    Returns X  such that CUMNOR(X)  =   P,  i.e., the  integral from -
!    infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P
!
!  Parameters:
!
!     P --> The probability whose normal deviate is sought.
!                    P is double precision
!
!     Q --> 1-P
!                    P is double precision
!
!  Method:
!
!     The  rational   function   on  page 95    of Kennedy  and  Gentle,
!     Statistical Computing, Marcel Dekker, NY , 1980 is used as a start
!     value for the Newton method of finding roots.
!
  integer, parameter :: maxit=100
  double precision, parameter :: eps=1.0d-13
  double precision, parameter :: r2pi = 0.3989422804014326d0
!
  double precision ccum
  double precision cum
  double precision dinvnr
  double precision dx
  integer i
  double precision p
  double precision pp
  double precision q
  double precision strtx
  double precision stvaln
  double precision x
  double precision xcur
!
  pp = min ( p, q )
  strtx = stvaln ( pp )
  xcur = strtx
!
!  Newton iterations.
!
  do i = 1, maxit

    call cumnor ( xcur, cum, ccum )
    dx = ( cum - pp ) / ( r2pi * exp ( - 0.5 * xcur**2 ) )
    xcur = xcur - dx

    if ( abs ( dx / xcur ) < eps ) then
      if ( p <= q ) then
        dinvnr = xcur
      else
        dinvnr = - xcur
      end if
      return
    end if

  end do

  if ( p <= q ) then
    dinvnr = strtx
  else
    dinvnr = - strtx
  end if

  return
end
subroutine d_swap ( x, y )
!
!*******************************************************************************
!
!! D_SWAP swaps two double precision values.
!
!
!  Modified:
!
!    01 May 2000
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input/output, double precision X, Y.  On output, the values of X and
!    Y have been interchanged.
!
  double precision x
  double precision y
  double precision z
!
  z = x
  x = y
  y = z

  return
end
subroutine dinvr ( status, x, fx, qleft, qhi )
!
!*******************************************************************************
!
!! DINVR bounds the zero of the function and invokes DZROR.
!
!
!  Function:
!
!    Bounds the    function  and  invokes  ZROR   to perform the   zero
!    finding.  STINVR  must  have   been  called  before this   routine
!    in order to set its parameters.
!
!  Parameters:
!
!     STATUS <--> At the beginning of a zero finding problem, STATUS
!                 should be set to 0 and INVR invoked.  (The value
!                 of parameters other than X will be ignored on this cal
!
!                 When INVR needs the function evaluated, it will set
!                 STATUS to 1 and return.  The value of the function
!                 should be set in FX and INVR again called without
!                 changing any of its other parameters.
!
!                 When INVR has finished without error, it will return
!                 with STATUS 0.  In that case X is approximately a root
!                 of F(X).
!
!                 If INVR cannot bound the function, it returns status
!                 -1 and sets QLEFT and QHI.
!                         integer STATUS
!
!     X <-- The value of X at which F(X) is to be evaluated.
!                         double precision X
!
!     FX --> The value of F(X) calculated when INVR returns with
!            STATUS = 1.
!                         double precision FX
!
!     QLEFT <-- Defined only if QMFINV returns .FALSE.  In that
!          case it is .TRUE. If the stepping search terminated
!          unsucessfully at SMALL.  If it is .FALSE. the search
!          terminated unsucessfully at BIG.
!                    QLEFT is LOGICAL
!
!     QHI <-- Defined only if QMFINV returns .FALSE.  In that
!          case it is .TRUE. if F(X) .GT. Y at the termination
!          of the search and .FALSE. if F(X) .LT. Y at the
!          termination of the search.
!                    QHI is LOGICAL
!
  double precision absstp
  double precision abstol
  double precision big
  double precision fbig
  double precision fsmall
  double precision fx
  integer i99999
  logical qbdd
  logical qhi
  logical qleft
  logical qup
  logical qxmon
  double precision relstp
  double precision reltol
  integer status
  double precision x
  double precision zabsst
  double precision zabsto
  double precision zbig
  double precision zrelst
  double precision zrelto
  double precision zsmall
  double precision zstpmu
  double precision zz

  double precision small,step,stpmul,xhi,xlb,xlo,xsave,xub,yy,zx,zy
  logical qcond,qdum1,qdum2,qincr,qlim
!
  save
!
  qxmon(zx,zy,zz) = zx <= zy .and. zy <= zz
!
  if (status>0) go to 310

  qcond = .not. qxmon(small,x,big)

  if ( qcond ) then
    write ( *, * ) ' '
    write ( *, * ) ' small, x, big not monotone in invr'
    stop
  end if

  xsave = x
!
!  See that SMALL and BIG bound the zero and set QINCR
!
  x = small
!
!  GET-function-VALUE
!
  assign 10 to i99999
  go to 300

   10 fsmall = fx
  x = big
!
!  GET-function-VALUE
!
  assign 20 to i99999
  go to 300

   20 fbig = fx
  qincr = fbig > fsmall

  if ( qincr ) then

    if ( fsmall > 0.0d0 ) then
      status = -1
      qleft = .true.
      qhi = .true.
      return
    end if

    if ( fbig < 0.0d0 ) then
      status = -1
      qleft = .false.
      qhi = .false.
      return
    end if

  else

    if ( fsmall < 0.0d0 ) then
      status = -1
      qleft = .true.
      qhi = .false.
      return
    end if

    if ( fbig > 0.0d0 ) then
      status = -1
      qleft = .false.
      qhi = .true.
      return
    end if

  end if

  x = xsave
  step = max ( absstp, relstp * abs ( x ) )
!
!  YY = F(X) - Y
!  GET-function-VALUE
!
  assign 90 to i99999
  go to 300

   90 yy = fx

  if ( yy == 0.0d0 ) then
    status = 0
    return
  end if

  100 qup = (qincr .and. (yy<0.0d0)) .or. (.not.qincr .and. (yy>0.0d0))
!
!  HANDLE CASE IN WHICH WE MUST STEP HIGHER
!
  if (.not. qup ) go to 170
  xlb = xsave
  xub = min(xlb+step,big)
  go to 120

  110 if (qcond) go to 150
!
!  YY = F(XUB) - Y
!
  120 x = xub
!
!  GET-function-VALUE
!
  assign 130 to i99999
  go to 300

  130 yy = fx
  qbdd = (qincr .and. (yy>=0.0d0)) .or. (.not.qincr .and. (yy<=0.0d0))
  qlim = xub >= big
  qcond = qbdd .or. qlim

  if (qcond) go to 140
  step = stpmul*step
  xlb = xub
  xub = min(xlb+step,big)

  140 continue

  go to 110

  150 continue

  if ( qlim .and. .not. qbdd ) then
    status = -1
    qleft = .false.
    qhi = .not. qincr
    x = big
    return
  end if

  160 go to 240
!
!  HANDLE CASE IN WHICH WE MUST STEP LOWER
!
  170 xub = xsave
  xlb = max(xub-step,small)
  go to 190

  180 if (qcond) go to 220
!
!  YY = F(XLB) - Y
!
  190 x = xlb
!
!  GET-function-VALUE
!
  assign 200 to i99999
  go to 300

  200 yy = fx
  qbdd = (qincr .and. (yy<=0.0d0)) .or.(.not.qincr .and. (yy>=0.0d0))
  qlim = xlb <= small
  qcond = qbdd .or. qlim
  if (qcond) go to 210
  step = stpmul*step
  xub = xlb
  xlb = max ( xub-step, small )
  210 go to 180

  220 if (.not. (qlim.and..not.qbdd)) go to 230
  status = -1
  qleft = .true.
  qhi = qincr
  x = small
  return

  230 continue
  240 continue

  call dstzr ( xlb, xub, abstol, reltol )
!
!     IF WE REACH HERE, XLB AND XUB BOUND THE ZERO OF F.
!
  status = 0
  go to 260

  250 if (.not. (status==1)) go to 290
  260 continue

  call dzror ( status, x, fx, xlo, xhi, qdum1, qdum2 )
  if (.not. (status==1)) go to 280
!
!     GET-function-VALUE
!
  assign 270 to i99999
  go to 300

  270 continue
  280 go to 250

  290 x = xlo
  status = 0
  return

entry dstinv ( zsmall, zbig, zabsst, zrelst, zstpmu, zabsto, zrelto )
!
!*******************************************************************************
!
!! DSTINV SeT INverse finder - Reverse Communication
!
!
!  Function
!
!     Given a monotone function F finds X
!     such that F(X) = Y.  Uses Reverse communication -- see invr.
!     This routine sets quantities needed by INVR.
!
!     F must be a monotone function, the results of QMFINV are
!     otherwise undefined.  QINCR must be .TRUE. if F is non-
!     decreasing and .FALSE. if F is non-increasing.
!
!     QMFINV will return .TRUE. if and only if F(SMALL) and
!     F(BIG) bracket Y, i. e.,
!          QINCR is .TRUE. and F(SMALL).LE.Y.LE.F(BIG) or
!          QINCR is .FALSE. and F(BIG).LE.Y.LE.F(SMALL)
!
!     if QMFINV returns .TRUE., then the X returned satisfies
!     the following condition.  let
!               TOL(X) = MAX(ABSTOL,RELTOL*ABS(X))
!     then if QINCR is .TRUE.,
!          F(X-TOL(X)) .LE. Y .LE. F(X+TOL(X))
!     and if QINCR is .FALSE.
!          F(X-TOL(X)) .GE. Y .GE. F(X+TOL(X))
!
!  Parameters:
!
!     SMALL --> The left endpoint of the interval to be
!          searched for a solution.
!                    SMALL is double precision
!
!     BIG --> The right endpoint of the interval to be
!          searched for a solution.
!                    BIG is double precision
!
!     ABSSTP, RELSTP --> The initial step size in the search
!          is MAX(ABSSTP,RELSTP*ABS(X)). See algorithm.
!                    ABSSTP is double precision
!                    RELSTP is double precision
!
!     STPMUL --> When a step doesn't bound the zero, the step
!                size is multiplied by STPMUL and another step
!                taken.  A popular value is 2.0
!                    double precision STPMUL
!
!     ABSTOL, RELTOL --> Two numbers that determine the accuracy
!          of the solution.  See function for a precise definition.
!                    ABSTOL is double precision
!                    RELTOL is double precision
!
!
!                              Method
!
!
!     Compares F(X) with Y for the input value of X then uses QINCR
!     to determine whether to step left or right to bound the
!     desired x.  the initial step size is
!          MAX(ABSSTP,RELSTP*ABS(S)) for the input value of X.
!     Iteratively steps right or left until it bounds X.
!     At each step which doesn't bound X, the step size is doubled.
!     The routine is careful never to step beyond SMALL or BIG.  If
!     it hasn't bounded X at SMALL or BIG, QMFINV returns .FALSE.
!     after setting QLEFT and QHI.
!
!     If X is successfully bounded then Algorithm R of the paper
!     'Two Efficient Algorithms with Guaranteed Convergence for
!     Finding a Zero of a Function' by J. C. P. Bus and
!     T. J. Dekker in ACM Transactions on Mathematical
!     Software, Volume 1, No. 4 page 330 (DEC. '75) is employed
!     to find the zero of the function F(X)-Y. This is routine
!     QRZERO.
!
  small = zsmall
  big = zbig
  absstp = zabsst
  relstp = zrelst
  stpmul = zstpmu
  abstol = zabsto
  reltol = zrelto
  return
!
!     TO GET-function-VALUE
!
  300 status = 1
  return

  310 continue
  go to i99999

end
function dlanor ( x )
!
!*******************************************************************************
!
!! DLANOR evaluates the logarithm of the asymptotic Normal.
!
!
!  Function
!
!    Computes the logarithm of the cumulative normal distribution
!    from abs ( x ) to infinity for abs ( x ) >= 5.
!
!  Reference:
!
!    Abramowitz and Stegun,  
!    Handbook of Mathematical Functions 
!    1966, Formula 26.2.12.
!
!  Parameters:
!
!      X --> Value at which cumulative normal to be evaluated
!                     double precision X
!
!
!                              Method
!
!      The relative error at X = 5 is about 0.5E-5.
!
!
!                              Note
!
!      ABS(X) must be >= 5 else there is an error stop.
!
  double precision, parameter :: dlsqpi = 0.91893853320467274177d0
!
  double precision alnrel
  double precision approx
  double precision, save, dimension ( 0:11 ) :: coef = (/ &
    -1.0d0,  3.0d0,  -15.0d0,  105.0d0,  -945.0d0,  10395.0d0, &
    -135135.0d0,  2027025.0d0,  -34459425.0d0,  654729075.0d0, &
    -13749310575d0,  316234143225.0d0 /)
  double precision correc
  double precision eval_pol
  double precision dlanor
  double precision x
  double precision xx
  double precision xx2
!
  xx = abs ( x )

  if ( abs ( x ) < 5.0d0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'DLANOR - Fatal error!'
    write ( *, * ) '  Argument too small.'
  end if

  approx = - dlsqpi - 0.5d0 * x**2 - log ( abs ( x ) )

  xx2 = xx * xx
  correc = eval_pol ( coef, 11, 1.0d0 / xx2 ) / xx2
  correc = alnrel ( correc )

  dlanor = approx + correc

  return
end
function dpmpar ( i )
!
!*******************************************************************************
!
!! DPMPAR provides double precision machine constants.
!
!
!  Author:
!
!    Alfred H Morris, Jr,
!    Naval Surface Weapons Center,
!    Dahlgren, Virginia.
!
!     It is assumed that the argument
!     I is an integer having one of the values 1, 2, or 3. if the
!     single precision arithmetic being used has M base B digits and
!     its smallest and largest exponents are EMIN and EMAX, then
!
!        DPMPAR(1) = B**(1 - M), the machine precision,
!
!        DPMPAR(2) = B**(EMIN - 1), the smallest magnitude,
!
!        DPMPAR(3) = B**EMAX*(1 - B**(-M)), the largest magnitude.
!
  double precision b
  double precision binv
  double precision bm1
  double precision dpmpar
  integer i

  double precision w,z
  integer emax,emin,ibeta,m
  integer ipmpar
!
  if ( i == 1 ) then

    b = ipmpar(4)
    m = ipmpar(8)
    dpmpar = b**(1-m)

  else if ( i == 2 ) then

    b = ipmpar(4)
    emin = ipmpar(9)
    binv = 1.0 / b
    w = b**( emin + 2 )
    dpmpar = (( w * binv ) * binv ) * binv

  else if ( i == 3 ) then

    ibeta = ipmpar(4)
    m = ipmpar(8)
    emax = ipmpar(10)

    b = ibeta
    bm1 = ibeta - 1
    z = b**( m - 1 )
    w = (( z - 1.0 ) * b + bm1 ) / ( b * z )

    z = b**( emax - 2 )
    dpmpar = (( w * z ) * b ) * b

  else

    write ( *, * ) ' '
    write ( *, * ) 'DPMPAR - Fatal error!'
    write ( *, * ) '  Illegal input argument I = ', i
    stop

  end if

  return
end
function dstrem ( z )
!
!*******************************************************************************
!
!! DSTREM computes the Sterling remainder Log(Gamma(Z)) - Sterling(Z).
!
!
!  Function:
!
!    Returns Log(Gamma(Z)) - Sterling(Z)  where   Sterling(Z)  is
!    Sterling's Approximation to Log(Gamma(Z))
!
!    Sterling(Z) = LOG( SQRT( 2*PI ) ) + ( Z-0.5 ) * LOG( Z ) - Z
!
!  Method:
!
!    If Z >= 6 uses 9 terms of series in Bernoulli numbers
!    (Values calculated using Maple)
!    Otherwise computes difference explicitly
!
!  Parameters:
!
!     Z --> Value at which Sterling remainder calculated
!           Must be positive.
!                  double precision Z
!
  double precision, parameter :: hln2pi = 0.91893853320467274178d0
  integer, parameter :: ncoef = 9
!
  double precision, parameter, dimension ( 0:ncoef ) :: coef = (/ &
    0.0d0, &
    0.0833333333333333333333333333333d0, &
    -0.00277777777777777777777777777778d0, &
    0.000793650793650793650793650793651d0, &
    -0.000595238095238095238095238095238d0, &
    0.000841750841750841750841750841751d0, &
    -0.00191752691752691752691752691753d0, &
    0.00641025641025641025641025641026d0, &
    -0.0295506535947712418300653594771d0, &
    0.179644372368830573164938490016d0 /)
  double precision eval_pol
  double precision dlngam
  double precision dstrem
  double precision sterl
  double precision z
!
  if ( z <= 0.0d0 ) then
    write ( *, * ) ' '
    write ( *, * ) 'DTSTREM - Fatal error!'
    write ( *, * ) '  Zero or negative argument.'
    stop
  end if

  if ( z > 6.0d0 ) then
    dstrem = eval_pol ( coef, ncoef, 1.0d0/z**2 ) * z
  else
    sterl = hln2pi + ( z - 0.5d0 ) * log ( z ) - z
!    dstrem = dlngam ( z ) - sterl
  end if

  return
end
function dt1 ( p, q, df )
!
!*******************************************************************************
!
!! DT1 computes an approximate inverse of the cumulative T distribution.
!
!
!  Discussion:
!
!    Returns  the  inverse   of  the T   distribution   function, i.e.,
!    the integral from 0 to INVT of the T density is P. This is an
!    initial approximation
!
!  Parameters:
!
!     P --> The p-value whose inverse from the T distribution is
!          desired.
!                    P is double precision
!
!     Q --> 1-P.
!                    Q is double precision
!
!     DF --> Degrees of freedom of the T distribution.
!                    DF is double precision
!
!
  double precision coef(0:4,4)
  double precision, parameter, dimension ( 4 ) :: denom = (/ &
    4.0d0, 96.0d0, 384.0d0, 92160.0d0 /)
  double precision denpow
  double precision eval_pol
  double precision df
  double precision dinvnr
  double precision dt1
  integer i
  integer, parameter, dimension ( 4 ) :: ideg = (/ 1, 2, 3, 4 /)
  double precision p
  double precision q
  double precision sum
  double precision term
  double precision x
  double precision xp
  double precision xx
!
  data (coef(i,1),i=0,4)/1.0d0,1.0d0,3*0.0d0/
  data (coef(i,2),i=0,4)/3.0d0,16.0d0,5.0d0,2*0.0d0/
  data (coef(i,3),i=0,4)/-15.0d0,17.0d0,19.0d0,3.0d0,0.0d0/
  data (coef(i,4),i=0,4)/-945.0d0,-1920.0d0,1482.0d0,776.0d0,79.0d0/
!
  x = abs ( dinvnr ( p, q ) )
  xx = x * x

  sum = x
  denpow = 1.0d0
  do i = 1, 4
    term = eval_pol ( coef(0,i), ideg(i), xx ) * x
    denpow = denpow * df
    sum = sum + term / ( denpow * denom(i) )
  end do

  if ( p >= 0.5d0 ) then
    xp = sum
  else
    xp = - sum
  end if

  dt1 = xp

  return
end
subroutine dzror ( status, x, fx, xlo, xhi, qleft, qhi )
!
!*******************************************************************************
!
!! DZROR seeks a zero of a function, using reverse communication.
!
!
!  Function
!
!    Performs the zero finding.  STZROR must have been called before
!    this routine in order to set its parameters.
!
!  Parameters:
!
!
!     STATUS <--> At the beginning of a zero finding problem, STATUS
!                 should be set to 0 and ZROR invoked.  (The value
!                 of other parameters will be ignored on this call.)
!
!                 When ZROR needs the function evaluated, it will set
!                 STATUS to 1 and return.  The value of the function
!                 should be set in FX and ZROR again called without
!                 changing any of its other parameters.
!
!                 When ZROR has finished without error, it will return
!                 with STATUS 0.  In that case (XLO,XHI) bound the answe
!
!                 If ZROR finds an error (which implies that F(XLO)-Y an
!                 F(XHI)-Y have the same sign, it returns STATUS -1.  In
!                 this case, XLO and XHI are undefined.
!                         integer STATUS
!
!     X <-- The value of X at which F(X) is to be evaluated.
!                         double precision X
!
!     FX --> The value of F(X) calculated when ZROR returns with
!            STATUS = 1.
!                         double precision FX
!
!     XLO <-- When ZROR returns with STATUS = 0, XLO bounds the
!             inverval in X containing the solution below.
!                         double precision XLO
!
!     XHI <-- When ZROR returns with STATUS = 0, XHI bounds the
!             inverval in X containing the solution above.
!                         double precision XHI
!
!     QLEFT <-- .TRUE. if the stepping search terminated unsucessfully
!                at XLO.  If it is .FALSE. the search terminated
!                unsucessfully at XHI.
!                    QLEFT is LOGICAL
!
!     QHI <-- .TRUE. if F(X) .GT. Y at the termination of the
!              search and .FALSE. if F(X) .LT. Y at the
!              termination of the search.
!                    QHI is LOGICAL
!
  double precision a
  double precision abstol
  double precision b
  double precision c
  double precision d
  integer ext
  double precision fa
  double precision fb
  double precision fc
  double precision fd
  double precision fda
  double precision fdb
  logical first
  double precision ftol
  double precision fx
  integer i99999
  double precision m
  double precision mb
  double precision p
  double precision q
  logical qhi
  logical qleft
  logical qrzero
  double precision reltol
  integer status
  double precision tol
  double precision w
  double precision x
  double precision xhi
  double precision xlo
  double precision, save :: xxhi = 0.0
  double precision, save :: xxlo = 0.0
  double precision zabstl
  double precision zreltl
  double precision zx
  double precision zxhi
  double precision zxlo
!
  save
!
  ftol(zx) = 0.5d0 * max ( abstol, reltol * abs ( zx ) )
!
  if (status>0) go to 280
  xlo = xxlo
  xhi = xxhi
  b = xlo
  x = xlo
!
!     GET-function-VALUE
!
  assign 10 to i99999
  go to 270

   10 fb = fx
  xlo = xhi
  a = xlo
  x = xlo
!
!     GET-function-VALUE
!
  assign 20 to i99999
  go to 270
!
!     Check that F(ZXLO) < 0 < F(ZXHI)  or
!                F(ZXLO) > 0 > F(ZXHI)
!
   20 if (.not. (fb<0.0d0)) go to 40
  if (.not. (fx<0.0d0)) go to 30
  status = -1
  qleft = fx < fb
  qhi = .false.
  return

   30 continue
   40 if (.not. (fb>0.0d0)) go to 60
  if (.not. (fx>0.0d0)) go to 50
  status = -1
  qleft = fx > fb
  qhi = .true.
  return

   50 continue
   60 fa = fx

  first = .true.
   70 c = a
  fc = fa
  ext = 0
   80 if (.not. ( abs ( fc ) < abs ( fb ) ) ) go to 100
  if (.not. (c/=a)) go to 90
  d = a
  fd = fa
   90 a = b
  fa = fb
  xlo = c
  b = xlo
  fb = fc
  c = a
  fc = fa
  100 tol = ftol(xlo)
  m = (c+b)*.5d0
  mb = m - b
  if (.not. ( abs ( mb ) > tol ) ) go to 240
  if (.not. (ext>3)) go to 110
  w = mb
  go to 190

  110 tol = sign(tol,mb)
  p = (b-a)*fb

  if ( first ) then
    q = fa - fb
    first = .false.
  else
    fdb = ( fd - fb ) / ( d - b )
    fda = ( fd - fa ) / ( d - a )
    p = fda * p
    q = fdb * fa - fda * fb
  end if

  130 if (.not. (p<0.0d0)) go to 140
  p = -p
  q = -q
  140 if (ext==3) p = p*2.0d0
  if (.not. ((p*1.0d0)==0.0d0.or.p<= (q*tol))) go to 150
  w = tol
  go to 180

  150 if (.not. (p< (mb*q))) go to 160
  w = p/q
  go to 170

  160 w = mb
  170 continue
  180 continue
  190 d = a
  fd = fa
  a = b
  fa = fb
  b = b + w
  xlo = b
  x = xlo
!
!     GET-function-VALUE
!
  assign 200 to i99999
  go to 270

  200 continue

  fb = fx

  if ( fc * fb >= 0.0d0 ) then

    go to 70

  else

    if ( w == mb ) then
      ext = 0
    else
      ext = ext + 1
    end if

    go to 80

  end if

  240 xhi = c
  qrzero = (fc>=0.0d0 .and. fb<=0.0d0) .or. (fc<0.0d0 .and. fb>=0.0d0)

  if ( qrzero ) then
    status = 0
  else
    status = -1
  end if

  return

  entry dstzr ( zxlo, zxhi, zabstl, zreltl )
!
!*******************************************************************************
!
!! DSTZR - SeT ZeRo finder - Reverse communication version
!
!
!  Function:
!
!    Sets quantities needed by ZROR.  The function of ZROR
!    and the quantities set is given here.
!
!    Given a function F, find XLO such that F(XLO) = 0.
!
!     Input condition. F is a double precision function of a single
!     double precision argument and XLO and XHI are such that
!          F(XLO)*F(XHI)  .LE.  0.0
!
!     If the input condition is met, QRZERO returns .TRUE.
!     and output values of XLO and XHI satisfy the following
!          F(XLO)*F(XHI)  .LE. 0.
!          ABS(F(XLO)  .LE. ABS(F(XHI)
!          ABS(XLO-XHI)  .LE. TOL(X)
!     where
!          TOL(X) = MAX(ABSTOL,RELTOL*ABS(X))
!
!     If this algorithm does not find XLO and XHI satisfying
!     these conditions then QRZERO returns .FALSE.  This
!     implies that the input condition was not met.
!
!  Parameters:
!
!
!     XLO --> The left endpoint of the interval to be
!           searched for a solution.
!                    XLO is double precision
!
!     XHI --> The right endpoint of the interval to be
!           for a solution.
!                    XHI is double precision
!
!     ABSTOL, RELTOL --> Two numbers that determine the accuracy
!                      of the solution.  See function for a
!                      precise definition.
!                    ABSTOL is double precision
!                    RELTOL is double precision
!
!
!                              Method
!
!     Algorithm R of the paper 'Two Efficient Algorithms with
!     Guaranteed Convergence for Finding a Zero of a Function'
!     by J. C. P. Bus and T. J. Dekker in ACM Transactions on
!     Mathematical Software, Volume 1, no. 4 page 330
!     (Dec. '75) is employed to find the zero of F(X)-Y.
!
  xxlo = zxlo
  xxhi = zxhi
  abstol = zabstl
  reltol = zreltl
  return
!
!     TO GET-function-VALUE
!
  270 status = 1
  return

  280 continue
  go to i99999

end
function erf ( x )
!
!*******************************************************************************
!
!! ERF evaluates the error function.
!
!
  double precision, parameter :: c = 0.564189583547756d0
!
  double precision, parameter, dimension ( 5 ) :: a = (/ &
    0.771058495001320d-04, &
    -0.133733772997339d-02, &
    0.323076579225834d-01, &
    0.479137145607681d-01, &
    0.128379167095513d+00 /)

  double precision ax
  double precision, parameter, dimension ( 3 ) :: b = (/ &
    0.301048631703895d-02, &
    0.538971687740286d-01, &
    0.375795757275549d+00 /)
  double precision bot
  double precision erf
  double precision p(8)
  double precision q(8)
  double precision r(5)
  double precision, parameter, dimension ( 4 ) :: s = (/ &
    9.41537750555460d+01, &
    1.87114811799590d+02, &
    9.90191814623914d+01, &
    1.80124575948747d+01 /)
  double precision t
  double precision top
  double precision x
  double precision x2
!
  data p(1)/-1.36864857382717d-07/,p(2)/5.64195517478974d-01/
  data p(3)/7.21175825088309d+00/,p(4)/4.31622272220567d+01/
  data p(5)/1.52989285046940d+02/,p(6)/3.39320816734344d+02/
  data p(7)/4.51918953711873d+02/,p(8)/3.00459261020162d+02/
  data q(1)/1.00000000000000d+00/,q(2)/1.27827273196294d+01/
  data q(3)/7.70001529352295d+01/,q(4)/2.77585444743988d+02/
  data q(5)/6.38980264465631d+02/,q(6)/9.31354094850610d+02/
  data q(7)/7.90950925327898d+02/,q(8)/3.00459260956983d+02/
  data r(1)/2.10144126479064d+00/,r(2)/2.62370141675169d+01/
  data r(3)/2.13688200555087d+01/,r(4)/4.65807828718470d+00/
  data r(5)/2.82094791773523d-01/
!
  ax = abs ( x )

  if ( ax <= 0.5d0 ) then

    t = x * x

    top = (((( a(1) * t + a(2) ) * t + a(3) ) * t + a(4) ) * t &
      + a(5) ) + 1.0d0

    bot = (( b(1) * t + b(2) ) * t + b(3) ) * t + 1.0d0
    erf = ax * ( top / bot )

  else if ( ax <= 4.0d0 ) then

    top = (((((( p(1) * ax + p(2) ) * ax + p(3) ) * ax + p(4) ) * ax &
      + p(5) ) * ax + p(6) ) * ax + p(7) ) * ax  + p(8)

    bot = (((((( q(1) * ax + q(2) ) * ax + q(3) ) * ax + q(4) ) * ax &
      + q(5) ) * ax + q(6) ) * ax + q(7) ) * ax + q(8)

    erf = 0.5d0 + ( 0.5d0 - exp ( - x * x ) * top / bot )

  else if ( ax < 5.8d0 ) then

    x2 = x * x
    t = 1.0d0 / x2

    top = ((( r(1) * t + r(2) ) * t + r(3) ) * t + r(4) ) * t + r(5)

    bot = ((( s(1) * t + s(2) ) * t + s(3) ) * t + s(4) ) * t + 1.0d0

    erf = ( c - top / ( x2 * bot )) / ax
    erf = 0.5d0 + ( 0.5d0 - exp ( - x2 ) * erf )

  else

    erf = 1.0

  end if

  if ( x < 0.0 ) then
    erf = - erf
  end if

  return
end
function erfc1 ( ind, x )
!
!*******************************************************************************
!
!! ERFC1 evaluates the complementary error function.
!
!
!  Modified:
!
!    09 December 1999
!
!  Parameters:
!
!    Input, integer IND, chooses the scaling.
!    If IND is nonzero, then the value returned has been multiplied by
!    EXP(X*X).
!
!    Input, double precision X, the argument of the function.
!
!    Output, double precision ERFC1, the value of the complementary error function.
!
  double precision, parameter :: c = 0.564189583547756d0
!
  double precision a(5)
  double precision ax
  double precision b(3)
  double precision bot
  double precision e
  double precision erfc1
  double precision exparg
  integer ind
  double precision p(8)
  double precision q(8)
  double precision r(5)
  double precision s(4)
  double precision t
  double precision top
  double precision w
  double precision x
!
  data a(1)/.771058495001320d-04/
  data a(2)/-.133733772997339d-02/
  data a(3)/.323076579225834d-01/
  data a(4)/.479137145607681d-01/
  data a(5)/.128379167095513d+00/
  data b(1)/.301048631703895d-02/
  data b(2)/.538971687740286d-01/
  data b(3)/.375795757275549d+00/
  data p(1)/-1.36864857382717d-07/,p(2)/5.64195517478974d-01/
  data p(3)/7.21175825088309d+00/,p(4)/4.31622272220567d+01/
  data p(5)/1.52989285046940d+02/,p(6)/3.39320816734344d+02/
  data p(7)/4.51918953711873d+02/,p(8)/3.00459261020162d+02/
  data q(1)/1.00000000000000d+00/,q(2)/1.27827273196294d+01/
  data q(3)/7.70001529352295d+01/,q(4)/2.77585444743988d+02/
  data q(5)/6.38980264465631d+02/,q(6)/9.31354094850610d+02/
  data q(7)/7.90950925327898d+02/,q(8)/3.00459260956983d+02/
  data r(1)/2.10144126479064d+00/,r(2)/2.62370141675169d+01/
  data r(3)/2.13688200555087d+01/,r(4)/4.65807828718470d+00/
  data r(5)/2.82094791773523d-01/
  data s(1)/9.41537750555460d+01/,s(2)/1.87114811799590d+02/
  data s(3)/9.90191814623914d+01/,s(4)/1.80124575948747d+01/
!
!  ABS(X) .LE. 0.5
!
  ax = abs ( x )

  if ( ax <= 0.5d0 ) then

    t = x * x
!
!  ERROR?  Did we lose a factor of T multiplying A(5)???
!  JVB
!
    top = (((( a(1) * t + a(2) ) * t + a(3) ) * t + a(4) ) * t + a(5) ) + 1.0d0

    bot = (( b(1) * t + b(2) ) * t + b(3) ) * t + 1.0d0

    erfc1 = 0.5d0 + ( 0.5d0 - x * ( top / bot ) )

    if ( ind /= 0 ) then
      erfc1 = exp ( t ) * erfc1
    end if

    return

  end if
!
!  0.5 .LT. ABS(X) .LE. 4
!
  if ( ax <= 4.0d0 ) then

    top = (((((( p(1) * ax + p(2)) * ax + p(3)) * ax + p(4)) * ax &
      + p(5)) * ax + p(6)) * ax + p(7)) * ax + p(8)

    bot = (((((( q(1) * ax + q(2)) * ax + q(3)) * ax + q(4)) * ax &
      + q(5)) * ax + q(6)) * ax + q(7)) * ax + q(8)

    erfc1 = top / bot
!
!  ABS(X) .GT. 4
!
  else

    if ( x <= - 5.6d0 ) then

      if ( ind == 0 ) then
        erfc1 = 2.0d0
      else
        erfc1 = 2.0d0 * exp ( x * x )
      end if

      return

    end if

    if ( ind == 0 ) then

      if ( x > 100.0d0 ) then
        erfc1 = 0.0
        return
      end if

      if ( x * x > - exparg ( 1 ) ) then
        erfc1 = 0.0
        return
      end if

    end if

    t = ( 1.0d0 / x ) **2

    top = ((( r(1) * t + r(2) ) * t + r(3) ) * t + r(4) ) * t + r(5)

    bot = ((( s(1) * t + s(2) ) * t + s(3) ) * t + s(4) ) * t + 1.0d0

    erfc1 = ( c - t * top / bot ) / ax

  end if
!
!  Final assembly.
!
  if ( ind /= 0 ) then

    if ( x < 0.0d0 ) then
      erfc1 = 2.0d0 * exp ( x * x ) - erfc1
    end if

  else

    w = x * x
    t = w
    e = w - t
    erfc1 = (( 0.5d0 + ( 0.5d0 - e )) * exp ( - t ) ) * erfc1

    if ( x < 0.0d0 ) then
      erfc1 = 2.0d0 - erfc1
    end if

  end if

  return
end
function esum ( mu, x )
!
!*******************************************************************************
!
!! ESUM evaluates EXP(MU + X).
!
!
  double precision esum
  integer mu
  double precision w
  double precision x
!
  if ( x <= 0.0d0 ) then
    if ( mu >= 0 ) then
      w = mu + x
      if ( w <= 0.0d0 ) then
        esum = exp ( w )
        return
      end if
    end if
  else if ( x > 0.0 ) then
    if ( mu <= 0 ) then
      w = mu + x
      if ( w >= 0.0d0 ) then
        esum = exp ( w )
        return
      end if
    end if
  end if

  w = mu
  esum = exp ( w ) * exp ( x )

  return
end
function exparg ( l )
!
!*******************************************************************************
!
!! EXPARG returns the largest or smallest legal argument for EXP.
!
!
!  Discussion:
!
!    Only an approximate limit for the argument of EXP is desired.
!
!  Modified:
!
!    09 December 1999
!
!  Parameters:
!
!    Input, integer L, indicates which limit is desired.
!    If L = 0, then the largest positive argument for EXP is desired.
!    Otherwise, the largest negative argument for EXP for which the
!    result is nonzero is desired.
!
  integer b
  double precision exparg
  integer ipmpar
  integer l
  double precision lnb
  integer m
!
!  Get the arithmetic base.
!
  b = ipmpar(4)
!
!  Compute the logarithm of the arithmetic base.
!
  if ( b == 2 ) then
    lnb = 0.69314718055995d0
  else if ( b == 8 ) then
    lnb = 2.0794415416798d0
  else if ( b == 16 ) then
    lnb = 2.7725887222398d0
  else
    lnb = log ( dble ( b ) )
  end if

  if ( l /= 0 ) then
    m = ipmpar(9) - 1
    exparg = 0.99999d0 * ( m * lnb )
  else
    m = ipmpar(10)
    exparg = 0.99999d0 * ( m * lnb )
  end if

  return
end
function fpser ( a, b, x, eps )
!
!*******************************************************************************
!
!! FPSER evaluates IX(A,B)
!
!          FOR B .LT. MIN(EPS,EPS*A) AND X .LE. 0.5.
!
!                  SET  FPSER = X**A
!
!
  double precision a
  double precision an
  double precision b
  double precision c
  double precision eps
  double precision exparg
  double precision fpser
  double precision s
  double precision t
  double precision tol
  double precision x
!
  fpser = 1.0d0

  if ( a > 1.0d-3 * eps ) then
    fpser = 0.0d0
    t = a * log ( x )
    if ( t < exparg ( 1 ) ) then
      return
    end if
    fpser = exp ( t )
  end if
!
!  1/B(A,B) = B
!
  fpser = ( b / a ) * fpser
  tol = eps / a
  an = a + 1.0d0
  t = x
  s = t / an

   20 continue

  an = an + 1.0d0
  t = x * t
  c = t / an
  s = s + c
  if ( abs ( c ) > tol ) then
    go to 20
  end if

  fpser = fpser * ( 1.0d0 + a * s )

  return
end
function gam1 ( a )
!
!*******************************************************************************
!
!! GAM1 computes 1 / GAMMA(A+1) - 1 for -0.5 <= A <= 1.5
!
!
  double precision a
  double precision bot
  double precision d
  double precision gam1
  double precision, parameter, dimension ( 7 ) :: p = (/ &
     0.577215664901533d+00, &
    -0.409078193005776d+00, &
    -0.230975380857675d+00, &
     0.597275330452234d-01, &
     0.766968181649490d-02, &
    -0.514889771323592d-02, &
     0.589597428611429d-03 /)
  double precision q(5)
  double precision r(9)
  double precision, parameter :: s1 = 0.273076135303957d+00
  double precision, parameter :: s2 = 0.559398236957378d-01
  double precision t
  double precision top
  double precision w
!
  data q(1)/.100000000000000d+01/,q(2)/.427569613095214d+00/
  data q(3)/.158451672430138d+00/,q(4)/.261132021441447d-01/
  data q(5)/.423244297896961d-02/
  data r(1)/-.422784335098468d+00/,r(2)/-.771330383816272d+00/
  data r(3)/-.244757765222226d+00/,r(4)/.118378989872749d+00/
  data r(5)/.930357293360349d-03/,r(6)/-.118290993445146d-01/
  data r(7)/.223047661158249d-02/,r(8)/.266505979058923d-03/
  data r(9)/-.132674909766242d-03/
!
  d = a - 0.5d0

  if ( d > 0.0d0 ) then
    t = d - 0.5d0
  else
    t = a
  end if

  if ( t == 0.0 ) then

    gam1 = 0.0d0

  else if ( t > 0.0 ) then

    top = ((((( p(7) * t + p(6) ) * t + p(5) ) * t + p(4) ) * t &
      + p(3) ) * t + p(2) ) * t + p(1)

    bot = ((( q(5) * t + q(4) ) * t + q(3) ) * t + q(2) ) * t + 1.0d0

    w = top / bot

    if ( d <= 0.0d0 ) then
      gam1 = a * w
    else
      gam1 = ( t / a ) * (( w - 0.5d0 ) - 0.5d0 )
    end if

  else if ( t < 0.0 ) then

    top = ((((((( r(9) * t + r(8)) * t + r(7)) * t + r(6)) * t &
      + r(5)) * t + r(4)) * t + r(3)) * t + r(2)) * t + r(1)

    bot = ( s2 * t + s1 ) * t + 1.0d0
    w = top / bot

    if ( d <= 0.0d0 ) then
      gam1 = a * (( w + 0.5d0 ) + 0.5d0 )
    else
      gam1 = t * w / a
    end if

  end if

  return
end
subroutine gamma_inc_inv ( a, x, x0, p, q, ierr )
!
!*******************************************************************************
!
!! GAMMA_INC_INV computes the inverse incomplete gamma ratio function.
!
!
!  Discussion:
!
!    The routine is given positive A, and nonnegative P and Q where P + Q = 1.
!    The value X is computed with the property that P(A,X) = P and Q(A,X) = Q.  
!    Schroder iteration is employed.  The routine attempts to compute X
!    to 10 significant digits if this is possible for the particular computer 
!    arithmetic being used.
!
!  Author:
!
!    Alfred H Morris, Jr,
!    Naval Surface Weapons Center,
!    Dahlgren, Virginia.
!
!  Parameters:
!
!    Input, double precision A, the parameter in the incomplete gamma
!    ratio.  A must be positive.
!
!    Output, double precision X, the computed point for which the
!    incomplete gamma functions have the values P and Q.
!
!    Input, double precision X0, an optional initial approximation
!    for the solution X.  If the user does not want to supply an
!    initial approximation, then X0 should be set to 0, or a negative
!    value.
!
!    Input, double precision P, Q, the values of the incomplete gamma
!    functions, for which the corresponding argument is desired.
!
!    Output, integer IERR, error flag.
!
!    0, the solution was obtained. Iteration was not used.
!    K > 0, The solution was obtained. IERR iterations were performed.
!    -2, A .LE. 0
!    -3, No solution was obtained. The ratio Q/A is too large.
!    -4, P + Q .NE. 1
!    -6, 20 iterations were performed. The most recent value obtained 
!        for X is given.  This cannot occur if X0 .LE. 0.
!    -7, Iteration failed. No value is given for X.
!        This may occur when X is approximately 0.
!    -8, A value for X has been obtained, but the routine is not certain
!        of its accuracy.  Iteration cannot be performed in this
!        case. If X0 .LE. 0, this can occur only when P or Q is 
!        approximately 0. If X0 is positive then this can occur when A is
!        exceedingly close to X and A is extremely large (say A .GE. 1.E20).
!
  double precision, parameter :: a0 = 3.31125922108741d0
  double precision, parameter :: a1 = 11.6616720288968d0
  double precision, parameter :: a2 = 4.28342155967104d0
  double precision, parameter :: a3 = 0.213623493715853d0
  double precision, parameter :: b1 = 6.61053765625462d0
  double precision, parameter :: b2 = 6.40691597760039d0
  double precision, parameter :: b3 = 1.27364489782223d0
  double precision, parameter :: b4 = .036117081018842d0
  double precision, parameter :: c = 0.577215664901533d0
  double precision, parameter :: ln10 = 2.302585d0
  double precision, parameter :: tol = 1.0d-5
!
  double precision a
  double precision alnrel
  double precision am1
  double precision amax
  double precision amin(2)
  double precision ap1
  double precision ap2
  double precision ap3
  double precision apn
  double precision b
  double precision bmin(2)
  double precision c1
  double precision c2
  double precision c3
  double precision c4
  double precision c5
  double precision d
  double precision dmin(2)
  double precision dpmpar
  double precision e
  double precision e2
  double precision emin(2)
  double precision eps
  double precision eps0(2)
  double precision g
  double precision gamma_log
  double precision gamma_ln1
  double precision gamma
  double precision h
  integer ierr
  integer iop
  double precision p
  double precision pn
  double precision q
  double precision qg
  double precision qn
  double precision r
  double precision rcomp
  double precision rta
  double precision s
  double precision s2
  double precision sum
  double precision t
  double precision u
  double precision w
  double precision x
  double precision x0
  double precision xmax
  double precision xmin
  double precision xn
  double precision y
  double precision z
!
  data eps0(1)/1.d-10/
  data eps0(2)/1.d-08/
  data amin(1)/500.0d0/
  data amin(2)/100.0d0/
  data bmin(1)/1.d-28/
  data bmin(2)/1.d-13/
  data dmin(1)/1.d-06/
  data dmin(2)/1.d-04/
  data emin(1)/2.d-03/
  data emin(2)/6.d-03/
!
!  E, XMIN, and XMAX are machine dependent constants.
!  E is the smallest number for which 1.0 + E .GT. 1.0.
!  XMIN is the smallest positive number.
!  XMAX is the largest positive number.
!
  e = epsilon ( 1.0d0 )
  xmin = dpmpar(2)
  xmax = dpmpar(3)

  x = 0.0d0

  if ( a <= 0.0d0 ) then
    ierr = - 2
    return
  end if

  t = dble ( p ) + dble ( q ) - 1.d0

  if ( abs ( t ) > e ) then
    ierr = - 4
    return
  end if

  ierr = 0

  if ( p == 0.0d0 ) then
    return
  end if

  if ( q == 0.0d0 ) then
    x = xmax
    return
  end if

  if ( a == 1.0d0 ) then
    if ( q >= 0.9d0 ) then
      x = - alnrel ( - p )
    else
      x = - log ( q )
    end if
    return
  end if

  e2 = 2.0d0 * e
  amax = 0.4d-10 / ( e * e )

  if ( e > 1.0d-10 ) then
    iop = 2
  else
    iop = 1
  end if

  eps = eps0(iop)
  xn = x0

  if ( x0 > 0.0d0 ) go to 160
!
!  Selection of the initial approximation XN OF X when A .LT. 1
!
  if ( a > 1.0d0 ) go to 80

  g = gamma ( a + 1.0d0 )
  qg = q * g

  if ( qg == 0.0d0 ) go to 360

  b = qg / a

  if ( qg > 0.6d0 * a ) go to 40

  if ( a < 0.30d0 .and. b >= 0.35d0 ) then
    t = exp ( - ( b + c ))
    u = t * exp ( t )
    xn = t * exp ( u )
    go to 160
  end if

  if ( b >= 0.45d0 ) go to 40
  if ( b == 0.0d0 ) go to 360
  y = - log ( b ) 
  s = 0.5d0 + ( 0.5d0 - a )
  z = log ( y ) 
  t = y - s * z

  if ( b >= 0.15d0 ) then
    xn = y - s * log ( t ) - log ( 1.0d0 + s / ( t + 1.0d0 ))
    go to 220
  end if

  if ( b > 0.01d0 ) then
    u = (( t + 2.0d0 * ( 3.0d0 - a ) ) * t + ( 2.0d0 - a ) * ( 3.0d0 - a )) / &
      (( t + ( 5.0d0 - a )) * t + 2.0d0 )
    xn = y - s * log ( t ) - log ( u )
    go to 220
  end if

   30 continue

  c1 = - s * z
  c2 = - s * ( 1.0d0 + c1 )

  c3 = s * (( 0.5d0 * c1 + ( 2.0d0 - a ) ) * c1 + ( 2.5d0 - 1.5d0 * a ) )

  c4 = - s * (((c1/3.0d0 + ( 2.5d0 - 1.5d0*a ) ) *c1 + ((a-6.0d0)*a+7.0d0)) &
    * c1 + ((11.0d0*a-46)*a+47.0d0)/6.0d0 )

  c5 = - s * (((( - c1 / 4.0d0 + ( 11.0d0 * a - 17.0d0 ) / 6.0d0 ) * c1 &
     + ((-3.0d0*a+13.0d0)*a-13.0d0)) * c1 &
     + 0.5d0* (((2.0d0*a-25.0d0)*a+72.0d0)*a-61.0d0)) * c1 &
     + (((25.0d0*a-195.0d0)*a+477.0d0)*a-379.0d0) / 12.0d0 )

  xn = (((( c5 / y + c4 ) / y + c3 ) / y + c2 ) / y + c1 ) + y

  if ( a > 1.0d0 ) then
    go to 220
  end if

  if ( b > bmin(iop) ) then
    go to 220
  end if

  x = xn
  return

   40 continue

  if ( b * q <= 1.0d-8 ) then
    xn = exp ( - ( q / a + c ))
  else if ( p > 0.9d0 ) then
    xn = exp (( alnrel ( - q ) + gamma_ln1 ( a )) / a )
  else
    xn = exp ( log ( p * g ) / a )
  end if

  if ( xn == 0.0d0 ) then
    ierr = - 3
    return
  end if

  t = 0.5d0 + ( 0.5d0 - xn / ( a + 1.0d0 ))
  xn = xn / t
  go to 160
!
!  Selection of the initial approximation XN of X when A > 1.
!
   80 continue

  if ( q > 0.5d0 ) then
    w = log ( p )
  else
    w = log ( q )
  end if

  t = sqrt ( - 2.0d0 * w )

  s = t - ((( a3 * t + a2 ) * t + a1 ) * t + a0 ) / (((( &
    b4 * t + b3 ) * t + b2 ) * t + b1 ) * t + 1.0d0 )

  if ( q > 0.5d0 ) then
    s = - s
  end if

  rta = sqrt ( a )
  s2 = s * s

  xn = a + s * rta + ( s2 - 1.0d0 ) / 3.0d0 + s * (s2-7.0d0) / (36.0d0*rta) &
     - ( (3.0d0*s2+7.0d0)*s2-16.0d0)/ (810.0d0*a) &
     + s * ((9.0d0*s2+256.0d0)*s2-433.0d0)/ (38880.0d0*a*rta)

  xn = max ( xn, 0.0d0 )

  if ( a >= amin(iop) ) then

    x = xn
    d = 0.5d0 + ( 0.5d0 - x / a )

    if ( abs ( d ) <= dmin(iop) ) then
      return
    end if

  end if

  110 continue

  if ( p <= 0.5d0 ) then
    go to 130
  end if

  if ( xn < 3.0d0 * a ) then
    go to 220
  end if

  y = - ( w + gamma_log ( a ) )
  d = max ( 2.0d0, a * ( a - 1.0d0 ) )

  if ( y >= ln10 * d ) then
    s = 1.0d0 - a
    z = log ( y )
    go to 30
  end if

  120 continue

  t = a - 1.0d0
  xn = y + t * log ( xn ) - alnrel ( - t / ( xn + 1.0d0 ) )
  xn = y + t * log ( xn ) - alnrel ( - t / ( xn + 1.0d0 ) )
  go to 220

  130 continue

  ap1 = a + 1.0d0

  if ( xn > 0.70d0 * ap1 ) then
    go to 170
  end if

  w = w + gamma_log ( ap1 )

  if ( xn <= 0.15d0 * ap1 ) then
    ap2 = a + 2.0d0
    ap3 = a + 3.0d0
    x = exp (( w + x ) / a )
    x = exp (( w + x - log ( 1.0d0+ ( x / ap1 ) * ( 1.0d0 + x / ap2 )))/ a )
    x = exp (( w + x - log ( 1.0d0+ ( x / ap1 ) * ( 1.0d0 + x / ap2 ))) / a )
    x = exp (( w + x - log ( 1.0d0+ ( x / ap1 ) * ( 1.0d0 + ( x / ap2 ) &
      * ( 1.0d0 + x / ap3 )))) / a )
    xn = x

    if ( xn <= 1.0d-2 * ap1 ) then
      if ( xn <= emin(iop) * ap1 ) then
        return
      end if
      go to 170
    end if

  end if

  apn = ap1
  t = xn / apn
  sum = 1.0d0 + t

  150 continue

  apn = apn + 1.0d0
  t = t * ( xn / apn )
  sum = sum + t
  if ( t > 1.0d-4 ) then
    go to 150
  end if

  t = w - log ( sum )
  xn = exp (( xn + t ) / a )
  xn = xn * ( 1.0d0 - ( a * log ( xn ) - xn - t ) / ( a - xn ))
  go to 170
!
!  Schroder iteration using P.
!
  160 continue

  if ( p > 0.5d0 ) then
    go to 220
  end if

  170 continue

  if ( p <= 1.0d10 * xmin ) then
    go to 350
  end if

  am1 = ( a - 0.5d0 ) - 0.5d0

  180 continue

  if ( a > amax ) then
    d = 0.5d0 + ( 0.5d0 - xn / a )
    if ( abs ( d ) <= e2 ) then
      go to 350
    end if
  end if

  190 continue

  if ( ierr >= 20 ) then
    ierr = -6
    return
  end if

  ierr = ierr + 1
  call gamma_inc ( a, xn, pn, qn, 0 )

  if ( pn == 0.0d0 .or. qn == 0.0d0 ) then
    go to 350
  end if

  r = rcomp ( a, xn )

  if ( r == 0.0d0 ) then
    go to 350
  end if

  t = ( pn - p ) / r
  w = 0.5d0 * ( am1 - xn )

  if ( abs ( t ) <= 0.1d0 .and. abs ( w * t ) <= 0.1d0 ) then
    go to 200
  end if

  x = xn * ( 1.0d0 - t )

  if ( x <= 0.0d0 ) then
    go to 340
  end if

  d = abs ( t )
  go to 210

  200 continue

  h = t * ( 1.0d0 + w * t )
  x = xn * ( 1.0d0 - h )

  if ( x <= 0.0d0 ) then
    go to 340
  end if

  if ( abs ( w ) >= 1.0d0 .and. abs ( w ) * t * t <= eps ) then
    return
  end if

  d = abs ( h )

  210 continue

  xn = x

  if ( d <= tol ) then

    if ( d <= eps ) then
      return
    end if

    if ( abs ( p - pn ) <= tol * p ) then
      return
    end if

  end if

  go to 180
!
!  Schroder iteration using Q.
!
  220 continue

  if ( q <= 1.0d10 * xmin ) then
    go to 350
  end if

  am1 = ( a - 0.5d0 ) - 0.5d0

  230 continue

  if ( a > amax ) then
    d = 0.5d0 + ( 0.5d0 - xn / a )
    if ( abs ( d ) <= e2 ) then
      go to 350
    end if
  end if

  if ( ierr >= 20 ) then
    ierr = -6
    return
  end if

  ierr = ierr + 1
  call gamma_inc ( a, xn, pn, qn, 0 )

  if ( pn == 0.0d0 .or. qn == 0.0d0 ) then
    go to 350
  end if

  r = rcomp ( a, xn )

  if ( r == 0.0d0 ) then
    go to 350
  end if

  t = ( q - qn ) / r
  w = 0.5d0 * ( am1 - xn )

  if ( abs ( t ) <= 0.1d0 .and. abs ( w * t ) <= 0.1d0 ) then
    go to 250
  end if

  x = xn * ( 1.0d0 - t )

  if ( x <= 0.0d0 ) then
    go to 340
  end if

  d = abs ( t )
  go to 260

  250 continue

  h = t * ( 1.0d0 + w * t )
  x = xn * ( 1.0d0 - h )

  if ( x <= 0.0d0 ) then
    ierr = - 7
    return
  end if

  if ( abs ( w ) >= 1.0d0 .and. abs ( w ) * t * t <= eps ) then
    return
  end if

  d = abs ( h )

  260 continue

  xn = x

  if ( d > tol ) then
    go to 230
  end if

  if ( d <= eps ) then
    return
  end if

  if ( abs ( q - qn ) <= tol * q ) then
    return
  end if

  go to 230
!
!  SPECIAL CASES
!
  340 ierr = -7
  return
!
  350 x = xn
  ierr = -8
  return
!
  360 x = xmax
  ierr = -8

  return
end
function gamma_log ( a )
!
!*******************************************************************************
!
!! GAMMA_LOG evaluates LN(GAMMA(A)) for positive A.
!
!
!  Reference:
!
!    A R DiDinato and A H Morris,
!    Algorithm 708: 
!    Significant Digit Computation of the Incomplete Beta Function Ratios,
!    ACM Transactions on Mathematical Software,
!    Volume 18, 1993, pages 360-373.
!
!  Author:
!
!    Alfred H Morris, Jr,
!    Naval Surface Weapons Center,
!    Dahlgren, Virginia.
!
  double precision, parameter :: c0 =   0.833333333333333d-01
  double precision, parameter :: c1 = - 0.277777777760991d-02
  double precision, parameter :: c2 =   0.793650666825390d-03
  double precision, parameter :: c3 = - 0.595202931351870d-03
  double precision, parameter :: c4 =   0.837308034031215d-03
  double precision, parameter :: c5 = - 0.165322962780713d-02
  double precision, parameter :: d  =   0.418938533204673d0
!
  double precision a
  double precision gamma_log
  double precision gamma_ln1
  integer i
  integer n
  double precision t
  double precision w
!
  if ( a <= 0.8d0 ) then

    gamma_log = gamma_ln1 ( a ) - log ( a )

  else if ( a <= 2.25d0 ) then

    t = ( a - 0.5d0 ) - 0.5d0
    gamma_log = gamma_ln1 ( t )

  else if ( a < 10.0d0 ) then

    n = a - 1.25d0
    t = a
    w = 1.0d0
    do i = 1, n
      t = t - 1.0d0
      w = t * w
    end do

    gamma_log = gamma_ln1 ( t - 1.0d0 ) + log ( w )

  else

    t = ( 1.0d0 / a )**2

    w = ((((( c5 * t + c4 ) * t + c3 ) * t + c2 ) * t + c1 ) * t + c0 ) / a

    gamma_log = ( d + w ) + ( a - 0.5d0 ) * ( log ( a ) - 1.0d0 )

  end if

  return
end
function gamma_ln1 ( a )
!
!*******************************************************************************
!
!! GAMMA_LN1 evaluates LN(GAMMA(1 + A)) FOR -0.2 .LE. A .LE. 1.25
!
!
  double precision, parameter :: p0 =   0.577215664901533d+00
  double precision, parameter :: p1 =   0.844203922187225d+00
  double precision, parameter :: p2 = - 0.168860593646662d+00
  double precision, parameter :: p3 = - 0.780427615533591d+00
  double precision, parameter :: p4 = - 0.402055799310489d+00
  double precision, parameter :: p5 = - 0.673562214325671d-01
  double precision, parameter :: p6 = - 0.271935708322958d-02
  double precision, parameter :: q1 =   0.288743195473681d+01
  double precision, parameter :: q2 =   0.312755088914843d+01
  double precision, parameter :: q3 =   0.156875193295039d+01
  double precision, parameter :: q4 =   0.361951990101499d+00
  double precision, parameter :: q5 =   0.325038868253937d-01
  double precision, parameter :: q6 =   0.667465618796164d-03
!
  double precision a
  double precision bot
  double precision gamma_ln1
  double precision, parameter :: r0 = 0.422784335098467d+00
  double precision, parameter :: r1 = 0.848044614534529d+00
  double precision, parameter :: r2 = 0.565221050691933d+00
  double precision, parameter :: r3 = 0.156513060486551d+00
  double precision, parameter :: r4 = 0.170502484022650d-01
  double precision, parameter :: r5 = 0.497958207639485d-03
  double precision, parameter :: s1 = 0.124313399877507d+01
  double precision, parameter :: s2 = 0.548042109832463d+00
  double precision, parameter :: s3 = 0.101552187439830d+00
  double precision, parameter :: s4 = 0.713309612391000d-02
  double precision, parameter :: s5 = 0.116165475989616d-03
  double precision top
  double precision w
  double precision x
!
  if ( a < 0.6d0 ) then

    top = (((((( p6 * a + p5 ) * a + p4 ) * a + p3 ) * a + p2 ) * a &
      + p1 ) * a + p0 ) 

    bot = (((((( q6 * a + q5 ) * a + q4 ) * a + q3 ) * a + q2 ) * a &
      + q1 ) * a + 1.0d0 )

    gamma_ln1 = - a * ( top / bot )

  else

    x = ( a - 0.5d0 ) - 0.5d0

    top = ((((( r5 * x + r4 ) * x + r3 ) * x + r2 ) * x + r1 ) * x + r0 ) 

    bot = ((((( s5 * x + s4 ) * x + s3 ) * x + s2 ) * x + s1 ) * x + 1.0d0 )

    gamma_ln1 = x * ( top / bot )

  end if

  return
end
function gamma ( a )
!
!*******************************************************************************
!
!! GAMMA evaluates the gamma function.
!
!
!     GAMMA(A) IS ASSIGNED THE VALUE 0 WHEN THE GAMMA function CANNOT
!     BE COMPUTED.
!
!  Author:
!
!    Alfred H Morris, Jr,
!    Naval Surface Weapons Center,
!    Dahlgren, Virginia.
!
  double precision, parameter :: d = 0.41893853320467274178d0
  double precision, parameter :: pi = 3.1415926535898d0
!
  double precision a
  double precision bot
  double precision dpmpar
  double precision exparg
  double precision g
  double precision gamma
  integer i
  integer j
  double precision lnx
  integer m
  integer n
  double precision p(7)
  double precision q(7)
  double precision, parameter :: r1 = 0.820756370353826d-03
  double precision, parameter :: r2 = -0.595156336428591d-03
  double precision, parameter :: r3 = 0.793650663183693d-03
  double precision, parameter :: r4 = -0.277777777770481d-02
  double precision, parameter :: r5 = 0.833333333333333d-01
  double precision s
  double precision t
  double precision top
  double precision w
  double precision x
  double precision z
!
  data p(1)/.539637273585445d-03/,p(2)/.261939260042690d-02/
  data p(3)/.204493667594920d-01/,p(4)/.730981088720487d-01/
  data p(5)/.279648642639792d+00/,p(6)/.553413866010467d+00/
  data p(7)/1.0d0/
  data q(1)/-.832979206704073d-03/,q(2)/.470059485860584d-02/
  data q(3)/.225211131035340d-01/,q(4)/-.170458969313360d+00/
  data q(5)/-.567902761974940d-01/,q(6)/.113062953091122d+01/
  data q(7)/1.0d0/
!
  gamma = 0.0d0
  x = a

  if ( abs ( a ) < 15.0d0 ) then
!
!  Evaluation of GAMMA(A) for ABS(A) .LT. 15
!
    t = 1.0d0
    m = int ( a ) - 1
!
!  Let T be the product of A-J when A .GE. 2
!
    if ( m >= 0 ) then

      do j = 1, m
        x = x - 1.0d0
        t = x * t
      end do

      x = x - 1.0d0
!
!  Let T be the product of A+J WHEN A .LT. 1
!
    else

      t = a

      if ( a <= 0.0d0 ) then

        m = - m - 1

        do j = 1, m
          x = x + 1.0d0
          t = x * t
        end do

        x = ( x + 0.5d0 ) + 0.5d0
        t = x * t
        if ( t == 0.0d0 ) then
          return
        end if

      end if
!
!  THE FOLLOWING CODE CHECKS IF 1/T CAN OVERFLOW. THIS
!  CODE MAY BE OMITTED IF DESIRED.
!
      if ( abs ( t ) < 1.0d-30 ) then
        if ( abs ( t ) * dpmpar(3) > 1.0001d0 ) then
          gamma = 1.0d0 / t
        end if
        return
      end if

    end if
!
!  COMPUTE GAMMA(1 + X) FOR  0 .LE. X .LT. 1
!
    top = p(1)
    bot = q(1)
    do i = 2, 7
      top = top * x + p(i)
      bot = bot * x + q(i)
    end do

    gamma = top / bot
!
!  Termination.
!
    if ( a >= 1.0d0 ) then
      gamma = gamma * t
    else
      gamma = gamma / t
    end if
!
!  Evaluation of GAMMA(A) FOR ABS(A) .GE. 15
!
  else

    if ( abs ( a ) >= 1000.0 ) then
      return
    end if

    if ( a <= 0.0d0 ) then

      x = - a
      n = x
      t = x - n

      if ( t > 0.9d0 ) then
        t = 1.0d0 - t
      end if

      s = sin ( pi * t ) / pi

      if ( mod ( n, 2 ) == 0 ) then
        s = - s
      end if

      if ( s == 0.0d0 ) then
        return
      end if

    end if
!
!  Compute the modified asymptotic sum.
!
    t = 1.0d0 / ( x * x )

    g = (((( r1 * t + r2 ) * t + r3 ) * t + r4 ) * t + r5 ) / x

    lnx = log ( x )
!
!  Final assembly.
!
    z = x
    g = ( d + g ) + ( z - 0.5d0 ) * ( lnx - 1.d0 )
    w = g
    t = g - dble ( w )
  
    if ( w > 0.99999d0 * exparg ( 0 ) ) then
      return
    end if

    gamma = exp ( w )* ( 1.0d0 + t )

    if ( a < 0.0d0 ) then
      gamma = ( 1.0d0 / ( gamma * s ) ) / x
    end if

  end if

  return
end
subroutine gamma_inc ( a, x, ans, qans, ind )
!
!*******************************************************************************
!
!! GAMMA_INC evaluates the incomplete gamma ratio functions P(A,X) and Q(A,X).
!
!
!  Author:
!
!    Alfred H Morris, Jr,
!    Naval Surface Weapons Center,
!    Dahlgren, Virginia.
!
!  Parameters:
!
!    Input, double precision A, X, the arguments of the incomplete
!    gamma ratio.  A and X must be nonnegative.  A and X cannot
!    both be zero.
!
!    Output, double precision ANS, QANS.  On normal output,
!    ANS = P(A,X) and QANS = Q(A,X).  However, ANS is set to 2 if
!    A or X is negative, or both are 0, or when the answer is
!    computationally indeterminate because A is extremely large
!    and X is very close to A.
!
!    Input, integer IND, indicates the accuracy request:
!    0, as much accuracy as possible.
!    1, to within 1 unit of the 6-th significant digit, 
!    otherwise, to within 1 unit of the 3rd significant digit.
!
  double precision, parameter :: alog10 = 2.30258509299405d0
!
  double precision a
  double precision a2n
  double precision a2nm1
  double precision acc
  double precision acc0(3)
  double precision am0
  double precision amn
  double precision an
  double precision an0
  double precision ans
  double precision apn
  double precision b2n
  double precision b2nm1
  double precision big(3)
  double precision c
  double precision c0
  double precision c1
  double precision c2
  double precision c3
  double precision c4
  double precision c5
  double precision c6
  double precision cma
  double precision d0(13)
  double precision d1(12)
  double precision d2(10)
  double precision d3(8)
  double precision d4(6)
  double precision d5(4)
  double precision d6(2)
  double precision d10
  double precision d20
  double precision d30
  double precision d40
  double precision d50
  double precision d60
  double precision d70
  double precision dpmpar
  double precision e
  double precision e0
  double precision e00(3)
  double precision erf
  double precision erfc1
  double precision g
  double precision gam1
  double precision gamma
  double precision h
  integer i
  integer ind
  integer iop
  double precision j
  double precision l
  integer m
  integer n
  integer n_max
  double precision qans
  double precision r
  double precision rexp
  double precision rlog
  double precision, parameter :: rt2pin = 0.398942280401433d0
  double precision rta
  double precision, parameter :: rtpi = 1.77245385090552d0
  double precision rtx
  double precision s
  double precision sum
  double precision t
  double precision t1
  double precision tol
  double precision twoa
  double precision u
  double precision w
  double precision wk(20)
  double precision x
  double precision x0
  double precision x00(3)
  double precision y
  double precision z
!
!     ALOG10 = LN(10)
!     RT2PIN = 1/SQRT(2*PI)
!     RTPI   = SQRT(PI)
!
  data acc0(1)/5.d-15/,acc0(2)/5.d-7/,acc0(3)/5.d-4/
  data big(1)/20.0d0/,big(2)/14.0d0/,big(3)/10.0d0/
  data e00(1)/.25d-3/,e00(2)/.25d-1/,e00(3)/.14d0/
  data x00(1)/31.0d0/,x00(2)/17.0d0/,x00(3)/9.7d0/
  data d0(1)/.833333333333333d-01/
  data d0(2)/-.148148148148148d-01/
  data d0(3)/.115740740740741d-02/,d0(4)/.352733686067019d-03/
  data d0(5)/-.178755144032922d-03/,d0(6)/.391926317852244d-04/
  data d0(7)/-.218544851067999d-05/,d0(8)/-.185406221071516d-05/
  data d0(9)/.829671134095309d-06/,d0(10)/-.176659527368261d-06/
  data d0(11)/.670785354340150d-08/,d0(12)/.102618097842403d-07/
  data d0(13)/-.438203601845335d-08/
  data d10/-.185185185185185d-02/,d1(1)/-.347222222222222d-02/
  data d1(2)/.264550264550265d-02/,d1(3)/-.990226337448560d-03/
  data d1(4)/.205761316872428d-03/,d1(5)/-.401877572016461d-06/
  data d1(6)/-.180985503344900d-04/,d1(7)/.764916091608111d-05/
  data d1(8)/-.161209008945634d-05/,d1(9)/.464712780280743d-08/
  data d1(10)/.137863344691572d-06/,d1(11)/-.575254560351770d-07/
  data d1(12)/.119516285997781d-07/
  data d20/.413359788359788d-02/,d2(1)/-.268132716049383d-02/
  data d2(2)/.771604938271605d-03/,d2(3)/.200938786008230d-05/
  data d2(4)/-.107366532263652d-03/,d2(5)/.529234488291201d-04/
  data d2(6)/-.127606351886187d-04/,d2(7)/.342357873409614d-07/
  data d2(8)/.137219573090629d-05/,d2(9)/-.629899213838006d-06/
  data d2(10)/.142806142060642d-06/
  data d30/.649434156378601d-03/,d3(1)/.229472093621399d-03/
  data d3(2)/-.469189494395256d-03/,d3(3)/.267720632062839d-03/
  data d3(4)/-.756180167188398d-04/,d3(5)/-.239650511386730d-06/
  data d3(6)/.110826541153473d-04/,d3(7)/-.567495282699160d-05/
  data d3(8)/.142309007324359d-05/
  data d40/-.861888290916712d-03/,d4(1)/.784039221720067d-03/
  data d4(2)/-.299072480303190d-03/,d4(3)/-.146384525788434d-05/
  data d4(4)/.664149821546512d-04/,d4(5)/-.396836504717943d-04/
  data d4(6)/.113757269706784d-04/
  data d50/-.336798553366358d-03/,d5(1)/-.697281375836586d-04/
  data d5(2)/.277275324495939d-03/,d5(3)/-.199325705161888d-03/
  data d5(4)/.679778047793721d-04/
  data d60/.531307936463992d-03/,d6(1)/-.592166437353694d-03/
  data d6(2)/.270878209671804d-03/
  data d70 / .344367606892378d-03/
!
!  E IS A MACHINE DEPendENT CONSTANT. E IS THE SMALLEST
!  FLOATING POINT NUMBER FOR WHICH 1.0 + E .GT. 1.0 .
!
  e = epsilon ( 1.0d0 )

  if ( a < 0.0d0 .or. x < 0.0d0 ) go to 430
  if ( a == 0.0d0 .and. x == 0.0d0 ) go to 430
  if ( a * x == 0.0d0 ) go to 420

  iop = ind + 1
  if ( iop /= 1 .and. iop /= 2 ) iop = 3
  acc = max ( acc0(iop), e )
  e0 = e00(iop)
  x0 = x00(iop)
!
!  SELECT THE APPROPRIATE ALGORITHM
!
  if ( a >= 1.0d0 ) go to 10
  if ( a == 0.5d0 ) go to 390
  if ( x < 1.1d0 ) go to 160
  t1 = a * log ( x ) - x
  u = a * exp ( t1 )
  if ( u == 0.0d0 ) go to 380
  r = u * ( 1.0d0 + gam1 ( a ) )
  go to 250

   10 continue

  if ( a >= big(iop) ) go to 30

  if ( a > x .or. x >= x0 ) go to 20

  twoa = a + a
  m = int ( twoa )

  if ( twoa == dble ( m ) ) then
    i = m / 2
    if ( a == dble ( i ) ) go to 210
    go to 220
  end if

   20 continue

  t1 = a * log ( x ) - x
  r = exp ( t1 ) / gamma ( a )
  go to 40

   30 continue

  l = x / a
  if ( l == 0.0d0 ) go to 370
  s = 0.5d0 + ( 0.5d0 - l )
  z = rlog ( l )
  if ( z >= 700.0d0 / a ) go to 410
  y = a * z
  rta = sqrt ( a )
  if ( abs ( s ) <= e0 / rta ) go to 330
  if ( abs ( s ) <= 0.4d0 ) go to 270

  t = ( 1.0d0 / a )**2
  t1 = (((0.75d0*t-1.0d0)*t+3.5d0)*t-105.0d0)/ (a*1260.0d0)
  t1 = t1 - y
  r = rt2pin * rta * exp ( t1 )

40    continue

  if (r==0.0d0) go to 420

  if ( x <= max ( a, alog10 ) ) go to 50
  if ( x < x0 ) go to 250
  go to 100
!
!  TAYLOR SERIES FOR P/R.
!
50    continue

  apn = a + 1.0d0
  t = x / apn
  wk(1) = t
  do n = 2, 20
    apn = apn + 1.0d0
    t = t * ( x / apn )
    if ( t <= 1.0d-3 ) then
      go to 70
    end if
    wk(n) = t
  end do

  n = 20

70    continue

  sum = t

  tol = 0.5d0 * acc

80    continue

  apn = apn + 1.0d0
  t = t * ( x / apn )
  sum = sum + t
  if ( t > tol ) then
    go to 80
  end if

  n_max = n - 1
  do m = 1, n_max
    n = n - 1
    sum = sum + wk(n)
  end do

  ans = ( r / a ) * ( 1.0d0 + sum )
  qans = 0.5d0 + ( 0.5d0 - ans )
  return
!
!  Asymptotic expansion.
!
  100 continue

  amn = a - 1.0d0
  t = amn / x
  wk(1) = t
  do n = 2, 20
    amn = amn - 1.0d0
    t = t * ( amn / x )
    if ( abs ( t ) <= 1.0d-3 ) then
      go to 120
    end if
    wk(n) = t
  end do

  n = 20

  120 continue

  sum = t

  130 continue

  if ( abs ( t ) > acc ) then
    amn = amn - 1.0d0
    t = t * ( amn / x )
    sum = sum + t
    go to 130
  end if

  140 continue

  n_max = n - 1
  do m = 1, n_max
    n = n - 1
    sum = sum + wk(n)
  end do
  qans = ( r / x ) * ( 1.0d0 + sum )
  ans = 0.5d0 + ( 0.5d0 - qans )
  return
!
!  Taylor series for P(A,X)/X**A
!
  160 continue

  an = 3.0d0
  c = x
  sum = x / ( a + 3.0d0 )
  tol = 3.0d0 * acc / ( a + 1.0d0 )

  170 continue

  an = an + 1.0d0
  c = - c * ( x / an )
  t = c / ( a + an )
  sum = sum + t
  if ( abs ( t ) > tol ) then
    go to 170
  end if

  j = a * x * (( sum / 6.0d0 - 0.5d0 / ( a + 2.0d0 )) * x + 1.0d0 &
    / ( a + 1.0d0 ))

  z = a * log ( x )
  h = gam1 ( a )
  g = 1.0d0 + h

  if ( x < 0.25d0 ) go to 180
  if ( a < x / 2.59d0 ) go to 200
  go to 190

  180 continue

  if ( z > -0.13394d0 ) go to 200

  190 continue

  w = exp ( z )
  ans = w * g * ( 0.5d0 + ( 0.5d0 - j ))
  qans = 0.5d0 + ( 0.5d0 - ans )
  return

200   continue

  l = rexp ( z )
  w = 0.5d0 + ( 0.5d0 + l )
  qans = ( w * j - l ) * g - h
  if ( qans < 0.0d0 ) then
    go to 380
  end if
  ans = 0.5d0 + ( 0.5d0 - qans )
  return
!
!  Finite sums for Q when A .GE. 1 AND 2*A IS AN integer
!
210   continue

  sum = exp ( - x )
  t = sum
  n = 1
  c = 0.0d0
  go to 230

220   continue

  rtx = sqrt ( x )
  sum = erfc1 ( 0, rtx )
  t = exp ( - x ) / ( rtpi * rtx )
  n = 0
  c = - 0.5d0

  230 continue

  if ( n /= i ) then
    n = n + 1
    c = c + 1.0d0
    t = ( x * t ) / c
    sum = sum + t
    go to 230
  end if

  240 continue

  qans = sum
  ans = 0.5d0 + ( 0.5d0 - qans )
  return
!
!  Continued fraction expansion.
!
250   continue

  tol = max ( 5.0d0 * e, acc )
  a2nm1 = 1.0d0
  a2n = 1.0d0
  b2nm1 = x
  b2n = x + ( 1.0d0 - a )
  c = 1.0d0

260   continue

  a2nm1 = x * a2n + c * a2nm1
  b2nm1 = x * b2n + c * b2nm1
  am0 = a2nm1 / b2nm1
  c = c + 1.0d0
  cma = c - a
  a2n = a2nm1 + cma * a2n
  b2n = b2nm1 + cma * b2n
  an0 = a2n / b2n

  if ( abs ( an0 - am0 ) >= tol * an0 ) then
    go to 260
  end if

  qans = r * an0
  ans = 0.5d0 + ( 0.5d0 - qans )
  return
!
!  General Temme expansion.
!
270   continue

  if ( abs ( s ) <= 2.0d0 * e .and. a*e*e>3.28d-3) go to 430

  c = exp ( - y )
  w = 0.5d0 * erfc1 ( 1, sqrt ( y ) )
  u = 1.0d0 / a
  z = sqrt ( z + z )
  if ( l < 1.0d0 ) then
    z = - z
  end if

  if ( iop < 2 ) then

    if ( abs ( s ) <= 1.0d-3 ) go to 340

    c0 = (((((((((((( d0(13)*z + d0(12))*z + d0(11))*z + d0(10))*z &
      + d0(9))*z + d0(8))*z + d0(7))*z + d0(6))*z + d0(5))*z &
      + d0(4))*z + d0(3))*z + d0(2))*z + d0(1))*z - 1.0 / 3.0

    c1 = ((((((((((( d1(12)*z + d1(11))*z + d1(10))*z + d1(9))*z &
      + d1(8))*z + d1(7))*z + d1(6))*z + d1(5))*z + d1(4))*z &
      + d1(3))*z + d1(2))*z + d1(1))*z + d10

    c2 = ((((((((( d2(10)*z + d2(9))*z + d2(8))*z + d2(7))*z &
       + d2(6))*z + d2(5))*z + d2(4))*z + d2(3))*z + d2(2))*z + d2(1))*z &
       + d20

    c3 = ((((((( d3(8)*z+d3(7))*z+d3(6))*z+d3(5))*z+d3(4))*z+d3(3))*z+ &
      d3(2))*z+d3(1))*z + d30

    c4 = ((((( d4(6)*z+d4(5))*z+d4(4))*z+d4(3))*z+d4(2))*z+d4(1))*z + d40

    c5 = (((d5(4)*z+d5(3))*z+d5(2))*z+d5(1))*z + d50

    c6 = (d6(2)*z+d6(1))*z + d60

    t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u + c0

  else if ( iop == 2 ) then

    c0 = (((((d0(6)*z+d0(5))*z+d0(4))*z +d0(3))*z+d0(2))*z+d0(1))*z - 1.0 / 3.0
    c1 = (((d1(4)*z+d1(3))*z+d1(2))*z+d1(1))*z + d10
    c2 = d2(1) * z + d20
    t = ( c2 * u + c1 ) * u + c0

  else if ( iop > 2 ) then

    t = (( d0(3) * z + d0(2) ) * z + d0(1) ) * z - 1.0 / 3.0

  end if

310   continue

  if ( l >= 1.0d0 ) then
    qans = c * ( w + rt2pin * t / rta )
    ans = 0.5d0 + ( 0.5d0 - qans )
  else
    ans = c * ( w - rt2pin * t / rta )
    qans = 0.5d0 + ( 0.5d0 - ans )
  end if

  return
!
!  Temme expansion for L = 1
!
  330 continue

  if (a*e*e>3.28d-3) go to 430
  c = 0.5d0 + ( 0.5d0 - y )
  w = (0.5d0-sqrt(y)* (0.5d0+ (0.5d0-y/3.0d0))/rtpi)/c
  u = 1.0d0 / a
  z = sqrt( z + z )

  if ( l < 1.0d0 ) then
    z = -z
  end if

  if (iop-2) 340,350,360

  340 continue

  c0 = ((((((d0(7)*z+d0(6))*z+d0(5))*z+d0(4))*z +d0(3))*z+d0(2))*z+ &
     d0(1))*z - 1.0 / 3.0

  c1 = (((((d1(6)*z+d1(5))*z+d1(4))*z+d1(3))*z+d1(2))*z+d1(1))*z + d10

  c2 = ((((d2(5)*z+d2(4))*z+d2(3))*z+d2(2))*z+d2(1))*z + d20

  c3 = (((d3(4)*z+d3(3))*z+d3(2))*z+d3(1))*z + d30

  c4 = (d4(2)*z+d4(1))*z + d40
  c5 = (d5(2)*z+d5(1))*z + d50
  c6 = d6(1)*z + d60
  t = ((((((d70*u+c6)*u+c5)*u+c4)*u+c3)*u+c2)*u+c1)*u + c0
  go to 310

  350 continue

  c0 = ( d0(2) * z + d0(1) ) * z - 1.0 / 3.0
  c1 = d1(1) * z + d10
  t = ( d20 * u + c1 ) * u + c0
  go to 310

  360 continue

  t = d0(1) * z - 1.0 / 3.0
  go to 310
!
!  Special cases
!
  370 continue

  ans = 0.0d0
  qans = 1.0d0
  return

  380 continue

  ans = 1.0d0
  qans = 0.0d0
  return

  390 continue

  if ( x < 0.25d0 ) then
    ans = erf ( sqrt ( x ) )
    qans = 0.5d0 + ( 0.5d0 - ans )
  else
    qans = erfc1 ( 0, sqrt ( x ) )
    ans = 0.5d0 + ( 0.5d0 - qans )
  end if

  return

  410 continue

  if ( abs ( s ) <= 2.0d0 * e ) go to 430

  420 continue

  if ( x <= a ) go to 370
  go to 380
!
!  Error return
!
  430 ans = 2.0d0

  return
end
subroutine gamma_rat1 ( a, x, r, p, q, eps )
!
!*******************************************************************************
!
!! GAMMA_RAT1 evaluates the incomplete gamma ratio functions P(A,X) and Q(A,X).
!
!
!     It is assumed that A .LE. 1.  EPS is the tolerance to be used.
!     The input argument R has the value E**(-X)*X**A/GAMMA(A).
!
  double precision a
  double precision a2n
  double precision a2nm1
  double precision am0
  double precision an
  double precision an0
  double precision b2n
  double precision b2nm1
  double precision c
  double precision cma
  double precision eps
  double precision erf
  double precision erfc1
  double precision g
  double precision gam1
  double precision h
  double precision j
  double precision l
  double precision p
  double precision q
  double precision r
  double precision rexp
  double precision sum
  double precision t
  double precision tol
  double precision w
  double precision x
  double precision z
!
  if ( a * x == 0.0d0 ) then

    if ( x <= a ) then
      p = 0.0d0
      q = 1.0d0
    else
      p = 1.0d0
      q = 0.0d0
    end if

    return
  end if

  if ( a == 0.5d0 ) then

    if ( x < 0.25d0 ) then
      p = erf ( sqrt ( x ) )
      q = 0.5d0 + ( 0.5d0 - p )
    else
      q = erfc1 ( 0, sqrt ( x ) )
      p = 0.5d0 + ( 0.5d0 - q )
    end if

    return

  end if

  if (x<1.1d0) go to 10
  go to 60
!
!  Taylor series for P(A,X)/X**A
!
   10 continue

  an = 3.0d0
  c = x
  sum = x / ( a + 3.0d0 )
  tol = 0.1d0 * eps / ( a + 1.0d0 )

   20 continue

  an = an + 1.0d0
  c = - c * ( x / an )
  t = c / ( a + an )
  sum = sum + t
  if ( abs ( t ) > tol ) then
    go to 20
  end if

  j = a * x * ((sum/6.0d0-0.5d0/ (a+2.0d0))*x+1.0d0/ (a+1.0d0))

  z = a * log ( x )
  h = gam1 ( a )
  g = 1.0d0 + h

  if ( x < 0.25d0 ) go to 30

  if ( a < x / 2.59d0 ) then
    go to 50
  else
    go to 40
  end if

   30 continue

  if ( z > - 0.13394d0 ) go to 50

   40 continue

  w = exp ( z )
  p = w * g * ( 0.5d0 + ( 0.5d0 - j ))
  q = 0.5d0 + ( 0.5d0 - p )
  return

   50 continue

  l = rexp ( z )
  w = 0.5d0 + ( 0.5d0 + l )
  q = ( w * j - l ) * g - h

  if ( q < 0.0d0 ) then
    p = 1.0
    q = 0.0
  else
    p = 0.5d0 + ( 0.5d0 - q )
  end if

  return
!
!  Continued fraction expansion.
!
   60 continue

  a2nm1 = 1.0d0
  a2n = 1.0d0
  b2nm1 = x
  b2n = x + ( 1.0d0 - a )
  c = 1.0d0

   70 continue

  a2nm1 = x * a2n + c * a2nm1
  b2nm1 = x * b2n + c * b2nm1
  am0 = a2nm1 / b2nm1
  c = c + 1.0d0
  cma = c - a
  a2n = a2nm1 + cma * a2n
  b2n = b2nm1 + cma * b2n
  an0 = a2n / b2n
  if ( abs ( an0 - am0 ) >= eps * an0 ) then
    go to 70
  end if

  q = r * an0
  p = 0.5d0 + ( 0.5d0 - q )

  return
end
function gsumln ( a, b )
!
!*******************************************************************************
!
!! GSUMLN evaluates the function LN(GAMMA(A + B)).
!
!
!  GSUMLN is used FOR 1 .LE. A .LE. 2  AND  1 .LE. B .LE. 2
!
!
  double precision a
  double precision alnrel
  double precision b
  double precision gamma_ln1
  double precision gsumln
  double precision x
!
  x = dble ( a ) + dble ( b ) - 2.d0

  if ( x <= 0.25d0 ) then
    gsumln = gamma_ln1 ( 1.0d0 + x )
  else if ( x <= 1.25d0 ) then
    gsumln = gamma_ln1 ( x ) + alnrel ( x )
  else
    gsumln = gamma_ln1 ( x - 1.0d0 ) + log ( x * ( 1.0d0 + x ))
  end if

  return
end
function ipmpar ( i )
!
!*******************************************************************************
!
!! IPMPAR returns integer machine constants. 
!
!
!  It is assumed that the argument I IS AN integer
!     HAVING ONE OF THE VALUES 1-10. IPMPAR(I) HAS THE VALUE:
!
!  integerS.
!
!     Assume integers are represented in the N-digit, base-A form
!
!               SIGN ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )
!
!               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,N-1.
!
!     IPMPAR(1) = A, THE BASE.
!
!     IPMPAR(2) = N, THE NUMBER OF BASE-A DIGITS.
!
!     IPMPAR(3) = A**N - 1, THE LARGEST MAGNITUDE.
!
!  FLOATING-POINT NUMBERS.
!
!     IT IS ASSUMED THAT THE SINGLE AND double precision FLOATING
!     POINT ARITHMETICS HAVE THE SAME BASE, SAY B, AND THAT THE
!     NONZERO NUMBERS ARE REPRESENTED IN THE FORM
!
!               SIGN (B**E) * (X(1)/B + ... + X(M)/B**M)
!
!               WHERE X(I) = 0,1,...,B-1 FOR I=1,...,M,
!               X(1) .GE. 1, AND EMIN .LE. E .LE. EMAX.
!
!     IPMPAR(4) = B, THE BASE.
!
!  SINGLE-PRECISION
!
!     IPMPAR(5) = M, THE NUMBER OF BASE-B DIGITS.
!
!     IPMPAR(6) = EMIN, THE SMALLEST EXPONENT E.
!
!     IPMPAR(7) = EMAX, THE LARGEST EXPONENT E.
!
!  DOUBLE-PRECISION
!
!     IPMPAR(8) = M, THE NUMBER OF BASE-B DIGITS.
!
!     IPMPAR(9) = EMIN, THE SMALLEST EXPONENT E.
!
!     IPMPAR(10) = EMAX, THE LARGEST EXPONENT E.
!
!     TO DEFINE THIS function FOR THE COMPUTER BEING USED, ACTIVATE
!     THE data STATMENTS FOR THE COMPUTER BY REMOVING THE C FROM
!     COLUMN 1. (ALL THE OTHER data STATEMENTS SHOULD HAVE C IN
!     COLUMN 1.)
!
!     IPMPAR IS AN ADAPTATION OF THE function I1MACH, WRITTEN BY
!     P.A. FOX, A.D. HALL, AND N.L. SCHRYER (BELL LABORATORIES).
!     IPMPAR WAS FORMED BY A.H. MORRIS (NSWC). THE CONSTANTS ARE
!     FROM BELL LABORATORIES, NSWC, AND OTHER SOURCES.
!
  integer i
  integer imach(10)
  integer ipmpar
!
!     MACHINE CONSTANTS FOR AMDAHL MACHINES.
!
!     data imach( 1) /   2 /
!     data imach( 2) /  31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /  16 /
!     data imach( 5) /   6 /
!     data imach( 6) / -64 /
!     data imach( 7) /  63 /
!     data imach( 8) /  14 /
!     data imach( 9) / -64 /
!     data imach(10) /  63 /
!
!     Machine constants for the AT&T 3B SERIES, AT&T
!     PC 7300, AND AT&T 6300.
!
!     data imach( 1) /     2 /
!     data imach( 2) /    31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /     2 /
!     data imach( 5) /    24 /
!     data imach( 6) /  -125 /
!     data imach( 7) /   128 /
!     data imach( 8) /    53 /
!     data imach( 9) / -1021 /
!     data imach(10) /  1024 /
!
!     Machine constants for the BURROUGHS 1700 SYSTEM.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   33 /
!     data imach( 3) / 8589934591 /
!     data imach( 4) /    2 /
!     data imach( 5) /   24 /
!     data imach( 6) / -256 /
!     data imach( 7) /  255 /
!     data imach( 8) /   60 /
!     data imach( 9) / -256 /
!     data imach(10) /  255 /
!
!     Machine constants for the BURROUGHS 5700 SYSTEM.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   39 /
!     data imach( 3) / 549755813887 /
!     data imach( 4) /    8 /
!     data imach( 5) /   13 /
!     data imach( 6) /  -50 /
!     data imach( 7) /   76 /
!     data imach( 8) /   26 /
!     data imach( 9) /  -50 /
!     data imach(10) /   76 /
!
!     Machine constants for the BURROUGHS 6700/7700 SYSTEMS.
!
!     data imach( 1) /      2 /
!     data imach( 2) /     39 /
!     data imach( 3) / 549755813887 /
!     data imach( 4) /      8 /
!     data imach( 5) /     13 /
!     data imach( 6) /    -50 /
!     data imach( 7) /     76 /
!     data imach( 8) /     26 /
!     data imach( 9) / -32754 /
!     data imach(10) /  32780 /
!
!     Machine constants for the CDC 6000/7000 SERIES
!     60 BIT ARITHMETIC, AND THE CDC CYBER 995 64 BIT
!     ARITHMETIC (NOS OPERATING SYSTEM).
!
!     data imach( 1) /    2 /
!     data imach( 2) /   48 /
!     data imach( 3) / 281474976710655 /
!     data imach( 4) /    2 /
!     data imach( 5) /   48 /
!     data imach( 6) / -974 /
!     data imach( 7) / 1070 /
!     data imach( 8) /   95 /
!     data imach( 9) / -926 /
!     data imach(10) / 1070 /
!
!     Machine constants for the CDC CYBER 995 64 BIT
!     ARITHMETIC (NOS/VE OPERATING SYSTEM).
!
!     data imach( 1) /     2 /
!     data imach( 2) /    63 /
!     data imach( 3) / 9223372036854775807 /
!     data imach( 4) /     2 /
!     data imach( 5) /    48 /
!     data imach( 6) / -4096 /
!     data imach( 7) /  4095 /
!     data imach( 8) /    96 /
!     data imach( 9) / -4096 /
!     data imach(10) /  4095 /
!
!     Machine constants for the CRAY 1, XMP, 2, AND 3.
!
!     data imach( 1) /     2 /
!     data imach( 2) /    63 /
!     data imach( 3) / 9223372036854775807 /
!     data imach( 4) /     2 /
!     data imach( 5) /    47 /
!     data imach( 6) / -8189 /
!     data imach( 7) /  8190 /
!     data imach( 8) /    94 /
!     data imach( 9) / -8099 /
!     data imach(10) /  8190 /
!
!     Machine constants for the data GENERAL ECLIPSE S/200.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   15 /
!     data imach( 3) / 32767 /
!     data imach( 4) /   16 /
!     data imach( 5) /    6 /
!     data imach( 6) /  -64 /
!     data imach( 7) /   63 /
!     data imach( 8) /   14 /
!     data imach( 9) /  -64 /
!     data imach(10) /   63 /
!
!     Machine constants for the HARRIS 220.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   23 /
!     data imach( 3) / 8388607 /
!     data imach( 4) /    2 /
!     data imach( 5) /   23 /
!     data imach( 6) / -127 /
!     data imach( 7) /  127 /
!     data imach( 8) /   38 /
!     data imach( 9) / -127 /
!     data imach(10) /  127 /
!
!     Machine constants for the HONEYWELL 600/6000
!     AND DPS 8/70 SERIES.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   35 /
!     data imach( 3) / 34359738367 /
!     data imach( 4) /    2 /
!     data imach( 5) /   27 /
!     data imach( 6) / -127 /
!     data imach( 7) /  127 /
!     data imach( 8) /   63 /
!     data imach( 9) / -127 /
!     data imach(10) /  127 /
!
!     Machine constants for the HP 2100
!     3 WORD double precision OPTION WITH FTN4
!
!     data imach( 1) /    2 /
!     data imach( 2) /   15 /
!     data imach( 3) / 32767 /
!     data imach( 4) /    2 /
!     data imach( 5) /   23 /
!     data imach( 6) / -128 /
!     data imach( 7) /  127 /
!     data imach( 8) /   39 /
!     data imach( 9) / -128 /
!     data imach(10) /  127 /
!
!     Machine constants for the HP 2100
!     4 WORD double precision OPTION WITH FTN4
!
!     data imach( 1) /    2 /
!     data imach( 2) /   15 /
!     data imach( 3) / 32767 /
!     data imach( 4) /    2 /
!     data imach( 5) /   23 /
!     data imach( 6) / -128 /
!     data imach( 7) /  127 /
!     data imach( 8) /   55 /
!     data imach( 9) / -128 /
!     data imach(10) /  127 /
!
!     Machine constants for the HP 9000.
!
!     data imach( 1) /     2 /
!     data imach( 2) /    31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /     2 /
!     data imach( 5) /    24 /
!     data imach( 6) /  -126 /
!     data imach( 7) /   128 /
!     data imach( 8) /    53 /
!     data imach( 9) / -1021 /
!     data imach(10) /  1024 /
!
!     Machine constants for the IBM 360/370 SERIES,
!     THE ICL 2900, THE ITEL AS/6, THE XEROX SIGMA
!     5/7/9 AND THE SEL SYSTEMS 85/86.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /   16 /
!     data imach( 5) /    6 /
!     data imach( 6) /  -64 /
!     data imach( 7) /   63 /
!     data imach( 8) /   14 /
!     data imach( 9) /  -64 /
!     data imach(10) /   63 /
!
!     Machine constants for the IBM PC.
!
!      data imach(1)/2/
!      data imach(2)/31/
!      data imach(3)/2147483647/
!      data imach(4)/2/
!      data imach(5)/24/
!      data imach(6)/-125/
!      data imach(7)/128/
!      data imach(8)/53/
!      data imach(9)/-1021/
!      data imach(10)/1024/
!
!     Machine constants for the MACINTOSH II - ABSOFT
!     MACFORTRAN II.
!
!     data imach( 1) /     2 /
!     data imach( 2) /    31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /     2 /
!     data imach( 5) /    24 /
!     data imach( 6) /  -125 /
!     data imach( 7) /   128 /
!     data imach( 8) /    53 /
!     data imach( 9) / -1021 /
!     data imach(10) /  1024 /
!
!     Machine constants for the MICROVAX - VMS FORTRAN.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /    2 /
!     data imach( 5) /   24 /
!     data imach( 6) / -127 /
!     data imach( 7) /  127 /
!     data imach( 8) /   56 /
!     data imach( 9) / -127 /
!     data imach(10) /  127 /
!
!     Machine constants for the PDP-10 (KA PROCESSOR).
!
!     data imach( 1) /    2 /
!     data imach( 2) /   35 /
!     data imach( 3) / 34359738367 /
!     data imach( 4) /    2 /
!     data imach( 5) /   27 /
!     data imach( 6) / -128 /
!     data imach( 7) /  127 /
!     data imach( 8) /   54 /
!     data imach( 9) / -101 /
!     data imach(10) /  127 /
!
!     Machine constants for the PDP-10 (KI PROCESSOR).
!
!     data imach( 1) /    2 /
!     data imach( 2) /   35 /
!     data imach( 3) / 34359738367 /
!     data imach( 4) /    2 /
!     data imach( 5) /   27 /
!     data imach( 6) / -128 /
!     data imach( 7) /  127 /
!     data imach( 8) /   62 /
!     data imach( 9) / -128 /
!     data imach(10) /  127 /
!
!     Machine constants for the PDP-11 FORTRAN SUPPORTING
!     32-BIT integer ARITHMETIC.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /    2 /
!     data imach( 5) /   24 /
!     data imach( 6) / -127 /
!     data imach( 7) /  127 /
!     data imach( 8) /   56 /
!     data imach( 9) / -127 /
!     data imach(10) /  127 /
!
!     Machine constants for the SEQUENT BALANCE 8000.
!
!     data imach( 1) /     2 /
!     data imach( 2) /    31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /     2 /
!     data imach( 5) /    24 /
!     data imach( 6) /  -125 /
!     data imach( 7) /   128 /
!     data imach( 8) /    53 /
!     data imach( 9) / -1021 /
!     data imach(10) /  1024 /
!
!     Machine constants for the SILICON GRAPHICS IRIS-4D
!     SERIES (MIPS R3000 PROCESSOR).
!
!     data imach( 1) /     2 /
!     data imach( 2) /    31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /     2 /
!     data imach( 5) /    24 /
!     data imach( 6) /  -125 /
!     data imach( 7) /   128 /
!     data imach( 8) /    53 /
!     data imach( 9) / -1021 /
!     data imach(10) /  1024 /
!
!     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
!     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
!     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
!
  data imach( 1) /     2 /
  data imach( 2) /    31 /
  data imach( 3) / 2147483647 /
  data imach( 4) /     2 /
  data imach( 5) /    24 /
  data imach( 6) /  -125 /
  data imach( 7) /   128 /
  data imach( 8) /    53 /
  data imach( 9) / -1021 /
  data imach(10) /  1024 /
!
!     Machine constants for the UNIVAC 1100 SERIES.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   35 /
!     data imach( 3) / 34359738367 /
!     data imach( 4) /    2 /
!     data imach( 5) /   27 /
!     data imach( 6) / -128 /
!     data imach( 7) /  127 /
!     data imach( 8) /   60 /
!     data imach( 9) /-1024 /
!     data imach(10) / 1023 /
!
!     Machine constants for the VAX 11/780.
!
!     data imach( 1) /    2 /
!     data imach( 2) /   31 /
!     data imach( 3) / 2147483647 /
!     data imach( 4) /    2 /
!     data imach( 5) /   24 /
!     data imach( 6) / -127 /
!     data imach( 7) /  127 /
!     data imach( 8) /   56 /
!     data imach( 9) / -127 /
!     data imach(10) /  127 /
!
  ipmpar = imach(i)
  return

end
function psi ( xx )
!
!*******************************************************************************
!
!! PSI evaluates the psi or digamma function, d/dx ln(gamma(x)).
!
!
!  Discussion:
!
!    The main computation involves evaluation of rational Chebyshev
!    approximations published in Math. Comp. 27, 123-127(1973) by
!    Cody, Strecok and Thacher.
!
!    PSI was written at Argonne National Laboratory for the FUNPACK
!    package.  PSI was modified by
!    A. H. Morris (NSWC).
!
!  Parameters:
!
!    Input, real XX, the argument of the psi function.
!
!    Output, real PSI, the value of the psi function.  PSI is assigned
!    the value 0 when the psi function is undefined.
!
  double precision, parameter :: dx0 = 1.461632144968362341262659542325721325d0
  double precision, parameter :: piov4 = 0.785398163397448d0
!
  double precision aug
  double precision den
  double precision dpmpar
  integer i
  integer ipmpar
  integer m
  integer n
  integer nq
  double precision, parameter, dimension ( 7 ) :: p1 = (/ &
   0.895385022981970d-02, &
   0.477762828042627d+01, &
   0.142441585084029d+03, &
   0.118645200713425d+04, &
   0.363351846806499d+04, &
   0.413810161269013d+04, &
   0.130560269827897d+04/)

  double precision p2(4)
  double precision psi
  double precision q1(6)
  double precision q2(4)
  double precision sgn
  double precision upper
  double precision w
  double precision x
  double precision xmax1
  double precision xmx0
  double precision xsmall
  double precision xx
  double precision z
!
!  Coefficients for rational approximation of
!  PSI(X) / (X - X0),  0.5 .LE. X .LE. 3.0
!

  data  q1(1)/.448452573429826d+02/
  data  q1(2)/.520752771467162d+03/
  data  q1(3)/.221000799247830d+04/
  data  q1(4)/.364127349079381d+04/
  data  q1(5)/.190831076596300d+04/
  data  q1(6)/.691091682714533d-05/
!
!  Coefficients for rational approximation of
!  PSI(X) - LN(X) + 1 / (2*X),  X .GT. 3.0
!
  data p2(1)/-.212940445131011d+01/
  data p2(2)/-.701677227766759d+01/
  data p2(3)/-.448616543918019d+01/
  data p2(4)/-.648157123766197d+00/

  data q2(1) / .322703493791143d+02/
  data q2(2) / .892920700481861d+02/
  data q2(3) / .546117738103215d+02/
  data q2(4) / .777788548522962d+01/
!
!  XMAX1 is the largest positive floating point constant with entirely 
!  integer representation.  It is also used as negative of lower bound 
!  on acceptable negative arguments and as the positive argument beyond which
!  psi may be represented as LOG(X).
!
  xmax1 = dble ( ipmpar(3) )
  xmax1 = min ( xmax1, 1.0d0 / epsilon ( 1.0d0 ) )
!
!  XSMALL is the absolute argument below which PI*COTAN(PI*X)
!  may be represented by 1/X.
!
  xsmall = 1.0d-9

  x = xx
  aug = 0.0d0

  if ( x == 0.0d0 ) then
    psi = 0.0
    return
  end if
!
!  X < 0.5,  Use reflection formula PSI(1-X) = PSI(X) + PI * COTAN(PI*X)
!
  if ( x < 0.5d0 ) then
!
!  0 < ABS(X) <= XSMALL.  USE 1/X AS A SUBSTITUTE FOR  PI*COTAN(PI*X)
!
    if ( abs ( x ) <= xsmall ) then
      aug = - 1.0d0 / x
      go to 40
    end if
!
!  Reduction of argument for cotan
!
    w = - x
    sgn = piov4

    if ( w <= 0.0d0 ) then
      w = - w
      sgn = - sgn
    end if
!
!  Make an error exit if X .LE. -XMAX1
!
    if ( w >= xmax1 ) then
      psi = 0.0
      return
    end if

    nq = int ( w )
    w = w - dble ( nq )
    nq = int ( w * 4.0d0 )
    w = 4.0d0 * ( w - dble ( nq ) * 0.25d0 )
!
!  W is now related to the fractional part of  4.0 * X.
!  Adjust argument to correspond to values in first
!  quadrant and determine sign
!
    n = nq / 2
    if ( n + n /= nq ) then
      w = 1.0d0 - w
    end if

    z = piov4 * w
    m = n / 2

    if ( m + m /= n ) then
      sgn = - sgn
    end if
!
!  Determine final value for -PI*COTAN(PI*X).
!
    n = ( nq + 1 ) / 2
    m = n / 2
    m = m + m

    if ( m == n ) then

      if ( z == 0.0d0 ) then
        psi = 0.0
       return
    end  if

      aug = 4.0 * sgn * ( cos(z) / sin(z) )

    else

      aug = 4.0 * sgn * ( sin(z) / cos(z) )

    end if

   40   continue

    x = 1.0d0 - x

  end if
!
!  0.5 <= X <= 3.0
!
  if ( x <= 3.0d0 ) then

    den = x
    upper = p1(1) * x

    do i = 1, 5
      den = ( den + q1(i) ) * x
      upper = ( upper + p1(i+1) ) * x
    end do

    den = ( upper + p1(7) ) / ( den + q1(6) )
    xmx0 = dble ( x ) - dx0
    psi = den * xmx0 + aug
!
!  3.0 < X < XMAX1
!
  else if ( x < xmax1 ) then

    w = 1.0d0 / x**2
    den = w
    upper = p2(1) * w

    do i = 1, 3
      den = ( den + q2(i) ) * w
      upper = ( upper + p2(i+1) ) * w
    end do

    aug = upper / ( den + q2(4) ) - 0.5d0 / x + aug
    psi = aug + log ( x )
!
!  XMAX1 <= X
!
  else

    psi = aug + log ( x )

  end if

  return
end
function rcomp ( a, x )
!
!*******************************************************************************
!
!! RCOMP evaluates EXP(-X)*X**A/GAMMA(A)
!
!
!     RT2PIN = 1/SQRT(2*PI)
! 
  double precision, parameter :: rt2pin = 0.398942280401433d0
!
  double precision a
  double precision gam1
  double precision gamma
  double precision rcomp
  double precision rlog
  double precision t
  double precision t1
  double precision u
  double precision x
!
  if ( a < 20.0d0 ) then

    t = a * log ( x ) - x

    if ( a < 1.0d0 ) then
      rcomp = ( a * exp ( t ) ) * ( 1.0d0 + gam1(a) )
    else
      rcomp = exp(t) / gamma(a)
    end if

  else

    u = x / a

    if ( u == 0.0 ) then
      rcomp = 0.0
    else
      t = ( 1.0d0 / a )**2
      t1 = ((( 0.75d0 * t - 1.0d0 ) * t + 3.5d0 ) * t - 105.0d0 ) &
        / ( a * 1260.0d0 )
      t1 = t1 - a * rlog ( u )
      rcomp = rt2pin * sqrt ( a ) * exp ( t1 )
    end if

  end if

  return
end
function rexp ( x )
!
!*******************************************************************************
!
!! REXP evaluates the function EXP(X) - 1.
!
!
!  Modified:
!
!    09 December 1999
!
!  Parameters:
!
!    Input, double precision X, the argument of the function.
!
!    Output, double precision REXP, the value of EXP(X)-1.
!
  double precision, parameter :: p1 =   0.914041914819518d-09
  double precision, parameter :: p2 =   0.238082361044469d-01
  double precision, parameter :: q1 = - 0.499999999085958d+00
  double precision, parameter :: q2 =   0.107141568980644d+00
  double precision, parameter :: q3 = - 0.119041179760821d-01
  double precision, parameter :: q4 =   0.595130811860248d-03
!
  double precision rexp
  double precision w
  double precision x
!
  if ( abs ( x ) <= 0.15d0 ) then

    rexp = x * ((( p2 * x + p1 ) * x + 1.0d0 ) / (((( q4 * x + q3 ) * x + q2 ) &
      * x + q1 ) * x + 1.0d0 ))

  else

    w = exp ( x )

    if ( x <= 0.0d0 ) then
      rexp = ( w - 0.5d0 ) - 0.5d0
    else
      rexp = w * ( 0.5d0 + ( 0.5d0 - 1.0d0 / w ) )
    end if

  end if

  return
end
function rlog ( x )
!
!*******************************************************************************
!
!! RLOG computes  X - 1 - LN(X).
!
!
!  Modified:
!
!    09 December 1999
!
!  Parameters:
!
!    Input, double precision X, the argument of the function.
!
!    Output, double precision RLOG, the value of the function.
!
  double precision, parameter :: a  =   0.566749439387324d-01
  double precision, parameter :: b  =   0.456512608815524d-01
  double precision, parameter :: p0 =   0.333333333333333d+00
  double precision, parameter :: p1 = - 0.224696413112536d+00
  double precision, parameter :: p2 =   0.620886815375787d-02
  double precision, parameter :: q1 = - 0.127408923933623d+01
  double precision, parameter :: q2 =   0.354508718369557d+00
!
  double precision r
  double precision rlog
  double precision t
  double precision u
  double precision w
  double precision w1
  double precision x
!
  if ( x < 0.61d0 ) then

    r = ( x - 0.5d0 ) - 0.5d0
    rlog = r - log ( x )

  else if ( x < 1.57d0 ) then

    if ( x < 0.82d0 ) then

      u = dble ( x ) - 0.7d0
      u = u / 0.7d0
      w1 = a - u * 0.3d0

    else if ( x < 1.18d0 ) then

      u = ( x - 0.5d0 ) - 0.5d0
      w1 = 0.0d0

    else if ( x < 1.57d0 ) then

      u = 0.75d0 * dble ( x ) - 1.d0
      w1 = b + u / 3.0d0

    end if

    r = u/ ( u + 2.0d0 )
    t = r * r
    w = (( p2 * t + p1 ) * t + p0 )/ ( ( q2 * t + q1 ) * t + 1.0d0 )
    rlog = 2.0d0 * t * ( 1.0d0 / ( 1.0d0 - r ) - r * w ) + w1

  else if ( x >= 1.57d0 ) then

    r = ( x - 0.5d0 ) - 0.5d0
    rlog = r - log ( x )

  end if

  return
end
function rlog1 ( x )
!
!*******************************************************************************
!
!! RLOG1 evaluates the function X - LN(1 + X)
!
!
  double precision, parameter :: a = 0.566749439387324d-01
  double precision, parameter :: b = 0.456512608815524d-01
  double precision, parameter :: p0 = 0.333333333333333d+00
  double precision, parameter :: p1 = -0.224696413112536d+00
  double precision, parameter :: p2 = 0.620886815375787d-02
  double precision, parameter :: q1 = -0.127408923933623d+01
  double precision, parameter :: q2 = 0.354508718369557d+00
!
  double precision h
  double precision r
  double precision rlog1
  double precision t
  double precision w
  double precision w1
  double precision x
!
  if (x<-0.39d0 .or. x>0.57d0) go to 40
  if (x<-0.18d0) go to 10
  if (x>0.18d0) go to 20
!
!  Argument reduction.
!
  h = x
  w1 = 0.0d0
  go to 30

10    continue

  h = dble ( x ) + 0.3d0
  h = h / 0.7d0
  w1 = a - h * 0.3d0
  go to 30

   20 h = 0.75d0 * dble ( x ) - 0.25d0
  w1 = b + h / 3.0d0
!
!  Series expansion.
!
   30 r = h / ( h + 2.0d0 )
  t = r * r
  w = (( p2 * t + p1 ) * t + p0 ) / (( q2 * t + q1 ) * t + 1.0d0 )
  rlog1 = 2.0d0 * t * ( 1.0d0 / ( 1.0d0 - r ) - r * w ) + w1
  return

   40 w = ( x + 0.5d0 ) + 0.5d0
  rlog1 = x - log ( w )

  return
end
function stvaln ( p )
!
!*******************************************************************************
!
!! STVALN provides starting values for the inverse of the normal distribution.
!
!
!  Discussion:
!
!    The routine returns X such that 
!      CUMNOR(X) = P,  
!    that is,
!      the integral from
!    -infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P
!
!  Reference:
!
!    Kennedy and Gentle,
!    Statistical Computing, 
!    Marcel Dekker, NY, 1980, page 95.
!
!  Parameters:
!
!     Input, double precision P, the probability whose normal deviate 
!     is sought.
!
!     Output, double precision STVALN, the normal deviate whose probability
!     is P.
!
  double precision eval_pol
  double precision p
  double precision sign
  double precision stvaln

  double precision, parameter, dimension(0:4) :: xden = (/ &
    0.993484626060d-1, &
    0.588581570495d0, &
    0.531103462366d0, &
    0.103537752850d0, &
    0.38560700634d-2 /)

  double precision, parameter, dimension(0:4) :: xnum = (/ &
    -0.322232431088d0, &
    -1.000000000000d0, &
    -0.342242088547d0, &
    -0.204231210245d-1, &
    -0.453642210148d-4 /)

  double precision y
  double precision z
!
  if ( p <= 0.5d0 ) then

    sign = -1.0d0
    z = p

  else

    sign = 1.0d0
    z = 1.0d0 - p

  end if

  y = sqrt ( - 2.0d0 * log ( z ) )
  stvaln = y + eval_pol ( xnum, 4, y ) / eval_pol ( xden, 4, y )
  stvaln = sign * stvaln

  return
end

C==========================================================================================
C Many subroutine/functions about random number generator,
C many kinds of distribution functions and cumulated distribution 
C functions contain to a module
C==========================================================================================
C The lists below are the generation of random numbers and several kinds of distributions
C most variables were declared to single precision
C==========================================================================================
C real function    genchi(df) 			=> generate random value of chi_square variable
C real function    genf(dfn,dfd) 		=> generate random deviate from the F-distr.
C real function    gengam(a,r) 			=> generate random deviate from gamma distr.
C subroutine       genmn(parm,x,work) 		=> generate multivariate normal deviate
C real function    gennor(av,sd) 		=> generate random deviate from a normal distr.
C real function    genunf(low,high) 		=> generate uniform real between low and high
C subroutine       getcgn(g) 			=> get generator
C subroutine       setcgn(g) 			=> set generator
C subroutine       getsd(iseed1,iseed2) 	=> get seeds(return two integer seeds)
C integer function ignlgi() 			=> generate large integer
C integer function ignuin(low,high) 		=> generate uniform integer
C subroutine       initgn(isdtyp) 		=> initialize current generator
C subroutine       inrgcm() 			=> initialize random number generator common
C integer function lennob(string) 		=> length not counting trailing blanks
C integer function mltmod(a,s,m) 		=> returns (a*s) mod m
C subroutine       phrtsd(phrase,seed1,seed2) 	=> phrase to seed2
C logical function qrgnin() 			=> Q random generators initialized?
C logical function qrgnsn(qvalue) 		=> Q random generators set whether initialzed 
C real function    ranf() 			=> random number generator as uniform)
C subroutine       setall(iseed1,iseed2) 	=> set all random number generators
C subroutine       rgnqsd(qssd) 		=> random number generator query seed set?
C subroutine       setgmn(meanv,covm,p,parm) 	=> set generate multivariate normal random deviate
C subroutine       setsd(iseed1,iseed2) 	=> set seed of current generator
C real function    sgamma(a) 			=> standard gamma distr.
C real function    sexpo()			=> standard exponential distr.
C real function    snorm() 			=> standard normal distibution
C subroutine       spofa(a,lda,n,info) 		=> spofa factors a real symmetric PDM
C real function    sdot(n,sx,incx,sy,incy) 	=> dependent to spofa

C==========================================================================================
C The lists below are the cumulative functions for several kinds of distr.
C most variables were declared to double precision
C==========================================================================================
C subroutine       cdfnor(which,p,q,x,mean,sd,status,bound) => CDF for normal distr.
C subroutine       cumnor(arg,result,ccum) 	=> cumulative normal distr.
C double function  devpl(a,n,x) 		=> evaluate a polynomial at x
C double function  dinvnr(p,q) 			=> normal distr. inverse
C integer function ipmpar(i) 			=> provide the integer machine constants
C double function  spmpar(i) 			=> provide the single precision machine constants
C double function  stvaln(p) 			=> start value for Newton-Raphon calc. of Norm. inverse


      module ranlib
!     private
!     public::genchi,genmn,gennor,genunf,getsd,ignuin,phrtsd,ranf,setall
!     public::setgmn,setsd,snorm,cdfnor

      CONTAINS


      REAL FUNCTION genchi(df)
C**********************************************************************
C
C     REAL FUNCTION GENCHI( DF )
C                Generate random value of CHIsquare variable
C
C
C                              Function
C
C
C     Generates random deviate from the distribution of a chisquare
C     with DF degrees of freedom random variable.
C
C
C                              Arguments
C
C
C     DF --> Degrees of freedom of the chisquare
C            (Must be positive)
C                         REAL DF
C
C
C                              Method
C
C
C     Uses relation between chisquare and gamma.
C
C**********************************************************************
C     .. Scalar Arguments ..
      REAL df
C     ..
C     .. External Functions ..
C     REAL gengam
C     EXTERNAL gengam
C     ..
C     .. Executable Statements ..
      IF (.NOT. (df.LE.0.0)) GO TO 10
      WRITE (*,*) 'DF <= 0 in GENCHI - ABORT'
      WRITE (*,*) 'Value of DF: ',df
      STOP 'DF <= 0 in GENCHI - ABORT'

   10 genchi = 2.0*gengam(1.0,df/2.0)
      RETURN

      END FUNCTION


      REAL FUNCTION gengam(a,r)
C**********************************************************************
C
C     REAL FUNCTION GENGAM( A, R )
C           GENerates random deviates from GAMma distribution
C
C
C                              Function
C
C
C     Generates random deviates from the gamma distribution whose
C     density is
C          (A**R)/Gamma(R) * X**(R-1) * Exp(-A*X)
C
C
C                              Arguments
C
C
C     A --> Location parameter of Gamma distribution
C                              REAL A
C
C     R --> Shape parameter of Gamma distribution
C                              REAL R
C
C
C                              Method
C
C
C     Renames SGAMMA from TOMS as slightly modified by BWB to use RANF
C     instead of SUNIF.
C
C     For details see:
C               (Case R >= 1.0)
C               Ahrens, J.H. and Dieter, U.
C               Generating Gamma Variates by a
C               Modified Rejection Technique.
C               Comm. ACM, 25,1 (Jan. 1982), 47 - 54.
C     Algorithm GD
C
C               (Case 0.0 <= R <= 1.0)
C               Ahrens, J.H. and Dieter, U.
C               Computer Methods for Sampling from Gamma,
C               Beta, Poisson and Binomial Distributions.
C               Computing, 12 (1974), 223-246/
C     Adapted algorithm GS.
C
C**********************************************************************
C     .. Scalar Arguments ..
      REAL a,r
C     ..
C     .. External Functions ..
C     REAL sgamma
C     EXTERNAL sgamma
C     ..
C     .. Executable Statements ..
      gengam = sgamma(r)
      gengam = gengam/a
      RETURN

      END FUNCTION
  
      SUBROUTINE genmtn(x,var,np,y)
      REAL x(np),var(np,np),y(np),parm(np*(np+3)/2+1),work(np)
      parm=0.0
      call setgmn(x,var,np,parm)
      call genmn(parm,y,work)
      return
      end subroutine

      SUBROUTINE genmn(parm,x,work)
C**********************************************************************
C
C     SUBROUTINE GENMN(PARM,X,WORK)
C              GENerate Multivariate Normal random deviate
C
C
C                              Arguments
C
C
C     PARM --> Parameters needed to generate multivariate normal
C               deviates (MEANV and Cholesky decomposition of
C               COVM). Set by a previous call to SETGMN.
C               1 : 1                - size of deviate, P
C               2 : P + 1            - mean vector
C               P+2 : P*(P+3)/2 + 1  - upper half of cholesky
C                                       decomposition of cov matrix
C                                             REAL PARM(*)
C
C     X    <-- Vector deviate generated.
C                                             REAL X(P)
C
C     WORK <--> Scratch array
C                                             REAL WORK(P)
C
C
C                              Method
C
C
C     1) Generate P independent standard normal deviates - Ei ~ N(0,1)
C
C     2) Using Cholesky decomposition find A s.t. trans(A)*A = COVM
C
C     3) trans(A)E + MEANV ~ N(MEANV,COVM)
C
C**********************************************************************
C     .. Array Arguments ..
      REAL parm(*),work(*),x(*)
C     ..
C     .. Local Scalars ..
      REAL ae
      INTEGER i,icount,j,p
C     ..
C     .. External Functions ..
C     REAL snorm
C     EXTERNAL snorm
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC int
C     ..
C     .. Executable Statements ..
      p = int(parm(1))
C
C     Generate P independent normal deviates - WORK ~ N(0,1)
C
      DO 10,i = 1,p
          work(i) = snorm()
   10 CONTINUE
      DO 30,i = 1,p
C
C     PARM (P+2 : P*(P+3)/2 + 1) contains A, the Cholesky
C      decomposition of the desired covariance matrix.
C          trans(A)(1,1) = PARM(P+2)
C          trans(A)(2,1) = PARM(P+3)
C          trans(A)(2,2) = PARM(P+2+P)
C          trans(A)(3,1) = PARM(P+4)
C          trans(A)(3,2) = PARM(P+3+P)
C          trans(A)(3,3) = PARM(P+2-1+2P)  ...
C
C     trans(A)*WORK + MEANV ~ N(MEANV,COVM)
C
          icount = 0
          ae = 0.0
          DO 20,j = 1,i
              icount = icount + j - 1
              ae = ae + parm(i+ (j-1)*p-icount+p+1)*work(j)
   20     CONTINUE
          x(i) = ae + parm(i+1)
   30 CONTINUE
      RETURN
C
      END SUBROUTINE


      REAL FUNCTION gennor(av,sd)
C**********************************************************************
C
C     REAL FUNCTION GENNOR( AV, SD )
C
C         GENerate random deviate from a NORmal distribution
C
C
C                              Function
C
C
C     Generates a single random deviate from a normal distribution
C     with mean, AV, and standard deviation, SD.
C
C
C                              Arguments
C
C
C     AV --> Mean of the normal distribution.
C                              REAL AV
C
C     SD --> Standard deviation of the normal distribution.
C                              REAL SD
C
C     GENNOR <-- Generated normal deviate.
C                              REAL GENNOR
C
C
C                              Method
C
C
C     Renames SNORM from TOMS as slightly modified by BWB to use RANF
C     instead of SUNIF.
C
C     For details see:
C               Ahrens, J.H. and Dieter, U.
C               Extensions of Forsythe's Method for Random
C               Sampling from the Normal Distribution.
C               Math. Comput., 27,124 (Oct. 1973), 927 - 937.
C
C
C**********************************************************************
C     .. Scalar Arguments ..
      REAL av,sd
C     ..
C     .. External Functions ..
C     REAL snorm
C     EXTERNAL snorm
C     ..
C     .. Executable Statements ..
      gennor = sd*snorm() + av
      RETURN

      END FUNCTION


      REAL FUNCTION genunf(low,high)
C**********************************************************************
C
C     REAL FUNCTION GENUNF( LOW, HIGH )
C
C               GeNerate Uniform Real between LOW and HIGH
C
C
C                              Function
C
C
C     Generates a real uniformly distributed between LOW and HIGH.
C
C
C                              Arguments
C
C
C     LOW --> Low bound (exclusive) on real value to be generated
C                         REAL LOW
C
C     HIGH --> High bound (exclusive) on real value to be generated
C                         REAL HIGH
C
C**********************************************************************
C     .. Scalar Arguments ..
      REAL high,low
C     ..
C     .. External Functions ..
C     REAL ranf
C     EXTERNAL ranf
C     ..
C     .. Executable Statements ..
      IF (.NOT. (low.GT.high)) GO TO 10
      WRITE (*,*) 'LOW > HIGH in GENUNF: LOW ',low,' HIGH: ',high
      WRITE (*,*) 'Abort'
      STOP 'LOW > High in GENUNF - Abort'

   10 genunf = low + (high-low)*ranf()

      RETURN

      END FUNCTION

      SUBROUTINE getcgn(g)
      INTEGER g
C**********************************************************************
C
C      SUBROUTINE GETCGN(G)
C                         Get GeNerator
C
C     Returns in G the number of the current random number generator
C
C
C                              Arguments
C
C
C     G <-- Number of the current random number generator (1..32)
C                    INTEGER G
C
C**********************************************************************
C
      INTEGER curntg,numg
      SAVE curntg
      PARAMETER (numg=32)
      DATA curntg/1/
C
      g = curntg
      RETURN

      ENTRY setcgn(g)
C**********************************************************************
C
C     SUBROUTINE SETCGN( G )
C                      Set GeNerator
C
C     Sets  the  current  generator to G.    All references to a generat
C     are to the current generator.
C
C
C                              Arguments
C
C
C     G --> Number of the current random number generator (1..32)
C                    INTEGER G
C
C**********************************************************************
C
C     Abort if generator number out of range
C
      IF (.NOT. (g.LT.0.OR.g.GT.numg)) GO TO 10
      WRITE (*,*) ' Generator number out of range in SETCGN:',
     +  ' Legal range is 1 to ',numg,' -- ABORT!'
      STOP ' Generator number out of range in SETCGN'

   10 curntg = g
      RETURN

      END SUBROUTINE

      SUBROUTINE getsd(iseed1,iseed2)
C**********************************************************************
C
C     SUBROUTINE GETSD(G,ISEED1,ISEED2)
C               GET SeeD
C
C     Returns the value of two integer seeds of the current generator
C
C     This  is   a  transcription from  Pascal   to  Fortran  of routine
C     Get_State from the paper
C
C     L'Ecuyer, P. and  Cote,  S. "Implementing a Random Number  Package
C     with   Splitting Facilities."  ACM  Transactions   on Mathematical
C     Software, 17:98-111 (1991)
C
C
C                              Arguments
C
C
C
C     ISEED1 <- First integer seed of generator G
C                                   INTEGER ISEED1
C
C     ISEED2 <- Second integer seed of generator G
C                                   INTEGER ISEED1
C
C**********************************************************************
C     .. Parameters ..
      INTEGER numg
      PARAMETER (numg=32)
C     ..
C     .. Scalar Arguments ..
      INTEGER iseed1,iseed2
C     ..
C     .. Scalars in Common ..
      INTEGER a1,a1vw,a1w,a2,a2vw,a2w,m1,m2
C     ..
C     .. Arrays in Common ..
      INTEGER cg1(numg),cg2(numg),ig1(numg),ig2(numg),lg1(numg),
     +        lg2(numg)
      LOGICAL qanti(numg)
C     ..
C     .. Local Scalars ..
      INTEGER g
C     ..
C     .. External Functions ..
C     LOGICAL qrgnin
C     EXTERNAL qrgnin
C     ..
C     .. External Subroutines ..
C     EXTERNAL getcgn
C     ..
C     .. Common blocks ..
      COMMON /globe/m1,m2,a1,a2,a1w,a2w,a1vw,a2vw,ig1,ig2,lg1,lg2,cg1,
     +       cg2,qanti
C     ..
C     .. Save statement ..
      SAVE /globe/
C     ..
C     .. Executable Statements ..
C     Abort unless random number generator initialized
      IF (qrgnin()) GO TO 10
      WRITE (*,*) ' GETSD called before random number generator ',
     +  ' initialized -- abort!'
      STOP ' GETSD called before random number generator initialized'

   10 CALL getcgn(g)
      iseed1 = cg1(g)
      iseed2 = cg2(g)
      RETURN

      END SUBROUTINE


      INTEGER FUNCTION ignlgi()
C**********************************************************************
C
C     INTEGER FUNCTION IGNLGI()
C               GeNerate LarGe Integer
C
C     Returns a random integer following a uniform distribution over
C     (1, 2147483562) using the current generator.
C
C     This is a transcription from Pascal to Fortran of routine
C     Random from the paper
C
C     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
C     with Splitting Facilities." ACM Transactions on Mathematical
C     Software, 17:98-111 (1991)
C
C**********************************************************************
C     .. Parameters ..
      INTEGER numg
      PARAMETER (numg=32)
C     ..
C     .. Scalars in Common ..
      INTEGER a1,a1vw,a1w,a2,a2vw,a2w,m1,m2
C     ..
C     .. Arrays in Common ..
      INTEGER cg1(numg),cg2(numg),ig1(numg),ig2(numg),lg1(numg),
     +        lg2(numg)
      LOGICAL qanti(numg)
C     ..
C     .. Local Scalars ..
      INTEGER curntg,k,s1,s2,z
      LOGICAL qqssd
C     ..
C     .. External Functions ..
C     LOGICAL qrgnin
C     EXTERNAL qrgnin
C     ..
C     .. External Subroutines ..
C     EXTERNAL getcgn,inrgcm,rgnqsd,setall
C     ..
C     .. Common blocks ..
      COMMON /globe/m1,m2,a1,a2,a1w,a2w,a1vw,a2vw,ig1,ig2,lg1,lg2,cg1,
     +       cg2,qanti
C     ..
C     .. Save statement ..
      SAVE /globe/
C     ..
C     .. Executable Statements ..
C
C     IF THE RANDOM NUMBER PACKAGE HAS NOT BEEN INITIALIZED YET, DO SO.
C     IT CAN BE INITIALIZED IN ONE OF TWO WAYS : 1) THE FIRST CALL TO
C     THIS ROUTINE  2) A CALL TO SETALL.
C
      IF (.NOT. (qrgnin())) CALL inrgcm()
      CALL rgnqsd(qqssd)
      IF (.NOT. (qqssd)) CALL setall(1234567890,123456789)
C
C     Get Current Generator
C
      CALL getcgn(curntg)
      s1 = cg1(curntg)
      s2 = cg2(curntg)
      k = s1/53668
      s1 = a1* (s1-k*53668) - k*12211
      IF (s1.LT.0) s1 = s1 + m1
      k = s2/52774
      s2 = a2* (s2-k*52774) - k*3791
      IF (s2.LT.0) s2 = s2 + m2
      cg1(curntg) = s1
      cg2(curntg) = s2
      z = s1 - s2
      IF (z.LT.1) z = z + m1 - 1
      IF (qanti(curntg)) z = m1 - z
      ignlgi = z
      RETURN

      END FUNCTION


      INTEGER FUNCTION ignuin(low,high)
C**********************************************************************
C
C     INTEGER FUNCTION IGNUIN( LOW, HIGH )
C
C               GeNerate Uniform INteger
C
C
C                              Function
C
C
C     Generates an integer uniformly distributed between LOW and HIGH.
C
C
C                              Arguments
C
C
C     LOW --> Low bound (inclusive) on integer value to be generated
C                         INTEGER LOW
C
C     HIGH --> High bound (inclusive) on integer value to be generated
C                         INTEGER HIGH
C
C
C                              Note
C
C
C     If (HIGH-LOW) > 2,147,483,561 prints error message on * unit and
C     stops the program.
C
C**********************************************************************

C     IGNLGI generates integers between 1 and 2147483562
C     MAXNUM is 1 less than maximum generable value
C     .. Parameters ..
      INTEGER maxnum
      PARAMETER (maxnum=2147483561)
      CHARACTER*(*) err1,err2
      PARAMETER (err1='LOW > HIGH in IGNUIN',
     +          err2=' ( HIGH - LOW ) > 2,147,483,561 in IGNUIN')
C     ..
C     .. Scalar Arguments ..
      INTEGER high,low
C     ..
C     .. Local Scalars ..
      INTEGER err,ign,maxnow,range,ranp1
C     ..
C     .. External Functions ..
C     INTEGER ignlgi
C     EXTERNAL ignlgi
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC mod
C     ..
C     .. Executable Statements ..
      IF (.NOT. (low.GT.high)) GO TO 10
      err = 1
C      ABORT-PROGRAM
      GO TO 80

   10 range = high - low
      IF (.NOT. (range.GT.maxnum)) GO TO 20
      err = 2
C      ABORT-PROGRAM
      GO TO 80

   20 IF (.NOT. (low.EQ.high)) GO TO 30
      ignuin = low
      RETURN

C     Number to be generated should be in range 0..RANGE
C     Set MAXNOW so that the number of integers in 0..MAXNOW is an
C     integral multiple of the number in 0..RANGE

   30 ranp1 = range + 1
      maxnow = (maxnum/ranp1)*ranp1
   40 ign = ignlgi() - 1
      IF (.NOT. (ign.LE.maxnow)) GO TO 50
      ignuin = low + mod(ign,ranp1)
      RETURN

   50 GO TO 40

   60 CONTINUE
   70 CONTINUE
   80 IF (.NOT. (err.EQ.1)) GO TO 90
      WRITE (*,*) err1
      GO TO 100

C     TO ABORT-PROGRAM
   90 WRITE (*,*) err2
  100 WRITE (*,*) ' LOW: ',low,' HIGH: ',high
      WRITE (*,*) ' Abort on Fatal ERROR'
      IF (.NOT. (err.EQ.1)) then
        print*,' ( HIGH - LOW ) > 2,147,483,561 in IGNUIN'
      else
        print*, 'LOW > HIGH in IGNUIN'
      endif

  120 END FUNCTION

      SUBROUTINE initgn(isdtyp)
C**********************************************************************
C
C     SUBROUTINE INITGN(ISDTYP)
C          INIT-ialize current G-e-N-erator
C
C     Reinitializes the state of the current generator
C
C     This is a transcription from Pascal to Fortran of routine
C     Init_Generator from the paper
C
C     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
C     with Splitting Facilities." ACM Transactions on Mathematical
C     Software, 17:98-111 (1991)
C
C
C                              Arguments
C
C
C     ISDTYP -> The state to which the generator is to be set
C
C          ISDTYP = -1  => sets the seeds to their initial value
C          ISDTYP =  0  => sets the seeds to the first value of
C                          the current block
C          ISDTYP =  1  => sets the seeds to the first value of
C                          the next block
C
C                                   INTEGER ISDTYP
C
C**********************************************************************
C     .. Parameters ..
      INTEGER numg
      PARAMETER (numg=32)
C     ..
C     .. Scalar Arguments ..
      INTEGER isdtyp
C     ..
C     .. Scalars in Common ..
      INTEGER a1,a1vw,a1w,a2,a2vw,a2w,m1,m2
C     ..
C     .. Arrays in Common ..
      INTEGER cg1(numg),cg2(numg),ig1(numg),ig2(numg),lg1(numg),
     +        lg2(numg)
      LOGICAL qanti(numg)
C     ..
C     .. Local Scalars ..
      INTEGER g
C     ..
C     .. External Functions ..
C     LOGICAL qrgnin
C     INTEGER mltmod
C     EXTERNAL qrgnin,mltmod
C     ..
C     .. External Subroutines ..
C     EXTERNAL getcgn
C     ..
C     .. Common blocks ..
      COMMON /globe/m1,m2,a1,a2,a1w,a2w,a1vw,a2vw,ig1,ig2,lg1,lg2,cg1,
     +       cg2,qanti
C     ..
C     .. Save statement ..
      SAVE /globe/
C     ..
C     .. Executable Statements ..
C     Abort unless random number generator initialized
      IF (qrgnin()) GO TO 10
      WRITE (*,*) ' INITGN called before random number generator ',
     +  ' initialized -- abort!'
      STOP ' INITGN called before random number generator initialized'

   10 CALL getcgn(g)
      IF ((-1).NE. (isdtyp)) GO TO 20
      lg1(g) = ig1(g)
      lg2(g) = ig2(g)
      GO TO 50

   20 IF ((0).NE. (isdtyp)) GO TO 30
      CONTINUE
      GO TO 50
C     do nothing
   30 IF ((1).NE. (isdtyp)) GO TO 40
      lg1(g) = mltmod(a1w,lg1(g),m1)
      lg2(g) = mltmod(a2w,lg2(g),m2)
      GO TO 50

   40 STOP 'ISDTYP NOT IN RANGE'

   50 cg1(g) = lg1(g)
      cg2(g) = lg2(g)
      RETURN

      END SUBROUTINE

      SUBROUTINE inrgcm()
C**********************************************************************
C
C     SUBROUTINE INRGCM()
C          INitialize Random number Generator CoMmon
C
C
C                              Function
C
C
C     Initializes common area  for random number  generator.  This saves
C     the  nuisance  of  a  BLOCK DATA  routine  and the  difficulty  of
C     assuring that the routine is loaded with the other routines.
C
C**********************************************************************
C     .. Parameters ..
      INTEGER numg
      PARAMETER (numg=32)
C     ..
C     .. Scalars in Common ..
      INTEGER a1,a1vw,a1w,a2,a2vw,a2w,m1,m2
C     ..
C     .. Arrays in Common ..
      INTEGER cg1(numg),cg2(numg),ig1(numg),ig2(numg),lg1(numg),
     +        lg2(numg)
      LOGICAL qanti(numg)
C     ..
C     .. Local Scalars ..
      INTEGER i
      LOGICAL qdum
C     ..
C     .. External Functions ..
C     LOGICAL qrgnsn
C     EXTERNAL qrgnsn
C     ..
C     .. Common blocks ..
      COMMON /globe/m1,m2,a1,a2,a1w,a2w,a1vw,a2vw,ig1,ig2,lg1,lg2,cg1,
     +       cg2,qanti
C     ..
C     .. Save statement ..
      SAVE /globe/
C     ..
C     .. Executable Statements ..
C     V=20;                            W=30;
C
C     A1W = MOD(A1**(2**W),M1)         A2W = MOD(A2**(2**W),M2)
C     A1VW = MOD(A1**(2**(V+W)),M1)    A2VW = MOD(A2**(2**(V+W)),M2)
C
C   If V or W is changed A1W, A2W, A1VW, and A2VW need to be recomputed.
C    An efficient way to precompute a**(2*j) MOD m is to start with
C    a and square it j times modulo m using the function MLTMOD.
C
      m1 = 2147483563
      m2 = 2147483399
      a1 = 40014
      a2 = 40692
      a1w = 1033780774
      a2w = 1494757890
      a1vw = 2082007225
      a2vw = 784306273
      DO 10,i = 1,numg
          qanti(i) = .FALSE.
   10 CONTINUE
C
C     Tell the world that common has been initialized
C
      qdum = qrgnsn(.TRUE.)
      RETURN

      END SUBROUTINE

      INTEGER FUNCTION lennob(string)
      IMPLICIT INTEGER (a-p,r-z),LOGICAL (q)
C**********************************************************************
C
C     INTEGER FUNCTION LENNOB( STRING )
C                LENgth NOt counting trailing Blanks
C
C
C                              Function
C
C
C     Returns the length of STRING up to and including the last
C     non-blank character.
C
C
C                              Arguments
C
C
C     STRING --> String whose length not counting trailing blanks
C                is returned.
C
C**********************************************************************
      CHARACTER*(*) string

      end = len(string)
      DO 20,i = end,1,-1
          IF (.NOT. (string(i:i).NE.' ')) GO TO 10
          lennob = i
          RETURN

   10     CONTINUE
   20 CONTINUE
      lennob = 0
      RETURN

      END FUNCTION

      INTEGER FUNCTION mltmod(a,s,m)
C**********************************************************************
C
C     INTEGER FUNCTION MLTMOD(A,S,M)
C
C                    Returns (A*S) MOD M
C
C     This is a transcription from Pascal to Fortran of routine
C     MULtMod_Decompos from the paper
C
C     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
C     with Splitting Facilities." ACM Transactions on Mathematical
C     Software, 17:98-111 (1991)
C
C
C                              Arguments
C
C
C     A, S, M  -->
C                         INTEGER A,S,M
C
C**********************************************************************
C     .. Parameters ..
      INTEGER h
      PARAMETER (h=32768)
C     ..
C     .. Scalar Arguments ..
      INTEGER a,m,s
C     ..
C     .. Local Scalars ..
      INTEGER a0,a1,k,p,q,qh,rh
C     ..
C     .. Executable Statements ..
C
C     H = 2**((b-2)/2) where b = 32 because we are using a 32 bit
C      machine. On a different machine recompute H
C
      IF (.NOT. (a.LE.0.OR.a.GE.m.OR.s.LE.0.OR.s.GE.m)) GO TO 10
      WRITE (*,*) ' A, M, S out of order in MLTMOD - ABORT!'
      WRITE (*,*) ' A = ',a,' S = ',s,' M = ',m
      WRITE (*,*) ' MLTMOD requires: 0 < A < M; 0 < S < M'
      STOP ' A, M, S out of order in MLTMOD - ABORT!'

   10 IF (.NOT. (a.LT.h)) GO TO 20
      a0 = a
      p = 0
      GO TO 120

   20 a1 = a/h
      a0 = a - h*a1
      qh = m/h
      rh = m - h*qh
      IF (.NOT. (a1.GE.h)) GO TO 50
      a1 = a1 - h
      k = s/qh
      p = h* (s-k*qh) - k*rh
   30 IF (.NOT. (p.LT.0)) GO TO 40
      p = p + m
      GO TO 30

   40 GO TO 60

   50 p = 0
C
C     P = (A2*S*H)MOD M
C
   60 IF (.NOT. (a1.NE.0)) GO TO 90
      q = m/a1
      k = s/q
      p = p - k* (m-a1*q)
      IF (p.GT.0) p = p - m
      p = p + a1* (s-k*q)
   70 IF (.NOT. (p.LT.0)) GO TO 80
      p = p + m
      GO TO 70

   80 CONTINUE
   90 k = p/qh
C
C     P = ((A2*H + A1)*S)MOD M
C
      p = h* (p-k*qh) - k*rh
  100 IF (.NOT. (p.LT.0)) GO TO 110
      p = p + m
      GO TO 100

  110 CONTINUE
  120 IF (.NOT. (a0.NE.0)) GO TO 150
C
C     P = ((A2*H + A1)*H*S)MOD M
C
      q = m/a0
      k = s/q
      p = p - k* (m-a0*q)
      IF (p.GT.0) p = p - m
      p = p + a0* (s-k*q)
  130 IF (.NOT. (p.LT.0)) GO TO 140
      p = p + m
      GO TO 130

  140 CONTINUE
  150 mltmod = p
C
      RETURN

      END FUNCTION

      SUBROUTINE phrtsd(phrase,seed1,seed2)
C**********************************************************************
C
C     SUBROUTINE PHRTSD( PHRASE, SEED1, SEED2 )
C               PHRase To SeeDs
C
C
C                              Function
C
C
C     Uses a phrase (character string) to generate two seeds for the RGN
C     random number generator.
C
C
C                              Arguments
C
C
C     PHRASE --> Phrase to be used for random number generation
C                         CHARACTER*(*) PHRASE
C
C     SEED1 <-- First seed for RGN generator
C                         INTEGER SEED1
C
C     SEED2 <-- Second seed for RGN generator
C                         INTEGER SEED2
C
C
C                              Note
C
C
C     Trailing blanks are eliminated before the seeds are generated.
C
C     Generated seed values will fall in the range 1..2^30
C     (1..1,073,741,824)
C
C**********************************************************************
C     .. Parameters ..
      CHARACTER*(*) table
      PARAMETER (table='abcdefghijklmnopqrstuvwxyz'//
     +          'ABCDEFGHIJKLMNOPQRSTUVWXYZ'//'0123456789'//
     +          '!@#$%^&*()_+[];:''"<>?,./')
      INTEGER twop30
      PARAMETER (twop30=1073741824)
C     ..
C     .. Scalar Arguments ..
      INTEGER seed1,seed2
      CHARACTER phrase* (*)
C     ..
C     .. Local Scalars ..
      INTEGER i,ichr,j,lphr
C     ..
C     .. Local Arrays ..
      INTEGER shift(0:4),values(5)
C     ..
C     .. External Functions ..
C     INTEGER lennob
C     EXTERNAL lennob
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC index,mod
C     ..
C     .. Data statements ..
      DATA shift/1,64,4096,262144,16777216/
C     ..
C     .. Executable Statements ..
      seed1 = 1234567890
      seed2 = 123456789
      lphr = lennob(phrase)
      IF (lphr.LT.1) RETURN
      DO 30,i = 1,lphr
          ichr = mod(index(table,phrase(i:i)),64)
          IF (ichr.EQ.0) ichr = 63
          DO 10,j = 1,5
              values(j) = ichr - j
              IF (values(j).LT.1) values(j) = values(j) + 63
   10     CONTINUE
          DO 20,j = 1,5
              seed1 = mod(seed1+shift(j-1)*values(j),twop30)
              seed2 = mod(seed2+shift(j-1)*values(6-j),twop30)
   20     CONTINUE
   30 CONTINUE
      RETURN

      END SUBROUTINE

      LOGICAL FUNCTION qrgnin()
C**********************************************************************
C
C     LOGICAL FUNCTION QRGNIN()
C               Q Random GeNerators INitialized?
C
C     A trivial routine to determine whether or not the random
C     number generator has been initialized.  Returns .TRUE. if
C     it has, else .FALSE.
C
C**********************************************************************
C     .. Scalar Arguments ..
      LOGICAL qvalue
C     ..
C     .. Local Scalars ..
      LOGICAL qinit
C     ..
C     .. Entry Points ..
      LOGICAL qrgnsn
C     ..
C     .. Save statement ..
      SAVE qinit
C     ..
C     .. Data statements ..
      DATA qinit/.FALSE./
C     ..
C     .. Executable Statements ..
      qrgnin = qinit
      RETURN

      ENTRY qrgnsn(qvalue)
C**********************************************************************
C
C     LOGICAL FUNCTION QRGNSN( QVALUE )
C               Q Random GeNerators Set whether iNitialized
C
C     Sets state of whether random number generator is initialized
C     to QVALUE.
C
C     This routine is actually an entry in QRGNIN, hence it is a
C     logical function.  It returns the (meaningless) value .TRUE.
C
C**********************************************************************
      qinit = qvalue
      qrgnsn = .TRUE.
      RETURN

      END FUNCTION

      REAL FUNCTION ranf()
C**********************************************************************
C
C     REAL FUNCTION RANF()
C                RANDom number generator as a Function
C
C     Returns a random floating point number from a uniform distribution
C     over 0 - 1 (endpoints of this interval are not returned) using the
C     current generator
C
C     This is a transcription from Pascal to Fortran of routine
C     Uniform_01 from the paper
C
C     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
C     with Splitting Facilities." ACM Transactions on Mathematical
C     Software, 17:98-111 (1991)
C
C**********************************************************************
C     .. External Functions ..
C     INTEGER ignlgi
C     EXTERNAL ignlgi
C     ..
C     .. Executable Statements ..
C
C     4.656613057E-10 is 1/M1  M1 is set in a data statement in IGNLGI
C      and is currently 2147483563. If M1 changes, change this also.
C
      ranf = ignlgi()*4.656613057E-10
      RETURN

      END FUNCTION

      SUBROUTINE setall(iseed1,iseed2)
C**********************************************************************
C
C      SUBROUTINE SETALL(ISEED1,ISEED2)
C               SET ALL random number generators
C
C     Sets the initial seed of generator 1 to ISEED1 and ISEED2. The
C     initial seeds of the other generators are set accordingly, and
C     all generators states are set to these seeds.
C
C     This is a transcription from Pascal to Fortran of routine
C     Set_Initial_Seed from the paper
C
C     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
C     with Splitting Facilities." ACM Transactions on Mathematical
C     Software, 17:98-111 (1991)
C
C
C                              Arguments
C
C
C     ISEED1 -> First of two integer seeds
C                                   INTEGER ISEED1
C
C     ISEED2 -> Second of two integer seeds
C                                   INTEGER ISEED1
C
C**********************************************************************
C     .. Parameters ..
      INTEGER numg
      PARAMETER (numg=32)
C     ..
C     .. Scalar Arguments ..
      INTEGER iseed1,iseed2
      LOGICAL qssd
C     ..
C     .. Scalars in Common ..
      INTEGER a1,a1vw,a1w,a2,a2vw,a2w,m1,m2
C     ..
C     .. Arrays in Common ..
      INTEGER cg1(numg),cg2(numg),ig1(numg),ig2(numg),lg1(numg),
     +        lg2(numg)
      LOGICAL qanti(numg)
C     ..
C     .. Local Scalars ..
      INTEGER g,ocgn
      LOGICAL qqssd
C     ..
C     .. External Functions ..
C     INTEGER mltmod
C     LOGICAL qrgnin
C     EXTERNAL mltmod,qrgnin
C     ..
C     .. External Subroutines ..
C     EXTERNAL getcgn,initgn,inrgcm,setcgn
C     ..
C     .. Common blocks ..
      COMMON /globe/m1,m2,a1,a2,a1w,a2w,a1vw,a2vw,ig1,ig2,lg1,lg2,cg1,
     +       cg2,qanti
C     ..
C     .. Save statement ..
      SAVE /globe/,qqssd
C     ..
C     .. Data statements ..
      DATA qqssd/.FALSE./
C     ..
C     .. Executable Statements ..
C
C     TELL IGNLGI, THE ACTUAL NUMBER GENERATOR, THAT THIS ROUTINE
C      HAS BEEN CALLED.
C
      qqssd = .TRUE.
      CALL getcgn(ocgn)
C
C     Initialize Common Block if Necessary
C
      IF (.NOT. (qrgnin())) CALL inrgcm()
      ig1(1) = iseed1
      ig2(1) = iseed2
      CALL initgn(-1)
      DO 10,g = 2,numg
          ig1(g) = mltmod(a1vw,ig1(g-1),m1)
          ig2(g) = mltmod(a2vw,ig2(g-1),m2)
          CALL setcgn(g)
          CALL initgn(-1)
   10 CONTINUE
      CALL setcgn(ocgn)
      RETURN

      ENTRY rgnqsd(qssd)
C**********************************************************************
C
C     SUBROUTINE RGNQSD
C                    Random Number Generator Query SeeD set?
C
C     Returns (LOGICAL) QSSD as .TRUE. if SETALL has been invoked,
C     otherwise returns .FALSE.
C
C**********************************************************************
      qssd = qqssd
      RETURN

      END SUBROUTINE


      SUBROUTINE setgmn(meanv,covm,p,parm)
C**********************************************************************
C
C     SUBROUTINE SETGMN( MEANV, COVM, P, PARM)
C            SET Generate Multivariate Normal random deviate
C
C
C                              Function
C
C
C      Places P, MEANV, and the Cholesky factoriztion of COVM
C      in GENMN.
C
C
C                              Arguments
C
C
C     MEANV --> Mean vector of multivariate normal distribution.
C                                        REAL MEANV(P)
C
C     COVM   <--> (Input) Covariance   matrix    of  the  multivariate
C                 normal distribution
C                 (Output) Destroyed on output
C                                        REAL COVM(P,P)
C
C     P     --> Dimension of the normal, or length of MEANV.
C                                        INTEGER P
C
C     PARM <-- Array of parameters needed to generate multivariate norma
C                deviates (P, MEANV and Cholesky decomposition of
C                COVM).
C                1 : 1                - P
C                2 : P + 1            - MEANV
C                P+2 : P*(P+3)/2 + 1  - Cholesky decomposition of COVM
C                                             REAL PARM(P*(P+3)/2 + 1)
C
C**********************************************************************
C     .. Scalar Arguments ..
      INTEGER p
C     ..
C     .. Array Arguments ..
      REAL covm(p,p),meanv(p),parm(p* (p+3)/2+1)
C     ..
C     .. Local Scalars ..
      INTEGER i,icount,info,j
C     ..
C     .. External Subroutines ..
C     EXTERNAL spofa
C     ..
C     .. Executable Statements ..
C
C
C     TEST THE INPUT
C
      IF (.NOT. (p.LE.0)) GO TO 10
      WRITE (*,*) 'P nonpositive in SETGMN'
      WRITE (*,*) 'Value of P: ',p
      STOP 'P nonpositive in SETGMN'

   10 parm(1) = p
C
C     PUT P AND MEANV INTO PARM
C
      DO 20,i = 2,p + 1
          parm(i) = meanv(i-1)
   20 CONTINUE
C
C      Cholesky decomposition to find A s.t. trans(A)*(A) = COVM
C
      CALL spofa(covm,p,p,info)
      IF (.NOT. (info.NE.0)) GO TO 30
      WRITE (*,*) ' COVM not positive definite in SETGMN'
      STOP ' COVM not positive definite in SETGMN'

   30 icount = p + 1
C
C     PUT UPPER HALF OF A, WHICH IS NOW THE CHOLESKY FACTOR, INTO PARM
C          COVM(1,1) = PARM(P+2)
C          COVM(1,2) = PARM(P+3)
C                    :
C          COVM(1,P) = PARM(2P+1)
C          COVM(2,2) = PARM(2P+2)  ...
C
      DO 50,i = 1,p
          DO 40,j = i,p
              icount = icount + 1
              parm(icount) = covm(i,j)
   40     CONTINUE
   50 CONTINUE
      RETURN
C
      END SUBROUTINE

      SUBROUTINE setsd(iseed1,iseed2)
C**********************************************************************
C
C     SUBROUTINE SETSD(ISEED1,ISEED2)
C               SET S-ee-D of current generator
C
C     Resets the initial  seed of  the current  generator to  ISEED1 and
C     ISEED2. The seeds of the other generators remain unchanged.
C
C     This is a transcription from Pascal to Fortran of routine
C     Set_Seed from the paper
C
C     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
C     with Splitting Facilities." ACM Transactions on Mathematical
C     Software, 17:98-111 (1991)
C
C
C                              Arguments
C
C
C     ISEED1 -> First integer seed
C                                   INTEGER ISEED1
C
C     ISEED2 -> Second integer seed
C                                   INTEGER ISEED1
C
C**********************************************************************
C     .. Parameters ..
      INTEGER numg
      PARAMETER (numg=32)
C     ..
C     .. Scalar Arguments ..
      INTEGER iseed1,iseed2
C     ..
C     .. Scalars in Common ..
      INTEGER a1,a1vw,a1w,a2,a2vw,a2w,m1,m2
C     ..
C     .. Arrays in Common ..
      INTEGER cg1(numg),cg2(numg),ig1(numg),ig2(numg),lg1(numg),
     +        lg2(numg)
      LOGICAL qanti(numg)
C     ..
C     .. Local Scalars ..
      INTEGER g
C     ..
C     .. External Functions ..
C     LOGICAL qrgnin
C     EXTERNAL qrgnin
C     ..
C     .. External Subroutines ..
C     EXTERNAL getcgn,initgn
C     ..
C     .. Common blocks ..
      COMMON /globe/m1,m2,a1,a2,a1w,a2w,a1vw,a2vw,ig1,ig2,lg1,lg2,cg1,
     +       cg2,qanti
C     ..
C     .. Save statement ..
      SAVE /globe/
C     ..
C     .. Executable Statements ..
C     Abort unless random number generator initialized
      IF (qrgnin()) GO TO 10
      WRITE (*,*) ' SETSD called before random number generator ',
     +  ' initialized -- abort!'
      STOP ' SETSD called before random number generator initialized'

   10 CALL getcgn(g)
      ig1(g) = iseed1
      ig2(g) = iseed2
      CALL initgn(-1)
      RETURN

      END SUBROUTINE


      REAL FUNCTION sgamma(a)
C**********************************************************************C
C                                                                      C
C                                                                      C
C     (STANDARD-)  G A M M A  DISTRIBUTION                             C
C                                                                      C
C                                                                      C
C**********************************************************************C
C**********************************************************************C
C                                                                      C
C               PARAMETER  A >= 1.0  !                                 C
C                                                                      C
C**********************************************************************C
C                                                                      C
C     FOR DETAILS SEE:                                                 C
C                                                                      C
C               AHRENS, J.H. AND DIETER, U.                            C
C               GENERATING GAMMA VARIATES BY A                         C
C               MODIFIED REJECTION TECHNIQUE.                          C
C               COMM. ACM, 25,1 (JAN. 1982), 47 - 54.                  C
C                                                                      C
C     STEP NUMBERS CORRESPOND TO ALGORITHM 'GD' IN THE ABOVE PAPER     C
C                                 (STRAIGHTFORWARD IMPLEMENTATION)     C
C                                                                      C
C     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   C
C     SUNIF.  The argument IR thus goes away.                          C
C                                                                      C
C**********************************************************************C
C                                                                      C
C               PARAMETER  0.0 < A < 1.0  !                            C
C                                                                      C
C**********************************************************************C
C                                                                      C
C     FOR DETAILS SEE:                                                 C
C                                                                      C
C               AHRENS, J.H. AND DIETER, U.                            C
C               COMPUTER METHODS FOR SAMPLING FROM GAMMA,              C
C               BETA, POISSON AND BINOMIAL DISTRIBUTIONS.              C
C               COMPUTING, 12 (1974), 223 - 246.                       C
C                                                                      C
C     (ADAPTED IMPLEMENTATION OF ALGORITHM 'GS' IN THE ABOVE PAPER)    C
C                                                                      C
C**********************************************************************C
C
C
C     INPUT: A =PARAMETER (MEAN) OF THE STANDARD GAMMA DISTRIBUTION
C     OUTPUT: SGAMMA = SAMPLE FROM THE GAMMA-(A)-DISTRIBUTION
C
C     COEFFICIENTS Q(K) - FOR Q0 = SUM(Q(K)*A**(-K))
C     COEFFICIENTS A(K) - FOR Q = Q0+(T*T/2)*SUM(A(K)*V**K)
C     COEFFICIENTS E(K) - FOR EXP(Q)-1 = SUM(E(K)*Q**K)
C
      DATA q1,q2,q3,q4,q5,q6,q7/.04166669,.02083148,.00801191,.00144121,
     +     -.00007388,.00024511,.00024240/
      DATA a1,a2,a3,a4,a5,a6,a7/.3333333,-.2500030,.2000062,-.1662921,
     +     .1423657,-.1367177,.1233795/
      DATA e1,e2,e3,e4,e5/1.,.4999897,.1668290,.0407753,.0102930/
C
C     PREVIOUS A PRE-SET TO ZERO - AA IS A', AAA IS A"
C     SQRT32 IS THE SQUAREROOT OF 32 = 5.656854249492380
C
      DATA aa/0.0/,aaa/0.0/,sqrt32/5.656854/
C
C     SAVE STATEMENTS
C
      SAVE aa,aaa,s2,s,d,q0,b,si,c
C
      IF (a.EQ.aa) GO TO 10
      IF (a.LT.1.0) GO TO 140
C
C     STEP  1:  RECALCULATIONS OF S2,S,D IF A HAS CHANGED
C
      aa = a
      s2 = a - 0.5
      s = sqrt(s2)
      d = sqrt32 - 12.0*s
C
C     STEP  2:  T=STANDARD NORMAL DEVIATE,
C               X=(S,1/2)-NORMAL DEVIATE.
C               IMMEDIATE ACCEPTANCE (I)
C
   10 t = snorm()
      x = s + 0.5*t
      sgamma = x*x
      IF (t.GE.0.0) RETURN
C
C     STEP  3:  U= 0,1 -UNIFORM SAMPLE. SQUEEZE ACCEPTANCE (S)
C
      u = ranf()
      IF (d*u.LE.t*t*t) RETURN
C
C     STEP  4:  RECALCULATIONS OF Q0,B,SI,C IF NECESSARY
C
      IF (a.EQ.aaa) GO TO 40
      aaa = a
      r = 1.0/a
      q0 = ((((((q7*r+q6)*r+q5)*r+q4)*r+q3)*r+q2)*r+q1)*r
C
C               APPROXIMATION DEPENDING ON SIZE OF PARAMETER A
C               THE CONSTANTS IN THE EXPRESSIONS FOR B, SI AND
C               C WERE ESTABLISHED BY NUMERICAL EXPERIMENTS
C
      IF (a.LE.3.686) GO TO 30
      IF (a.LE.13.022) GO TO 20
C
C               CASE 3:  A .GT. 13.022
C
      b = 1.77
      si = .75
      c = .1515/s
      GO TO 40
C
C               CASE 2:  3.686 .LT. A .LE. 13.022
C
   20 b = 1.654 + .0076*s2
      si = 1.68/s + .275
      c = .062/s + .024
      GO TO 40
C
C               CASE 1:  A .LE. 3.686
C
   30 b = .463 + s + .178*s2
      si = 1.235
      c = .195/s - .079 + .16*s
C
C     STEP  5:  NO QUOTIENT TEST IF X NOT POSITIVE
C
   40 IF (x.LE.0.0) GO TO 70
C
C     STEP  6:  CALCULATION OF V AND QUOTIENT Q
C
      v = t/ (s+s)
      IF (abs(v).LE.0.25) GO TO 50
      q = q0 - s*t + 0.25*t*t + (s2+s2)*alog(1.0+v)
      GO TO 60

   50 q = q0 + 0.5*t*t* ((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v
C
C     STEP  7:  QUOTIENT ACCEPTANCE (Q)
C
   60 IF (alog(1.0-u).LE.q) RETURN
C
C     STEP  8:  E=STANDARD EXPONENTIAL DEVIATE
C               U= 0,1 -UNIFORM DEVIATE
C               T=(B,SI)-DOUBLE EXPONENTIAL (LAPLACE) SAMPLE
C
   70 e = sexpo()
      u = ranf()
      u = u + u - 1.0
      t = b + sign(si*e,u)
      IF (.NOT. (u.GE.0.0)) GO TO 80
      t = b + si*e
      GO TO 90

   80 t = b - si*e

C
C     STEP  9:  REJECTION IF T .LT. TAU(1) = -.71874483771719
C
   90 IF (t.LT. (-.7187449)) GO TO 70
C
C     STEP 10:  CALCULATION OF V AND QUOTIENT Q
C
      v = t/ (s+s)
      IF (abs(v).LE.0.25) GO TO 100
      q = q0 - s*t + 0.25*t*t + (s2+s2)*alog(1.0+v)
      GO TO 110

  100 q = q0 + 0.5*t*t* ((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v
C
C     STEP 11:  HAT ACCEPTANCE (H) (IF Q NOT POSITIVE GO TO STEP 8)
C
  110 IF (q.LE.0.0) GO TO 70
      IF (q.LE.0.5) GO TO 120
      w = exp(q) - 1.0
      GO TO 130

  120 w = ((((e5*q+e4)*q+e3)*q+e2)*q+e1)*q
C
C               IF T IS REJECTED, SAMPLE AGAIN AT STEP 8
C
  130 IF (c*abs(u).GT.w*exp(e-0.5*t*t)) GO TO 70
      x = s + 0.5*t
      sgamma = x*x
      RETURN
C
C     ALTERNATE METHOD FOR PARAMETERS A BELOW 1  (.3678794=EXP(-1.))
C
  140 aa = 0.0
      b = 1.0 + .3678794*a
  150 p = b*ranf()
      IF (p.GE.1.0) GO TO 160
      sgamma = exp(alog(p)/a)
      IF (sexpo().LT.sgamma) GO TO 150
      RETURN

  160 sgamma = -alog((b-p)/a)
      IF (sexpo().LT. (1.0-a)*alog(sgamma)) GO TO 150
      RETURN

      END FUNCTION


      REAL FUNCTION sexpo()
C**********************************************************************C
C                                                                      C
C                                                                      C
C     (STANDARD-)  E X P O N E N T I A L   DISTRIBUTION                C
C                                                                      C
C                                                                      C
C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     FOR DETAILS SEE:                                                 C
C                                                                      C
C               AHRENS, J.H. AND DIETER, U.                            C
C               COMPUTER METHODS FOR SAMPLING FROM THE                 C
C               EXPONENTIAL AND NORMAL DISTRIBUTIONS.                  C
C               COMM. ACM, 15,10 (OCT. 1972), 873 - 882.               C
C                                                                      C
C     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM       C
C     'SA' IN THE ABOVE PAPER (SLIGHTLY MODIFIED IMPLEMENTATION)       C
C                                                                      C
C     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   C
C     SUNIF.  The argument IR thus goes away.                          C
C                                                                      C
C**********************************************************************C
C
      DIMENSION q(8)
      EQUIVALENCE (q(1),q1)
C
C     Q(N) = SUM(ALOG(2.0)**K/K!)    K=1,..,N ,      THE HIGHEST N
C     (HERE 8) IS DETERMINED BY Q(N)=1.0 WITHIN STANDARD PRECISION
C
      DATA q/.6931472,.9333737,.9888778,.9984959,.9998293,.9999833,
     +     .9999986,.9999999/
C
   10 a = 0.0
      u = ranf()
      GO TO 30

   20 a = a + q1
   30 u = u + u
      IF (u.LE.1.0) GO TO 20
   40 u = u - 1.0
      IF (u.GT.q1) GO TO 60
   50 sexpo = a + u
      RETURN

   60 i = 1
      ustar = ranf()
      umin = ustar
   70 ustar = ranf()
      IF (ustar.LT.umin) umin = ustar
   80 i = i + 1
      IF (u.GT.q(i)) GO TO 70
   90 sexpo = a + umin*q1
      RETURN

      END FUNCTION


      REAL FUNCTION snorm()
C**********************************************************************C
C                                                                      C
C                                                                      C
C     (STANDARD-)  N O R M A L  DISTRIBUTION                           C
C                                                                      C
C                                                                      C
C**********************************************************************C
C**********************************************************************C
C                                                                      C
C     FOR DETAILS SEE:                                                 C
C                                                                      C
C               AHRENS, J.H. AND DIETER, U.                            C
C               EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM             C
C               SAMPLING FROM THE NORMAL DISTRIBUTION.                 C
C               MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.          C
C                                                                      C
C     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'  C
C     (M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)  C
C                                                                      C
C     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of   C
C     SUNIF.  The argument IR thus goes away.                          C
C                                                                      C
C**********************************************************************C
C
      DIMENSION a(32),d(31),t(31),h(31)
C
C     THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
C     H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
C
      DATA a/0.0,.3917609E-1,.7841241E-1,.1177699,.1573107,.1970991,
     +     .2372021,.2776904,.3186394,.3601299,.4022501,.4450965,
     +     .4887764,.5334097,.5791322,.6260990,.6744898,.7245144,
     +     .7764218,.8305109,.8871466,.9467818,1.009990,1.077516,
     +     1.150349,1.229859,1.318011,1.417797,1.534121,1.675940,
     +     1.862732,2.153875/
      DATA d/5*0.0,.2636843,.2425085,.2255674,.2116342,.1999243,
     +     .1899108,.1812252,.1736014,.1668419,.1607967,.1553497,
     +     .1504094,.1459026,.1417700,.1379632,.1344418,.1311722,
     +     .1281260,.1252791,.1226109,.1201036,.1177417,.1155119,
     +     .1134023,.1114027,.1095039/
      DATA t/.7673828E-3,.2306870E-2,.3860618E-2,.5438454E-2,
     +     .7050699E-2,.8708396E-2,.1042357E-1,.1220953E-1,.1408125E-1,
     +     .1605579E-1,.1815290E-1,.2039573E-1,.2281177E-1,.2543407E-1,
     +     .2830296E-1,.3146822E-1,.3499233E-1,.3895483E-1,.4345878E-1,
     +     .4864035E-1,.5468334E-1,.6184222E-1,.7047983E-1,.8113195E-1,
     +     .9462444E-1,.1123001,.1364980,.1716886,.2276241,.3304980,
     +     .5847031/
      DATA h/.3920617E-1,.3932705E-1,.3950999E-1,.3975703E-1,
     +     .4007093E-1,.4045533E-1,.4091481E-1,.4145507E-1,.4208311E-1,
     +     .4280748E-1,.4363863E-1,.4458932E-1,.4567523E-1,.4691571E-1,
     +     .4833487E-1,.4996298E-1,.5183859E-1,.5401138E-1,.5654656E-1,
     +     .5953130E-1,.6308489E-1,.6737503E-1,.7264544E-1,.7926471E-1,
     +     .8781922E-1,.9930398E-1,.1155599,.1404344,.1836142,.2790016,
     +     .7010474/
C
   10 u = ranf()
      s = 0.0
      IF (u.GT.0.5) s = 1.0
      u = u + u - s
   20 u = 32.0*u
      i = int(u)
      IF (i.EQ.32) i = 31
      IF (i.EQ.0) GO TO 100
C
C                                START CENTER
C
   30 ustar = u - float(i)
      aa = a(i)
   40 IF (ustar.LE.t(i)) GO TO 60
      w = (ustar-t(i))*h(i)
C
C                                EXIT   (BOTH CASES)
C
   50 y = aa + w
      snorm = y
      IF (s.EQ.1.0) snorm = -y
      RETURN
C
C                                CENTER CONTINUED
C
   60 u = ranf()
      w = u* (a(i+1)-aa)
      tt = (0.5*w+aa)*w
      GO TO 80

   70 tt = u
      ustar = ranf()
   80 IF (ustar.GT.tt) GO TO 50
   90 u = ranf()
      IF (ustar.GE.u) GO TO 70
      ustar = ranf()
      GO TO 40
C
C                                START TAIL
C
  100 i = 6
      aa = a(32)
      GO TO 120

  110 aa = aa + d(i)
      i = i + 1
  120 u = u + u
      IF (u.LT.1.0) GO TO 110
  130 u = u - 1.0
  140 w = u*d(i)
      tt = (0.5*w+aa)*w
      GO TO 160

  150 tt = u
  160 ustar = ranf()
      IF (ustar.GT.tt) GO TO 50
  170 u = ranf()
      IF (ustar.GE.u) GO TO 150
      u = ranf()
      GO TO 140

      END FUNCTION


      SUBROUTINE spofa(a,lda,n,info)
      INTEGER lda,n,info
      REAL a(lda,1)
C
C     SPOFA FACTORS A REAL SYMMETRIC POSITIVE DEFINITE MATRIX.
C
C     SPOFA IS USUALLY CALLED BY SPOCO, BUT IT CAN BE CALLED
C     DIRECTLY WITH A SAVING IN TIME IF  RCOND  IS NOT NEEDED.
C     (TIME FOR SPOCO) = (1 + 18/N)*(TIME FOR SPOFA) .
C
C     ON ENTRY
C
C        A       REAL(LDA, N)
C                THE SYMMETRIC MATRIX TO BE FACTORED.  ONLY THE
C                DIAGONAL AND UPPER TRIANGLE ARE USED.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C     ON RETURN
C
C        A       AN UPPER TRIANGULAR MATRIX  R  SO THAT  A = TRANS(R)*R
C                WHERE  TRANS(R)  IS THE TRANSPOSE.
C                THE STRICT LOWER TRIANGLE IS UNALTERED.
C                IF  INFO .NE. 0 , THE FACTORIZATION IS NOT COMPLETE.
C
C        INFO    INTEGER
C                = 0  FOR NORMAL RETURN.
C                = K  SIGNALS AN ERROR CONDITION.  THE LEADING MINOR
C                     OF ORDER  K  IS NOT POSITIVE DEFINITE.
C
C     LINPACK.  THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS SDOT
C     FORTRAN SQRT
C
C     INTERNAL VARIABLES
C
C     REAL sdot
      REAL s,t
      INTEGER j,jm1,k
C     BEGIN BLOCK WITH ...EXITS TO 40
C
C
      DO 30 j = 1,n
          info = j
          s = 0.0E0
          jm1 = j - 1
          IF (jm1.LT.1) GO TO 20
          DO 10 k = 1,jm1
              t = a(k,j) - sdot(k-1,a(1,k),1,a(1,j),1)
              t = t/a(k,k)
              a(k,j) = t
              s = s + t*t
   10     CONTINUE
   20     CONTINUE
          s = a(j,j) - s
C     ......EXIT
          IF (s.LE.0.0E0) GO TO 40
          a(j,j) = sqrt(s)
   30 CONTINUE
      info = 0
   40 CONTINUE
      RETURN

      END SUBROUTINE

      REAL FUNCTION sdot(n,sx,incx,sy,incy)
      REAL sx(1),sy(1),stemp
      INTEGER i,incx,incy,ix,iy,m,mp1,n

      stemp = 0.0E0
      sdot = 0.0E0
      IF (n.LE.0) RETURN
      IF (incx.EQ.1 .AND. incy.EQ.1) GO TO 20
      ix = 1
      iy = 1
      IF (incx.LT.0) ix = (-n+1)*incx + 1
      IF (incy.LT.0) iy = (-n+1)*incy + 1
      DO 10 i = 1,n
          stemp = stemp + sx(ix)*sy(iy)
          ix = ix + incx
          iy = iy + incy
   10 CONTINUE
      sdot = stemp
      RETURN

   20 m = mod(n,5)
      IF (m.EQ.0) GO TO 40
      DO 30 i = 1,m
          stemp = stemp + sx(i)*sy(i)
   30 CONTINUE
      IF (n.LT.5) GO TO 60
   40 mp1 = m + 1
      DO 50 i = mp1,n,5
          stemp = stemp + sx(i)*sy(i) + sx(i+1)*sy(i+1) +
     +            sx(i+2)*sy(i+2) + sx(i+3)*sy(i+3) + sx(i+4)*sy(i+4)
   50 CONTINUE
   60 sdot = stemp
      RETURN

      END FUNCTION

C These below were about CDF


      SUBROUTINE cdfnor(which,p,q,x,mean,sd,status,bound)
C**********************************************************************
C
C      SUBROUTINE CDFNOR( WHICH, P, Q, X, MEAN, SD, STATUS, BOUND )
C               Cumulative Distribution Function
C               NORmal distribution
C
C
C                              Function
C
C
C     Calculates any one parameter of the normal
C     distribution given values for the others.
C
C
C                              Arguments
C
C
C     WHICH  --> Integer indicating  which of the  next  parameter
C     values is to be calculated using values  of the others.
C     Legal range: 1..4
C               iwhich = 1 : Calculate P and Q from X,MEAN and SD
C               iwhich = 2 : Calculate X from P,Q,MEAN and SD
C               iwhich = 3 : Calculate MEAN from P,Q,X and SD
C               iwhich = 4 : Calculate SD from P,Q,X and MEAN
C                    INTEGER WHICH
C
C     P <--> The integral from -infinity to X of the normal density.
C            Input range: (0,1].
C                    DOUBLE PRECISION P
C
C     Q <--> 1-P.
C            Input range: (0, 1].
C            P + Q = 1.0.
C                    DOUBLE PRECISION Q
C
C     X < --> Upper limit of integration of the normal-density.
C             Input range: ( -infinity, +infinity)
C                    DOUBLE PRECISION X
C
C     MEAN <--> The mean of the normal density.
C               Input range: (-infinity, +infinity)
C                    DOUBLE PRECISION MEAN
C
C     SD <--> Standard Deviation of the normal density.
C             Input range: (0, +infinity).
C                    DOUBLE PRECISION SD
C
C     STATUS <-- 0 if calculation completed correctly
C               -I if input parameter number I is out of range
C                1 if answer appears to be lower than lowest
C                  search bound
C                2 if answer appears to be higher than greatest
C                  search bound
C                3 if P + Q .ne. 1
C                    INTEGER STATUS
C
C     BOUND <-- Undefined if STATUS is 0
C
C               Bound exceeded by parameter number I if STATUS
C               is negative.
C
C               Lower search bound if STATUS is 1.
C
C               Upper search bound if STATUS is 2.
C
C
C                              Method
C
C
C
C
C     A slightly modified version of ANORM from
C
C     Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN
C     Package of Special Function Routines and Test Drivers"
C     acm Transactions on Mathematical Software. 19, 22-32.
C
C     is used to calulate the  cumulative standard normal distribution.
C
C     The rational functions from pages  90-95  of Kennedy and Gentle,
C     Statistical  Computing,  Marcel  Dekker, NY,  1980 are  used  as
C     starting values to Newton's Iterations which compute the inverse
C     standard normal.  Therefore no  searches  are necessary for  any
C     parameter.
C
C     For X < -15, the asymptotic expansion for the normal is used  as
C     the starting value in finding the inverse standard normal.
C     This is formula 26.2.12 of Abramowitz and Stegun.
C
C
C                              Note
C
C
C      The normal density is proportional to
C      exp( - 0.5 * (( X - MEAN)/SD)**2)
C
C
C**********************************************************************
C     .. Parameters ..
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION bound,mean,p,q,sd,x
      INTEGER status,which
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION z,pq
C     ..
C     .. External Functions ..
C     DOUBLE PRECISION dinvnr,spmpar
C     EXTERNAL dinvnr,spmpar
C     ..
C     .. External Subroutines ..
C     EXTERNAL cumnor
C     ..
C     .. Executable Statements ..
C
C     Check arguments
C
      status = 0
      IF (.NOT. ((which.LT.1).OR. (which.GT.4))) GO TO 30
      IF (.NOT. (which.LT.1)) GO TO 10
      bound = 1.0D0
      GO TO 20

   10 bound = 4.0D0
   20 status = -1
      RETURN

   30 IF (which.EQ.1) GO TO 70
C
C     P
C
      IF (.NOT. ((p.LE.0.0D0).OR. (p.GT.1.0D0))) GO TO 60
      IF (.NOT. (p.LE.0.0D0)) GO TO 40
      bound = 0.0D0
      GO TO 50

   40 bound = 1.0D0
   50 status = -2
      RETURN

   60 CONTINUE
   70 IF (which.EQ.1) GO TO 110
C
C     Q
C
      IF (.NOT. ((q.LE.0.0D0).OR. (q.GT.1.0D0))) GO TO 100
      IF (.NOT. (q.LE.0.0D0)) GO TO 80
      bound = 0.0D0
      GO TO 90

   80 bound = 1.0D0
   90 status = -3
      RETURN

  100 CONTINUE
  110 IF (which.EQ.1) GO TO 150
C
C     P + Q
C
      pq = p + q
      IF (.NOT. (abs(((pq)-0.5D0)-0.5D0).GT.
     +    (3.0D0*spmpar(1)))) GO TO 140
      IF (.NOT. (pq.LT.0.0D0)) GO TO 120
      bound = 0.0D0
      GO TO 130

  120 bound = 1.0D0
  130 status = 3
      RETURN

  140 CONTINUE
  150 IF (which.EQ.4) GO TO 170
C
C     SD
C
      IF (.NOT. (sd.LE.0.0D0)) GO TO 160
      bound = 0.0D0
      status = -6
      RETURN

  160 CONTINUE
C
C     Calculate ANSWERS
C
  170 IF ((1).EQ. (which)) THEN
C
C     Computing P
C
          z = (x-mean)/sd
          CALL cumnor(z,p,q)

      ELSE IF ((2).EQ. (which)) THEN
C
C     Computing X
C
          z = dinvnr(p,q)
          x = sd*z + mean

      ELSE IF ((3).EQ. (which)) THEN
C
C     Computing the MEAN
C
          z = dinvnr(p,q)
          mean = x - sd*z

      ELSE IF ((4).EQ. (which)) THEN
C
C     Computing SD
C
          z = dinvnr(p,q)
          sd = (x-mean)/z
      END IF

      RETURN

      END SUBROUTINE


      SUBROUTINE cumnor(arg,result,ccum)
C**********************************************************************
C
C     SUBROUINE CUMNOR(X,RESULT,CCUM)
C
C
C                              Function
C
C
C     Computes the cumulative  of    the  normal   distribution,   i.e.,
C     the integral from -infinity to x of
C          (1/sqrt(2*pi)) exp(-u*u/2) du
C
C     X --> Upper limit of integration.
C                                        X is DOUBLE PRECISION
C
C     RESULT <-- Cumulative normal distribution.
C                                        RESULT is DOUBLE PRECISION
C
C     CCUM <-- Compliment of Cumulative normal distribution.
C                                        CCUM is DOUBLE PRECISION
C
C
C     Renaming of function ANORM from:
C
C     Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN
C     Package of Special Function Routines and Test Drivers"
C     acm Transactions on Mathematical Software. 19, 22-32.
C
C     with slight modifications to return ccum and to deal with
C     machine constants.
C
C**********************************************************************
C
C
C Original Comments:
C------------------------------------------------------------------
C
C This function evaluates the normal distribution function:
C
C                              / x
C                     1       |       -t*t/2
C          P(x) = ----------- |      e       dt
C                 sqrt(2 pi)  |
C                             /-oo
C
C   The main computation evaluates near-minimax approximations
C   derived from those in "Rational Chebyshev approximations for
C   the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
C   This transportable program uses rational functions that
C   theoretically approximate the normal distribution function to
C   at least 18 significant decimal digits.  The accuracy achieved
C   depends on the arithmetic system, the compiler, the intrinsic
C   functions, and proper selection of the machine-dependent
C   constants.
C
C*******************************************************************
C*******************************************************************
C
C Explanation of machine-dependent constants.
C
C   MIN   = smallest machine representable number.
C
C   EPS   = argument below which anorm(x) may be represented by
C           0.5  and above which  x*x  will not underflow.
C           A conservative value is the largest machine number X
C           such that   1.0 + X = 1.0   to machine precision.
C*******************************************************************
C*******************************************************************
C
C Error returns
C
C  The program returns  ANORM = 0     for  ARG .LE. XLOW.
C
C
C Intrinsic functions required are:
C
C     ABS, AINT, EXP
C
C
C  Author: W. J. Cody
C          Mathematics and Computer Science Division
C          Argonne National Laboratory
C          Argonne, IL 60439
C
C  Latest modification: March 15, 1992
C
C------------------------------------------------------------------
      INTEGER i
      DOUBLE PRECISION a,arg,b,c,d,del,eps,half,p,one,q,result,sixten,
     +                 temp,sqrpi,thrsh,root32,x,xden,xnum,y,xsq,zero,
     +                 min,ccum
      DIMENSION a(5),b(4),c(9),d(8),p(6),q(5)
C------------------------------------------------------------------
C  External Function
C------------------------------------------------------------------
C     DOUBLE PRECISION spmpar
C     EXTERNAL spmpar
C------------------------------------------------------------------
C  Mathematical constants
C
C  SQRPI = 1 / sqrt(2*pi), ROOT32 = sqrt(32), and
C  THRSH is the argument for which anorm = 0.75.
C------------------------------------------------------------------
      DATA one,half,zero,sixten/1.0D0,0.5D0,0.0D0,1.60D0/,
     +     sqrpi/3.9894228040143267794D-1/,thrsh/0.66291D0/,
     +     root32/5.656854248D0/
C------------------------------------------------------------------
C  Coefficients for approximation in first interval
C------------------------------------------------------------------
      DATA a/2.2352520354606839287D00,1.6102823106855587881D02,
     +     1.0676894854603709582D03,1.8154981253343561249D04,
     +     6.5682337918207449113D-2/
      DATA b/4.7202581904688241870D01,9.7609855173777669322D02,
     +     1.0260932208618978205D04,4.5507789335026729956D04/
C------------------------------------------------------------------
C  Coefficients for approximation in second interval
C------------------------------------------------------------------
      DATA c/3.9894151208813466764D-1,8.8831497943883759412D00,
     +     9.3506656132177855979D01,5.9727027639480026226D02,
     +     2.4945375852903726711D03,6.8481904505362823326D03,
     +     1.1602651437647350124D04,9.8427148383839780218D03,
     +     1.0765576773720192317D-8/
      DATA d/2.2266688044328115691D01,2.3538790178262499861D02,
     +     1.5193775994075548050D03,6.4855582982667607550D03,
     +     1.8615571640885098091D04,3.4900952721145977266D04,
     +     3.8912003286093271411D04,1.9685429676859990727D04/
C------------------------------------------------------------------
C  Coefficients for approximation in third interval
C------------------------------------------------------------------
      DATA p/2.1589853405795699D-1,1.274011611602473639D-1,
     +     2.2235277870649807D-2,1.421619193227893466D-3,
     +     2.9112874951168792D-5,2.307344176494017303D-2/
      DATA q/1.28426009614491121D00,4.68238212480865118D-1,
     +     6.59881378689285515D-2,3.78239633202758244D-3,
     +     7.29751555083966205D-5/
C------------------------------------------------------------------
C  Machine dependent constants
C------------------------------------------------------------------
      eps = spmpar(1)*0.5D0
      min = spmpar(2)
C------------------------------------------------------------------
      x = arg
      y = abs(x)
      IF (y.LE.thrsh) THEN
C------------------------------------------------------------------
C  Evaluate  anorm  for  |X| <= 0.66291
C------------------------------------------------------------------
          xsq = zero
          IF (y.GT.eps) xsq = x*x
          xnum = a(5)*xsq
          xden = xsq
          DO 10 i = 1,3
              xnum = (xnum+a(i))*xsq
              xden = (xden+b(i))*xsq
   10     CONTINUE
          result = x* (xnum+a(4))/ (xden+b(4))
          temp = result
          result = half + temp
          ccum = half - temp
C------------------------------------------------------------------
C  Evaluate  anorm  for 0.66291 <= |X| <= sqrt(32)
C------------------------------------------------------------------
      ELSE IF (y.LE.root32) THEN
          xnum = c(9)*y
          xden = y
          DO 20 i = 1,7
              xnum = (xnum+c(i))*y
              xden = (xden+d(i))*y
   20     CONTINUE
          result = (xnum+c(8))/ (xden+d(8))
          xsq = aint(y*sixten)/sixten
          del = (y-xsq)* (y+xsq)
          result = exp(-xsq*xsq*half)*exp(-del*half)*result
          ccum = one - result
          IF (x.GT.zero) THEN
              temp = result
              result = ccum
              ccum = temp
          END IF
C------------------------------------------------------------------
C  Evaluate  anorm  for |X| > sqrt(32)
C------------------------------------------------------------------
      ELSE
          result = zero
          xsq = one/ (x*x)
          xnum = p(6)*xsq
          xden = xsq
          DO 30 i = 1,4
              xnum = (xnum+p(i))*xsq
              xden = (xden+q(i))*xsq
   30     CONTINUE
          result = xsq* (xnum+p(5))/ (xden+q(5))
          result = (sqrpi-result)/y
          xsq = aint(x*sixten)/sixten
          del = (x-xsq)* (x+xsq)
          result = exp(-xsq*xsq*half)*exp(-del*half)*result
          ccum = one - result
          IF (x.GT.zero) THEN
              temp = result
              result = ccum
              ccum = temp
          END IF

      END IF

      IF (result.LT.min) result = 0.0D0
      IF (ccum.LT.min) ccum = 0.0D0
C------------------------------------------------------------------
C  Fix up for negative argument, erf, etc.
C------------------------------------------------------------------
C----------Last card of ANORM ----------

      END SUBROUTINE


      DOUBLE PRECISION FUNCTION devlpl(a,n,x)
C**********************************************************************
C
C     DOUBLE PRECISION FUNCTION DEVLPL(A,N,X)
C              Double precision EVALuate a PoLynomial at X
C
C
C                              Function
C
C
C     returns
C          A(1) + A(2)*X + ... + A(N)*X**(N-1)
C
C
C                              Arguments
C
C
C     A --> Array of coefficients of the polynomial.
C                                        A is DOUBLE PRECISION(N)
C
C     N --> Length of A, also degree of polynomial - 1.
C                                        N is INTEGER
C
C     X --> Point at which the polynomial is to be evaluated.
C                                        X is DOUBLE PRECISION
C
C**********************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION x
      INTEGER n
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION a(n)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION term
      INTEGER i
C     ..
C     .. Executable Statements ..
      term = a(n)
      DO 10,i = n - 1,1,-1
          term = a(i) + term*x
   10 CONTINUE
      devlpl = term
      RETURN

      END FUNCTION


      DOUBLE PRECISION FUNCTION dinvnr(p,q)
C**********************************************************************
C
C     DOUBLE PRECISION FUNCTION DINVNR(P,Q)
C     Double precision NoRmal distribution INVerse
C
C
C                              Function
C
C
C     Returns X  such that CUMNOR(X)  =   P,  i.e., the  integral from -
C     infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P
C
C
C                              Arguments
C
C
C     P --> The probability whose normal deviate is sought.
C                    P is DOUBLE PRECISION
C
C     Q --> 1-P
C                    P is DOUBLE PRECISION
C
C
C                              Method
C
C
C     The  rational   function   on  page 95    of Kennedy  and  Gentle,
C     Statistical Computing, Marcel Dekker, NY , 1980 is used as a start
C     value for the Newton method of finding roots.
C
C
C                              Note
C
C
C     If P or Q .lt. machine EPS returns +/- DINVNR(EPS)
C
C**********************************************************************
C     .. Parameters ..
      INTEGER maxit
      PARAMETER (maxit=100)
      DOUBLE PRECISION eps
      PARAMETER (eps=1.0D-13)
      DOUBLE PRECISION r2pi
      PARAMETER (r2pi=0.3989422804014326D0)
      DOUBLE PRECISION nhalf
      PARAMETER (nhalf=-0.5D0)
C     ..
C     .. Scalar Arguments ..
      DOUBLE PRECISION p,q
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION strtx,xcur,cum,ccum,pp,dx
      INTEGER i
      LOGICAL qporq
C     ..
C     .. External Functions ..
C     DOUBLE PRECISION stvaln
C     EXTERNAL stvaln
C     ..
C     .. External Subroutines ..
C     EXTERNAL cumnor
C     ..
C     .. Statement Functions ..
      DOUBLE PRECISION dennor,x

      dennor(x) = r2pi*exp(nhalf*x*x)
C     ..
C     .. Executable Statements ..
C
C     FIND MINIMUM OF P AND Q
C
      qporq = p .LE. q
      IF (.NOT. (qporq)) GO TO 10
      pp = p
      GO TO 20

   10 pp = q
C
C     INITIALIZATION STEP
C
   20 strtx = stvaln(pp)
      xcur = strtx
C
C     NEWTON INTERATIONS
C
      DO 30,i = 1,maxit
          CALL cumnor(xcur,cum,ccum)
          dx = (cum-pp)/dennor(xcur)
          xcur = xcur - dx
          IF (abs(dx/xcur).LT.eps) GO TO 40
   30 CONTINUE
      dinvnr = strtx
C
C     IF WE GET HERE, NEWTON HAS FAILED
C
      IF (.NOT.qporq) dinvnr = -dinvnr
      RETURN
C
C     IF WE GET HERE, NEWTON HAS SUCCEDED
C
   40 dinvnr = xcur
      IF (.NOT.qporq) dinvnr = -dinvnr
      RETURN

      END FUNCTION


      INTEGER FUNCTION ipmpar(i)
C-----------------------------------------------------------------------
C
C     IPMPAR PROVIDES THE INTEGER MACHINE CONSTANTS FOR THE COMPUTER
C     THAT IS USED. IT IS ASSUMED THAT THE ARGUMENT I IS AN INTEGER
C     HAVING ONE OF THE VALUES 1-10. IPMPAR(I) HAS THE VALUE ...
C
C  INTEGERS.
C
C     ASSUME INTEGERS ARE REPRESENTED IN THE N-DIGIT, BASE-A FORM
C
C               SIGN ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )
C
C               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,N-1.
C
C     IPMPAR(1) = A, THE BASE.
C
C     IPMPAR(2) = N, THE NUMBER OF BASE-A DIGITS.
C
C     IPMPAR(3) = A**N - 1, THE LARGEST MAGNITUDE.
C
C  FLOATING-POINT NUMBERS.
C
C     IT IS ASSUMED THAT THE SINGLE AND DOUBLE PRECISION FLOATING
C     POINT ARITHMETICS HAVE THE SAME BASE, SAY B, AND THAT THE
C     NONZERO NUMBERS ARE REPRESENTED IN THE FORM
C
C               SIGN (B**E) * (X(1)/B + ... + X(M)/B**M)
C
C               WHERE X(I) = 0,1,...,B-1 FOR I=1,...,M,
C               X(1) .GE. 1, AND EMIN .LE. E .LE. EMAX.
C
C     IPMPAR(4) = B, THE BASE.
C
C  SINGLE-PRECISION
C
C     IPMPAR(5) = M, THE NUMBER OF BASE-B DIGITS.
C
C     IPMPAR(6) = EMIN, THE SMALLEST EXPONENT E.
C
C     IPMPAR(7) = EMAX, THE LARGEST EXPONENT E.
C
C  DOUBLE-PRECISION
C
C     IPMPAR(8) = M, THE NUMBER OF BASE-B DIGITS.
C
C     IPMPAR(9) = EMIN, THE SMALLEST EXPONENT E.
C
C     IPMPAR(10) = EMAX, THE LARGEST EXPONENT E.
C
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER i
C     ..
C     .. Local Arrays ..
      INTEGER imach(10)
C     ..
C     .. Data statements ..
C
C
C     MACHINE CONSTANTS FOR THE HP 9000.
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -126 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR THE IBM PC.
C
C      DATA imach(1)/2/
C      DATA imach(2)/31/
C      DATA imach(3)/2147483647/
C      DATA imach(4)/2/
C      DATA imach(5)/24/
C      DATA imach(6)/-125/
C      DATA imach(7)/128/
C      DATA imach(8)/53/
C      DATA imach(9)/-1021/
C      DATA imach(10)/1024/
C
C     MACHINE CONSTANTS FOR THE MACINTOSH II - ABSOFT
C     MACFORTRAN II.
C
C     DATA IMACH( 1) /     2 /
C     DATA IMACH( 2) /    31 /
C     DATA IMACH( 3) / 2147483647 /
C     DATA IMACH( 4) /     2 /
C     DATA IMACH( 5) /    24 /
C     DATA IMACH( 6) /  -125 /
C     DATA IMACH( 7) /   128 /
C     DATA IMACH( 8) /    53 /
C     DATA IMACH( 9) / -1021 /
C     DATA IMACH(10) /  1024 /
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
C     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
C     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
C
      imach(1)=2
      imach(2)=31
      imach(3)=2147483647
      imach(4)=2
      imach(5)=24
      imach(6)=-125
      imach(7)=128
      imach(8)=53
      imach(9)=-1021
      imach(10)=1024
C
C
      ipmpar = imach(i)
      RETURN

      END FUNCTION


      DOUBLE PRECISION FUNCTION spmpar(i)
C-----------------------------------------------------------------------
C
C     SPMPAR PROVIDES THE SINGLE PRECISION MACHINE CONSTANTS FOR
C     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
C     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
C     SINGLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
C     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN
C
C        SPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,
C
C        SPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,
C
C        SPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.
C
C-----------------------------------------------------------------------
C     WRITTEN BY
C        ALFRED H. MORRIS, JR.
C        NAVAL SURFACE WARFARE CENTER
C        DAHLGREN VIRGINIA
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C     MODIFIED BY BARRY W. BROWN TO RETURN DOUBLE PRECISION MACHINE
C     CONSTANTS FOR THE COMPUTER BEING USED.  THIS MODIFICATION WAS
C     MADE AS PART OF CONVERTING BRATIO TO DOUBLE PRECISION
C-----------------------------------------------------------------------
C     .. Scalar Arguments ..
      INTEGER i
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION b,binv,bm1,one,w,z
      INTEGER emax,emin,ibeta,m
C     ..
C     .. External Functions ..
C     INTEGER ipmpar
C     EXTERNAL ipmpar
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dble
C     ..
C     .. Executable Statements ..
C
      IF (i.GT.1) GO TO 10
      b = ipmpar(4)
      m = ipmpar(8)
      spmpar = b** (1-m)
      RETURN
C
   10 IF (i.GT.2) GO TO 20
      b = ipmpar(4)
      emin = ipmpar(9)
      one = dble(1)
      binv = one/b
      w = b** (emin+2)
      spmpar = ((w*binv)*binv)*binv
      RETURN
C
   20 ibeta = ipmpar(4)
      m = ipmpar(8)
      emax = ipmpar(10)
C
      b = ibeta
      bm1 = ibeta - 1
      one = dble(1)
      z = b** (m-1)
      w = ((z-one)*b+bm1)/ (b*z)
C
      z = b** (emax-2)
      spmpar = ((w*z)*b)*b
      RETURN

      END FUNCTION

      DOUBLE PRECISION FUNCTION stvaln(p)
C
C**********************************************************************
C
C     DOUBLE PRECISION FUNCTION STVALN(P)
C                    STarting VALue for Neton-Raphon
C                calculation of Normal distribution Inverse
C
C
C                              Function
C
C
C     Returns X  such that CUMNOR(X)  =   P,  i.e., the  integral from -
C     infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P
C
C
C                              Arguments
C
C
C     P --> The probability whose normal deviate is sought.
C                    P is DOUBLE PRECISION
C
C
C                              Method
C
C
C     The  rational   function   on  page 95    of Kennedy  and  Gentle,
C     Statistical Computing, Marcel Dekker, NY , 1980.
C
C**********************************************************************
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION p
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION sign,y,z
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION xden(5),xnum(5)
C     ..
C     .. External Functions ..
C     DOUBLE PRECISION devlpl
C     EXTERNAL devlpl
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC dble,log,sqrt
C     ..
C     .. Data statements ..
      DATA xnum/-0.322232431088D0,-1.000000000000D0,-0.342242088547D0,
     +     -0.204231210245D-1,-0.453642210148D-4/
      DATA xden/0.993484626060D-1,0.588581570495D0,0.531103462366D0,
     +     0.103537752850D0,0.38560700634D-2/
C     ..
C     .. Executable Statements ..
      IF (.NOT. (p.LE.0.5D0)) GO TO 10
      sign = -1.0D0
      z = p
      GO TO 20

   10 sign = 1.0D0
      z = 1.0D0 - p
   20 y = sqrt(-2.0D0*log(z))
      stvaln = y + devlpl(xnum,5,y)/devlpl(xden,5,y)
      stvaln = sign*stvaln
      RETURN

      END FUNCTION



      end module

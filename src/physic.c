#include <tgmath.h>
#include <stdio.h>
#include "def.h"
#include "bunch.h"
#include "physic.h"

#define N            19
#define TM           6

#define SQRT_PI      1.7724538509055158819
#define SQRT_POW_1_5 5.5683279968317078712

#define TM_POW_2     36.0
#define A0           0.29540897515091929515
#define FACT_AN      -0.068538919452009433586

/*
static inline
_Complex double sin_na(const unsigned n, const _Complex double a)
{
  if(n == 0)      return 0.0;
  else if(n == 1) return sin(a);
  else            return 2.0*cos(a)*sin_na(n-1, a) - sin_na(n-2, a);
}
*/

static inline
void sin_na(const unsigned n, const _Complex double a, _Complex double * res)
{
/* TODO fix compilation
  const _Complex double cosa = cos(a);
  res[0] = 0.0;
  res[1] = sin(a);
  unsigned i;
  for(i = 2; i <= n; i++)
  {
    res[i] = 2.0*cosa*res[i-1] - res[i-2];
  }
*/
}

int complex_errorfunc_wofz_abrarov_quine_approx9(const double x, const double y, double w[2])
{
/*
  const double a0 = SQRT_PI / TM;
  const double fact_an = -1.0*M_PI*M_PI/(4.0*TM_POW_2);
*/
  
//  const _Complex double exp1 = exp(-pow(x + 1.0i*y, 2));

  const _Complex double exp1 = (exp(-x*x)*exp(-1.0i*x*y)*exp(-y*y));
  
  _Complex double factor = 1.0;
  
  _Complex double sum = 0.0;
  _Complex double sin_na_res[N+1];
  sin_na(N, (M_PI/TM) * (- 1.0i * x + y), sin_na_res);
  unsigned n;
  for(n = 1; n <= N; n++)
  {
   // sum += (exp(FACT_AN*n*n)/n) * sin_na_res[n];
  }
  factor -= A0 * sum * 2.0 * TM / SQRT_POW_1_5;
  
  factor -= A0*(-1.0i * x + y) / SQRT_PI;
  
  _Complex double w_tmp = exp1 * factor;
  w[0] = creal(w_tmp);
  w[1] = cimag(w_tmp);
  
  return 1;
}


int complex_errorfunc_wofz_algo680(xi, yi, uu)
  const double xi, yi;
  double uu[2];
{
  int    i, j, n;
  int    iflag;
  int    kapn, nu, np1;
  int    ilogicA, ilogicB;
  const double factor   = 1.12837916709551257388;
  const double rmaxreal = 0.5e154;
  const double rmaxexp  = 708.503061461606;
  const double rmaxgoni = 3.53711887601422e08;
  double xabs, yabs, x, y;
  double qrho, xabsq, xquad, yquad;
  double xsum, ysum, xaux, u1, v1, u2, v2, w1, daux;
  double h, h2, qlamda, rx, ry, sx, sy;
  double tx, ty, c;
/*********************************************************************
C      ALGORITHM 680, COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 16, NO. 1, PP. 47.
C
C  GIVEN A COMPLEX NUMBER z = (xi, yi), THIS SUBROUTINE COMPUTES C  THE VALUE OF THE FADDEEVA-FUNCTION w(z) = exp(-z**2)*erfc(-i*z), C  WHERE ERFC 

IS THE COMPLEX COMPLEMENTARY 

ERROR-FUNCTION AND i C  MEANS SQRT(-1). C  THE ACCURACY OF THE ALGORITHM FOR Z IN THE 1ST AND 2ND QUADRANT C  IS 14 SIGNIFICANT DIGITS; IN THE 3RD 

AND 4TH IT IS 13 SIGNIFICANT C  DIGITS 

OUTSIDE A CIRCULAR REGION WITH RADIUS 0.126 AROUND A ZERO C  OF THE FUNCTION. C  ALL REAL VARIABLES IN THE PROGRAM ARE DOUBLE PRECISION. C C C  THE 

CODE CONTAINS A FEW COMPILER-DEPENDENT 

PARAMETERS :
C     RMAXREAL = THE MAXIMUM VALUE OF RMAXREAL EQUALS THE ROOT OF
C                RMAX = THE LARGEST NUMBER WHICH CAN STILL BE
C                IMPLEMENTED ON THE COMPUTER IN DOUBLE PRECISION
C                FLOATING-POINT ARITHMETIC
C     RMAXEXP  = LN(RMAX) - LN(2)
C     RMAXGONI = THE LARGEST POSSIBLE ARGUMENT OF A DOUBLE PRECISION
C                GONIOMETRIC FUNCTION (DCOS, DSIN, ...)
C  THE REASON WHY THESE PARAMETERS ARE NEEDED AS THEY ARE DEFINED WILL C  BE EXPLAINED IN THE CODE BY MEANS OF COMMENTS C C C  PARAMETER LIST
C     xi     = REAL      PART OF z
C     yi     = IMAGINARY PART OF z
C     u      = REAL      PART OF w(z)
C     v      = IMAGINARY PART OF w(z)
C     iflag  = AN ERROR FLAG INDICATING WHETHER OVERFLOW WILL
C              OCCUR OR NOT; TYPE LOGICAL;
C              THE VALUES OF THIS VARIABLE HAVE THE FOLLOWING
C              MEANING :
C              FLAG=.FALSE. : NO ERROR CONDITION
C              FLAG=.TRUE.  : OVERFLOW WILL OCCUR, THE ROUTINE
C                             BECOMES INACTIVE
C  xi, yi      ARE THE INPUT-PARAMETERS
C  u, u, iflag ARE THE OUTPUT-PARAMETERS
C
C  FURTHERMORE THE PARAMETER FACTOR EQUALS 2/SQRT(PI)
C
C  THE ROUTINE IS NOT UNDERFLOW-PROTECTED BUT ANY VARIABLE CAN BE C  PUT TO 0 UPON UNDERFLOW; C C  REFERENCE - GPM POPPE, CMJ WIJERS; MORE EFFICIENT 

COMPUTATION OF C  THE COMPLEX 

ERROR-FUNCTION, ACM TRANS. MATH. SOFTWARE. C
***************************************************************************
* <mbtrack>                                 Last modified: August 2007
***************************************************************************/
  
      iflag   = 1;  /*** if iflag = 1, it means that there is no error ****/
      ilogicA = 0;

      xabs = (fabs)(xi);
      yabs = (fabs)(yi);
      x    = xabs/6.3;
      y    = yabs/4.4;

/*****************************************************************************
C
C     THE FOLLOWING IF-STATEMENT PROTECTS
C     QRHO = (X**2 + Y**2) AGAINST OVERFLOW
C ******************************************************************************/
      if(xabs > rmaxreal || yabs > rmaxreal) goto pt_100;

      qrho = x*x + y*y;

      xabsq = xabs*xabs;
      xquad = xabsq - yabs*yabs;
      yquad = 2.0*xabs*yabs;
 
      if(qrho < 0.085264) {
         ilogicA = 1;
        

         /********************************************************************
         C  IF (QRHO.LT.0.085264D0) THEN THE FADDEEVA-FUNCTION IS EVALUATED
         C  USING A POWER-SERIES (ABRAMOWITZ/STEGUN, EQUATION (7.1.5), P.297)
         C  N IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
         C  ACCURACY
         *********************************************************************/

         qrho  = (1.0 - 0.85*y)*sqrt(qrho);
         n     = (int)(6.0 + 72.0*qrho);
         j     = 2*n + 1;
         xsum  = 1.0/((double) j);
         ysum  = 0.0;
         for(i=n; i>=1; i--) {
            j = j - 2;
            xaux = (xsum*xquad - ysum*yquad)/((double) i);
            ysum = (xsum*yquad + ysum*xquad)/((double) i);
            xsum = xaux + 1.0/((double) j);
         }
         u1   = -factor*(xsum*yabs + ysum*xabs) + 1.0;
         v1   =  factor*(xsum*xabs - ysum*yabs);
         daux = exp(-xquad);
         u2   = daux*cos(yquad);
         v2   =-daux*sin(yquad);

         uu[0]   = u1*u2 - v1*v2;
         uu[1]   = u1*v2 + v1*u2;

      } else {

         /**********************************************************************
         C  IF (QRHO.GT.1.O) THEN W(Z) IS EVALUATED USING THE LAPLACE
         C  CONTINUED FRACTION
         C  NU IS THE MINIMUM NUMBER OF TERMS NEEDED TO OBTAIN THE REQUIRED
         C  ACCURACY
         C
         C  IF ((QRHO.GT.0.085264D0).AND.(QRHO.LT.1.0)) THEN W(Z) IS EVALUATED
         C  BY A TRUNCATED TAYLOR EXPANSION, WHERE THE LAPLACE CONTINUED FRACTION
         C  IS USED TO CALCULATE THE DERIVATIVES OF W(Z)
         C  KAPN IS THE MINIMUM NUMBER OF TERMS IN THE TAYLOR EXPANSION NEEDED
         C  TO OBTAIN THE REQUIRED ACCURACY
         C  NU IS THE MINIMUM NUMBER OF TERMS OF THE CONTINUED FRACTION NEEDED
         C  TO CALCULATE THE DERIVATIVES WITH THE REQUIRED ACCURACY
         ***********************************************************************/

        if(qrho > 1.0) {
           h    = 0.0;
           kapn = 0;
           qrho = sqrt(qrho);
           nu   = (int) (3.0 + 1442.0/(26.0*qrho + 77.0));
        } else {
           qrho = (1.0 - y)*sqrt(1.0 - qrho);
           h    = 1.88*qrho;
           h2   = 2.0*h;
           kapn = (int)(7.0  + 34.0*qrho);
           nu   = (int)(16.0 + 26.0*qrho);
        }

        if(h > 0.0) {
	   ilogicB = 1; 
	   qlamda  = pow(h2, kapn);
	} else {
	   ilogicB = 0;
	}     

        rx = 0.0;
        ry = 0.0;
        sx = 0.0;
        sy = 0.0;

        for(n=nu; n>=0; n--) {
           np1 = n + 1;
           tx  = yabs + h + np1*rx;
           ty  = xabs - np1*ry;
           c   = 0.5/(tx*tx + ty*ty);
           rx  = c*tx;
           ry  = c*ty;
           if(ilogicB && (n <= kapn)) {
              tx = qlamda + sx;
              sx = rx*tx - ry*sy;
              sy = ry*tx + rx*sy;
              qlamda = qlamda/h2;
           }
        }

        if(h == 0.0) {
           uu[0] = factor*rx;
           uu[1] = factor*ry;
        } else {
           uu[0] = factor*sx;
           uu[1] = factor*sy;
        }

        if(yabs == 0.0) uu[0] = exp(-xabs*xabs);

     }


/*****************************************************************************
C
C  EVALUATION OF W(Z) IN THE OTHER QUADRANTS
C ******************************************************************************/
      
      if(yi < 0.0) {
           
         if(ilogicA) {
            u2 = 2.0*u2;
            v2 = 2.0*v2;
         } else {
            xquad = -xquad;
            /******************************************************************
            C         THE FOLLOWING IF-STATEMENT PROTECTS 2*EXP(-Z**2)
            C         AGAINST OVERFLOW
            *******************************************************************/
            if((yquad > rmaxgoni) || (xquad > rmaxexp)) goto pt_100;
            w1 = 2.0*exp(xquad);
            u2 = w1*cos(yquad);
            v2 =-w1*sin(yquad);
         }


         uu[0] = u2 - uu[0];
         uu[1] = v2 - uu[1];
         if(xi > 0.0) uu[1] = -uu[1];
      } else {
         if(xi < 0.0) uu[1] = -uu[1];
      }

      return iflag;

      pt_100:;
      iflag =-1; 
/*********************************************************************************/
      return iflag;
}

void reglin(const double * x, const double * y, const unsigned n, double * a, double * b)
{
  double xavr = 0.0;
  double yavr = 0.0;
  unsigned i;
  for(i = 0; i < n; i++)
  {
    xavr += x[i];
    yavr += y[i];
  }
  xavr = xavr / ((double) n);
  yavr = yavr / ((double) n);
  
  double k1 = 0.0;
  double k2 = 0.0;
  double tmp;
  for(i = 0; i < n; i++)
  {
    k1 += (x[i] - xavr) * (y[i] - yavr);
    tmp = x[i] - xavr;
    k2 += tmp*tmp;
  }
  *a = k1/k2;
  *b = yavr - (k1/k2) * xavr;
}

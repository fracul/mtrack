/**
 * @file
 * Physic module
 */

#ifndef MBTRACK_PHYSIC_H
#define MBTRACK_PHYSIC_H

/**
 * We use algorithm 680 to computer complex error function
 */
#define complex_errorfunc_wofz(x, y, u) complex_errorfunc_wofz_algo680(x, y, u)

/**
 * Complex error function
 *
 * Algorithm approximation 9 by Abrarov and Quine
 * Publised in Applied Mathematics and Computation, No. 218, Pages 1894-1902
 *
 * @param x
 *   Real part of z
 * @param y
 *   Imaginary part of z
 * @param w
 *   Result w(z)
 */
int complex_errorfunc_wofz_abrarov_quine_approx9(const double x, const double y, double w[2]);


/**
 * GIVEN A COMPLEX NUMBER z = (xi, yi), THIS SUBROUTINE COMPUTES THE VALUE
 * OF THE FADDEEVA-FUNCTION w(z) = exp(-z**2)*erfc(-i*z), WHERE ERFC IS THE
 * COMPLEX COMPLEMENTARY ERROR-FUNCTION AND i MEANS SQRT(-1).
 *
 * @param xi
 *   Real part of z
 * @param yi
 *   Imaginary part of z
 * @param uu
 *   Result w(z)
 *
 * ALGORITHM 680, COLLECTED ALGORITHMS FROM ACM.
 * THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
 * VOL. 16, NO. 1, PP. 47.
 */
int complex_errorfunc_wofz_algo680(const double xi, const double yi, double uu[2]);

/**
 * Linear regression, using least squares method
 *
 * @param x
 *   Array of n double values
 * @param y
 *   Array of n double values
 * @param n
 *   Array size
 * @param a
 *   Coefficient a where y(i) ~ a * x(i) + b
 * @param b
 *   Coefficient b where y(i) ~ a * x(i) + b
 */
void reglin(const double * x, const double * y, const unsigned n, double * a, double * b);

#endif /* MBTRACK_PHYSIC_H */

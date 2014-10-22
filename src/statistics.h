/**
 * @file
 * Module for statistics
 */

#ifndef STATISTICS_H
#define STATISTICS_H

#include <stdio.h>
#include "bunch.h"

inline vector_t
vector_new(const double xtau, const double x, const double z);

inline vector_t
vector_add(const vector_t a, const vector_t b);

inline vector_t
vector_mul_scalar(const vector_t a, const double k);

inline vector_t
vector_square(const vector_t a);

/**
 * Computes the variance from sum(x_i) and sum(x_i^2) over N turns:
 * Var = 1/(N-1) * {sum(x_i^2) - 1/N * sum(x_i)^2}
 * @param a Sum of position values: sum(x_i)
 * @param b_sqr Sum of square of positon values: sum(x_i)^2
 * @param N Sum over N values
 */
vector_t
vector_variance(const vector_t a, const vector_t b_sqr, const long int N);

/**
 * Compute moments of a weak bunch
 * See http://reference.wolfram.com/mathematica/ref/Moment.html
 */
void moments(const weak_bunch_t * bunch, const plane_t plane,
             unsigned int orders, double * moments);

/**
 * Compute moments > 2 of a weak bunch, takes 1 and 2 from bunch stats
 * See http://reference.wolfram.com/mathematica/ref/Moment.html
 */
void moments_quick(const weak_bunch_t * bunch, const plane_t plane,
             unsigned int orders, double * moments);

/**
 * Compute central moments of a weak bunch
 * See http://reference.wolfram.com/mathematica/ref/CentralMoment.html
 */
void central_moments(const weak_bunch_t * bunch, const plane_t plane,
                     const double mean, unsigned int orders, double * moments);

/**
 * Compute root mean square
 *
 * @param in
 *    Input double array
 * @param n
 *    Input array size
 * @param avr
 *    Average of input
 */
double rms(const double * in, const size_t n, const double avr);

/**
 * Compute weighted root mean square
 *
 * @param in
 *    Input double array
 * @param weights
 *    Weights of input
 * @param n
 *    Input array size
 * @param avr
 *    Average of input (WEIGHTED average)
 */
double weighted_rms(const double * in, const int * weights, const size_t n, const double avr);

/**
 * Generate an histogram from an array of double values
 *
 * @param input
 *    Input double array
 * @param n
 *    Size of the input array
 * @param middle
 *    Middle of the histogram
 * @param m
 *    Number of bins in the histogram
 * @param width
 *    Width of one bin
 * @param hist
 *    Generated histogram, array of size m (result variable)
 */
int histogram_generate(const double * input, const unsigned n,
                       const double middle, const unsigned m, const double width,
                       unsigned * hist);

/**
 * Generated a 2D histogram from two arrays of doubles values
 */
int histogram_2d_generate(const double * input1, const double * input2, const unsigned n,
                          const double middle1, const unsigned m1, const double width1,
                          const double middle2, const unsigned m2, const double width2,
                          unsigned * hist);

/**
 * Generate two 1D histograms, and one 2D histogram, and write results
 * in two separate files
 *
 * @param input1
 *    Input for X axis values
 * @param input2
 *    Input for Y axis values
 * @param n
 *    Size of the input arrays
 * @param m
 *    Number of histogram bins along 1 axis
 * @param filename_1d
 *    Filename for 1D histograms writing
 * @param filename_2d
 *    Filename for 2D histograms writing
 */
int histogram_write_1d_2d(const double * input1, const double * input2,
                          const unsigned n, const unsigned m,
                          const double middle1, const double width1,
                          const double middle2, const double width2,
                          const char filename_1d[FILENAME_MAX],
                          const char filename_2d[FILENAME_MAX]);

/**
 * Generate histograms for each of 6 input, print them in a file
 *
 * @param input1
 *    Input double array #1
 * @param middle1
 *    Middle of the histogram #1
 * @param sigma1
 *    Standard deviation of input #1
 * @param n
 *    Input array size
 * @param filename
 *    Filename for histograms output
 */
int particules_histograms(const double * input1, const double middle1, const double sigma1,
                          const double * input2, const double middle2, const double sigma2,
                          const double * input3, const double middle3, const double sigma3,
                          const double * input4, const double middle4, const double sigma4,
                          const double * input5, const double middle5, const double sigma5,
                          const double * input6, const double middle6, const double sigma6,
                          const unsigned n,
                          const char filename[FILENAME_MAX]);


/**
 * Calculate average and RMS position for one plane
 *
 * @param bunch Bunch (result variable)
 * @param plane Plane to consider
 *
 */
void
weak_bunch_mean_rms(weak_bunch_t * bunch, unsigned int plane);

/**
 * Calculate actual bunch statistics (average and RMS position)
 *
 * @param bunch Bunch (result variable)
 * @param track Tracking configuration to indicate tracked planes
 */
void
weak_bunch_calc_statistics(weak_bunch_t * bunch, const tracking_t *track);


/**
 * Calculate actual ampinv
 *
 * @param bunch Bunch (result variable)
 * @param plane Plane to consider
 */
void
weak_bunch_calc_ampinv(weak_bunch_t * bunch,  const tracking_t * track);

/**
 * Add actual stast to sums
 *
 * @param allstats Summed stats (result variable)
 * @param addstats Stats to add
 */
void
weak_bunch_add_statistics(bunch_stats_t * allstats, const bunch_stats_t * addstats);


/**
 * Update bunch statistics (average over turns and or bunches)
 *
 * @param bunch Bunch (result variable)
 * @param Nstat averaging over Nstat data-sets
 */
void
weak_bunch_update_statistics(bunch_stats_t * bstats, const long int Nstat);

/**
 * Reset statistics' sums to zero
 *
 * @param bstats Bunch statistics (result variable)
 */
void
weak_bunch_reset_statistics(bunch_stats_t * bstats);

#endif /* STATISTICS_H */

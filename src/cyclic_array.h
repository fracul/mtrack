#ifndef MBTRACK_CYCLIC_ARRAY_H
#define MBTRACK_CYCLIC_ARRAY_H



#include <stdlib.h>
#include <stdbool.h>

/**
 * \brief A two dimension double array (double[n][m]), cyclic around n.
 */
typedef struct cyclic_array
{
  double * arr; /**< (n x m)-Array containing distribution history, cyclic around n */
  unsigned n; /**< Bunch statistics history for the last n buckets: n = (turns considered) * Nbunches */
  unsigned m; /**< Moments stored up to mth order, m = MOMENTS_ORDERS */
}
cyclic_array_t;

/**
 * Create a new cyclic array (allocate memory)
 */
bool cyclic_array_create(cyclic_array_t * carray, unsigned n, unsigned m);

/**
 * Destroy a cyclic array (free memory)
 */
void cyclic_array_destroy(cyclic_array_t * carray);

/**
 * Get the 1D array from array[n][]
 * 
 * _________________________________
 *      o->      o->  ....
 * ---------------------------------
 *      |        |
 * n:   kb      kb+1 ...
 * 
 */
const double * cyclic_array_get(int n, const cyclic_array_t * carray);

/**
 * Access the 1D array from array[n][]
 * 
 */
double * cyclic_array_access(int n, cyclic_array_t * carray);

/**
 * Set the 1d array array[n][] values
 */
void cyclic_array_set(unsigned n, double * inp, cyclic_array_t * carray);

/**
 * Fill the n 1D array with the input 1D array
 * WARNING : not sure if cyclic_array_fill works
 */
void cyclic_array_fill(double * inp, cyclic_array_t * carray, unsigned int kb, unsigned int Nturn);

#endif /* MBTRACK_CYCLIC_ARRAY_H */

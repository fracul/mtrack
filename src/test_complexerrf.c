#define _GNU_SOURCE
#include <fenv.h>
#include <signal.h>

#include <math.h>
#include <stdlib.h>
#include "fbi.h"
#include "test.h"

/* Global variables */
ring_t ring;
tracking_t track;
grid_t fbii_grid;

/* #define complex_errorfunc_wofz(x, y, u) complex_errorfunc_wofz_abrarov_quine_approx9(x, y, u) */

int main(int argc, char ** argv)
{
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_UNDERFLOW | FE_OVERFLOW);
  int n = 1;
  if(argc == 2)
  {
    n = atoi(argv[1]);
  }
  
  test_init(true);

  int i;
  for(i = 0; i < n; i++)
  {
    double w[2];

    complex_errorfunc_wofz(0.0, 0.0, w);
    ASSERT_EQUAL_FLOAT(w[0], 1.0, 1.0e-8);
    ASSERT_EQUAL_FLOAT(w[1], 0.0, 1.0e-8);
    
    complex_errorfunc_wofz(0.01, 0.01, w);
    ASSERT_EQUAL_FLOAT(w[0], 0.9887176929549547, 1.0e-8);
    ASSERT_EQUAL_FLOAT(w[1], 0.01108529605747726, 1.0e-8);
    
    complex_errorfunc_wofz(0.5, 0.5, w);
    ASSERT_EQUAL_FLOAT(w[0], 0.5331567079121750, 1.0e-8);
    ASSERT_EQUAL_FLOAT(w[1], 0.2304882313844584, 1.0e-8);
    
    complex_errorfunc_wofz(0.5, 1.0, w);
    ASSERT_EQUAL_FLOAT(w[0], 0.391234, 1.0e-5);
    ASSERT_EQUAL_FLOAT(w[1], 0.127202, 1.0e-5);
    
    complex_errorfunc_wofz(1.0, 1.0, w);
    ASSERT_EQUAL_FLOAT(w[0], 0.3047442052569128, 1.0e-8);
    ASSERT_EQUAL_FLOAT(w[1], 0.2082189382028316, 1.0e-8);
    
    complex_errorfunc_wofz(2.5, 2.5, w);
    ASSERT_EQUAL_FLOAT(w[0], 0.1167371250446503, 1.0e-8);
    ASSERT_EQUAL_FLOAT(w[1], 0.1079085859964814, 1.0e-8);
    
    complex_errorfunc_wofz(15.0, 15.0, w);
    ASSERT_EQUAL_FLOAT(w[0], 0.01882714532513676, 1.0e-8);
    ASSERT_EQUAL_FLOAT(w[1], 0.01878535427799565, 1.0e-8);

  }
  
  return test_backend();
}

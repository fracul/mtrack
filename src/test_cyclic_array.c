#include <stdlib.h>
#include <math.h>
#include "test.h"
#include "cyclic_array.h"

int main()
{
  cyclic_array_t carray;
  unsigned int in, im;
  double * input = NULL;

  test_init(true);

  ASSERT_TRUE(cyclic_array_create(&carray, 123, 6));
  ASSERT_EQUAL(carray.n, 123);
  ASSERT_EQUAL(carray.m, 6);

  input = (double *) malloc(carray.m * sizeof(double));
  ASSERT_NOT_EQUAL(input, NULL);

  /* Test random set */
  srand(0x0123);
  for(in = 0; in < carray.n; in++)
  {
    for(im = 0; im < carray.m; im++)
    {
      input[im] = rand();
    }
    cyclic_array_set(in, input, &carray);
  }
  
  const double * a1 = cyclic_array_get(3, &carray);
  const double * a2 = cyclic_array_get(3 + 2 * carray.n, &carray);
  const double * a3 = cyclic_array_get(3 - 7 * carray.n, &carray);
  for(im = 0; im < carray.m; im++)
  {
    ASSERT_EQUAL(a1[im], a2[im]);
    ASSERT_EQUAL(a1[im], a3[im]);
  }

  /* Test random filling */
  srand(0x0654);
  for(im = 0; im < carray.m; im++)
  {
    input[im] = rand();
  }
  cyclic_array_fill(input, &carray);
  
  for(in = 0; in < carray.n; in++)
  {
    const double * a = cyclic_array_get(in, &carray);
    
    for(im = 0; im < carray.m; im++)
    {
      ASSERT_EQUAL(a[im], input[im]);
    }
  }

  cyclic_array_destroy(&carray);

  return test_backend();
}

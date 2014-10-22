#include <stdio.h>
#include "test.h"
#include "physic.h"

#define N 43

int main()
{
  double x[N];
  double y[N];
  double a, b;
  double aref = 2.5;
  double bref = 4.0;
  int i;
  
  test_init(true);
  
  for(i = 0; i < N; i++)
  {
    x[i] = i;
    y[i] = aref * x[i] + bref;
  }
  
  reglin(x, y, N, &a, &b);
    
  ASSERT_EQUAL_FLOAT(a, aref, 1.0e-8);
  ASSERT_EQUAL_FLOAT(b, bref, 1.0e-8);
  
  return test_backend();
}

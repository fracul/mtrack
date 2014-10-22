#include <stdio.h>
#include "test.h"
#include "statistics.h"

#define N 1000
#define M 10

int main()
{
  double xx[N];
  unsigned i = 0;
  double sum = 0.0;
  for(i = 0; i < N; i++)
  {
    xx[i] = i;
    sum += xx[i];
  }
  double average = sum/((double) N);
  unsigned hist[M];
  if(histogram_generate(xx, N, average, M, ((double) N)/((double) M), hist) < 0)
  {
    fprintf(stderr, "histogram_generate failed\n");
    return -1;
  }
  else
  {
    test_init(true);
    
    for(i = 0; i < M; i++)
    {
      ASSERT_EQUAL_FLOAT(hist[i], ((double) N)/((double) M), 1.0e-8);
    }
    
    return test_backend();
  }
}

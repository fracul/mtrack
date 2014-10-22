#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include "cyclic_array.h"

bool cyclic_array_create(cyclic_array_t * carray, unsigned n, unsigned m)
{
  carray->arr = (double *) calloc(n * m, sizeof(double));
  carray->n = n;
  carray->m = m;
  return carray->arr != NULL;
}

void cyclic_array_destroy(cyclic_array_t * carray)
{
  free(carray->arr);
  carray->n = 0;
  carray->m = 0;
}

const double * cyclic_array_get(int in, const cyclic_array_t * carray)
{
  const int n = (int) carray->n;
  const int in_mod = in < 0 ? (in%n + n)%n : in%n;
  return &(carray->arr[in_mod * carray->m]);
}


double * cyclic_array_access(int in, cyclic_array_t * carray)
{
  const int n = (int) carray->n;
  const int in_mod = in < 0 ? (in%n + n)%n : in%n;
  return &(carray->arr[in_mod * carray->m]);
}

void cyclic_array_set(unsigned n, double * inp, cyclic_array_t * carray)
{
  unsigned key = (n%carray->n)*carray->m;
  unsigned im;
  for(im = 0; im < carray->m; im++)
    carray->arr[key + im] = inp[im];
}


/* Fill cyclic_array with start moments */
void cyclic_array_fill(double * inp, cyclic_array_t * carray, unsigned int kb, unsigned int Nturn)
{
  /* printf( "\n WARNING : not sure if cyclic_array_fill works !\n"); */
  unsigned it, im;
  unsigned int h = carray->n / Nturn; 
  for(it = 1; it < Nturn; it++)
  {
      for(im = 0; im < carray->m; im++)
      {
        carray->arr[kb + it*h*carray->m + im] = inp[im];

      }
  }
}

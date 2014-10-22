#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "test.h"
#include "input.h"
#include "physic.h"
#include "def.h"
#include "types.h"

double tau_to_minus1(const plane_t plane, const ring_t ring)
{
  const double beta = ring.beta1[plane];
  const double w0 = ring.w0;
  const double I = ring.Ib*FMILLI; /* mA -> A */
  const double E = ring.E0*FGIGA; /* GeV -> eV */
  const double R = ring.R;
  const double beff3 = pow(ring.beff[2], 3);
  const double rho = 2.8e-08; /* TODO read from input */
  double Q = 0.0;
  if(plane == HOR)      Q = ring.QH0;
  else if(plane == VER) Q = ring.QV0;
  const double delta_Qbeta = fmod(Q, 1.0); /* fractional part of Q */
  
  return (beta*w0*I / (4.0*M_PI*E))
         * (R/beff3)
         * sqrt(2.0 * rho / ((1.0-delta_Qbeta)*w0*EPSILON0) );
}

int main(int argc, char ** argv)
{
  if(argc != 3)
  {
  fprintf(stderr, "ERROR ! \n\nUsage:\ntest_rw input_file output_file\n");
  return -1;
  }
  
  /* Read inputs from file */
  ring_t ring;
  tracking_t track;
  bunch_macroparticle_model_t macrop_model;
  e_beam_t ebeam;
  
  if(!read_input(argv[1], &ring, &track, &macrop_model, &(ebeam.distrib)))
  {
    fprintf(stderr, "ERROR: read_input failed\n");
    return -1;
  }
  
  FILE * fp = fopen(argv[2], "r");
  if(fp == NULL)
  {
    fprintf(stderr, "ERROR: failed opening %s\n", argv[2]);
  }
  
  /* Setup parameters */
  setup_ring_parameters(&ring);
  macrop_model_setup_parameters(ring, track, &macrop_model);
  
  double * x = (double *) malloc((1 + track.NrevTot/track.NrevMon) * sizeof(double));
  double * y = (double *) malloc((1 + track.NrevTot/track.NrevMon) * sizeof(double));
  double a;
  double b;
  int i = 0;
  test_init(true);
  
  while(!feof(fp))
  {
    int rev;
    double ampinv1;
    double ampinv2;
    if(fscanf(fp, "%d %lf %lf", &rev, &ampinv1, &ampinv2) == 3)
    {
      x[i] = rev*ring.Lc/C_LIGHT;
      y[i] = log(ampinv1);
      i++;
    }
  }
  fclose(fp);
  reglin(x, y, i, &a, &b);
  free(x);
  free(y);
  
  double taum1H = tau_to_minus1(HOR, ring)*2.0;
  double taum1V = tau_to_minus1(VER, ring)*2.0;
  
  printf("Analytic estimation\n  tau[HOR] = %g ms\n  tau[VER] = %g ms\n", FKILO*(1.0/taum1H), FKILO*(1.0/taum1V));
  printf("Simulation (%d measure points)\n  tau[?input?] = %g ms\n", i,  FKILO*(1.0/a));

  return test_backend();
}

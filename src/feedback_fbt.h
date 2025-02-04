#ifndef MBTRACK_FEEDBACK_FBT
#define MBTRACK_FEEDBACK_FBT

#include "cyclic_array.h"
#include "def.h"

typedef struct fbt
{
  int tap;
  int downsampling;
  double gain;
  double tune;
  double phase;
  double * coeffs_FIR;
  cyclic_array_t * offset_history;
  int filt_type;
  plane_t plane;
} fbt_feedback_t;

bool fbt_init(fbt_feedback_t * fbt);

double fbt_kick(fbt_feedback_t * fbt,int rev);

void fbt_calc_coeffs_Dimtel(fbt_feedback_t * fbt);

void fbt_calc_coeffs_Spring8(fbt_feedback_t * fbt);

void fbt_calc_coeffs_energy_sensing(fbt_feedback_t * fbt);
#endif

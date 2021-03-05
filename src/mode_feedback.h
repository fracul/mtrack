#ifndef MODE_FEEDBACK_H
#define MODE_FEEDBACK_H

#include "types.h"
#include "bunch.h"

/*typedef struct mode_feedback
{
  int mode;
  double gain;
  int diff_delay;
  double output_phase;
  int averaging;
} mode_feedback_t;*/

void
modefb_get_phidiff(mode_feedback_t * mfb,cyclic_array_t * all_moments,const int bpos, ring_t * ring, e_beam_t * ebeam, double * phasediff);

#endif // MODE_FEEDBACK_H

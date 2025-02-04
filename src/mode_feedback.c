#include <math.h>
#include "cyclic_array.h"
#include "types.h"
#include "bunch.h"
#include "mode_feedback.h"
#include "def.h"

void modefb_get_phidiff(mode_feedback_t * mfb,cyclic_array_t * all_moments,const int bpos, ring_t * ring, e_beam_t * ebeam, double * phasediff) {

  const int h = ring->Nharm;
  unsigned int n;
  unsigned int m;
  int k, in_mod;
  const double * moments;
  double mode_proj [2] = {0, 0};

  for (m=0; m<2; m++) {
    for (n=0; n<ring->Nharm*mfb->averaging; n++){
      k = bpos+n+m*ring->Nharm*mfb->diff_delay;
      in_mod = k < 0 ? (k%h + h)%h: k%h;
      if (ebeam->nfFill[in_mod]) {
	moments = cyclic_array_get(k, all_moments);
	mode_proj[m] += cos(2*M_PI*n*mfb->mode/ring->h)*moments[2];
      }
    }
  }

  *phasediff = mfb->gain*(mode_proj[0]-mode_proj[1])*2*M_PI*ring->frf/mfb->averaging/ring->h+mfb->output_phase;
}

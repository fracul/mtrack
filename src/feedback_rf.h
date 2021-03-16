#ifndef MBTRACK_FEEDBACK_RF
#define MBTRACK_FEEDBACK_RF

#include "cyclic_array.h"

typedef struct rffeedback
{
  double vrf_design;
  double phi0_design;
  int lr_resonator;
  int len_average;
  int count;
  cyclic_array_t * voltage_history;
} rf_feedback_t;

void
rffb_get_vrf_phi(rf_feedback_t * RfFeedback, double vreal, double vimag);

void
rffb_init(rf_feedback_t * RfFeedback);

void
rffb_calc_mean_voltage_phase(rf_feedback_t * RfFeedback, double * vrfbar, double * phibar);

#endif

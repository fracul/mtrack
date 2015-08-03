#include <math.h>
#include "feedback_rf.h"

void
rffb_get_vrf_phi(rf_feedback_t * RfFeedback,double vreal,double vimag) {
  double phi = atan2(vreal,vimag);
  double vrf = sqrt(vreal*vreal+vimag*vimag);
  double v_point[2] = {vrf,phi};
  cyclic_array_set(RfFeedback->count,v_point,RfFeedback->voltage_history);
  RfFeedback->count++;
}

void
rffb_init(rf_feedback_t * RfFeedback) {
  RfFeedback->voltage_history = (cyclic_array_t *) malloc(1*sizeof(cyclic_array_t));
  cyclic_array_create(RfFeedback->voltage_history,RfFeedback->len_average,2);
  RfFeedback->count = 0;
}

void
rffb_calc_mean_voltage_phase(rf_feedback_t * RfFeedback, double * vrfbar, double * phibar) {
  //double * vrf = cyclic_array_get(0,RfFeedback->voltage_history);
  //double * phi = cyclic_array_get(1,RfFeedback->voltage_history);
  double vrf_sum = 0;
  double phi_sum = 0;
  int m;
  for (m=0; m<RfFeedback->len_average; m++) {
    double * vrfphi = cyclic_array_get(m,RfFeedback->voltage_history);
    vrf_sum += vrfphi[0];
    phi_sum += vrfphi[1];
  }
  *vrfbar = vrf_sum/RfFeedback->len_average;
  *phibar = phi_sum/RfFeedback->len_average;
}

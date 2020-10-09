#include "feedback_fbt.h"
#include "coeffs_Sp8.h"
#include <math.h>

bool
fbt_init(fbt_feedback_t * fbt) {
  fbt->offset_history = (cyclic_array_t *)  malloc(1*sizeof(cyclic_array_t));
  cyclic_array_create(fbt->offset_history,fbt->downsampling*fbt->tap,1);
  //calculate filter coefficients
  fbt->coeffs_FIR = (double *) malloc(fbt->tap*sizeof(double));
  if (fbt->filt_type==1) {
    fbt_calc_coeffs_Dimtel(fbt);
  }
  else if (fbt->filt_type==2) {
    fbt_calc_coeffs_Spring8(fbt);
  }
  else if (fbt->filt_type==3) {
    fbt_calc_coeffs_energy_sensing(fbt);
  }
  else return false;

  return true;
}

void
fbt_calc_coeffs_Dimtel(fbt_feedback_t * fbt) {
  int t;
  double norm = 0;
  for (t=0; t<fbt->tap; t++) {
    double coeff_tmp = sin(2.0*M_PI*fbt->tune*t+fbt->phase*M_PI/180.0);
    norm += coeff_tmp;
    fbt->coeffs_FIR[t] = coeff_tmp; 
  }
  norm = norm/fbt->tap;

  if (fbt->tap>1) {
    for (t=0; t<fbt->tap; t++) {
      fbt->coeffs_FIR[t] = fbt->coeffs_FIR[t]-norm;
    }
  }
}

void
fbt_calc_coeffs_Spring8(fbt_feedback_t * fbt){
  int t;
  double depha = fbt->phase*M_PI/180.0;
  double * par1, * par2, * par3, * par4; 
  double norm_1 = 0;
  double norm_2 = 0;
  FIR_coeff(fbt->tap,fbt->tune,&par1,&par2,&par3,&par4);
  for (t=0; t<fbt->tap; t++) {
    double coeff_tmp =  sin(depha)*(par1[t]-par3[t])-cos(depha)*(par2[t]-par4[t])
      +sin(depha)*depha*par4[t]+cos(depha)*depha*par3[t];
    norm_1 += coeff_tmp*cos(depha*t);
    norm_2 += coeff_tmp*sin(depha*t);
    fbt->coeffs_FIR[t] = coeff_tmp;
  }
  for (t=0; t<fbt->tap; t++) {
    fbt->coeffs_FIR[t] = fbt->coeffs_FIR[t]/cabs(norm_1+I*norm_2);
  }
}

void 
fbt_calc_coeffs_energy_sensing(fbt_feedback_t * fbt) {
  int t;
  for (t=0; t<fbt->tap; t++) fbt->coeffs_FIR[t] = cos(fbt->phase*M_PI/180.0);
}

double
fbt_kick(fbt_feedback_t * fbt, int rev) {
  int t;
  double fbt_kick = 0;
  for (t=0;t<fbt->tap; t++) {
    //int index = rev-(fbt->tap-t-1)*fbt->downsampling;
    //int index = rev-t*fbt->downsampling;
    int index = (rev/fbt->downsampling)*fbt->downsampling-t*fbt->downsampling;
    fbt_kick += fbt->coeffs_FIR[t]*cyclic_array_get(index,fbt->offset_history)[0];
  }
  
  return fbt->gain*fbt_kick;
}

/**
 * @file
 * Module particle tunes
 */



#ifndef MBTRACK_TUNE_H
#define MBTRACK_TUNE_H

#include <stdbool.h>
#include <stdlib.h>
#include "types.h"
#include "bunch.h"
#include "transform_weak.h"


typedef struct incoherent_tune_model
{
  unsigned N;
  unsigned Nt;  /**< Length of tune data: Nt = N-th power of two */
  unsigned Ntp; /**< Number of bins / particles incoherent tune is determined at */
  unsigned kb;  /**< Bunch number reference partcles are in */
  unsigned Nstart; /**< Particle number to start from */  
  unsigned Nmax; /**< with zero-padding */
  
  double * datax; /**< Longitudinal pos. of Ntp reference particles over 2^Nt turns */
  double * dataxprime;
  
  double * tunes;  /**< Incoherent synchrotron tune of Ntp reference paricles */ 
  
  unsigned * pk_numbr; // [Bins]
  double * pk_value; // [Bins][Peaks-> max 10]
  double * pk_ampli; // [Bins][Peaks-> max 10]
}
tuneI_t;



int 
incoherent_tune_eval(tuneI_t * synctune, const long unsigned NrevTune,  FILE * fp);


int 
incoherent_tune_fill(tuneI_t * tune, const weak_bunch_t * bunch, unsigned long int rev);

typedef struct coherent_tune_model
{
  int TunePlane_C[3]; /**< Long. tune and transverse tunes */
  double UL;
  double UH;
  double UV;
  int Nt;  /**< Length of tune data: Nt-th power of two */

  double * xt_tune_C;
  double * xx_tune_C;
  double * xz_tune_C;
  
  double * freq_data;
  
  double * tunes;  /**< Incoherent synchrotron tune of Ntp reference paricles */
  
  unsigned * pk_numbr; // [Bins]
  double * pk_value; // [Bins][Peaks-> max 10]
  double * pk_ampli; // [Bins][Peaks-> max 10]
}
tuneC_t;


// bool
// incoherent_tune_init(tuneI_t * tune);
// 
// int 
// incoherent_tune_fill(tuneI_t * tune, const weak_bunch_t * bunch, unsigned int rev);
// 
// int 
// incoherent_tune_eval(tuneI_t * synctune, const unsigned NrevTune,  FILE * fp);
// 
// int 
// incoherent_tune_reset(const tuneI_t * tune, weak_bunch_t * bunch, const selffield_model_t * SelfFieldModel);
// 
// void
// tune_weak_incoherent_destroy(tuneI_t * tune);



bool
coherent_tune_init(tuneC_t * tune);

int 
identify_tunepeaks(tuneC_t * tune, double * cofu, const unsigned Pl, const long unsigned NrevTune,  FILE * fp);

int
acquire_tunedata_coherent(tuneC_t * tune, int kn, weak_bunch_t * bunch, plane_t plane, selffield_model_t * SelfFieldModel);


void
tune_weak_model_destroy(tuneC_t * tune);


#endif /* MBTRACK_TUNE_H */
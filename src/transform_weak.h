/**
 * @file
 * All transformations for weak bunch: Optical, Resonator selffield, RW longrange, Resonator longrange
 */

#ifndef TRANSFORM_WEAK_H
#define TRANSFORM_WEAK_H

#include <stdbool.h>
#include <stdlib.h>
#include "types.h"
#include "bunch.h"


/**
 * Optical transformation (for all planes) of a bunch (weak model)
 * 
 * @param bunch Bunch which is transformed
 * @param Pl Plane in which the transformation takes place
 */
int
transform_weak_bunch_optic(weak_bunch_t * bunch, 
                           const bunch_macroparticle_model_t * bunchModel, int iseed, int kb, const ring_t * ring);

/**
 * Initialize SelfField model, allocate memory for Greens function and set to zero
 * 
 * @param SelfFieldModel Parameters of the Selffield, to be initialized (result variable)
 */
bool
selffield_model_init(selffield_model_t * SelfFieldModel, double Nsigma, weak_bunch_t * bunch);

/**
 * Free dynamically allocated memory in SelfField model
 * 
 * @param SelfFieldModel Parameters of the Selffield, to be deleted
 */
void
selffield_model_destroy(selffield_model_t * SelfFieldModel);


/**
 * Calculates the Greens Function for all resonator impedance in all (tracked) planes
 * 
 * @param SelfFieldModel Binned selffield of ResImp resonator (result variable)
 */
int
construct_greensfunc_resonator(selffield_model_t * SelfFieldModel,
                     const ring_t * ring);


/**
 * Calculates the Greens Function contribution of the RW impedance in all (tracked) planes
 * 
 * @param SelfFieldModel Binned selffield including RW contribution (result variable)
 */
int
construct_greensfunc_RW(selffield_model_t * SelfFieldModel, 
                               const ring_t * ring);


/**
 * Transformation of electron distribution due to resonator + RW selffield
 * 
 * @param bunch Bunch which is transformed (result variable)
 * @param SelfFieldModel Binned selffield of all resonators and RW
 * @param Pl Plane in which the transformation takes place
 * @param rev Actual turn number
 * 
 */
int
transform_weak_bunch_selffield(weak_bunch_t * bunch,
			       const selffield_model_t SelfFieldModel, 
			       const cyclic_array_t * all_moments, const ring_t *ring,
			       long unsigned int rev, FILE * fp, double scan_val, int kb, 
			       e_beam_t * ebeam, double * phasor_end, double * phasor_end_HOR, double * phasor_end_VER,
			       double * fnp_ring, double * fnp_HOR, double * fnp_VER);


/**
 * Long-range resistive wall bunch transformation (weak model)
 */
int
transform_weak_bunch_RW_longrange(int in, unsigned int kb, const plane_t plane,
                         const e_beam_t ebeam, const ring_t * ring,
                         const bunch_macroparticle_model_t MPmodel,
                         const bunch_CM_history_strong_t * bCMhist,
                         weak_bunch_t * bunch);

int
transform_weak_bunch_RW_longrange_cyclic(const int in, const int bpos,
                                  cyclic_array_t * all_moments,
                                  const plane_t plane,
                                  const bunch_macroparticle_model_t * MPmodel,
                                  e_beam_t * ebeam,
                                  weak_bunch_t * bunch,
                                  ring_t * ring,
                                  tracking_t * track);

/**
 * Stores actual particle distribution in bins of all bunches
 * @param fnp_ring result variable
 */
int
fnp_ring_update(ring_t * ring, double * fnp_ring, const selffield_model_t * SelfFieldModel, weak_bunch_t * bunch, int kb, int plane);


/**
 * Initializes the phasor scheme and computes initial values
 * @param phasor_end result variable
 */
void
wake_phasor_init(ring_t * ring, double * fnp_ring, const selffield_model_t * SelfFieldModel, weak_bunch_t * bunch, int kb, double * phasor_end, e_beam_t * ebeam, int plane);


/**
 * Computes sum of phasors due to actual particle distribution
 * @param lr_wake result variable
 */
void
construct_wake_phasor(double * lr_wake, double * phasor_end, 
                       const selffield_model_t * SelfFieldModel,
		      const ring_t * ring, int kb, const double * fnp, e_beam_t * ebeam, weak_bunch_t * bunch, int plane);

void
fnp_ring_destroy(double * fnp_ring);


/**
 * Resistive wall contribution for first bin of LON greensfunction
 * @param wert result variable
 * @param uval averaged RW wake over [0; uval]
 */
int
RW_table_LON(double uval, double * wert);

int
RW_table_TRANS(double uval, double * wert);



#endif /* TRANSFORM_WEAK_H */

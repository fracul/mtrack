/**
 * @file
 * Manager and worker headers
 */


#include "types.h"
#include "bunch.h"
#include "fbi.h"

/**
 * Manager for weak-weak model
 */
void manager_weakweak(ring_t ring, const tracking_t track,
                      e_beam_t ebeam,
                      const bunch_macroparticle_model_t bunchModel, 
                      const selffield_model_t * SelfFieldModel);

/**
 * Worker for weak-weak model
 */
void worker_weakweak(ring_t ring, const tracking_t track,
                     e_beam_t ebeam,
                     const bunch_macroparticle_model_t bunchModel0,
                     selffield_model_t SelfFieldModel);
/**
 * Manager for weak-strong model
 */
void manager_weakstrong(const ring_t ring, const tracking_t track,
                        e_beam_t ebeam,
                        const bunch_macroparticle_model_t macrop_model);
/**
 * Worker for weak-strong model
 */
void worker_weakstrong(const ring_t ring, const tracking_t track,
                       e_beam_t ebeam,
                       const bunch_macroparticle_model_t macrop_model);


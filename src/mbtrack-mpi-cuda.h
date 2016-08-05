/**
 * @file
 * Manager and worker headers that use cuda
 * Uses two mpi nodes
 * worker node runs simulation on all the bunches
 * manager node saves output
 * if cuda is available worker offloads simulations to GPU
 */

#include "types.h"
#include "bunch.h"
#include "fbi.h"

/**
 * Manager for weak-weak model that uses cuda
 */
void manager_weakweak_cuda(ring_t ring, const tracking_t track,
			   e_beam_t ebeam,
			   const bunch_macroparticle_model_t bunchModel, 
			   selffield_model_t SelfFieldModel);

/**
 * Worker for weak-weak model that uses cuda
 */
void worker_weakweak_cuda(ring_t ring, const tracking_t track,
			  e_beam_t ebeam,
			  const bunch_macroparticle_model_t bunchModel0,
			  selffield_model_t SelfFieldModel);

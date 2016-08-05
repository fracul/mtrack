/**
 * @file
 * transform_weak_bunch_optic running on GPU using CUDA
 */

#ifndef TRANSFORM_WEAK_CUDA_H
#define TRANSFORM_WEAK_CUDA_H

#include "types.h"
#include "bunch.h"


/**
 * Optical transformation (for all planes) of a bunch (weak model)
 * 
 * @param bunch Bunch which is transformed
 * @param Pl Plane in which the transformation takes place
 */
int
transform_weak_bunch_optic_cuda(weak_bunch_t * bunch, const ring_t * ring,
				const bunch_macroparticle_model_t * bunchModel, int iseed, int kb);

void 
transfer_ring(const ring_t * ring, const bunch_macroparticle_model_t * bunchModel, const tracking_t * track, const selffield_model_t * SelfFieldModel, int Np, int seed);

void 
transfer_bunch(weak_bunch_t * bunch, int kb);

void
setup_cuda(int nbunches, int np, int ncell);

void 
allocate_bunch(weak_bunch_t * bunch, int kb);

void 
free_cuda();

void 
free_bunch(weak_bunch_t * bunch, int kb);

void 
transfer_bunch_to_device(weak_bunch_t * bunch, int kb);

void
transfer_bunch_from_device(weak_bunch_t * bunch, int kb);

int transform_weak_bunch_selffield_cuda(weak_bunch_t * bunch, const selffield_model_t SelfFieldModel,
					int kb, FILE *fp, int rev, double scan_val);


#endif

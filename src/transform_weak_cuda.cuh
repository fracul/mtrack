/**
 * @file
 * transform_weak_bunch_optic running on GPU using CUDA
 */

#ifndef TRANSFORM_WEAK_CUDA_H
#define TRANSFORM_WEAK_CUDA_H



#include "types.h"
#include "bunch.h"


//#include <thrust/device_vector.h>
//#include <thrust/sort.h>
//#include <thrust/count.h>

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
transfer_ebeam(const e_beam_t *ebeam);

void 
transfer_phasor(int nbunches, int resonators, double *phasor_end);

void
setup_cuda(int nbunches, int np, int ncell, double *fnp_ring);

void 
allocate_bunch(weak_bunch_t * bunch, int kb);

void 
sync_device();

void 
free_cuda();

void 
free_bunch(weak_bunch_t * bunch, int kb);

void 
transfer_bunch_to_device(weak_bunch_t * bunch, int kb);

void
transfer_bunch_from_device(weak_bunch_t * bunch, int kb);

void 
fnp_ring_update_cuda(const weak_bunch_t * bunch, const selffield_model_t SelfFieldModel,
		     int kb);

int
construct_wake_phasor_cuda(int Nbunch, int Np, int Ncell, int kb, int resonators);

int 
transform_weak_bunch_selffield_cuda(weak_bunch_t * bunch, const selffield_model_t SelfFieldModel,
				    int kb, FILE *fp, int rev, double scan_val, int resonators);

void
transform_bunch_RW_longrange_cyclic_cuda(const int in, const int bpos, const int plane,
					   const bunch_macroparticle_model_t * MPmodel,
					   ring_t *ring, tracking_t * track, 
					   int kb, int Np, cyclic_array_t *dipole_RW);
  
void
initialize_cyclic_array_cuda(cyclic_array_t dipole_RW);

void
update_cyclic_array_cuda(cyclic_array_t *dipole_RW);

void 
free_cyclic_array_cuda(cyclic_array_t *dipole_RW);


#endif

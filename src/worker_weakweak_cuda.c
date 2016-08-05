#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include "confmpi.h"
#include "input.h"
#include "types.h"
#include "bunch.h"
#include "fbi.h"
#include "statistics.h"
#include "transform_weak.h"
#include "tracking.h"
#include "cyclic_array.h"
#include "tune.h"
#include <time.h>
#include <sys/time.h>

#ifdef MBTRACK_CUDA
#include "transform_weak_cuda.cuh"
#endif

void worker_weakweak_cuda(ring_t ring, const tracking_t track,
			  e_beam_t ebeam,
			  const bunch_macroparticle_model_t bunchModel0,
			  selffield_model_t SelfFieldModel)
{

//TODO: receive pointers to device memory from manager and take over data transfer and output
//for cuda this can be done usinc IPC to sedn device memory address
//for host? do we even need host version of two processes?

printf("Worker process doing nothing!\n");

}

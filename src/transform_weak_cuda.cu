extern "C" {
#include "transform_weak_cuda.cuh"
}

#include <stdio.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>

#define BLOCK_SIZE 128

/*
 *TODO: handling of constant memory: dont transfer all of the ring just the variables needed? 
 * this applies to all the structs holding conf parameters (ring, bunchModel, tracking etc.).
 * some pros and const to consider:
 * + better distinction between __constant__ memory that can be used for static parameters
 *   and global memroy that will need to be used for dynamic arrays in the struct
 * + not all the variables are needed on the device
 * - copy each variable seperatelly using cudaMemcpyToSymbol will result in lots of copies
 *   may be slower and will make the code more cluttered
 */
__constant__ ring_t dring;
__constant__ bunch_macroparticle_model_t dbunchModel;
__constant__ tracking_t dtrack;
__constant__ selffield_model_t dSelfFieldModel;

__device__ active_HC_t *dactive_HC;

__device__ int dcellmin;
__device__ int dcellmax;

//array to hold particle in a cell
int *dmapcell;

//device arrays to hold greens functions
double *dSelfFieldGl;
double *dSelfFieldGh;
double *dSelfFieldGv;

//temporary arrays to hold values for potential calcualtions
//array size is ncell * nbunch, so each bunch gets a seperate storage
double *dfnp;
double *ddipoleV;
double *ddipoleH;
double *dGL1;
double *dGlambdaV;
double *dGlambdaH;
double *dlr_wake;

curandState *d_cudaRndStates;
particle_t **d_particles;
active_HC_t *hactive_HC;

__global__ void kernelInitRandoms(curandState *state, int size, int seed) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < size)
    curand_init(seed, idx, 0, &state[idx]);
}

__global__ void kernelTransformWeakBunchOptic(particle_t *particles, curandState *rndState, 
					      int np, int kb) 
{

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int tid = threadIdx.x;

  //load active_HC in shared memory
  extern __shared__ active_HC_t sactive_HC[];
  while (tid < dring.active_HC_size) {
    sactive_HC[tid] = dactive_HC[tid];
    tid += blockDim.x;
  }

  if (idx < np) {

    //TODO: Do this once and store in constant memory?
    /* Longitudinal trafo, ALWAYS */
    double rffocus;
    double dringT0 = dring.T0;
    const double damp = 1.0 - dring.De;
    const double excite = 2.0 * dbunchModel.slope_sgm.v[LON] * sqrt(dringT0/dring.taue);
    const double quantumV_pos = dbunchModel.pos_sgm.v[VER] * sqrt(2.0 * dringT0 / dring.tauz);
    const double quantumV_slope = dbunchModel.slope_sgm.v[VER] * sqrt(2.0 * dringT0 / dring.tauz);
    const double quantumH_pos = dbunchModel.pos_sgm.v[HOR] * sqrt(2.0 * dringT0 / dring.taux);
    const double quantumH_slope = dbunchModel.slope_sgm.v[HOR] * sqrt(2.0 * dringT0 / dring.taux);
    
    double eps_const = dring.Vrf0 / (dring.E0 * FKILO);

    curandState state = rndState[idx];
    particle_t particle0 = particles[idx];
    particle_t particle = particle0;

    double xeps_gainj;

    if (dtrack.EnableIdealHC == 0) {
      xeps_gainj = eps_const * sin(dring.wrf * particle0.pos.xtau + dring.phai0);
    
      //adding potential of active HCs
      for (int j = 0; j < dring.active_HC_size; j++) {
	active_HC_t aHC = sactive_HC[j];

	xeps_gainj += aHC.Vpeak/(dring.E0 * FGIGA) * sin(dring.wrf * aHC.nHC * particle0.pos.xtau + aHC.phi_aHC + aHC.nHC * kb * 2 * M_PI);
      }
    } else {
      xeps_gainj = eps_const * (sin(dring.wrf * particle0.pos.xtau + dring.phai0) + dring.HC_k * sin(dring.m_aHC * dring.wrf * particle0.pos.xtau + dring.m_aHC * dring.phi_n));
    }

    particle.slope.xtau = particle0.slope.xtau + xeps_gainj - dring.Urad;

    if (dtrack.EnableQuantum) {
      particle.slope.xtau = particle.slope.xtau * damp + excite * curand_normal_double(&state);
      rffocus = (1.0 + particle0.slope.xtau)/ (1.0 + particle0.slope.xtau + xeps_gainj);
    }
    particle.pos.xtau = particle0.pos.xtau - particle.slope.xtau * dringT0*dring.ac;


    //********************************************************//
    if (dtrack.TrackPlane[VER]) {
      double PsiVj, cosVj, sinVj, amv11j, amv21j, amv12j, amv22j;
      double PsiV0   = 2.0*M_PI*dring.QV0;

      PsiVj  = PsiV0 * (1.0 + dring.Gziz * particle0.slope.xtau);
      cosVj  = cos(PsiVj);
      sinVj  = sin(PsiVj);
      amv11j = cosVj + dring.alpha1[VER]*sinVj,   amv12j =         dring.beta1[VER]*sinVj;
      amv21j =       - dring.gamma1[VER]*sinVj,   amv22j = cosVj - dring.alpha1[VER]*sinVj;

      particle.pos.z = amv11j*particle0.pos.z + amv12j*particle0.slope.z;
      particle.slope.z = amv21j*particle0.pos.z + amv22j*particle0.slope.z;

      if(dtrack.EnableQuantum) {
	particle.pos.z = particle.pos.z + quantumV_pos * curand_normal_double(&state);
	particle.slope.z = rffocus * particle.slope.z + quantumV_slope * curand_normal_double(&state); 
      }

    }

    if (dtrack.TrackPlane[HOR]) {
      double PsiHj, cosHj, sinHj, amh11j, amh12j, amh13j, amh21j, amh22j, amh23j;
      double PsiH0   = 2.0*M_PI*dring.QH0; 

      PsiHj  = PsiH0 * (1.0 + dring.Gzix * particle0.slope.xtau);
      cosHj  = cos(PsiHj);  
      sinHj = sin(PsiHj);          
      amh11j = cosHj + dring.alpha1[HOR]*sinHj;
      amh12j =         dring.beta1[HOR]*sinHj;
      amh13j = (1.0 - amh11j)*dring.dispH1 - amh12j*dring.disppH1;
      amh21j =       - dring.gamma1[HOR]*sinHj;
      amh22j = cosHj - dring.alpha1[HOR]*sinHj;
      amh23j =       -amh21j *dring.dispH1 + (1.0 - amh22j)*dring.disppH1;
      
      particle.pos.x = amh11j*particle0.pos.x + amh12j*particle0.slope.x  + amh13j*particle.slope.xtau;
      particle.slope.x = amh21j*particle0.pos.x + amh22j*particle0.slope.x + amh23j*particle.slope.xtau;
          
      if(dtrack.EnableQuantum) {
	particle.pos.x = particle.pos.x + quantumH_pos * curand_normal_double(&state);
	particle.slope.x = rffocus * particle.slope.x + quantumH_slope * curand_normal_double(&state); 
      } 
    }

    rndState[idx] = state;    
    particles[idx] = particle;
  }

}

__global__ void kernelAssignBin(particle_t *particles, int *mapcell, int np, int ncell)
{

  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  if (idx < np) {
    //load particle from global memory
    particle_t particle = particles[idx];

    double sgm_xtau = dSelfFieldModel.sigma_tau; /* [s] */
    double Tau = -particle.pos.xtau / sgm_xtau;

    int icell = (int) ((Tau + dSelfFieldModel.Nsigma)/dSelfFieldModel.dT + 1);
    //TODO: N_trash_low and N_trash_high should added here, will require 2 more attomicAdds
    if (icell < 0) icell = 0;
    if (icell >= ncell - 1) icell = ncell - 1;

    mapcell[idx] = icell;
  }
}

__device__ double atomicAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
                                          (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed, 
                        __double_as_longlong(val + 
                        __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}

__global__ void kernelCountBinAtomic(particle_t *particles, double *dipoleV, double *dipoleH,
				     double *fnp, int *mapcell, int np, int ncell, int nblock)
{

  int idx = blockIdx.x * blockDim.x + threadIdx.x;

  //create arrays in shared memory
  extern __shared__ double smem[];
  double *s_fnp = (double*)smem;
  double *s_dipoleV = (double*)&smem[ncell];
  double *s_dipoleH = (double*)&smem[2*ncell];

  //zero shared memory
  for (int id = threadIdx.x; id < ncell; id += blockDim.x) {
    s_fnp[id] = 0.0;
    s_dipoleV[id] = 0.0;
    s_dipoleH[id] = 0.0;
  }

  __syncthreads();

  //calculate local histogram 
  if (idx < np) {
    int map = mapcell[idx];
    particle_t p = particles[idx];
    
    atomicAdd(&s_fnp[map], 1.0);
    atomicAdd(&s_dipoleV[map], p.pos.z);
    atomicAdd(&s_dipoleH[map], p.pos.x);
  }

  __syncthreads();

  //load local hist to global memory
  for (int id = threadIdx.x; id < ncell; id += blockDim.x) {
    atomicAdd(&fnp[id], s_fnp[id]);
    atomicAdd(&dipoleV[id], s_dipoleV[id]);
    atomicAdd(&dipoleH[id], s_dipoleH[id]);
  }

}


//one block per bin, each block 
__global__ void kernelCountBin(particle_t *particles, double *dipoleV, double *dipoleH, double *fnp, 
			       int *mapcell, int np, int ncell)
{
  
  int bin = blockIdx.x;
  int blockSize = blockDim.x;
  int tid = threadIdx.x;


  //create arrays in shared memory
  extern __shared__ double smem[];
  double *s_fnp = (double*)smem;
  double *s_dipoleV = (double*)&smem[blockSize];
  double *s_dipoleH = (double*)&smem[2*blockSize];

  //set shared memory to 0.0
  for (int id = threadIdx.x; id < blockSize; id += blockSize) {
    s_fnp[id] = 0.0;
    s_dipoleV[id] = 0.0;
    s_dipoleH[id] = 0.0;
  }

  //make sure global and shared memory is set in all the threads
  __syncthreads();

  //each block loops trough np with a step of blockDim.x
  for (int id = threadIdx.x; id < np; id += blockSize) {
    //load map
    int map = mapcell[id];

    //if bin of the particle == blockIdx.x load particle from global memory
    //and calculate local fnp, dipoleV, dipoleH in shared memory
    if (bin == map) {
      particle_t p = particles[id];
      s_fnp[tid] += 1.0;
      s_dipoleV[tid] += p.pos.z;
      s_dipoleH[tid] += p.pos.x;
    }
  }

  __syncthreads();

  //when particle loop is done perform parallel reduction to get total value per bin
  int id = threadIdx.x;
  for (int step = blockSize >> 1; step > 0; step = step >> 1) {
    if (id < step) {
      s_fnp[id] += s_fnp[id + step];
      s_dipoleV[id] += s_dipoleV[id + step];
      s_dipoleH[id] += s_dipoleH[id + step];
    }
  }

  __syncthreads();

  //write out to global memory
  if (tid == 0) {
    fnp[bin] = s_fnp[0];
    dipoleV[bin] = s_dipoleV[0];
    dipoleH[bin] = s_dipoleH[0];
  }
  

}
//thrust::find to search for icellmin and icellmax

//launch one block with multiple threads, load fnp in shared memory and then thread 0 loops trough it
//will be slow, but ncell is small, so it will be small % of the whole
//simulation, not really worth using thrust
__global__ void kernelGetMinMax(double *fnp) {

  extern __shared__ double s_fnp[];

  int tid = threadIdx.x;
  int Ncell = dSelfFieldModel.Ncell;

  for (int id = threadIdx.x; id < Ncell; id += blockDim.x)
    s_fnp[id] = fnp[id];

  if (tid == 0) {
    int icellmin = Ncell - 1;
    int icellmax = 0;

    for (int icell = 0; icell < Ncell - 1; icell++) {
      if (s_fnp[icell] > 0) {
	icellmin = icell;
	break;
      }
    }
    for (int icell = Ncell - 1; icell > 0; icell--) {
      if (s_fnp[icell] > 0) {
	icellmax = icell;
	break;
      }
    }
    dcellmin = icellmin;
    dcellmax = icellmax;
  }
}

//from icellmin + 1, to icellmax + 1, before that zero GL1
__global__ void kernelWakePotential(int ncell, int icellmin, int icellmax, double *GL1, double *fnp,
				    double *SelfFieldGl)
{
  int icell = blockIdx.x * blockDim.x + threadIdx.x + icellmin + 1;

  //load fnp and SelfField.Gl to shared memory
  extern __shared__ double smem[];
  double * s_fnp = (double*)smem;
  double * s_Gl = (double*)&smem[ncell];

  for(int tid = threadIdx.x; tid < ncell; tid += blockDim.x) 
  {
    s_fnp[tid] = fnp[tid];
    s_Gl[tid] = SelfFieldGl[tid];
  }

  //make sure all threads complete shared memory load before continue
  __syncthreads();

  if (icell < icellmax + 1) {
    double l_GL1 = 0.0;

    /***  Contribution of BBR impedances  ***/
    double Gn = 0.0;
    for (int k = icell - 1; k >= icellmin; k--) {
      Gn += s_fnp[k] * s_Gl[icell - k];
    }
    l_GL1 = Gn + s_fnp[icell] * s_Gl[0];

    /***  Contribution of purely resistive impedances  ***/
    if (dSelfFieldModel.RsisL > 0)
    {
      double sgm_xtau = dSelfFieldModel.sigma_tau; /* [s] */
      double ntmoy = (0.5 * s_fnp[icell-1] + s_fnp[icell] + 0.5 * s_fnp[icell+1]) / 2.0;
      l_GL1 += dSelfFieldModel.RsisL * ntmoy/(sgm_xtau * dSelfFieldModel.dT);
    }

    /***  Contribution of purely inductive impedances  ***/
    if (dSelfFieldModel.aindL > 0)
    {
      double sgm_xtau = dSelfFieldModel.sigma_tau; /* [s] */
      double sgmatau2 = sgm_xtau * sgm_xtau; //for purely inductive
      double dTau2    = dSelfFieldModel.dT * dSelfFieldModel.dT; //for purely inductive
      double ntmoy0 = (0.25 * s_fnp[icell-3] + 0.5 * s_fnp[icell-2] + s_fnp[icell-1] + 0.5 * s_fnp[icell]  + 0.25 * s_fnp[icell+1]) / 2.5;
      double ntmoy1 = (0.25 * s_fnp[icell-1] + 0.5 * s_fnp[icell]   + s_fnp[icell+1] + 0.5 * s_fnp[icell+2] + 0.25 * s_fnp[icell+3]) / 2.5;
      l_GL1 += (ntmoy1 - ntmoy0) / (2.0 * sgmatau2 * dTau2) * dSelfFieldModel.aindL;
    }

    //write l_GL1 to global memory
    GL1[icell] = l_GL1;
  }

}

__global__ void kernelWakePotentialPlanesVH(double *Glambda, double *Gdipole, double *SelfFieldGl, 
					   int ncell) 
{

  int icell = blockIdx.x * blockDim.x + threadIdx.x + 1;

  //load fnp and SelfField.Gl to shared memory
  extern __shared__ double smem[];
  double * s_dipole = (double*)smem;
  double * s_Gl = (double*)&smem[ncell];
  
  for(int tid = threadIdx.x; tid < ncell; tid += blockDim.x) 
  {
    s_dipole[tid] = Gdipole[tid];
    s_Gl[tid] = SelfFieldGl[tid];
  }

  //make sure all threads complete shared memory load before continue
  __syncthreads();

  if (icell < ncell) {
    double Gn = 0.0;
    for (int k = icell - 1; k >= 0; k--)
      Gn += s_dipole[k] * s_Gl[icell - k];

    Glambda[icell] = Gn;
  }

}

__global__ void kernelWakePotentialEffect(particle_t *particles, int *mapcell, double *GL1, double *GlambdaV, double *GlambdaH, double bunchIb, int Np, int Ncell) {

  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  int tid = threadIdx.x;
  
  //load GL1 to shared memory
  extern __shared__ double smem[];
  double *sGL1 = (double*)smem;
  double *sGlambdaV = (double*)&smem[Ncell];
  double *sGlambdaH = (double*)&smem[2*Ncell];
  while (tid < Ncell) {
    sGL1[tid] = GL1[tid];
    sGlambdaV[tid] = GlambdaV[tid];
    sGlambdaH[tid] = GlambdaH[tid];
    tid += blockDim.x;
  }

  __syncthreads();

  if (idx < Np) {
    int map = mapcell[idx];

    double factG = -dring.T0 * bunchIb / (Np * dring.E0 * FGIGA);
    particles[idx].slope.xtau += factG * sGL1[map];

    if (dSelfFieldModel.PlaneV > 0)
      particles[idx].slope.z += -factG * sGlambdaV[map];

    if (dSelfFieldModel.PlaneH > 0)
      particles[idx].slope.x += -factG * sGlambdaH[map];
  }

}


 void transfer_ring(const ring_t * ring, const bunch_macroparticle_model_t * bunchModel, 
		    const tracking_t * track, const selffield_model_t * SelfFieldModel, 
		    int Np, int seed)
{
  
  cudaError_t e1, e2, e3;

  //setup random numbers
  int threads = BLOCK_SIZE;
  int blocks = Np / threads + 1;
  kernelInitRandoms<<<blocks, threads>>>(d_cudaRndStates, Np, seed);
  e1 = cudaGetLastError();
  if (e1 != cudaSuccess)
    fprintf(stderr, "Error ! CUDA init randoms failed\n");

  //copy bunch model to device constant memory
  e1 = cudaMemcpyToSymbol( dbunchModel, bunchModel, sizeof(bunch_macroparticle_model_t) );
  if (e1 != cudaSuccess)
    fprintf(stderr, "Error ! CUDA bunchModel to device failed\n");

  //copy ring to device constant memory
  e1 = cudaMemcpyToSymbol( dring, ring, sizeof(ring_t) );
  if (e1 != cudaSuccess)
    fprintf(stderr, "Error ! CUDA ring to device failed\n");

  //allocate memory for dactive_HC in global memory
  e1 = cudaMalloc((void**) &hactive_HC, sizeof(active_HC_t) * ring->active_HC_size);
  e2 = cudaMemcpy(hactive_HC, ring->active_HC, ring->active_HC_size * sizeof(active_HC_t), 
	     cudaMemcpyHostToDevice);
  e3 = cudaMemcpyToSymbol(dactive_HC, &hactive_HC, sizeof(hactive_HC));
  if (e1 != cudaSuccess)
    fprintf(stderr, "Error ! CUDA ring to device malloc failed\n");
  if (e2 != cudaSuccess)
    fprintf(stderr, "Error ! CUDA ring to device memcpy failed\n");
  if (e3 != cudaSuccess)
    fprintf(stderr, "Error ! CUDA ring to device memcpytosymbol failed\n");

  //TODO: copy ring_resonators, longrange_resonators, lr_wake, active_HC to device memory

  //copy track to device constant memory
  e1 = cudaMemcpyToSymbol( dtrack, track, sizeof(tracking_t) );
  if (e1 != cudaSuccess)
    fprintf(stderr, "Error ! CUDA track to device failed\n");

  //copy SelfFieldModel to device memroy
  e1 = cudaMemcpyToSymbol( dSelfFieldModel, SelfFieldModel, sizeof(selffield_model) );
  if (e1 != cudaSuccess)
    fprintf(stderr, "Error ! CUDA SelfFieldModel to device failed\n");

  //copy greens functions to GPU
  e1 = cudaMemcpy( dSelfFieldGl, SelfFieldModel->Gl, sizeof(double) * SelfFieldModel->Ncell,
		   cudaMemcpyHostToDevice);
  e2 = cudaMemcpy( dSelfFieldGh, SelfFieldModel->Gh, sizeof(double) * SelfFieldModel->Ncell,
		   cudaMemcpyHostToDevice);
  e3 = cudaMemcpy( dSelfFieldGv, SelfFieldModel->Gv, sizeof(double) * SelfFieldModel->Ncell,
		   cudaMemcpyHostToDevice);
  if (e1 != cudaSuccess || e2 != cudaSuccess || e3 != cudaSuccess)
    fprintf(stderr, "Error ! CUDA copy Greens functions failed!\n");

}

void setup_cuda(int nbunches, int np, int ncell) {
  int ndevices = 0;
  cudaGetDeviceCount(&ndevices);
  cudaSetDevice(ndevices - 1);

  //create device pointer for each bunch
  d_particles = new particle_t*[nbunches];

  //allocate device memory for random numbers
  cudaError_t e1 = cudaMalloc( (void**) &d_cudaRndStates, sizeof(curandState) * np );

  //allocate memory for Greens functions
  e1 = cudaMalloc( (void**) &dSelfFieldGl, sizeof(double) * ncell );
  e1 = cudaMalloc( (void**) &dSelfFieldGv, sizeof(double) * ncell );
  e1 = cudaMalloc( (void**) &dSelfFieldGh, sizeof(double) * ncell );

  //allocate memory for each bunch for temporary arrays
  cudaMalloc( (void**) &dfnp, sizeof(double) * ncell * nbunches);
  cudaMalloc( (void**) &ddipoleV, sizeof(double) * ncell * nbunches);
  cudaMalloc( (void**) &ddipoleH, sizeof(double) * ncell * nbunches);
  cudaMalloc( (void**) &dGL1, sizeof(double) * ncell * nbunches);
  cudaMalloc( (void**) &dGlambdaV, sizeof(double) * ncell * nbunches);
  cudaMalloc( (void**) &dGlambdaH, sizeof(double) * ncell * nbunches);
  cudaMalloc( (void**) &dlr_wake, sizeof(double) * ncell * nbunches);

  
  cudaMemset(dfnp, 0, sizeof(double) * ncell * nbunches);
  cudaMemset(ddipoleV, 0, sizeof(double) * ncell * nbunches);
  cudaMemset(ddipoleH, 0, sizeof(double) * ncell * nbunches);
  cudaMemset(dGL1, 0, sizeof(double) * ncell * nbunches);
  cudaMemset(dGlambdaV, 0, sizeof(double) * ncell * nbunches);
  cudaMemset(dGlambdaH, 0, sizeof(double) * ncell * nbunches);
  cudaMemset(dlr_wake, 0, sizeof(double) * ncell * nbunches);
  
  cudaMalloc( (void**) &dmapcell, sizeof(int) * np);

}

void allocate_bunch(weak_bunch_t * bunch, int kb) {
  
  //allocate memory on the GPU for particle array
  cudaError_t e1 = cudaMalloc( (void**) &d_particles[kb], sizeof(particle_t) * bunch->Np);

  if (e1 != cudaSuccess)
    fprintf(stderr, "Error ! CUDA memory allocation failed!\n");

  /*
   * TODO: for 480 bunches might be a bad idea to page lock memory, but then we loose data 
   * transfer speed speed and true async possibilities. This needs some testing how much
   * host memory we can lock before we run into trouble. (If we limit the output to files and
   * move all the simulation to GPU this might not be so important)
   */
  //page lock bunch->particles for faster host->device->host transfers
  e1 = cudaHostRegister(bunch->particles, sizeof(particle_t) * bunch->Np, cudaHostRegisterPortable);
  if (e1 != cudaSuccess)
    fprintf(stderr, "Error ! CUDA host register memory failed!\n");

}

void free_cuda() {
  //free device memory
  cudaFree(d_cudaRndStates);
  cudaFree(hactive_HC);
  cudaFree(dSelfFieldGl);
  cudaFree(dSelfFieldGh);
  cudaFree(dSelfFieldGv);

  cudaFree(dfnp);
  cudaFree(ddipoleV);
  cudaFree(ddipoleH);
  cudaFree(dGL1);
  cudaFree(dGlambdaV);
  cudaFree(dGlambdaH);
  cudaFree(dlr_wake);

  cudaFree(dmapcell);
}

void free_bunch(weak_bunch_t * bunch, int kb) {
  //free device memory
  cudaFree(d_particles[kb]);

  //unregister page-locked memory from
  cudaHostUnregister(bunch->particles);
}

void transfer_bunch_to_device(weak_bunch_t * bunch, int kb) {

  cudaError_t e;
  e = cudaMemcpy(d_particles[kb], bunch->particles, bunch->Np * sizeof(particle_t), 
		 cudaMemcpyHostToDevice); 
  if (e != cudaSuccess)
    fprintf(stderr, "Error ! CUDA bunch to device failed\n");
}

void transfer_bunch_from_device(weak_bunch_t * bunch, int kb) {
  cudaError_t e;
  e = cudaMemcpy(bunch->particles, d_particles[kb], bunch->Np * sizeof(particle_t), 
		 cudaMemcpyDeviceToHost); 
  if (e != cudaSuccess)
    fprintf(stderr, "Error ! CUDA bunch from device failed\n");
}

int
transform_weak_bunch_optic_cuda(weak_bunch_t * bunch, const ring_t * ring, 
				const bunch_macroparticle_model_t * bunchModel, int iseed, int kb)
{

  int threads = BLOCK_SIZE;
  int blocks = bunch->Np / threads + 1;

  int shared_size = ring->active_HC_size * sizeof(active_HC_t);

  kernelTransformWeakBunchOptic<<<blocks, threads, shared_size>>>(d_particles[kb], 
								  d_cudaRndStates, 
								  bunch->Np,
								  kb);

  cudaError_t e1 = cudaGetLastError();
  if (e1 != cudaSuccess) {
    fprintf(stderr, "Error ! CUDA exec transform weak bunch optic failed\n");
    return -1;
  }
  

  return 1;
}

int
transform_weak_bunch_selffield_cuda(weak_bunch_t * bunch, const selffield_model_t SelfFieldModel,
				    int kb, FILE *fp, int rev, double scan_val) 
{

  int offset = kb * SelfFieldModel.Ncell;
  int bytes = SelfFieldModel.Ncell * sizeof(double);

  //reset temporary GPU memory
  cudaMemset(&dfnp[offset], 0, bytes);
  cudaMemset(&ddipoleV[offset], 0, bytes);
  cudaMemset(&ddipoleH[offset], 0, bytes);
  cudaMemset(&dGL1[offset], 0, bytes);
  cudaMemset(&dGlambdaV[offset], 0, bytes);
  cudaMemset(&dGlambdaH[offset], 0, bytes);
  cudaMemset(&dlr_wake[offset], 0, bytes);

  //assign each particle a bin
  int smem_size = 0;
  int threads = 128;
  int blocks = bunch->Np / threads + 1;
  kernelAssignBin<<<blocks, threads>>>(d_particles[kb], dmapcell, bunch->Np, SelfFieldModel.Ncell);

  //count particles per bin and sum H and V positions inside a bin
  threads = 1024;
  blocks = bunch->Np / threads + 1;
  smem_size = 3 * SelfFieldModel.Ncell * sizeof(double);
  cudaMemset(&dfnp[kb * SelfFieldModel.Ncell], 0, SelfFieldModel.Ncell * sizeof(double));
  cudaMemset(&ddipoleV[kb * SelfFieldModel.Ncell], 0, SelfFieldModel.Ncell * sizeof(double));
  cudaMemset(&ddipoleH[kb * SelfFieldModel.Ncell], 0, SelfFieldModel.Ncell * sizeof(double));
  kernelCountBinAtomic<<<blocks, threads, smem_size >>>(d_particles[kb], 
							&ddipoleV[kb * SelfFieldModel.Ncell], 
							&ddipoleH[kb * SelfFieldModel.Ncell], 
							&dfnp[kb * SelfFieldModel.Ncell], 
							dmapcell, bunch->Np, SelfFieldModel.Ncell,
							bunch->Np / blocks);
  

  //find min and max bins with particles
  threads = 128;
  blocks = 1;
  //smem_size = threads * sizeof(double);
  smem_size = SelfFieldModel.Ncell * sizeof(double);
  kernelGetMinMax<<<blocks, threads, smem_size>>>(&dfnp[kb * SelfFieldModel.Ncell]);
 
  int imin, imax;
  cudaMemcpyFromSymbol(&imin, dcellmin, sizeof(int)); 
  cudaMemcpyFromSymbol(&imax, dcellmax, sizeof(int));

  //get the effects of wake potentials on LON plane
  //TODO: use streams to increase device utilization by launching all three kernels
  //simultaniously
  threads = 128;
  blocks = (imax - imin) / threads + 1;
  smem_size = 2 * SelfFieldModel.Ncell * sizeof(double);
  kernelWakePotential<<<blocks, threads, smem_size>>>(SelfFieldModel.Ncell, imin, imax, 
						      &dGL1[kb * SelfFieldModel.Ncell],
						      &dfnp[kb * SelfFieldModel.Ncell], 
						      dSelfFieldGl);

  //calc the effects of wake potentials on VER and HOR planes
  threads = 128;
  blocks = (SelfFieldModel.Ncell - 1) / threads + 1;
  smem_size = 2 * SelfFieldModel.Ncell * sizeof(double);
  if (SelfFieldModel.PlaneV > 0)
    kernelWakePotentialPlanesVH<<<blocks, threads, smem_size>>>(&dGlambdaV[kb * SelfFieldModel.Ncell], 
								&ddipoleV[kb * SelfFieldModel.Ncell],
								dSelfFieldGl,
								SelfFieldModel.Ncell);
  if (SelfFieldModel.PlaneH > 0)
    kernelWakePotentialPlanesVH<<<blocks, threads, smem_size>>>(&dGlambdaH[kb * SelfFieldModel.Ncell], 
								&ddipoleH[kb * SelfFieldModel.Ncell],
								dSelfFieldGl,
								SelfFieldModel.Ncell);

  //apply effects of wake potentials to particles
  threads = 128;
  blocks = bunch->Np / threads + 1;
  smem_size = 3 * SelfFieldModel.Ncell * sizeof(double);
  kernelWakePotentialEffect<<<blocks, threads, smem_size>>>(d_particles[kb], dmapcell, 
							    &dGL1[kb * SelfFieldModel.Ncell],
							    &dGlambdaV[kb * SelfFieldModel.Ncell],
							    &dGlambdaH[kb * SelfFieldModel.Ncell],
							    bunch->Ib, bunch->Np,
							    SelfFieldModel.Ncell);

  cudaError_t e1 = cudaGetLastError();
  if (e1 != cudaSuccess)
    fprintf(stderr, "Error ! CUDA error transform weak bunch selffield!\n");


  //log potentials
  if (bunch->kb_out == 1 && rev%track.NrevMon == 0) {
    double fnp[SelfFieldModel.Ncell], GL1[SelfFieldModel.Ncell];
    double dipoleV[SelfFieldModel.Ncell], GlambdaV[SelfFieldModel.Ncell];
    double dipoleH[SelfFieldModel.Ncell], GlambdaH[SelfFieldModel.Ncell];
    double lr_wake[SelfFieldModel.Ncell];

    memset(fnp, 0.0, sizeof(fnp));
    memset(dipoleV, 0.0, sizeof(dipoleV));
    memset(dipoleH, 0.0, sizeof(dipoleH));
    memset(GL1, 0.0, sizeof(GL1));
    memset(GlambdaV, 0.0, sizeof(GlambdaV));
    memset(GlambdaH, 0.0, sizeof(GlambdaH));
    memset(lr_wake, 0.0, sizeof(lr_wake));   

    cudaMemcpy(fnp, &dfnp[offset], bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(GL1, &dGL1[offset], bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(dipoleV, &ddipoleV[offset], bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(GlambdaV, &dGlambdaV[offset], bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(dipoleH, &ddipoleH[offset], bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(GlambdaH, &dGlambdaH[offset], bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(lr_wake, &dlr_wake[offset], bytes, cudaMemcpyDeviceToHost);
    
    int icell;  
    double sgm_xtau = SelfFieldModel.sigma_tau; /* [s] */
    double factG = -ring.T0 * bunch->Ib / (bunch->Np * ring.E0 * FGIGA); /* Wakes are in [V/C] or [V/Cm] -> per unit length & for a test charge with Q = 1 C */

    double rftest, taucell;
    const double dTau = SelfFieldModel.dT;
    double active, FacI, ideal;
    double eps_const = ring.Vrf0 / (ring.E0 * FKILO);
    int j;
    active_HC_t * aHC;
    if(track.EnableIdealHC == 0)
      FacI = 0.;
    else
      FacI = eps_const *  ring.HC_k;

    for (icell = 0; icell < SelfFieldModel.Ncell; icell++) {
      taucell = -(-SelfFieldModel.Nsigma + dTau * (icell - 0.5))* sgm_xtau;
      rftest = eps_const * sin(ring.wrf * taucell + ring.phai0) - ring.Urad;
      // adding potential of active HCs
      active = 0.0;
      for(j = 0; j < ring.active_HC_size; j++) {
	aHC = &(ring.active_HC[j]);
	active += aHC->Vpeak/(ring.E0 * FGIGA) * sin(ring.wrf * aHC->nHC * taucell + aHC->phi_aHC + aHC->nHC * kb * 2 * M_PI);
      }
      ideal = FacI * sin(ring.m_aHC * ring.wrf * taucell + ring.m_aHC * ring.phi_n);
      fprintf(fp, "\n %ld   %e   %e   %e    %e   %e   %e   %e   %e",rev, scan_val, taucell,  factG * GL1[icell], lr_wake[icell], rftest, active, ideal, fnp[icell]);
      if(SelfFieldModel.PlaneV > 0)
	fprintf(fp, "   %e   %e", -factG * GlambdaV[icell], dipoleV[icell]);
      if(SelfFieldModel.PlaneH > 0)
	fprintf(fp, "   %e   %e", -factG * GlambdaH[icell], dipoleH[icell]);
    }
    fprintf(fp, "\n");

    cudaError_t e2 = cudaGetLastError();
    if (e2 != cudaSuccess)
      fprintf(stderr, "Error ! CUDA error logging potentials!\n");
  }

  return 1;
}

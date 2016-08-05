#include <stdio.h>
#include <stdbool.h>
#include <mpi.h>
#include "confmpi.h"
#include "input.h"
#include "bunch.h"

#ifdef DEBUG_FLOATING_POINT
#define _GNU_SOURCE
#include <fenv.h>
#include <signal.h>
#endif

#define MSG_SIZE 128

#include "mbtrack-mpi-cuda.h"

/* Global variables */
ring_t ring;
tracking_t track;
grid_t fbii_grid;

int main(int argc, char ** argv)
{

  char mbtrack_inputfile[FILENAME_MAX];
  if(argc == 3) /* Usage: ./executable [input_file] [work_dir] */
  {
    snprintf(mbtrack_inputfile, FILENAME_MAX, "%s", argv[1]); /* Argument 1 */
    snprintf(track.work_path, FILENAME_MAX, "%s", argv[2]); /* Argument 2 */
  }
  else if(argc == 2) /* Usage: ./executable [input_file] */
  {
    snprintf(mbtrack_inputfile, FILENAME_MAX, "%s", argv[1]); /* Argument 1 */
    snprintf(track.work_path, FILENAME_MAX, "work/"); /* Using default value */
  }
  else
  {
    fprintf(stderr, "ERROR ! \n\nUsage:\nmbtrack-mpi inputfile.conf [/path/to/work]\n");
    return -1;
  }
  
  bunch_macroparticle_model_t macrop_model;
  selffield_model_t SelfFieldModel;
  e_beam_t ebeam;

  /* MPI Initialisation */
  ctx_t comm;
  if(!confmpi_init(&argc, &argv, &comm))
  {
    fprintf(stderr, "ERROR: confmpi_init failed\n");
    MPI_Finalize();
    return -1;
  }

#ifdef DEBUG_FLOATING_POINT
  /* Enable floating point exception */
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_UNDERFLOW | FE_OVERFLOW);
#endif
  
  /* Read inputs from file */
  if(!read_input(mbtrack_inputfile, &ring, &track, &ebeam, &macrop_model, &SelfFieldModel))
  {
    fprintf(stderr, "ERROR: read_input failed\n");
    MPI_Finalize();
    return -1;
  }
  
  /* setup parameters */
  setup_ring_parameters(&ring);
  macrop_model_setup_parameters(ring, track, &macrop_model);
  
  /* FBII related parameters */
  fbii_grid_setup_parameters(&fbii_grid, &ebeam, macrop_model);
  
  /* setup ebeam filling */
  if(!e_beam_setup(&track, &ring, &ebeam))
  {
    fprintf(stderr, "ERROR: e_beam_setup failed\n");
    MPI_Finalize();
    return -1;
  }

  /* Check number of processes */
  if (comm.world_size != 2) 
  {
    fprintf(stderr, "Process number must be 2; processes = %d).\n", comm.world_size);
    MPI_Finalize();
    return -1;
  }

  switch (comm.proc_id)
  {
    /* manager specific code */
  case MANAGER_RANK:
    
    setup_tracking_parameters(&ring, &track, &SelfFieldModel, &macrop_model);
    
    /* Check parameters, stop if errors */
    if(!check_parameters(&ring, &track, &ebeam))
    {
      MPI_Abort(MPI_COMM_WORLD, -1);
      return -1;
    }
    #ifdef DEBUG_INPUT
    else
    {
      FILE * log_input = confmpi_log_fopen("input");
      if(log_input == NULL)
      {
	WARNING("confmpi_log_fopen(\"input\")");
      }
      else
      {
	fprint_ring(log_input, ring);
	fprint_tracking(log_input, track);
	fprintf_bunch_macroparticle_model(log_input, macrop_model, track.TrackPlane);
        
	fprint_parameters(log_input, ring, macrop_model);
	fprint_grid(log_input, fbii_grid);
	fprint_e_beam(log_input, ring, ebeam);
          
	if(ring.resonators_size > 0)
	  fprint_resonators(log_input, ring, SelfFieldModel, macrop_model);
        
	fclose(log_input);
      }
    }
    #endif
	
    if(track.BunchModel == MODEL_WEAK)
    {
      manager_weakweak_cuda(ring, track, ebeam, macrop_model, SelfFieldModel);
      ring_destroy(&ring);
      tracking_destroy(&track);
    }
    else if(track.BunchModel == MODEL_STRONG) 
    {
      fprintf(stderr, "Model STRONG is no implemented to run on GPU.\n");
      MPI_Finalize();
      return -1;
    }
        
    break;
    /* end of manager specific code */
    
    /* worker specific code */    
  default:

    if(track.BunchModel == MODEL_WEAK)
    {
      setup_tracking_parameters(&ring, &track, &SelfFieldModel, &macrop_model);
      worker_weakweak_cuda(ring, track, ebeam, macrop_model, SelfFieldModel);
      ring_destroy(&ring);
      tracking_destroy(&track);
    }
    else if(track.BunchModel == MODEL_STRONG)
    {
      fprintf(stderr, "Model STRONG is no implemented to run on GPU.\n");
      MPI_Finalize();
      return -1;
    }
    break;
    /* end of worker specific code */   
  }
  
  MPI_Finalize();
  return 0;
}

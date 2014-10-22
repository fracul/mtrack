#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include "confmpi.h"
#include "types.h"

/* Global variable */
extern tracking_t track;

bool confmpi_init(int * argc, char *** argv, ctx_t * c)
{
  MPI_Init(argc, argv);
  if(!confmpi_ctx_init(c))
  {
    return false;
  }
  return true;
}

bool confmpi_ctx_init(ctx_t * c)
{
  int namelen;
  int * univ_size_p;
  int flag;
  if(MPI_Comm_rank(MPI_COMM_WORLD, &(c->proc_id)) != MPI_SUCCESS)
    return false;
  if(MPI_Get_processor_name((char *) &(c->proc_name), &namelen) != MPI_SUCCESS)
    return false;
  if(MPI_Comm_size(MPI_COMM_WORLD, &(c->world_size)) != MPI_SUCCESS)
    return false;
  if(MPI_Attr_get(MPI_COMM_WORLD, MPI_UNIVERSE_SIZE, &univ_size_p, &flag)
     != MPI_SUCCESS)
     return false;
  if (!flag)
  { 
     printf("This MPI does not support UNIVERSE_SIZE.\n"); 
     return false;
  }
  c->univ_size = *univ_size_p;
  return true;
}

FILE * confmpi_log_fopen(char * filelog_prefix)
{
  if(filelog_prefix == NULL)
  {
    return NULL;
  }
  else
  {
    ctx_t comm;
    confmpi_ctx_init(&comm);
    char log_filename[FILENAME_MAX];
    snprintf(log_filename, FILENAME_MAX, "%s/%s_%d.log", track.work_path, filelog_prefix, comm.proc_id);
    FILE * fp = fopen(log_filename, "a+");
    
    /* Increase buffer size to 64 kio */
    if(fp != NULL)
      setvbuf(fp, NULL, _IOFBF, 65536);
    
    return fp;
  }
}

/**
 * @file
 * MPI related definitions
 */

#ifndef MBTRACK_CONFMPI_H
#define MBTRACK_CONFMPI_H

#include <mpi.h>
#include <stdio.h>
#include <stdbool.h>

#ifdef MPI_MAX_PROCESSOR_NAME
#define MAX_PROCESSOR_NAME MPI_MAX_PROCESSOR_NAME
#else
#define MAX_PROCESSOR_NAME 128
#endif

#define MBTRACK_TAG 4

#define MANAGER_RANK 0

typedef struct ctx
{
  int proc_id;
  char proc_name[MAX_PROCESSOR_NAME];
  int world_size;
  int univ_size;
} ctx_t;

/**
 * Initialize MPI
 */
bool confmpi_init(int * argc, char *** argv, ctx_t * c);

/**
 * Initialize ctx_t structure with host information
 */
bool confmpi_ctx_init(ctx_t * c);

/**
 * Open a log/debug file specific to the process (ex : prefix_0.log)
 */
FILE * confmpi_log_fopen(char * filelog_prefix);

#endif /* MBTRACK_CONFMPI_H */

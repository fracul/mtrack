#include <stdio.h>
#include "stats.h"
#include "types.h"
#include "def.h"

int e_beam_rms_size_fprint(FILE * fp, const e_beam_stats_t * stats, const ring_t ring)
{
  if(fp == NULL || stats == NULL)
  {
    return -1;
  }
  else
  {
    int count = 0;
    
    unsigned int kb;
    for(kb = 0; kb < ring.Nharm; kb++)
    {
      /*** Vertical rms beamsize written out in [mm] ***/
      const int r = fprintf(fp, "%-12.3lf ", FKILO * stats->zz_rms[kb]);
      if(r < 0) return -1;
      else count += r;
    }
    
    const int r = fprintf(fp, "\n");
    if(r < 0) return -1;
    else count += r;
    
    return count;
  }
}

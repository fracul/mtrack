#include <stdlib.h>
#include <math.h>
#include "types.h"
#include "bunch.h"
#include "tracking.h"
#include "transform_weak.h"


int
tracking_scan_step_manager(long int rev, double * scan_val, const tracking_t * track, ring_t * ring, e_beam_t * ebeam, weak_bunch_t * bunches)
{
  if(track->scan == 1)
  {
    double fac;
    ring->Iring += track->scan_step;  
    *scan_val = ring->Iring;
    unsigned int i;
    for(i = 0; i < ring->Nharm; i++)
      if(ebeam->nfFill[i])
      {
        if(track->EnableDiffCurr)
        {
          if(i%2 == 0) fac = track->current_ratio;
          else fac = 1.0 - track->current_ratio;              
          ring->Ibunch[i] += 2 * fac * track->scan_step /((double) ebeam->Nbunch);
          bunches[i].Ib += 2 * fac * track->scan_step /(((double) ebeam->Nbunch) * FKILO);

        }
        else
        {
          ring->Ibunch[i] += track->scan_step /((double) ebeam->Nbunch);
          bunches[i].Ib += track->scan_step /(((double) ebeam->Nbunch) * FKILO);
        }
      }
  }

  else if (track->scan == 2) {
    ring->taue += track->scan_step;
    ring->aexpe = 1.0/ring->taue;
    ring->De    = 2 * ring->aexpe * ring->T0;
    
    *scan_val = ring->taue;
  }
  
  return 0;
}



int
tracking_scan_step_worker(long int rev, double * scan_val, const tracking_t * track, ring_t * ring, e_beam_t * ebeam, weak_bunch_t * bunch)
{
  if(track->scan == 1)
  {
    ring->Iring += track->scan_step;
    *scan_val = ring->Iring;
    bunch->Ib += track->scan_step /(((double) ebeam->Nbunch) * FKILO);
    unsigned int i;
    for(i = 0; i < ring->Nharm; i++)
      if(ebeam->nfFill[i])
        ring->Ibunch[i] += track->scan_step /((double)ebeam->Nbunch);
  }

  else if (track->scan == 2) {
    ring->taue += track->scan_step;
    ring->aexpe = 1.0/ring->taue;
    ring->De    = 2 * ring->aexpe * ring->T0;
    
    *scan_val = ring->taue;
  }
  
  return 0;
}




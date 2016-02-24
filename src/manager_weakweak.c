#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "confmpi.h"
#include "input.h"
#include "types.h"
#include "bunch.h"
#include "fbi.h"
#include "statistics.h"
#include "transform_weak.h"
#include "tracking.h"
#include <time.h>

void manager_weakweak(ring_t ring, const tracking_t track,
                      e_beam_t ebeam,
                      const bunch_macroparticle_model_t bunchModel, const selffield_model_t * SelfFieldModel)
{
  time_t now;
  struct tm * time_now;
  unsigned int kb;
  long int rev;
  
  MPI_Status status;
  
  bunch_CM_history_weak_t CMhist;
  bunch_CM_history_weak_init(&CMhist, (track.NrevTot / track.NrevMon * ebeam.Nbunch));
  bunch_stats_t allstats;
  
  const int Nscan = track.NrevTot / track.NrevMon + track.Nscan;
  int m = 0;
  double scan_val_hist[Nscan];
  double scan_val = track.scan_start;
  unsigned int * branks = (unsigned int *) malloc(ring.Nharm * sizeof(int));
  
  printf("Starting job \"%s\"\n\n", track.jobtitle);  
  printf("Generating initial distribution...\n");
  
  /* TODO FBII like in generate_eCMS: 
   * Creation of local transverse grid points around each xxCM0[icell]...
   * TODO necessary in manager ? maybe in worker, not in manager ?
   */
  
  /* Create bunches, allocate memory */
  weak_bunch_t * bunches = (weak_bunch_t *) malloc(ebeam.Nbunch * sizeof(weak_bunch_t));
  
  for(kb = 0; kb < ebeam.Nbunch; kb++)
  {
    /* Special bunch current configuration happens HERE, by setting Ib */
    unsigned int Np = bunchModel.Np; /* WARNING assuming equal number of particles */
    double bunch_Ib, fac;
    
    if(track.EnableDiffCurr)
    {
      if(kb%2 == 0) fac = track.current_ratio;
      else fac = 1.0 - track.current_ratio;      
      bunch_Ib = 2 * fac * ring.Iring/(((double) ebeam.Nbunch) * FKILO); 
    }
    else
      bunch_Ib = ring.Iring/(((double) ebeam.Nbunch) * FKILO); /* WARNING assuming equal curent in bunches */
      
      MPI_Send(&Np, 1, MPI_UNSIGNED, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(&bunch_Ib, 1, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
    ring.Ibunch[kb] = bunch_Ib * FKILO;
    
    /* If FBII enabled, allocating particles in master */
    if(track.EnableFBII)
    {
      /* Create bunch, allocate Np particles */
      if(!weak_bunch_create(&(bunches[kb]), &track, Np, bunch_Ib, kb, true))
        ERROR("weak_bunch_create", return);
    }
    else /* If FBII disabled, not allocating particles in master */
    {
      /* Create bunch, allocate 0 particles, only statistics */
      if(!weak_bunch_create(&(bunches[kb]), &track, Np, bunch_Ib, kb, false))
        ERROR("weak_bunch_create", return);
    }
  }
  
  /* Broadcast to workers */
  MPI_Bcast((void *) &(ebeam.Nbunch), 1, MPI_INT, MANAGER_RANK, MPI_COMM_WORLD);
  MPI_Bcast(ring.Ibunch, 500, MPI_DOUBLE, MANAGER_RANK, MPI_COMM_WORLD);
  
  unsigned int bnum = 0;
  for(kb = 0; kb < ebeam.Nbunch; kb++)
  {
    if (ebeam.nfFill[kb]) {
      /* Send bucket number to each worker */
      bnum++;
      branks[kb] = bnum;
      MPI_Send(&kb, 1, MPI_INT, bnum, MBTRACK_TAG, MPI_COMM_WORLD);
    }
  }
  
  /* Recieve, calculate and output initial distribution statistics over all bunches */
  char filename_bstats[FILENAME_MAX] = "";
  snprintf(filename_bstats, FILENAME_MAX, "%s/mean_bunch_all.dat", track.work_path);
  FILE * bstats_fp = fopen(filename_bstats, "w+");
  if(bstats_fp == NULL)
    ERROR("fopen_stat_all", return);
  weak_bunch_reset_statistics(&allstats);
  for(kb = 0; kb < ebeam.Nbunch; kb++)
  {
    weak_bunch_t * bunch = &(bunches[kb]);
    bunch_stats_t * bstats = &(bunch->stats);
    MPI_Recv(&(bstats->pos), 3, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(&(bstats->pos_sigma), 3, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(&(bstats->slope), 3, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(&(bstats->slope_sigma), 3, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    
    weak_bunch_add_statistics(&allstats, bstats);
  }
  
  weak_bunch_update_statistics(&allstats, ebeam.Nbunch);
  fprintf(bstats_fp, " # bunch statistics, AVEraged over all bunches & turns, Ib = %g A ;\n", bunches[0].Ib);
  fprintf(bstats_fp, " # turn - scan_var - AVE_CM_LON - SIG_CM_LON - AVE_delE_LON - SIG_delE_LON -  AVE_sigmaL_LON - SIG_sigmaL_LON - AVE_sigmaE_LON - SIG_sigmaE_LON... (_HOR) ... (_VER)\n # ");
  bunch_stats_fprintf(bstats_fp, -100, &allstats, 0.0);
  fprintf(bstats_fp, "\n");
  weak_bunch_reset_statistics(&allstats);
  
  
  if(ring.longrange_resonators_size > 0)
  {
    unsigned int i;
    double ttrash[SelfFieldModel->Ncell];
    for(i = 0; i < ring.Nharm; i++)
      if(ebeam.nfFill[i])           
      {
        /* Broadcast from workers to workers, manager does not nead the info -> trash */
        MPI_Bcast(ttrash, SelfFieldModel->Ncell, MPI_DOUBLE, branks[i], MPI_COMM_WORLD);        
      }
  }
  
  if(track.EnableRW_long > 0)
  {
    unsigned int i;
    double trash[2];
    for(i = 0; i < ring.Nharm; i++)
      if(ebeam.nfFill[i])
      {
        /* Broadcast from workers to workers, manager does not nead the info -> trash */
        MPI_Bcast(trash, 2, MPI_DOUBLE, branks[i], MPI_COMM_WORLD);
      }
  }
  
  printf("Starting multi bunch tracking...\n");
  
  /* Counting turns for statistics */ 
  long int Nstat = 0;
  /* ----------------------------------------------------------------------------------- */
  /* Tracking for NrevTot turns */
  for(rev = 0; rev < track.NrevTot; rev++)
  {   
    /* Terminal output */
    if(rev%5000 == 0)
    {
      if(rev+4999 >= track.NrevTot)
        printf("Transformation for turn %3ld to %3ld", rev, track.NrevTot-1);
      else
        printf("Transformation for turn %3ld to %3ld", rev, rev+4999);
      fflush(stdout);
    }
    /* Workers do selffield effect */    
    /* Workers do optics & HC trafo */
    
    /*------------------------------------------------------------------------------
     *          Fast Beam-Ion Interaction: ((Introduced in August 2007))
     *----------------------------------------------------------------------------*/
    
    /* FBII initialization TODO (should avoid reset if UNIFORM filling ?) */
    /* FBII communication TODO */
    /* FBII communication TODO */
    /* FBII writeout TODO */    
    /* TODO
     *    if(writeoutVCMmotionsPM2(track.work_path, rev, bunches, &ebeam) < 0)
     *    {
     *      fprintf(stderr, "Error: writeoutVCMmotionsPM2 failed\n");
  }
  */    
    /* TODO writeout_rms_Vbeamsize */
    
    
    if(track.EnableRW_long > 0)
    {
      unsigned int i;
      double trash[2];
      for(i = 0; i < ring.Nharm; i++)
        if(ebeam.nfFill[i])
        {
          /* Broadcast from workers to workers, manager does not nead the info -> trash */
          MPI_Bcast(trash, 2, MPI_DOUBLE, branks[i], MPI_COMM_WORLD);
        }
    }    
    
    /* Recieve statistics from each worker, update allstats and RW history */
    for(kb = 0; kb < ebeam.Nbunch; kb++)
    {
      weak_bunch_t * bunch = &(bunches[kb]);
      bunch_stats_t * bstats = &(bunch->stats);
      MPI_Recv(&(bstats->pos), 3, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
      MPI_Recv(&(bstats->pos_sigma), 3, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
      MPI_Recv(&(bstats->slope), 3, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
      MPI_Recv(&(bstats->slope_sigma), 3, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
      
      /* Every NrevMon turn, recieve more statistics */
      if((rev+1)%(track.NrevMon) == 0 && rev+1>=track.NrevOutputStart+(rev/track.NrevScan)*track.NrevScan)// && track.EnableAmpinv_out)
      {
        const long int irevmon = rev/track.NrevMon;
        MPI_Recv(&(CMhist.ampinv[irevmon * ebeam.Nbunch + kb]), 3, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&(CMhist.ampinv_cm[irevmon * ebeam.Nbunch + kb]), 3, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);       
      }
      weak_bunch_add_statistics(&allstats, bstats);
    }
    Nstat++;
    
    if(ring.longrange_resonators_size > 0)
    {
      unsigned int i;
      double ttrash[SelfFieldModel->Ncell];
      for(i = 0; i < ring.Nharm; i++)
        if(ebeam.nfFill[i])
        {
          /* Broadcast from workers to workers, manager does not nead the info -> trash */
          MPI_Bcast(ttrash, SelfFieldModel->Ncell, MPI_DOUBLE, branks[i], MPI_COMM_WORLD);
        }
    }
    
    if((Nstat % track.NrevMon == 0 || (rev+1) % track.NrevScan == 0) && rev+1>=track.NrevOutputStart+(rev/track.NrevScan)*track.NrevScan)
    { /* Outputfile ampinv and bunch stats */
    Nstat *= ebeam.Nbunch;
    weak_bunch_update_statistics(&allstats, Nstat);
    scan_val_hist[m] = scan_val;
    m++;
    if(track.EnableAmpinv_out) bunch_weak_writeout_ampinv(rev, rev/track.NrevMon, &CMhist, ebeam);
    bunch_stats_fprintf(bstats_fp, rev, &allstats, scan_val);
    weak_bunch_reset_statistics(&allstats);
    Nstat = 0;
    }
    
    /* Update scan value if rev % NrevScan */
    if((rev+1) % track.NrevScan == 0 && rev > 0)
      tracking_scan_step_manager(rev, &scan_val, &track, &ring, &ebeam, bunches);
    
    if(rev%500 == 499)
      printf("."); fflush(stdout);
    if(rev%5000 == 4999 || rev == track.NrevTot-1)
    {    
      time(&now);
      time_now = localtime(&now);
      printf(" done (at time: %d:%d:%d)\n", time_now->tm_hour, time_now->tm_min, time_now->tm_sec);
    }
    
  } /* end tracking ---------------------------------------------------------------------------*/
  
  /* Recieve, calculate and output end distribution statistics over all bunches */
  if(track.scan == 0)
  {
    char filename_bstats_end[FILENAME_MAX] = "";
    snprintf(filename_bstats_end, FILENAME_MAX, "%s/mean_bunch_end.dat", track.work_path);
    FILE * bstats_end_fp = fopen(filename_bstats_end, "w+");
    if(bstats_end_fp == NULL)
      ERROR("fopen_stat_end", return);
    fprintf(bstats_end_fp, " # bunch statistics end of run, turns = %i; Ib = %g A ;\n", rev, bunches[0].Ib);
    fprintf(bstats_end_fp, " # bunch_num    CM    bunch_length    energy_dev    energy_spread;\n");
    for(kb = 0; kb < ebeam.Nbunch; kb++)
    {
      weak_bunch_t * bunch = &(bunches[kb]);
      bunch_stats_t * bstats = &(bunch->stats);
      MPI_Recv(&(bstats->pos), 3, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
      MPI_Recv(&(bstats->pos_sigma), 3, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
      MPI_Recv(&(bstats->slope), 3, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
      MPI_Recv(&(bstats->slope_sigma), 3, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
      
      fprintf(bstats_end_fp, " %d     %e    %e    %e    %e\n", kb, bstats->pos.xtau, bstats->pos_sigma.xtau, bstats->slope.xtau, bstats->slope_sigma.xtau);
    }
    fclose(bstats_end_fp);
  }
  
  /* Writeout ampinv history */
  //if(track.EnableAmpinv_out)
  //{
  if(track.TrackPlane[HOR])
    weak_bunch_writeout_mean_ampinv(track.NrevTot, track.NrevMon, &CMhist, ebeam, HOR, scan_val_hist);
  if(track.TrackPlane[VER])
    weak_bunch_writeout_mean_ampinv(track.NrevTot, track.NrevMon, &CMhist, ebeam, VER, scan_val_hist);
  //}
  
  /* free memory */
  for(kb = 0; kb < ebeam.Nbunch; kb++)
  {
    weak_bunch_destroy(&(bunches[kb]));
  }
  free(bunches);
  fclose(bstats_fp);
  bunch_CM_history_weak_destroy(&CMhist);
}

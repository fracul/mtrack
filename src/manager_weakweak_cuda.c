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

void manager_weakweak_cuda(ring_t ring, const tracking_t track,
			   e_beam_t ebeam, const bunch_macroparticle_model_t bunchModel, 
			   selffield_model_t SelfFieldModel)
{

  time_t now;
  struct tm * time_now;
  unsigned int kb;
  long int rev;

  struct timeval time_start, time_end;
  gettimeofday(&time_start, NULL);

  cyclic_array_t moments_history;

  MPI_Status status;

  bunch_CM_history_weak_t CMhist;
  bunch_CM_history_weak_init(&CMhist, (track.NrevTot / track.NrevMon * ebeam.Nbunch));
  bunch_stats_t allstats;

  int m = 0;
  int mm = 0;
  const int Nscan = track.NrevTot / track.NrevMon + track.Nscan;
  double scan_val_hist[Nscan];
  double scan_val = track.scan_start;
  
  printf("Starting job \"%s\"\n\n", track.jobtitle);  
  printf("Generating initial distribution...\n");

  /* Create bunches, allocate memory */
  weak_bunch_t * bunches = (weak_bunch_t *) malloc(ebeam.Nbunch * sizeof(weak_bunch_t));

  /* WARNING assuming equal number of particles */
  unsigned int Np = bunchModel.Np; /* WARNING assuming equal number of particles */

  /* Special bunch current configuration happens HERE, by setting Ib */
  for(kb = 0; kb < ebeam.Nbunch; kb++)
  {
    double bunch_Ib, fac;

    if(track.EnableDiffCurr)
    {
      if(kb%2 == 0) fac = track.current_ratio;
      else fac = 1.0 - track.current_ratio;      
      bunch_Ib = 2 * fac * ring.Iring/(((double) ebeam.Nbunch) * FKILO); 
    }
    else
      bunch_Ib = ring.Iring/(((double) ebeam.Nbunch) * FKILO); 
     /* WARNING assuming equal curent in bunches */

    ring.Ibunch[kb] = bunch_Ib * FKILO;

    /* Create bunch, allocate memory */
    if(!weak_bunch_create(&bunches[kb], &track, Np, bunch_Ib, kb, true))
      ERROR("weak_bunch_create", return);
    /* Generate bunch distribution, with offsets */
    if(!weak_generate_bunch_distribution(&bunches[kb], bunchModel, track.TrackPlane))
      ERROR("weak_bunch_generate_bunch_distribution", return);  
    /* Initialize selffield model, and allocate memory */
    if(!selffield_model_init(&SelfFieldModel, bunchModel.pos_sgm.xtau, &bunches[kb]))
      ERROR("selffield_model_init", return);
  }

  printf("Number of particles: %d\n", Np);

  /* Construct Greensfunctions for all given resonators */
  if(ring.resonators_size > 0)
    construct_greensfunc_resonator(&SelfFieldModel, &ring);

  /* Adding RW wake to SelfFieldModel */
  if(track.EnableRW_short)
    construct_greensfunc_RW(&SelfFieldModel, &ring);  

  FILE *wakes_fp = NULL;
  char filename_wakes[FILENAME_MAX] = "";
  snprintf(filename_wakes, FILENAME_MAX, "%s/wakes.dat", track.work_path);  
  wakes_fp = fopen(filename_wakes, "w+");
  fprintf(wakes_fp, "# bin (%.2e s)   total wake L     - V -      - H -\n", 
	  SelfFieldModel.dT*SelfFieldModel.sigma_tau);

  int nc;
  for (nc = 0; nc < SelfFieldModel.Ncell; nc++)
  {
    fprintf(wakes_fp, "%d  %e   %e   %e\n", nc, SelfFieldModel.Gl[nc], SelfFieldModel.Gv[nc], 
	    SelfFieldModel.Gh[nc]);
  }

  /* Update statistics with this new distribution */
  weak_bunch_reset_statistics(&allstats);
  for (kb = 0; kb < ebeam.Nbunch; kb++)
  {
    weak_bunch_calc_statistics(&bunches[kb], &track);
    weak_bunch_add_statistics(&allstats, &bunches[kb].stats);
  }
  weak_bunch_update_statistics(&allstats, ebeam.Nbunch);

  char filename_bstats[FILENAME_MAX] = "";
  snprintf(filename_bstats, FILENAME_MAX, "%s/mean_bunch_all.dat", track.work_path);
  FILE * bstats_fp_all = fopen(filename_bstats, "w+");
  fprintf(bstats_fp_all, " # bunch statistics, AVEraged over all bunches & turns, Ib = %g A ;\n", 
	  bunches[0].Ib);
  fprintf(bstats_fp_all, " # turn - scan_var - AVE_CM_LON - SIG_CM_LON - AVE_delE_LON - SIG_delE_LON -  AVE_sigmaL_LON - SIG_sigmaL_LON - AVE_sigmaE_LON - SIG_sigmaE_LON... (_HOR) ... (_VER)\n # ");
  bunch_stats_fprintf(bstats_fp_all, -100, &allstats, 0.0);
  fprintf(bstats_fp_all, "\n");
  weak_bunch_reset_statistics(&allstats);

  /* initialise cyclic array */
  int fnp_size = (int)SelfFieldModel.Ncell * ring.Nharm;
  double * fnp_ring;
  int phasor_size = ebeam.Nbunch*2*ring.longrange_resonators_size;
  double phasor_end[phasor_size]; // real & imag. part of phasor at end of turn
  if (ring.longrange_resonators_size > 0)
  {
    int Nbin = SelfFieldModel.Ncell;
    fnp_ring = (double *) calloc(fnp_size, sizeof(double));

    for (kb = 0; kb < ebeam.Nbunch; kb++)
    {
      if(!fnp_ring_update(&ring, fnp_ring, &SelfFieldModel, &bunches[kb], kb))
	ERROR("fnp_ring_init", return);

    }

    for (kb = 0; kb < ebeam.Nbunch; kb++)
    {
      int offset = kb*2*ring.longrange_resonators_size;
      wake_phasor_init(&ring, fnp_ring, &SelfFieldModel, &bunches[kb], 
		       kb, &phasor_end[offset], &ebeam);
    }
  }

  /* initialise cyclic array */
  cyclic_array_t dipole_RW;
  if (track.EnableRW_long > 0)
  {
    if (!cyclic_array_create(&dipole_RW, (track.Nmlt+2)*ring.Nharm, 2))
      ERROR("cyclic_array_create", return);
    
    
    double bunch_dipole[2];
    for (kb = 0; kb < ebeam.Nbunch; kb++)
    {
      bunch_dipole[0] = bunches[kb].stats.pos.x;
      bunch_dipole[1] = bunches[kb].stats.pos.z;
      cyclic_array_set(kb, &bunch_dipole[0], &dipole_RW);
    }

    unsigned int i, j;
    for(i = 0; i < ring.Nharm; i++) {
      if(ebeam.nfFill[i]) {
        double * mom = cyclic_array_access(i, &dipole_RW);
        for(j = 1; j < track.Nmlt; j++)
          cyclic_array_set((j*ring.Nharm + i), mom, &dipole_RW);
      }
    }
  }
  printf("Starting multi bunch tracking...\n");

  /* Initialize model for FBI interaction */
  fbii_model_t fbii_model;
  fbii_model_init(&fbii_model);
  fbii_model.factBE = sqrt(2*(bunchModel.sgm_xx*bunchModel.sgm_xx-bunchModel.sgm_zz*bunchModel.sgm_zz))/FKILO;  

  FILE *bstats_fp[ebeam.Nbunch];
  FILE *trafo_fp[ebeam.Nbunch];
  for (kb = 0; kb < ebeam.Nbunch; kb++) 
  {
    bstats_fp[kb] = NULL;
    trafo_fp[kb] = NULL;
  }

  /* Open output files */
  for (kb = 0; kb < ebeam.Nbunch; kb++) 
  {

    if (bunches[kb].kb_out == 1) {
      char filename_bstats[FILENAME_MAX] = "";
      snprintf(filename_bstats, FILENAME_MAX, "%s/mean_bunch_%d.dat", track.work_path, kb);  
      bstats_fp[kb] = fopen(filename_bstats, "w+");
      if(bstats_fp[kb] == NULL)
	ERROR("fopen_stat_bunch0", return);
      worker_weak_stat_output(bstats_fp[kb], &bunches[kb].stats, bunches[kb].Ib);
      weak_bunch_reset_statistics(&bunches[kb].stats);
  
      if (track.NrevPotentialsOut<=track.NrevTot) {
	char filename_trafo[FILENAME_MAX] = "";
	snprintf(filename_trafo, FILENAME_MAX, "%s/potentials_bunch_%d.dat", track.work_path, kb);  
	trafo_fp[kb] = fopen(filename_trafo, "w+");
	if(trafo_fp[kb] == NULL)
	  ERROR("fopen_potentials_bunch0", return);
	fprintf(trafo_fp[kb], " # potentials and distributions, Ib = %g A ;\n", bunches[kb].Ib);
	fprintf(trafo_fp[kb], " # rev  scan_val         tau         LON_wake_pot       HC_pot        rf_pot      active_cav      ideal_HC        bunch_shape");
	if(SelfFieldModel.PlaneV > 0)
	  fprintf(trafo_fp[kb], "     VER_wake_pot     VER_dipole_mom");
	if(SelfFieldModel.PlaneH > 0)
	  fprintf(trafo_fp[kb], "     HOR_wake_pot     HOR_dipole_mom");
	fprintf(trafo_fp[kb], "\n");
      }
    }
  }

  /* Counting turns for statistics */
  int long Nstat = 0;
  int iseed = bunchModel.iseed[LON];
  //setup device memory, init random numbers, page lock host memory
  setup_cuda(ebeam.Nbunch, Np, SelfFieldModel.Ncell, fnp_ring);
  transfer_ring(&ring, &bunchModel, &track, &SelfFieldModel, Np, iseed);
  transfer_ebeam(&ebeam);
  transfer_phasor(ebeam.Nbunch, ring.longrange_resonators_size, phasor_end);

  if (track.EnableRW_long > 0)
    initialize_cyclic_array_cuda(dipole_RW);

  for (kb = 0; kb < ebeam.Nbunch; kb++) {
    allocate_bunch(&bunches[kb], kb);
    transfer_bunch_to_device(&bunches[kb], kb);
  }

  printf("DEBUG: tracking %d bunches for %d turns\n", ebeam.Nbunch, track.NrevTot);

  /* ----------------------------------------------------------------------------------------- */
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

    //transfer bunch to device
    //for (kb = 0; kb < ebeam.Nbunch; kb++)
    //  transfer_bunch_to_device(&bunches[kb], kb);

    if (ring.longrange_resonators_size > 0) {
      for (kb = 0; kb < ebeam.Nbunch; kb++)
	fnp_ring_update_cuda(&bunches[kb], SelfFieldModel, kb);

      for (kb = 0; kb < ebeam.Nbunch; kb++)
	construct_wake_phasor_cuda(Np, SelfFieldModel.Ncell, kb, ring.longrange_resonators_size);
    }

    for (kb = 0; kb < ebeam.Nbunch; kb++)
    {
      /* Resonator selffield transformation includes effect of harmonic cavity and RW */
      /* PlaneL etc give information if a resonator in this plane is given AND if plane is tracked */ 
      if(SelfFieldModel.PlaneL + SelfFieldModel.PlaneV + SelfFieldModel.PlaneH 
	 + ring.longrange_resonators_size > 0)
      {
	transform_weak_bunch_selffield_cuda(&bunches[kb], SelfFieldModel, kb, trafo_fp[kb], 
					     rev, scan_val, ring.longrange_resonators_size);
      }

      /* Optics transformaiton including quantum excitation & radiation 
       * damping and active or passive HC.
       * Long. plane is always tracked, hor. and vert. has to be enabled */
      if (transform_weak_bunch_optic_cuda(&bunches[kb], &ring, &bunchModel, iseed+rev, kb) < 0) {
	fprintf(stderr, "Error ! transform_optic_LON_CUDA failed\n");   
	MPI_Abort(MPI_COMM_WORLD, 1);
      }
      

      //transfer bunch back from device
      transfer_bunch_from_device(&bunches[kb], kb);
    }
    sync_device();
    /* end tracking for one bunch */

    /* Update the statistics every turn */
    Nstat++;
    for (kb = 0; kb < ebeam.Nbunch; kb++) {
      bunch_stats_t * bstats = &(bunches[kb].stats);

      weak_bunch_calc_statistics(&bunches[kb], &track);
      
      if(bunches[kb].kb_out == 1)
      {
	if(Nstat % track.NrevMon == 0 || (rev+1) % track.NrevScan == 0)  
	{ /* Outputfile ampinv and bunch stats */
	  weak_bunch_update_statistics(bstats, Nstat);
	  bunch_stats_fprintf(bstats_fp[kb], rev, bstats, scan_val);
	  weak_bunch_reset_statistics(bstats);
	}
      }
    
      if (rev%(track.NrevMon) == 0)
	weak_bunch_calc_ampinv(&bunches[kb], &track);
      weak_bunch_add_statistics(&allstats, &bunches[kb].stats);
    }

    /* TODO update cyclic array with actual statistics */
    //update the cyclic array on the CPU side, since the statistics are on the CPU
    
    if (track.EnableRW_long) {
      for (kb = 0; kb < ebeam.Nbunch; kb++) {
	bunch_stats_t * bstats = &(bunches[kb].stats);
	double *mom = cyclic_array_access(-rev*ring.Nharm + kb, &dipole_RW);
	mom[0] = bstats->pos.x;
	mom[1] = bstats->pos.z;
      }
      update_cyclic_array_cuda(&dipole_RW);
    }
    

    /* TODO include interbunch (long-range resistive-wall interactions */
    if (track.EnableRW_long) {
      for (kb = 0; kb < ebeam.Nbunch; kb++) {
	const int bpos = -rev*ring.Nharm + kb;  //Actual position of bunch kb in cyclic array
	if (track.TrackPlane[HOR])
	  transform_bunch_RW_longrange_cyclic_cuda(m, bpos, HOR, &bunchModel, &ring, &track, kb, bunches[kb].Np, &dipole_RW);
	if (track.TrackPlane[VER])
	  transform_bunch_RW_longrange_cyclic_cuda(m, bpos, VER, &bunchModel, &ring, &track, kb, bunches[kb].Np, &dipole_RW);

	//transfer bunch back from device
	transfer_bunch_from_device(&bunches[kb], kb);
      }
      sync_device();
      m++;
    }

    /* write out statistics */
    if(Nstat % track.NrevMon == 0 || (rev+1) % track.NrevScan == 0)
    { /* Outputfile ampinv and bunch stats */
      Nstat *= ebeam.Nbunch;
      weak_bunch_update_statistics(&allstats, Nstat);
      scan_val_hist[mm] = scan_val;
      mm++;
      if(track.EnableAmpinv_out) 
	bunch_weak_writeout_ampinv(rev, rev/track.NrevMon, &CMhist, ebeam);
      bunch_stats_fprintf(bstats_fp_all, rev, &allstats, scan_val);
      weak_bunch_reset_statistics(&allstats);
      Nstat = 0;
    }

    if((rev+1) % track.NrevScan == 0) {
      tracking_scan_step_worker(rev, &scan_val, &track, &ring, &ebeam, &bunches[kb]);
      m = 0;
    }

    if(rev%500 == 499)
      printf("."); fflush(stdout);
    if(rev%5000 == 4999 || rev == track.NrevTot-1)
    {    
      time(&now);
      time_now = localtime(&now);
      printf(" done (at time: %d:%d:%d)\n", time_now->tm_hour, time_now->tm_min, time_now->tm_sec);
    }

  }
  /* end tracking --------------------------------------------------------------- */

  if(track.scan == 0)
  {
    char filename_bstats_end[FILENAME_MAX] = "";
    snprintf(filename_bstats_end, FILENAME_MAX, "%s/mean_bunch_end.dat", track.work_path);
    FILE * bstats_end_fp = fopen(filename_bstats_end, "w+");
    if(bstats_end_fp == NULL)
      ERROR("fopen_stat_end", return);
    fprintf(bstats_end_fp, " # bunch statistics end of run, turns = %i; Ib = %g A ;\n", rev, bunches[0].Ib);
    fprintf(bstats_end_fp, " # bunch_num    CM    bunch_length    energy_dev    energy_spread;\n");


    for (kb = 0; kb < ebeam.Nbunch; kb++) {

      if(ebeam.nfFill[kb] && (bunches[kb].N_trash_low != 0 || bunches[kb].N_trash_high != 0))
      {
	printf("\n WARNING - Binning for selffields bunch #%d:\n %d particles below first bin, %d particles above last bin\n", kb, bunches[kb].N_trash_low, bunches[kb].N_trash_high);
      }

      bunch_stats_t * bstats = &(bunches[kb].stats);
      weak_bunch_calc_statistics(&bunches[kb], &track);
      fprintf(bstats_end_fp, " %d     %e    %e    %e    %e\n", kb, bstats->pos.xtau, bstats->pos_sigma.xtau, bstats->slope.xtau, bstats->slope_sigma.xtau);
    }

  }
  
  free_cuda();
  for (kb = 0; kb < ebeam.Nbunch; kb++)
    free_bunch(&bunches[kb], kb);
  selffield_model_destroy(&SelfFieldModel);
  for (kb = 0; kb < ebeam.Nbunch; kb++)
    weak_bunch_destroy(&bunches[kb]);
  if(ring.longrange_resonators_size > 0)
  {    
    free(fnp_ring);
    fnp_ring = NULL;
  }

  if(track.EnableRW_long > 0) {
    free_cyclic_array_cuda(&dipole_RW);
    cyclic_array_destroy(&dipole_RW);
  }

  gettimeofday(&time_end, NULL);

  double timervalue = ( (time_end.tv_sec - time_start.tv_sec) * 1000000 + 
			(time_end.tv_usec - time_start.tv_usec)) * 1e-6;

  printf("Total execution time: %f\n", timervalue);

}

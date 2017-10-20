#include <math.h>
#include <mpi.h>
#include "confmpi.h"
#include "types.h"
#include "bunch.h"
#include "fbi.h"
#include "input.h"
#include "cyclic_array.h"
#include "statistics.h"
#include "transform_weak.h"
#include "tune.h"
#include "tracking.h"
#include "feedback_rf.h"


/* Global variables */
extern grid_t fbii_grid;

void worker_weakweak(ring_t ring, const tracking_t track, e_beam_t ebeam,
                     const bunch_macroparticle_model_t bunchModel0,
                     selffield_model_t SelfFieldModel)
{
  long int rev;
  
  MPI_Status status;
  
  bunch_macroparticle_model_t bunchModel = bunchModel0;
  weak_bunch_t bunch;
  bunch_stats_t * bstats = &(bunch.stats);
  FILE * bstats_fp = NULL;
  FILE * trafo_fp = NULL;
  FILE * wakes_fp = NULL;
  cyclic_array_t moments_history;
  cyclic_array_t dipole_RW;
  double phasor_end[2*ring.longrange_resonators_size[LON]]; // real & imag. part of phasor at end of turn
  double phasor_end_HOR[2*ring.longrange_resonators_size[HOR]]; // real & imag. part of phasor at end of turn
  double phasor_end_VER[2*ring.longrange_resonators_size[VER]]; // real & imag. part of phasor at end of turn
  int lr_res_sizetot = ring.longrange_resonators_size[LON] + ring.longrange_resonators_size[HOR] + ring.longrange_resonators_size[VER];
  double * fnp_ring; // # macroparticles per bin of last turn, for passive HC phasor
  double * fnp_HOR; // # dipole moment in bins of last turn, for horizontal HOMs
  double * fnp_VER; // # dipole moment in bins of last turn, for vertical HOMs
  double scan_val = track.scan_start;
  unsigned int * branks = (unsigned int *) malloc(ring.Nharm*sizeof(int));
  unsigned int bnum=0;
  unsigned int i;
  //double ring.rf_feedback->phi0_design, ring.rf_feedback->vrf_design;
  for (i=0; i<ring.Nharm; i++)
      if (ebeam.nfFill[i])
      {
	bnum++;
	branks[i]=bnum;
      }
  
  fbii_model_t fbii_model;
  /*ion_beam_local_t ionbeam;*/ /* TODO Ion beam for single bunch transformation */
  
  /* Recieve bunch parameters */
  unsigned int Np;
  double bunch_I;
  MPI_Recv(&Np, 1, MPI_UNSIGNED, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
  MPI_Recv(&bunch_I, 1, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
  bunchModel.Np = Np;
  bunchModel.fNp = (double) Np;
  
  MPI_Bcast((void *) &(ebeam.Nbunch), 1, MPI_INT, MANAGER_RANK, MPI_COMM_WORLD);
  MPI_Bcast(ring.Ibunch, ring.Nharm, MPI_DOUBLE, MANAGER_RANK, MPI_COMM_WORLD);
  
  /* Recieve bunch number from manager */
  int kb;
  MPI_Recv(&kb, 1, MPI_INT, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
  
  /* Create bunch, allocate memory */
  if(!weak_bunch_create(&bunch, &track, Np, bunch_I, kb, true))
    ERROR("weak_bunch_create", return);
  /* Generate bunch distribution, with offsets */
  if(!weak_generate_bunch_distribution(&bunch, bunchModel, track.TrackPlane))
    ERROR("weak_bunch_generate_bunch_distribution", return);  
  /* Initialize selffield model, and allocate memory */
  if(!selffield_model_init(&SelfFieldModel, bunchModel.pos_sgm.xtau, &bunch))
    ERROR("selffield_model_init", return);
  /* Construct Greensfunctions for all given resonators */
  
  if(ring.resonators_size > 0)
    construct_greensfunc_resonator(&SelfFieldModel, &ring);
  
  /* Adding RW wake to SelfFieldModel */
  if(track.EnableRW_short)
    construct_greensfunc_RW(&SelfFieldModel, &ring);

  if(bunch.kb_out == 1)
  {
    char filename_wakes[FILENAME_MAX] = "";
    snprintf(filename_wakes, FILENAME_MAX, "%s/wakes.dat", track.work_path);  
    wakes_fp = fopen(filename_wakes, "w+");
    fprintf(wakes_fp, "# bin (%.2e s)   total wake L     - V -      - H -\n", SelfFieldModel.dT*SelfFieldModel.sigma_tau);
    int nc;
    for (nc = 0; nc < SelfFieldModel.Ncell; nc++)
    {
      fprintf(wakes_fp, "%d  %e   %e   %e\n", nc, SelfFieldModel.Gl[nc], SelfFieldModel.Gv[nc], SelfFieldModel.Gh[nc]);
    }
  }
  
  /* Update statistics with this new distribution and send to manager */
  weak_bunch_calc_statistics(&bunch, &track);
  
  MPI_Send(&(bstats->pos), 3, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
  MPI_Send(&(bstats->pos_sigma), 3, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
  MPI_Send(&(bstats->slope), 3, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
  MPI_Send(&(bstats->slope_sigma), 3, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD); 
  
  /* initialise cyclic array*/
  int Nbin = SelfFieldModel.Ncell;
  if(ring.longrange_resonators_size[LON] > 0)
  {
    if (ring.has_rf_feedback)
      rffb_init(ring.rf_feedback);
  
    int Nbin = SelfFieldModel.Ncell;
    fnp_ring =
    (double *) calloc((int)Nbin * ring.Nharm, sizeof(double));
    
    if(!fnp_ring_update(&ring, fnp_ring, &SelfFieldModel, &bunch, kb, LON))
      ERROR("fnp_ring_init", return);
    
    unsigned int i;    
    for(i = 0; i < ring.Nharm; i++)
        if(ebeam.nfFill[i])
	{
	  /* if kb = i: bunch_moments are sent, if not recieved and stored fnp_ring */
	  MPI_Bcast(&fnp_ring[i*Nbin], Nbin, MPI_DOUBLE, branks[i], MPI_COMM_WORLD);
	}
      
    wake_phasor_init(&ring, fnp_ring, &SelfFieldModel, &bunch, kb, phasor_end, &ebeam, LON);

    /* feedback: calculate current voltage and phase and try to restore to desired value after longrange wake initialisation*/
    if (ring.has_rf_feedback) {
      double vmbar, phibar;
      rffb_calc_mean_voltage_phase(ring.rf_feedback,&vmbar,&phibar);
      vmbar = vmbar*ring.E0*FKILO;
      //ring.rf_feedback->phi0_design = 1*ring.phai0;
      //ring.rf_feedback->vrf_design = 1*ring.Vrf0;
      double tmp_numrtor = ring.rf_feedback->vrf_design*sin(ring.rf_feedback->phi0_design)-vmbar*sin(phibar);
      //ring.phai0 = atan2(tmp_numrtor*ring.wrf,ring.rf_feedback->vrf_design*ring.wrf*cos(ring.rf_feedback->phi0_design)-vmbar*ring.longrange_resonators[0].wr*cos(phibar));
      ring.phai0 = atan2(tmp_numrtor,ring.rf_feedback->vrf_design*cos(ring.rf_feedback->phi0_design)-vmbar*cos(phibar));
      ring.Vrf0 = tmp_numrtor/sin(ring.phai0);
    }
  }
  if(track.TrackPlane[HOR] && ring.longrange_resonators_size[HOR] > 0)
  {
    fnp_HOR =
    (double *) calloc((int)Nbin * ring.Nharm, sizeof(double));
    
    if(!fnp_ring_update(&ring, fnp_HOR, &SelfFieldModel, &bunch, kb, HOR))
      ERROR("fnp_ring_init", return);
    
    unsigned int i;    
    for(i = 0; i < ring.Nharm; i++)
        if(ebeam.nfFill[i])
	{
	  /* if kb = i: bunch_moments are sent, if not recieved and stored fnp_ring */
	  MPI_Bcast(&fnp_HOR[i*Nbin], Nbin, MPI_DOUBLE, branks[i], MPI_COMM_WORLD);
	}
      
    wake_phasor_init(&ring, fnp_HOR, &SelfFieldModel, &bunch, kb, phasor_end_HOR, &ebeam, HOR);
  }
  if(track.TrackPlane[VER] && ring.longrange_resonators_size[VER] > 0)
  {
    fnp_VER =
    (double *) calloc((int)Nbin * ring.Nharm, sizeof(double));
    
    if(!fnp_ring_update(&ring, fnp_VER, &SelfFieldModel, &bunch, kb, VER))
      ERROR("fnp_ring_init", return);
    
    unsigned int i;    
    for(i = 0; i < ring.Nharm; i++)
        if(ebeam.nfFill[i])
	{
	  /* if kb = i: bunch_moments are sent, if not recieved and stored fnp_ring */
	  MPI_Bcast(&fnp_VER[i*Nbin], Nbin, MPI_DOUBLE, branks[i], MPI_COMM_WORLD);
	}
      
    wake_phasor_init(&ring, fnp_VER, &SelfFieldModel, &bunch, kb, phasor_end_VER, &ebeam, VER);
  }
  
  
  /* initialise cyclic array*/
  if(track.EnableRW_long > 0)
  {
    /* analog to longrange_resonators */
    if(!cyclic_array_create(&dipole_RW, (track.Nmlt+2)*ring.Nharm, 2))
      ERROR("cyclic_array_create", return);
    double bunch_dipole[2] = {bstats->pos.x, bstats->pos.z};
    cyclic_array_set(kb, &bunch_dipole[0], &dipole_RW);
    unsigned int i, j;
    for(i = 0; i < ring.Nharm; i++)
      if(ebeam.nfFill[i])
      {
        double * mom = cyclic_array_access(i, &dipole_RW);
        MPI_Bcast(mom, dipole_RW.m, MPI_DOUBLE, branks[i], MPI_COMM_WORLD);
        for(j = 1; j < track.Nmlt; j++)
          cyclic_array_set((j*ring.Nharm + i), mom, &dipole_RW);
      }
  }  
  
  /* Initialize model for FBI interaction */
  fbii_model_init(&fbii_model);
  fbii_model.factBE = sqrt(2*(bunchModel.sgm_xx*bunchModel.sgm_xx-bunchModel.sgm_zz*bunchModel.sgm_zz))/FKILO;  
  
  if(bunch.kb_out == 1)  /* Flag for output */ 
  {
    /* Open output files */
    char filename_bstats[FILENAME_MAX] = "";
    snprintf(filename_bstats, FILENAME_MAX, "%s/mean_bunch_%d.dat", track.work_path, kb);  
    bstats_fp = fopen(filename_bstats, "w+");
    if(bstats_fp == NULL)
      ERROR("fopen_stat_bunch0", return);
    worker_weak_stat_output(bstats_fp, bstats, bunch.Ib);
    weak_bunch_reset_statistics(bstats);
    
    if (track.NrevPotentialsOut<=track.NrevTot) {
      char filename_trafo[FILENAME_MAX] = "";
      snprintf(filename_trafo, FILENAME_MAX, "%s/potentials_bunch_%d.dat", track.work_path, kb);  
      trafo_fp = fopen(filename_trafo, "w+");
      if(trafo_fp == NULL)
	ERROR("fopen_potentials_bunch0", return);
      fprintf(trafo_fp, " # potentials and distributions, Ib = %g A ;\n", bunch.Ib);
      fprintf(trafo_fp, " # rev  scan_val         tau         LON_wake_pot       HC_pot        rf_pot      active_cav      ideal_HC        bunch_shape");
      if(SelfFieldModel.PlaneV > 0)
	fprintf(trafo_fp, "     VER_wake_pot     VER_dipole_mom");
      if(SelfFieldModel.PlaneH > 0)
	fprintf(trafo_fp, "     HOR_wake_pot     HOR_dipole_mom");
      fprintf(trafo_fp, "\n");
    }
  }
  
  /* Counting turns for statistics */
  int long Nstat = 0;
  int m = 0;
  int iseed = bunchModel0.iseed[LON] + bunch.kb;
  
  /* ----------------------------------------------------------------------------------------- */
  /* Tracking for NrevTot turns */
  for(rev = 0; rev < track.NrevTot; rev++)
  {      
    /* Resonator selffield transformation includes effect of harmonic cavity and RW */
    /* PlaneL etc give information if a resonator in this plane is given AND if plane is tracked */    
    
    if(SelfFieldModel.PlaneL + SelfFieldModel.PlaneV + SelfFieldModel.PlaneH + lr_res_sizetot > 0)
      transform_weak_bunch_selffield(&bunch, SelfFieldModel, &moments_history, &ring, rev, trafo_fp, scan_val, kb, &ebeam, phasor_end, phasor_end_HOR, phasor_end_VER, fnp_ring, fnp_HOR, fnp_VER);

    /* RF feedback: calculate and apply turn by turn correction to RF parameters */
    if (ring.has_rf_feedback) {
      double vmbar, phibar;
      rffb_calc_mean_voltage_phase(ring.rf_feedback,&vmbar,&phibar);
      vmbar = vmbar*ring.E0*FKILO;
      double tmp_numrtor = ring.rf_feedback->vrf_design*sin(ring.rf_feedback->phi0_design)-vmbar*sin(phibar);
      ring.phai0 = atan2(tmp_numrtor,ring.rf_feedback->vrf_design*cos(ring.rf_feedback->phi0_design)-vmbar*cos(phibar));
      //ring.phai0 = atan2(tmp_numrtor*ring.wrf,ring.rf_feedback->vrf_design*ring.wrf*cos(ring.rf_feedback->phi0_design)-vmbar*ring.longrange_resonators[0].wr*cos(phibar));
      ring.Vrf0 = tmp_numrtor/sin(ring.phai0);
    }
    
    /* Optics transformaiton including quantum excitation & radiation damping and active or passive HC*/
    /* Long. plane is always tracked, hor. and vert. has to be enabled */
    if(transform_weak_bunch_optic(&bunch, &bunchModel, iseed+rev, kb, &ring) < 0)   
      fprintf(stderr, "Error ! transform_optic_LON failed\n");
    
    /* FBII bunch TODO */
    /* FBII Communication TODO */    
    /* FBII train TODO (different approach) */    
    /* FBII Communication TODO */
    
    /* Update the statistics every turn */
    weak_bunch_calc_statistics(&bunch, &track);
    
    if(bunch.kb_out == 1)
    {
      Nstat++;    
      if (Nstat % track.NrevMon == 0 || (rev+1) % track.NrevScan == 0)
      { /* Outputfile ampinv and bunch stats */
      weak_bunch_update_statistics(bstats, Nstat);
      bunch_stats_fprintf(bstats_fp, rev, bstats, scan_val);
      weak_bunch_reset_statistics(bstats);
      Nstat = 0;
      }
    }
    
    /* update cyclic array with actual statistics */
    if(track.EnableRW_long > 0)
    {
      double * mom = cyclic_array_access(-rev*ring.Nharm + kb, &dipole_RW);
      mom[0] = bstats->pos.x;
      mom[1] = bstats->pos.z;
      unsigned int i;
      for(i = 0; i < ring.Nharm; i++)
        if(ebeam.nfFill[i])
        {
          double * mom = cyclic_array_access(-rev*ring.Nharm + i, &dipole_RW);        
          MPI_Bcast(mom, dipole_RW.m, MPI_DOUBLE, branks[i], MPI_COMM_WORLD);
        }
    }
    
    /* Send statistics to manager, every turn */
    MPI_Send(&(bstats->pos), 3, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(&(bstats->pos_sigma), 3, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(&(bstats->slope), 3, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(&(bstats->slope_sigma), 3, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    
    /* Calc and send more statistics to manager, every NrevMon turn if ampinv is outputed */
    if((rev+1)%(track.NrevMon) == 0)// && track.EnableAmpinv_out)
    { 
      weak_bunch_calc_ampinv(&bunch, &track);
      MPI_Send(&(bstats->ampinv), 3, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
      MPI_Send(&(bstats->ampinv_cm), 3, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    }    
    
    if(track.EnableRW_long)
    {
      /* Now include interbunch (long-range resistive-wall) interactions
       * by receiving the CM information of all filled bunches from the master
       * We shall keep this original structure in view of including the RW effects
       * in future.
       */       
      if(track.TrackPlane[HOR])
      {
        const int bpos = -rev*ring.Nharm + kb;  //Actual position of bunch kb in cyclic array
        if(transform_weak_bunch_RW_longrange_cyclic(m, bpos, &dipole_RW, HOR, &bunchModel, &ebeam, &bunch, &ring, &track) < 0)
          fprintf(stderr, "Error ! transformSbunch_longrange HOR failed\n");
      }
      if(track.TrackPlane[VER])
      {
        const int bpos = -rev*ring.Nharm + kb;  //Actual position of bunch kb in cyclic array
        if(transform_weak_bunch_RW_longrange_cyclic(m, bpos, &dipole_RW, VER, &bunchModel, &ebeam, &bunch, &ring, &track) < 0)
          fprintf(stderr, "Error ! transformSbunch_longrange VER failed\n");
      }      
      m++;
    }
    
    if(ring.longrange_resonators_size[LON] > 0)
    {
      if(!fnp_ring_update(&ring, fnp_ring, &SelfFieldModel, &bunch, kb, LON))
        ERROR("fnp_ring_update", return);
      
      unsigned int i;      
      for(i = 0; i < ring.Nharm; i++)
        if(ebeam.nfFill[i])
          /* if kb = i: bunch_moments are sent, if not recieved and stored fnp_ring */
          MPI_Bcast(&fnp_ring[i*SelfFieldModel.Ncell], SelfFieldModel.Ncell, MPI_DOUBLE, branks[i], MPI_COMM_WORLD);
    }
    if(ring.longrange_resonators_size[HOR] > 0 && track.TrackPlane[HOR])
    {
      if(!fnp_ring_update(&ring, fnp_HOR, &SelfFieldModel, &bunch, kb, HOR))
        ERROR("fnp_ring_update", return);
      
      for(i = 0; i < ring.Nharm; i++)
        if(ebeam.nfFill[i])
          /* if kb = i: bunch_moments are sent, if not recieved and stored fnp_ring */
          MPI_Bcast(&fnp_HOR[i*SelfFieldModel.Ncell], SelfFieldModel.Ncell, MPI_DOUBLE, branks[i], MPI_COMM_WORLD);
    }
    if(ring.longrange_resonators_size[VER] > 0 && track.TrackPlane[VER])
    {
      if(!fnp_ring_update(&ring, fnp_VER, &SelfFieldModel, &bunch, kb, VER))
        ERROR("fnp_ring_update", return);
      
      for(i = 0; i < ring.Nharm; i++)
        if(ebeam.nfFill[i])
          /* if kb = i: bunch_moments are sent, if not recieved and stored fnp_ring */
          MPI_Bcast(&fnp_VER[i*SelfFieldModel.Ncell], SelfFieldModel.Ncell, MPI_DOUBLE, branks[i], MPI_COMM_WORLD);
    }
    
    /* Update scan value if rev % NrevScan */
    if((rev+1) % track.NrevScan == 0)
    {
      tracking_scan_step_worker(rev, &scan_val, &track, &ring, &ebeam, &bunch);
      m = 0;
    }
  }
  /*  end tracking ------------------------------------------------------------------------------------------------- */
  /* Update statistics with this new distribution and send to manager */
  if(track.scan == 0)
  {
    weak_bunch_calc_statistics(&bunch, &track);
    
    MPI_Send(&(bstats->pos), 3, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(&(bstats->pos_sigma), 3, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(&(bstats->slope), 3, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(&(bstats->slope_sigma), 3, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD); 
  }
  
  if(bunch.kb_out == 1)
  {
  /* Writeout distribution of first bunch */
    weak_writeout_bunch_distribution(&bunch, bunchModel, LON, bunchModel.nGen[LON]);
    if(track.TrackPlane[HOR])
      weak_writeout_bunch_distribution(&bunch, bunchModel, HOR, bunchModel.nGen[HOR]);
    if(track.TrackPlane[VER]) 
      weak_writeout_bunch_distribution(&bunch, bunchModel, VER, bunchModel.nGen[VER]);

    if(track.EnableRW_long)
    {
      int i, j;
      char filename_hist_dipole[FILENAME_MAX] = "";
      snprintf(filename_hist_dipole, FILENAME_MAX, "%s/dipoles_history.dat", track.work_path);
      FILE * hist_dipole_fp = fopen(filename_hist_dipole, "w+");
      if(hist_dipole_fp == NULL)
        ERROR("fopen_stat_end", return);
      fprintf(hist_dipole_fp, " # history of dipole moments at rev = %i; Ib = %g A ;\n", rev, bunch.Ib);
      fprintf(hist_dipole_fp, " # Nhist  kb   rev   dipole_HOR    dipole_VER\n");
      for(j = 0; j < track.Nmlt; j++)
      {
        for(i = 0; i < ring.Nharm; i++)
        {
          double * mom = cyclic_array_get((j*ring.Nharm + i), &dipole_RW);
          fprintf(hist_dipole_fp, " %d    %d    %d           %e    %e \n", j*ring.Nharm + i, i, j, mom[0], mom[1]);
        }
        
      }
      fclose(hist_dipole_fp);
    }
    fclose(bstats_fp);
    if (trafo_fp!=NULL) fclose(trafo_fp); 
  }
  
  if(ebeam.nfFill[kb] && (bunch.N_trash_low != 0 || bunch.N_trash_high != 0))
  {
    printf("\n WARNING - Binning for selffields bunch #%d:\n %d particles below first bin,    %d particles above last bin\n", kb, bunch.N_trash_low, bunch.N_trash_high);
  }
  
  selffield_model_destroy(&SelfFieldModel);
  weak_bunch_destroy(&bunch);
  if(ring.longrange_resonators_size[LON] > 0)
  {    
    free(fnp_ring);
    fnp_ring = NULL;
  }
  if(ring.longrange_resonators_size[HOR] > 0 && track.TrackPlane[HOR])
  {    
    free(fnp_HOR);
    fnp_ring = NULL;
  }
  if(ring.longrange_resonators_size[VER] > 0 && track.TrackPlane[VER])
  {    
    free(fnp_VER);
    fnp_ring = NULL;
  }
  
  if(track.EnableRW_long > 0)
    cyclic_array_destroy(&dipole_RW);
  
}

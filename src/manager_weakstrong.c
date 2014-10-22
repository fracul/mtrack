#include <stdio.h>
#include "confmpi.h"
#include "types.h"
#include "input.h"
#include "fbi.h"

void manager_weakstrong(const ring_t ring, const tracking_t track,
                        e_beam_t ebeam,
                        const bunch_macroparticle_model_t macrop_model)
{
  bunch_CM_history_strong_t CMhist;
  ion_beam_t ionbeam;
  
  MPI_Status status;
  double xxCM, xpCM, zzCM, zpCM;
  
  bunch_strong_t bunch0; /* TODO init bunch0 */
  bunch_strong_t bunches[500]; /* TODO dynamic allocation */
  
  printf("Starting job \"%s\"\n\n", track.jobtitle);

  /*
   * Generate electron CM distribution
   */
  
  bunch_CM_history_strong_init(&CMhist, (int)(track.NrevTot / track.NrevMon * ebeam.Nbunch));
  /* Distribution parameters in standard case */
  /* used only in weak model
  const double maxTau   = 7.0;
  const double dTau     = 0.1;
  */
  /* Distribution parameters in FBII case */
  /* TODO move to input
  const double maxTauFBI = 7.0;
  const double dTauFBI   = 2*maxTauFBI/((double) (2*ebeam.distrib.Nslc+1)); 
  */
  
  printf("Generating initial distribution...\n");
  
  /* generate distribution for .... */
  int writeout = 1; /* Only the master will writeout distribution files */
  if(strong_generate_eCMs(macrop_model, &ebeam, &bunch0, writeout) < 0)
  {
    fprintf(stderr, "ERROR: generate_eCMs failed\n");
    return;
  }
  
  /* generate general distribution for multibunch tracking */
  if(generate_initial_general_distribution2(bunches, &bunch0,
                                            macrop_model, ring,
                                            ebeam) < 0)
  {
    fprintf(stderr, "Error ! generate_initial_general_distribution2 failed\n");
    return;
  }
  
  /*
   * Multi-bunch tracking
   */
  
  const int slices = 2 * ebeam.distrib.Nslc + 1;
  
  unsigned int kb;
  for(kb = 0; kb < ring.Nharm; kb++)
  {
    /* WARNING: assuming all bunches use bunch0's distribution */
    if(ebeam.nfFill[kb])
      MPI_Send(bunch0.npop, slices, MPI_INT, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
  }
  
  /* update and sort cm data TODO? */
  
  printf("Starting multi bunch tracking...\n");
  
  /* Multi Bunch tracking */
  long int rev;
  for(rev = 0; rev < track.NrevTot; rev++)
  {
    if(rev%100 == 0)
    {
      if(rev+99 >= track.NrevTot)
        printf("Transformation for turn %3ld to %3ld", rev, track.NrevTot-1);
      else
        printf("Transformation for turn %3ld to %3ld", rev, rev+99);
      fflush(stdout);
    }
    
    /* Send filled bunches, worker will start one-turn transformation */
    int ib = 0; /* Counts the number of active (filled) bunches */
    for(kb = 0; kb < ring.Nharm; kb++)
      if(ebeam.nfFill[kb])
      {
        MPI_Send((void *) &(ebeam.Nbunch), 1, MPI_INT, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
        /*MPI_Send(&(rev), 1, MPI_INT, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);*/
        MPI_Send(&kb, 1, MPI_INT, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
        MPI_Send(&ib, 1, MPI_INT, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
        MPI_Send((void *) ring.Ibunch, 500, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
        MPI_Send(bunches[kb].xtau, slices, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
        MPI_Send(bunches[kb].xeps, slices, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
        MPI_Send(bunches[kb].xx, slices, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
        MPI_Send(bunches[kb].xp, slices, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
        MPI_Send(bunches[kb].zz, slices, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
        MPI_Send(bunches[kb].zp, slices, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
        ib++;
      }
    
    /* Recieve result of one-turn transformation (short-range) */
    ib = 0; /* Counts the number of active (filled) bunches */
    for(kb = 0; kb < ring.Nharm; kb++)
      if(ebeam.nfFill[kb])
      {
        int kbn, ibn;
        MPI_Recv(&kbn, 1, MPI_INT, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&ibn, 1, MPI_INT, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(bunches[kb].xtau, slices, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status); 
        MPI_Recv(bunches[kb].xeps, slices, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status); 
        MPI_Recv(bunches[kb].xx, slices, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status); 
        MPI_Recv(bunches[kb].xp, slices, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status); 
        MPI_Recv(bunches[kb].zz, slices, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status); 
        MPI_Recv(bunches[kb].zp, slices, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status); 
        MPI_Recv(&xxCM, 1, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&xpCM, 1, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&zzCM, 1, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&zpCM, 1, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
        
        bunch_CM_history_append(&CMhist, kb, xxCM, xpCM, zzCM, zpCM);
        
        if(rev%(track.NrevMon) == 0) /* Every NrevMon turn */
        {
          const int irevmon = rev/track.NrevMon;
          MPI_Recv(&(CMhist.cm[irevmon * ebeam.Nbunch + kb]), 3, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
          MPI_Recv(&(CMhist.rms[irevmon * ebeam.Nbunch + kb]), 3, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
          MPI_Recv(&(CMhist.cm_slope[irevmon * ebeam.Nbunch + kb]), 3, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
          MPI_Recv(&(CMhist.rms_slope[irevmon * ebeam.Nbunch + kb]), 3, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
          MPI_Recv(&(CMhist.ampinv[irevmon * ebeam.Nbunch + kb]), 3, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
          MPI_Recv(&(CMhist.ampinv_cm[irevmon * ebeam.Nbunch + kb]), 3, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
        }
	ib++;
      }
      if(rev%(track.NrevMon) == 0)
        bunch_strong_writeout_ampinv(rev, rev/track.NrevMon, &CMhist, ebeam);

/*** Let all filled bunches now include interbunch (long-range resistive-wall) interactions     ***/
    /*** Sort first the obtained CM data xxCM etc  ***/
    if(sortCMdata(rev, &CMhist) < 0)
      fprintf(stderr, "Error ! sortCMdata failed\n");
    
    ib = 0;  /*** Counts the number of active (filled) bunches ***/
    for(kb = 0; kb < ring.Nharm; kb++)
      if(ebeam.nfFill[kb])
      {
        MPI_Send(&kb, 1, MPI_INT, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
        MPI_Send(&ib, 1, MPI_INT, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
/* Not necessary to get nfFill[] and Ibunch[] again
         pvm_pkdouble(Ibunch, 500, 1);
         pvm_pkint(nfFill, 500, 1);
*/
        MPI_Send(&(CMhist.dataVV), 5000, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
        MPI_Send(&(CMhist.dataVP), 5000, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
        MPI_Send(&(CMhist.dataXX), 5000, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
        MPI_Send(&(CMhist.dataXP), 5000, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
        ib++;     
      }

/*** Wait and receive the results of interbunch (long-range resistive-wall) interactions        ***/
     ib = 0;  /*** Counts the number of active (filled) bunches ***/
    for(kb = 0; kb < ring.Nharm; kb++)
    { 
      if(ebeam.nfFill[kb])
      {
        int kbn, ibn;
        MPI_Recv(&kbn, 1, MPI_INT, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(&ibn, 1, MPI_INT, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
        MPI_Recv(bunches[kb].xtau, slices, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status); 
        MPI_Recv(bunches[kb].xeps, slices, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status); 
        MPI_Recv(bunches[kb].xx, slices, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status); 
        MPI_Recv(bunches[kb].xp, slices, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status); 
        MPI_Recv(bunches[kb].zz, slices, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status); 
        MPI_Recv(bunches[kb].zp, slices, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
        ib++;     
      } 
    }
    
/*   ----------------------------------------------------------------------------------------------*/
/*          Fast Beam-Ion Interaction: ((Introduced in August 2007))                               */
/*   ----------------------------------------------------------------------------------------------*/

    /* Counts the number of active (filled) bunches */
    ib = 0;
    /* Initialise the number of overall macro-ions to zero (for the moment) */
    ionbeam.model.nFio = 0;
    /* Bunch address up to which the total macro-ion distribution is created */
    ionbeam.model.kbFio = 0;
    
    for(kb=0; kb<ring.Nharm; kb++)
    if(ebeam.nfFill[kb])
    {
      int kbn, ibn;
       
      /* Send to worker */
      MPI_Send(&kb, 1, MPI_INT, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
      MPI_Send(&ib, 1, MPI_INT, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
      MPI_Send(&(ionbeam.model.nFio), 1, MPI_INT, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
      MPI_Send(&(ionbeam.model.kbFio), 1, MPI_INT, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
      MPI_Send(ionbeam.NriFio, ionbeam.model.nFio, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
      MPI_Send(ionbeam.yyFio, ionbeam.model.nFio, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
      MPI_Send(ionbeam.ydFio, ionbeam.model.nFio, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
      MPI_Send(ionbeam.xxFio, ionbeam.model.nFio, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);
      MPI_Send(ionbeam.xdFio, ionbeam.model.nFio, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD);

      /* Recieve the results (FBII) ***/
      MPI_Recv(&kbn, 1, MPI_INT, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
      MPI_Recv(&ibn, 1, MPI_INT, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
      MPI_Recv(&(ionbeam.model.nFio), 1, MPI_INT, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
      MPI_Recv(&(ionbeam.model.kbFio), 1, MPI_INT, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
      MPI_Recv(ionbeam.NriFio, ionbeam.model.nFio, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
      MPI_Recv(ionbeam.yyFio, ionbeam.model.nFio, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
      MPI_Recv(ionbeam.ydFio, ionbeam.model.nFio, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
      MPI_Recv(ionbeam.xxFio, ionbeam.model.nFio, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
      MPI_Recv(ionbeam.xdFio, ionbeam.model.nFio, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
      MPI_Recv(bunches[kb].xtau, slices, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
      MPI_Recv(bunches[kb].xeps, slices, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
      MPI_Recv(bunches[kb].xx, slices, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
      MPI_Recv(bunches[kb].xp, slices, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
      MPI_Recv(bunches[kb].zz, slices, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
      MPI_Recv(bunches[kb].zp, slices, MPI_DOUBLE, kb+1, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    
      if(track.EnableFBII
         && writeoutCMImotions(track.work_path, rev, ib, &ionbeam, &ebeam, macrop_model) < 0)
      {
        fprintf(stderr, "Error: writeoutCMImotions failed\n");
      }
      
      ib++;
    }
    
    if(writeoutVCMmotionsPM2(track.work_path, rev, bunches, &ebeam) < 0)
    {
      fprintf(stderr, "Error: writeoutVCMmotionsPM2 failed\n");
    }
     
    /* TODO writeout_rms_Vbeamsize */
    
    if(rev%10 == 9)
      printf("."); fflush(stdout);
    if(rev%100 == 99 || rev == track.NrevTot-1)
      printf(" done\n");
  }
  
  if(track.TrackPlane[HOR])
    bunch_writeout_mean_ampinv_pos(track.NrevTot, track.NrevMon, &CMhist, ebeam, HOR);
  if(track.TrackPlane[VER])
    bunch_writeout_mean_ampinv_pos(track.NrevTot, track.NrevMon, &CMhist, ebeam, VER);
  
  bunch_CM_history_strong_destroy(&CMhist);
  
  
/* TODO DEBUG CPUTIME
*/



}

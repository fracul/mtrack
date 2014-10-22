#include <math.h>
#include <stdio.h>
#include "confmpi.h"
#include "types.h"
#include "fbi.h"
#include "rw.h"

/* Global variables */
extern ring_t ring;
extern tracking_t track;
extern grid_t fbii_grid;

void worker_weakstrong(const ring_t ring, const tracking_t track,
                       e_beam_t ebeam,
                       const bunch_macroparticle_model_t macrop_model)
{
  bunch_strong_t bunch0;
  bunch_strong_t bunch;
  
  bunch_CM_history_strong_t CMhist;
  ion_beam_t ionbeam_mbtrack;
  
  /*
   * Generate electron CM distribution
   */
  
  bunch_CM_history_strong_init(&CMhist, (int)(track.NrevTot / track.NrevMon * ebeam.Nbunch));

  /* Only the master will writeout distribution files */
  int writeout = 0;
  if(strong_generate_eCMs(macrop_model, &ebeam, &bunch0, writeout) < 0)
  {
    fprintf(stderr, "ERROR: generate_eCMs failed\n");
    return;
  }

  int icell;
  const int slices = 2 * ebeam.distrib.Nslc + 1;
  
  /* Ion beam for single bunch transformation */
  ion_beam_local_t ionbeam;
  
  /* Used to define ionbeam.mapslice2i[icell][msx+offst][msz+offst] */
  int ic = 0;         
  int msx, msz;
  
  /*ion_model_t * ionmodel = &(ionbeam_mbtrack.model);*/
  const int offst = fbii_grid.offst;
                                             
  for(icell=0; icell<slices; icell++)
  {
    for(msx=-offst; msx<offst+1; msx++)
    {
      for(msz=-offst; msz<offst+1; msz++)
      {
        ionbeam.mapslice2i[icell][msx+offst][msz+offst] = ic; 
        ic++;
      }   
    }   
  }
  
  /* Model for FBI interaction */
  fbii_model_t fbii_model;
  fbii_model_init(&fbii_model);
  fbii_model.factBE = sqrt(2*(macrop_model.sgm_xx*macrop_model.sgm_xx-macrop_model.sgm_zz*macrop_model.sgm_zz))/FKILO;
  
  MPI_Status status;
  
  /*
   * Single bunch tranformation
   */
  
  /* Recieve particles population from manager */
  MPI_Recv(bunch0.npop, slices, MPI_INT, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
  
  int rev;
  for(rev = 0; rev < track.NrevTot; rev++)
  {
    int kb = 0;
    int ib = 0;
    
    /* Recieve from manager */
    MPI_Recv((void *) &(ebeam.Nbunch), 1, MPI_INT, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    /*MPI_Recv(&rev, 1, MPI_INT, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);*/
    MPI_Recv(&kb, 1, MPI_INT, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(&ib, 1, MPI_INT, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv((void *) ring.Ibunch, 500, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(bunch0.xtau, slices, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(bunch0.xeps, slices, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(bunch0.xx, slices, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(bunch0.xp, slices, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(bunch0.zz, slices, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(bunch0.zp, slices, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    
    /* Update CM */
    bunch = bunch0;
    
    transformSbunch_shortrange(rev, kb, ebeam, macrop_model, &bunch0, &bunch);
    
    /* FBI for this bunch */
    if(track.EnableFBII
       && fbi_weakstrong_bunch(macrop_model, &(ionbeam_mbtrack.model), fbii_model, ebeam, &ionbeam, &bunch, kb) < 0)
    {
      fprintf(stderr, "Error ! FBI_bunch failed\n");
    }
    
    /* Send results to manager */
    MPI_Send(&kb, 1, MPI_INT, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(&ib, 1, MPI_INT, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(bunch.xtau, slices, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(bunch.xeps, slices, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(bunch.xx, slices, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(bunch.xp, slices, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(bunch.zz, slices, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(bunch.zp, slices, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(&(bunch.stats.pos.x), 1, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(&(bunch.stats.slope.x), 1, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(&(bunch.stats.pos.z), 1, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(&(bunch.stats.slope.z), 1, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    
    if(rev%(track.NrevMon) == 0) /* Every NrevMon turn */
    {
      bunch_strong_distribution_statistics(&bunch, slices);
      MPI_Send(&(bunch.stats.pos), 3, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
      MPI_Send(&(bunch.stats.pos_sigma), 3, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
      MPI_Send(&(bunch.stats.slope), 3, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
      MPI_Send(&(bunch.stats.slope_sigma), 3, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
      MPI_Send(&(bunch.stats.ampinv), 3, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
      MPI_Send(&(bunch.stats.ampinv_cm), 3, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    }
    
/*** Now include interbunch (long-range resistive-wall) interactions            ***/
/*** by receiving the CM information of all filled bunches from the master      ***/
/*** We shall keep this original structure in view of including the RW effects  ***/
/*** in future.                                                                 ***/

/* TODO check again necessity for in and jn
    pvm_upkint(&in, 1, 1);
    pvm_upkint(&jn, 1, 1);
*/
    MPI_Recv(&kb, 1, MPI_INT, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(&ib, 1, MPI_INT, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
/* not necessary to re-send
    pvm_upkdouble(Ibunch, 500, 1);
    pvm_upkint(nfFill, 500, 1);
*/
    MPI_Recv(&(CMhist.dataVV), 5000, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(&(CMhist.dataVP), 5000, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(&(CMhist.dataXX), 5000, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(&(CMhist.dataXP), 5000, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    
    if(resortCMdata(rev, &CMhist) < 0)
      fprintf(stderr, "Error ! resortCMdata failed\n");

    if(transformSbunch_optic(rev, kb, ebeam, &bunch) < 0)
      fprintf(stderr, "Error ! transformSbunch_longrange failed\n");

    if(track.EnableRW_long)
    {
      if(track.TrackPlane[HOR]
         && transformSbunch_longrange(rev, kb, HOR, ebeam, macrop_model, &CMhist, &bunch) < 0)
      fprintf(stderr, "Error ! transformSbunch_longrange LON failed\n");

      if(track.TrackPlane[VER]
         && transformSbunch_longrange(rev, kb, VER, ebeam, macrop_model, &CMhist, &bunch) < 0)
      fprintf(stderr, "Error ! transformSbunch_longrange VER failed\n");
    }

    bunch_strong_updateCM(&bunch, 2 * ebeam.distrib.Nslc + 1, HOR);
    bunch_strong_updateCM(&bunch, 2 * ebeam.distrib.Nslc + 1, VER);


    /* Send the results back to the manager program (long) ***/
    MPI_Send(&kb, 1, MPI_INT, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(&ib, 1, MPI_INT, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(bunch.xtau, slices, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(bunch.xeps, slices, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(bunch.xx, slices, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(bunch.xp, slices, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(bunch.zz, slices, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(bunch.zp, slices, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
     
/*   ----------------------------------------------------------------------------------------------*/
/*          Fast Beam-Ion Interaction: ((Introduced in August 2007))                               */
/*   ----------------------------------------------------------------------------------------------*/

    /* Recieve from manager */
    MPI_Recv(&kb, 1, MPI_INT, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(&ib, 1, MPI_INT, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(&(ionbeam_mbtrack.model.nFio), 1, MPI_INT, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(&(ionbeam_mbtrack.model.kbFio), 1, MPI_INT, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(ionbeam_mbtrack.NriFio, ionbeam_mbtrack.model.nFio, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(ionbeam_mbtrack.yyFio, ionbeam_mbtrack.model.nFio, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(ionbeam_mbtrack.ydFio, ionbeam_mbtrack.model.nFio, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(ionbeam_mbtrack.xxFio, ionbeam_mbtrack.model.nFio, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    MPI_Recv(ionbeam_mbtrack.xdFio, ionbeam_mbtrack.model.nFio, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD, &status);
    
#ifdef DEBUG_LOG
FILE * idebug = confmpi_log_fopen("ion");
fprintf(idebug, "\nkb = %d ; ib = %d ; deltT = %g\n", kb, ib, ebeam.distrib.deltT);
fprintf(idebug, "yyi0[0] = %g ; ydi0[0] = %g\n", ionbeam.yyi0[0], ionbeam.ydi0[0]);
fprintf(idebug, "BEFORE: nFio = %d ; \n", ionbeam_mbtrack.model.nFio);
if(ionbeam_mbtrack.model.nFio > 0)
{
  fprintf(idebug, "NriFio[0] = %g\n", ionbeam_mbtrack.NriFio[0]);
  fprintf(idebug, "yyFio[0] = %g ; ydFio[0] = %g\n", ionbeam_mbtrack.yyFio[0], ionbeam_mbtrack.ydFio[0]);
  fprintf(idebug, "xxFio[0] = %g ; xdFio[0] = %g\n", ionbeam_mbtrack.xxFio[0], ionbeam_mbtrack.xdFio[0]);
}
#endif
    if(track.EnableFBII
       && fbi_weakstrong_train(macrop_model, fbii_model, ebeam, &ionbeam, &ionbeam_mbtrack, &bunch, ib, kb) < 0)
    {
      fprintf(stderr, "Error ! FBI_train failed\n");
    }
#ifdef DEBUG_LOG
fprintf(idebug, "AFTER: nFio = %d\n", ionbeam_mbtrack.model.nFio);
if(ionbeam_mbtrack.model.nFio > 0)
{
fprintf(idebug, "yyFio[0] = %g ; ydFio[0] = %g\n", ionbeam_mbtrack.yyFio[0], ionbeam_mbtrack.ydFio[0]);
fprintf(idebug, "xxFio[0] = %g ; xdFio[0] = %g\n", ionbeam_mbtrack.xxFio[0], ionbeam_mbtrack.xdFio[0]);
}
fclose(idebug);
#endif
    

    /* Send the results back to the master program (FBII) ***/
    MPI_Send(&kb, 1, MPI_INT, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(&ib, 1, MPI_INT, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(&(ionbeam_mbtrack.model.nFio), 1, MPI_INT, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(&(ionbeam_mbtrack.model.kbFio), 1, MPI_INT, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(ionbeam_mbtrack.NriFio, ionbeam_mbtrack.model.nFio, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(ionbeam_mbtrack.yyFio, ionbeam_mbtrack.model.nFio, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(ionbeam_mbtrack.ydFio, ionbeam_mbtrack.model.nFio, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(ionbeam_mbtrack.xxFio, ionbeam_mbtrack.model.nFio, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(ionbeam_mbtrack.xdFio, ionbeam_mbtrack.model.nFio, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(bunch.xtau, slices, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(bunch.xeps, slices, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(bunch.xx, slices, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(bunch.xp, slices, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(bunch.zz, slices, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    MPI_Send(bunch.zp, slices, MPI_DOUBLE, MANAGER_RANK, MBTRACK_TAG, MPI_COMM_WORLD);
    
    
  }
  
  bunch_CM_history_strong_destroy(&CMhist);
}

#include "confmpi.h"
#include "bunch.h"

/* Global variables */
extern ring_t ring;
extern tracking_t track;

#include "fbi.h"

/* Global variable */
extern grid_t fbii_grid;

int generate_initial_general_distribution2(bunch_strong_t * bunches,
                                           const bunch_strong_t * bunch0,
                                           const bunch_macroparticle_model_t MPmodel,
                                           const ring_t ring,
                                           const e_beam_t ebeam)
{ 
/*---------------------------------------------------------------------------------------------*/
/*             Preparation of bunch distributions for all bunches in the ring                  */
/*             prior to performing the multibunch tracking                                     */
/*                                                                                             */     
/*      - The single bunch distribution will be identical for all single bunches,              */
/*        and is given by "generate_initial_sbunch_distribution2()"                            */
/*---------------------------------------------------------------------------------------------*/

  unsigned int kb;
  int icell;
  double  xxCMkb, xpCMkb, zzCMkb, zpCMkb;
  int     iseedN = 571030;
  
  for(kb=0; kb < ring.Nharm; kb++)
  {
    bunch_strong_t * bunch = &(bunches[kb]);
    if(ebeam.nfFill[kb])
    {
       xxCMkb = MPmodel.xxCM_offset * FMILLI * c_fnorm(iseedN);
       xpCMkb = MPmodel.xpCM_offset * FMILLI * c_fnorm(iseedN); 
       zzCMkb = MPmodel.zzCM_offset * FMILLI * c_fnorm(iseedN); 
       zpCMkb = MPmodel.zpCM_offset * FMILLI * c_fnorm(iseedN); 
       for(icell=0; icell < 2 * ebeam.distrib.Nslc + 1; icell++)
       {
          bunch->xtau[icell] = bunch0->xtau[icell];
          bunch->xeps[icell] = bunch0->xeps[icell];
          bunch->xx[icell]   = bunch0->xx[icell] + xxCMkb;
          bunch->xp[icell]   = bunch0->xp[icell] + xpCMkb;
          bunch->zz[icell]   = bunch0->zz[icell] + zzCMkb;
          bunch->zp[icell]   = bunch0->zp[icell] + zpCMkb;
       }   
    }
  }

  return 1;
}

/*----------------------------------------------------------------------------*/

int
strong_generate_eCMs(MPmodel, ebeam, bunch0, writeout)
  const bunch_macroparticle_model_t MPmodel;
  e_beam_t * ebeam;
  bunch_strong_t * bunch0;
  int writeout;
{
  int    i, jp, icell, k, l, ip;
  int    msx, msz, nsx, nsz;
  double Tau;
  double xx0[30000], xp0[30000], zz0[30000], zp0[30000], xtau0[30000], xeps0[30000];
/*  int dbg0 = 1;
  int dbg1 = 1;*/
  
  /* Distribution of macro-electrons (weak) */
  bunch_weak_distribution_t * weakdistrib = &(ebeam->w_distrib);
  
  /* Equivalent distribution of cells (string) */
  bunch_strong_distribution_t * distrib = &(ebeam->distrib);
  
/*----------------------------------------------------------------------------*/

  i = generate_initial_sbunch_distribution2(MPmodel.nGen[LON], writeout,
                                            LON, xtau0, xeps0,
                                            ring, MPmodel);
  if(i < 0) return -1;
  i = generate_initial_sbunch_distribution2(MPmodel.nGen[HOR], writeout,
                                            HOR, xx0, xp0,
                                            ring, MPmodel);
  if(i < 0) return -1;
  i = generate_initial_sbunch_distribution2(MPmodel.nGen[VER], writeout,
                                            VER, zz0, zp0,
                                            ring, MPmodel);
  if(i < 0) return -1;

/*------------------------------------------------------------------------------
 *  Define the number of macro-electrons { npop[icell] } and the initial
 * center of mass coordinates
 * { xxCM0[icell], xpCM0[icell], yyCM0[icell], ypCM0[icell] }
 * for e-bunch slices icell's  (icell = 0, ..., 2*Nslc)
 *----------------------------------------------------------------------------*/

  for(icell=0; icell < NCELL_MAX; icell++) {
     bunch0->xx[icell] = 0.0;     /*** xxCM0[icell]: Centre of mass H position of electrons in cell "icell"       ***/
     bunch0->xp[icell] = 0.0;     /*** xpCM0[icell]: Centre of mass H slope of electrons in cell "icell"          ***/
     bunch0->zz[icell] = 0.0;     /*** zzCM0[icell]: Centre of mass V position of electrons in cell "icell"       ***/
     bunch0->zp[icell] = 0.0;     /*** zzCM0[icell]: Centre of mass V slope of electrons in cell "icell"          ***/
     /*** npop[icell] : Total number of macro electrons in cell "icell"              ***/
     bunch0->npop[icell]  = 0;
  }
  
/*------------------------------------------------------------------------------
 * Stacking Np macro electrons into 2*Nslc+1 cells (Nslc = 0,1, ...)
 * and deducing xxCM0[icell], xpCM0[icell], yyCM0[icell] and ypCM0[icell]
 * in each slice:
 *----------------------------------------------------------------------------*/

#ifdef DEBUG_LOG
  FILE * log_distrib = confmpi_log_fopen("distrib");
  if(log_distrib == NULL)
  {
    fprintf(stderr, "Error ! Cannot open distrib log for writing\n");
    return -1;
  }
#endif
  for(jp=0; jp<MPmodel.Np; jp++)
  {
     Tau   =-xtau0[jp]/MPmodel.sgmatau;
     icell = (int) ((weakdistrib->maxTau + Tau)/weakdistrib->dTau);
     if(icell < 1)        icell = 0;
     if(icell > 2*distrib->Nslc+1) icell = 2*distrib->Nslc+1;
#ifdef DEBUG_LOG
    fprintf(log_distrib, "jp = %3d ; xtau0[%3d] = %g ; Tau = %g\n", jp, jp, xtau0[jp], Tau);
    fprintf(log_distrib, "sgmatau = %g\n", MPmodel.sgmatau);
    fprintf(log_distrib, "icell = %4d\n", icell);
#endif
     weakdistrib->mapcell[jp] = icell;
     weakdistrib->mapcell2jp[icell][bunch0->npop[icell]] = jp;
     (bunch0->npop[icell])++;
     bunch0->xx[icell] = bunch0->xx[icell] + xx0[jp];
     bunch0->xp[icell] = bunch0->xp[icell] + xp0[jp];
     bunch0->zz[icell] = bunch0->zz[icell] + zz0[jp];
     bunch0->zp[icell] = bunch0->zp[icell] + zp0[jp];
  }

  for(icell=0; icell<2*distrib->Nslc+1; icell++)
  {
     if(bunch0->npop[icell]) {
        bunch0->xx[icell] = bunch0->xx[icell]/bunch0->npop[icell];
        bunch0->xp[icell] = bunch0->xp[icell]/bunch0->npop[icell];
        bunch0->zz[icell] = bunch0->zz[icell]/bunch0->npop[icell];
        bunch0->zp[icell] = bunch0->zp[icell]/bunch0->npop[icell];
     }	

#ifdef DEBUG_LOG
    fprintf(log_distrib, "icell: %3d   npop: %5d  xxCM0: %g    xpCM0: %g    yyCM0: %g    ypCM0: %g\n", icell, bunch0->npop[icell], bunch0->xx[icell], bunch0->xp[icell], bunch0->zz[icell], bunch0->zp[icell]);
#endif
     
  }
  
#ifdef DEBUG_LOG
  fclose(log_distrib);
#endif
  
/*  
  if(dbg1) if((i = write2outfile(fpo,1," <FBI_bunch> I have terminated defining electron CMs in <<
   generate_eCMs >> ...")) < 0);
*/

  /*----------------------------------------------------------------------------
   * Creation of local transverse grid points around each xxCM0[icell]
   * and zzCM0[icell] and allocating electron density on the grid points: 
   *--------------------------------------------------------------------------*/

  /* Local constants variables, for better performance */
  const int offst = fbii_grid.offst;
  const double agLOCx = fbii_grid.agLOCx;
  const double agLOCz = fbii_grid.agLOCz;
 
  for(icell=0; icell<2*distrib->Nslc+1; icell++) {
     for(k=0; k<50; k++) for(l=0; l<50; l++) distrib->Nofe[icell][k][l] = 0;  
     for(ip=0; ip < bunch0->npop[icell]; ip++) {
        jp = weakdistrib->mapcell2jp[icell][ip];

        nsz = (int)(2.0*(zz0[jp] - bunch0->zz[icell])/agLOCz);
        if(nsz >= 0)    msz = (nsz+1)/2; else msz = (nsz-1)/2; 
        if(msz > offst) msz = offst;   if(msz < -offst) msz = -offst;

       	nsx = (int)(2.0*(xx0[jp] - bunch0->xx[icell])/agLOCx);
        if(nsx >= 0)    msx = (nsx+1)/2; else msx = (nsx-1)/2;
        if(msx > offst) msx = offst;   if(msx < -offst) msx = -offst;

        distrib->Nofe[icell][msx+offst][msz+offst]++;   

        /*** Nofe[icell][msx+offst][msz+offst]: Number of macro electrons on the local grid "(msx,msz)" of the slice "icell"  ***/
        /*** (msx(z) = -offst, -offst+1,...,-1, 0, 1,..., offst-1, offst)                                                    ***/
     }  /*** End of loop over ip ***/
  }   
#ifdef DEBUG_LOG
{
  FILE * fp = confmpi_log_fopen("nofe");
  int ip;
  if(fp != NULL)
  {
    fprintf(fp, "zzCM0[0] = %4.4g ; agLOCz = %4.4g\n", bunch0->zz[0], agLOCz);
    fprintf(fp, "xxCM0[0] = %4.4g ; agLOCx = %4.4g\n", bunch0->xx[0], agLOCx);
    for(ip=0; ip < bunch0->npop[icell]; ip++)
    {
      const int jp = weakdistrib->mapcell2jp[0][ip];
      fprintf(fp, "zz0[%d] = %4.4g ; ", jp, zz0[jp]);
      fprintf(fp, "xx0[%d] = %4.4g\n", jp, xx0[jp]);
    }
    fprintf_matrix(fp, 50, 50, distrib->Nofe[0]);
  }
  fclose(fp);
}
#endif

  return 0;
}

/*----------------------------------------------------------------------------*/

void bunch_strong_updateCM(bunch_strong_t * bunch, int slices, plane_t plane)
{
  double uuCM = 0.0;
  double upCM = 0.0;
  int npopCM = 0;
  int icell;
  switch(plane)
  {
    case LON:
    {
      for(icell = 0; icell < slices; icell++)
      {
        uuCM += bunch->xtau[icell]*((double) bunch->npop[icell]);
        upCM += bunch->xeps[icell]*((double) bunch->npop[icell]);
        npopCM += bunch->npop[icell];
      }   
      bunch->stats.pos.xtau = uuCM/((double) npopCM);
      bunch->stats.slope.xtau = upCM/((double) npopCM);
    };
    break;
    case HOR:
    {
      for(icell = 0; icell < slices; icell++)
      {
        uuCM += bunch->xx[icell]*((double) bunch->npop[icell]);
        upCM += bunch->xp[icell]*((double) bunch->npop[icell]);
        npopCM += bunch->npop[icell];
      }   
      bunch->stats.pos.x = uuCM/((double) npopCM);
      bunch->stats.slope.x = upCM/((double) npopCM);
    };
    break;
    case VER:
    {
      for(icell = 0; icell < slices; icell++)
      {
        uuCM += bunch->zz[icell]*((double) bunch->npop[icell]);
        upCM += bunch->zp[icell]*((double) bunch->npop[icell]);
        npopCM += bunch->npop[icell];
      }   
      bunch->stats.pos.z = uuCM/((double) npopCM);
      bunch->stats.slope.z = upCM/((double) npopCM);
    };
    break;
  }
}

/*----------------------------------------------------------------------------*/



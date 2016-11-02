#include <stdlib.h>
#include <math.h>
#include "def.h"
#include "types.h"
#include "bunch.h"
#include "fbi.h"
#include "statistics.h"
#include "confmpi.h"

/* Global variables */
extern ring_t ring;
extern tracking_t track;

void bunch_CM_history_weak_init(bunch_CM_history_weak_t * bCM, int Nout)
{
  bCM->ampinv = (vector_t *) calloc((Nout + 2), sizeof(vector_t));
  bCM->ampinv_cm = (vector_t *) calloc((Nout + 2), sizeof(vector_t));
  
  bCM->cm = (vector_t *) calloc((Nout + 2), sizeof(vector_t));
  bCM->cm_slope = (vector_t *) calloc((Nout + 2), sizeof(vector_t));
  bCM->rms = (vector_t *) calloc((Nout + 2), sizeof(vector_t));
  bCM->rms_slope = (vector_t *) calloc((Nout + 2), sizeof(vector_t));
}


void bunch_CM_history_strong_init(bunch_CM_history_strong_t * bCM, int Nout)
{
  int i, j;
  
  for(i=0; i<500; i++) for(j=0; j<30; j++)
  {
    bCM->xxhist[i][j] = 0.0;
    bCM->xphist[i][j] = 0.0;
    bCM->zzhist[i][j] = 0.0;
    bCM->zphist[i][j] = 0.0;
  }
  
  for(i=0; i<15000; i++)
  {
    bCM->dataXX[i] = 0.0;
    bCM->dataXP[i] = 0.0;
    bCM->dataVV[i] = 0.0;
    bCM->dataVP[i] = 0.0;
  }
  
  bCM->ampinv = (vector_t *) calloc((Nout + 2), sizeof(vector_t));
  bCM->ampinv_cm = (vector_t *) calloc((Nout + 2), sizeof(vector_t));
  
  bCM->cm = (vector_t *) calloc((Nout + 2), sizeof(vector_t));
  bCM->cm_slope = (vector_t *) calloc((Nout + 2), sizeof(vector_t));
  bCM->rms = (vector_t *) calloc((Nout + 2), sizeof(vector_t));
  bCM->rms_slope = (vector_t *) calloc((Nout + 2), sizeof(vector_t));
}



void bunch_CM_history_weak_destroy(bunch_CM_history_weak_t * bCM)
{
  free(bCM->ampinv);
  bCM->ampinv = NULL;
  free(bCM->ampinv_cm);
  bCM->ampinv_cm = NULL;
  
  free(bCM->cm);
  free(bCM->cm_slope);
  free(bCM->rms);
  free(bCM->rms_slope);
  bCM->cm = NULL;
  bCM->cm_slope = NULL;
  bCM->rms = NULL;
  bCM->rms_slope = NULL;
}


void bunch_CM_history_strong_destroy(bunch_CM_history_strong_t * bCM)
{
  free(bCM->ampinv);
  bCM->ampinv = NULL;
  free(bCM->ampinv_cm);
  bCM->ampinv_cm = NULL;
  
  free(bCM->cm);
  free(bCM->cm_slope);
  free(bCM->rms);
  free(bCM->rms_slope);
  bCM->cm = NULL;
  bCM->cm_slope = NULL;
  bCM->rms = NULL;
  bCM->rms_slope = NULL;
}


void
bunch_CM_history_append(bunch_CM_history_strong_t * bCM,
                        int kb, double xxCM, double xpCM,
                        double zzCM, double zpCM)
{
  int i;

  /***  cf) NmultiT = 3 means that we consider up to 2 previous turns  ***/
  for(i = 0; i < track.NmultiT; i++) {
     int j = track.NmultiT-i-2;
     if(j >= 0) {
        bCM->xxhist[kb][j+1] = bCM->xxhist[kb][j];
        bCM->xphist[kb][j+1] = bCM->xphist[kb][j];
        bCM->zzhist[kb][j+1] = bCM->zzhist[kb][j];
        bCM->zphist[kb][j+1] = bCM->zphist[kb][j];
     } else {
        bCM->xxhist[kb][0]   = xxCM;
        bCM->xphist[kb][0]   = xpCM;
        bCM->zzhist[kb][0]   = zzCM;
        bCM->zphist[kb][0]   = zpCM;
     }
  }
}

int sortCMdata(in, bCM)
  int in;
  bunch_CM_history_strong_t * bCM;
{ 
  int j, k, jconst;
  unsigned int jb;
  /*
  FILE   *fcp;
  int i, jb, j, k, jconst;
  char cc;
  int dbg0 = 0;
  int dbg1 = 1;
  */

  for(j = 0; j < track.NmultiT; j++) {
     jconst = j*ring.Nharm;
     for(jb = 0; jb < ring.Nharm; jb++) {
        k = jb + jconst;
        bCM->dataXX[k] = bCM->xxhist[jb][j];
        bCM->dataXP[k] = bCM->xphist[jb][j];
        bCM->dataVV[k] = bCM->zzhist[jb][j];
        bCM->dataVP[k] = bCM->zphist[jb][j];
     }
     // printf("\n j = %d, jb = %d \n", j, jb);
  }
  
// printf("\n k = %d \n", k);
/* TODO DEBUG
  if(dbg0  && in<2) {
     erase140(dum140A); erase140(dum140B);
     sprintf(dum140A,"sortCMdata.in%1d", in);
     strmrg(path_name, dum140A, dum140B); 
     if((fcp=fopen(dum140B,"w")) == NULL) {
         printf("\n file = %s does not exist!  Execution terminated\n", dum140B);
         exit(-1);
     } 
     for(j=0; j<NmultiT; j++) {
        jconst = j*Nharm;
        for(jb=0; jb<Nharm; jb++) {
           k = jb + jconst;
           erase80(dum140A); 
           // sprintf(dum140A,"xtauCM[%3d][%1d] = %12.5e,   xxCM[%3d][%1d] = %12.5e,   zzCM[%3d][%1d] = %12.5e", 
           // jb,j,xtauCM[jb][j],jb,j,xxCM[jb][j],jb,j,zzCM[jb][j]);
           putc140N(fcp, dum140A);
        }
     }
     fclose(fcp);
  }
*/
 
/* TODO DEBUG
  if(dbg1  && in>=9 && in<13) {
     erase140(dum140A); erase140(dum140B);
     sprintf(dum140A,"sortCMdata.in%1d", in);
     strmrg(path_name, dum140A, dum140B); 
     if((fcp=fopen(dum140B,"w")) == NULL) {
         printf("\n file = %s does not exist!  Execution terminated\n", dum140B);
         exit(-1);
     } 
     for(jb=0; jb<Nharm; jb++) {
        for(j=0; j<NmultiT; j++) {
           erase140(dum140A);  
           sprintf(dum140A,"zCM[%3d,%1d]=%10.3e, ",  jb,j,zzCM[jb][j]);  
           i = 0;
           cc = dum140A[i];
           while( cc != eow ) {
              putc(cc, fcp);
              i++;  cc = dum140A[i];
           }
        }
        cc = crt;
        putc(cc, fcp);   
     }
     fclose(fcp);
  }
 
*/
  return 0;
}

int resortCMdata(int in, bunch_CM_history_strong_t * CMhist)
{ 
  unsigned int jb;
  int j, k;

  for(j = 0; j < track.NmultiT; j++)
  {
    int jconst = j*ring.Nharm;
    for(jb=0; jb<ring.Nharm; jb++)
    {
      k = jb + jconst;
      CMhist->xxhist[jb][j]   = CMhist->dataXX[k];
      CMhist->xphist[jb][j]   = CMhist->dataXP[k];
      CMhist->zzhist[jb][j]   = CMhist->dataVV[k];
      CMhist->zphist[jb][j]   = CMhist->dataVP[k];
    }
  }
  return 1;
}



/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/






int generate_initial_sbunch_distribution2(int nGen, int writeout, plane_t iplane,
                                          double uu0[30000], double up0[30000],
                                          const ring_t ring,
                                          const bunch_macroparticle_model_t MPmodel)
{ 
  int    idummy, jp;
  double xdummy, ydummy, xx0d, xp0d, zz0d, zp0d;
  int    modeHT;
  double fmodeHT;
  /*int amodeHT; 
  double tauhat, amptau;*/
  double psi0k, argA, argB, cosA, sinA, cosB, sinB;
/*  int    dbg0 = 0;
  int    dbg1 = 0;*/
  int    mode_excitation = 0;     /*** Introduced on 17 July 2006 ***/

/*---------------------------------------------------------------------------------------------*/
/*                      Treatment in case there is only one particle:                          */
/*  [NB] This routine was copied and added to sub_sbtranf.c from sub_mbtrack.c by adjusting    */
/*       the array dimensions for xtau0 etc.                                                   */
/*---------------------------------------------------------------------------------------------*/
  if(MPmodel.Np == 1) {
     if(iplane == LON) {
        uu0[0] = MPmodel.xtauCM_offset * FNANO;
        up0[0] = MPmodel.xepsCM_offset;
     }   
     if(iplane == HOR) {
        uu0[0] = MPmodel.xxCM_offset * FMILLI;
        up0[0] = MPmodel.xpCM_offset * FMILLI;
     }   
     if(iplane == VER) {
        uu0[0] = MPmodel.zzCM_offset * FMILLI;
        up0[0] = MPmodel.zpCM_offset * FMILLI;
     }   
     return 1;
  } 
/*---------------------------------------------------------------------------------------------*/
/*  Generation of Bunch Distribution:                                                          */
/*---------------------------------------------------------------------------------------------*/
  if(nGen) {
     /* if nGen is nonzero, calculate the distribution, instead of reading a file */

     if(iplane == LON) {
        for(jp=0; jp < MPmodel.Np; jp++) {  
           uu0[jp] = (MPmodel.xtauCM_offset + MPmodel.sgm_xtau*c_fnorm(MPmodel.iseed[LON])) * FNANO;
           up0[jp] = MPmodel.xepsCM_offset + MPmodel.sgm_xeps*c_fnorm(MPmodel.iseed[LON]);
        }
     }
     if(iplane == HOR) {
        for(jp=0; jp < MPmodel.Np; jp++) {  
           xx0d    = MPmodel.sgm_xx*c_fnorm(MPmodel.iseed[HOR]);
           xp0d    = MPmodel.sgm_xp*c_fnorm(MPmodel.iseed[HOR]);
           uu0[jp] = xx0d * FMILLI;
           up0[jp] = xp0d * FMILLI;
        }
     }
     if(iplane == VER) {
        if(mode_excitation) {
           modeHT  =-1;	/* TODO get from input */
           /*amodeHT = fabs(modeHT); */
           fmodeHT = (double) modeHT; 
           for(jp=0; jp < MPmodel.Np; jp++) {  
              /*tauhat = sqrt(uu0[jp]*uu0[jp] + ring.fc12*up0[jp]*up0[jp])/MPmodel.sgmatau;
              amptau = pow(tauhat, amodeHT);*/
              psi0k  = atan(ring.fc1*up0[jp]/uu0[jp]);
              if(psi0k < 0.0  && up0[jp] > 0.0) psi0k = psi0k + M_PI;
              if(psi0k > 0.0  && up0[jp] < 0.0) psi0k = psi0k + M_PI;
              argA   = ring.wgziV*uu0[jp];
              argB   = fmodeHT*psi0k;
              cosA   = cos(argA);
              sinA   = sin(argA);
              cosB   = cos(argB);
              sinB   = sin(argB);
              uu0[jp] = (MPmodel.zzCM_offset + MPmodel.sgm_zz*(cosA*cosB + sinA*sinB)) * FMILLI;
              up0[jp] = (MPmodel.zpCM_offset - MPmodel.sgm_zz/ring.beta1[VER]*(sinA*cosB - cosA*sinB)) * FMILLI;
           }
/*           if(dbg1) exit(0);*/
        } else {
           for(jp=0; jp < MPmodel.Np; jp++) {  
              zz0d    = MPmodel.sgm_zz*c_fnorm(MPmodel.iseed[VER]);
              zp0d    = MPmodel.sgm_zp*c_fnorm(MPmodel.iseed[VER]);
              uu0[jp] = zz0d * FMILLI;
              up0[jp] = zp0d * FMILLI;	
          }
        }
     }
     
     /* if writeout, write the distribution in a file */
     if(writeout)
     {
       FILE * fp = fopen(MPmodel.gendst_datafile[iplane], "w");
       if(fp == NULL)
       {
         fprintf(stderr, "Error ! Cannot open %s\n", MPmodel.gendst_datafile[iplane]);
         return -1;
       }
       
       if(iplane == LON)
       {
         for(jp = 0; jp < MPmodel.Np; jp++)
         {
           fprintf(fp, "%5d %15.12g %15.12g\n", jp,uu0[jp]/FNANO, up0[jp]);
         }
       }
       else
       {
         for(jp = 0; jp < MPmodel.Np; jp++)
         {
           fprintf(fp, "%5d %15.12g %15.12g\n", jp,uu0[jp]/FMILLI, up0[jp]/FMILLI);
         }
       }
       fclose(fp);
     }

  } else {
     /*** if nGen is zero, read the distribution from an existing file   ***/  
     FILE * fpp = fopen(MPmodel.gendst_datafile[iplane], "r");
     if(fpp == NULL)
     {
       fprintf(stderr, "Error ! Cannot open %s\n", MPmodel.gendst_datafile[iplane]);
       return -1;
     }
     
     for(jp=0; jp < MPmodel.Np; jp++)
     {      
        if(fscanf(fpp,"%d%lf%lf",&idummy,&xdummy,&ydummy) != 3)
          return -1;
/*        if(dbg0) printf("\n %5d   %10.5lf   %10.5lf", idummy,xdummy,ydummy);*/

        if(iplane == LON) {
           uu0[jp] = xdummy * FNANO;
           up0[jp] = ydummy;
        }   
        if(iplane == HOR) {
           uu0[jp] = xdummy * FMILLI;
           up0[jp] = ydummy * FMILLI;
        }   
        if(iplane == VER) {
           uu0[jp] = xdummy * FMILLI;
           up0[jp] = ydummy * FMILLI;
        }   
     }  
     fclose(fpp);

  }   

/*---------------------------------------------------------------------------------------------*/
  return 1;
}


//static int iseed0 = 0;

double c_fnorm(const int iseed)
{ 
  int  i;
  int truncate = 0;
/*  int dbg0 = 0;*/
  double s, u, v;

  if(iseed != iseed0) {
     srand(iseed);
     iseed0 = iseed;
  }

  pointA:;
  s = 0.0;
  
    //printf("\nRandMax = %e\n", RAND_MAX);
  
  for(i=0; i<12; i++) {
     u = ((double) rand() )/ ((double) RAND_MAX);
     v = u - 0.5;
     s = s + v;
/*     if(dbg0) printf("\n i=%2d: u = %10.7lf, v = %10.7lf, s = %10.5lf", i,u,v,s);*/
  }
  if(truncate && fabs(s) > 1.0) goto pointA;   

  return s;
}

double box_muller(const int iseed)
{
  float x1, x2, w, y1;
  static float y2;
  static int use_last = 0;

  if(iseed != iseed0) {
     srand(iseed);
     iseed0 = iseed;
  }

  if (use_last) {
    y1 = y2;
    use_last = 0;
  }
  else {
    do {
      x1 = 2.0*((double) rand())/((double) RAND_MAX)-1.0;
      x2 = 2.0*((double) rand())/((double) RAND_MAX)-1.0;
      w = x1*x1 + x2*x2;
    } while ( w >= 1.0 );
    
    w = sqrt((-2.0*log(w))/w);
    y1 = x1*w;
    y2 = x2*w;
    use_last = 1;
  }

  return y1;
}

void fprintf_matrix(FILE * fp, const unsigned lines, const unsigned columns, int matrix[][50])
{
  unsigned l, c;
  if(fp == NULL) return;
  for(l = 0; l < lines; l++)
  {
    for(c = 0; c < columns; c++)
    {
      if(c%10 == 0)
      {
        if(c != 0)
        {
          fprintf(fp, "\n");
        }
        if(c+10 >= columns)
        {
          fprintf(fp, "Line %u, Column %u to %u\n", l, c, columns-1);
        }
        else
        {
          fprintf(fp, "Line %u, Column %u to %u\n", l, c, c+9);
        }
      }
      fprintf(fp, "%6d  ", matrix[l][c]);
    }
    fprintf(fp, "\n");
  }
}

static inline
double _ampinv(const bunch_strong_t * bunch, const int slices, const plane_t plane)
{
  int islice;
  double ampinv = 0.0;
  int npops = 0;
  switch(plane)
  {
    case LON:
    {
      for(islice = 0; islice < slices; islice++)
      {
        ampinv += bunch->npop[islice] *
                ( ring.gamma1[plane] * pow(bunch->xtau[islice], 2)
                  + 2.0 * ring.alpha1[plane] * bunch->xtau[islice] * bunch->xeps[islice]
                  + ring.beta1[plane] * pow(bunch->xeps[islice], 2) );
        npops += bunch->npop[islice];
      }
      return ampinv/((double) npops);
    } break;
    case HOR:
    {
      for(islice = 0; islice < slices; islice++)
      {
        ampinv += bunch->npop[islice] *
               ( ring.gamma1[plane] * pow(bunch->xx[islice], 2)
                 + 2.0 * ring.alpha1[plane] * bunch->xx[islice] * bunch->xp[islice]
                 + ring.beta1[plane] * pow(bunch->xp[islice], 2) );
        npops += bunch->npop[islice];
      }
      return ampinv/((double) npops);
    } break;
    case VER:
    {
      for(islice = 0; islice < slices; islice++)
      {
        ampinv += bunch->npop[islice] *
                ( ring.gamma1[plane] * pow(bunch->zz[islice], 2)
                  + 2.0 * ring.alpha1[plane] * bunch->zz[islice] * bunch->zp[islice]
                  + ring.beta1[plane] * pow(bunch->zp[islice], 2) );
        npops += bunch->npop[islice];
      }
      return ampinv/((double) npops);
    } break;
  }
  return -1.0;
}


int bunch_strong_distribution_statistics(bunch_strong_t * bunch, int slices)
{
  bunch_stats_t * stats = &(bunch->stats);
  if(track.TrackPlane[LON])
  {
    bunch_strong_updateCM(bunch, slices, LON);
    stats->pos_sigma.xtau = weighted_rms(bunch->xtau, bunch->npop, slices, stats->pos.xtau);
    stats->slope_sigma.xtau = weighted_rms(bunch->xeps, bunch->npop, slices, stats->slope.xtau);
    stats->ampinv.xtau = _ampinv(bunch, slices, LON);
    stats->ampinv_cm.xtau = ring.gamma1[LON] * pow(stats->pos.xtau, 2)
                             + 2.0*ring.alpha1[LON] * stats->pos.xtau * stats->slope.xtau
                             + ring.beta1[LON] * pow(stats->slope.xtau, 2);
  }
  if(track.TrackPlane[HOR])
  {
    bunch_strong_updateCM(bunch, slices, HOR);
    stats->pos_sigma.x = weighted_rms(bunch->xx, bunch->npop, slices, stats->pos.x);
    stats->slope_sigma.x = weighted_rms(bunch->xp, bunch->npop, slices, stats->slope.x);
    stats->ampinv.x = _ampinv(bunch, slices, HOR);
    stats->ampinv_cm.x = ring.gamma1[HOR] * pow(stats->pos.x, 2)
                          + 2.0*ring.alpha1[HOR] * stats->pos.x * stats->slope.x
                          + ring.beta1[HOR] * pow(stats->slope.x, 2);
  }
  if(track.TrackPlane[VER])
  {
    bunch_strong_updateCM(bunch, slices, VER);
    stats->pos_sigma.z = weighted_rms(bunch->zz, bunch->npop, slices, stats->pos.z);
    stats->slope_sigma.z = weighted_rms(bunch->zp, bunch->npop, slices, stats->slope.z);
    stats->ampinv.z = _ampinv(bunch, slices, VER);
    stats->ampinv_cm.z = ring.gamma1[VER] * pow(stats->pos.z, 2)
                          + 2.0*ring.alpha1[VER] * stats->pos.z * stats->slope.z
                          + ring.beta1[VER] * pow(stats->slope.z, 2);
  }
  return 1;
}

void bunch_fprint_ampinv(FILE * fp, const long int rev, const vector_t * ampinv, const vector_t * ampinv_cm, const int Nharm, const plane_t plane, const e_beam_t ebeam)
{
  int jb;
  fprintf(fp,"# turn: %-5ld\n", rev+1);
  fprintf(fp,"# bunch   ampinv_cm   ampinv\n");
  switch(plane)
  {
    case LON:
    {
      for(jb=0; jb<Nharm; jb++)
      if(ebeam.nfFill[jb])
      {
        fprintf(fp," %3d %15.8e %15.8e\n", jb, ampinv_cm[jb].xtau, ampinv[jb].xtau);
      }
    } break;
    case HOR:
    {
      for(jb=0; jb<Nharm; jb++)
      if(ebeam.nfFill[jb])
      {
        fprintf(fp," %3d %15.8e %15.8e\n", jb, ampinv_cm[jb].x, ampinv[jb].x);
      }
    } break;
    case VER:
    {
      for(jb=0; jb<Nharm; jb++)
      if(ebeam.nfFill[jb])
      {
        fprintf(fp," %3d %15.8e %15.8e\n", jb, ampinv_cm[jb].z, ampinv[jb].z);
      }
    } break;
  }
}

/*----------------------------------------------------------------------------*/

int bunch_weak_writeout_ampinv(const long int rev, const long int irevmon,
                          const bunch_CM_history_weak_t * CMhist,
                          const e_beam_t ebeam)
{
  int ret = 1;
  if(track.TrackPlane[LON])
  {
    char filename[FILENAME_MAX] = "";
    snprintf(filename, FILENAME_MAX, "%s/ampinv_LON_%ld.dat", track.work_path, rev+1);
    FILE * fp = fopen(filename, "w");
    if(fp == NULL)
    {
      ret = -1;
    }
    else
    {
      bunch_fprint_ampinv(fp, rev, &CMhist->ampinv[irevmon], &CMhist->ampinv_cm[irevmon], ring.Nharm, LON, ebeam);
      fclose(fp);
    }
  }
  if(track.TrackPlane[HOR])
  {
    char filename[FILENAME_MAX] = "";
    snprintf(filename, FILENAME_MAX, "%s/ampinv_HOR_%ld.dat", track.work_path, rev+1);
    FILE * fp = fopen(filename, "w");
    if(fp == NULL)
    {
      ret = -1;
    }
    else
    {
      bunch_fprint_ampinv(fp, rev, &CMhist->ampinv[irevmon], &CMhist->ampinv_cm[irevmon], ring.Nharm, HOR, ebeam);
      fclose(fp);
    }
  }
  if(track.TrackPlane[VER])
  {
    char filename[FILENAME_MAX] = "";
    snprintf(filename, FILENAME_MAX, "%s/ampinv_VER_%ld.dat", track.work_path, rev+1);
    FILE * fp = fopen(filename, "w");
    if(fp == NULL)
    {
      ret = -1;
    }
    else
    {
      bunch_fprint_ampinv(fp, rev, &CMhist->ampinv[irevmon], &CMhist->ampinv_cm[irevmon], ring.Nharm, VER, ebeam);
      fclose(fp);
    }
  }
  
  return ret;
}


/*----------------------------------------------------------------------------*/

int bunch_strong_writeout_ampinv(const long int rev, const long int irevmon,
                          const bunch_CM_history_strong_t * CMhist,
                          const e_beam_t ebeam)
{
  int ret = 1;
  if(track.TrackPlane[LON])
  {
    char filename[FILENAME_MAX] = "";
    snprintf(filename, FILENAME_MAX, "%s/ampinv_LON_%ld.dat", track.work_path, rev+1);
    FILE * fp = fopen(filename, "w");
    if(fp == NULL)
    {
      ret = -1;
    }
    else
    {
      bunch_fprint_ampinv(fp, rev, &CMhist->ampinv[irevmon], &CMhist->ampinv_cm[irevmon], ring.Nharm, LON, ebeam);
      fclose(fp);
    }
  }
  if(track.TrackPlane[HOR])
  {
    char filename[FILENAME_MAX] = "";
    snprintf(filename, FILENAME_MAX, "%s/ampinv_HOR_%ld.dat", track.work_path, rev+1);
    FILE * fp = fopen(filename, "w");
    if(fp == NULL)
    {
      ret = -1;
    }
    else
    {
      bunch_fprint_ampinv(fp, rev, &CMhist->ampinv[irevmon], &CMhist->ampinv_cm[irevmon], ring.Nharm, HOR, ebeam);
      fclose(fp);
    }
  }
  if(track.TrackPlane[VER])
  {
    char filename[FILENAME_MAX] = "";
    snprintf(filename, FILENAME_MAX, "%s/ampinv_VER_%ld.dat", track.work_path, rev+1);
    FILE * fp = fopen(filename, "w");
    if(fp == NULL)
    {
      ret = -1;
    }
    else
    {
      bunch_fprint_ampinv(fp, rev, &CMhist->ampinv[irevmon], &CMhist->ampinv_cm[irevmon], ring.Nharm, VER, ebeam);
      fclose(fp);
    }
  }
  
  return ret;
}

/*----------------------------------------------------------------------------*/



int
bunch_writeout_mean_ampinv_pos(const long int NrevTot, const long int NrevMon,
                               const bunch_CM_history_strong_t * CMhist,
                               const e_beam_t ebeam, const plane_t plane)
{
    char filename[FILENAME_MAX] = "";
    char filename_pos[FILENAME_MAX] = "";
    switch(plane)
    {
      case HOR:
        snprintf(filename, FILENAME_MAX, "%s/ampinv_HOR_mean.dat", track.work_path);
        snprintf(filename_pos, FILENAME_MAX, "%s/cmpos_HOR_mean.dat", track.work_path);
      break;
      case VER:
        snprintf(filename, FILENAME_MAX, "%s/ampinv_VER_mean.dat", track.work_path);
        snprintf(filename_pos, FILENAME_MAX, "%s/cmpos_VER_mean.dat", track.work_path);
      break;
      default:
        return -1;
      break;
    }
    FILE * fp = fopen(filename, "w");
    FILE * fp_pos = fopen(filename_pos, "w");
    if(fp != NULL && fp_pos != NULL)
    {
      long int rev;
      for(rev = 0; rev < NrevTot; rev += NrevMon)
      {
        long int irevmon = rev/NrevMon;
        unsigned int ibunch;
        double sum_pos = 0.0;
        double sum_ampinv = 0.0;
        double sum_ampinv_cm = 0.0;
        int Nbunch = 0;
        if(plane == HOR)
          for(ibunch = 0; ibunch < ring.Nharm; ibunch++)
          if(ebeam.nfFill[ibunch])
          {
            sum_pos += CMhist->cm[irevmon * ebeam.Nbunch + ibunch].x;
            sum_ampinv += CMhist->ampinv[irevmon * ebeam.Nbunch + ibunch].x;
            sum_ampinv_cm += CMhist->ampinv_cm[irevmon * ebeam.Nbunch + ibunch].x;
            Nbunch++;
          }
        if(plane == VER)
          for(ibunch = 0; ibunch < ring.Nharm; ibunch++)
          if(ebeam.nfFill[ibunch])
          {
            sum_pos += CMhist->cm[irevmon * ebeam.Nbunch + ibunch].z;
            sum_ampinv += CMhist->ampinv[irevmon * ebeam.Nbunch + ibunch].z;
            sum_ampinv_cm += CMhist->ampinv_cm[irevmon * ebeam.Nbunch + ibunch].z;
            Nbunch++;
          }
        double fNbunch = (double) Nbunch;
        fprintf(fp," %3ld %15.8e %15.8e\n", rev, sum_ampinv/fNbunch, sum_ampinv_cm/fNbunch);
        fprintf(fp_pos," %3ld %15.8e\n", rev, sum_pos/fNbunch);
      }
      fclose(fp);
      fclose(fp_pos);
      return 1;
    }
    else
    {
      return -1;
    }
}

/*----------------------------------------------------------------------------*/


void
worker_weak_stat_output(FILE * fp, const bunch_stats_t * bstats, const double Ib)
{
    fprintf(fp, " # bunch statistics, averaged over turns, Ib = %g A ;\n", Ib);
    fprintf(fp, " # turn - scan_var - AVE_CM_LON - SIG_CM_LON - AVE_delE_LON - SIG_delE_LON -  AVE_sigmaL_LON - SIG_sigmaL_LON - AVE_sigmaE_LON - SIG_sigmaE_LON... (_HOR) ... (_VER)\n # ");
    fprintf(fp, "beg      -/-    LON    %.4e       -/-        %.4e       -/-        %.4e      -/-        %.4e        -/- ",
           bstats->pos.xtau, bstats->slope.xtau, bstats->pos_sigma.xtau, bstats->slope_sigma.xtau);
    if(track.TrackPlane[HOR])
           fprintf(fp, "    HOR     %.4e       -/-        %.4e       -/-        %.4e      -/-        %.4e        -/- ",
           bstats->pos.x, bstats->slope.x, bstats->pos_sigma.x, bstats->slope_sigma.x);
    if(track.TrackPlane[VER])
           fprintf(fp, "    VER     %.4e       -/-       %.4e       -/-        %.4e      -/-        %.4e       -/- ",
          bstats->pos.z, bstats->slope.z, bstats->pos_sigma.z, bstats->slope_sigma.z);
    fprintf(fp, "\n");
}




int
bunch_stats_fprintf(FILE * fp, long int rev, bunch_stats_t * stats, double scan)
{
  fprintf(fp, "%4ld   %.4e   LON    %.4e   %.4e   %.4e   %.4e   %.4e   %.4e   %.4e   %.4e", rev + 1, scan, stats->sum_pos.xtau, sqrt(stats->sum_pos_sqr.xtau), stats->sum_slope.xtau, sqrt(stats->sum_slope_sqr.xtau), stats->sum_pos_sigma.xtau, sqrt(stats->sum_pos_sigma_sqr.xtau), stats->sum_slope_sigma.xtau, sqrt(stats->sum_slope_sigma_sqr.xtau));
   
  if(track.TrackPlane[HOR])
    fprintf(fp, "  HOR    %.4e   %.4e   %.4e   %.4e   %.4e   %.4e   %.4e   %.4e",
            stats->sum_pos.x, sqrt(stats->sum_pos_sqr.x), stats->sum_slope.x, sqrt(stats->sum_slope_sqr.x), stats->sum_pos_sigma.x, sqrt(stats->sum_pos_sigma_sqr.x), stats->sum_slope_sigma.x, sqrt(stats->sum_slope_sigma_sqr.x));
  if(track.TrackPlane[VER])
    fprintf(fp, "  VER    %.4e   %.4e   %.4e   %.4e   %.4e   %.4e   %.4e   %.4e",
            stats->sum_pos.z, sqrt(stats->sum_pos_sqr.z), stats->sum_slope.z, sqrt(stats->sum_slope_sqr.z), stats->sum_pos_sigma.z, sqrt(stats->sum_pos_sigma_sqr.z), stats->sum_slope_sigma.z, sqrt(stats->sum_slope_sigma_sqr.z));
   fprintf(fp, "\n");
  return 0;
}

/*----------------------------------------------------------------------------*/




bool
fbii_grid_setup_parameters(grid_t * grid,
                           e_beam_t * ebeam,
                           const bunch_macroparticle_model_t MPmodel)
{

  /* Adopt the number taken by RZ (Raubenheimer/Zimmermann), 12 on each side */
  grid->Ngrid     = 5;
  /* Number of division of transverse sigma */
  grid->Ndivs     =  3;
  /* Offset: Note that Ngrid should always be odd */
  if((grid->Ngrid)%2 == 0)
  {
    fprintf(stderr, "Ngrid is not an odd number\n");
    return false;
  }
  grid->offst     = (grid->Ngrid-1)/2;
  /* agLOCz: Local vertical grid size */
  grid->agLOCz    = MPmodel.sgm_zz/((double) grid->Ndivs)/FKILO;
  /* agLOCx: Local horizontal grid size */
  grid->agLOCx    = MPmodel.sgm_xx/((double) grid->Ndivs)/FKILO;
  /* agABSz: Absolute vertical grid size */
  grid->agABSz    = 2.0*grid->agLOCz;
  /* agABSx: Absolute horizontal grid size   */
  grid->agABSx    = 2.0*grid->agLOCx;

  /* Distribution of macro-electrons (weak) */
  //bunch_weak_distribution_t * weakdistrib = &(ebeam->w_distrib);
  
  /* Initialize weak distribution parameters */
  /* for wake potential, it was removed
  ncellmax = ((int) (2.0*weakdistrib->maxTau/weakdistrib->dTau)) + 1;
  */
  
  /* Equivalent distribution of cells (string) */
  bunch_strong_distribution_t * distrib = &(ebeam->distrib);
  
  /* Initialize strong distribution parameters */
  distrib->deltT     = distrib->dTau * MPmodel.sgmatau;
  
  return true;
}

void
grid_generate(grid_t * grid,
              e_beam_t * ebeam,
              bunch_strong_t * sbunch0,
              const weak_bunch_t * wbunch0,
              const bunch_macroparticle_model_t MPmodel)
{
  unsigned int jp;
  int icell, k, l, ip;
  int msx, msz, nsx, nsz;
  double Tau;
  
  /* Distribution of macro-electrons (weak) */
  bunch_weak_distribution_t * weakdistrib = &(ebeam->w_distrib);
  
  /* Equivalent distribution of cells (string) */
  bunch_strong_distribution_t * distrib = &(ebeam->distrib);

/*------------------------------------------------------------------------------
 *  Define the number of macro-electrons { npop[icell] } and the initial
 * center of mass coordinates
 * { xxCM0[icell], xpCM0[icell], yyCM0[icell], ypCM0[icell] }
 * for e-bunch slices icell's  (icell = 0, ..., 2*Nslc)
 *----------------------------------------------------------------------------*/

  for(icell = 0; icell < NCELL_MAX; icell++)
  {
     sbunch0->xx[icell] = 0.0;     /*** xxCM0[icell]: Centre of mass H position of electrons in cell "icell"       ***/
     sbunch0->xp[icell] = 0.0;     /*** xpCM0[icell]: Centre of mass H slope of electrons in cell "icell"          ***/
     sbunch0->zz[icell] = 0.0;     /*** zzCM0[icell]: Centre of mass V position of electrons in cell "icell"       ***/
     sbunch0->zp[icell] = 0.0;     /*** zzCM0[icell]: Centre of mass V slope of electrons in cell "icell"          ***/
     /*** npop[icell] : Total number of macro electrons in cell "icell"              ***/
     sbunch0->npop[icell]  = 0;
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
    WARNING("log_distrib")
  }
#endif
  for(jp = 0; jp < wbunch0->Np; jp++)
  {
    const particle_t * particle = &(wbunch0->particles[jp]);
    Tau = - particle->pos.xtau / MPmodel.sgmatau;
    icell = (int) ((weakdistrib->maxTau + Tau)/weakdistrib->dTau);
    if(icell < 1)        icell = 0;
    if(icell > 2*distrib->Nslc+1) icell = 2*distrib->Nslc+1;
#ifdef DEBUG_LOG
    fprintf(log_distrib, "jp = %3d ; xtau0[%3d] = %g ; Tau = %g\n", jp, jp, particle->pos.xtau, Tau);
    fprintf(log_distrib, "sgmatau = %g\n", MPmodel.sgmatau);
    fprintf(log_distrib, "icell = %4d\n", icell);
#endif
    weakdistrib->mapcell[jp] = icell;
    weakdistrib->mapcell2jp[icell][sbunch0->npop[icell]] = jp;
    (sbunch0->npop[icell])++;
    sbunch0->xx[icell] = sbunch0->xx[icell] + particle->pos.x;
    sbunch0->xp[icell] = sbunch0->xp[icell] + particle->slope.x;
    sbunch0->zz[icell] = sbunch0->zz[icell] + particle->pos.z;
    sbunch0->zp[icell] = sbunch0->zp[icell] + particle->slope.z;
  }

  for(icell = 0; icell < 2*distrib->Nslc+1; icell++)
  {
     if(sbunch0->npop[icell])
     {
        sbunch0->xx[icell] = sbunch0->xx[icell]/sbunch0->npop[icell];
        sbunch0->xp[icell] = sbunch0->xp[icell]/sbunch0->npop[icell];
        sbunch0->zz[icell] = sbunch0->zz[icell]/sbunch0->npop[icell];
        sbunch0->zp[icell] = sbunch0->zp[icell]/sbunch0->npop[icell];
     }	

#ifdef DEBUG_LOG
    fprintf(log_distrib, "icell: %3d   npop: %5d  xxCM0: %g    xpCM0: %g    yyCM0: %g    ypCM0: %g\n", icell, sbunch0->npop[icell], sbunch0->xx[icell], sbunch0->xp[icell], sbunch0->zz[icell], sbunch0->zp[icell]);
#endif
     
  }
  
#ifdef DEBUG_LOG
  fclose(log_distrib);
#endif
  
  /*----------------------------------------------------------------------------
   * Creation of local transverse grid points around each xxCM0[icell]
   * and zzCM0[icell] and allocating electron density on the grid points: 
   *--------------------------------------------------------------------------*/

  /* Local constants variables, for better performance */
  const int offst = grid->offst;
  const double agLOCx = grid->agLOCx;
  const double agLOCz = grid->agLOCz;
 
  for(icell=0; icell<2*distrib->Nslc+1; icell++) {
    for(k=0; k<50; k++) for(l=0; l<50; l++) distrib->Nofe[icell][k][l] = 0;  
    for(ip=0; ip < sbunch0->npop[icell]; ip++)
    {
      jp = weakdistrib->mapcell2jp[icell][ip];
      const particle_t * particle = &(wbunch0->particles[jp]);

      nsz = (int)(2.0*(particle->pos.z - sbunch0->zz[icell])/agLOCz);
      if(nsz >= 0)    msz = (nsz+1)/2; else msz = (nsz-1)/2; 
      if(msz > offst) msz = offst;   if(msz < -offst) msz = -offst;

      nsx = (int)(2.0*(particle->pos.x - sbunch0->xx[icell])/agLOCx);
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
    fprintf(fp, "zzCM0[0] = %4.4g ; agLOCz = %4.4g\n", sbunch0->zz[0], agLOCz);
    fprintf(fp, "xxCM0[0] = %4.4g ; agLOCx = %4.4g\n", sbunch0->xx[0], agLOCx);
    for(ip=0; ip < sbunch0->npop[icell]; ip++)
    {
      const int jp = weakdistrib->mapcell2jp[0][ip];
      const particle_t particle = wbunch0->particles[jp];
      fprintf(fp, "zz0[%d] = %4.4g ; ", jp, particle.pos.z);
      fprintf(fp, "xx0[%d] = %4.4g\n", jp, particle.pos.x);
    }
    fprintf_matrix(fp, 50, 50, distrib->Nofe[0]);
  }
  fclose(fp);
}
#endif
}

int fprint_grid(FILE * fp, const grid_t grid)
{
  if(fp != NULL)
  {
    fprintf(fp, "  ===================================================================\n");
    fprintf(fp, "    Grid:                                                            \n");
    fprintf(fp, "  ===================================================================\n");
    fprintf(fp, "  Ngrid      = %d\n", grid.Ngrid);
    fprintf(fp, "  Ndivs      = %d\n", grid.Ndivs);
    fprintf(fp, "  offst      = %d\n", grid.offst);
    fprintf(fp, "  agLOCx     = %12.8g  ;  agLOCz     = %12.8g\n", grid.agLOCx, grid.agLOCz);
    fprintf(fp, "  agABSx     = %12.8g  ;  agABSz     = %12.8g\n", grid.agABSx, grid.agABSz);
    fprintf(fp, "\n");
    return 0;
  }
  else
  {
    return -1;
  }
}

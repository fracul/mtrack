#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tune.h"
#include "def.h"
#include "types.h"
#include "bunch.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

/********************* Coherent tunes from sbtrack not yet tested **********/
/***************************************************************************/
bool
coherent_tune_init(tuneC_t * tune)
{
  tune->Nt = 9;
  tune->TunePlane_C[LON] = 1;
  tune->TunePlane_C[VER] = 1;
  tune->TunePlane_C[HOR] = 1;
  const unsigned NpeakMax = 10;
 if(tune->Nt > 0)
  {
    int N = pow(2,tune->Nt);
    
    if(tune->TunePlane_C[LON])
    {
      tune->xt_tune_C = (double *) malloc(N * sizeof(double));
      if ( tune->xt_tune_C == NULL)
        return false;
    }
    if(tune->TunePlane_C[HOR])
    {
      tune->xx_tune_C = (double *) malloc(N * sizeof(double));
      if ( tune->xx_tune_C == NULL)
        return false;
    }
    if(tune->TunePlane_C[VER])
    {   
      tune->xz_tune_C = (double *) malloc(N * sizeof(double));
      if ( tune->xz_tune_C == NULL)
        return false;
    }

    tune->UL = 0.0;
    tune->UH = 0.0;
    tune->UV = 0.0;
    
    tune->tunes = (double *) malloc(N * sizeof(double));
    if (tune->tunes == NULL)
      return false;
    tune->pk_numbr = (unsigned *) malloc(N * sizeof(unsigned));
    if (tune->pk_numbr == NULL)
      return false;
    tune->pk_value = (double *) malloc(N * NpeakMax * sizeof(double));
    if (tune->pk_value == NULL)
      return false;
    tune->pk_ampli = (double *) malloc(N * NpeakMax * sizeof(double));
    if (tune->pk_ampli == NULL)
      return false;
    
  }
  return true;
}




int
acquire_tunedata_coherent(tuneC_t * tune, int kn, weak_bunch_t * bunch, plane_t plane, selffield_model_t * SelfFieldModel)
{ 
  double sgm_xtau = SelfFieldModel->sigma_tau;  /**< Initial given bunch length [s] */
/*---------------------------------------------------------------------------------------------*/ 
  double dT = sgm_xtau/3.0; /***  Defines the Time Window [s] to Filter the Signal with frev  ***/
/*---------------------------------------------------------------------------------------------*/

  tune->UL = 0.0;
  tune->UH = 0.0;
  tune->UV = 0.0;
  tune->xt_tune_C[kn] = 0.0;
  tune->xx_tune_C[kn] = 0.0;
  tune->xz_tune_C[kn] = 0.0;  
  unsigned int fnc = 0;
  unsigned int jp;

  // Collects both transverse tune datas if either HOR or VER is tracked
  switch (plane)
  {
    case LON:
    {
      for (jp = 0; jp < bunch->Np; jp++)
      {
        particle_t * particle = &(bunch->particles[jp]);
        tune->UL += particle->pos.xtau;
      }
      tune->xt_tune_C[kn] = tune->UL/bunch->Np;
    }
    
    case HOR:   
    {   
      for (jp = 0; jp < bunch->Np; jp++)
      {
        particle_t * particle = &(bunch->particles[jp]);
        if(0.0 < particle->pos.xtau && particle->pos.xtau < dT) 
        {
          tune->UH  += particle->pos.x;   
          fnc++;
        }
      }
      if(fnc > 0) 
      {
        tune->xx_tune_C[kn]   = tune->UH/fnc;
      }
    }
    
    case VER:   
    {   
      for (jp = 0; jp < bunch->Np; jp++)
        {
        particle_t * particle = &(bunch->particles[jp]);
        if(0.0 < particle->pos.xtau && particle->pos.xtau < dT) 
        {
          tune->UV  += particle->pos.z;  
          fnc++;
        }
      }
      if(fnc > 0) 
      {
        tune->xz_tune_C[kn]   = tune->UV/fnc;
      }
    }
  }

/*---------------------------------------------------------------------------------------------*/
  return 1;
}     


/*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$*/

int get_tunepeaks_coherent(tuneC_t * tune, plane_t pl, FILE * fp)
{ 
//  FILE   *fcp;

  int N = pow(2,tune->Nt);

  int    i, j, k;
  int    jxmax, jzmax, jsmax;
  int    NptsTune;
//  int    dbg0 = 0;
  double mumin, arg, co, sn;
  double cax, cbx, caz, cbz, cas, cbs;
  double cofxmax, cofzmax, cofsmax;
  double eps = 1.0e-10;
  double rjxpk, rjzpk, rjspk;
  double cxcen, cxlow, cxhigh;
  double czcen, czlow, czhigh;
  double cscen, cslow, cshigh;
  double tunxpk, tunzpk, tunspk;
  double fNrevTune;
  
  int gsl_fft_real_radix2_transform();
  double Amax;
  double ndata[2*N], ampli[N];

  fNrevTune = (double) N;
  mumin     = 2 * M_PI/fNrevTune;
  NptsTune  = N/2;
  
  
  switch(pl)
  {
    case LON:
      // TODO: not yet tested
    {
      double * cofs = (double *) calloc(N, sizeof(double));
      double * cfsn = (double *) calloc(N, sizeof(double));
      
      for(j=0; j<=NptsTune; j++) 
      {
        cas = cbs = 0.0;
        for(k=0; k<N; k++) 
        {
          arg = ((double) (j*k))*mumin;
          co  = cos(arg);
          sn  = sin(arg);
          cas = cas + co*tune->xt_tune_C[k];
          cbs = cbs + sn*tune->xt_tune_C[k];
        }
        cofs[j] = sqrt(cas*cas + cbs*cbs)/fNrevTune;
        /***printf("\n cofz[%3d] = %12.5e", j,cofz[j]);***/
      }
      cofs[0] = cofs[0]/2.0;
      cofsmax = 0.0;
      jsmax   = 0;
      for(j=1; j<=NptsTune; j++) 
      {
        if(cofs[j] > cofsmax) 
        {
        jsmax = j;
        cofsmax = cofs[j];
        } 
      }
      /*-----------------------------------------------------------------------*/
      if(jsmax == 0 || jsmax == NptsTune) 
      {
        rjspk = (double) jsmax;
        goto pointAs;
      }
      cscen = cofs[jsmax];
      if(fabs(cscen) < eps) 
      {
        rjspk = 0.0;
        goto pointAs;
      }
      cslow  = cofs[jsmax-1];
      cshigh = cofs[jsmax+1];
      rjspk  = jsmax + (cshigh - cslow)/(2.0*cscen - cslow - cshigh)/2.0;
      
      pointAs:;
      tunspk = rjspk/fNrevTune;
      
      printf("\n tunspk = %10.5lf", tunspk);
      printf("\n");
     
      /*---------------------------------------------------------------------------------*/
      /* Use normalised tune spectra to identify other peaks:                            */
      /*---------------------------------------------------------------------------------*/
      for(j=1; j<=NptsTune; j++) cfsn[j] = cofs[j]/cofsmax;

      identify_tunepeaks(tune, cfsn, pl, N, fp);  

      free(cofs);
      free(cfsn);
      break; 
    }
    
    
    case HOR:
    {  
      for(k=0; k<N; k++) 
      {
        ndata[k] = tune->xx_tune_C[k];
        ndata[k+N] = 0.0;
      }
      gsl_fft_real_radix2_transform(&ndata[0], 1, 2*N);
      Amax = 0.0; /* Normalize tune data */
      for (j = 0; j < N; j++)
      {
        ampli[j] = sqrt(ndata[j] * ndata[j] + ndata[2*N-j] * ndata[2*N-j]);
        if(ampli[j] > Amax) Amax = ampli[j];
      }
      for (j = 0; j < N; j++) ampli[j] /= Amax;
      
      identify_tunepeaks(tune, ampli, pl, N, fp);

      break;
    }
    
    
    case VER:
    {
      for(k=0; k<N; k++) 
      {
        ndata[k] = tune->xz_tune_C[k];
        ndata[k+N] = 0.0;
      }
      gsl_fft_real_radix2_transform(&ndata[0], 1, 2*N);
      Amax = 0; /* Normalize tune data */
      for (j = 0; j < N; j++)
      {
        ampli[j] = sqrt(ndata[j] * ndata[j] + ndata[2*N-j] * ndata[2*N-j]);
        if(ampli[j] > Amax) Amax = ampli[j];
      }
      for (j = 0; j < N; j++) 
      {
        ampli[j] /= Amax;
      //printf("\n  %d  %e  %e", j, ampli[j], (double)j/N); 
      }
      identify_tunepeaks(tune, ampli, pl, N, fp);

      break;
    }
      
      
//       double * cofz = (double *) calloc(N, sizeof(double));
//       double * cfzn = (double *) calloc(N, sizeof(double));
//       for(j=0; j<=NptsTune; j++) 
//       {
//         caz = cbz = 0.0;
//         for(k=0; k<N; k++) 
//         {
//           arg = ((double) (j*k))*mumin;
//           co  = cos(arg);
//           sn  = sin(arg);
//           caz += co*tune->xz_tune_C[k];
//           cbz += sn*tune->xz_tune_C[k];
//         }
//         cofz[j] = sqrt(caz*caz + cbz*cbz)/fNrevTune;
// 
//       }
//       cofz[0] = cofz[0]/2.0;
//       cofzmax = 0.0;
//       jzmax = 0;
//       for(j=1; j<=NptsTune; j++) 
//       {
//         if(cofz[j] > cofzmax) 
//         {
//           jzmax = j;
//           cofzmax = cofz[j];
//         }
//       }
// 
//       if(jzmax == 0 || jzmax == NptsTune) 
//       {
//         rjzpk = (double) jzmax;
//       }
//       else
//       {
//         czcen = cofz[jzmax];
//         if(fabs(czcen) < eps) 
//         {
//           rjzpk = 0.0;
//         }
//         else
//         {
//           czlow  = cofz[jzmax-1];
//           czhigh = cofz[jzmax+1];
//           rjzpk  = jzmax + (czhigh - czlow)/(2.0*czcen - czlow - czhigh)/2.0;
//         }
//       }
// 
//       tunzpk = rjzpk/fNrevTune;
//       
//       printf("\n tunzpk = %10.5lf", tunzpk);
//       printf("\n");
// 
//       for(j=1; j<=NptsTune; j++) 
//       {
//         cfzn[j] = cofz[j]/cofzmax;
//         printf("\n  %d  %e  %e", j, cfzn[j], (double)j/N); 
//       }
//       
//       identify_tunepeaks(tune, cfzn, pl, N, fp);  
//             
//       free(cofz);
//       free(cfzn);
//       
//       break;
//     }
  }


// /*---------------------------------------------------------------------------------*/
//   if(write2file) {
//      erase80(dum140A); erase80(dum140B);
//      if(Ib <  10.0) sprintf(dum140A,"tunespectra_%4.2lfmA", Ib);
//      if(Ib >= 10.0) sprintf(dum140A,"tunespectra_%5.2lfmA", Ib);
//      strmrg(path_name, dum140A, dum140B); 
//      if((fcp=fopen(dum140B,"w")) == NULL) {
//          printf("\n file = %s does not exist!  Execution terminated\n", dum140B);
//          exit(-1);
//      } 
//      for(j=1; j<=NptsTune; j++) {
//         erase80(dum140A); 
//         sprintf(dum140A,"%7.5lf   %6.4lf   %6.4lf   %6.4lf", ((double)j)/fNrevTune,cfxn[j],cfzn[j],cfsn[j]);  
//         putc140N(fcp, dum140A);
//      }
//      fclose(fcp);
//   }
/*---------------------------------------------------------------------------------*/
  return 1;
} 



// Data-array: cofu (Use normalised tune spectra to identify other peaks)

int 
identify_tunepeaks(tuneC_t * tune, double * cofu, const unsigned Pl, const long unsigned NrevTune,  FILE * fp)
{  
  const unsigned NpeakMax = 10;
  
  int  i, j, k, jn;
  int  jpk, Npeak, ipeak, jmx[NpeakMax], pcase[NpeakMax];
  int  NptsTune;
  double Vpeak, fNrevTune;    
  double tune_ampmax;    
  double rjpk, cucen, culow, cuhigh;
  double tuneval[NpeakMax], tuneamp[NpeakMax];

    
/*--------------------------------------------------------------------------------------------*/
  //unsigned int NpeakMax = 10;     /*** Maximum number of peaks to count (in case there are so many)      ***/
  double ThreshPS = 0.2;    /*** Threshold for the identification of peaks                         ***/
/*--------------------------------------------------------------------------------------------*/
  fNrevTune = (double) NrevTune; // NrevTune     : Number of revolutions used to analyse the tunes
  NptsTune  = NrevTune/2;
  tune->pk_numbr[Pl] = 0;
  for(i=0; i<NpeakMax; i++) {
     jmx[i] = 0;
     tune->pk_value[Pl * NpeakMax + i] = tune->pk_ampli[Pl * NpeakMax + i] = 0.0;
  }   
  Npeak = 0;   /*** Number of identified peaks (counts up to NpeakMax)  ***/
  j     = 1;
  /***if(iplane == ver) for(j=0; j<NptsTune-1; j++) printf("\n cofu[%3d] = %12.5e", j,cofu[j]);***/
  while(j < NptsTune-1) 
  {
     ipeak = 0;
     if(cofu[j] > cofu[j-1] + ThreshPS) 
     {
        if(cofu[j] > cofu[j+1] + ThreshPS)  
        {  
           ipeak = 1;   
           jpk = j;  
           jn = j+2; 
           
           //output:
           Vpeak = (double) jpk/fNrevTune;
           pcase[Npeak] = 1;
        } 
        else if(cofu[j] < cofu[j+1] - ThreshPS)  
        {  
           ipeak = 0;             jn = j+1;  
        } 
        else if(fabs(cofu[j] - cofu[j+1]) < ThreshPS) 
        {
           k = j+2;
           while(fabs(cofu[j]-cofu[k]) < ThreshPS  &&  k < NptsTune-1) k++;
           if(cofu[j] > cofu[k]+ThreshPS) 
           {  
              ipeak = 1;   jpk = (j+k)/2;  
              pcase[Npeak] = 2;
              //output:
              Vpeak = (double) jpk/fNrevTune;
           }
           if(cofu[j] < cofu[k]+ThreshPS) ipeak = 0;
           jn = k;
        }
        j = jn;

        if(ipeak) 
        {
           cucen  = cofu[jpk];
           culow  = cofu[jpk-1];
           cuhigh = cofu[jpk+1];
           rjpk   = jpk + (cuhigh - culow)/(2.0*cucen - culow - cuhigh)/2.0;
           Vpeak  = rjpk/fNrevTune;
           tune->pk_value[Pl * NpeakMax + Npeak] = Vpeak;
           tune->pk_ampli[Pl * NpeakMax + Npeak] = cofu[jpk];
           Npeak++;
        }
     } 
     else j++;
     if(Npeak == NpeakMax) 
     {
       j = NptsTune;
       printf("\n  maximum peak number (= %d) of tune determination reached", NpeakMax);
       printf("\n  at tune %5.3lf of %5.3lf.", Vpeak, (NptsTune-1)/fNrevTune);
     }
  }  
  
  int rev = 1;
  //output:
     for(i=0; i<Npeak; i++) {
        fprintf(fp,"\n %d   %d    %d    %d    %e    %e", rev, Pl, i, pcase[i], tune->pk_value[Pl * NpeakMax + i], tune->pk_ampli[Pl * NpeakMax + i]);
     }
     fprintf(fp,"\n");

/*---------------------------------------------------------------------------------*/
  return 1;
}      







void
tune_weak_model_destroy(tuneC_t * tune)
{
    //TODO: NOT YET TESTED 
    
  free(tune->xt_tune_C);
  tune->xt_tune_C = NULL;
  free(tune->xx_tune_C);
  tune->xx_tune_C = NULL;
  free(tune->xz_tune_C);
  tune->xz_tune_C = NULL;

  tune->UL = 0.0;
  tune->UL = 0.0;
  tune->UL = 0.0;
  
  free(tune->tunes);
  tune->tunes = NULL;
  
  free(tune->pk_numbr);
  tune->pk_numbr = NULL;
  
  free(tune->pk_value);
  tune->pk_value = NULL;
  
  free(tune->pk_ampli);
  tune->pk_ampli = NULL;
  
}


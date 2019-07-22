#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <errno.h>
#include "bunch.h"
#include "statistics.h"

bool
weak_bunch_create(weak_bunch_t * bunch, const tracking_t * track, const unsigned Np, const double Ib,
                  const int kb, bool allocate)
{
  if(bunch == NULL)
  {
    errno = EFAULT;
    return false;
  }

  bunch->Ib = Ib;
  bunch->Np = Np;
  bunch->kb = kb;
  bunch->kb_out = 0;     
  bunch->particles = NULL;
  
  bunch->stats.ampinv = vector_new(0.0, 0.0, 0.0);
  bunch->stats.ampinv_cm = vector_new(0.0, 0.0, 0.0);
  
  bunch->stats.pos = vector_new(0.0, 0.0, 0.0);
  bunch->stats.pos_sigma = vector_new(0.0, 0.0, 0.0);
  bunch->stats.slope = vector_new(0.0, 0.0, 0.0);
  bunch->stats.slope_sigma = vector_new(0.0, 0.0, 0.0);
  
  bunch->stats.sum_pos = vector_new(0.0, 0.0, 0.0);
  bunch->stats.sum_pos_sigma = vector_new(0.0, 0.0, 0.0);
  bunch->stats.sum_slope = vector_new(0.0, 0.0, 0.0);
  bunch->stats.sum_slope_sigma = vector_new(0.0, 0.0, 0.0);
  
  bunch->stats.sum_pos_sqr = vector_new(0.0, 0.0, 0.0);
  bunch->stats.sum_pos_sigma_sqr = vector_new(0.0, 0.0, 0.0);
  bunch->stats.sum_slope_sqr = vector_new(0.0, 0.0, 0.0);
  bunch->stats.sum_slope_sigma_sqr = vector_new(0.0, 0.0, 0.0);
  
   int i;
   for(i = 0; i < track->Nbunch_out; i++)  
      if(kb == track->bunch_out[i]) bunch->kb_out = 1;
  
  if(bunch->Np == 0)
  {
    bunch->qp = 0.0;
    return true;
  }
  else
  {
    bunch->qp = bunch->Ib * ring.T0 / Np;
    if(allocate)
    {
      bunch->particles = (particle_t *) malloc(bunch->Np * sizeof(particle_t));
      return bunch->particles != NULL;
    }
    else
    {
      return true;
    }
  }
}

bool
weak_bunch_destroy(weak_bunch_t * bunch)
{
  if(bunch == NULL)
    return false;
  
  /* Free used memory */
  if(bunch->particles != NULL) free(bunch->particles);
  
  /* Reset values */
  bunch->Np = 0;
  bunch->qp = 0.0;
  bunch->kb = 0;
  bunch->particles = NULL;
  
  return true;
}


static void
weak_bunch_normal_distribution(weak_bunch_t * bunch,
                               const plane_t plane,
                               const bunch_macroparticle_model_t bunchModel)
{
  double cbphase = 2*M_PI/ring.Nharm*bunchModel.modeCB[plane]*bunch->kb;
  const double pos_offset = bunchModel.pos_offset.v[plane]*sin(cbphase);
  const double pos_sgm = bunchModel.pos_sgm.v[plane];
  const double slope_offset = bunchModel.slope_offset.v[plane]*cos(cbphase);
  const double slope_sgm = bunchModel.slope_sgm.v[plane];
  const double correlate = bunchModel.correlate.v[plane];
  const int iseed = bunchModel.iseed[plane] + bunch->kb; /* (+ bunch->kb) to have different distr. for every bunch */

  unsigned jp;
  for(jp = 0; jp < bunch->Np; jp++)
  {
    particle_t * particle = &(bunch->particles[jp]);
    double tmp_pos = pos_sgm * c_fnorm(iseed);
    particle->pos.v[plane] = pos_offset + tmp_pos;
    particle->slope.v[plane] = slope_offset + slope_sgm * c_fnorm(iseed) + correlate * tmp_pos;
  }
}


static void
weak_bunch_excitation_distribution(weak_bunch_t * bunch,
                                   const bunch_macroparticle_model_t bunchModel)
{
  const double fmodeHT = (double) bunchModel.modeHT; 
  /*const int amodeHT = fabs(modeHT);
  const double sgmatau = bunchModel.sgm_xtau/FGIGA;
  double tauhat, amptau;*/
  double psi0k;
  double argA, argB, cosA, sinA, cosB, sinB;
  
  weak_bunch_mean_rms(bunch, LON);
  double xtauCM_offset = bunch->stats.pos.v[LON];
  double xtaupCM_offset = bunch->stats.slope.v[LON];
  
  const double zzCM_offset = bunchModel.zzCM_offset;
  const double sgm_zz = bunchModel.sgm_zz;
  const double zpCM_offset = bunchModel.zpCM_offset;

  int iseed = bunchModel.iseed[VER];
  if (iseed!=iseed0) {
     srand(iseed);
     iseed0 = iseed;
  }
  double phi_beta = 0*((double) rand())/((double) RAND_MAX)*2*M_PI;
  
  unsigned jp;
  double norm_fac;
  double norm_constant;
  if (track.EnableIdealHC==1) {
    norm_constant = 1.85407*ring.ac*ring.T0/(M_PI*M_PI)/ring.Nharm/ring.wso/sqrt((ring.m_aHC*ring.m_aHC-1)/6.0);
  }
  for(jp = 0; jp < bunch->Np; jp++)
  {
    particle_t * particle = &(bunch->particles[jp]);
    double slopextau = particle->slope.xtau-xtaupCM_offset;
    double posxtau = particle->pos.xtau-xtauCM_offset;
    double posxtau2 = posxtau*posxtau;
    /*tauhat = sqrt(particle->pos.xtau*particle->pos.xtau
                  + ring.fc12*particle->slope.xtau*particle->slope.xtau) / sgmatau;
    amptau = pow(tauhat, amodeHT);*/

    if (track.EnableIdealHC==1)
      norm_fac = norm_constant/sqrt((posxtau2+sqrt(posxtau2*posxtau2+4*slopextau*slopextau*norm_constant*norm_constant))/2.0);
    else norm_fac = ring.fc1;


    psi0k  = atan(norm_fac*slopextau/(particle->pos.xtau-xtauCM_offset));
    if(psi0k < 0.0  && slopextau > 0.0) psi0k = psi0k + M_PI;
    if(psi0k > 0.0  && slopextau < 0.0) psi0k = psi0k + M_PI;
    argA   = 0 * ring.wgziV * particle->pos.xtau;
    argB   = fmodeHT*psi0k+phi_beta;
    cosA   = cos(argA);
    sinA   = sin(argA);
    cosB   = cos(argB);
    sinB   = sin(argB);
    particle->pos.z = (zzCM_offset + sqrt(2)*sgm_zz*(cosA*cosB + sinA*sinB))/FKILO;
    particle->slope.z = (zpCM_offset - sqrt(2)*sgm_zz/ring.beta1[VER]*(sinA*cosB - cosA*sinB) 
			 + sqrt(2)*sgm_zz*ring.alpha1[VER]/ring.beta1[VER]*(cosA*cosB + sinA*sinB))/FKILO;

    /*fprintf(fp, "%10d  %10.5lf  %10.5lf", jp,zz0[0][jp]*FKILO, zp0[0][jp]*FKILO);*/
  }
}

int
weak_generate_bunch_distribution(weak_bunch_t * bunch,
                                 const bunch_macroparticle_model_t bunchModel,
                                 const int TrackPlane[3])
{
  /* TODO get from input */
  //int mode_excitation = 0;     /*** Introduced on 17 July 2006 ***/
  //int modeHT = -1;
  //IMPLEMENTED AS PART OF bunchModel **** 28 August 2014 ****
  
  if(bunch->Np == 0)
    /* Nothing to do */
    return 1;
    
  else if(bunch->Np == 1)
  {
    double cbphase = 2*M_PI/ring.Nharm*bunchModel.modeCB[LON]*bunch->kb;
    bunch->particles[0].pos.v[LON] = bunchModel.pos_offset.v[LON]*sin(cbphase);
    bunch->particles[0].slope.v[LON] = bunchModel.slope_offset.v[LON]*cos(cbphase); 
    
     if(TrackPlane[HOR])
     {
        cbphase = 2*M_PI/ring.Nharm*bunchModel.modeCB[HOR]*bunch->kb;
        bunch->particles[0].pos.v[HOR] = bunchModel.pos_offset.v[HOR]*sin(cbphase);
        bunch->particles[0].slope.v[HOR] = bunchModel.slope_offset.v[HOR]*cos(cbphase); 
     }
     
    if(TrackPlane[VER])
     {
        cbphase = 2*M_PI/ring.Nharm*bunchModel.modeCB[VER]*bunch->kb;
        bunch->particles[0].pos.v[VER] = bunchModel.pos_offset.v[VER]*sin(cbphase);
        bunch->particles[0].slope.v[VER] = bunchModel.slope_offset.v[VER]*cos(cbphase); 
     }
    
    return 1;
  }
  
  else
  {
    if(bunchModel.nGen[LON]) 
     /* if nGen is nonzero, calculate the distribution, instead of reading a file */
        weak_bunch_normal_distribution(bunch, LON, bunchModel);
    else
    { /*** if nGen is zero, read the distribution from an existing file   ***/
      FILE * fpL = fopen(bunchModel.gendst_datafile[LON], "r");
      if(fpL == NULL)
      {
        fprintf(stderr, "Error ! Cannot open %s\n", bunchModel.gendst_datafile[LON]);
        return -1;
      }
      unsigned int idummy, jp;
      double xdummy, ydummy;
      for(jp=0; jp < bunchModel.Np; jp++)
      {      
        particle_t * particle = &(bunch->particles[jp]);
        if(fscanf(fpL,"%d %lf %lf\n",&idummy,&xdummy,&ydummy) != 3)
          return -1;
        //printf("\n %5d   %e   %e", idummy,xdummy,ydummy);
        particle->pos.v[LON] = xdummy * FNANO;
        particle->slope.v[LON] = ydummy;
      }
      fclose(fpL);
    }  
    
    if(TrackPlane[HOR])
    {
      if(bunchModel.nGen[HOR])
          weak_bunch_normal_distribution(bunch, HOR, bunchModel);
      else
      {
        FILE * fpH = fopen(bunchModel.gendst_datafile[HOR], "r");
        if(fpH == NULL)
        {
          fprintf(stderr, "Error ! Cannot open %s\n", bunchModel.gendst_datafile[HOR]);
          return -1;
        }
        unsigned int idummy, jp;
        double xdummy, ydummy;
        for(jp=0; jp < bunchModel.Np; jp++)
        {      
          particle_t * particle = &(bunch->particles[jp]);
          if(fscanf(fpH,"%d %lf %lf\n",&idummy,&xdummy,&ydummy) != 3)
            return -1;
          particle->pos.v[HOR] = xdummy * FMILLI;
          particle->slope.v[HOR] = ydummy * FMILLI;
        } 
        fclose(fpH);
      } 
    }
      
    if(TrackPlane[VER])
    {
      if(bunchModel.nGen[VER])
      {
        if(bunchModel.mode_excitation)
          weak_bunch_excitation_distribution(bunch, bunchModel);
        else
          weak_bunch_normal_distribution(bunch, VER, bunchModel);
      }
    
      else
      {
        FILE * fpV = fopen(bunchModel.gendst_datafile[VER], "r");
        if(fpV == NULL)
        {
          fprintf(stderr, "Error ! Cannot open %s\n", bunchModel.gendst_datafile[VER]);
          return -1;
        }    
        unsigned int idummy, jp;
        double xdummy, ydummy;
        for(jp=0; jp < bunchModel.Np; jp++)
        {      
          particle_t * particle = &(bunch->particles[jp]);
          if(fscanf(fpV,"%d %lf %lf\n",&idummy,&xdummy,&ydummy) != 3)
            return -1;
          particle->pos.v[VER] = xdummy * FMILLI;
          particle->slope.v[VER] = ydummy * FMILLI;
        }
        fclose(fpV); 
      } 
    }
    return 1;
  }
}


bool
weak_writeout_bunch_distribution(const weak_bunch_t * bunch,
                                 const bunch_macroparticle_model_t bunchModel,
                                 const plane_t iplane, int generated)
{
  unsigned int jp;
  char distfname[FILENAME_MAX];

  if (generated==0) {
    unsigned int ns = 0;
    unsigned int st = 0;
    while (bunchModel.gendst_datafile[iplane][ns]!='\0') {
      if (bunchModel.gendst_datafile[iplane][ns]=='/') st = 1*ns;
      ns += 1;
    }
    strcpy(distfname,track.work_path);
    strcat(distfname,bunchModel.gendst_datafile[iplane]+st);
  }
  else strcpy(distfname,bunchModel.gendst_datafile[iplane]);

  FILE * fp = fopen(distfname, "w");
  if(fp == NULL)
  {
    fprintf(stderr, "Error ! Cannot open %s\n", distfname);
    return false;
  }
     
  double prefix_pos[3] = {FGIGA, FKILO, FKILO};   
  double prefix_slope[3] = {1.0, FKILO, FKILO};
     
  for(jp = 0; jp < bunch->Np; jp++)
  {
    particle_t * particle = &(bunch->particles[jp]);
    fprintf(fp, "%5d %15.12g %15.12g\n", jp, particle->pos.v[iplane]*prefix_pos[iplane], particle->slope.v[iplane]*prefix_slope[iplane]);
  }
  fclose(fp);
  return true;
}


void
weak_bunch_update_current(weak_bunch_t * bunch, double Inew)
{
  bunch->Ib = Inew;
  if(bunch->Np != 0) bunch->qp = bunch->Ib * ring.T0 / bunch->Np;
}

int
weak_bunch_writeout_tbtbbb(const long int NrevTot, const long int NrevMon,
			   const bunch_CM_history_weak_t * CMhist,
			   const e_beam_t ebeam, const plane_t plane, const double * scan_val_hist)
{
  char filename[FILENAME_MAX] = "";
    switch(plane)
    {
      case HOR:
        snprintf(filename, FILENAME_MAX, "%s/bunchbybunch_HOR.dat", track.work_path);        
      break;
      case VER:
        snprintf(filename, FILENAME_MAX, "%s/bunchbybunch_VER.dat", track.work_path);
      break;
      case LON:
        snprintf(filename, FILENAME_MAX, "%s/bunchbybunch_LON.dat", track.work_path);
      break;
      default:
        return -1;
      break;
    }
    FILE * fp = fopen(filename,"w");
    if (fp!=NULL)
    {
      fprintf(fp,"# turn    position of bunch 1,2,3...\n");
      long int rev;
      int m = 0;
      for (rev=NrevMon; rev<NrevTot+1; rev+=NrevMon)
      {
	long int irevmon = rev/NrevMon-1;
	unsigned int ibunch;
	fprintf(fp," %3ld %15.8e", rev, scan_val_hist[m]);
	for (ibunch=0; ibunch<ebeam.Nbunch; ibunch++) {
	  if (plane==HOR) fprintf(fp," %15.8e",CMhist->cm[irevmon*ebeam.Nbunch+ibunch].x);
	  if (plane==VER) fprintf(fp," %15.8e",CMhist->cm[irevmon*ebeam.Nbunch+ibunch].z);
	  if (plane==LON) fprintf(fp," %15.8e",CMhist->cm[irevmon*ebeam.Nbunch+ibunch].xtau);
	}
	fprintf(fp,"\n");
      }
      fclose(fp);
      return 1;
    }
    else {
      return -1;
    }
}

int
weak_bunch_writeout_mean_ampinv(const long int NrevTot, const long int NrevMon,
                               const bunch_CM_history_weak_t * CMhist,
                               const e_beam_t ebeam, const plane_t plane, const double * scan_val_hist)
{
    char filename[FILENAME_MAX] = "";
    switch(plane)
    {
      case HOR:
        snprintf(filename, FILENAME_MAX, "%s/ampinv_HOR_mean.dat", track.work_path);        
      break;
      case VER:
        snprintf(filename, FILENAME_MAX, "%s/ampinv_VER_mean.dat", track.work_path);
      break;
      default:
        return -1;
      break;
    }
    FILE * fp = fopen(filename, "w");
    if(fp != NULL)
    {
      fprintf(fp,"# turn   AVE_ampinv    AVE_ampinv_cm   current \n");
      long int rev;
      int m = 0;
      for(rev = 0; rev < NrevTot; rev += NrevMon)
      {
        long int irevmon = rev/NrevMon;
        unsigned int ibunch;
        double sum_ampinv = 0.0;
        double sum_ampinv_cm = 0.0;
        int Nbunch = 0;
        if(plane == HOR)
          for(ibunch = 0; ibunch < ring.Nharm; ibunch++)
          if(ebeam.nfFill[ibunch])
          {
            sum_ampinv += CMhist->ampinv[irevmon * ebeam.Nbunch + ibunch].x;
            sum_ampinv_cm += CMhist->ampinv_cm[irevmon * ebeam.Nbunch + ibunch].x;
            Nbunch++;
          }
        if(plane == VER)
          for(ibunch = 0; ibunch < ring.Nharm; ibunch++)
          if(ebeam.nfFill[ibunch])
          {
            sum_ampinv += CMhist->ampinv[irevmon * ebeam.Nbunch + ibunch].z;
            sum_ampinv_cm += CMhist->ampinv_cm[irevmon * ebeam.Nbunch + ibunch].z;
            Nbunch++;
          }
        double fNbunch = (double) Nbunch;
        fprintf(fp," %3ld %15.8e %15.8e  %15.8e\n", rev, sum_ampinv/fNbunch, sum_ampinv_cm/fNbunch, scan_val_hist[m]);
        m++;
      }
      fclose(fp);
      return 1;
    }
    else
    {
      return -1;
    }
}


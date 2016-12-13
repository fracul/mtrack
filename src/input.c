#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <math.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include "input.h"
#include "types.h"
#include "def.h"
#include "transform_weak.h"
#include "test.h"

static const char * const fillings_str[] =
{
  /* 0 */ "UNIFORM",
  /* 1 */ "ONETHIRD",
  /* 2 */ "TWOTHIRD",
  /* 3 */ "SINGLE",
  /* 4 */ "SINGLE2",
  /* 5 */ "SINGLE3",
  /* 6 */ "THREEFOURTH",
  /* 7 */ "ONEFOURTH",
  /* 8 */ "SINGLE10",
  /* 9 */ "SINGLE20",
  /* 10 */ "SINGLE30",
  /* 11 */ "EVENHALF",
  /* 12 */ "EVENFOURTH"
};

bool scan_filling(char * str, filling_t * filling)
{
  int i;
  for(i = 0; i < 13; i++)
  {
    if(strcmp(fillings_str[i], str) == 0)
    {
      *filling = i;
      return true;
    }
  }
  return false;
}

static const char * const models_str[] =
{
  /* 0 */ "WEAK",
  /* 1 */ "STRONG"
};

bool scan_model(char * str, model_t * model)
{
  int i;
  for(i = 0; i < 2; i++)
  {
    if(strcmp(models_str[i], str) == 0)
    {
      *model = i;
      return true;
    }
  }
  return false;
}

bool scan_plane(char * str, plane_t * plane)
{
  if(strcmp(str, "LON") == 0)
  {
    *plane = LON;
    return true;
  }
  else  if(strcmp(str, "HOR") == 0)
  {
    *plane = HOR;
    return true;
  }
  else  if(strcmp(str, "VER") == 0)
  {
    *plane = VER;
    return true;
  }
  else
  {
    return false;
  }
}

void c_strng(FILE * fp, char * str, const size_t strlen)
{
  char c0, c1;
  unsigned int  i;
  char blank = 0x20;
  char eow = 0x00;
  char crt = 0x0a;
       
  c0 = getc(fp);
  while(c0 == 'C' || c0 == 'c') {
     do  c1 = getc(fp); while ( c1 != crt ); 
     c0 = getc(fp);
  }   

  for(i=0; i < strlen; i++) str[i] = blank;
  i = 0;
  while( c0 != crt ) {
        str[i] = c0;
        i++;
        c0=getc(fp);
  }
  str[i] = eow;
}

static
bool get_extension(const char * filename, char * extension)
{
  int l = strlen(filename);
  int i;
  for(i = l-1; i >= 0; i--)
  {
    if(filename[i] == '.') break;
  }
  if(i < 0)
  {
    return false;
  }
  else
  {
    int j;
    i++;
    for(j = i; j <= l; j++)
      extension[j-i] = filename[j];
    return true;
  }
}

bool read_input(char filename[FILENAME_MAX], ring_t * ring, tracking_t * track,
                e_beam_t * ebeam,
                bunch_macroparticle_model_t * macrop_model,
                selffield_model_t * SelfFieldModel)
{
  bunch_strong_distribution_t * bunchStrongD = &(ebeam->distrib);
  /* Copy the input filename */
  strncpy(track->input_filename, filename, FILENAME_MAX);
  
  FILE * fp = fopen(filename, "r");
  if(fp == NULL)
  {
    fprintf(stderr, "ERROR: Cannot open file %s for reading.\n", filename);
    return false;
  }
  else
  {
    char extension[FILENAME_MAX] = "";
    if(!get_extension(filename, extension))
    {
      fprintf(stderr, "ERROR: Couldn't get input file extension.\n");
      fclose(fp);
      return false;
    }
    else
    {
      bool r = false;
      
      /* Default values */
      
      bunchStrongD->Nslc = 0;
      ring->fYokoyaV = 1.0; ring->fYokoyaH = 1.0; /* Circular ring */
      ring->QH0 = 18.20; ring->QV0 = 10.30; /* SOLEIL specific val. */
      track->EnableFBII = 0;
      track->EnableQuantum = 1;
      track->EnableAmpinv_out = 0;
      track->EnableRW_short = 0; 
      track->EnableRW_long = 0; 
      track->filling = SINGLE3;
      track->triggerRW = 0;
      track->BunchModel = MODEL_WEAK;
      /* Distribution parameters in standard case */
      /* used in weak model */
      ebeam->w_distrib.maxTau = 7.0;
      ebeam->w_distrib.dTau = 0.3;
      ebeam->distrib.maxTau = 7.0; /* not used ? */
      ebeam->distrib.dTau = 0.3; /* not used ? */
      
      /* Read input file */
      
      if(strcmp(extension, "inp") == 0)
        r = read_input_file(fp, ring, track, macrop_model, bunchStrongD);
      else if(strcmp(extension, "conf") == 0)
        r = read_conf_file(fp, ring, track, macrop_model, bunchStrongD, SelfFieldModel);
      else
        fprintf(stderr, "ERROR: Unknown input file extension.\n");
      
      fclose(fp);
      return r;
    }
  }
}

bool read_input_file(FILE * fp, ring_t * ring, tracking_t * track,
                    bunch_macroparticle_model_t * macrop_model,
                    bunch_strong_distribution_t * bunchStrongD)
{
  char str80[81];
  if(fp == NULL)
  {
      return false;
  }

  /* Beam and Machine Parameters */

  c_strng(fp, str80, 80);
  snprintf(track->jobtitle, 80, "%s", str80);
  c_strng(fp, str80, 80);
  sscanf(str80,"%lf", &(ring->E0));
  c_strng(fp, str80, 80);
  sscanf(str80,"%lf%lf", &(ring->Lc), &(ring->rho0));
  c_strng(fp, str80, 80);
  sscanf(str80,"%lf%lf%lf", &(ring->beta1[HOR]), &(ring->beta1[VER]),
                            &(ring->dispH1));
  c_strng(fp, str80, 80);
  sscanf(str80,"%lf%lf%lf", &(ring->alpha1[HOR]), &(ring->alpha1[VER]),
                           &(ring->disppH1));
  c_strng(fp, str80, 80);
  sscanf(str80,"%lf%lf", &(ring->emittanceH), &(ring->couplbeta));
  c_strng(fp, str80, 80);
  sscanf(str80,"%lf%lf", &(ring->ac),&(ring->Je));
  c_strng(fp, str80, 80);
  sscanf(str80,"%lf%lf", &(ring->Vrf0), &(ring->frf));
  c_strng(fp, str80, 80);
  sscanf(str80,"%lf", &(ring->Iring));
  c_strng(fp, str80, 80);
  sscanf(str80,"%lf%lf", &(ring->ApertureH), &(ring->ApertureV));
  c_strng(fp, str80, 80);
  sscanf(str80,"%lf%lf", &(ring->Gzix), &(ring->Gziz));

  /* Tracking Parameters */

  c_strng(fp, str80, 80);
  sscanf(str80,"%d%d%d", &(track->TrackPlane[0]), &(track->TrackPlane[1]),
                         &(track->TrackPlane[2]));
  c_strng(fp, str80, 80);
  sscanf(str80,"%ld%d%ld", &(track->NrevTot), &(track->Nmlt), &(track->NrevMon));
  track->NmultiT = track->Nmlt+2; 

  /* Beam Distribution Parameters */
  
  c_strng(fp, str80, 80);
  sscanf(str80, "%d", &(macrop_model->Np));
  macrop_model->fNp = (double) macrop_model->Np; 
  c_strng(fp, str80, 80);
  sscanf(str80, "%d%d%d%d%d%d", &(macrop_model->nGen[LON]),
         &(macrop_model->iseed[LON]), &(macrop_model->nGen[HOR]),
         &(macrop_model->iseed[HOR]), &(macrop_model->nGen[VER]),
         &(macrop_model->iseed[VER]));
  /* iseed0 = 0; */
  c_strng(fp, str80, 80);
  sscanf(str80, "%s", macrop_model->gendst_datafile[LON]);
  c_strng(fp, str80, 80);
  sscanf(str80, "%s", macrop_model->gendst_datafile[HOR]);
  c_strng(fp, str80, 80);
  sscanf(str80, "%s", macrop_model->gendst_datafile[VER]);
  c_strng(fp, str80, 80);
  sscanf(str80, "%lf%lf%lf%lf", &(macrop_model->xtauCM_offset),
         &(macrop_model->xepsCM_offset), &(macrop_model->sgm_xtau),
         &(macrop_model->sgm_xeps));
  c_strng(fp, str80, 80);
  sscanf(str80, "%lf%lf%lf%lf", &(macrop_model->xxCM_offset),
         &(macrop_model->xpCM_offset), &(macrop_model->sgm_xx),
         &(macrop_model->sgm_xp));
  c_strng(fp, str80, 80);
  sscanf(str80, "%lf%lf%lf%lf", &(macrop_model->zzCM_offset),
         &(macrop_model->zpCM_offset), &(macrop_model->sgm_zz),
         &(macrop_model->sgm_zp));
  
  return true;
}

bool read_conf_file(FILE * fp, ring_t * ring, tracking_t * track,
                    bunch_macroparticle_model_t * macrop_model,
                    bunch_strong_distribution_t * bunchStrongD,
                    selffield_model_t * SelfFieldModel)
{
  int i;
  char section[32] = "";
  char filling_str[32] = "";
  char model_str[32] = "";
  char plane_str[32] = "";
  if(fp == NULL)
  {
      return false;
  }
  
  /* [tracking] */
  snprintf(section, 32, "tracking");
  /* Print errors */
  config_get_fprint = true;

  if(!config_get_str(fp, section, "jobtitle", track->jobtitle))
    return false;
  
  if(!config_get_str(fp, section, "filling", filling_str)
     || !scan_filling(filling_str, &(track->filling)))
  {
    fprintf(stderr, "ERROR: Unknown filling pattern \"%s\"\n", filling_str);
    return false;
  }
  
  if(!config_get_str(fp, section, "BunchModel", model_str)
     || !scan_model(model_str, &(track->BunchModel)) )
    return false;
                    
  if(!config_get_int(fp, section, "TrackPlane_LON", &(track->TrackPlane[0])))
    return false;
  if(track->TrackPlane[0] == 0)
  {
    fprintf(stderr, "ERROR: LON tracking switched off, has to be on ALWAYS!\n");
    return false;
  }
  if(!config_get_int(fp, section, "TrackPlane_HOR", &(track->TrackPlane[1])))
    return false;
  if(!config_get_int(fp, section, "TrackPlane_VER", &(track->TrackPlane[2])))
    return false;

  if(!config_get_long_int(fp, section, "NrevTot", &(track->NrevTot)))
    return false;
  if(!config_get_long_int(fp, section, "NrevMon", &(track->NrevMon)))
    return false;

  if(!config_get_int(fp, section, "EnableRW_short", &(track->EnableRW_short)))
    return false;
  if(!config_get_int(fp, section, "EnableRW_long", &(track->EnableRW_long)))
    return false;
  if(!config_get_int(fp, section, "EnableFBII", &(track->EnableFBII)))
    return false;
 if(!config_get_int(fp, section, "EnableQuantum", &(track->EnableQuantum)))
    return false; 
  if(!config_get_int(fp, section, "EnableIdealHC", &(track->EnableIdealHC)))
    return false;
  
  config_get_fprint = false;
  if(!config_get_long_int(fp, section, "NrevPotentialsOut", &(track->NrevPotentialsOut)))
    track->NrevPotentialsOut = 0;
  if(!config_get_long_int(fp, section, "EnableRW_short_LON", &(track->EnableRW_short_LON)))
    track->EnableRW_short_LON = 1;
  if(!config_get_int(fp, section, "triggerRW", &(track->triggerRW)))
    track->triggerRW = 0;
  if(!config_get_int(fp, section, "Nmlt", &(track->Nmlt)) && track->EnableRW_long)
    track->Nmlt = 0;
  track->NmultiT = track->Nmlt+2;
  if(!config_get_int(fp, section, "EnableAmpinv_out", &(track->EnableAmpinv_out)))
    track->EnableAmpinv_out = 0;

  int enable_resonators;
  if (!config_get_int(fp, section, "EnableResonators", &(enable_resonators)))
    enable_resonators=1;
  
  if(config_get_int(fp, section, "EnableDiffCurr", &(track->EnableDiffCurr)))
  {
    config_get_fprint = true;
    if(!config_get_double(fp, section, "current_ratio", &(track->current_ratio)))
      return false;
    if(track->current_ratio < 0. || track->current_ratio > 1.) 
    {
      fprintf(stderr, "ERROR: current_ration not between 0 and 1.\n");
      return false;
    }
  }

  /* [ring] - Beam and Machine Parameters */
  snprintf(section, 32, "ring");
  /* Print errors */
  config_get_fprint = true;

  if(!config_get_double(fp, section, "E0", &(ring->E0)))
    return false;
  if(!config_get_double(fp, section, "Lc", &(ring->Lc)))
    return false;
  if(!config_get_double(fp, section, "rho0", &(ring->rho0)))
    return false;
  if(!config_get_double(fp, section, "betaH", &(ring->beta1[HOR])))
    return false;
  if(!config_get_double(fp, section, "betaV", &(ring->beta1[VER])))
    return false;
  if(!config_get_double(fp, section, "dispH", &(ring->dispH1)))
    return false;
  if(!config_get_double(fp, section, "alphaH", &(ring->alpha1[HOR])))
    return false;
  if(!config_get_double(fp, section, "alphaV", &(ring->alpha1[VER])))
    return false;
  if(!config_get_double(fp, section, "disppH", &(ring->disppH1)))
    return false;
  if(!config_get_double(fp, section, "ac", &(ring->ac)))
    return false;
  if(!config_get_double(fp, section, "Je", &(ring->Je)))
    return false;
  if(!config_get_double(fp, section, "Vrf0", &(ring->Vrf0)))
    return false;
  if(!config_get_double(fp, section, "frf", &(ring->frf)))
    return false;
  if(!config_get_double(fp, section, "Iring", &(ring->Iring)))
    return false;
  if(!config_get_double(fp, section, "ApertureH", &(ring->ApertureH))
     | !config_get_double(fp, section, "ApertureV", &(ring->ApertureV)))
    return false;
  if(!config_get_double(fp, section, "QH0", &(ring->QH0))
       | !config_get_double(fp, section, "QV0", &(ring->QV0)))
       return false;
    
  double tmp;
  /* Gzix -> normalized Chromaticity;  Cx not normalized */
  config_get_fprint = false; /* no error printing */
  if(!config_get_int(fp, section, "mIdeal", &(ring->m_aHC)))
    ring->m_aHC = 3;
  if(!config_get_double(fp, section, "Gzix", &(ring->Gzix)))
  {
    config_get_fprint = true; /* activate error printing */
    if(config_get_double(fp, section, "Cx", &tmp))
      ring->Gzix = tmp / ring->QH0;
    else
      return false;
  }
  config_get_fprint = false; /* no error printing */  
  if(!config_get_double(fp, section, "Gziz", &(ring->Gziz)))
  {
    config_get_fprint = true; /* activate error printing */
    if(config_get_double(fp, section, "Cz", &tmp))
      ring->Gziz = tmp / ring->QV0;
    else
      return false;
  }

  /* couplbeta is the emittance-ratio */
  
  if(!config_get_double(fp, section, "emittanceH", &(ring->emittanceH)))
    return false;
  config_get_fprint = false; /* no error printing */  
  if(!config_get_double(fp, section, "couplbeta", &(ring->couplbeta)))
  {
    config_get_fprint = true; /* activate error printing */
    if(config_get_double(fp, section, "emittanceV", &tmp))
      ring->couplbeta = tmp / ring->emittanceH;
    else
      return false;
  }
    
  
  /* Wall resistivity needed if RW enabled */
  if(track->EnableRW_short == 1 || track->EnableRW_long == 1)
  {
    if(!config_get_double(fp, section, "rhorw", &(ring->rhorw)))
      return false;
    if(!config_get_double(fp, section, "beffL1", &(ring->beffL[0])))
      return false;
    if(!config_get_double(fp, section, "beffL2", &(ring->beffL[1])))
      return false;
    if(track->TrackPlane[HOR])
    {
      if(!config_get_double(fp, section, "beffH3", &(ring->beffH[0])))
        return false;
      if(!config_get_double(fp, section, "beffH4", &(ring->beffH[1])))
        return false;
    }
    if(track->TrackPlane[VER])
    {
      if(!config_get_double(fp, section, "beffV3", &(ring->beffV[0])))
        return false;
      if(!config_get_double(fp, section, "beffV4", &(ring->beffV[1])))
        return false;
    }
  }
  else
  {
    ring->rhorw = 0.0;
    ring->beffL[0] = 0.0;
    ring->beffL[1] = 0.0;
    ring->beffH[0] = 0.0;
    ring->beffH[1] = 0.0;
    ring->beffV[0] = 0.0;
    ring->beffV[1] = 0.0;
  }
  
  if(track->EnableIdealHC != 1)
  {
    ring->m_aHC = 1;
  }
  
  /* [ring] Optional parameters */
  config_get_fprint = false; /* no error printing */
  config_get_double(fp, section, "fYokoyaV", &(ring->fYokoyaV));
  config_get_double(fp, section, "fYokoyaH", &(ring->fYokoyaH));
  
  /* Enabeling the option to give U0 and taue instead of Je */
  /* if both are given, Urad and Je are calculated therfrom */
  ring->U0 = -1.0;
  ring->taue = -1.0;
  if(config_get_double(fp, section, "U0", &tmp))
    ring->U0 = tmp * FKILO; /* [keV] -> [eV] */
  if(config_get_double(fp, section, "taue", &tmp))
    ring->taue = tmp * FMILLI; /* [ms] -> [s] */
  
  config_get_fprint = true; /* activate error printing */


  /* [distribution] - Beam Distribution Parameters */
  snprintf(section, 32, "distribution");
  /* Print errors */
  config_get_fprint = true;
  
  if(!config_get_int(fp, section, "Nslc", &(bunchStrongD->Nslc)))
    return false;
  
  if(!config_get_int(fp, section, "Np", &(macrop_model->Np)))
    return false;
  macrop_model->fNp = (double) macrop_model->Np;

  config_get_fprint = false;
  if (!config_get_int(fp,section,"modeHT", &(macrop_model->modeHT)))
    macrop_model->mode_excitation=false;
  else macrop_model->mode_excitation=true;
  config_get_fprint = true;

  if(!config_get_int(fp, section, "nGen_LON", &(macrop_model->nGen[LON]))
     | !config_get_int(fp, section, "nGen_HOR", &(macrop_model->nGen[HOR]))
     | !config_get_int(fp, section, "nGen_VER", &(macrop_model->nGen[VER])))
    return false;

  config_get_fprint = false;
  int tmp_seed=time(NULL); /*If a seed is not given, use current time*/
  if (!config_get_int(fp, section, "iseed_LON", &(macrop_model->iseed[LON])))
    macrop_model->iseed[LON] = tmp_seed;
  if (!config_get_int(fp, section, "iseed_HOR", &(macrop_model->iseed[HOR])))
    macrop_model->iseed[HOR] = tmp_seed;
  if (!config_get_int(fp, section, "iseed_VER", &(macrop_model->iseed[VER])))
    macrop_model->iseed[VER] = tmp_seed;
  config_get_fprint = true;
  /* iseed0 = 0; */
    
  if(!config_get_str(fp, section, "gendst_datafile_LON", macrop_model->gendst_datafile[LON])
     | !config_get_str(fp, section, "gendst_datafile_HOR", macrop_model->gendst_datafile[HOR])
     | !config_get_str(fp, section, "gendst_datafile_VER", macrop_model->gendst_datafile[VER]))
    return false;
  
  if(!config_get_double(fp, section, "xtauCM_offset", &(macrop_model->xtauCM_offset))
     | !config_get_double(fp, section, "xepsCM_offset", &(macrop_model->xepsCM_offset))
     | !config_get_double(fp, section, "sgm_xtau", &(macrop_model->sgm_xtau))
     | !config_get_double(fp, section, "sgm_xeps", &(macrop_model->sgm_xeps)))
    return false;

  if(!config_get_double(fp, section, "xxCM_offset", &(macrop_model->xxCM_offset))
     | !config_get_double(fp, section, "xpCM_offset", &(macrop_model->xpCM_offset))
     | !config_get_double(fp, section, "sgm_xx", &(macrop_model->sgm_xx))
     | !config_get_double(fp, section, "sgm_xp", &(macrop_model->sgm_xp)))
    return false;
  
  if(!config_get_double(fp, section, "zzCM_offset", &(macrop_model->zzCM_offset))
     | !config_get_double(fp, section, "zpCM_offset", &(macrop_model->zpCM_offset))
     | !config_get_double(fp, section, "sgm_zz", &(macrop_model->sgm_zz))
     | !config_get_double(fp, section, "sgm_zp", &(macrop_model->sgm_zp)))
    return false;
  
  /*
   * [resonator_1], [resonator_2], ...
   */
  ring->resonators = (resonator_t *) malloc(1 * sizeof(resonator_t));
  config_get_fprint = true;
  
  /* first resonator */
  i = 1;
  snprintf(section, 32, "resonator_%d", i);
  while(config_have_section(fp, section) && enable_resonators!=0)
  {    
    ring->resonators = (resonator_t *) realloc(ring->resonators, i * sizeof(resonator_t));
    resonator_t * resonator = &(ring->resonators[i-1]);
    
    if(!config_get_double(fp, section, "Rs", &(resonator->Rs)))
      return false;
    
    if(!config_get_double(fp, section, "fres", &tmp))
      return false;
    resonator->wr = (tmp*FGIGA)*2.0*M_PI; /* [GHz] -> [Hz] */
    
    if(!config_get_double(fp, section, "Qfactor", &(resonator->Qfactor)))
      return false;
    
    if(resonator->Qfactor == 0.5)
    fprintf(stderr, "WARNING: resonator %d has Qfactor = 0.5 and will be ignored! \n", i);
    
  if(!config_get_str(fp, section, "plane", plane_str)
     || !scan_plane(plane_str, &(resonator->plane)) )
    return false;

    
    /* next resonator */
    i++;
    snprintf(section, 32, "resonator_%d", i);
  }
  ring->resonators_size = i-1;
  
  /*
   * [harmonic_cavity_1], [harmonic_cavity_2], ...
   */
  ring->longrange_resonators = (LR_resonator_t *) malloc(1 * sizeof(LR_resonator_t));
  i = 1;
   /* first HC */
  snprintf(section, 32, "harmonic_cavity_%d", i);
  while(config_have_section(fp, section))
  {
    ring->longrange_resonators = (LR_resonator_t *) realloc(ring->longrange_resonators, i * sizeof(LR_resonator_t));
    LR_resonator_t * lr_resonator = &(ring->longrange_resonators[i-1]);
    
    lr_resonator->facc = NULL;
    lr_resonator->facs = NULL;

    if(!config_get_int(fp, section, "m", &(lr_resonator->m)))
      return false;

    if(!config_get_double(fp, section, "Rs", &(lr_resonator->Rs)))
      return false;
    
    if(!config_get_double(fp, section, "detune", &(lr_resonator->detune)))
      return false;
    
    if(!config_get_double(fp, section, "Qfactor", &(lr_resonator->Qfactor)))
      return false;
    
    /* Long range resonator always act in longitudinal axis */
    
    /* next resonator */
    i++;
    snprintf(section, 32, "harmonic_cavity_%d", i);
  }
  ring->longrange_resonators_size = i-1;

  /*
   * [rf_feedback]
   */
  ring->rf_feedback = (rf_feedback_t *) malloc(1 * sizeof(rf_feedback_t));
  rf_feedback_t * rffb = ring->rf_feedback;
  snprintf(section, 32, "rf_feedback");
  if (config_have_section(fp, section)) {
    ring->has_rf_feedback = 1;
    if (!config_get_int(fp, section, "resonator", &(rffb->lr_resonator)))
      return false;

    if (!config_get_int(fp, section, "averaging_length", &(rffb->len_average)))
      return false;

  }
  
   /*
   * [active_HC_1], [active_HC_2], ...
   */
  ring->active_HC = (active_HC_t *) malloc(1 * sizeof(active_HC_t));
  config_get_fprint = true;  
  /* first active HC */
  i = 1;
  snprintf(section, 32, "active_HC_%d", i);
  while(config_have_section(fp, section))
  {    
    ring->active_HC = (active_HC_t *) realloc(ring->active_HC, i * sizeof(active_HC_t));
    active_HC_t * aHC = &(ring->active_HC[i-1]);
    
    if(!config_get_double(fp, section, "Vpeak", &(aHC->Vpeak)))
      return false;
    
    if(!config_get_double(fp, section, "nHC", &(aHC->nHC)))
      return false;
    
    if(!config_get_double(fp, section, "phi_aHC", &(aHC->phi_aHC)))
      return false;

    /* next active HC */
    i++;
    snprintf(section, 32, "active_HC_%d", i);
  }
  ring->active_HC_size = i-1;


  /*
   * [selffield]
   */
  snprintf(section, 32, "selffield");
  
  SelfFieldModel->RsisL = 0.0;
  SelfFieldModel->aindL = 0.0;
  SelfFieldModel->RsisV = 0.0;
  SelfFieldModel->aindV = 0.0;
  SelfFieldModel->RsisH = 0.0;
  SelfFieldModel->aindH = 0.0;
  config_get_fprint = false; /* no error printing */
  config_get_double(fp, section, "RsisL", &(SelfFieldModel->RsisL)); 
  if(config_get_double(fp, section, "aindL", &(tmp)))
    SelfFieldModel->aindL = tmp / (2.0 * M_PI * FGIGA);  
  config_get_double(fp, section, "RsisV", &(SelfFieldModel->RsisV)); 
  if(config_get_double(fp, section, "aindV", &(tmp)))
    SelfFieldModel->aindV = tmp / (2.0 * M_PI * FGIGA);
    config_get_double(fp, section, "RsisH", &(SelfFieldModel->RsisH)); 
  if(config_get_double(fp, section, "aindH", &(tmp)))
    SelfFieldModel->aindH = tmp / (2.0 * M_PI * FGIGA);
  
  config_get_fprint = true;  
  if(ring->resonators_size > 0 || track->EnableRW_short == 1 || ring->longrange_resonators_size > 0)
  {
    if(!config_have_section(fp, section))
      return false;
    
    if(!config_get_double(fp, section, "Nsigma", &(SelfFieldModel->Nsigma))
       | !config_get_int(fp, section, "Ncell", &(SelfFieldModel->Ncell)) )
      return false;
  }
  else
  {
    SelfFieldModel->Nsigma = 1.0;
    SelfFieldModel->Ncell = 0;
  }
  
   /*
   * [scan]
   */
  track->scan = 0;
  track->Nscan = 1;
  track->NrevScan = track->NrevTot + 1;
  track->scan_start = 0.0;
  track->scan_step = 0.0;
  snprintf(section, 32, "scan");
  config_get_fprint = true;
  if(config_have_section(fp, section))
  {
    if(!config_get_int(fp, section, "scan_type", &(track->scan)))
      return false;
    
    if(track->scan != 0)
    {       
      if(!config_get_long_int(fp, section, "Nscan", &(track->Nscan)))
        return false;
    
      if(!config_get_long_int(fp, section, "NrevScan", &(track->NrevScan)))
        return false;
    
      if(track->scan == 1)
      {
        if(!config_get_double(fp, section, "scan_start", &(track->scan_start)))
          return false;
      
        if(!config_get_double(fp, section, "scan_step", &(track->scan_step)))
          return false;
      }
      
      if(track->scan == 4)
      {
        if(!config_get_double(fp, section, "scan_start", &(track->scan_start)))
          return false;
      }

      track->NrevTot = track->Nscan * track->NrevScan;
    }
  }
     
  /* At this point everything went well */
  return true;
}

int config_scan_section(const char * line, const char * section_name)
{
  if(line[0] == '#' || line[0] == ';') /* Commented line */
  {
    return -1;
  }
  else
  {
    unsigned int i = 0;
    unsigned int begin = 0;
    unsigned int end = 0;
    
    /* Pass spaces and tabs at the beginning */
    while(line[i] != '\0' && (line[i] == ' ' || line[i] == '\t'))
      i++;
    
    /* Look for opening bracket */
    if(line[i] != '[') return -1;
    
    i++;
    begin = i;
    
    /* Look for end of work */
    while(line[i] != '\0'
          && line[i] != ' ' && line[i] != '\t')
      i++;
    
    /* Look for closing bracket */
    if(line[i-1] != ']') return -1;
    
    end = i-2;
    
    /* Check section name length */
    if(strlen(section_name) != end-begin+1)
    {
      return 0;
    }
    /* Match section name */
    if(strncmp(line + begin, section_name, end-begin+1) == 0)
    {
      return 1;
    }
    else
    {
      return 0;
    }
  }
}

int config_scan(const char * line, const char * variable_name, const char * format, va_list args)
{
  int r = 0;
  int i = 0;
  
  /* Pass spaces and tabs at the beginning */
  while(line[i] != '\0' && (line[i] == ' ' || line[i] == '\t'))
    i++;
  
  if(line[i] == '#' || line[i] == ';') /* Commented line */
  {
    r = 0;
  }
  else
  {
    unsigned int begin = i;
    unsigned int end = 0;
    
    /* Look for end of identifier */
    while(line[i] != '\0'
          && line[i] != ' ' && line[i] != '\t'
          && line[i] != ':' && line[i] != '=')
     i++;
    end = i-1;
   
    /* Pass spaces and tabs after the identifier */
    while(line[i] != '\0' && (line[i] == ' ' || line[i] == '\t'))
      i++;
    
    /* Look for assign operator */
    if(line[i] != ':' && line[i] != '=')
    {
      r = 0;
    }
    else
    {
      i++; /* Get pass assign operator */
      
      /* Pass spaces and tabs after the assign operator */
      while(line[i] != '\0' && (line[i] == ' ' || line[i] == '\t'))
        i++;
      
      /* Check variable length */
      if(strlen(variable_name) != end-begin+1)
      {
        r = 0;
      }
      /* Check variable name */
      else if(strncmp(line + begin, variable_name, end-begin+1) == 0)
      {
        if(strcmp(format, "%s") == 0)
        {
          strcpy(va_arg(args, char *), line + i);
          r = 1;
        }
        else
        {
          r = vsscanf(line + i, format, args);
        }
      }
      else
      {
        r = 0;
      }
    }
  }
  return r;
}

bool config_get_fprint = false; /* Global variable */

bool config_have_section(FILE * fp, const char * section)
{
  int count = 0;
  
  if(fp == NULL)
  {
    return false;
  }
  else
  {
    char line[256];
    fseek(fp, 0, SEEK_SET);
    while(!feof(fp))
    {
      if(fgets(line, 256, fp) == NULL)
      {
        break;
      }
      else
      {
        int len = strlen(line);
        if(line[len-1] == '\n') line[len-1] = '\0';
        
        if(config_scan_section(line, section) > 0)
        {
          count++;
        }
      }
    }
  }
  
  if(count > 1 && config_get_fprint)
    fprintf(stderr, "ERROR: Section [%s] is not unique !\n", section);
  
  return count == 1;
}

int config_get(FILE * fp, const char * section, const char * identifier, const char * format, ...)
{
  int r = 0;
  
  va_list args;
  va_start(args, format);
  
  if(fp == NULL || format == NULL)
  {
    r = -1;
  }
  else
  {
    char line[256];
    bool in_section = false;
    fseek(fp, 0, SEEK_SET);
    while(!feof(fp))
    {
      if(fgets(line, 256, fp) == NULL)
      {
        r = -1;
      }
      else
      {
        int len = strlen(line);
        if(line[len-1] == '\n') line[len-1] = '\0';
        
        int section_result = config_scan_section(line, section);
        if(section_result == 0)
        {
          in_section = false;
        }
        else if(section_result > 0)
        {
          in_section = true;
        }
        else if(in_section)
        {
          int id_result = config_scan(line, identifier, format, args);
          if(id_result > 0)
          {
            r = id_result;
            break;
          }
        }
      }
    }
  }
  
  va_end(args);
  
  if(r <= 0 && config_get_fprint)
    fprintf(stderr, "ERROR: Could not get identifier %s in section [%s]\n", identifier, section);
  return r;
}

inline
bool config_get_int(FILE * fp, const char * section, const char * identifier, int * r)
{
  if(config_get(fp, section, identifier, "%d", r) == 1)
    return true;
  else
    return false;
}

inline
bool config_get_long_int(FILE * fp, const char * section, const char * identifier, long int * r)
{
  if(config_get(fp, section, identifier, "%d", r) == 1)
    return true;
  else
    return false;
}

inline
bool config_get_double(FILE * fp, const char * section, const char * identifier, double * r)
{
  if(config_get(fp, section, identifier, "%lf", r) == 1)
    return true;
  else
    return false;
}

inline
bool config_get_str(FILE * fp, const char * section, const char * identifier,
                    char * r)
{
  if(config_get(fp, section, identifier, "%s", r) == 1)
    return true;
  else
    return false;
}

bool check_parameters(ring_t * ring, tracking_t * track, e_beam_t * ebeam)
{
  /* Warnings (only print a message) */
  test_init(true);
  ASSERT_TRUE(ring->beffL[0] < 1);
  ASSERT_TRUE(ring->beffL[1] < 1);
  ASSERT_TRUE(ring->beffH[0] < 1);
  ASSERT_TRUE(ring->beffH[1] < 1);
  ASSERT_TRUE(ring->beffV[0] < 1);
  ASSERT_TRUE(ring->beffV[1] < 1);
  if(test_results.errors > 0)
  {
    fprintf(stderr, "WARNING: %d possible errors\n", test_results.errors);
  }
  
  /* Errors (print a warning + end as a FAILURE) */
  test_init(true);

//  int tmp;  
//   tmp = (track->NrevTot - 1)/track->NrevMon;
//   ASSERT(tmp <= 100,
//          "Monitor turns limit reached (%d > 100), try increase NrevMon.", tmp);
  
  ASSERT(2*ebeam->distrib.Nslc+1 <= NCELL_MAX,
         "Cell limit reached (%d > %d), try decrease Nslc.",
         2*ebeam->distrib.Nslc+1, NCELL_MAX);
  
  if(test_results.errors > 0)
  {
    fprintf(stderr, "FAILURE: %d input errors\n", test_results.errors);
    return false;
  }
  else
  {
    return true;
  }
}

int fprint_ring(FILE * fp, const ring_t ring)
{
  if(fp == NULL)
  {
    return -1;
  }
  else
  {
    fprintf(fp, "  ===================================================================\n");
    fprintf(fp, "    Beam and Machine Parameters, inputs:                                     \n");
    fprintf(fp, "  ===================================================================\n");
    fprintf(fp, "  E0      = %8.2lf [GeV]\n", ring.E0);
    fprintf(fp, "  Lc      = %8.2lf [m]\n",   ring.Lc);
    fprintf(fp, "  rho0    = %8.2lf [m]\n",   ring.rho0);
    fprintf(fp, "  betaH   = %8.2lf [m],          betaV     = %8.2lf,          dispH    = %8.2lf\n",
            ring.beta1[HOR], ring.beta1[VER], ring.dispH1);
    fprintf(fp, "  alphaH  = %8.2lf [m],          alphaV    = %8.2lf,          disppH   = %8.2lf\n",
            ring.alpha1[HOR], ring.alpha1[VER], ring.disppH1);
    fprintf(fp, "  emitt   = %8.2lf [nm],         coupling  = %9.4lf\n",
            ring.emittanceH, ring.couplbeta);
    fprintf(fp, "  ac      = %9.4e,            Je        = %9.4lf           De*1e3   = %9.4lf\n", 
            ring.ac, ring.Je, ring.De*1.e3);
    fprintf(fp, "  Vrf0    = %8.2lf [MV],         frf       = %9.3lf [MHz]\n",
            ring.Vrf0, ring.frf);
    fprintf(fp, "  Iring   = %8.2lf [mA]\n",
            ring.Iring);
    fprintf(fp, "  AperH   = %8.2lf [mm],         AperV     = %8.2lf [mm]\n",
            ring.ApertureH, ring.ApertureV);
    fprintf(fp, "  Gzix    = %9.4lf,              Gziz      = %9.4lf\n",
            ring.Gzix, ring.Gziz);
    fprintf(fp, "  rhorw   = %.3e [Ohm m]\n", ring.rhorw);
    fprintf(fp, "  beffL(1)   = %.3e [m],        beffL(2)   = %.3e [m]\n", ring.beffL[0], ring.beffL[1]);
    fprintf(fp, "  beffH(3)   = %.3e [m],        beffH(4)   = %.3e [m]\n", ring.beffH[0], ring.beffH[1]);
    fprintf(fp, "  beffV(3)   = %.3e [m],        beffV(4)   = %.3e [m]\n", ring.beffV[0], ring.beffV[1]);
    fprintf(fp, "\n");
    return 0;
  }
}

int fprint_tracking(FILE * fp, const tracking_t track)
{
  if(fp == NULL)
  {
    return -1;
  }
  else
  {
    fprintf(fp, "  ===================================================================\n");
    fprintf(fp, "    Tracking Parameters:                                             \n");
    fprintf(fp, "  ===================================================================\n");
    fprintf(fp, "  TrackPlane[LON] = %1d,         TrackPlane[HOR] = %1d,         TrackPlane[VER] = %1d\n",
           track.TrackPlane[LON], track.TrackPlane[HOR], track.TrackPlane[VER]);
    fprintf(fp, "  NrevTot         = %6ld,    NrevMon         = %6ld\n",
            track.NrevTot, track.NrevMon);
    if(track.scan != 0)
    {
      if(track.scan == 1)
          fprintf(fp, "  CURRENT scan:    Istart       = %8.2lf [mA],    Istep         = %8.2lf [mA]\n",
            track.scan_start, track.scan_step);
      if(track.scan == 2)
          fprintf(fp, "  CHROMA scan:    Gstart         = %9.4lf,             Gstep         = %9.4lf\n",
            track.scan_start, track.scan_step);    
      if(track.scan == 3)
          fprintf(fp, "  QFACTOR scan:\n");    
      if(track.scan == 4)
          fprintf(fp, "  QFACTOR scan:   start      * = %8.2lf  -> down\n",
            track.scan_start);
      fprintf(fp, "  NrevScan         = %6ld,    Nscan            = %6ld\n",
            track.NrevScan, track.Nscan);
    }
    fprintf(fp, "  Quantum excit. & rad. damp. enabled:     %d\n", track.EnableQuantum);
    fprintf(fp, "  Resistive wall self-force enabled:       %d", track.EnableRW_short);
    if (track.EnableRW_short)
    {
      fprintf(fp, ",  s0 = %.2e m (dTau = %d*s0) \n", track.s0, (int)track.mult);
    }
    fprintf(fp, "\n  Resistive wall longrange enabled:       %d, trigger: %d", track.EnableRW_long, track.triggerRW);
    if (track.EnableRW_long)
    {
      fprintf(fp, "  Nmlt = %d turns", track.Nmlt);
    }
    fprintf(fp, "\n  Fast beam ion interaction enabled:   %d", track.EnableFBII);
    if (ring.m_aHC != 1)
    {
      fprintf(fp, "\n  Ideal HC:   m = %d,    phi_n =  %.2e,   k =  %.2e", ring.m_aHC, ring.phi_n, ring.HC_k);
    }
    fprintf(fp, "\n");
    if(track.EnableDiffCurr)
      fprintf(fp, "  Diff. currents for even and odd bunches enabled, ratio = %.2e \n", track.current_ratio);
    fprintf(fp, "\n");
    
    return 0;
  }
}


int fprintf_bunch_macroparticle_model(FILE * fp, const bunch_macroparticle_model_t model, const int TrackPlane[3])
{
  if(fp == NULL)
  {
    return -1;
  }
  else
  {
    fprintf(fp, "  ===================================================================\n");
    fprintf(fp, "    Beam Distribution Parameters:                                    \n");
    fprintf(fp, "  ===================================================================\n");
    fprintf(fp, "  Np        = %5d\n", model.Np);
    if(model.nGen[LON]) 
    {
      fprintf(fp, "  xtauCM_offset = %8.6lf [nsec],    xepsCM_offset  = %8.6lf\n",
            model.pos_offset.xtau * FGIGA, model.slope_offset.xtau);
      fprintf(fp, "  sgm_xtau      = %8.6lf [nsec],    sgm_xeps       = %8.6lf\n",
            model.pos_sgm.xtau * FGIGA, model.slope_sgm.xtau);
    }
    else
      fprintf(fp, "  LON distribution read from file:   %s\n",
            model.gendst_datafile[LON]);
    if(TrackPlane[HOR])
    {
      if(model.nGen[HOR]) 
      {
        fprintf(fp, "  xxCM_offset   = %8.6lf [mm],      xpCM_offset    = %8.6lf [mrad]\n",
            model.pos_offset.x * FKILO, model.slope_offset.x * FKILO);
        fprintf(fp, "  sgm_xx        = %8.6lf [mm],      sgm_xp         = %8.6lf [mrad]\n",
            model.pos_sgm.x * FKILO, model.slope_sgm.x * FKILO);
      }
      else
        fprintf(fp, "  HOR distribution read from file:   %s\n",
            model.gendst_datafile[HOR]);
    }
    if(TrackPlane[VER]) 
    {
      if(model.nGen[VER]) 
      {
        fprintf(fp, "  zzCM_offset   = %8.6lf [mm],      zpCM_offset    = %8.6lf [mrad]\n",
            model.pos_offset.z * FKILO, model.slope_offset.z * FKILO);
        fprintf(fp, "  sgm_zz        = %8.6lf [mm],      sgm_zp         = %8.6lf [mrad]\n",
            model.pos_sgm.z * FKILO, model.slope_sgm.z * FKILO);
      }    
      else
        fprintf(fp, "  VER distribution read from file:   %s\n",
            model.gendst_datafile[VER]);
    }
    fprintf(fp, "\n");
    return 0;
  }
}

bool e_beam_setup(tracking_t * track, ring_t * ring, e_beam_t * ebeam)
{
  int i;
  unsigned int kb;
  
  int filling = track->filling;

  if(ring == NULL)
  {
    fprintf(stderr, "ring is NULL\n");
    return false;
  }
  else if(ebeam == NULL)
  {
    fprintf(stderr, "ebeam is NULL\n");
    return false;
  }
  
  for(i=0; i<1000; i++)
  {
    ebeam->nfFill[i] = 0;
  }
  
  /* Beam filling */
  /* !! No fills with two gaps possible, filled bunches are expected to come one after the other !!! */
  /* otherwise problems in the ampinv statistics ! */
  if(filling == uniform)     
  {
    for(kb=0; kb < ring->Nharm; kb++) ebeam->nfFill[kb] = 1;
    track->Nbunch_out = 1;
    track->bunch_out = (int *) malloc(track->Nbunch_out * sizeof(int));
    track->bunch_out[0] = 0;
  }
  else if(filling == onefourth)   
  {
    for(kb=0; kb < ring->Nharm/4; kb++)   ebeam->nfFill[kb] = 1;
    track->Nbunch_out = 3;
    track->bunch_out = (int *) malloc(track->Nbunch_out * sizeof(int));
    track->bunch_out[0] = 0;
    track->bunch_out[1] = (int)ring->Nharm/8;
    track->bunch_out[2] = (int)ring->Nharm/4 - 1;    
  }
  else if(filling == onethird)    
  { 
    for(kb=0; kb < ring->Nharm/3; kb++)   ebeam->nfFill[kb] = 1;
    track->Nbunch_out = 3;
    track->bunch_out = (int *) malloc(track->Nbunch_out * sizeof(int));
    track->bunch_out[0] = 0;
    track->bunch_out[1] = (int)ring->Nharm/6;
    track->bunch_out[2] = (int)ring->Nharm/3 - 1; 
  }
  else if(filling == twothird)    
  { 
    for(kb=0; kb < 2*ring->Nharm/3; kb++) ebeam->nfFill[kb] = 1;
    track->Nbunch_out = 3;
    track->bunch_out = (int *) malloc(track->Nbunch_out * sizeof(int));
    track->bunch_out[0] = 0;
    track->bunch_out[1] = (int)ring->Nharm/3;
    track->bunch_out[2] = (int)2*ring->Nharm/3 - 1; 
  }
  else if(filling == threefourth) 
  {
    for(kb=0; kb < 3*ring->Nharm/4; kb++) ebeam->nfFill[kb] = 1;
    track->Nbunch_out = 3;
    track->bunch_out = (int *) malloc(track->Nbunch_out * sizeof(int));
    track->bunch_out[0] = 0;
    track->bunch_out[1] = (int)3*ring->Nharm/8;
    track->bunch_out[2] = (int)3*ring->Nharm/4 - 1; 
  }
  else if(filling == single)      
  {
    for(kb=0; kb<1; kb++)         ebeam->nfFill[kb] = 1;
    track->Nbunch_out = 1;
    track->bunch_out = (int *) malloc(track->Nbunch_out * sizeof(int));
    track->bunch_out[0] = 0;
  }
  else if(filling == single2)     
  {
    for(kb=0; kb<2; kb++)         ebeam->nfFill[kb] = 1;
    track->Nbunch_out = 2;
    track->bunch_out = (int *) malloc(track->Nbunch_out * sizeof(int));
    track->bunch_out[0] = 0;
    track->bunch_out[1] = 1;
  }
  else if(filling == single3)     
  {
    for(kb=0; kb<3; kb++)         ebeam->nfFill[kb] = 1;
    track->Nbunch_out = 3;
    track->bunch_out = (int *) malloc(track->Nbunch_out * sizeof(int));
    track->bunch_out[0] = 0;
    track->bunch_out[1] = 1;
    track->bunch_out[2] = 2; 
  }
  else if(filling == single10)    
  {
    for(kb=0; kb<10; kb++)        ebeam->nfFill[kb] = 1;
    track->Nbunch_out = 3;
    track->bunch_out = (int *) malloc(track->Nbunch_out * sizeof(int));
    track->bunch_out[0] = 0;
    track->bunch_out[1] = 4;
    track->bunch_out[2] = 9; 
  }
  
  else if(filling == single20)    
  {
    for(kb=0; kb<20; kb++)        ebeam->nfFill[kb] = 1;
    track->Nbunch_out = 3;
    track->bunch_out = (int *) malloc(track->Nbunch_out * sizeof(int));
    track->bunch_out[0] = 0;
    track->bunch_out[1] = 9;
    track->bunch_out[2] = 19; 
  }
  else if(filling == single30)  
  {
    for(kb=0; kb<30; kb++)        ebeam->nfFill[kb] = 1;
    track->Nbunch_out = 3;
    track->bunch_out = (int *) malloc(track->Nbunch_out * sizeof(int));
    track->bunch_out[0] = 0;
    track->bunch_out[1] = 14;
    track->bunch_out[2] = 29; 
  }
  else if(filling == evenhalf)
  {
    for (kb=0; kb<ring->Nharm/2; kb++) ebeam->nfFill[2*kb] = 1;
    track->Nbunch_out = 1;
    track->bunch_out = (int *) malloc(track->Nbunch_out * sizeof(int));
    track->bunch_out[0] = 0;
  }
  else if(filling == evenfourth)
  {
    for (kb=0; kb<ring->Nharm/4; kb++) ebeam->nfFill[4*kb] = 1;
    track->Nbunch_out = 1;
    track->bunch_out = (int *) malloc(track->Nbunch_out * sizeof(int));
    track->bunch_out[0] = 0;
  }
  else return false;
  
  ebeam->Nbunch = 0;
  for(kb=0; kb<ring->Nharm; kb++) if(ebeam->nfFill[kb]) (ebeam->Nbunch)++;
  
  return true;
}

int fprint_e_beam(FILE * fp, const ring_t ring, const e_beam_t ebeam)
{
  if(fp == NULL)
  {
    return -1;
  }
  else
  {
    unsigned int kb = 0;
    fprintf(fp, "  ===================================================================\n");
    fprintf(fp, "    Electron Beam Filling                                            \n");
    fprintf(fp, "  ===================================================================\n");
    fprintf(fp, "%5d   %5d\n", kb, ebeam.nfFill[kb]);
    int prec = ebeam.nfFill[kb];
    for(kb=1; kb < ring.Nharm; kb++)
    {
      if(ebeam.nfFill[kb] != prec)
      {
        fprintf(fp, "  ...   %5d\n", prec);
        fprintf(fp, "%5u   %5d\n", kb-1, prec);
        fprintf(fp, "%5u   %5d\n", kb, ebeam.nfFill[kb]);
        prec = ebeam.nfFill[kb];
      }
    }
    if(ebeam.nfFill[kb] == prec)
    {
      fprintf(fp, "  ...   %5d\n", prec);
      fprintf(fp, "%5u   %5d\n", kb, ebeam.nfFill[kb]);
    }
    return 0;
  }
}

int e_beam_print(const ring_t ring, const e_beam_t ebeam)
{
  return fprint_e_beam(stdout, ring, ebeam);
}


bool setup_ring_parameters(ring_t * ring)
{ 
  /* physic constants */
  
  const double Cg     = 8.85e-05;    /***  [meter GeV-3] ***/
  
  /* Parameters */

  ring->T0      = ring->Lc / C_LIGHT; /*SI*/
  ring->R       = ring->Lc / (2.0*M_PI); /*SI*/
  ring->w0      = (2.0*M_PI) / ring->T0; /*SI*/
  ring->Gamma   = 1957.0*ring->E0; /* E0 in GeV needed */  /*SI*/
  ring->Gamma2  = ring->Gamma * ring->Gamma; /*SI*/
  ring->Gamma3  = ring->Gamma * ring->Gamma2; /*SI*/
  ring->wrf     = 2.0*M_PI * ring->frf * FMEGA;  /*SI*/
  //ring->h       = ((int) (ring->wrf / ring->w0));  
  ring->h       = ((int) (ring->frf * FMEGA * ring->Lc / C_LIGHT + 0.5));  
  ring->Nharm   = ring->h;
  
  /*** U0 : Energy Loss per Turn [keV]  ***/
  /*** Urad : Radiation Loss Term used in the tracking ***/
  /*** Longitudinal Damping Time etc. ***/ 
  if (ring->U0 >= 0 && ring->taue >= 0) 
  {
    /* if (ring->Je > 0) 
      fprintf(stderr, "WARNING: Je is not take from input but calculated from given U0 and taue !\n"); */
    ring->Urad    = ring->U0 / (ring->E0 * FGIGA);
    if (ring->Urad!=0) {
      if (ring->Je < 0)
	ring->Je      = 2 * ring->T0 / (ring->taue * ring->Urad);
      if (ring->rho0 < 0)
	ring->rho0    = Cg * ring->E0 * ring->E0 * ring->E0 * ring->E0 * FGIGA / ring->U0;
    }
    ring->De      = 2*ring->T0/(ring->taue);
    ring->tauC    = ring->Je * ring->taue;
    ring->aexpe   = 1.0 / ring->taue;
  }
  else
  {
    ring->U0      = Cg * ring->E0 * ring->E0 * ring->E0 * ring->E0 * FGIGA / ring->rho0; /* E0 in GeV needed */
    ring->Urad    = ring->U0 / (ring->E0 * FGIGA);
    ring->tauC    = 4.0 * M_PI * ring->R * ring->rho0 / (Cg*C_LIGHT) / pow(ring->E0, 3); /* E0 in GeV needed */
    ring->taue    = ring->tauC / ring->Je;
    ring->aexpe   = 1.0 / ring->taue;
    ring->De      = 2.0 * ring->aexpe * ring->T0;
  }
  
  ring->q       = ring->U0 / ring->Vrf0 / FMEGA ;
  ring->Fq      = 2.0*(sqrt(1.0 /ring->q / ring->q - 1.0) - acos(ring->q));
  ring->epsMax  = ring->U0 / (M_PI * ring->ac * ring->h * ring->E0 * FGIGA) * ring->Fq;
  ring->epsMax  = sqrt(ring->epsMax);  
  ring->Vrfp    = FMEGA * ring->Vrf0 * ring->wrf * sqrt(1.0 - ring->q * ring->q);
  ring->wso     = sqrt(fabs(ring->ac) * ring->Vrfp / ( ring->T0 * ring->E0 * FGIGA));
  ring->wso2    = ring->wso * ring->wso;
  ring->fso     = ring->wso / (2.0 * M_PI);
  ring->Qso     = ring->fso * ring->T0;
  ring->wgziH   = ring->QH0 * ring->w0 * ring->Gzix / ring->ac;
  ring->wgziV   = ring->QV0 * ring->w0 * ring->Gziz / ring->ac;


  ring->sn0     = ring->q;
  ring->phai0   = asin( ring->q );
  ring->fc1     = ring->ac / ring->wso;
  ring->fc12    = pow(ring->fc1, 2);

  /*** Transverse Damping Times etc. ***/
  
  ring->Jz      = 1.0;
  ring->Jx      = 4.0 - ring->Jz - ring->Je;
  ring->taux    = ring->tauC / ring->Jx;
  ring->tauz    = ring->tauC / ring->Jz;
  
  /* Optical parameters */
  
  ring->gamma1[HOR] = (1.0 + ring->alpha1[HOR] * ring->alpha1[HOR])/ring->beta1[HOR];
  ring->gamma1[VER] = (1.0 + ring->alpha1[VER] * ring->alpha1[VER])/ring->beta1[VER];
  
    /* Longrange resonator parameters & constants */
  int i;
  double tmp;
  unsigned Ntmp = 0;
  double mult =1.;
  double genphase = 0;
  for(i = 1; i <= ring->longrange_resonators_size; i++)
  {
    LR_resonator_t * lr_resonator = &(ring->longrange_resonators[i-1]);
    lr_resonator->wr = lr_resonator->m * ring->wrf + lr_resonator->detune * 2 * M_PI;
    if (lr_resonator->m>1) mult *= (double)lr_resonator->m*lr_resonator->m / (lr_resonator->m*lr_resonator->m - 1.);  
    else if (!ring->has_rf_feedback) genphase = atan(lr_resonator->Qfactor * ( lr_resonator->wr/ring->wrf - ring->wrf/lr_resonator->wr ));
    tmp = (lr_resonator->wr * 0.5 / lr_resonator->Qfactor);
    lr_resonator->Nturn = (unsigned) (log(2) * 10 / tmp / ring->T0) + 1;
    lr_resonator->Nbu = lr_resonator->Nturn * ring->h;
    if(lr_resonator->Nbu < 2)
    {
      fprintf(stderr, "Given longrange resonator is shortrange!\n");
      return false;
    }
    if (lr_resonator->Nbu > Ntmp)
      Ntmp = lr_resonator->Nbu;
  }
  ring->phai0 = asin(mult * ring->q) - genphase/2.0;
  ring->Nbumax = Ntmp;
  ring->lr_order = 6;

  if (ring->ac<0) ring->phai0 = M_PI-ring->phai0;
   
  return true;
}


bool setup_tracking_parameters(ring_t * ring, tracking_t * track, selffield_model_t * SelfFieldModel, bunch_macroparticle_model_t * macrop_model)
{
  if(track->EnableIdealHC == 1)
  {
    ring->phai0 = asin(ring->m_aHC*ring->m_aHC /(ring->m_aHC*ring->m_aHC - 1.) * ring->q);
    ring->phi_n = atan(tan(ring->phai0) / (double) ring->m_aHC) / (double) ring->m_aHC;
    if (ring->ac<0) {
      ring->phai0 = M_PI-ring->phai0;
      ring->phi_n = M_PI-ring->phi_n;
    }
    ring->HC_k = -cos(ring->phai0) / (ring->m_aHC * cos(ring->m_aHC * ring->phi_n));  
  }
  if(track->scan == 1)
    { //TODO: Chroma scan
      ring->Iring = track->scan_start;
    }
  if(track->scan == 3)
    {
      if(ring->longrange_resonators_size < 1)
      {
        fprintf(stderr, "Qfactor scan but no longrange resonator given!\n");
        return false;
      }
      track->scan_start = 1./track->Nscan;
      track->scan_step = 1./track->Nscan;
      unsigned int i;
      for(i = 0; i < ring->longrange_resonators_size; i++)
        ring->longrange_resonators[i].Qfactor = ring->longrange_resonators[i].Qfactor / track->Nscan;
    }
    if(track->scan == 4)
    {
      if(ring->longrange_resonators_size < 1)
      {
        fprintf(stderr, "Qfactor scan but no longrange resonator given!\n");
        return false;
      }
      track->scan_step = (1. - track->scan_start) / (track->Nscan - 1);
      unsigned int i;
      for(i = 0; i < ring->longrange_resonators_size; i++)
        ring->longrange_resonators[i].Qfactor = ring->longrange_resonators[i].Qfactor * track->scan_start;
    }
    
    double b = ring->beffL[1];
    if(ring->beffH[1] > b && track->TrackPlane[1] == 1) b = ring->beffH[1];
    else if(ring->beffV[1] > b && track->TrackPlane[2] == 1) b = ring->beffV[1];

    if (track->EnableRW_short)
    { // RW ala Bane

      track->s0 = pow((2 * b * b * ring->rhorw / Z_0), 1.0/3.0);
      double dT = 2.0 * SelfFieldModel->Nsigma / ((double) SelfFieldModel->Ncell - 1);
      double rel_bin_with = dT * macrop_model->pos_sgm.xtau * C_LIGHT;
      track->mult = rel_bin_with / track->s0;
        if (track->mult < 20.)
        {
          fprintf(stderr, "Selffield Model has too small binwidth for RW treatment!\n");
          return false;
        }
    }
    
    if(track->EnableRW_long && track->Nmlt == 0)
    {
      double tmp;
      tmp = ring->Lc * sqrt(ring->rhorw) / pow(b, 3);
      track->Nmlt = 1 + ((int) 500e3 / tmp); /* Factor changed from 333 to 500 since the max b from all planes is taken */
      track->NmultiT = track->Nmlt+2;
      if (track->Nmlt == 0)
        {
          fprintf(stderr, "Given RW is to shortrange, RW_long will be switched off!\n");
          track->EnableRW_long = 0;
        }
    }
  
  return true;
}


bool macrop_model_setup_parameters(const ring_t ring, const tracking_t tracking, 
                           bunch_macroparticle_model_t * macrop_model)
{
  /* physic constants */
  
  const double Cq     = 3.84e-13;    /***  [meter]       ***/
  
  /* beam parameters */
  
  if(tracking.TrackPlane[LON] && macrop_model->sgm_xeps < 0.0)
    macrop_model->sgm_xeps = sqrt(Cq * ring.Gamma2 / (ring.Je * ring.rho0));
  if(tracking.TrackPlane[LON] && macrop_model->sgm_xtau < 0.0)
    macrop_model->sgm_xtau = FGIGA * fabs(ring.ac) / ring.wso * macrop_model->sgm_xeps;

  if(macrop_model->sgm_xx < 0.0)
    macrop_model->sgm_xx = FKILO*sqrt(ring.emittanceH/FGIGA *  ring.beta1[HOR]);
  if(macrop_model->sgm_xp < 0.0)
    macrop_model->sgm_xp = FKILO*sqrt(ring.emittanceH/FGIGA * ring.gamma1[HOR]);
  if(macrop_model->sgm_zz < 0.0)
    macrop_model->sgm_zz = FKILO*sqrt(ring.emittanceH/FGIGA *  ring.beta1[VER] * ring.couplbeta);
  if(macrop_model->sgm_zp < 0.0)
     macrop_model->sgm_zp = FKILO*sqrt(ring.emittanceH/FGIGA * ring.gamma1[VER] * ring.couplbeta);
  
  macrop_model->sgmatau = macrop_model->sgm_xtau/FGIGA;
  
        /* Offsets and sgms */
    macrop_model->pos_offset.xtau = macrop_model->xtauCM_offset * FNANO;
    macrop_model->pos_offset.x = macrop_model->xxCM_offset * FMILLI;
    macrop_model->pos_offset.z = macrop_model->zzCM_offset * FMILLI;

    macrop_model->slope_offset.xtau = macrop_model->xepsCM_offset;
    macrop_model->slope_offset.x = macrop_model->xpCM_offset * FMILLI;
    macrop_model->slope_offset.z = macrop_model->zpCM_offset * FMILLI;
      
    macrop_model->pos_sgm.xtau = macrop_model->sgm_xtau * FNANO;
    macrop_model->pos_sgm.x = macrop_model->sgm_xx * FMILLI;
    macrop_model->pos_sgm.z = macrop_model->sgm_zz * FMILLI;

    macrop_model->slope_sgm.xtau = macrop_model->sgm_xeps;
    macrop_model->slope_sgm.x = macrop_model->sgm_xp * FMILLI;
    macrop_model->slope_sgm.z = macrop_model->sgm_zp * FMILLI;
  
  return true;
}

int fprint_parameters(FILE * fp, const ring_t ring,
                      const bunch_macroparticle_model_t macrop_model)
{
  fprintf(fp, "  ===================================================================\n");
  fprintf(fp, "    Summary ring parameters, determined from inputs                  \n");
  fprintf(fp, "  ===================================================================\n");
  
  fprintf(fp, "\n  rel. gamma =  %9.lf,      T0 =   %9.lf [nsec]", ring.Gamma, ring.T0 * FGIGA);
  fprintf(fp, "\n  av. radius =  %9.2lf [m],  Urad =   %9.3lf  [MeV]", ring.R, ring.U0/FMEGA);
//fprintf(fp, "\n  Element Name at the Transformation Point: %s", ge[ib0]->name);
  fprintf(fp, "\n  beta1[u]  [m]   (u = hor, ver):  %9.5lf   %9.5lf", ring.beta1[HOR], ring.beta1[VER]);
  fprintf(fp, "\n  alpha1[u]       (u = hor, ver):  %9.5lf   %9.5lf", ring.alpha1[HOR], ring.alpha1[VER]);
  fprintf(fp, "\n  gamma1[u] [m-1] (u = hor, ver):  %9.5lf   %9.5lf", ring.gamma1[HOR], ring.gamma1[VER]);
  fprintf(fp, "\n  sgm_xtau = %8.6lf [ns],   sgm_xeps = %8.6lf", macrop_model.sgm_xtau, macrop_model.sgm_xeps);
  fprintf(fp, "\n  sgm_xx   = %8.6lf [mm],   sgm_xp   = %8.6lf [mrad]", macrop_model.sgm_xx, macrop_model.sgm_xp);
  fprintf(fp, "\n  sgm_zz   = %8.6lf [mm],   sgm_zp   = %8.6lf [mrad]", macrop_model.sgm_zz, macrop_model.sgm_zp);
  fprintf(fp, "\n  taue  = %8.4lf [msec],   T0/taue = %9.6lf", FKILO*ring.taue, ring.T0/ring.taue);
  fprintf(fp, "\n  taux  = %8.4lf [msec],   T0/taux = %9.6lf", FKILO*ring.taux, ring.T0/ring.taux);
  fprintf(fp, "\n  tauz  = %8.4lf [msec],   T0/tauz = %9.6lf", FKILO*ring.tauz, ring.T0/ring.tauz);
  fprintf(fp, "\n  QH = %8.3lf,   QV = %8.3lf", ring.QH0, ring.QV0);
  fprintf(fp, "\n  Qso   = %9.5lf,  fso = %9.5lf [kHz]", ring.Qso, ring.fso/FKILO);
  fprintf(fp, "\n  U0    = %9.5lf [keV],   phi0 = %9.5lf [rad]", ring.U0/FKILO, ring.phai0);
  fprintf(fp, "\n  wgziH = %9.5lf [GHz]", ring.wgziH/FGIGA);
  fprintf(fp, "\n  wgziV = %9.5lf [GHz]", ring.wgziV/FGIGA);
  fprintf(fp, "\n  fc1   = %12.5e [sec]", ring.fc1);
  fprintf(fp, "\n  Nharm = %5u,    (diff int - double: %e)", ring.Nharm, ring.Nharm - ring.frf * FMEGA * ring.Lc / C_LIGHT);
  fprintf(fp, "\n");
  fprintf(fp, "\n");
  return 0;
}



int
fprint_resonators(FILE * fp, const ring_t ring, const selffield_model_t SelfFieldModel, const bunch_macroparticle_model_t macrop_model)
{
     unsigned k;
     double dT = 2.0 * macrop_model.sgm_xtau * (double) SelfFieldModel.Nsigma / ((double) SelfFieldModel.Ncell - 1.);
     fprintf(fp, "\n  ===================================================================");
     fprintf(fp, "\n    Ring Impedance Inputs:                                           ");
     fprintf(fp, "\n  ===================================================================");
     fprintf(fp, "\n  Pure resistive / pure inductive impedance: ");
     fprintf(fp, "\n  ResistL0      = %8.2lf [ohm],     aindL0      = %8.2lf [ohm*GHz-1]",SelfFieldModel.RsisL ,SelfFieldModel.aindL*2.0*M_PI*FGIGA);
     fprintf(fp, "\n  Total number of BBRs:    Nbbr = %2d \n", ring.resonators_size);
     if(ring.resonators_size > 0)
     {
       fprintf(fp, "\n    n    R(n)[kohm]     Q(n)      fr(n)[GHz]     plane(n)");
       fprintf(fp, "\n  ===================================================================");
     }
     for(k = 0; k < ring.resonators_size; k++) 
     {
        const resonator_t resonator = ring.resonators[k];
        fprintf(fp, "\n  %2d     %8.2lf    %8.2lf     %8.2lf       %d",  k+1, resonator.Rs/FKILO, resonator.Qfactor, resonator.wr/(2*M_PI*FGIGA), resonator.plane);  
     }
       fprintf(fp, "\n  Selffield:\n  Ncell = %2d,     Nsigma = %.2f,   binwidth = %.2e [ns]\n", SelfFieldModel.Ncell, SelfFieldModel.Nsigma, dT);   
  fprintf(fp, " ===================================================================\n");
  fprintf(fp, "   Harmonic Cavity Input: \n");
  fprintf(fp, " ===================================================================\n");
  fprintf(fp, "  Passive HCs, number of LR resonators: %u\n", ring.longrange_resonators_size);
  if (ring.longrange_resonators_size > 0)
    fprintf(fp, "  Nbu max    = %d = %d turns,  order:  %u\n", ring.Nbumax, (int)ring.Nbumax/ring.Nharm, ring.lr_order);
  for(k = 0; k < ring.longrange_resonators_size; k++) 
  {
    const LR_resonator_t resonator = ring.longrange_resonators[k];
    fprintf(fp, "  R [%d]     = %8.2lf [Mohm],         Q [%d]      = %8.2lf, \n",
            k, resonator.Rs/FMEGA, k, resonator.Qfactor);
    fprintf(fp, "  wr [%d]    = %8.2lf [GHz rad],     detuning: %8.2lf [kHz]\n",
            k, resonator.wr/FGIGA, resonator.detune/1000.);  
    fprintf(fp, "  Nbu [%d]    = %d = %d turns,          harm.: %d\n",
            k, resonator.Nbu, resonator.Nturn, resonator.m);
  } 
  fprintf(fp, "  Active HCs, number of cavities: %u\n", ring.active_HC_size);
  for(k = 0; k < ring.active_HC_size; k++) 
  {
    const active_HC_t aHC = ring.active_HC[k];
    fprintf(fp, "  Vpeak [%d] = %8.2lf [MV],  nHC [%d] = %8.2lf,  phase [%d] = %8.2lf [rad] \n",
            k, aHC.Vpeak/FMEGA, k, aHC.nHC, k, aHC.phi_aHC);
  }
  return 0; 
}


int ring_destroy(ring_t * ring)
{
  int i;
  for(i = 0; i < ring->resonators_size; i++)
  {  
    free(ring->resonators);
    ring->resonators = NULL;
  }
  
  for(i = 0; i < ring->longrange_resonators_size; i++)
  {
    LR_resonator_t * ResonImp = &ring->longrange_resonators[i];
    
    ResonImp->Rs = 0;
    ResonImp->wr = 0;
    ResonImp->Qfactor = 0;
    ResonImp->Nbu = 0;
    
    free(ResonImp->facc);
    ResonImp->facc = NULL;
    
    free(ResonImp->facs);
    ResonImp->facs = NULL;
  }  
  
  if(ring->longrange_resonators_size > 0)
  {
    free(ring->lr_wake);
    ring->lr_wake = NULL; 
  }
    
  return 1;
}


int
tracking_destroy(tracking_t * track)
{
  free(track->bunch_out);
  return 1;
}






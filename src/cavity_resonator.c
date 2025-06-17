//   mbtrack
//   Copyright (C) 2000-2019, Synchrotron SOLEIL
//   Copyright (C) 2018-2019, High Energy Accelerator Research Organisation, KEK
//
//   This program is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "def.h"
#include "types.h"
#include "input.h"
#include "transform_weak.h"
#include "cavity_resonator.h"

/* Global variables */
//extern ring_t ring;
extern tracking_t track;

static const char * const cavScan_str[]=
{
  /* 0 */ "detuneRAD",
  /* 1 */ "fr"
};

bool scanCavType(char * str, cavScanType_t * scanType)
{
  int i;
  for(i = 0; i < 2; i++)
  {
    if(strcmp(cavScan_str[i], str) == 0)
    {
      *scanType = i;
      return true;
    }
  }
  return false;
}


bool read_conf_cav_res(FILE * fp, ring_t * ring, tracking_t * track)
{
  /*
   * [cavity_resonator_1], [cavity_resonator_2], ... , 09/07/2018, created by Naoto Yamamoto
   */
  int i,itmp;
  double tmp,frev;
  char section[32] = "";

  ring->cavity_resonators = (CAVITY_resonator_t *) malloc(1 * sizeof(CAVITY_resonator_t));
  i = 1;
  ring->cav_resonators_drfFB_active = 0 ;
  ring->cavity_resonator_main = 0;
  snprintf(section, 32, "cavity_resonator_%d", i);
  while(config_have_section(fp, section))
  {
    ring->cavity_resonators = (CAVITY_resonator_t *) realloc(ring->cavity_resonators, i * sizeof(CAVITY_resonator_t));
    CAVITY_resonator_t * cav_resonator = &(ring->cavity_resonators[i-1]);
    /////////////////////////////////////////////////////////
    // cavity setup 
    /////////////////////////////////////////////////////////
    if(!config_get_double(fp, section, "m", &(cav_resonator->m)) || cav_resonator->m < 0)
    {
      fprintf(stderr, "m (m > 0) is not found. \n");
      return false;
    }
    //cav_resonator->mInt = (int)(cav_resonator->m + 0.5);
    frev = ring->frf * FMEGA / ((int) (ring->frf * FMEGA * ring->Lc / C_LIGHT + 0.5));//f0
    cav_resonator->wref = 2* M_PI * (int)(cav_resonator->m * ring->frf * FMEGA / frev + 0.5) * frev ;
    if(!config_get_double(fp, section, "Rs", &(cav_resonator->Rs)))  return false;
    if(!config_get_double(fp, section, "Qload", &(cav_resonator->Qload))) return false;
    if(!config_get_double(fp, section, "Qzero", &(cav_resonator->Qzero))) return false;

    /////////////////////////////////////////////////////////
    // operation parameter setup
    ////////////////////////////////////////////////////////
    config_get_fprint = false; /* no error printing */
    // Optional
    cav_resonator->switchOnTurn=0;
    cav_resonator->switchOffTurn=track->NrevTot+1;
    if(config_get_double(fp, section, "switchOnTurn", &tmp)) cav_resonator->switchOnTurn = (int)tmp;
    if(config_get_double(fp, section, "switchOffTurn", &tmp)) 
	    if(tmp > cav_resonator->switchOnTurn ) cav_resonator->switchOffTurn = (int)tmp;
    
    cav_resonator->mode = 0; // default is without FB.
    config_get_int(fp, section, "FBmode", &(cav_resonator->mode));
    ///////////////////////////////////////////////////////
    // For automatic generator paramter
    cav_resonator->VcAP[0] = 0;
    cav_resonator->VcAP[1] = 0;
    if(config_get_double(fp, section, "Vc", &tmp))
	cav_resonator->VcAP[0] = tmp / (ring->E0 * FGIGA);
    if (cav_resonator->VcAP[0] > 0)	
    {
      if(!config_get_double(fp, section, "syncPhase", &(cav_resonator->VcAP[1])))
      {
    	fprintf(stderr, "If Vc > 0, syncPhase is mandatory. \n");
 	return false;
      }
    }
    /////////////////////////////////////////////////////////
    // For detuning angle
    cav_resonator->FBparaIniSet = 0; //190820 for cavity feedback NY
    cav_resonator->detune = -10;
    if(config_get_double(fp, section, "detuneRad", &(cav_resonator->detune)))
    {
      cav_resonator->wr = 0.5 * cav_resonator->wref *(tan(cav_resonator->detune) / cav_resonator->Qload+sqrt( tan(cav_resonator->detune) / cav_resonator->Qload * tan(cav_resonator->detune) / cav_resonator->Qload + 4 ));
    }
    else if(config_get_double(fp, section, "detuneHz", &tmp))
    {
      cav_resonator->wr = cav_resonator->m * 2.0 * M_PI * ring->frf * FMEGA +  2.0 * M_PI * tmp;
      cav_resonator->detune = atan(cav_resonator->Qload * (cav_resonator->wr / cav_resonator->wref - cav_resonator->wref / cav_resonator->wr));
    }
    if(cav_resonator->VcAP[0] == 0 && cav_resonator->detune == -10 && fabs((int)(cav_resonator->m*2)-2*cav_resonator->m) < 1e-20)
    {
      fprintf(stderr, "If Vc = 0 and m is (half)integer, detuneHz or detuneRad is mandatory. \n");
      return false;
    }
    ////////////////////////////////////////////////////////////
    // For detuning offset, this parameter is used in auto generator parameter calulation mode
    cav_resonator->detuneOffset = 0.0;
    cav_resonator->wrOffset = 0.0;
    if(config_get_double(fp, section, "detuneOffsetHz", &tmp))
      cav_resonator->wrOffset = tmp * 2.0 * M_PI; // [Hz rad]
    else //if(cav_resonator->wrOffset == 0.0) 
    {
      if(config_get_double(fp, section, "detuneOffsetDeg", &tmp))
         cav_resonator->detuneOffset = tmp / 180 * M_PI;
      config_get_double(fp, section, "detuneOffsetRad", &(cav_resonator->detuneOffset));
    }
    /////////////////////////////////////////////////////////////
    // For manual (active) generator parameter
    cav_resonator->geneVolt = 1;// default is active cavity
    cav_resonator->geneVoltRes[0] = 0;
    cav_resonator->genePhase = 0;
    if(config_get_double(fp, section, "geneVoltRes[0]", &tmp)) 
    {
      cav_resonator->geneVoltRes[0] = tmp / (ring->E0 * FGIGA);
      cav_resonator->geneVolt = fabs(cav_resonator->geneVoltRes[0] * cos(cav_resonator->detune));
      if(cav_resonator->geneVolt > 0)
      {
	if(!config_get_double(fp, section, "genePhase", &(cav_resonator->genePhase)))
	{
          fprintf(stderr, "If geneVolt > 0, genePhase is mandatory.\n");
      	  return false;
	}
      }
    } else if(config_get_double(fp, section, "geneVolt", &tmp)) 
    {
      cav_resonator->geneVolt = tmp / (ring->E0 * FGIGA);
      cav_resonator->geneVoltRes[0] = fabs(cav_resonator->geneVolt / cos(cav_resonator->detune));
      if(cav_resonator->geneVolt > 0)
      {
	if(!config_get_double(fp, section, "genePhase", &(cav_resonator->genePhase)))
	{
          fprintf(stderr, "If geneVolt > 0, genePhase is mandatory.\n");
      	  return false;
	}
      }
    }
    cav_resonator->EnableCavInst = 1;
    track->EnableCavReCalc = 0;
    config_get_int(fp, section, "EnableCavInst", &(cav_resonator->EnableCavInst));
    /////////////////////////////////////////////////////////////////
    // For cavity feedback  
    cav_resonator->fbGain = 1.0;
    cavityRFVoltage_feedback_t * cavFB = &(cav_resonator->cavVolFB);
    cavFB->PIcontrol=0;
    cavFB->loopDelay=0;
    cavFB->limitRange = 0.0000;
    cavFB->avgNum=1.0;
    cavFB->sample=4.0;
    cavFB->switchOnTurn=0.0;
    cav_resonator->cavVolFB.switchOffTurn= cav_resonator->switchOffTurn;
    //config_get_int(fp, section, "FBPIcontrol", &itmp);
    if(config_get_double(fp, section, "FBAvgNum", &tmp))  if(tmp > 0) cavFB->avgNum=(int)tmp;

    for (int i=0;i<3;i++) 
    {
      cavFB->gainProp[i]=0.1;
      cavFB->gainInte[i]=1.0;
    }
    for (int i=0;i<2;i++)
    {
      cavFB->integralMemory[i]=0;
      cavFB->integralMemoryFac[i]=1/ring->frf / FMEGA;//1.0/frev/ring->Nharm;//1.0/frev;//*cavFB->avgNum; // 1/2^20
    }

    config_get_int(fp, section, "FBPIcontrol", &(cavFB->PIcontrol));
    config_get_int(fp, section, "FBSampleInterval", &(cavFB->sample));
    if(config_get_double(fp, section, "FBloopDelay", &tmp)) 
	    if(tmp > 0) cavFB->loopDelay=(int)tmp;
    //if(config_get_int(fp, section, "FBPIcontrol", &itmp)) if(itmp > 0) cavFB->PIcontrol = 1;
    if(config_get_double(fp, section, "FBPgain", &(cavFB->gainProp[0])))
    {
       cavFB->gainProp[1] = cavFB->gainProp[0];
    }else{
      config_get_double(fp, section, "FBPgainAmp", &(cavFB->gainProp[0]));
      config_get_double(fp, section, "FBPgainPhase", &(cavFB->gainProp[1]));
    }
    if(config_get_double(fp, section, "FBIgain", &(cavFB->gainInte[0])))
    {
       cavFB->gainInte[1] = cavFB->gainInte[0];
    }else{
      config_get_double(fp, section, "FBIgainAmp",   &(cavFB->gainInte[0]));
      config_get_double(fp, section, "FBIgainPhase", &(cavFB->gainInte[1]));
    }
    config_get_double(fp, section, "FBPgainTuner", &(cavFB->gainProp[2]));
    if(config_get_double(fp, section, "FBIntegralFactorAmp", &tmp)) cavFB->integralMemoryFac[0] = tmp;
    if(config_get_double(fp, section, "FBIntegralFactorPhase", &tmp)) cavFB->integralMemoryFac[1] = tmp;

    if(cavFB->PIcontrol==0)
    {
      config_get_double(fp, section, "FBgain", &(cavFB->gainProp[0]));
      if(!config_get_double(fp, section, "FBgainP", &(cavFB->gainProp[1]))) cavFB->gainProp[1] = cavFB->gainProp[0];
      if(!config_get_double(fp, section, "FBgainT", &(cavFB->gainProp[2]))) cavFB->gainProp[2] = cavFB->gainProp[0];
    }

    cavFB->tuneRange =  0.05; // 19/06/20 NY
    cavFB->tuneDeadRange = M_PI/180 * 0.1; //0.1 deg
    config_get_double(fp, section, "tuneRange", &(cavFB->tuneRange));
    cav_resonator->fbMon = track->NrevMon; 
    config_get_double(fp, section, "FBlimit", &(cavFB->limitRange));
    if(config_get_double(fp, section, "FBswitchOnTurn", &tmp)) 
	    if(tmp > 0) cavFB->switchOnTurn=(int)tmp;
    if(config_get_double(fp, section, "FBswitchOffTurn", &tmp)) 
	    if(tmp > cavFB->switchOnTurn) cavFB->switchOffTurn=(int)tmp;
  
    // For Fourier componet of beam load
    cav_resonator->FourierAmp = 1.0;
    cav_resonator->FourierPhase = 0.0;
    config_get_double(fp, section, "FourierAmp", &(cav_resonator->FourierAmp));
    config_get_double(fp, section, "FourierPhase", &(cav_resonator->FourierPhase));

    // For cavity voltage compensation 190809 N.Yamamoto
    cav_resonator->geneModfreqMrev=0;
    config_get_int(fp, section, "geneModfreqMrev", &(cav_resonator->geneModfreqMrev));
    for(itmp=0;itmp<5;itmp++) cav_resonator->geneModPhaseOffsetM[itmp]=0;
    for(itmp=0;itmp<5;itmp++) cav_resonator->geneModVolt[itmp]=0;
    if(config_get_int(fp, section, "geneModPhaseOffsetM", &itmp)) cav_resonator->geneModPhaseOffsetM[0] = itmp;
    if(config_get_int(fp, section, "geneModPhaseOffsetM_1", &itmp)) cav_resonator->geneModPhaseOffsetM[0] = itmp;
    if(config_get_int(fp, section, "geneModPhaseOffsetM_2", &itmp)) cav_resonator->geneModPhaseOffsetM[1] = itmp;
    if(config_get_int(fp, section, "geneModPhaseOffsetM_3", &itmp)) cav_resonator->geneModPhaseOffsetM[2] = itmp;
    if(config_get_double(fp, section, "geneModVolt", &tmp)) cav_resonator->geneModVolt[0] = tmp / (ring->E0 * FGIGA);
    if(config_get_double(fp, section, "geneModVolt_1", &tmp)) cav_resonator->geneModVolt[0] = tmp / (ring->E0 * FGIGA);
    if(config_get_double(fp, section, "geneModVolt_2", &tmp)) cav_resonator->geneModVolt[1] = tmp / (ring->E0 * FGIGA);
    if(config_get_double(fp, section, "geneModVolt_3", &tmp)) cav_resonator->geneModVolt[2] = tmp / (ring->E0 * FGIGA);

   
    // direct RF feedback
    cav_resonator->directRF_feedback.fbAmp=0.0;
    cav_resonator->directRF_feedback.fbPhase=0.0;
    cav_resonator->directRF_feedback.gain=0.0;
    cav_resonator->directRF_feedback.phaseShift=0;
    cav_resonator->directRF_feedback.switchOnTurn=0;
    cav_resonator->directRF_feedback.switchOffTurn=track->NrevTot;
    cav_resonator->directRF_feedback.loopDelay=0.0;
    cav_resonator->directRF_feedback.AmpRatio=0.5;
    cav_resonator->directRF_feedback.Vg[0]=0;
    cav_resonator->directRF_feedback.Vg[1]=0;
    if(cav_resonator->VcAP[0] > 0 && config_get_double(fp, section, "directRF_feedback.gain", &tmp)) {
      cav_resonator->directRF_feedback.gain = tmp ;
      ring->cav_resonators_drfFB_active = 1 ;
    }
    if(config_get_double(fp, section, "directRF_feedback.phaseShift", &tmp)) cav_resonator->directRF_feedback.phaseShift = tmp;
    if(config_get_double(fp, section, "directRF_feedback.phaseShiftDeg", &tmp)) cav_resonator->directRF_feedback.phaseShift = tmp / 180 * M_PI;
    if(config_get_double(fp, section, "directRF_feedback.switchOnTurn", &tmp)) cav_resonator->directRF_feedback.switchOnTurn = (int)tmp;
    if(config_get_double(fp, section, "directRF_feedback.switchOffTurn", &tmp)) 
	    if(tmp > cav_resonator->directRF_feedback.switchOnTurn ) cav_resonator->directRF_feedback.switchOffTurn = (int)tmp;
    if(cav_resonator->directRF_feedback.gain  <= 0)  cav_resonator->directRF_feedback.switchOffTurn = -1;
    if(config_get_double(fp, section, "directRF_feedback.loopDelay", &tmp)) 
	    if(tmp > 0) cav_resonator->directRF_feedback.loopDelay=(int)tmp;
    if(config_get_double(fp, section, "directRF_feedback.AmpRatio", &tmp)) 
	    if(tmp > 0) cav_resonator->directRF_feedback.AmpRatio = tmp;
    if(config_get_double(fp, section, "directRF_feedback.AmpRatiodB", &tmp)) 
	    if(tmp < 0) cav_resonator->directRF_feedback.AmpRatio = pow(10,tmp/20.0);

    //scan; 2021/04/02 N.Yamamoto
    cav_resonator->scan = (cavScan_t *) malloc(1 * sizeof(cavScan_t));
    int scan_num=1;
    char cavScan_strTmp[32]="";
    snprintf(section, 32, "cavity_resonator_%d.scan_%d", i,scan_num);
    while(config_have_section(fp, section))
    {
      cav_resonator->scan = (cavScan_t *) realloc(cav_resonator->scan, scan_num * sizeof(cavScan_t));
      cavScan_t * pscan = &(cav_resonator->scan[scan_num-1]);

      //if(!config_get_int(fp, section, "CAVscan_type", &(pscan->type)))  return false;
      if(!config_get_str(fp, section, "CAVscan_type",cavScan_strTmp) || !scanCavType(cavScan_strTmp,&(pscan->type) )) {
	fprintf(stderr, "ERROR: Unknown scan parameter \"%s\" in section \"%s\" \n", cavScan_strTmp,section);
        return false;
      }
      if(!config_get_double(fp, section, "CAVscan_step", &(pscan->step))) {
	fprintf(stderr, "ERROR: No scan_step parameter (CAVscan_step = XXX) \"%s\" in section \"%s\" \n", cavScan_strTmp,section);
	return false;
      }
      if(config_get_double(fp, section, "CAVscan_nrev", &tmp))  pscan->nrev = (int)tmp;
      if(config_get_double(fp, section, "CAVscan_start", &tmp))  pscan->start = tmp;
      else 
      {
	switch (pscan->type)
	{
	  case detuneRAD: // detune
	    pscan->start =  cav_resonator->detune;
	    //printf("### %lf\n",pscan->start);
	    break;
	  case fr: // wr , NOTICE: It does not work if m is not integer or half integer
	    pscan->start =  cav_resonator->wr;
	    //printf("### %lf\n",pscan->start);
	    break;
	  default:
	    printf("ERROR: CAVscan type = %d is not defined yet!\n",pscan->type);
	    return false;
	}
      }
      scan_num++;
      snprintf(section, 32, "cavity_resonator_%d.scan_%d", i,scan_num);
    }
    cav_resonator->scan_size=scan_num-1;

    if (cav_resonator->geneVolt > 0 && cav_resonator->m == 1)
      ring->cavity_resonator_main = i;
    
    config_get_fprint = true; /* error printing */
    /* next resonator */
    i++;
    snprintf(section, 32, "cavity_resonator_%d", i);
  }
  ring->cavity_resonators_size = i-1;

  return true;
}

bool setup_cavity_parameters(ring_t * ring)
{
  unsigned int i,k;
  double cavbeta,tmp;
  double Vbr; // Beam loading Voltage at resonance
  double Vb[2],optVg[2];// complex voltage, 0: real, 1: imaginary
  unsigned Ntmp = 0;
  for(i = 1; i <= ring->cavity_resonators_size; i++)
  {
    CAVITY_resonator_t * cav_res = &(ring->cavity_resonators[i-1]);
    /* setup the cavity paramter to construct wake phasor (PB. Wilson definision, but Rs=0.5Vc^2/Pc) */
    cavbeta = cav_res->Qzero / cav_res->Qload - 1.0; // cavity coupling, using only in this loop
    if(cavbeta < 0)// error check
    {
      fprintf(stderr, "Qzero should be equal to or larger than Qload!\n");
      return false;
    }
    Vbr= 2 * ring->Iring * cav_res->Rs * cav_res->FourierAmp / (1 + cavbeta) *FMILLI / (ring->E0 * FGIGA); // [V]

    ///////////////////////////////////////////////////////////
    // For active cavity, calculation for geneVolt and genePhase 
    if(cav_res->VcAP[0] > 0 )// if Vc is given, 
    {
      cav_res->Vc[0] = cav_res->VcAP[0] * cos(cav_res->VcAP[1]);// set Vc => VcAP; maybe only once
      cav_res->Vc[1] = cav_res->VcAP[0] * sin(cav_res->VcAP[1]);
      if(cav_res->detune == -10) // if detune is not given 
      {// Calculation for Optimum operation condition of a cavity 
 	cav_res->detune = atan(- Vbr/cav_res->VcAP[0] * sin(cav_res->VcAP[1]));

        cav_res->wr = 0.5 * cav_res->wref *(tan(cav_res->detune) / cav_res->Qload 
        +sqrt( tan(cav_res->detune) / cav_res->Qload * tan(cav_res->detune) / cav_res->Qload + 4 ));
        tmp = ring->frf * FMEGA / ((int) (ring->frf * FMEGA * ring->Lc / C_LIGHT + 0.5));//f0
        cav_res->wref = 2* M_PI * (int)(cav_res->m * ring->frf * FMEGA / tmp + 0.5) * tmp ;
      }else{ cav_res->FBparaIniSet=1;}

      ///////////////////////////////////////////////////////////
      // calculation for detuneOffset
      if(cav_res->wrOffset != 0.0) 
      {
        cav_res->wr += cav_res->wrOffset;
        tmp = ring->frf * FMEGA / ((int) (ring->frf * FMEGA * ring->Lc / C_LIGHT + 0.5));//f0
        cav_res->wref = 2* M_PI * (int)(cav_res->m * ring->frf * FMEGA / tmp + 0.5) * tmp ;
        cav_res->detuneOffset = -cav_res->detune;// previous detune angle
	cav_res->detune = atan(cav_res->Qload * (cav_res->wr / cav_res->wref - cav_res->wref / cav_res->wr));
        cav_res->detuneOffset += cav_res->detune; 
        cav_res->wr -= cav_res->wrOffset;
	cav_res->detune = atan(cav_res->Qload * (cav_res->wr / cav_res->wref - cav_res->wref / cav_res->wr));

	//calc detuneOffset =  Setdetune - optDetune
	cav_res->detuneOffset = cav_res->detune - atan(- Vbr/cav_res->VcAP[0] * sin(cav_res->VcAP[1]));//2020/12/07 NY for scan
      } else if (cav_res->detuneOffset != 0.0)// calculation for wrOffset
      {
        cav_res->detune += cav_res->detuneOffset;
        if(fabs(cav_res->detune) > 0.5 * M_PI)
          fprintf(stderr, "Detune angle is out of range (-pi/2 < psi < pi/2 ).\nPlease use the detuneOffsetHz, if you need large frequency offset.\n");
        cav_res->wrOffset = -cav_res->wr;
        cav_res->wr = 0.5 * cav_res->wref *(tan(cav_res->detune) / cav_res->Qload 
          +sqrt( tan(cav_res->detune) / cav_res->Qload * tan(cav_res->detune) / cav_res->Qload + 4 ));
        cav_res->wrOffset += cav_res->wr;
  
        cav_res->detune -= cav_res->detuneOffset;
        cav_res->wr = 0.5 * cav_res->wref *(tan(cav_res->detune) / cav_res->Qload 
         +sqrt( tan(cav_res->detune) / cav_res->Qload * tan(cav_res->detune) / cav_res->Qload + 4 ));
      }
      //////////////////////////////////////////////////////////////
      // Vb calculation, assuming uniform fill
      Vb[0] = Vbr*cos(cav_res->detune)*cos(cav_res->detune + M_PI - cav_res->FourierPhase);// real(Vb) [V]
      Vb[1] = Vbr*cos(cav_res->detune)*sin(cav_res->detune + M_PI - cav_res->FourierPhase);// imag(Vb) [V]
      // initial vg calculation from target Vc,
      optVg[0] = 	cav_res->Vc[0] - Vb[0];		
      optVg[1] = 	cav_res->Vc[1] - Vb[1];
      if(cav_res->geneVolt > 0 && cav_res->FBparaIniSet==0) // not passive mode
      {
	cav_res->geneVolt = sqrt(optVg[0] * optVg[0] + optVg[1] * optVg[1]);		
	cav_res->geneVoltRes[0] = fabs(cav_res->geneVolt / cos (cav_res->detune));		
	cav_res->genePhase = atan2(optVg[1],  optVg[0]);
	if(fabs(cav_res->VcAP[1]) > 0.5 * M_PI)
	  cav_res->genePhase =  - 2 * cav_res->VcAP[1] + cav_res->genePhase - 2*M_PI;
      }
    } else if(fabs((int)(cav_res->m*2)-2*cav_res->m) > 1e-20) //assuming HOM
    {// m is neither integer nor half-integer (m has a information for the detune angle)
      cav_res->wr = cav_res->m * 2.0 * M_PI * ring->frf * FMEGA;
      cav_res->detune = atan(cav_res->Qload * (cav_res->wr / cav_res->wref - cav_res->wref / cav_res->wr));
      if(cav_res->EnableCavInst!=1)
      {
        cav_res->EnableCavInst = 1;
        cav_res->wref = 0;// calculate in wr frame, not reference frame
        cav_res->detune = 0;
      }
    }
    // calculate in reference frame, not wr frame //210303, 210405
    //if(cav_res->VcAP[0] == 0) cav_res->wref = 0;
    cav_res->delOmega = cav_res->wr - cav_res->wref;

    // for cavity compensation. assuming the use of fundamental cavity
    for(k=0;k<cav_res->geneModfreqMrev;k++) 
      cav_res->wMod[k] = 2.0 * M_PI * ring->frf * FMEGA / ring->h * (k+1); 

    ////////////////////////////////////////////////////////////////////////////////////////////
    // calculation for cavity response parameters, filling rate, loss factor & pre-calculation number
    cav_res->fillrate = cav_res->wr * 0.5 / cav_res->Qload; 
    cav_res->lossfactor = cav_res->wr * 0.5 * cav_res->Rs / cav_res->Qzero;

    cav_res->Nturn = (unsigned) (log(2) * 10 / cav_res->fillrate / ring->T0) + 1;
    cav_res->Nbu = cav_res->Nturn * ring->h;
    if(cav_res->Nbu < 2)
    {
      fprintf(stderr, "Given longrange resonator is shortrange!\n");
      return false;
    }
    if (cav_res->Nbu > Ntmp) Ntmp = cav_res->Nbu;
    ////////////////////////////////////////////////////////////////////////////////////////////
  }
  ring->phai0   = acos( ring->q );// cosine difine for cavity resonator
  ring->Nbumax = Ntmp;
  ring->lr_order = 6;
  
  return true;
}

void
calc_resonatorVb(int frag,ring_t * ring, double * fnp_ring, const selffield_model_t * SelfFieldModel, e_beam_t * ebeam, double bunchNp,CAVITY_resonator_t * cav_res, double detune, double *vbout)
	// CAVITY_resonator_t * cav_res, double detune 
	//  frag = 0 -> reference frame
{
  int i, j, k, m;
  double  C[2];
  double progress[2],progressTb[2],progressEb[2],V_old[2]={0.0,0.0}, V_new[2]={0.0,0.0},avgVb[2]={0.0,0.0};
   
  // Determines the phasor after all bunches have Nturn-times passed // 11/07/2018 Naoto Yamamoto
    //CAVITY_resonator_t * cav_res = &(ring->cavity_resonators[l]);   
  int Nbin = SelfFieldModel->Ncell;
  int Nturns = (int)(ring->Nbumax / ring->h); // 10 damping times
  double dTau = SelfFieldModel->dT * SelfFieldModel->sigma_tau; // Bin-width
  double tbucket = 1.0/ring->frf/FMEGA; // = ring->T0 / ring->h distance between two bucket
  double tbunch = (tbucket - Nbin * dTau); // distance between two bunches
  double Vb[2],absVb,Vbr=0;
  double dtmp=- tbunch - dTau * (Nbin * 0.5 );//- 0.5);// interval to synchronus time from 

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  double wr, delOmega,fillrate,lossfactor,fac;
  //wr = 0.5 * cav_res->wref * ( tan(detune) / cav_res->Qload +sqrt( tan(detune) / cav_res->Qload * tan(detune) / cav_res->Qload + 4 ));
  if(frag == 1) wr =cav_res->wr;
  else wr = 0.5 * cav_res->wref * ( tan(detune) / cav_res->Qload +sqrt( tan(detune) / cav_res->Qload * tan(detune) / cav_res->Qload + 4 ));
  delOmega = wr - cav_res->wref*(1-frag);
  fillrate = wr * 0.5 / cav_res->Qload;// wr relates to filling time & lossfactor.
  lossfactor = wr * 0.5 * cav_res->Rs / cav_res->Qzero;
  fac = 2.0 * lossfactor * ring->T0 / ring->E0 / FGIGA / bunchNp * FMILLI;
  ///////////////////////////////////////////////////////////////////////////////////////////
	  
  C[0] = -fillrate;
  C[1] = delOmega;
  // Decay and oscillation of the phasors during the passage of 1 bin with length dTau
  progress[0] = exp(C[0] * dTau) * cos(C[1] * dTau);
  progress[1] = exp(C[0] * dTau) * sin(C[1] * dTau);
  progressTb[0] = exp(C[0] * tbunch) * cos(C[1] * tbunch); // Decay and rotation between one bucket
  progressTb[1] = exp(C[0] * tbunch) * sin(C[1] * tbunch);
  progressEb[0] = exp(C[0] * tbucket) * cos(C[1] * tbucket);// Decay and rotation of phasor during an empty bucket
  progressEb[1] = exp(C[0] * tbucket) * sin(C[1] * tbucket);      

  // ||||||||||||||||||____________||||||||||||||||||____________
  // <- dTau * Nbin  -><- tbunch  ->      |: bin
  // <======    tbucket     =======>      _: no bin 
  //          <------ dtmp -------->

  for(k=0; k<Nturns; k++)
  {
    for(m=0; m<ring->h; m++)
    {
      i = ring->h - m - 1;
      if(ebeam->nfFill[i] > 0.0)
      {
        for(j=0; j<Nbin; j++)
        {        
          V_new[0] = V_old[0] * progress[0] - V_old[1] * progress[1] - ring->Ibunch[i] * fnp_ring[i*Nbin + j] * fac;
          V_new[1] = (V_old[0] * progress[1] + V_old[1] * progress[0]);
          V_old[0] = V_new[0];
          V_old[1] = V_new[1]; 
        }   // Decay and oscillation of the phasors inbetween two bunches
        V_new[0] = V_old[0] * progressTb[0] - V_old[1] * progressTb[1];
        V_new[1] = V_old[0] * progressTb[1] + V_old[1] * progressTb[0];
        V_old[0] = V_new[0];
        V_old[1] = V_new[1];
      } else
      { // Decay and rotation of phasor during the passage of an empty buncket
        V_new[0] = V_old[0] * progressEb[0] - V_old[1] * progressEb[1];
        V_new[1] = V_old[0] * progressEb[1] + V_old[1] * progressEb[0];
        V_old[0] = V_new[0];
        V_old[1] = V_new[1];
      }
      //if(cav_res->VcAP[0] > 0 && k == Nturns*3-1) // if Vc is given and last turn 
      if(k == Nturns-1) // if Vc is given and last turn 
      { 
        Vb[0] = V_old[0] * exp(C[0]*dtmp)*cos(C[1]*dtmp) - V_old[1] * exp(C[0]*dtmp)*sin(C[1]*dtmp);
        Vb[1] = V_old[0] * exp(C[0]*dtmp)*sin(C[1]*dtmp) + V_old[1] * exp(C[0]*dtmp)*cos(C[1]*dtmp);
        avgVb[0] += Vb[0]/ ring->h;	
        avgVb[1] += Vb[1]/ ring->h;	
      }
    }
  }

  absVb = sqrt(avgVb[0] * avgVb[0] + avgVb[1] * avgVb[1]);
  Vbr = absVb /cos(atan(avgVb[1]/avgVb[0]));// Vbr (incl. fill pattern in ref. frame)
  vbout[0]=avgVb[0];
  vbout[1]=avgVb[1];
  if (frag == 1){
    vbout[0]=V_new[0];
    vbout[1]=V_new[1];
  }
  //printf("\n%e %e %e %e\n",wr,lossfactor,avgVb[0],avgVb[1]);
  return;
}

void
calc_resonatorVg(CAVITY_resonator_t * cav_res, double *calcVb)
{
    //CAVITY_resonator_t * cav_res = &(ring->cavity_resonators[l]);   
  double optVg[2]; 

  optVg[0] =  cav_res->Vc[0] - calcVb[0];// + cav_res->directRF_feedback.Vg[0];
  optVg[1] =  cav_res->Vc[1] - calcVb[1];// + cav_res->directRF_feedback.Vg[1];
  cav_res->genePhase = atan2(optVg[1],optVg[0]);

  cav_res->geneVolt = sqrt(optVg[0] * optVg[0] + optVg[1] * optVg[1]);
  cav_res->geneVoltRes[0] = fabs(cav_res->geneVolt / cos (cav_res->detune));		
  return;
}

void
set_resonatorDetune(CAVITY_resonator_t * cav_res, double detune)
{
  // optimum conditon only
  cav_res->detune = detune;

  cav_res->wr = 0.5 * cav_res->wref * ( tan(cav_res->detune) / cav_res->Qload
            +sqrt( tan(cav_res->detune) / cav_res->Qload * tan(cav_res->detune) / cav_res->Qload + 4 ));
  cav_res->delOmega = cav_res->wr - cav_res->wref;
  cav_res->fillrate = cav_res->wr * 0.5 / cav_res->Qload;// wr relates to filling time & lossfactor.
  cav_res->lossfactor = cav_res->wr * 0.5 * cav_res->Rs / cav_res->Qzero;

  // if directRF feedback exists, Vg cannot be calculated from only beam voltage.
  if(cav_res->directRF_feedback.Vg[0] == 0){
    cav_res->geneVoltRes[0] = fabs(cav_res->geneVolt / cos (cav_res->detune));		
  }  

  return;
}


void
set_resonatorDetuneOffset(CAVITY_resonator_t * cav_res)
{
  cav_res->detune += cav_res->detuneOffset;
  if(fabs(cav_res->detune) > 0.5 * M_PI)
    fprintf(stderr, "Detune angle is out of range (-pi/2 < psi < pi/2 ).\nPlease use the detuneOffsetHz, if you need large frequency offset.\n");

  set_resonatorDetune(cav_res, cav_res->detune);

  return;
}


void
BKwake_phasor_Recalc(const ring_t * ring, CAVITY_resonator_t * cav_res, int l,const double * fnp_ring, const selffield_model_t * SelfFieldModel, 
		  weak_bunch_t * bunch, int kb, double * phasor_end, e_beam_t * ebeam)
{
  int i, j, k, m;
  double  C[2]; 
  double progress[2],V_old[2], V_new[2]={0.0,0.0};
   
  int Nbin = SelfFieldModel->Ncell;
  int Nturns = (int)(ring->Nbumax / ring->h); // 10 damping times
  int ii,numCalcVbr = 1;
  double dTau = SelfFieldModel->dT * SelfFieldModel->sigma_tau; // Bin-width
  double tbucket = 1.0/ring->frf/FMEGA; // = ring->T0 / ring->h distance between two bucket
  double tbunch = (tbucket - Nbin * dTau); // distance between two bunches
  double fac;
  double avgVb[2],Vb[2],optVg[2],absVb,Vbr=0;
  double tmpDetune,dtmp;// time of bin center from first bin 
  V_old[0] = V_old[1] = 0.0;
  avgVb[0] = avgVb[1] = 0.0; // 1-turn average Vb 
  if(cav_res->VcAP[0] > 0) {
    if(cav_res->FBparaIniSet == 0) numCalcVbr = 1; 
  }
  for(ii=0;ii<numCalcVbr;ii++)
  {
    dtmp=- tbunch - dTau * (Nbin/2);
    fac = 2 * cav_res->lossfactor * ring->T0 / ring->E0 / FGIGA / bunch->Np * FMILLI;
    C[0] = -cav_res->fillrate;
    C[1] = cav_res->delOmega;
    // Decay and oscillation of the phasors during the passage of 1 bin with length dTau
    progress[0] = exp(C[0] * dTau) * cos(C[1] * dTau);
    progress[1] = exp(C[0] * dTau) * sin(C[1] * dTau);
     
    for(k=0; k<Nturns*3; k++)
    {
      for(m=0; m<ring->h; m++)
      {
        i = ring->h - m - 1;
        if(ebeam->nfFill[i] > 0.0)
        {
          for(j=0; j<Nbin; j++)      
          {        
            //V_new[0] = (V_old[0] * progress[0] - V_old[1] * progress[1]) - ring->Ibunch[i] * fnp_ring[i*Nbin + j] * fac;
	    V_new[0] = (V_old[0] - ring->Ibunch[i] * fnp_ring[i*Nbin + j] * fac) * progress[0] - V_old[1] * progress[1];
            V_new[1] = (V_old[0] * progress[1] + V_old[1] * progress[0]);

            V_old[0] = V_new[0];
            V_old[1] = V_new[1]; 
          }
          // Decay and oscillation of the phasors inbetween two bunches
          V_new[0] = V_old[0] * exp(C[0]*tbunch)*cos(C[1]*tbunch) - V_old[1] * exp(C[0]*tbunch)*sin(C[1]*tbunch);
          V_new[1] = V_old[0] * exp(C[0]*tbunch)*sin(C[1]*tbunch) + V_old[1] * exp(C[0]*tbunch)*cos(C[1]*tbunch);
           
          V_old[0] = V_new[0];
          V_old[1] = V_new[1];
        }   
        else
        { // Decay and rotation of phasor during the passage of an empty buncket
          V_new[0] = V_old[0] * exp(C[0] * tbucket)*cos(C[1] * tbucket) - V_old[1] * exp(C[0] * tbucket)*sin(C[1] * tbucket);
          V_new[1] = V_old[0] * exp(C[0] * tbucket)*sin(C[1] * tbucket) + V_old[1] * exp(C[0] * tbucket)*cos(C[1] * tbucket);
       	  V_old[0] = V_new[0];
          V_old[1] = V_new[1];
        }
        if(cav_res->VcAP[0] > 0 && k == Nturns*3-1) // if Vc is given and last turn 
        {
          Vb[0] = V_old[0] * exp(C[0]*dtmp)*cos(C[1]*dtmp) - V_old[1] * exp(C[0]*dtmp)*sin(C[1]*dtmp);
          Vb[1] = V_old[0] * exp(C[0]*dtmp)*sin(C[1]*dtmp) + V_old[1] * exp(C[0]*dtmp)*cos(C[1]*dtmp);
          avgVb[0] += Vb[0]/ ring->h;	
          avgVb[1] += Vb[1]/ ring->h;	
        }
      }
    }

    if(cav_res->VcAP[0] > 0 && cav_res->FBparaIniSet==0)
    {
      dtmp=0;
      absVb = sqrt(avgVb[0] * avgVb[0] + avgVb[1] * avgVb[1]);
      Vbr = absVb /cos(atan(avgVb[1]/avgVb[0]));// Vbr (incl. fill pattern)
      if(cav_res->geneVolt > 1e-20) //active cavity
      {//Vb -> Vg // Re-calculate Vg parameter considering with the transient beam loading
        optVg[0] =  cav_res->Vc[0] - avgVb[0];
        optVg[1] =  cav_res->Vc[1] - avgVb[1];
        dtmp = cav_res->genePhase;
        cav_res->geneVolt = sqrt(optVg[0] * optVg[0] + optVg[1] * optVg[1]);
	cav_res->geneVoltRes[0] = fabs(cav_res->geneVolt / cos (cav_res->detune));		
        cav_res->genePhase = atan2(optVg[1],optVg[0]);
       
        dtmp -=  cav_res->genePhase;// rad
      } else if(cav_res->geneVolt ==0) //passive cavity
      {
	dtmp = cav_res->wr; // Hz rad
        tmpDetune = atan(- Vbr/cav_res->VcAP[0] * sin(cav_res->VcAP[1]));
        cav_res->wr = 0.5 * cav_res->wref * ( tan(tmpDetune) / cav_res->Qload
             +sqrt( tan(tmpDetune) / cav_res->Qload * tan(tmpDetune) / cav_res->Qload + 4 ));
	dtmp -= cav_res->wr;
            
        //cav_res->wr -= dtmp;// wr should be also tuned. Hz rad
        //cav_res->detune -= dtmp; // rad
        cav_res->delOmega -= dtmp; // Hz rad
        cav_res->fillrate = cav_res->wr * 0.5 / cav_res->Qload;// wr relates to filling time & lossfactor.
        cav_res->lossfactor = cav_res->wr * 0.5 * cav_res->Rs / cav_res->Qzero;
        cav_res->detune = atan(cav_res->Qload * (cav_res->wr / cav_res->wref - cav_res->wref / cav_res->wr));
      }
      if(cav_res->directRF_feedback.gain > 0)
      {
	cav_res->directRF_feedback.Vg[0]=0;
	cav_res->directRF_feedback.Vg[1]=0;
	cav_res->directRF_feedback.fbAmp=0;
	cav_res->directRF_feedback.fbPhase=0;
      }

      if(ring->cav_resonators_cavFB_outBunch == kb)
      {
        fprintf(ring->cavity_resonators_cavFB_fp,"\n#Pre %d: Vbr = %.5e [V], angVb = %.6lf [rad], psi = %.6lf, Vg(Amp,Phase) = %.5e %.6lf",
          l,Vbr*(ring->E0 * FGIGA),atan(avgVb[1]/avgVb[0]),cav_res->detune,cav_res->geneVolt,cav_res->genePhase);	
        fprintf(ring->cavity_resonators_cavFB_fp,"\n#Pre %d: delta f  = %.6lf [Hz],Vb: %.10e %.10e",
          l, cav_res->delOmega/2/ M_PI,avgVb[0], avgVb[1]);
      }
      avgVb[0]=0;// 1-turn average Vb 
      avgVb[1]=0;
    }
  } 
  phasor_end[l*2] = V_new[0];// for tracking, the resulting (complex) voltage is succeeded to next turn.
  phasor_end[l*2 + 1] = V_new[1]; 
}

// not test
void
wake_phasor_Recalc(ring_t * ring, CAVITY_resonator_t * cav_res, int l,double * fnp_ring, const selffield_model_t * SelfFieldModel, 
		  weak_bunch_t * bunch, int kb, double * phasor_end, e_beam_t * ebeam)
{
  double calvb[2],optDetune, foutFrag = 0;

  if(cav_res->VcAP[0] > 0 )
  {
    if(cav_res->FBparaIniSet==0) //**** Cavity parameter Search & Set ***
    {
     //// 1. VBR => optimum detune
      calc_resonatorVb(0,ring, fnp_ring, SelfFieldModel,ebeam, bunch->Np, cav_res, 0, calvb); 
      optDetune = atan(- sqrt(calvb[0]*calvb[0]+calvb[1]*calvb[1])/cav_res->VcAP[0] * sin(cav_res->VcAP[1]));
      set_resonatorDetune(cav_res,optDetune); // genePhase = 0 for passive cavity.
      if(ring->cav_resonators_cavFB_outBunch == kb) 
        fprintf(ring->cavity_resonators_cavFB_fp,
       "## seq.1 ## : Search Optimum detuning from Beam loading voltage at resonance.\n# Re[Vb] = %e, Im[Vb] = %e, optDetune = %.5lf\n",
       calvb[0],calvb[1],optDetune); 
    
      //// 2. Optimum detune => Vg //   ==> (Comment out)Operating detune
      calc_resonatorVb(0,ring, fnp_ring, SelfFieldModel,ebeam, bunch->Np, cav_res, optDetune, calvb);
      if(cav_res->geneVolt > 1e-20) calc_resonatorVg(cav_res, calvb);
      if(ring->cav_resonators_cavFB_outBunch == kb) 
	fprintf(ring->cavity_resonators_cavFB_fp,
	"## seq.2 ## : Search generator voltage with obtained operating detuning.\n# Re[Vb] = %e, Im[Vb] = %e, PhaseVg = %.5lf, Detune = %.5lf\n",
	calvb[0],calvb[1],cav_res->genePhase,cav_res->detune); 

      ///// 3. Set detuning offset => Vg
      set_resonatorDetuneOffset(cav_res);
      calc_resonatorVb(0,ring, fnp_ring, SelfFieldModel,ebeam, bunch->Np, cav_res, cav_res->detune, calvb);
      if(cav_res->geneVolt > 1e-20) calc_resonatorVg(cav_res, calvb);
      if(ring->cav_resonators_cavFB_outBunch == kb) 
	fprintf(ring->cavity_resonators_cavFB_fp,
	"## seq.3 ## : Adding detuning offset.\n# Re[Vb] = %e, Im[Vb] = %e, PhaseVg = %.5lf, Detune = %.5lf\n",
	calvb[0],calvb[1],cav_res->genePhase,cav_res->detune); 
     
      if(cav_res->directRF_feedback.gain > 0) //**** Direct RF paramter Searcg & Set ***
      {
        if(cav_res->FBparaIniSet==1) calc_resonatorVb(0,ring, fnp_ring, SelfFieldModel,ebeam, bunch->Np, cav_res, cav_res->detune, calvb); 
        // direct RF feedback set up
        if(ring->cav_resonators_cavFB_outBunch == kb)
        {
          fprintf(ring->cavity_resonators_cavFB_fp,"## seq.3.5 ## : Calculating Direct RF parameter.\n");
          foutFrag = 1;
        }
        init_directRF_feedback(cav_res,calvb, ring->cavity_resonators_drfFB_fp,foutFrag,l);   
      }

     cav_res->delOmega = cav_res->wr;////////////
     if(ring->cav_resonators_cavFB_outBunch == kb) 
     fprintf(ring->cavity_resonators_cavFB_fp,"## seq.4 ## : Move to Real frame.\n"); 
    }
  }
  ///// 4. Calc pot voltage 
  calc_resonatorVb(1,ring, fnp_ring, SelfFieldModel,ebeam, bunch->Np, cav_res, cav_res->detune, calvb);
  if(ring->cav_resonators_cavFB_outBunch == kb) 
    fprintf(ring->cavity_resonators_cavFB_fp,"# Re[Vb] = %e, Im[Vb] = %e, AmpVg = %e, PhaseVg = %.5lf, Detune = %.5lf\n",
      //"## seq.4 ## : Move to Real frame.\n# Re[Vb] = %e, Im[Vb] = %e, AmpVg = %e, PhaseVg = %.5lf, Detune = %.5lf\n",
      calvb[0],calvb[1],cav_res->geneVolt,cav_res->genePhase,cav_res->detune); 
	
  phasor_end[l*2] = calvb[0];// for tracking, the resulting (complex) voltage is succeeded to next turn.
  phasor_end[l*2 + 1] = calvb[1]; 

  init_cavVol_feedback(&cav_res->cavVolFB,cav_res->VcAP); 

}
 

void
wake_phasor_initTR(ring_t * ring, double * fnp_ring, const selffield_model_t * SelfFieldModel, 
		  weak_bunch_t * bunch, int kb, double * phasor_end, e_beam_t * ebeam)
{
  int l;
   
  if(ring->cav_resonators_cavFB_outBunch == kb)
  { // for feedback
    if (ring->cavity_resonators_cavFB_fp == NULL) {
      char filename_fb[FILENAME_MAX] = "";
      snprintf(filename_fb, FILENAME_MAX, "%s/cavfeedback_bunch.dat", track.work_path);  
      ring->cavity_resonators_cavFB_fp = fopen(filename_fb, "w+");
      if(ring->cavity_resonators_cavFB_fp == NULL) ERROR("fopen_feedback_bunch", return);
    }
    fprintf(ring->cavity_resonators_cavFB_fp, "##### CAVfeedback, Ib = %g A\n", bunch->Ib);
  }
  if (ring->cav_resonators_drfFB_active == 1 && ring->cav_resonators_cavFB_outBunch == kb) //bunch->kb_out == 1)
  {
    if (ring->cavity_resonators_drfFB_fp == NULL) {
      char filename_dfb[FILENAME_MAX] = "";
      snprintf(filename_dfb, FILENAME_MAX, "%s/directRFfeedback_bunch.dat", track.work_path);  
      ring->cavity_resonators_drfFB_fp = fopen(filename_dfb, "w+");
      if(ring->cavity_resonators_drfFB_fp == NULL) ERROR("fopen_directRFfeedback_bunch", return);
      fprintf(ring->cavity_resonators_drfFB_fp, "# RFfeedback, Ib = %g A ;\n", bunch->Ib);
    }
  }

  // Calculate the phasor after all bunches have Nturn-times passed // 11/07/2018 Naoto Yamamoto
  for(l = 0; l < ring->cavity_resonators_size; l++)
  {   
    CAVITY_resonator_t * cav_res = &(ring->cavity_resonators[l]);   
    double calvb[2],optDetune, foutFrag = 0;

    // output cavity feedback 
    if(ring->cav_resonators_cavFB_outBunch == kb)
    {
      if(cav_res->VcAP[0] > 0)
	fprintf(ring->cavity_resonators_cavFB_fp,"### %d: FBmode = %d , targetVc: %.6e %.6e\n",l,cav_res->mode,cav_res->VcAP[0],cav_res->VcAP[1]);
      else
        fprintf(ring->cavity_resonators_cavFB_fp,"### %d: FBmode = %d\n",l,cav_res->mode);
    }

    if(cav_res->VcAP[0] > 0 )
    {
      if(cav_res->FBparaIniSet==0) //**** Cavity parameter Search & Set ***
      {
      //// 1. VBR => optimum detune
      calc_resonatorVb(0,ring, fnp_ring, SelfFieldModel,ebeam, bunch->Np, cav_res, 0, calvb); 
      optDetune = atan(- sqrt(calvb[0]*calvb[0]+calvb[1]*calvb[1])/cav_res->VcAP[0] * sin(cav_res->VcAP[1]));
      set_resonatorDetune(cav_res,optDetune); // genePhase = 0 for passive cavity.
      if(ring->cav_resonators_cavFB_outBunch == kb) 
        fprintf(ring->cavity_resonators_cavFB_fp,
       "## seq.1 ## : Search Optimum detuning from Beam loading voltage at resonance.\n# Re[Vb] = %e, Im[Vb] = %e, optDetune = %.5lf\n",
       calvb[0],calvb[1],optDetune); 
    
      //// 2. Optimum detune => Vg //   ==> (Comment out)Operating detune
      calc_resonatorVb(0,ring, fnp_ring, SelfFieldModel,ebeam, bunch->Np, cav_res, optDetune, calvb);
      if(cav_res->geneVolt > 1e-20) calc_resonatorVg(cav_res, calvb);
      if(ring->cav_resonators_cavFB_outBunch == kb) 
	fprintf(ring->cavity_resonators_cavFB_fp,
	"## seq.2 ## : Search generator voltage with obtained operating detuning.\n# Re[Vb] = %e, Im[Vb] = %e, PhaseVg = %.5lf, Detune = %.5lf\n",
	calvb[0],calvb[1],cav_res->genePhase,cav_res->detune); 

      ///// 3. Set detuning offset => Vg
      set_resonatorDetuneOffset(cav_res);
      calc_resonatorVb(0,ring, fnp_ring, SelfFieldModel,ebeam, bunch->Np, cav_res, cav_res->detune, calvb);
      if(cav_res->geneVolt > 1e-20) calc_resonatorVg(cav_res, calvb);
      if(ring->cav_resonators_cavFB_outBunch == kb) 
	fprintf(ring->cavity_resonators_cavFB_fp,
	"## seq.3 ## : Adding detuning offset.\n# Re[Vb] = %e, Im[Vb] = %e, PhaseVg = %.5lf, Detune = %.5lf\n",
	calvb[0],calvb[1],cav_res->genePhase,cav_res->detune); 
      }

      if(cav_res->directRF_feedback.gain > 0) //**** Direct RF paramter Searcg & Set ***
      {
        if(cav_res->FBparaIniSet==1) calc_resonatorVb(0,ring, fnp_ring, SelfFieldModel,ebeam, bunch->Np, cav_res, cav_res->detune, calvb); 

        // direct RF feedback set up
        if(ring->cav_resonators_cavFB_outBunch == kb){
          fprintf(ring->cavity_resonators_cavFB_fp,"## seq.3.5 ## : Calculating Direct RF parameter.\n");
          foutFrag = 1;
        }
        init_directRF_feedback(cav_res,calvb, ring->cavity_resonators_drfFB_fp,foutFrag,l);   
      }

      cav_res->delOmega = cav_res->wr;////////////
      if(ring->cav_resonators_cavFB_outBunch == kb) 
      fprintf(ring->cavity_resonators_cavFB_fp,"## seq.4 ## : Move to Real frame.\n"); 
    }

    ///// 4. Calc pot voltage 
    calc_resonatorVb(1,ring, fnp_ring, SelfFieldModel,ebeam, bunch->Np, cav_res, cav_res->detune, calvb);
    if(ring->cav_resonators_cavFB_outBunch == kb) 
      fprintf(ring->cavity_resonators_cavFB_fp,"# Re[Vb] = %e, Im[Vb] = %e, AmpVg = %e, PhaseVg = %.5lf, Detune = %.5lf\n",
      //"## seq.4 ## : Move to Real frame.\n# Re[Vb] = %e, Im[Vb] = %e, AmpVg = %e, PhaseVg = %.5lf, Detune = %.5lf\n",
      calvb[0],calvb[1],cav_res->geneVolt,cav_res->genePhase,cav_res->detune); 

	
    phasor_end[l*2] = calvb[0];// for tracking, the resulting (complex) voltage is succeeded to next turn.
    phasor_end[l*2 + 1] = calvb[1]; 

    init_cavVol_feedback(&cav_res->cavVolFB,cav_res->VcAP); 
    cav_res->geneVoltRes[1]= cav_res->genePhase - cav_res->detune;
    Vg2Ig(cav_res);
    if(ring->cav_resonators_cavFB_outBunch == kb) 
      fprintf(ring->cavity_resonators_cavFB_fp,"# Generator current; Amp = %.5e, Phase = %.5lf (Vgr phase = %.5lf)\n",
		      cav_res->IgAP[0],cav_res->IgAP[1],cav_res->geneVoltRes[1]);
    Ig2Vg(cav_res,ring);
    if(ring->cav_resonators_cavFB_outBunch == kb) 
      fprintf(ring->cavity_resonators_cavFB_fp,"# Generator voltage; Amp = %.5e, Phase = %.5lf,(Vgr phase = %.5lf)\n",
		      cav_res->geneVolt,cav_res->genePhase, cav_res->geneVoltRes[1]);
    cav_res->tunerOffset = cav_res->VcAP[1]-cav_res->geneVoltRes[1];
    if(ring->cav_resonators_cavFB_outBunch == kb) 
      fprintf(ring->cavity_resonators_cavFB_fp,"# TunerOffset = %.5lf, IntegralMemoryFac = %5e\n",cav_res->tunerOffset,cav_res->cavVolFB.integralMemoryFac[0]);
  }
  if(ring->cav_resonators_cavFB_outBunch == kb )
    fprintf(ring->cavity_resonators_cavFB_fp,"\n# Vg(Amp,Phase), detune wr, Vc(Amp,Phase) : ratio(Amp,Phase,tuner)");
}
 

 /**
  * 
  * @param lr_wake Stores wake function of all bins from several longrange resonators
  * @param phasor_end Stores phasors at the end of the last turn
  */
 void
 construct_wake_phasorTR(double * lr_wake, double * phasor_end, const selffield_model_t * SelfFieldModel,
                       const ring_t * ring, int kb, const double * fnp,e_beam_t * ebeam, weak_bunch_t * bunch)
 {
  int i, j, l, m,ii,scan_num;
  int foutFrag = 0;
  double progress[2],progressTb[2],progressEb[2],progressDt[2],V_old[2],V_new[2];
  double C[2];
  //double prog2beam[2],progb2beam[2];
   
  static unsigned long int turn=1;
  if(turn*ring->h >= ULONG_MAX) turn = 1;

  // Calculation the phasor voltage after all bunches have passed and stores // 2018/07/16. Naoto Yamamoto
  for(l = 0; l < ring->cavity_resonators_size; l++)
  {     
    CAVITY_resonator_t * cav_res = &(ring->cavity_resonators[l]);
    if(cav_res->switchOnTurn > turn ) continue;
    else if(turn > 1 && cav_res->switchOnTurn == turn && cav_res->mode >0 )  wake_phasor_Recalc(ring, cav_res, l, fnp, SelfFieldModel, bunch, kb, phasor_end, ebeam);
    if(cav_res->switchOffTurn <= turn ) continue;
    int Nbin = SelfFieldModel->Ncell;   
    double cavField,tauP,ptmp,tmpVc[2],dtmp,ttmp,avgVb[2],Vb[2];
    double dTau = SelfFieldModel->dT * SelfFieldModel->sigma_tau; // Bin-width [s]
    double tbucket = 1.0/ring->frf/FMEGA; // distance between two bucket
    double tbunch = (tbucket - Nbin * dTau); // distance between two bunches
    double fac = 2.0 * cav_res->lossfactor * ring->T0 / ring->E0 / FGIGA / bunch->Np * FMILLI;
    C[0] = -cav_res->fillrate;
    C[1] = cav_res->wr;//delOmega;// [Hz rad] 
    dtmp=- tbunch - dTau * (Nbin/2);
    ttmp=0;

    progress[0] = exp(C[0] * dTau) * cos(C[1] * dTau); // Decay and rotation of phasor during one bin
    progress[1] = exp(C[0] * dTau) * sin(C[1] * dTau);
    progressTb[0] = exp(C[0] * tbunch) * cos(C[1] * tbunch); // Decay and rotation between one bucket
    progressTb[1] = exp(C[0] * tbunch) * sin(C[1] * tbunch);
    progressEb[0] = exp(C[0] * tbucket) * cos(C[1] * tbucket);// Decay and rotation of phasor during an empty bucket
    progressEb[1] = exp(C[0] * tbucket) * sin(C[1] * tbucket);      
    progressDt[0] = exp(C[0] * dtmp) * cos(C[1] * dtmp); // Decay and rotation between one bucket
    progressDt[1] = exp(C[0] * dtmp) * sin(C[1] * dtmp);
    V_old[0] = phasor_end[l*2];// real voltage just before the passage of the target bunch. 
    V_old[1] = phasor_end[l*2 + 1];
    V_new[0] = phasor_end[l*2]; // real voltage just after the passage of the target bunch.
    V_new[1] = phasor_end[l*2 + 1]; 

    avgVb[0]=0;// 1-turn average Vb 
    avgVb[1]=0;
    for(m=0; m<ring->h; m++)
    {
      i = ring->h - m - 1;
      if(ebeam->nfFill[i] > 0.0)
      {
       for(j=0; j<Nbin; j++)
       {
	 V_new[0] = (V_old[0] * progress[0] - V_old[1] * progress[1] - ring->Ibunch[i] * fnp[i*Nbin + j] * fac); 
	 V_new[1] = (V_old[0] * progress[1] + V_old[1] * progress[0]);	 
	 
	 if(i == kb)
	 {// for calculation of interaction with beam, phase ref. to the beam in driving frequency is needed.
	   cavField = 0;
	   tauP = (-SelfFieldModel->Nsigma + SelfFieldModel->dT * (j - 0.5)) * SelfFieldModel->sigma_tau;// time [s] 
	   ptmp = cav_res->m * ring->wrf * ( tauP + tbucket * ( kb + turn*ring->h)); // [rad]
	   //////////////////////////////////////////////////////////////
	   // beam induced  component // 2019/09/12 N.Yamamoto
           cavField +=  V_new[0] + 0.5*ring->Ibunch[i]*fnp[i*Nbin+j]*fac;  //ptmp = 0;
	   ///////////////////////////////////////////////////////////////////////////////////
	   // generator component // 2019/09/12 N.Yamamoto
	   // Vg=geneVolt * exp(i*omega_rf*t)*exp(i*genePhase); t=0+ptmp;
	   //   =geneVolt * exp(i*ptmp)*exp(i*genePhase)
	   // real Vg = geneVolt * cos(genPhase + ptmp)
           cavField += cav_res->geneVolt * cos(cav_res->genePhase + ptmp);
	   // Modulated generator
	   for(ii=0; ii< cav_res->geneModfreqMrev; ii++)
           {// Mrev=1; 1 revolution; 1 cycle for harmonics h; 2*pi/h/Mrev*bucket // 2019/09/16
	     cavField += cav_res->geneModVolt[ii] * cos( 2*M_PI/ring->h*(ii+1)* (kb + cav_res->geneModPhaseOffsetM[ii])+ptmp + cav_res->genePhase);
	   }
	   //////Direct RF Feedback // 2020/04/13 N.Yamamoto
           if( cav_res->directRF_feedback.gain > 0){
             if (turn > cav_res->directRF_feedback.switchOnTurn && turn < cav_res->directRF_feedback.switchOffTurn) 
	         cavField += cav_res->directRF_feedback.fbAmp* cos(ptmp + cav_res->directRF_feedback.fbPhase);
	   }

	   lr_wake[j] += cavField;
         }

	 V_old[0] = V_new[0];
         V_old[1] = V_new[1];
       }
       // Deca and rotation of phasor between bunches
       V_new[0] = V_old[0] * progressTb[0] - V_old[1] * progressTb[1];
       V_new[1] = V_old[0] * progressTb[1] + V_old[1] * progressTb[0];

       V_old[0] = V_new[0];
       V_old[1] = V_new[1]; 
      }
      else  // nfFill[i] < 0.0, empty bucket
      { // Decay and rotation of phasor during the passage of an empty buncket
        V_new[0] = V_old[0] * progressEb[0] - V_old[1] * progressEb[1];
        V_new[1] = V_old[0] * progressEb[1] + V_old[1] * progressEb[0];

        V_old[0] = V_new[0];
        V_old[1] = V_new[1];
      }
      Vb[0] = V_old[0] * progressDt[0] - V_old[1] * progressDt[1];
      Vb[1] = V_old[0] * progressDt[1] + V_old[1] * progressDt[0];
      avgVb[0] += Vb[0]/ ring->h;	
      avgVb[1] += Vb[1]/ ring->h;	


      //high speed FB???
      if(cav_res->VcAP[0] > 0 && cav_res->cavVolFB.PIcontrol > 0 && cav_res->directRF_feedback.gain) 
	      highSpeedFB(cav_res, ring, l,turn, m, kb, Vb);
    }    
    phasor_end[l*2] = V_new[0];
    phasor_end[l*2 + 1] = V_new[1];   

    ///////////////////////////////////////////////////////////////////
    // Field error for Cavity feedback 
    tmpVc[0]=cav_res->geneVolt * cos(cav_res->genePhase)+avgVb[0] + cav_res->directRF_feedback.Vg[0];
    tmpVc[1]=cav_res->geneVolt * sin(cav_res->genePhase)+avgVb[1] + cav_res->directRF_feedback.Vg[1];
    cav_res->VcTrackAP[0] = sqrt(tmpVc[0] * tmpVc[0] + tmpVc[1] * tmpVc[1]);
    cav_res->VcTrackAP[1] = atan2(tmpVc[1],tmpVc[0]);
    put_cavVolFB_value(&cav_res->cavVolFB, cav_res->VcTrackAP);

    if(ring->cav_resonators_cavFB_outBunch == kb && turn%cav_res->fbMon == 0) // for debug
    {
      fprintf(ring->cavity_resonators_cavFB_fp,"\n%d %lu  %.5e %.6lf  %.6lf %.5e  %.4e %.5lf",l,turn,
      cav_res->geneVolt,cav_res->genePhase, cav_res->detune,cav_res->wr*0.5/M_PI,cav_res->VcTrackAP[0],cav_res->VcTrackAP[1]);
    } 
    if(cav_res->VcAP[0] > 0)
    { // Field error for Cavity feedback
      if(cav_res->cavVolFB.PIcontrol==0)
      {
        dtmp=cav_res->VcAP[0] / get_cavVolFB_value(&cav_res->cavVolFB,0);// dtmp > 1.0, Vc -> small
        ptmp=get_cavVolFB_value(&cav_res->cavVolFB,1) - cav_res->VcAP[1];// ptmp < 1.0, Phase of Vc -> small
        //ttmp= get_cavVolFB_value(&cav_res->cavVolFB,1) - cav_res->geneVoltRes[1] - cav_res->detuneOffset;
	ttmp = cav_res->genePhase - cav_res->detune + cav_res->detuneOffset - cav_res->VcAP[1];	
      
        if(fabs(1-dtmp) < cav_res->cavVolFB.limitRange) dtmp = 1.0; // dtmp=1.0 -> FB off; (1.0)^x = 1.0
        if(fabs(1-ptmp/cav_res->VcAP[1]) < cav_res->cavVolFB.limitRange) ptmp = 0.0;
        if(ring->cav_resonators_cavFB_outBunch == kb && turn%cav_res->fbMon == 0) fprintf(ring->cavity_resonators_cavFB_fp,"  %.6lf %.6lf %.6lf",dtmp,ptmp,ttmp); // for debug
        if(turn > cav_res->cavVolFB.switchOnTurn && turn < cav_res->cavVolFB.switchOffTurn) cavity_feedback(cav_res,dtmp,ptmp);// if Vc is given
      } else {
        dtmp= 1- tmpVc[0]/cav_res->Vc[0];// dtmp > 1.0, Vc -> high
        ptmp= 1- tmpVc[1]/cav_res->Vc[1];// dtmp > 1.0, Vc -> high
	//dtmp = sqrt(tmpVc[0]*tmpVc[0]+tmpVc[1]*tmpVc[1])/sqrt(cav_res->Vc[0]*cav_res->Vc[0]+cav_res->Vc[1]*cav_res->Vc[1]);
	//ptmp = atan2(tmpVc[1],tmpVc[0])-atan2(cav_res->Vc[1]/cav_res->Vc[0]);
        dtmp=cav_res->VcAP[0] / get_cavVolFB_value(&cav_res->cavVolFB,0);// dtmp > 1.0, Vc -> small
        ptmp=get_cavVolFB_value(&cav_res->cavVolFB,1) - cav_res->VcAP[1];// ptmp < 1.0, Phase of Vc -> small	
        ttmp= get_cavVolFB_value(&cav_res->cavVolFB,1) - cav_res->geneVoltRes[1] - cav_res->tunerOffset;

        if(turn > cav_res->cavVolFB.switchOnTurn && turn < cav_res->cavVolFB.switchOffTurn) cavity_PIfeedback(cav_res,ring,0,0,ttmp);// if Vc is given
        if(ring->cav_resonators_cavFB_outBunch == kb && turn%cav_res->fbMon == 0) 
		fprintf(ring->cavity_resonators_cavFB_fp,"  %.5e %.5e %.6lf %.5e",dtmp,ptmp,ttmp,cav_res->cavVolFB.gainInte[0]*cav_res->cavVolFB.integralMemory[0]); // for debug
		//fprintf(ring->cavity_resonators_cavFB_fp,"  %.5e %.5e %.6lf %.5e %.5e",dtmp,ptmp,ttmp,cav_res->geneVoltRes[1],cav_res->cavVolFB.gainInte[0]*cav_res->cavVolFB.integralMemory[0]); // for debug
      }

    }
    ///////////////////////////////////////////////////////////////////
    //////Direct RF Feedback // 2020/04/13 N.Yamamoto
    if(cav_res->directRF_feedback.gain > 0) 
    {
      if(ring->cav_resonators_cavFB_outBunch == kb && turn%cav_res->fbMon == 0 ) foutFrag = 1;
      //put_directRF_feedback(&cav_res->directRF_feedback, tmpVc, cav_res->Vc,cav_res->geneVolt,ring->cavity_resonators_drfFB_fp,foutFrag,l,turn);
      put_directRF_feedback(cav_res, ring->cavity_resonators_drfFB_fp,foutFrag,l,turn);
      foutFrag = 0;
    }

    //scan 2021/04/02 N.Yamamoto
    for(scan_num = 0; scan_num < cav_res->scan_size; scan_num++)
    {
      if(turn % cav_res->scan[scan_num].nrev == 0)
      {
	switch (cav_res->scan[scan_num].type)
	{
	  case detuneRAD:
            cav_res->detune += cav_res->scan[scan_num].step;
            cav_res->wr = 0.5 * cav_res->wref *(tan(cav_res->detune) / cav_res->Qload
              +sqrt( tan(cav_res->detune) / cav_res->Qload * tan(cav_res->detune) / cav_res->Qload + 4 ));
            cav_res->delOmega = cav_res->wr; 
	    break;
	  case fr:
	    cav_res->wr += 2.0 * M_PI * cav_res->scan[scan_num].step;
            if(cav_res->wref != 0) cav_res->detune = atan(cav_res->Qload * (cav_res->wr / cav_res->wref - cav_res->wref / cav_res->wr));
	    cav_res->delOmega = cav_res->wr;
	    break;
	}
      }
    }
  }
  turn++;
}

void highSpeedFB(CAVITY_resonator_t * cav_res,const ring_t * ring, int l,long int turn, int bucket,int kb, double * vbmon)
{
  int i;
  double vcmon[2],dtmp,ptmp;
  static double avgVb[2] = {0};

  // averaging vb between sampling period
  for(i=0; i < 2; i++) avgVb[i] += vbmon[i]/cav_res->cavVolFB.sample;

  if( bucket % cav_res->cavVolFB.sample != cav_res->cavVolFB.sample -1) return;

  vcmon[0]=cav_res->geneVolt * cos(cav_res->genePhase)+avgVb[0];//vbmon[0];// + cav_res->directRF_feedback.Vg[0];
  vcmon[1]=cav_res->geneVolt * sin(cav_res->genePhase)+avgVb[1];//vbmon[1];// + cav_res->directRF_feedback.Vg[1];

  if(cav_res->mode < 100) // for AmpPhase FB
  {
    dtmp = 1- vcmon[0]/cav_res->Vc[0];// dtmp > 1.0, Vc -> high
    ptmp = 1- vcmon[1]/cav_res->Vc[1];// dtmp > 1.0, Vc -> high
  } else {
    dtmp = 1- sqrt( vcmon[0] * vcmon[0] + vcmon[1] * vcmon[1]) / cav_res->VcAP[0];
    ptmp = 1- atan2(vcmon[1], vcmon[0]) / cav_res->VcAP[1];
  }

  //Loopdelay  
  cav_res->cavVolFB.loopD[cav_res->cavVolFB.loopIndex] = dtmp;
  cav_res->cavVolFB.loopD[cav_res->cavVolFB.loopIndex + (cav_res->cavVolFB.loopDelay + 1) ] = ptmp;


  i = cav_res->cavVolFB.loopIndex + 1;
  if(i >  cav_res->cavVolFB.loopDelay - 1) i = 0; 


  dtmp = cav_res->cavVolFB.loopD[i] ; 
  ptmp = cav_res->cavVolFB.loopD[i + cav_res->cavVolFB.loopDelay + 1]; 

  //if(kb == 0 && turn < 2) printf("LoopIndex = %d, i = %d, %e %e \n",cav_res->cavVolFB.loopIndex,i,dtmp,cav_res->cavVolFB.loopD[i]);
  if(turn > cav_res->cavVolFB.switchOnTurn && turn < cav_res->cavVolFB.switchOffTurn) cavity_PIfeedback(cav_res,ring,dtmp,ptmp,0);// if Vc is given

  cav_res->cavVolFB.loopIndex++;
  if (cav_res->cavVolFB.loopIndex == cav_res->cavVolFB.loopDelay) cav_res->cavVolFB.loopIndex = 0;
  
  
  for(i=0; i < 2; i++) avgVb[i]=0;

// for debug
//  double turnDouble = turn + (bucket+1)/ring->h -1;
//  double tmpAmp,tmpPhase;
//  tmpAmp = sqrt(vcmon[0] * vcmon[0] + vcmon[1] * vcmon[1]);
//  tmpPhase = atan2(vcmon[1],vcmon[0]);
//  if(ring->cav_resonators_cavFB_outBunch == kb && turn % cav_res->fbMon == 0) 
//    fprintf(ring->cavity_resonators_cavFB_fp,"\n%d %.2lf %.5e %.6lf %.6lf %.5e %.5e %.5e  %.5e %.5e 0 %.5e %.5e",
//      l,turnDouble,cav_res->geneVolt,cav_res->genePhase, cav_res->detune,cav_res->wr*0.5/M_PI,tmpAmp,tmpPhase,
//      dtmp,ptmp,cav_res->cavVolFB.gainProp[0]*dtmp,cav_res->cavVolFB.gainInte[0]*cav_res->cavVolFB.integralMemory[0]); // for debug
  
  
  return;
}

///////////////////////////////////////////////////
// Initialize direct RF feedback,
///////////////////////////////////////////////////
//// for DRFB before tracking
// test 2: PAC09 WE5PFP089 L.H.Chang, boundary :  drf->phase & drf->AmpRatio
void init_directRF_feedback(CAVITY_resonator_t * cav_res,double * vb, FILE *fp, int frag, int cav)
{
  int i;
  double newVg[2];
  directRF_feedback_t * drf = &(cav_res->directRF_feedback);

  //if(drf->loopDelay == 0) drf->loopDelay = 1;
  drf->loopV = NULL;
  drf->loopV = (double*)malloc(sizeof(double) * (drf->loopDelay + 1) * 2);
 
  drf->fbAmp = drf->gain * cav_res->VcAP[0];// * cos(cav_res->detune);
  drf->fbPhase = cav_res->VcAP[1] + drf-> phaseShift;// + cav_res->detune;
  drf->Vg[0] = drf->fbAmp * cos(drf->fbPhase);
  drf->Vg[1] = drf->fbAmp * sin(drf->fbPhase);
  newVg[0] = cav_res->Vc[0] - vb[0] - drf->Vg[0];
  newVg[1] = cav_res->Vc[1] - vb[1] - drf->Vg[1];
  cav_res->geneVolt = sqrt( newVg[0] * newVg[0] + newVg[1] * newVg[1]);
  cav_res->geneVoltRes[0] = fabs(cav_res->geneVolt / cos (cav_res->detune));
  cav_res->genePhase = atan2( newVg[1],newVg[0]);

  drf->loopIndex = 0;
  for(i =0; i < drf->loopDelay + 1 ;i++) {
    drf->loopV[i]               = cav_res->VcAP[0] * cos(cav_res->VcAP[1]);// real !!!
    drf->loopV[i + drf->loopDelay + 1] = cav_res->VcAP[0] * sin(cav_res->VcAP[1]);// imag. part !!!
  } 
  if(frag == 1)
  {
    fprintf(fp,"#%d: Gain = %.2e, phaseShift = %.6lf, loopDelay = %d, Average = %d, SwitchOnTurn = %d, \n## target: %.6e %.6e\n",
           cav,drf->gain,drf->phaseShift,drf->loopDelay,drf->loopDelay,drf->switchOnTurn,cav_res->VcAP[0],cav_res->VcAP[1]);
    fprintf(fp,"## MainVgAmp = %.5e, MainVgPhase = %.5lf \n",cav_res->geneVolt,cav_res->genePhase);
    fprintf(fp,"## drfFBVgAmp = %.5e, drfFBPhase = %.5lf \n",drf->fbAmp,drf->fbPhase);
   // fprintf(fp,"# cavNum turn loopDelay ReVg ImVg ReVc ImVc Vdrf/Vg dB");
    fprintf(fp,"# cavNum turn loopDelay AmpDRF PhaseDRF Vdrf/Vg dB");
  }
  return;
}

// for DRFB during tracking
void put_directRF_feedback(CAVITY_resonator_t * cav_res, FILE *fp,int frag, int cav, long int turn)
{
  int i;
  directRF_feedback_t * drf = &(cav_res->directRF_feedback);

  // loopV = [ AmpVc1 AmpVc2 ,,, PhaseVc1, PhaseVc2, ,,, ]
  for(i=0; i < 2; i++) drf->loopV[drf->loopIndex + (drf->loopDelay + 1) * i ] = cav_res->VcTrackAP[i];

  i = drf->loopIndex - 1;
  if(i < 0) i = drf->loopDelay; 

  drf->loopIndex++;
  if (drf->loopIndex > drf->loopDelay) drf->loopIndex = 0;

  drf->fbAmp = drf->gain * drf->loopV[i] ; //* cos(cav_res->detune);
  drf->fbPhase = drf->loopV[i + drf->loopDelay + 1] + drf-> phaseShift; //cav_res->detune + drf-> phaseShift;
  drf->Vg[0] = drf->fbAmp * cos(drf->fbPhase);
  drf->Vg[1] = drf->fbAmp * sin(drf->fbPhase);

  double fac,dB;
  fac = drf->fbAmp/ cav_res->geneVolt;
  dB = 20*log10(fac);

  if(frag == 1) fprintf(fp,"\n%d %ld %d %.5e %.5e %.5e %.5e", cav,turn,drf->loopIndex,drf->fbAmp,drf->fbPhase,fac,dB);
  //if(frag == 1) fprintf(fp,"\n%d %ld %d %.5e %.5e %.5e %.5e %.5e %.5e", cav,turn,drf->loopIndex,drf->Vg[0],drf->Vg[1],v[0],v[1],fac,dB);
  return;
}

///////////////////////////////////////////////////
// Initialize cavity voltage feedback,
///////////////////////////////////////////////////
void init_cavVol_feedback(cavityRFVoltage_feedback_t * cavFB, double * vAP)
{
  int i;

  if(cavFB->avgNum == 0) cavFB->avgNum = 1;
  cavFB->ampVc = NULL; 
  cavFB->ampVc   = (double*)malloc(sizeof(double) * cavFB->avgNum);
  cavFB->phaseVc == NULL;
  cavFB->phaseVc = (double*)malloc(sizeof(double) * cavFB->avgNum);
  
  cavFB->avgIndex = 0;
  for(i =0; i < cavFB->avgNum;i++) {
    cavFB->ampVc[i] = vAP[0];
    cavFB->phaseVc[i] = vAP[1];
  } 

  if(cavFB->PIcontrol == 0) return;
  cavFB->loopD = NULL;
  cavFB->loopD = (double*)malloc(sizeof(double) * (cavFB->loopDelay + 1) * 2);
  cavFB->loopIndex = 0;
  for(i =0; i < cavFB->loopDelay + 1 ;i++) {
    cavFB->loopD[i]               =  0;
    cavFB->loopD[i + cavFB->loopDelay + 1] = 0;
  } 

  return;
}

void put_cavVolFB_value(cavityRFVoltage_feedback_t * cavFB, double * v)
// xx = 0 -> amplitude,  1 -> phase
{
  cavFB->ampVc[cavFB->avgIndex]   = v[0];//sqrt(v[0] * v[0] + v[1] * v[1]);
  cavFB->phaseVc[cavFB->avgIndex] = v[1];//atan2(v[1],v[0]);
  cavFB->avgIndex++;
  if (cavFB->avgIndex == cavFB->avgNum) cavFB->avgIndex = 0;
  return;
}

double get_cavVolFB_value( cavityRFVoltage_feedback_t * cavFB, int xx)
// return averaged value
// xx = 0 -> ampVc,  1 -> phaseVc
{
  int i;
  double sum=0;
  if(xx==0) for(i = 0;i < cavFB->avgNum;i++) sum += cavFB->ampVc[i];
  if(xx==1) for(i = 0;i < cavFB->avgNum;i++) sum += cavFB->phaseVc[i];
  return sum/cavFB->avgNum;
}

// Vc feedback, Auto Gain Control for cavity voltage
void cavity_feedback(CAVITY_resonator_t * cav_res, double dtmp,double ptmp)
{
  double tuneRange = cav_res->cavVolFB.tuneRange;
  //double fbGain=cav_res->fbGain;
  //double fbGainP=cav_res->fbGainP;
  double fbGain=cav_res->cavVolFB.gainProp[0];
  double fbGainP=cav_res->cavVolFB.gainProp[1];
  cavityRFVoltage_feedback_t * cavFB = &(cav_res->cavVolFB);  

  switch (cav_res->mode) 
  {
    // 12 : Vg <= Vc (1), Psi <= Pc (2)
    case 12: // for conventional active cavity, geneVolt & detune angle 
    { // general conditioni, 2019/08/06 N.Yamamoto
      // note: There is an undershoot problem for cavity voltage. It makes the Satatic Robinson threshold lower. 

      if(fabs(ptmp-1)< tuneRange &&  cav_res->geneVolt > 1e-20) cav_res->geneVoltRes[0] = cav_res->geneVoltRes[0] * pow(dtmp,fbGain);
      
      //if(cav_res->detune > 0) ptmp = 1/ptmp;// detune > 0 -> wr < (n)wo, for HC
      dtmp = cav_res->detune * (pow(ptmp,fbGainP)-1);// 1/10000 gain for phase
      cav_res->detune += dtmp;
      cav_res->wr = 0.5 * cav_res->wref *(tan(cav_res->detune) / cav_res->Qload
           +sqrt( tan(cav_res->detune) / cav_res->Qload * tan(cav_res->detune) / cav_res->Qload + 4 ));
      cav_res->delOmega = cav_res->wr; //- cav_res->wref;
      cav_res->genePhase += dtmp;// genePhase should be also tuned.
      cav_res->geneVolt = fabs(cav_res->geneVoltRes[0] * cos(cav_res->detune));// geneVolt should be also changed.

      break;
    }
    // 13 : Vg <= Vc (1), Psi <= Pc (2)
    case 13: // constant detune angle, phase shifted by addition not multiplicatively
    { 
      if(fabs(ptmp)< tuneRange &&  cav_res->geneVolt > 1e-20) cav_res->geneVoltRes[0] = cav_res->geneVoltRes[0] * pow(dtmp,fbGain);
      
      cav_res->genePhase = cav_res->genePhase - ptmp*fbGainP;
      cav_res->geneVolt = fabs(cav_res->geneVoltRes[0] * cos(cav_res->detune));
      break;
    }
    case 14: // constant detune angle, phase shifted by addition not multiplicatively with integral loop
    { 
      cavFB->integralMemory[0] = cavFB->integralMemory[0]*pow(dtmp,cavFB->integralMemoryFac[0]);
      cavFB->integralMemory[1] += ptmp*cavFB->integralMemoryFac[1];
      if(fabs(ptmp)< tuneRange &&  cav_res->geneVolt > 1e-20)
	cav_res->geneVoltRes[0] = cav_res->geneVoltRes[0] * pow(dtmp,cavFB->gainProp[0])*pow(cavFB->integralMemory[0],cavFB->gainInte[0]);
      cav_res->genePhase = cav_res->genePhase - ptmp*cavFB->gainProp[1]-cavFB->integralMemory[1]*cavFB->gainInte[1];
      cav_res->geneVolt = fabs(cav_res->geneVoltRes[0] * cos(cav_res->detune));

      //cav_res->geneVolt = (cavFB->gainProp[0] * diff[0] + cavFB->gainInte[0] * cavFB->integralMemory[0]) * cav_res->IgAP[0];
      //cav_res->IgAP[1] = cav_res->IgAP[1] - (cavFB->gainProp[1] * diff[1] + cavFB->gainInte[1] * cavFB->integralMemory[1]);
      //cav_res->genePhase = cav_res->genePhase - ptmp*fbGainP;
      
      break;
    }
    case 15: // constant detune angle, phase shifted by addition not multiplicatively, +P frequency loop
    {
      double dttmp = 0;
      if(fabs(ptmp)< tuneRange &&  cav_res->geneVolt > 1e-20) {
	//dttmp = cav_res->geneVoltRes[1]-cav_res->VcAP[1]-cav_res->detuneOffset;
	dttmp = cav_res->genePhase - cav_res->detune + cav_res->detuneOffset - cav_res->VcAP[1];
	cav_res->geneVoltRes[0] = cav_res->geneVoltRes[0] * pow(dtmp,fbGain);
      }

      cav_res->detune = cav_res->detune + dttmp*fbGainP;
      cav_res->wr = 0.5 * cav_res->wref *(tan(cav_res->detune) / cav_res->Qload
           +sqrt( tan(cav_res->detune) / cav_res->Qload * tan(cav_res->detune) / cav_res->Qload + 4 ));      
      cav_res->genePhase = cav_res->genePhase - ptmp*fbGainP;
      cav_res->geneVolt = fabs(cav_res->geneVoltRes[0] * cos(cav_res->detune));
      break;
    }
    // 11 : Vg <= Vc (1), Psi <= Vc (1);
    case 11:// for active harmonic cavity?
    {
      if(cav_res->geneVolt > 1e-20) cav_res->geneVoltRes[0] = cav_res->geneVoltRes[0] * pow(dtmp,fbGain);

      cav_res->detune =  cav_res->detune * pow(dtmp,fbGainP);
      cav_res->wr = 0.5 * cav_res->wref *(tan(cav_res->detune) / cav_res->Qload
           +sqrt( tan(cav_res->detune) / cav_res->Qload * tan(cav_res->detune) / cav_res->Qload + 4 ));
      cav_res->delOmega = cav_res->wr ;//- cav_res->wref;
      cav_res->geneVolt = fabs(cav_res->geneVoltRes[0] * cos(cav_res->detune));// geneVolt should be also changed.
      break;
    }
    // 01 : Vg <= none (0), Psi <= Vc (1);
    case 1:// for passive cavity, detune angle only,
    {
      if(cav_res->detune > 0) dtmp = 1/dtmp;
      //cav_res->detune =  cav_res->detune * pow(dtmp,fbGainP); // for main cavity OK 19/06/21
      cav_res->detune =  cav_res->detune * pow(1/dtmp,fbGainP);
      cav_res->wr = 0.5 * cav_res->wref *(tan(cav_res->detune) / cav_res->Qload
           +sqrt( tan(cav_res->detune) / cav_res->Qload * tan(cav_res->detune) / cav_res->Qload + 4 ));
      //cav_res->delOmega = cav_res->wr ;//- cav_res->wref;
      cav_res->geneVolt = fabs(cav_res->geneVoltRes[0] * cos(cav_res->detune));// geneVolt should be also changed.
      break;
    }
    case 3:// for passive cavity, detune angle only, non-exponential
    {
      //if(cav_res->detune > 0) dtmp = 1/dtmp;
      //cav_res->detune =  cav_res->detune * pow(dtmp,fbGainP); // for main cavity OK 19/06/21
      cav_res->detune =  cav_res->detune - (dtmp - 1) * fbGainP / tan(cav_res->detune);
      cav_res->wr = 0.5 * cav_res->wref *(tan(cav_res->detune) / cav_res->Qload
           +sqrt( tan(cav_res->detune) / cav_res->Qload * tan(cav_res->detune) / cav_res->Qload + 4 ));
      //cav_res->delOmega = cav_res->wr ;//- cav_res->wref;
      cav_res->geneVolt = fabs(cav_res->geneVoltRes[0] * cos(cav_res->detune));// geneVolt should be also changed.
      break;
    }    
    // 10 : Vg <= Vc (1), Psi <= none (0);
    case 10: // This method works when the cavity is tuned (cavity phase is alomst eqaul to the target value). 
    { 
      if(fabs(ptmp-1)< tuneRange &&  cav_res->geneVolt > 1e-20) cav_res->geneVoltRes[0] = cav_res->geneVoltRes[0] * pow(dtmp,fbGain);
      cav_res->geneVolt = fabs(cav_res->geneVoltRes[0] * cos(cav_res->detune));// geneVolt should be also changed.
      break;
    }
    // 20 : Vg <= PhaseVc(2), Psi <= none (0);
    case 20: // This method works when the cavity is tuned (cavity phase is alomst eqaul to the target value). 
    { 
      if(fabs(ptmp-1)< 0.001 &&  cav_res->geneVolt > 1e-20) cav_res->geneVoltRes[0] = cav_res->geneVoltRes[0] * pow(ptmp,fbGain);
      cav_res->geneVolt = fabs(cav_res->geneVoltRes[0] * cos(cav_res->detune));// geneVolt should be also changed.
      break;
    }
    // 30 : Vg,PhaseVg <= Vc (1),PhaseVc, Psi <= none (0);
    case 30: //  
    { 
      if(cav_res->geneVolt > 1e-20) cav_res->geneVoltRes[0] = cav_res->geneVoltRes[0] * pow(dtmp,fbGain);
      cav_res->geneVolt = fabs(cav_res->geneVoltRes[0] * cos(cav_res->detune));// geneVolt should be also changed.
      // not work, genePhase should be also feedbacked ?
      break;
    }
    // 02 : Vg <= none (0), Psi <= PhaseVc (2);
    case 2:
    { // 2019/06/20 MC test -> OK

      //if(cav_res->detune > 0) ptmp = 1/ptmp;// detune > 0 -> wr < (n)wo, for HC

      dtmp = cav_res->detune * (pow(ptmp,fbGainP)-1);// 1/10000 gain for phase
      cav_res->detune += dtmp;
      cav_res->wr = 0.5 * cav_res->wref *(tan(cav_res->detune) / cav_res->Qload
       +sqrt( tan(cav_res->detune) / cav_res->Qload * tan(cav_res->detune) / cav_res->Qload + 4 ));
      cav_res->delOmega = cav_res->wr; //- cav_res->wref;
      cav_res->genePhase += dtmp;// genePhase should be also tuned.
      cav_res->geneVolt = fabs(cav_res->geneVoltRes[0] * cos(cav_res->detune));// geneVolt should be also changed.
      break;
    }
    // 22 : Vg <= Pc (0), Psi <= PhaseVc (2);
    case 22:
    { // NOT TEST
      if(cav_res->geneVolt > 1e-20) cav_res->geneVoltRes[0] = cav_res->geneVoltRes[0] * pow(ptmp,fbGain);

      dtmp = cav_res->detune * (pow(ptmp,fbGainP)-1);// 1/10000 gain for phase
      cav_res->detune += dtmp;
      cav_res->wr = 0.5 * cav_res->wref *(tan(cav_res->detune) / cav_res->Qload
         +sqrt( tan(cav_res->detune) / cav_res->Qload * tan(cav_res->detune) / cav_res->Qload + 4 ));
      cav_res->delOmega = cav_res->wr; //- cav_res->wref;
      cav_res->genePhase += dtmp;// genePhase should be also tuned.
      cav_res->geneVolt = fabs(cav_res->geneVoltRes[0] * cos(cav_res->detune));// geneVolt should be also changed.
      break;
    }
  }
        
  cav_res->fillrate = cav_res->wr * 0.5 / cav_res->Qload; 
  cav_res->lossfactor = cav_res->wr * 0.5 * cav_res->Rs / cav_res->Qzero;
  return;
}

// Vc feedback, Auto Gain Control for cavity voltage
// dtmp,ptmp : difference
//
//      ref(V,P)
//      |
// Vc   V
//  -->(-)-----------> (* Pgain)------------------------->(+) ---->//--->
//       diff        |                                     |          Vg or Ig
//                   ---->(* Igain)--->(sum)--->(/2^X)-----> 
//
//             offset
//               |
//               V
//  Pgen--(-)-->(-)-->(*Pgain)--->(+ detune)
//         |    
//  Pcav----     
//             
void cavity_PIfeedback(CAVITY_resonator_t * cav_res, const ring_t * ring, double dtmp,double ptmp,double ttmp)
  //dtmp = 1- vcmon[0]/cav_res->Vc[0];
  //ptmp = 1- vcmon[1]/cav_res->Vc[1];
{
  int i;
  double diff[3];
  double igIQ[2];
  double fac=1;//e-3;//1e-7;
  double tunerFac=1e-3;

  cavityRFVoltage_feedback_t * cavFB = &(cav_res->cavVolFB);

  diff[0] = dtmp * fac;
  diff[1] = ptmp * fac;
  diff[2] = ttmp; // tan(get_cavVolFB_value(&cav_res->cavVolFB,1) - cav_res->geneVoltRes[1] - cav_res->tunerOffset)
  igIQ[0]= cav_res->IgAP[0] * cos(cav_res->IgAP[1]);
  igIQ[1]= cav_res->IgAP[0] * sin(cav_res->IgAP[1]);

  if(fabs(ttmp) > cavFB->tuneRange) ttmp = 0.0;
  if(fabs(ttmp) < cavFB->tuneDeadRange) ttmp = 0.0;

  switch (cav_res->mode) 
  { 
    // 12 : Vg <= Vc (1), Psi <= Pc (2)
    case 12: // for conventional active cavity, geneVolt & detune angle 
    { // IQ based feedback general condition, 2022/05/16 N.Yamamoto

      // Genrerator PID FB  
      for (i=0;i<2;i++)
      {
        cavFB->integralMemory[i] += diff[i] * cavFB->integralMemoryFac[0];
        igIQ[i] += (cavFB->gainProp[0]  *  diff[i] + cavFB->gainInte[0] * cavFB->integralMemory[i]) * igIQ[i]; //cav_res->Vc[i];
      }
      cav_res->IgAP[0] = sqrt(igIQ[0]*igIQ[0] + igIQ[1]*igIQ[1]);
      cav_res->IgAP[1] = atan2(igIQ[1],igIQ[0]);
      Ig2VgIQ(cav_res,ring,igIQ);
      //Tuner 
      cav_res->detune -= cavFB->gainProp[2] * diff[2] * tunerFac;
      break;
    }
    case 112: // for conventional active cavity, geneVolt & detune angle 
    { // AP based feedback general condition, 2022/05/16 N.Yamamoto

      // Genrerator PID FB  
      for (i=0;i<2;i++)
      {
        cavFB->integralMemory[i] += diff[i] * cavFB->integralMemoryFac[i];
        cav_res->IgAP[i] += (cavFB->gainProp[i] * diff[i] + cavFB->gainInte[i] * cavFB->integralMemory[i]) * cav_res->IgAP[i];
      }

      Ig2Vg(cav_res,ring);
      //Tuner 
      cav_res->detune -= cavFB->gainProp[2] * diff[2] * tunerFac;
      break;
    }
    case 113: // for conventional active cavity, geneVolt & detune angle 
    { // AP based feedback general condition, 2022/05/16 N.Yamamoto

      // Genrerator PID FB  
      for (i=0;i<2;i++)
      {
        cavFB->integralMemory[i] += diff[i] * cavFB->integralMemoryFac[i];
      }
      cav_res->IgAP[0] = (cavFB->gainProp[0] * diff[0] + cavFB->gainInte[0] * cavFB->integralMemory[0]) * cav_res->IgAP[0];
      cav_res->IgAP[1] = cav_res->IgAP[1] - (cavFB->gainProp[1] * diff[1] + cavFB->gainInte[1] * cavFB->integralMemory[1]);

      Ig2Vg(cav_res,ring);
      //Tuner 
      //cav_res->detune -= cavFB->gainProp[2] * diff[2] * tunerFac;
      break;
    }
 
  }
  cav_res->wr = 0.5 * cav_res->wref *(tan(cav_res->detune) / cav_res->Qload
           +sqrt( tan(cav_res->detune) / cav_res->Qload * tan(cav_res->detune) / cav_res->Qload + 4 ));
  cav_res->delOmega = cav_res->wr; //- cav_res->wref;
  cav_res->fillrate = cav_res->wr * 0.5 / cav_res->Qload; 
  cav_res->lossfactor = cav_res->wr * 0.5 * cav_res->Rs / cav_res->Qzero;
  return;
}



void Vg2Ig(CAVITY_resonator_t * cav_res)
// kn * ig~ = Fillrate (1-itanPsi) * Vg~
{
  double igIQ[2],VgIQ[2];
  double fac = cav_res->fillrate/cav_res->lossfactor;// [1/s] / [Ohm]
  double tanPsi = tan(cav_res->detune);

  //ig = cav_res->fillrate/cav_res->lossfactor*(1-i*tan(cav_res->detune))*cav_res->geneVolt;
  VgIQ[0] = cav_res->geneVolt * cos(cav_res->genePhase);// genevolt: V / eV
  VgIQ[1] = cav_res->geneVolt * sin(cav_res->genePhase);

  igIQ[0] =  fac * (VgIQ[0] + VgIQ[1] * tanPsi);
  igIQ[1] =  fac * (-VgIQ[0] * tanPsi + VgIQ[1]);

  cav_res->IgAP[0] = sqrt(igIQ[0]*igIQ[0] + igIQ[1]*igIQ[1]);
  cav_res->IgAP[1] = atan2(igIQ[1],igIQ[0]);
  return;
}

void Ig2Vg(CAVITY_resonator_t * cav_res, const ring_t * ring)
{
  double v[]={0,0};
  Ig2VgIQ(cav_res,ring,v);
}

void Ig2VgIQ(CAVITY_resonator_t * cav_res, const ring_t * ring, double igIQ[2] )
{
  //double igIQ[2];
  double oldVgIQ[2],VgIQ[2];
  double Trf = ring->T0/ring->Nharm;
  double fac=cav_res->fillrate * Trf;
  double lossT = cav_res->lossfactor * Trf;
  double Vgfac = exp(-fac);
  double tanPsi = tan(cav_res->detune)*cav_res->fillrate* Trf;

  if(igIQ[0]+igIQ[1] == 0)
  {
    igIQ[0]= cav_res->IgAP[0] * cos(cav_res->IgAP[1]);
    igIQ[1]= cav_res->IgAP[0] * sin(cav_res->IgAP[1]);
  }
 
  oldVgIQ[0] = cav_res->geneVolt * cos(cav_res->genePhase);
  oldVgIQ[1] = cav_res->geneVolt * sin(cav_res->genePhase);

  VgIQ[0] = Vgfac * ( oldVgIQ[0] * cos(tanPsi) - oldVgIQ[1] * sin(tanPsi)) + lossT * igIQ[0];
  VgIQ[1] = Vgfac * ( oldVgIQ[0] * sin(tanPsi) + oldVgIQ[1] * cos(tanPsi)) + lossT * igIQ[1];

  //printf("\n %e, %e, %e\n",Vgfac, lossT, Trf);

  cav_res->geneVolt = sqrt(VgIQ[0]*VgIQ[0] + VgIQ[1]*VgIQ[1]);
  cav_res->genePhase = atan2( VgIQ[1],VgIQ[0]);// + cav_res->detune;
  cav_res->geneVoltRes[0]= cav_res->geneVolt / cos(cav_res->detune);
  cav_res->geneVoltRes[1]= cav_res->genePhase - cav_res->detune;//cav_res->IgAP[1];

  return;
}
 
int fprint_cavresonator(FILE *fp,const ring_t ring)
{
  int k;
  fprintf(fp, " ===================================================================\n");
  fprintf(fp, "   Cavity-type Resonators with Transient Beam Loading Effect Input: \n");
  fprintf(fp, " ===================================================================\n");
  fprintf(fp, "  Number of cavity-type resonators: %u\n", ring.cavity_resonators_size);
  if (ring.cavity_resonator_main > 0)
    fprintf(fp, "  Index of detected main cavity: %u\n", ring.cavity_resonator_main);
  if (ring.cavity_resonators_size > 0)
    fprintf(fp, "  Nbu max    = %d = %d turns,  order:  %u\n", ring.Nbumax, (int)ring.Nbumax/ring.Nharm, ring.lr_order);
  int scan_num;
  for(k = 0; k < ring.cavity_resonators_size; k++) 
  {
    const CAVITY_resonator_t resonator = ring.cavity_resonators[k];
    fprintf(fp, "  harm. [%d]     = %8.8lf,         R [%d]      = %8.2lf [MOhm], \n",
            k, resonator.m, k, resonator.Rs/FMEGA);
    fprintf(fp, "  Qzero [%d]     = %8.2lf ,        Qload [%d]  = %8.2lf, \n",
            k, resonator.Qzero, k, resonator.Qload);
    fprintf(fp, "  Cavity Voltage [%d]  = %8.8lf [MV], Phase [%d]     = %8.5lf [rad]\n",
            k, resonator.VcAP[0]/FMEGA*(ring.E0 * FGIGA),	k, resonator.VcAP[1]);
    //fprintf(fp, "  geneIndVol (Vg) [%d]     = %8.6lf [MW],    genePhase [%d]      = %8.5lf [rad], \n",
    //        k, resonator.geneVolt/FMEGA*(ring.E0 * FGIGA), k, resonator.genePhase);
    fprintf(fp, "  wr [%d]    = 2PI x %8.8lf [GHz rad], delta Omega = 2PI x %.6e [Hz rad]\n",
            k, (resonator.wr+resonator.wrOffset)/FGIGA*0.5/M_PI, (resonator.wr-resonator.wref+resonator.wrOffset)*0.5/M_PI);
    fprintf(fp, "  wref [%d]  = 2PI x %8.8lf [GHz rad]\n",
            k, resonator.wref*0.5/M_PI/FGIGA);  
    fprintf(fp, "  detune [%d] = %8.5lf [rad], offset = %8.5lf [rad] (%.2lf [deg])\n",
            k, resonator.detune+resonator.detuneOffset, resonator.detuneOffset,resonator.detuneOffset/M_PI*180);  
    fprintf(fp,"  Generator voltage at resonance [%d] = %e [V]\n",k,resonator.geneVoltRes[0]*(ring.E0 * FGIGA));
    fprintf(fp,"  geneVolt [%d] = %e [V], geneModVolt_1 = %e [V], _2 = %e [V], _3 = %e [V]\n",
            k, resonator.geneVolt*(ring.E0 * FGIGA), resonator.geneModVolt[0]*(ring.E0 * FGIGA),resonator.geneModVolt[1]*(ring.E0 * FGIGA),resonator.geneModVolt[2]*(ring.E0 * FGIGA));  
    fprintf(fp, "  genePhase [%d] = %8.5lf [rad], geneModPhaseOffsetM_1 = %d, _2 = %d, _3 = %d \n",
            k, resonator.genePhase, resonator.geneModPhaseOffsetM[0],resonator.geneModPhaseOffsetM[1],resonator.geneModPhaseOffsetM[2]);  
    fprintf(fp, "  geneModfreqMrev [%d] = %d, wMod_1 = 2PI x %9.8lf [GHz rad]\n",
            k, resonator.geneModfreqMrev, resonator.wMod[0]*0.5/M_PI/FGIGA);  
    fprintf(fp, "  Nbu [%d]    = %d RF buckets = %d turns\n", k, resonator.Nbu, resonator.Nturn);
    fprintf(fp, "  Cavity Switch On turn [%d]= %d,  Switch Off turn [%d] = %d\n",
            k, resonator.switchOnTurn, k, resonator.switchOffTurn);
    fprintf(fp, "  FBmode [%d] = %d,  FBPIcontrol [%d] = %d,  FBlimit [%d] = %.8lf,  FBaverage [%d] = %d turn\n",
            k, resonator.mode, k,resonator.cavVolFB.PIcontrol, k, resonator.cavVolFB.limitRange,k,resonator.cavVolFB.avgNum);
    if(resonator.cavVolFB.PIcontrol==0) 
      fprintf(fp, "  FBgain [%d] = %.4e,  FBgainP [%d] = %.4e\n",
            k, resonator.cavVolFB.gainProp[0], k, resonator.cavVolFB.gainProp[1]);
    else {
      fprintf(fp, "  FBLoopDelay [%d] = %d, FBSample [%d] = %d (Unit in bucket number)\n",
		      k, resonator.cavVolFB.loopDelay,k, resonator.cavVolFB.sample);
      if(resonator.mode > 100)
      {
        fprintf(fp, "  FBPgainAmp [%d] = %.5lf,  FBPgainPhase [%d] = %.5lf,  FBPgainTuner [%d] = %.5lf\n",
		      k, resonator.cavVolFB.gainProp[0], k, resonator.cavVolFB.gainProp[1], k, resonator.cavVolFB.gainProp[2]);
        fprintf(fp, "  FBIgainAmp [%d] = %.5lf,  FBIgainPhase [%d] = %.5lf,  FBIntegralFacorAmp [%d] = %.3e,  FBIntegralFactorPhase [%d] = %.3e\n",
		      k, resonator.cavVolFB.gainInte[0], k, resonator.cavVolFB.gainInte[1], k, resonator.cavVolFB.integralMemoryFac[0],k, resonator.cavVolFB.integralMemoryFac[1]);
      } else {
        fprintf(fp, "  FBPgain [%d] = %.5lf,  FBIgain [%d] = %.5lf,  FBPgainTuner [%d] = %.5lf\n",
		      k, resonator.cavVolFB.gainProp[0], k, resonator.cavVolFB.gainInte[0], k, resonator.cavVolFB.gainProp[2]);
      }
    }
    fprintf(fp, "  TuneRange [%d] = %8.5f, TuneDeadRange [%d] = %8.5f,  FourierAmp [%d]= %8.5f,  FourierPhase [%d] = %8.5f\n",
            k, resonator.cavVolFB.tuneRange, k, resonator.cavVolFB.tuneDeadRange, k, resonator.FourierAmp, k, resonator.FourierPhase);
    fprintf(fp, "  Cavity Voltage RF feedback : switchOnTurn [%d] = %d [Turn], switchOffTurn [%d] = %d [Turn]\n",
            k, resonator.cavVolFB.switchOnTurn, k, resonator.cavVolFB.switchOffTurn);
    fprintf(fp, "  Direct RF feedback : Gain [%d] = %8.5f, phaseShift [%d]= %8.5f [rad],  loop delay [%d] = %d [turn] \n",
            k, resonator.directRF_feedback.gain, k, resonator.directRF_feedback.phaseShift, k, resonator.directRF_feedback.loopDelay);
    fprintf(fp, "  Direct RF feedback : switchOnTurn [%d] = %d [Turn], switchOffTurn [%d] = %d [Turn]\n",
            k, resonator.directRF_feedback.switchOnTurn, k, resonator.directRF_feedback.switchOffTurn);
    //fprintf(fp, "  CavIndInst [%d] = %d (0:OFF 1:ON)\n\n", k, resonator.EnableCavInst);
    for(scan_num = 0; scan_num < resonator.scan_size; scan_num++) 
    {
      const cavScan_t pscan = resonator.scan[scan_num];
      fprintf(fp, "  CAVscan_%d_type [%d] = %s, CAVscan_%d_start [%d] = %8.5f, CAVscan_%d_step [%d] = %8.5f\n", 
		      scan_num+1, k, cavScan_str[pscan.type],scan_num+1, k, pscan.start,scan_num+1, k, pscan.step);
      fprintf(fp, "  CAVscan_%d_nrev [%d] = %d\n",scan_num+1, k, pscan.nrev);

    }
    fprintf(fp, "\n");
  } 

  return 0; 
}

int fprint_cavresonator_last(long int rev)
{
  int l;
  char filename_hist_rf[FILENAME_MAX] = "";
  snprintf(filename_hist_rf, FILENAME_MAX, "%s/input_0.log", track.work_path);
  FILE * fp = fopen(filename_hist_rf,"a");
  if(fp == NULL) ERROR("fopen_input_end", return -1);
  fprintf(fp, "\n  ===================================================================");
  fprintf(fp, "\n     Last parameters for Cavity-type Resonator at %ld turns",rev);
  fprintf(fp, "\n  ===================================================================");
  for(l = 0; l < ring.cavity_resonators_size; l++)
  {
    CAVITY_resonator_t * cav_res = &(ring.cavity_resonators[l]);
    
    //fprintf(fp,"\n  FBmode [%d]    = %d, FBgain [%d]    = %8.5lf",l,cav_res->mode,l,cav_res->cavVolFB.gainProp[0]);
    fprintf(fp,"\n  Target Cavity voltage [%d]  = %e [V], Phase [%d] = %8.7lf [rad]",
	l,cav_res->VcAP[0]*(ring.E0 * FGIGA),l,cav_res->VcAP[1]);
    fprintf(fp,"\n  Generator voltage at resonance [%d] = %e [V]",l,cav_res->geneVoltRes[0]*(ring.E0 * FGIGA));
    if(cav_res->directRF_feedback.gain > 0) fprintf(fp,"\n  Generator voltage [%d] = %e [V], phase [%d] = %8.7lf [rad]",l,
		    sqrt((cav_res->geneVolt*cos(cav_res->genePhase)+cav_res->directRF_feedback.Vg[0])*(cav_res->geneVolt*cos(cav_res->genePhase)+cav_res->directRF_feedback.Vg[0])+ 
			 (cav_res->geneVolt*sin(cav_res->genePhase)+cav_res->directRF_feedback.Vg[1])*(cav_res->geneVolt*sin(cav_res->genePhase)+cav_res->directRF_feedback.Vg[1]))*(ring.E0 * FGIGA),
		    l,atan2(cav_res->geneVolt*sin(cav_res->genePhase)+cav_res->directRF_feedback.Vg[1],cav_res->geneVolt*cos(cav_res->genePhase)+cav_res->directRF_feedback.Vg[0]));
    else fprintf(fp,"\n  Generator voltage [%d] = %e [V], phase [%d] = %8.7lf [rad]",l,cav_res->geneVolt*(ring.E0 * FGIGA),l,cav_res->genePhase);
    fprintf(fp,"\n  wr     [%d] = 2PI x %.8lf [GHz rad]",l,cav_res->wr*0.5/M_PI/FGIGA);
    fprintf(fp,"\n  detune [%d] = %8.7lf [rad]",l,cav_res->detune);
    fprintf(fp,"\n  Vc at last turn [%d]   = %e [V], phase [%d] = %8.7lf [rad]\n",
	l,cav_res->VcTrackAP[0]*(ring.E0 * FGIGA),l,cav_res->VcTrackAP[1]);
   
  }
  fclose(fp);
  return 0;
}

void  cavres_destory(ring_t * ring, const selffield_model_t * SelfFieldModel)
{
  int l;
  //int Nbin = SelfFieldModel->Ncell;
  if (ring->cavity_resonators_cavFB_fp!=NULL) fclose(ring->cavity_resonators_cavFB_fp);
  if (ring->cavity_resonators_drfFB_fp!=NULL) fclose(ring->cavity_resonators_drfFB_fp);
  for(l = 0; l < ring->cavity_resonators_size; l++)
  {
    CAVITY_resonator_t * cav_res = &(ring->cavity_resonators[l]);
    //if(cav_res->mode >0 && cav_res->cavVolFB.avgNum > 0){
      free(cav_res->cavVolFB.ampVc);
      free(cav_res->cavVolFB.phaseVc);
      cav_res->cavVolFB.ampVc = NULL;
      cav_res->cavVolFB.phaseVc = NULL;
    //}
  }
  return;
}

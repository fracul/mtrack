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

#ifdef MBTRACK_CAVR_H


bool read_conf_cav_res(FILE * fp, ring_t * ring, tracking_t * track);
bool setup_cavity_parameters(ring_t * ring);
int fprint_cavresonator(FILE *fp,const ring_t ring);
int fprint_cavresonator_last(long int rev);
void cavres_destory(ring_t * ring,const selffield_model_t * SelfFieldModel);

/**
 * Initializes the phasor scheme and computes initial values
 * @param phasor_end result variable
 */
void
wake_phasor_initTR(ring_t * ring, double * fnp_ring, const selffield_model_t * SelfFieldModel,
	       	weak_bunch_t * bunch, int kb, double * phasor_end, e_beam_t * ebeam);

/**
 * Computes sum of phasors due to actual particle distribution
 * @param lr_wake result variable
 */
void
construct_wake_phasorTR(double * lr_wake, double * phasor_end, 
		      const selffield_model_t * SelfFieldModel,
		      const ring_t * ring, int kb, const double * fnp, e_beam_t * ebeam, weak_bunch_t * bunch);

void
wake_phasor_Recalc(ring_t * ring,CAVITY_resonator_t * cav_res, int l,double * fnp_ring, const selffield_model_t * SelfFieldModel, 
		  weak_bunch_t * bunch, int kb, double * phasor_end, e_beam_t * ebeam);
//wake_phasor_Recalc(const ring_t * ring,CAVITY_resonator_t * cav_res, int l,const double * fnp_ring, const selffield_model_t * SelfFieldModel, 
//		  weak_bunch_t * bunch, int kb, double * phasor_end, e_beam_t * ebeam);

void calc_resonatorVb(int frag, ring_t * ring, double * fnp_ring, const selffield_model_t * SelfFieldModel, 
		e_beam_t * ebeam, double bunchNp,CAVITY_resonator_t * cav_res, double detune, double *vbout);
void calc_resonatorVg(CAVITY_resonator_t * cav_res, double *calcVb);
void set_resonatorDetune(CAVITY_resonator_t * cav_res, double detune);
void set_resonatorDetuneOffset(CAVITY_resonator_t * cav_res);

void Vg2Ig(CAVITY_resonator_t * cav_res);
void Ig2Vg(CAVITY_resonator_t * cav_res, const ring_t * ring);
void Ig2VgIQ(CAVITY_resonator_t * cav_res, const ring_t * ring, double igIQ[2] );

//// RF feedback
// cavityVoltage_feedback
void cavity_feedback(CAVITY_resonator_t * cav_res, double dtmp,double ptmp);
void cavity_PIfeedback(CAVITY_resonator_t * cav_res,const ring_t * ring, double dtmp,double ptmp,double ttmp);
void highSpeedFB(CAVITY_resonator_t * cav_res,const ring_t * ring, int l,long int turn, int bucket,int kb, double * vbmon);
void init_cavVol_feedback(cavityRFVoltage_feedback_t * cavFB, double * vAP);
void put_cavVolFB_value( cavityRFVoltage_feedback_t * cavFB, double * v);
double get_cavVolFB_value( cavityRFVoltage_feedback_t * cavFB, int xx);

//directRF_feedback
void init_directRF_feedback(CAVITY_resonator_t * cav_res,double *v, FILE *fp, int frag, int cav);
void put_directRF_feedback(CAVITY_resonator_t * cav_res, FILE *fp,int frag, int cav, long int turn);

#endif

#ifndef MBTRACK_CAVR_H
#define MBTRACK_CAVR_H


typedef enum { // check cavityresonator.c if the definition is changed.
  detuneRAD = 0,
  fr = 1
} cavScanType_t;


typedef struct cavScan
{ // 2020/04/14 N.Yamamoto
  //int type;
  	// 1 -> detuneRad
	// candidate: Vc, syncPhase, wr for HOM
  double start;/** Start and step value of the parameter */
  double step;/** Start and step value of the parameter */
  int nrev; /** Number of turns per scan step */
  //char scanType[80];
  cavScanType_t type;
}
cavScan_t;


typedef struct directRF_feedback
{ // 2020/04/14 N.Yamamoto
  double gain; /** Gain of the feedback, |Vg_fb| = gain * |Vg| */
  //double gain4mainVg; /** Gain of the main Vg feedback. This is used for tuning the AmpRatio */
  double phaseShift; /** phase shift of the feedback, ang(Vg_fb) = ang(moniterd Vc) + phaseShift */
  //int avgNum;
  int loopDelay;  /** loopDelay of the feedback, unit: Turn */
  int loopIndex;  /** Index of loop memory */
  int switchOnTurn; /** If turn > switchOnTurn, the feedback ON */
  int switchOffTurn; /** If turn > switchOffTurn, the feedback OFF */
  double *loopV;
  double fbAmp;
  double fbPhase;
  double Vg[2]; /** real[0] or imaginary[1] part of drfFB voltage */
  double AmpRatio;
  //double initialSetupAccuracy;
}
directRF_feedback_t;

typedef struct cavityRFVoltage_feedback
{ // 2022/05/16 N.Yamamoto
  int PIcontrol;
  int sample; // sampling interval indicated by bucket number 
  double gainProp[3];// 0: Amp, 1:phase & 2:tuner gain
  double gainInte[3];// 0: Amp, 1:phase & 2:tuner gain
  double integralMemory[2];// Integral values; 0: Amp, 1:phase
  double integralMemoryFac[2];// Integral factor; ex 1/2**10
  double limitRange; /**< Limit for RF feedback, if diff < limitRange, FB is canceled */
  double tuneRange; /**< geneVolt feedback works when cavity is tuned (cavity phase closes to the target value)*/
  double tuneDeadRange; /**< Tuner feedback does not work in lower value*/
  //double gain; /** Gain of the feedback, |Vg_fb| = gain * |Vg| */
  //double phaseShift; /** phase shift of the feedback, ang(Vg_fb) = ang(moniterd Vc) + phaseShift */
  int avgNum;  /** avgerage number of the feedback, unit: Turn */
  int avgIndex;  /** Index of loop memory */
  int loopIndex;  /** Index of loop memory */
  int loopDelay;  /** loopDelay of the feedback, unit: Turn */
  int switchOnTurn; /** If turn > switchOnTurn, the feedback ON */
  int switchOffTurn; /** If turn > switchOffTurn, the feedback OFF */
  double *ampVc;
  double *phaseVc;
  double *loopD;
}
cavityRFVoltage_feedback_t;


/**
 * \Cavity-type Resonator, can have SelfField effect and/or long-range Wake Field effect
 * \ 09/07/2018, created by Naoto Yamamoto
 */
typedef
struct CAVITY_resonator
{
  int mode; /**< Operation mode,  0:active, 1:passive (geneVolt=0) */
  double m; /**< higher harmonic: m*wrf = wr, don't has to be an integer */
  double Rs; /**< Impedance [ohm] in LON plane */
  double Vc[2]; /**< Target Cavity voltage [V], 0 : real & 1: imaginary */
  double VcAP[2]; /**< Target Cavity voltage [V], 0 : amplitude & 1: phase */
  double VcTrackAP[2]; /**< Tracked Cavity phase at last turn */
  double Vbr; /**< Beam loading at resonance */ 
  double delOmega; /**< cavity detuning [Hz rad]: wr = m*wrf + delOmega , for phasor calculation*/
  double detune; /**< cavity detuning [rad]*/
  double detuneOffset; /**< cavity detuning Offset [rad] */
  double Qload; /**< Loaded Quality factor = (Qzero)/(1+beta)*/
  double Qzero; /**< unlodade Quality factor */
  double genePhase;  /**< generator induced RF phase with respect to beam axis (see Wilson's phasor diagram) [rad] , for phasor calculation*/
  double geneVolt; /**< generator induced voltage [V] , for phasor calculation*/
  double geneVoltRes[2]; /**< generator induced voltage at resonance [V] , for phasor calculation; 0:Amp, 1:phase */
  double IgAP[2];/**< generator current; 0:Amp, 1:phase */
  double tunerOffset; /**< Tuner offset */
  double wr; /**< Angular resonance frequency [Hz rad] = 2*M_PI*freq */
  double wref; /**< Angular frequency of reference frame[Hz rad], nearest frequency of multiple revolution */
  double wrOffset; /**< Offset of Angular resonance frequency [Hz rad] */
  double fillrate;  /**< inverse of cavity filling time [1/s] = 0.5*wr/Qload) , for phasor calculation*/
  double lossfactor; /**< energy loss by a charge passing through a cavity = 0.5*wr*Rs/Qzero , for phasor calculation*/
  double fbGain; /**< Gain for RF feedback, factor of change angle to diffrence */
  double FourierAmp; /**< Fourier component for beam loading, to modify Vbr calculation  */
  double FourierPhase;
  unsigned int fbMon;
  unsigned int Nbu; /**< History of Nbu bunches is taken into account */
  unsigned int Nturn; /**< History of Nturn turns is taken into account */

  int EnableCavInst; /**< Cavity induced Instability is supressed. (Calculate in deltaOmega = wr,) 0:OFF, 1:ON */
  double fbGainP; /**< Gain for detune angle feedback,0*/

  //for compeansation 190809 N.Yamamoto
  int geneModfreqMrev;/**< modulated generator parameters for kicker cavity*/  
  double geneModVolt[5];/**< modulated generator parameters for kicker cavity*/
  double wMod[5];/**< modulated generator parameters for kicker cavity*/
  int geneModPhaseOffsetM[5];/**< modulated generator parameters for kicker cavity*/
  
  int FBparaIniSet; /**< For cavity feedback with initial detune offset */
  int switchOnTurn;
  int switchOffTurn;

  cavityRFVoltage_feedback_t cavVolFB;
  directRF_feedback_t directRF_feedback;

  unsigned scan_size;
  cavScan_t * scan;
}
CAVITY_resonator_t;

#endif /* MBTRACK_CAVR_H */

/**
 * @file
 * Accelerator ring and tracking parameters structures
 */

#ifndef MBTRACK_TYPES_H
#define MBTRACK_TYPES_H

#include <stdio.h>
#include "def.h"
#include "cyclic_array.h"

/**
 * \brief Resonator, can have SelfField effect and/or long-range Wake Field effect
 */
typedef struct resonator
{
  double Rs; /**< Impedance [ohm/m] in LON plane, or [ohm] */
  double wr; /**< Angular resonance frequency [rad/s] = 2*M_PI*freq */
  double Qfactor; /**< Quality factor */
  plane_t plane;  /**< Plane: LON, VER, HOR */
}
resonator_t;

typedef
struct LR_resonator
{
  double Rs; /**< Impedance [ohm] in LON plane */
  double wr; /**< Angular resonance frequency [rad/s] = 2*M_PI*freq */
  double Qfactor; /**< Quality factor */
  int m; /**< higher harmonic: m*wrf = wr */
  double detune; /**< HC detuning: wr = m*wrf + detune*2pi */
  unsigned int Nbu; /**< History of Nbu bunches is taken into account */
  unsigned int Nturn; /**< History of Nturn turns is taken into account */
  double * facc; /**< Cosine coefficient of wake field and it's derivatives */
  double * facs; /**< Sine coefficient of wake field and it's derivatives */
}
LR_resonator_t;

typedef struct active_HC
{
  double Vpeak; /** Peak voltage */ 
  double nHC; /** higher harmonic of rf cavity, don't has to be an integer */ 
  double phi_aHC; /** phase with respect to rf cavity of first bunch*/
}
active_HC_t;

/**
 * \brief Ring and machine parameters
 */
typedef struct ring
{
  unsigned int Nharm; /**< Harmonic number */
  double h; /**< Harmonic number, as non-truncated double */
  double E0; /**< Nominal energy [GeV] */
  double Gamma, Gamma2, Gamma3; /**< Relativistic gamma and its powers */
  double Iring; /**< Total ring current [mA] */
  double * Ibunch; /**< Bunch current [mA], WARNING: assuming equal current in bunches */
  double emittanceH; /**< Hor. emittance [nm] */
  double couplbeta; /**< emittance ration (coupling) */
  double Gzix, Gziz; /**< Normalised chromaticity */
  double beffL[2], beffV[2], beffH[2]; /**< Long. effetive beampipe radius for 1rd and 2th power , vert. and hor. effective beampipe radius for 3rd and 4th power */
  double fYokoyaV, fYokoyaH;
  double ApertureH, ApertureV;  /**< Half aperture in [mm] */
  double Vrf0; /**< Total Effective RF Voltage [MV]*/
  double frf; /**< RF Frequency [MHz] */
  double Je; /**< Long. damping partition, WARNING: Je is not take from input but calculated if U0 and taue are given */
  double Jx, Jz; /**< Transverse damping partition*/
  double rho0; /**< Radius of curv (bends) [m] used to calculate U0, not used if U0 is given  */
  double Lc; /**< Ring circumference [m] */
  double ac; /**< Momentum compaction */
  double R; /**< Average ring radius = Lc / 2pi [m] */
  double T0; /**< Revolution time [s] */
  double w0; /**< Angular rev. frequency */
  double U0; /**< Energy Loss per Turn [keV], calculated from E0 and rho0 if not given */
  double Urad;   /**< Relative loss per turn = U0/E0 */
  double Vrfp, q, sn0;
  double wrf; /**< Angular rf-frequency */
  double phai0; /**< Synchrounus phase */
  double De; /**< energy dependent loss */
  double aexpe; /**< Energy damping increment = 1/tauE and its powers */
  double fso, wso, wso2; /**< (Angular) synchrotron frequency and its powers */
  double Fq, epsMax;
  double fc1, fc12; /**< [s] */
  double taue; /**< [s] energy damping time, calculated from Je if not given */
  double taux, tauz; /** [s] transverse damping times */
  double tauC;
  int m_aHC; /** Harmonic of ACTIVE HC */
  double phi_n;
  double HC_k;
  
  double QH0, QV0, Qso; /**< Tunes, H: hor., V: vert., s: sync. */
  double wgziH, wgziV; /**< [Hz] */
  double beta1[3], alpha1[3], gamma1[3];
  double dispH1,   disppH1;
  
  double rhorw; /**< Wall resistivity */
   
  resonator_t * resonators; /**< List of resonators used for selffields */
  unsigned resonators_size; /**< Number of shortrange resonators for selffield */
  
  LR_resonator_t * longrange_resonators;  /**< Longrange resonator, acts over several turns */
  unsigned longrange_resonators_size; /** number of lr resonators */
  unsigned Nbumax; /** longest Nbu of all LR resonators -> length of cyclic_array */
  unsigned lr_order; /** Approximation of the resonator wake up to this order */
  double * lr_wake; /** Wake pot. at m*delta_bucket and its (order) drivatives of all lr
  resonators */
  
  active_HC_t * active_HC;
  unsigned active_HC_size;
}
ring_t;


/**
 * \brief Availables models of electron bunch
 */
typedef enum
{
  MODEL_WEAK = 0,
  MODEL_STRONG = 1
}
model_t;

/**
 * \brief Tracking parameters
 */
typedef struct tracking
{
  int TrackPlane[3];
  long int NrevTot; /**< Total number of revolutions */
  long int NrevMon; /**< Turn interval between beam diagnostics */
  long int NrevPotentialsOut; /**< Start turn for writing potentials to output */
  
  model_t BunchModel; /** Flag for bunch models, MODEL_WEAK = 0, MODEL_STRONG = 1 */
  
  int EnableRW_short; /**< Activate the Resistive Wall self-force */
  int EnableRW_short_LON; /**< Activate the Resistive Wall self-force in the longitudinal plane */
  int EnableRW_long; /**< Activate the Resistive Wall long range force (between different bunches) */
  int EnableFBII; /**< Activate the fast beam ion interaction */
  int EnableQuantum; /**< Activate quantume excitation and radiation damping */
  int EnableAmpinv_out; /**< Activate output of the amplitude-invariant every NrevMon turns */
  int EnableActiveHC; /**< Activate active harmonic cavity in optics transformation */
  int EnableIdealHC;
  int EnableDiffCurr;
  
  int Nmlt; /**< Multiturns for the long range force */
  int NmultiT; /**< MultiT = Nmlt+2 is used in the program to stock CM data for all values of Nmlt */
  int triggerRW; /**< Trigger resistive wall effect, by exaggerating its effect in the very firsts turns. YES (0) or NO (1) */
  
  filling_t filling; /**< General filling pattern */
  int * bunch_out;   /**< Bunch numbers for output files */
  int Nbunch_out;   /**< Number of output files */
  double current_ratio;  /**< If EnableDiffCurr, defines current ration Ib(even)/Ib(odd) */

  char jobtitle[80]; /**< Job title */
  char input_filename[FILENAME_MAX]; /**< Path to the input file */
  char work_path[FILENAME_MAX]; /**< Path of the working directory, where output will be writen */
  
  int scan; /**< Flag for scan options; 0: no scan, 1: (ring)current scan, 2: chroma scan, 3: Q of HC */
  long int Nscan; /**< Number of scan steps */
  long int NrevScan; /**< Number of turns per scan step, NrevTot = (Nscan + 1) * NrevScan */
  double scan_start; /**< Starting value for current / chroma */
  double scan_step; /**< Added every scan step to current / chroma */
  
  double s0; /**< Characteristic distance RW */
  double mult; /**< Selffield Model bin width is to mult * s0 */
}
tracking_t;


typedef
struct selffield_model
{
  double * Gh; /**< Transversal Greens function, (horizontal) */
  double * Gl; /**< Longitudinal Greens function */
  double * Gv; /**< Transversal Greens function, (vertical) */
  int Ncell;    /**< Number of Bins */
  double Nsigma; /**< Length of Binning as multiples of sigma_tau */
  double dT; /**< Delta Tau, bin width */
  double sigma_tau; /**< Initial given bunch length [s] */
  
  double RsisL, aindL; /** Purely resistive and purely inductive impedance (LON) */
  double RsisV, aindV; /** Purely resistive and purely inductive impedance (VER) */
  double RsisH, aindH; /** Purely resistive and purely inductive impedance (HOR) */
  
  int PlaneL; /**< To label if resonator in this plane exists AND plane is tracked */
  int PlaneH; /**< To label if resonator in this plane exists AND plane is tracked */
  int PlaneV; /**< To label if resonator in this plane exists AND plane is tracked */

  int ImportL; /**< To label if the import data in this plane exists AND plane is tracked **/
  double * importT;
  double * importW;
}
selffield_model_t;


/**
 * \brief Particle with position, slope
 */
typedef struct particle
{
  vector_t pos;
  vector_t slope;
}
particle_t;



/**
 * \brief Statistics on a bunch in SI units
 */
typedef struct bunch_stats
{
  vector_t pos; /**< Mean position of particles in bunch */
  vector_t slope; /**< Mean slope of particles in bunch */
  vector_t pos_sigma; /**< Root mean square value for position */
  vector_t slope_sigma; /**< Root mean square value for slope */

  // SNAPSHOT not averaged over time
  vector_t ampinv; /**< Invariant amplitude, on particules */
  vector_t ampinv_cm; /**< Invariant amplitude, on center of mass */
  
  /* To calculate the mean later*/
  vector_t sum_pos; /**< Sum of particle position over NrevMon turns */
  vector_t sum_slope; /**< Sum of particle slope over NrevMon turns */
  vector_t sum_pos_sigma; /**< Sum of particle sigma over NrevMon turns */
  vector_t sum_slope_sigma; /**< Sum of particle slope sigma over NrevMon turns */
  
  /* To calculate the variance later: 1/(n-1) * {sum(x_i^2) - 1/n * sum(x_i)^2})*/
  vector_t sum_pos_sqr; /**< Sum of (bunch position)^2 over NrevMon turns */
  vector_t sum_slope_sqr; /**< Sum of (bunch slope)^2 over NrevMon turns */
  vector_t sum_pos_sigma_sqr; /**< Sum of (bunch sigma)^2 over NrevMon turns */
  vector_t sum_slope_sigma_sqr; /**< Sum of (bunch slope sigma)^2 over NrevMon turns */
  
} bunch_stats_t;


/**
 * \brief Bunch weak model (general case, with macro-particles)
 */
typedef struct weak_bunch
{
  double Ib; /**< Bunch current [A] */
  unsigned Np; /**< Number of macro-particle in bunch */
  double qp; /**< Charge of each macro-particles [C] */
  
  particle_t * particles; /**< Array with position and slope of each macroparticle */
  
  bunch_stats_t stats; /**< Statistics on a bunch: mean, sigma, ampinv, average over NrevMon turns */
  
  int kb; /**< Bucket number */
  int kb_out;       /**< Activates statistic output if kb ==  bunch_out[i] */
  
  int N_trash_low; /**< Number of particles with pos. below selffield binning range */
  int N_trash_high; /**< Number of particles with pos. above selffield binning range */
  
  /* from sbtransf.c */
  /*
  double xx0[50], xp0[50], zz0[50], yp0[50], xtau0[50], xeps0[50];
  double xx[50], xp[50], zz[50], yp[50], xtau[50], xeps[50];
  */
}
weak_bunch_t;



#endif /* MBTRACK_TYPES_H */

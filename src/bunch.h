#ifndef MBTRACK_BUNCH_H
#define MBTRACK_BUNCH_H

/**
 * @file
 * Module bunch
 */

#include <stdio.h>
#include <stdbool.h>
#include "def.h"
#include "types.h"

/* Global variables */
extern ring_t ring;
extern tracking_t track;
static int iseed0=0;

/**
 * \brief Parameters for bunch macroparticule model
 */
typedef struct bunch_macroparticle_model
{
  int Np;
  double fNp;
  int iseed[3];
  int nGen[3];
  char gendst_datafile[3][FILENAME_MAX];

  double sgmatau;
  
  /* Offset [m] */
  vector_t pos_offset;
  vector_t slope_offset;
  
  /* Sgm [m] */
  vector_t pos_sgm;
  vector_t slope_sgm;

  int modeHT;
  bool mode_excitation;
  
  /* Offset and Sgm, with various units: [m] or [mm] or [nm] */
  /* Set by read-in, if negative: calculated by macrop_model_setup_parameters */
  double xtauCM_offset, xepsCM_offset, sgm_xtau,    sgm_xeps;
  double xxCM_offset,   xpCM_offset,   sgm_xx,      sgm_xp;
  double zzCM_offset,   zpCM_offset,   sgm_zz,      sgm_zp;
}
bunch_macroparticle_model_t;

/**
 * \brief Parameters for bunch strong (slice based) distribution 
 */
typedef struct bunch_strong_distribution
{
  int Nslc;
/*  int Ncellmax;*/
  double maxTau, dTau; /* USED ? */
  double deltT; /**< Width [sec] of the bunch slicing */
  int    Nofe[50][50][50]; /**< Number of Macro-Electrons in grid */
}
bunch_strong_distribution_t;

#define NCELL_MAX 50 /**< Maximum possible value of Ncell, constant equal to 50 (old icellmax) */

/**
 * \brief Parameters for bunch weak (macro-particle based) distribution
 * : Wake force is calculated for bins. Every macroparticle belongs to one bin -> stored in mapcell
 */
typedef struct bunch_weak_distribution
{
  double maxTau; /**< Wake force is calculated for bins, last bin at sigma_tau * maxTau */
  double dTau; /**< Wake force is calculated for bins, bin width as multiples of sigma_tau */
  int mapcell[10000]; /**< List of bin numbers for each macroparticle */
  int mapcell2jp[50][10000];
}
bunch_weak_distribution_t;

/**
 * \brief Electron beam and filling parameters
 */
typedef struct e_beam
{
  int    nfFill[500]; /**< Filling by bucket. 0: empty,  1: filled  */
  int    Nbunch; /**< Total number of filled bunches in the ring */
  
  bunch_strong_distribution_t distrib;
  
  bunch_weak_distribution_t w_distrib;
  
}
e_beam_t;


/**
 * Print actual bunch_stats
 */
int
bunch_stats_fprintf(FILE * fp, long int rev, bunch_stats_t * stats, double scan);

/**
 * Print bunch_stats header
 */
void
worker_weak_stat_output(FILE * fp, const bunch_stats_t * bstats, const double Ib);


/**
 * \brief Bunch strong model (sliced)
 */
typedef struct bunch_strong
{
  double xx[50]; /**< X-axis position array */
  double xp[50]; /**< X-axis slope array */
  double zz[50]; /**< Z-axis position array */
  double zp[50]; /**< Z-axis slope array */
  double xtau[50]; /**< Longitudinal position array */
  double xeps[50]; /**< Longitudinal slope array */
  int npop[50]; /**< Particules population */
  
  bunch_stats_t stats;
  
  /* from sbtransf.c */
  /*
  double xx0[50], xp0[50], zz0[50], yp0[50], xtau0[50], xeps0[50];
  double xx[50], xp[50], zz[50], yp[50], xtau[50], xeps[50];
  */
}
bunch_strong_t;


/**
 * \brief Center of mass information for every electron bunches, with history
 */
typedef struct bunch_CM_history_weak
{  
  vector_t * ampinv; /**< Invariant amplitude, on particules */
  vector_t * ampinv_cm; /**< Invariant amplitude, on center of mass */
  
  // only used for weakstrong:
  vector_t * cm; /**< Bunch center of mass */
  vector_t * cm_slope; /**< Bunch center of mass */
  vector_t * rms; /**< Root mean square value for position */
  vector_t * rms_slope; /**< Root mean square value for slope */
}
bunch_CM_history_weak_t;


typedef struct bunch_CM_history_strong
{
  // For RW longrange: 
  // ! history over more than 30 turns not possible, Nbunch not more than 500 possible ! //
  double zzhist[500][30]; /**< History of Z-axis positions */
  double zphist[500][30]; /**< History of Z-axis slopes */
  double dataVV[15000];
  double dataVP[15000];
  
  double xxhist[500][30]; /**< History of X-axis positions */
  double xphist[500][30]; /**< History of X-axis slopes */
  double dataXX[15000];
  double dataXP[15000];
  
// for ampinv (phase space ellipse) output:
// old:  
//   vector_t ampinv[100][500]; /** 100: max of rev/NrevMon; 500 max of bunches in the ring */
//   vector_t ampinv_cm[100][500]; /**< Invariant amplitude, on center of mass */
  
  vector_t * ampinv; /**< Invariant amplitude, on particules */
  vector_t * ampinv_cm; /**< Invariant amplitude, on center of mass */
  
  // only used for weakstrong:
  vector_t * cm; /**< Bunch center of mass */
  vector_t * cm_slope; /**< Bunch center of mass */
  vector_t * rms; /**< Root mean square value for position */
  vector_t * rms_slope; /**< Root mean square value for slope */
}
bunch_CM_history_strong_t;



/**
 * Create a bunch: initialize number and charge of particles,
 * and allocate a particle vector
 *
 * @param bunch Bunch to be initialized (result variable)
 * @param Np Number of particle in bunch
 * @param qp Charge of each particle
 * @param kb Bucket number
 * 
 */
bool
weak_bunch_create(weak_bunch_t * bunch, const tracking_t * track, const unsigned Np, const double Ib,
                  const int kb, bool allocate);

/**
 * Destroy a bunch: free memory used for particles, and reset values
 *
 * Note: call this after weak_bunch_create.
 */
bool
weak_bunch_destroy(weak_bunch_t * bunch);

/**
 * Generate distribution of a bunch
 *
 * Note: call this after weak_bunch_create.
 *
 * @param bunch Bunch (result variable)
 */
int
weak_generate_bunch_distribution(weak_bunch_t * bunch,
                                 const bunch_macroparticle_model_t bunchModel,
                                 const int TrackPlane[3]);

/**
 * Writeout weak bunch distribution to a file
 */
bool
weak_writeout_bunch_distribution(const weak_bunch_t * bunch,
                                 const bunch_macroparticle_model_t bunchModel,
                                 const plane_t iplane,int generated);


/**
 * Update bunch current
 *
 * @param bunch Bunch (result variable)
 * @param Inew new BUNCH current
 */
void
weak_bunch_update_current(weak_bunch_t * bunch, double Inew);


/**
 * Replaces bunch_writeout_mean_ampinv_pos,
 * Writes out only average ampinv over bunches, per turn, to a file
 */
int
weak_bunch_writeout_mean_ampinv(const long int NrevTot, const long int NrevMon,
                               const bunch_CM_history_weak_t * CMhist,
                               const e_beam_t ebeam, const plane_t plane, const double * scan_val_hist);

/**
 * Initialize structure bunch_CM_history_t
 *
 * @param bCM Bunches center of mass (result variable)
 */
void bunch_CM_history_weak_init(bunch_CM_history_weak_t * bCM, int Nout);
void bunch_CM_history_strong_init(bunch_CM_history_strong_t * bCM, int Nout);

void bunch_CM_history_weak_destroy(bunch_CM_history_weak_t * bCM);
void bunch_CM_history_strong_destroy(bunch_CM_history_strong_t * bCM);

/**
 * Push a new center of mass value in historic, removing the older one
 */
void
bunch_CM_history_append(bunch_CM_history_strong_t * bCM,
                        int kb,
                        double xxCM, double xpCM,
                        double zzCM, double zpCM);

/**
 * Update bunch_CM_history's data.. arrays, using ..hist 2D arrays.
 */
int sortCMdata(int in, bunch_CM_history_strong_t * bCM);

/**
 * Update bunch_CM_history's ..hist 2D arrays, using data.. arrays.
 * This does the opposite of sortCMdata().
 */
int resortCMdata(int in, bunch_CM_history_strong_t * CMhist);

/**
 * Write out ampinv at one turn to a file
 */
int bunch_weak_writeout_ampinv(const long int rev, const long int irevmon,
                          const bunch_CM_history_weak_t * CMhist,
                          const e_beam_t ebeam);
int bunch_strong_writeout_ampinv(const long int rev, const long int irevmon,
                          const bunch_CM_history_strong_t * CMhist,
                          const e_beam_t ebeam);

/**
 * Write out average ampinv over bunches, per turn, to a file
 */
int
bunch_weak_writeout_mean_ampinv_pos(const long int NrevTot, const long int NrevMon,
                               const bunch_CM_history_weak_t * CMhist,
                               const e_beam_t ebeam, const plane_t plane);
int
bunch_strong_writeout_mean_ampinv_pos(const long int NrevTot, const long int NrevMon,
                               const bunch_CM_history_strong_t * CMhist,
                               const e_beam_t ebeam, const plane_t plane);

/**
 * Update center of mass information of a bunch, for one plane
 */
void bunch_strong_updateCM(bunch_strong_t * bunch, int slices, plane_t plane);

/**
 * Sort information on center of mass
 */
int sortCMdata(int in, bunch_CM_history_strong_t * bCM);

/**
 * Compute statistics (average, RMS, and invariant amplitude) of a bunch
 *
 * @param bunch
 *    Bunch in strong model (result variable)
 * @param slices
 *    Number of slices in bunch
 */
int bunch_strong_distribution_statistics(bunch_strong_t * bunch, int slices);

/**
 * Generate electron bunch distribution, based on electron center of mass of
 * a number of longitudinal slices
 *
 * eCMS : electron center of mass
 *
 * @param ring Ring parameters
 * @param macrop_model Bunch Macroparticule model parameters
 * @param ebeam Electron beam (result variable)
 * @param bCM Bunches center of mass (result variable)
 * @param fbii_grid Grid for FBII interaction (result variable)
 */
int
strong_generate_eCMs(const bunch_macroparticle_model_t MPmodel,
                     e_beam_t * ebeam,
                     bunch_strong_t * bunch0,
                     int writeout);

/**
 * Generation electron bunch distribution, based on numerous macro particules
 * 
 * @param nGen
 * @param iplane
 * @param uu0
 * @param up0
 * @param ring Ring parameters
 * @param macrop_model Bunch Macroparticule model parameters
 * @return
 *    Always return 1
 */
int generate_initial_sbunch_distribution2(int nGen, int writeout, plane_t iplane,
                                          double uu0[30000], double up0[30000],
                                          const ring_t ring,
                                          const bunch_macroparticle_model_t MPmodel);

/**
 * Preparation of bunch distributions for all bunches in the ring
 * prior to performing the multibunch tracking.
 *
 * The single bunch distribution will be identical for all single bunches,
 * and is given by "generate_initial_sbunch_distribution2()"
 *
 * @return
 *   Always return 1
 */
 int generate_initial_general_distribution2(bunch_strong_t * bunches,
                                           const bunch_strong_t * bunch0,
                                           const bunch_macroparticle_model_t MPmodel,
                                           const ring_t ring,
                                           const e_beam_t ebeam);

/**
 * \brief Grid parameters. Used for ions-electrons interaction
 */
typedef struct grid
{
  int Ngrid;
  int Ndivs;
  int offst;
  
  /* grid size */
  double agLOCx, agLOCz, agABSx, agABSz;
} grid_t;

/**
 * Initialize grid parameters
 */
bool
fbii_grid_setup_parameters(grid_t * grid,
                           e_beam_t * ebeam,
                           const bunch_macroparticle_model_t MPmodel);

/**
 * Generate grid for beam-beam interaction
 */
void
grid_generate(grid_t * grid,
              e_beam_t * ebeam,
              bunch_strong_t * sbunch0,
              const weak_bunch_t * wbunch0,
              const bunch_macroparticle_model_t MPmodel);

/**
 * Print grid
 */
int fprint_grid(FILE *, const grid_t);

/**
 * Random number generations following normal distribution
 */
double c_fnorm(const int iseed);

/**
 * Random number generations following normal distribution, hopefully faster
 */
double box_muller(const int iseed);

/**
 * Print a 2D matrix
 */
void fprintf_matrix(FILE * fp, const unsigned lines, const unsigned columns, int matrix[][50]);

#endif /* MBTRACK_BUNCH_H */

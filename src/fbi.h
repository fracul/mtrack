#ifndef MBTRACK_FBI_H
#define MBTRACK_FBI_H

/**
 * @file
 * Module for FBII interation
 */

#include <stdbool.h>
#include <stdio.h>
/*#include "types.h"*/
#include "bunch.h"
#include "physic.h"

/**
 * \brief Fast Beam Interaction model
 */
typedef struct fbii_model
{
  double me, echarge, mu0, epsilon0, re, mpOme;
  
  double Pave, Pu20, Pm, dm, sigmam, Lu20, Ldmsgm;
  double Aion, factBE;
  double factF;
}
fbii_model_t;

/**
 * Initialize values for the Fast Beam Interaction model
 *
 * @param FBIImodel
 *    FBII model (result variable)
 */
bool fbii_model_init(fbii_model_t * FBIImodel);

/**
 * \brief TODO description
 */
typedef struct ion_model
{
  /*int     Nslc, mapslice2i[50][50][50];*/
  int    Nions;
  int    nFio, kbFio;
/*  int    Nofe[50][50][50];     */
/*  double maxTauFBI, dtauFBI; // moved to distribFBI, part of e_beam_t */
}
ion_model_t;

/**
 * \brief Describe the ion beam
 */
typedef struct ion_beam
{

  double NriFio[2500000];
  double yyFio[2500000], ydFio[2500000];
  double xxFio[2500000], xdFio[2500000];
  
  ion_model_t model;
  bunch_weak_distribution_t distrib; /* TODO setup */
}
ion_beam_t;

/**
 * \brief Describe the ion beam at the individual bunch level
 */
typedef struct ion_beam_local
{
  double Nri[5000];
  double Ey[5000], yyi0[5000], ydi0[5000];
  double Ex[5000], xxi0[5000], xdi0[5000];
  
  /* Map cells in three dimensional space to Ex/Ey index */
  int mapslice2i[50][50][50];
}
ion_beam_local_t;

/**
 * Generate ions created by collision between a single electron bunch and 
 * residual gas
 *
 * @param kb
 *   Bucket identifier (a bucket may contain no bunch)
 */
int fbi_weakstrong_bunch(const bunch_macroparticle_model_t MPmodel,
                         ion_model_t * ionmodel,
                         const fbii_model_t FBIImodel,
                         const e_beam_t ebeam, ion_beam_local_t * ionbeam,
                         bunch_strong_t * bCM, int kb);

/**
 * Simulate two-beam interaction (ions beam and electrons beam)
 * 
 * @param ib
 *   Bunch identifier
 * @param kb
 *   Bucket identifier (a bucket may contain no bunch)
 */
int fbi_weakstrong_train(const bunch_macroparticle_model_t MPmodel,
                         const fbii_model_t FBIImodel,
                         const e_beam_t ebeam,
                         ion_beam_local_t * ionlocalbeam, ion_beam_t * ionbeam,
                         bunch_strong_t * bCM,
                         int ib, int kb);

/**
 * Interpolate electric force between ions and electrons from already computed
 * electric field (on a grid).
 *
 * @params
 *   Grid cell identifier
 * @param xx
 *   xx represents the horizontal coordiate in the absolute axis
 *   at which, the E-field shall be calculated
 * @param yy
 *   xx represents the vertical coordiate in the absolute axis
 *   at which, the E-field shall be calculated
 * @param ionbeam
 *   Local ion beam structure
 * @param Etr
 *   Results
 */
int get_linearly_interpolated_Efields(int icell, double xx, double yy,
                                      const ion_beam_local_t * ionbeam,
                                      double Etr[2]);

/**
 * Write ions center of mass
 *
 * @param fpoH
 *   File for horizontal position informations
 * @param fpoV
 *   File for vertical position informations
 * @param in
 *   Turn number
 * @param ib
 *   Bunch identifier
 */
int writeoutCMImotions(const char * output_path, int in, int ib,
                       const ion_beam_t * ionbeam, const e_beam_t * ebeam,
                       const bunch_macroparticle_model_t MPmodel);

/**
 * Write motion information on centers of mass motion.
 * 
 * Use PMdata (post mortem data) format
 *
 * Aurelien Bence (15 June 2011)
 *
 * @param fpoV
 *   File for vertical motion informations
 * @param in
 *   Turn number
 */
int writeoutVCMmotionsPM2(const char * output_path, int in,
                          const bunch_strong_t * bunches,
                          const e_beam_t * ebeam);

#endif /* MBTRACK_FBI_H */

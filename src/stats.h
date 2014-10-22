/**
 * @file
 * Statistics module, TO BE REMOVED
 */

#ifndef MBTRACK_STATS_H
#define MBTRACK_STATS_H

#include "types.h"

/* TODO REMOVE THIS MODULE, use bunch_stats instead if necessary */

/**
 * \brief Statistics on the electron beam
 */
typedef struct e_beam_stats
{
  double uu_ave[3], uu_rms[3], up_ave[3], up_rms[3];
  double xtau_ave[500], xtau_rms[500], xeps_ave[500], xeps_rms[500];
  double zz_ave[500],   zz_rms[500],   zp_ave[500],   zp_rms[500];
  double xx_ave[500],   xx_rms[500],   xp_ave[500],   xp_rms[500];
  double ampinv[3], ampinvCM[3];
  double ainvH[500], ainvV[500], ainvL[500];
  double ainvCMH[500], ainvCMV[500], ainvCML[500];
}
e_beam_stats_t;

/**
 * \brief Print vertical rms beamsize written out in [mm]
 */
int e_beam_rms_size_fprint(FILE * fp, const e_beam_stats_t * stats, const ring_t ring);

#endif /* MBTRACK_STATS_H */

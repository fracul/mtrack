/**
 * @file
 * Resistive wall module
 */

#ifndef MBTRACK_RW_H
#define MBTRACK_RW_H

#include "types.h"
#include "bunch.h"

/**
 * Short-range bunch transformation (strong model)
 */
int transformSbunch_shortrange(int in, unsigned int kb,
                               const e_beam_t ebeam,
                               const bunch_macroparticle_model_t MPmodel,
                               bunch_strong_t * bunch0,
                               bunch_strong_t * bunch);

/**
 * Long-range resistive wall bunch transformation (strong model)
 */
int
transformSbunch_longrange(int in, unsigned int kb, const plane_t plane,
                          const e_beam_t ebeam,
                          const bunch_macroparticle_model_t MPmodel,
                          const bunch_CM_history_strong_t * bCMhist,
                          bunch_strong_t * bCM);

/**
 * Matrix based optical transformation of a bunch (strong model)
 */
int
transformSbunch_optic(int in, int kb,
                      const e_beam_t ebeam, bunch_strong_t * bCM);

#endif /* MBTRACK_RW_H */

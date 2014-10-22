/**
 * @file
 * Input module
 */
#ifndef MBTRACK_INPUT_H
#define MBTRACK_INPUT_H

#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include "types.h"
#include "fbi.h"
#include "transform_weak.h"

#define dbhor         0
#define dbver         1
#define dbboth_planes -1

/* Input format */

bool read_input(char filename[FILENAME_MAX], ring_t * ring, tracking_t * track,
                e_beam_t * ebeam,
                bunch_macroparticle_model_t * macrop_model,
                selffield_model_t * SelfFieldModel);

bool read_input_file(FILE * fp, ring_t * ring, tracking_t * track,
                    bunch_macroparticle_model_t * macrop_model,
                    bunch_strong_distribution_t * bunchStrongD);

bool read_conf_file(FILE * fp, ring_t * ring, tracking_t * track,
                    bunch_macroparticle_model_t * macrop_model,
                    bunch_strong_distribution_t * bunchStrongD,
                    selffield_model_t * SelfFieldModel);

extern bool config_get_fprint;

/**
 * Check if config have a section, and ensure that section is unique
 */
bool config_have_section(FILE * fp, const char * section);

inline
bool config_get_int(FILE * fp, const char * section, const char * identifier, int * r);

inline
bool config_get_long_int(FILE * fp, const char * section, const char * identifier, long int * r);

inline
bool config_get_double(FILE * fp, const char * section, const char * identifier, double * r);

inline
bool config_get_str(FILE * fp, const char * section, const char * identifier,
                    char * r);

/* Checking */

bool check_parameters(ring_t * ring, tracking_t * track, e_beam_t * ebeam);

/* Pretty printing */

int fprint_ring(FILE * fp, const ring_t ring);

int fprint_tracking(FILE * fp, const tracking_t track);

int fprintf_bunch_macroparticle_model(FILE * fp, const bunch_macroparticle_model_t model, const int TrackPlane[3]);

/* electron beam */

bool e_beam_setup(tracking_t * track, ring_t * ring, e_beam_t * ebeam);

int fprint_e_beam(FILE * fp, const ring_t ring, const e_beam_t ebeam);

int e_beam_print(const ring_t ring, const e_beam_t ebeam);

bool setup_ring_parameters(ring_t * ring);

bool setup_tracking_parameters(ring_t * ring, tracking_t * track, 
                               selffield_model_t * SelfFieldModel,
                               bunch_macroparticle_model_t * macrop_model);

bool macrop_model_setup_parameters(const ring_t ring, const tracking_t tracking, 
                           bunch_macroparticle_model_t * macrop_model);

int fprint_parameters(FILE * fp, const ring_t ring,
                      const bunch_macroparticle_model_t macrop_model);

int
fprint_resonators(FILE * fp, const ring_t ring, const selffield_model_t SelfFieldModel, const bunch_macroparticle_model_t macrop_model);


int ring_destroy(ring_t * ring);

int
tracking_destroy(tracking_t * track);



#endif /* MBTRACK_INPUT_H */

/**
 * @file
 * parameter scan while tracking
 */

#ifndef TRACKING_H
#define TRACKING_H

#include <stdbool.h>
#include <stdlib.h>
#include "types.h"
#include "bunch.h"

int
tracking_scan_step_manager(long int rev, double * scan_val, const tracking_t * track, ring_t * ring, e_beam_t * ebeam, weak_bunch_t * bunches);

int
tracking_scan_step_worker(long int rev, double * scan_val, const tracking_t * track, ring_t * ring, e_beam_t * ebeam, weak_bunch_t * bunch);

#endif /* TRACING_H */
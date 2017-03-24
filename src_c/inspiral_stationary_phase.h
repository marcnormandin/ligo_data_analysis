/*
 * stationary_phase.h
 *
 *  Created on: Mar 2, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_STATIONARY_PHASE_H_
#define SRC_C_STATIONARY_PHASE_H_

#include <stddef.h>
#include <gsl/gsl_complex.h>

#include "inspiral_chirp.h"
#include "strain.h"

typedef struct stationary_phase_s {
	size_t 			len;
	gsl_complex		*spa_0;
	gsl_complex		*spa_90;
} stationary_phase_t;

stationary_phase_t* SP_malloc(size_t size);

void SP_free(stationary_phase_t *sp);

double SP_g(double f_low, double f_high, chirp_time_t *chirp, strain_t *interp_strain);

void SP_compute(double coalesce_phase, double time_delay,
		chirp_time_t *chirp, strain_t *interp_strain,
		double f_low, double f_high,
		stationary_phase_t *out_sp);

void SP_save(char *filename, strain_t *interp_strain, stationary_phase_t *sp);

#endif /* SRC_C_STATIONARY_PHASE_H_ */

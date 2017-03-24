/*
 * antenna_patterns.h
 *
 *  Created on: Mar 2, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_ANTENNA_PATTERNS_H_
#define SRC_C_ANTENNA_PATTERNS_H_

#include "sky.h"
#include "detector.h"

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

/* GW Interferometer Antenna Patterns */
typedef struct antenna_patterns_s {
	/* polarization-independent antenna patterns */
	double u, v;

	/* polarization-dependent antenna patterns */
	double f_plus, f_cross;

} antenna_patterns_t;

typedef struct antenna_patterns_workspace_s {
	gsl_vector* n_hat;
	gsl_vector* ex_i;
	gsl_vector* ey_j;
	gsl_matrix* epsilon_plus;
	gsl_matrix* epsilon_cross;
	gsl_matrix* temp;
	gsl_matrix* wp_matrix;

} antenna_patterns_workspace_t;

antenna_patterns_workspace_t* antenna_patterns_workspace_alloc();
void antenna_patterns_workspace_free( antenna_patterns_workspace_t* workspace );

int antenna_patterns(detector_t *d, sky_t *sky, double polarization_angle,
		antenna_patterns_workspace_t *workspace, antenna_patterns_t *ant);

#endif /* SRC_C_ANTENNA_PATTERNS_H_ */

/*
 * antenna_patterns.h
 *
 *  Created on: Mar 2, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_ANTENNA_PATTERNS_H_
#define SRC_C_ANTENNA_PATTERNS_H_

#include "sky.h"

/* GW Interferometer Antenna Patterns */
typedef struct antenna_patterns_s {
	/* polarization-independent antenna patterns */
	double u, v;

	/* polarization-dependent antenna patterns */
	double f_plus, f_cross;

} antenna_patterns_t;

int antenna_patterns(const char *iid, sky_t *sky, double polarization_angle, antenna_patterns_t *ant);

#endif /* SRC_C_ANTENNA_PATTERNS_H_ */

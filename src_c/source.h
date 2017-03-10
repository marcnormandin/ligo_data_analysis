/*
 * source.h
 *
 *  Created on: Mar 2, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_SOURCE_H_
#define SRC_C_SOURCE_H_

#include "sky.h"

typedef struct {
	double 	m1;
	double 	m2;
	double 	time_of_arrival;
	sky_t	sky;

	double	polarization_angle;
	double	coalesce_phase;
	double	inclination_angle;

	double snr;
} source_t;

void Print_Source(source_t* source);

void Load_Source(source_t* source);

#endif /* SRC_C_SOURCE_H_ */

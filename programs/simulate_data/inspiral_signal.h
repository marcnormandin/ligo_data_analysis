/*
 * datagen.h
 *
 *  Created on: Feb 22, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_DATAGEN_H_
#define SRC_C_DATAGEN_H_

#include <gsl/gsl_rng.h>

#include "inspiral_signal.h"
#include "spectral_density.h"
#include "detector.h"
#include "detector_network.h"

#include "sky.h"

#include "strain.h"
typedef struct source_s {
	double 	m1;
	double 	m2;
	double 	time_of_arrival;
	sky_t	sky;

	double	polarization_angle;
	double	coalesce_phase;
	double	inclination_angle;

	double snr;
} source_t;


strain_half_fft_t* inspiral_template_half_fft(double f_low, double f_high, size_t num_time_samples, detector_t *det, source_t *source);

network_strain_half_fft_t* inspiral_template_unscaled(double f_low, double f_high, size_t num_time_samples, detector_network_t *net, source_t *source);

/* This is the main routine that generates the inspiral template. The memory must be freed. */
network_strain_half_fft_t* inspiral_template(double f_low, double f_high, size_t num_time_samples, detector_network_t *net, source_t *source);

void Source_print(source_t* source);

void Source_load_testsource(source_t* source);

#endif /* SRC_C_DATAGEN_H_ */

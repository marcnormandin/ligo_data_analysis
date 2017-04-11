/*
 * datagen.h
 *
 *  Created on: Feb 22, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_DATAGEN_H_
#define SRC_C_DATAGEN_H_

#include <gsl/gsl_rng.h>

#include "inspiral_source.h"
#include "signal.h"
#include "spectral_density.h"
#include "detector.h"
#include "detector_network.h"

typedef struct inspiral_signal_half_fft_s {
	/* length of the full fft = num_time_samples */
	size_t full_len;

	/* length of the template, which is only half of the full fft */
	size_t half_fft_len;
	gsl_complex *half_fft;

} inspiral_template_half_fft_t;

inspiral_template_half_fft_t* inspiral_template_half_fft_alloc(size_t num_time_samples);
void inspiral_signal_half_fft_free(inspiral_template_half_fft_t *signal);

inspiral_template_half_fft_t* inspiral_template_half_fft(double f_low, double f_high, size_t num_time_samples, detector_t *det, source_t *source);

inspiral_template_half_fft_t** inspiral_template_unscaled(double f_low, double f_high, size_t num_time_samples, detector_network_t *net, source_t *source);

/* This is the main routine that generates the inspiral template. The memory must be freed. */
inspiral_template_half_fft_t** inspiral_template(double f_low, double f_high, size_t num_time_samples, detector_network_t *net, source_t *source);

#endif /* SRC_C_DATAGEN_H_ */

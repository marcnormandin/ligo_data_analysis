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
#include "strain.h"
#include "detector.h"
#include "detector_network.h"

typedef struct inspiral_signal_half_fft_s {
	/* length of the full fft */
	size_t full_len;

	/* length of the template, which is only half of the full fft */
	size_t half_fft_len;
	gsl_complex *half_fft;

} inspiral_signal_half_fft_t;

inspiral_signal_half_fft_t* inspiral_signal_half_fft_alloc(size_t full_len);
void inspiral_signal_half_fft_free(inspiral_signal_half_fft_t *signal);

inspiral_signal_half_fft_t* inspiral_signal_half_fft(double f_low, double f_high, detector_t *det, strain_t *half_strain, source_t *source);

/* This is the main routine that simulates the data
   The memory must be freed. */
inspiral_signal_half_fft_t** simulate_inspiral(gsl_rng *rng, double f_low, double f_high, detector_network_t *net, strain_t *half_strain, source_t *source);

inspiral_signal_half_fft_t** simulate_inspiral_unscaled(double f_low, double f_high, detector_network_t *net, strain_t *half_strain, source_t *source);

#endif /* SRC_C_DATAGEN_H_ */

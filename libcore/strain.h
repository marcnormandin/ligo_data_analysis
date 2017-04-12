/*
 * strain.h
 *
 *  Created on: Apr 12, 2017
 *      Author: marcnormandin
 */

#ifndef LIBCORE_STRAIN_H_
#define LIBCORE_STRAIN_H_

#include <stddef.h>
#include <gsl/gsl_complex.h>

typedef struct strain_half_fft_s {
	/* length of the full fft = num_time_samples */
	size_t full_len;

	/* length of the template, which is only half of the full fft */
	size_t half_fft_len;
	gsl_complex *half_fft;

} strain_half_fft_t;

typedef struct network_strain_half_fft_s {
	size_t num_strains;
	size_t num_time_samples;

	strain_half_fft_t **strains;

} network_strain_half_fft_t;


strain_half_fft_t* strain_half_fft_alloc(size_t num_time_samples);
void strain_half_fft_free(strain_half_fft_t *strain);

network_strain_half_fft_t* network_strain_half_fft_alloc(size_t num_strains, size_t num_time_samples);
void network_strain_half_fft_free(network_strain_half_fft_t *strains);

#endif /* LIBCORE_STRAIN_H_ */

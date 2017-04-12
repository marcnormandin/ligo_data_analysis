/*
 * strain.c
 *
 *  Created on: Apr 12, 2017
 *      Author: marcnormandin
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_complex.h>

#include "sampling_system.h"
#include "strain.h"

strain_half_fft_t* strain_half_fft_alloc(size_t num_time_samples) {
	strain_half_fft_t *signal = (strain_half_fft_t*) malloc( sizeof(strain_half_fft_t) );
	if (signal == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory for strain_half_fft_t. Aborting.\n");
		abort();
	}

	signal->full_len = num_time_samples;
	signal->half_fft_len = SS_half_size(signal->full_len);

	signal->half_fft = (gsl_complex*) malloc( signal->half_fft_len * sizeof(gsl_complex) );
	if (signal->half_fft == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory for strain_half_fft_t.half_fft. Aborting.\n");
		abort();
	}

	return signal;
}

void strain_half_fft_free(strain_half_fft_t *strain) {
	assert(strain != NULL);
	free(strain->half_fft);
	free(strain);
	strain = NULL;
}

network_strain_half_fft_t* network_strain_half_fft_alloc(size_t num_strains, size_t num_time_samples) {
	size_t i;

	network_strain_half_fft_t *network_strain = (network_strain_half_fft_t*) malloc( sizeof(network_strain_half_fft_t) );
	if (network_strain == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory for network_strain_half_fft. Aborting.\n");
		abort();
	}

	network_strain->num_strains = num_strains;
	network_strain->num_time_samples = num_time_samples;

	network_strain->strains = (strain_half_fft_t**) malloc( network_strain->num_strains * sizeof(strain_half_fft_t*) );
	if (network_strain->strains == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory for the network_strain_half_fft.strains. Aborting.\n");
	}

	for (i = 0; i < network_strain->num_strains; i++) {
		network_strain->strains[i] = strain_half_fft_alloc( network_strain->num_time_samples );
	}

	return network_strain;
}

void network_strain_half_fft_free(network_strain_half_fft_t *network_strain) {
	assert(network_strain != NULL);

	size_t i;
	for (i = 0; i < network_strain->num_strains; i++) {
		strain_half_fft_free(network_strain->strains[i]);
	}
	free(network_strain);
	network_strain = NULL;
}

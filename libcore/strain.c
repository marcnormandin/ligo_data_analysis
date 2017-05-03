#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_fft_complex.h>

#include "sampling_system.h"
#include "strain.h"

strain_half_fft_t* strain_half_fft_alloc(size_t num_time_samples) {
	strain_half_fft_t *signal = (strain_half_fft_t*) malloc( sizeof(strain_half_fft_t) );
	if (signal == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory for strain_half_fft_t. Exiting.\n");
		exit(-1);
	}

	signal->full_len = num_time_samples;
	signal->half_fft_len = SS_half_size(signal->full_len);

	signal->half_fft = (gsl_complex*) malloc( signal->half_fft_len * sizeof(gsl_complex) );
	if (signal->half_fft == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory for strain_half_fft_t.half_fft. Exiting.\n");
		exit(-1);
	}

	return signal;
}

void strain_half_fft_free(strain_half_fft_t *strain) {
	assert(strain != NULL);
	free(strain->half_fft);
	free(strain);
	strain = NULL;
}

strain_full_fft_t* strain_full_fft_alloc(size_t num_time_samples) {
	strain_full_fft_t *signal = (strain_full_fft_t*) malloc( sizeof(strain_full_fft_t) );
	if (signal == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory for strain_full_fft_t. Exiting.\n");
		exit(-1);
	}

	signal->full_len = num_time_samples;
	signal->half_fft_len = SS_half_size(signal->full_len);

	signal->full_fft = (gsl_complex*) malloc( signal->full_len * sizeof(gsl_complex) );
	if (signal->full_fft == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory for strain_full_fft_t.full_fft. Exiting.\n");
		exit(-1);
	}

	return signal;
}

void strain_full_fft_free(strain_full_fft_t *strain) {
	assert(strain != NULL);
	free(strain->full_fft);
	free(strain);
}

strain_t* strain_alloc(size_t num_time_samples) {
	strain_t *strain = (strain_t*) malloc( sizeof(strain_t) );
	if (strain == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory for strain_t. Exiting.\n");
		exit(-1);
	}

	strain->num_time_samples = num_time_samples;

	strain->samples = (double*) malloc( strain->num_time_samples * sizeof(double) );
	if (strain->samples == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory for strain_t.samples. Exiting.\n");
		exit(-1);
	}

	return strain;
}

void strain_free( strain_t* strain ) {
	assert(strain != NULL);
	free(strain->samples);
	free(strain);
}

strain_full_fft_t* strain_half_to_full(strain_half_fft_t *one_sided) {
	strain_full_fft_t *two_sided = strain_full_fft_alloc( one_sided->full_len );

	/* Form the two-sided template so we can take the inverse FFT */
	SS_make_two_sided (one_sided->half_fft_len, one_sided->half_fft, two_sided->full_len, two_sided->full_fft);

	return two_sided;
}

strain_t* strain_full_fft_to_strain( strain_full_fft_t* fft) {
	size_t j;
	gsl_fft_complex_wavetable *fft_wavetable;
	gsl_fft_complex_workspace *fft_workspace;
	fft_wavetable = gsl_fft_complex_wavetable_alloc( fft->full_len );
	fft_workspace = gsl_fft_complex_workspace_alloc( fft->full_len );

	/* GSL needs an input linear array arranged as (real, complex) pairs in sequence. */
	double *template_ifft = (double*) malloc( 2 * fft->full_len * sizeof(double) );
	if (template_ifft == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory for the template_ifft. Aborting.\n");
		abort();
	}
	for (j = 0; j < fft->full_len; j++) {
		template_ifft[2*j + 0] = GSL_REAL( fft->full_fft[j] );
		template_ifft[2*j + 1] = GSL_IMAG( fft->full_fft[j] );
	}

	/* Compute the ifft */
	gsl_fft_complex_inverse( template_ifft, 1, fft->full_len, fft_wavetable, fft_workspace );

	strain_t *strain = strain_alloc(fft->full_len);

	for (j = 0; j < fft->full_len; j++) {
		/* Only save the real part */
		strain->samples[j] = template_ifft[2*j + 0];
	}

	free(template_ifft);
	gsl_fft_complex_workspace_free( fft_workspace );
	gsl_fft_complex_wavetable_free( fft_wavetable );

	return strain;
}

network_strain_half_fft_t* network_strain_half_fft_alloc(size_t num_strains, size_t num_time_samples) {
	size_t i;

	network_strain_half_fft_t *network_strain = (network_strain_half_fft_t*) malloc( sizeof(network_strain_half_fft_t) );
	if (network_strain == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory for network_strain_half_fft. Exiting.\n");
		exit(-1);
	}

	network_strain->num_strains = num_strains;
	network_strain->num_time_samples = num_time_samples;

	network_strain->strains = (strain_half_fft_t**) malloc( network_strain->num_strains * sizeof(strain_half_fft_t*) );
	if (network_strain->strains == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory for the network_strain_half_fft.strains. Exiting.\n");
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

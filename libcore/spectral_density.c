#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_interp.h>

#include "hdf5_file.h"
#include "inspiral_stationary_phase.h" /* needed for find_index */
#include "sampling_system.h"
#include "spectral_density.h"


asd_t* ASD_alloc( size_t len ) {
	asd_t *asd;

	asd = (asd_t*) malloc (sizeof(asd_t));
	if (asd == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory in ASD_malloc. Exiting.\n");
		exit(-1);
	}

	asd->len = len;

	asd->asd = (double*) malloc( len * sizeof(double) );
	if (asd->asd == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory in ASD_malloc. Exiting.\n");
		exit(-1);
	}

	asd->f = (double*) malloc( len * sizeof(double) );
	if (asd->f == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory in ASD_malloc. Exiting.\n");
		exit(-1);
	}

	return asd;
}

void ASD_free( asd_t *asd) {
	assert(asd != NULL);

	assert(asd->asd != NULL);
	free(asd->asd);
	asd->asd = NULL;

	assert(asd->f != NULL);
	free(asd->f);
	asd->f = NULL;

	free(asd);
}

void ASD_init_from_psd( psd_t *psd, asd_t *asd) {
	assert(psd != NULL);
	assert(psd->psd != NULL);
	assert(psd->f != NULL);
	assert(asd != NULL);
	assert(asd->asd != NULL);
	assert(asd->f != NULL);

	/* The allocated structures must have the same array lengths. */
	assert( psd->len == asd->len );

	size_t i, N;

	N = psd->len;

	/* ASD is the square-root of PSD */
	for ( i = 0; i < N; i++) {
		asd->asd[i] = sqrt( psd->psd[i] );
	}

	memcpy(asd->f, psd->f, psd->len * sizeof(double));

	switch (psd->type) {
	case PSD_ONE_SIDED: asd->type = ASD_ONE_SIDED; break;
	case PSD_TWO_SIDED: asd->type = ASD_TWO_SIDED; break;
	default:
		fprintf(stderr, "Error: PSD type is invalid. Can not convert to ASD. Exiting.\n");
		exit(-1);
	}
}

psd_t* PSD_alloc( size_t len ) {
	psd_t *psd;

	psd = (psd_t*) malloc (sizeof(psd_t));
	if (psd == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory in PSD_malloc. Exiting.\n");
		exit(-1);
	}

	psd->len = len;

	psd->psd = (double*) malloc( len * sizeof(double) );
	if (psd->psd == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory in PSD_malloc. Exiting.\n");
		exit(-1);
	}

	psd->f = (double*) malloc( len * sizeof(double) );
	if (psd->f == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory in PSD_malloc. Exiting.\n");
		exit(-1);
	}

	return psd;
}

void PSD_free( psd_t *psd) {
	assert(psd != NULL);

	assert(psd->psd != NULL);
	free(psd->psd);
	psd->psd = NULL;

	assert(psd->f != NULL);
	free(psd->f);
	psd->f = NULL;

	free(psd);
}

void PSD_init_from_asd( asd_t *asd, psd_t *psd) {
	assert(asd != NULL);
	assert(asd->asd != NULL);
	assert(asd->f != NULL);
	assert(psd != NULL);
	assert(psd->psd != NULL);
	assert(psd->f != NULL);

	/* The allocated structures must have the same array lengths. */
	assert( asd->len == psd->len );

	size_t i, N;

	N = asd->len;

	/* PSD is the square of ASD */
	for ( i = 0; i < N; i++) {
		psd->psd[i] = asd->asd[i] * asd->asd[i];
	}

	memcpy(psd->f, asd->f, N * sizeof(double));

	switch (asd->type) {
	case ASD_ONE_SIDED: psd->type = PSD_ONE_SIDED; break;
	case ASD_TWO_SIDED: psd->type = PSD_TWO_SIDED; break;
	default:
		fprintf(stderr, "Error: ASD type is invalid. Can not convert to PSD. Exiting.\n");
		exit(-1);
	}
}

psd_t* PSD_load( const char *hdf_filename ) {
	assert(hdf_filename != NULL);

	size_t len_psd;

	len_psd = hdf5_get_dataset_array_length( hdf_filename, "/psd/PSD" );

	psd_t* psd = PSD_alloc ( len_psd );

	hdf5_load_array( hdf_filename, "/psd/PSD", psd->psd );
	hdf5_load_array( hdf_filename, "/psd/Freq", psd->f );

	psd->type = PSD_ONE_SIDED;

	return psd;
}

void PSD_save( const char *hdf_filename, psd_t *psd ) {
	assert(hdf_filename != NULL);
	assert(psd != NULL);

	hdf5_create_group(hdf_filename, "/psd");

	hdf5_save_array( hdf_filename, "/psd", "PSD", psd->len, psd->psd );
	hdf5_save_array( hdf_filename, "/psd", "Freq", psd->len, psd->f );
}

void ASD_save( const char *hdf_filename, asd_t *asd) {
	assert(hdf_filename != NULL);
	assert(asd != NULL);

	hdf5_create_group(hdf_filename, "/asd");

	hdf5_save_array( hdf_filename, "/asd", "ASD", asd->len, asd->asd );
	hdf5_save_array( hdf_filename, "/asd", "Freq", asd->len, asd->f );
}

/* This takes a PSD that isn't specified uniformly over frequency and returns one that is. */
psd_t* PSD_nonuniform_to_uniform(psd_t *nonuniform, size_t num_time_samples, double sampling_frequency) {
	assert(nonuniform);
	assert(sampling_frequency >= 1.0);

	size_t j;
	size_t half_size = SS_half_size( num_time_samples );

	if (nonuniform->type != PSD_ONE_SIDED) {
		fprintf(stderr, "Error. PSD_nonuniform_to_uniform only works for one-sided PSDs. Exiting.\n");
		exit(-1);
	}

	/* Set the frequencies we want to use. */
	psd_t *psd = PSD_alloc( half_size );
	psd->type = PSD_ONE_SIDED;
	SS_frequency_array(sampling_frequency, num_time_samples, psd->len, psd->f);

	/* Interpolate the PSD to the desired frequencies. */
	gsl_interp* interp = gsl_interp_alloc(gsl_interp_linear, nonuniform->len);
	gsl_interp_accel* acc = gsl_interp_accel_alloc();
	gsl_interp_init(interp, nonuniform->f, nonuniform->psd, nonuniform->len);
	for (j = 0; j < psd->len; j++) {
		double asd_interpolated = gsl_interp_eval(interp, nonuniform->f, nonuniform->psd, psd->f[j], acc);
		psd->psd[j] = asd_interpolated;
	}
	gsl_interp_accel_free(acc);
	gsl_interp_free(interp);

	return psd;
}

void PSD_flatten_edges(double f_low, double f_high, psd_t *psd) {
	assert(psd != 0);
	assert(f_low >= 0.0);
	assert(f_high >= 0.0 && f_high >= f_low);
	size_t j;

	if (psd->type != PSD_ONE_SIDED) {
		fprintf(stderr, "Error. PSD_flatten_edges only works for one-sided PSDs. Exiting.\n");
		exit(-1);
	}

	size_t f_low_index = find_index_low(f_low, psd->len, psd->f);
	size_t f_high_index = find_index_high(f_high, psd->len, psd->f);
	if (f_low_index == -1 || f_high_index == -1) {
		fprintf(stderr, "Error. PSD_flatten_edges can not work. f_low_index or f_high_index were not found. Exiting.\n");
		exit(-1);
	}

	/* Flatten the ends */
	for (j = 0; j < psd->len; j++) {
		if (j <= f_low_index) {
			psd->psd[j] = psd->psd[f_low_index];
		} else if (j >= f_high_index) {
			psd->psd[j] = psd->psd[f_high_index];
		}
	}
}

psd_t* PSD_make_suitable_for_network_analysis(psd_t *nonuniform, size_t num_time_samples, double sampling_frequency, double f_low, double f_high) {
	psd_t *psd = PSD_nonuniform_to_uniform(nonuniform, num_time_samples, sampling_frequency);
	PSD_flatten_edges(f_low, f_high, psd);
	return psd;
}

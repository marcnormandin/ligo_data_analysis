#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "hdf5_file.h"
#include "spectral_density.h"


asd_t* ASD_malloc( size_t len ) {
	asd_t *asd;
	asd = (asd_t*) malloc (sizeof(asd_t));
	asd->len = len;
	asd->asd = (double*) malloc( len * sizeof(double) );
	asd->f = (double*) malloc( len * sizeof(double) );
	return asd;
}

void ASD_free( asd_t *asd) {
	free(asd->asd);
	free(asd->f);
	free(asd);
	asd = NULL;
}

void ASD_init_from_psd( psd_t *psd, asd_t *asd) {
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
		fprintf(stderr, "Error: PSD type is invalid. Can not convert to ASD. Aborting.\n");
		abort();
	}
}

psd_t* PSD_malloc( size_t len ) {
	psd_t *psd;
	psd = (psd_t*) malloc (sizeof(psd_t));
	psd->len = len;
	psd->psd = (double*) malloc( len * sizeof(double) );
	psd->f = (double*) malloc( len * sizeof(double) );
	return psd;
}

void PSD_free( psd_t *psd) {
	free(psd->psd);
	free(psd->f);
	free(psd);
	psd = NULL;
}

void PSD_init_from_asd( asd_t *asd, psd_t *psd) {
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
		fprintf(stderr, "Error: ASD type is invalid. Can not convert to PSD. Aborting.\n");
		abort();
	}
}

psd_t* PSD_load( const char *hdf_filename ) {
	size_t len_psd;

	len_psd = hdf5_get_dataset_array_length( hdf_filename, "/psd/PSD" );

	psd_t* psd = PSD_malloc ( len_psd );

	hdf5_load_array( hdf_filename, "/psd/PSD", psd->psd );
	hdf5_load_array( hdf_filename, "/psd/Freq", psd->f );

	psd->type = PSD_ONE_SIDED;

	return psd;
}

void PSD_save( const char *hdf_filename, psd_t *psd ) {
	hdf5_create_group(hdf_filename, "/psd");

	hdf5_save_array( hdf_filename, "/psd", "PSD", psd->len, psd->psd );
	hdf5_save_array( hdf_filename, "/psd", "Freq", psd->len, psd->f );
}

void ASD_save( const char *hdf_filename, asd_t *asd) {
	hdf5_create_group(hdf_filename, "/asd");

	hdf5_save_array( hdf_filename, "/asd", "ASD", asd->len, asd->asd );
	hdf5_save_array( hdf_filename, "/asd", "Freq", asd->len, asd->f );
}

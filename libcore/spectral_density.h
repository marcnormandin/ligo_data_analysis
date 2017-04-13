#ifndef SPECTRAL_DENSITY_H_
#define SPECTRAL_DENSITY_H_

#include <stddef.h>

typedef enum {
	PSD_ONE_SIDED,
	PSD_TWO_SIDED
} PSD_TYPE;

/* Power Spectral Density */
typedef struct psd_s {
	PSD_TYPE type;
	size_t len;
	double *psd;
	double *f;

} psd_t;

typedef enum {
	ASD_ONE_SIDED,
	ASD_TWO_SIDED
} ASD_TYPE;

/* Amplitude Spectral Density.
 * This is the square-root of the PSD. */
typedef struct asd_s {
	ASD_TYPE type;
	size_t len;
	double *asd;
	double *f;

} asd_t;

psd_t* PSD_malloc( size_t len );
void PSD_free( psd_t *psd);
void PSD_init_from_asd( asd_t *asd, psd_t *psd);

psd_t* PSD_load( const char *hdf_filename );
void PSD_save( const char *hdf_filename, psd_t *psd );

asd_t* ASD_malloc( size_t len );
void ASD_free( asd_t *asd);
void ASD_init_from_psd( psd_t *psd, asd_t *asd);


#endif /* SPECTRAL_DENSITY_H_ */

#ifndef STRAIN_H_
#define STRAIN_H_

#include <stddef.h>

#include <gsl/gsl_complex.h>

typedef struct strain_half_fft_s {
	/* length of the full fft = num_time_samples */
	size_t full_len;

	/* length of the template, which is only half of the full fft */
	size_t half_fft_len;
	gsl_complex *half_fft;

} strain_half_fft_t;

typedef struct strain_full_fft_s {
	size_t full_len;
	size_t half_fft_len;
	gsl_complex* full_fft;

} strain_full_fft_t;

typedef struct strain_s {
	size_t num_time_samples;

	double *samples;

} strain_t;

typedef struct network_strain_half_fft_s {
	size_t num_strains;
	size_t num_time_samples;

	strain_half_fft_t **strains;

} network_strain_half_fft_t;


strain_half_fft_t* strain_half_fft_alloc(size_t num_time_samples);
void strain_half_fft_free(strain_half_fft_t *strain);

strain_full_fft_t* strain_full_fft_alloc(size_t num_time_samples);
void strain_full_fft_free(strain_full_fft_t *strain);

strain_t* strain_alloc(size_t num_time_samples);
void strain_free( strain_t* strain );

strain_full_fft_t* strain_half_to_full(strain_half_fft_t *one_sided);
strain_t* strain_full_fft_to_strain( strain_full_fft_t* fft);

network_strain_half_fft_t* network_strain_half_fft_alloc(size_t num_strains, size_t num_time_samples);
void network_strain_half_fft_free(network_strain_half_fft_t *strains);

#endif /* STRAIN_H_ */

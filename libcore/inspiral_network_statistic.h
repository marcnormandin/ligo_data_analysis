#ifndef SRC_C_NETWORK_ANALYSIS_H_
#define SRC_C_NETWORK_ANALYSIS_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics_double.h>

#include "detector_antenna_patterns.h"
#include "detector_network.h"
#include "inspiral_chirp_time.h"
#include "inspiral_stationary_phase.h"
#include "spectral_density.h"
#include "strain.h"

#if defined (__cplusplus)
extern "C" {
#endif

/* There is one helper per detector */
typedef struct coherent_network_helper_s {
	double w_plus_input;
	double w_minus_input;
	gsl_complex *c_plus;
	gsl_complex *c_minus;

} coherent_network_helper_t;

void CN_template_chirp_time(double f_low, double chirp_time0, double chirp_time1_5, inspiral_chirp_time_t *ct);

coherent_network_helper_t* CN_helper_alloc(size_t num_time_samples);

void CN_helper_free( coherent_network_helper_t* helper);

typedef struct coherent_network_workspace_s {
	size_t num_time_samples;

	size_t num_helpers;
	coherent_network_helper_t **helpers;

	stationary_phase_workspace_t *sp_lookup;
	stationary_phase_t *sp;

	/* temporary array that is repeatedly used.
	 * Its size must be the same as the number of frequencies.
	 */
	gsl_complex *temp_array;

	gsl_complex **terms;
	double **fs;

	double *temp_ifft;
	gsl_fft_complex_wavetable *fft_wavetable;
	gsl_fft_complex_workspace *fft_workspace;

	detector_antenna_patterns_workspace_t *ap_workspace;
	detector_antenna_patterns_t *ap;

	/* g, normalization factor */
	double *normalization_factors;

} coherent_network_workspace_t;

coherent_network_workspace_t* CN_workspace_alloc(size_t num_time_samples, detector_network_t *net, size_t num_half_freq,
		double f_low, double f_high);

void CN_workspace_free( coherent_network_workspace_t *workspace );

void CN_do_work(size_t num_time_samples, size_t f_low_index, size_t f_high_index, gsl_complex *spa, asd_t *asd, gsl_complex *whitened_data, gsl_complex *temp, gsl_complex *out_c);

void CN_save(char* filename, size_t len, double* tmp_ifft);

void coherent_network_statistic(
		detector_network_t* net,
		double f_low,
		double f_high,
		inspiral_chirp_time_t *chirp,
		sky_t *sky,
		network_strain_half_fft_t *network_strain,
		coherent_network_workspace_t *workspace,
		double *out_network_css_value,
		int *out_network_css_index,
		char *out_network_css_filename);

#if defined (__cplusplus)
}
#endif

#endif /* SRC_C_NETWORK_ANALYSIS_H_ */

/*
 * network_analysis.h
 *
 *  Created on: Mar 2, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_NETWORK_ANALYSIS_H_
#define SRC_C_NETWORK_ANALYSIS_H_

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_fft_complex.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics_double.h>

#include "inspiral_stationary_phase.h"
#include "strain.h"
#include "strain_interpolate.h"

#include "inspiral_chirp.h"
#include "inspiral_template.h"
#include "spectral_density.h"
#include "detector_network.h"
#include "detector_antenna_patterns.h"
#include "signal.h"

/* Data that the network will search for an inspiral in */
/*
typedef struct coherent_network_data_s {
	size_t len_data;
	array_t *data;

} coherent_network_data_t;
*/

/* There is one helper per detector */
typedef struct coherent_network_helper_s {
	double w_plus_input;
	double w_minus_input;
	gsl_complex *c_plus;
	gsl_complex *c_minus;

} coherent_network_helper_t;

coherent_network_helper_t* CN_helper_malloc(size_t num_time_samples);

void CN_helper_free( coherent_network_helper_t* helper);

typedef struct coherent_network_workspace_s {
	size_t num_time_samples;

	size_t num_helpers;
	coherent_network_helper_t **helpers;

	stationary_phase_lookup_t *sp_lookup;
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

	antenna_patterns_workspace_t *ap_workspace;
	antenna_patterns_t *ap;

} coherent_network_workspace_t;

coherent_network_workspace_t* CN_workspace_malloc(size_t num_time_samples, size_t num_detectors, size_t num_half_freq,
		double f_low, double f_high, asd_t *asd);

void CN_workspace_free( coherent_network_workspace_t *workspace );

void do_work(size_t num_time_samples, gsl_complex *spa, asd_t *asd, gsl_complex *whitened_data, gsl_complex *temp, gsl_complex *out_c);

void CN_save(char* filename, size_t len, double* tmp_ifft);

void coherent_network_statistic(
		detector_network_t* net,
		double f_low,
		double f_high,
		chirp_time_t *chirp,
		sky_t *sky,
		double polarization_angle,
		inspiral_template_half_fft_t **signals,
		coherent_network_workspace_t *workspace,
		double *out_val);

#endif /* SRC_C_NETWORK_ANALYSIS_H_ */

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

#include "stationary_phase.h"
#include "strain.h"
#include "strain_interpolate.h"
#include "network_analysis.h"
#include "notused_templategen.h"

#include "chirp.h"
#include "strain.h"
#include "detector_network.h"
#include "signal.h"

// There is one helper per detector
typedef struct coherent_network_helper_s {
	double w_plus_input;
	double w_minus_input;
	gsl_complex *c_plus;
	gsl_complex *c_minus;

} coherent_network_helper_t;

coherent_network_helper_t* CN_helper_malloc(size_t num_frequencies);

void CN_helper_free( coherent_network_helper_t* helper);

typedef struct coherent_network_workspace_s {
	size_t num_helpers;
	coherent_network_helper_t **helpers;

	stationary_phase_t *sp;

	// temporary array that is repeatedly used. Its size must be the same as the number of frequencies.
	gsl_complex *temp_array;

	gsl_complex **terms;
	double **fs;


	double *temp_ifft;
	gsl_fft_complex_wavetable *fft_wavetable;
	gsl_fft_complex_workspace *fft_workspace;

} coherent_network_workspace_t;

coherent_network_workspace_t* CN_workspace_malloc(size_t num_detectors, size_t len_freq);

void CN_workspace_free( coherent_network_workspace_t *workspace );

void do_work(gsl_complex *spa, strain_t *regular_strain, gsl_complex *whitened_data, gsl_complex *temp, gsl_complex *out_c);

void CN_save(char* filename, size_t len, double* tmp_ifft);

void coherent_network_statistic(
		detector_network_t* net,
		strain_t *regular_strain,
		double f_low,
		double f_high,
		chirp_time_t *chirp,
		sky_t *sky,
		double polarization_angle,
		signal_t **signals,
		coherent_network_workspace_t *workspace,
		double *out_val);

#endif /* SRC_C_NETWORK_ANALYSIS_H_ */

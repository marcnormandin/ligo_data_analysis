/*
 * coherent_new_sec2.c
 *
 *  Created on: Feb 23, 2017
 *      Author: marcnormandin
 */

#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_statistics_double.h>

#include "strain.h"
#include "strain_interpolate.h"
#include "network_analysis.h"
#include "detector.h"
#include "detector_network.h"
#include "sampling_system.h"
#include "time_delay.h"

coherent_network_helper_t* CN_helper_malloc(size_t num_frequencies) {
	size_t s;
	coherent_network_helper_t *h;

	h = (coherent_network_helper_t*) malloc( sizeof(coherent_network_helper_t) );

	/* Fixme
	 * Check that this is always valid. It may depend on nyquist frequency being present
	 */
	s = 2*num_frequencies - 2;
	h->c_plus = (gsl_complex*) malloc( s * sizeof(gsl_complex) );
	h->c_minus = (gsl_complex*) malloc( s * sizeof(gsl_complex) );
	return h;
}

void CN_helper_free( coherent_network_helper_t* helper) {
	free(helper->c_plus);
	free(helper->c_minus);
	free(helper);
}

coherent_network_workspace_t* CN_workspace_malloc(size_t num_detectors, size_t len_freq) {
	coherent_network_workspace_t * work;
	size_t i;
	size_t len_terms;
	size_t s;

	work = (coherent_network_workspace_t*) malloc(sizeof(coherent_network_workspace_t));
	work->num_helpers = num_detectors;
	work->helpers = (coherent_network_helper_t**) malloc( work->num_helpers * sizeof(coherent_network_helper_t*));
	for (i = 0; i < work->num_helpers; i++) {
		work->helpers[i] = CN_helper_malloc( len_freq );
	}

	work->sp = SP_malloc( len_freq );

	work->temp_array = (gsl_complex*) malloc( len_freq * sizeof(gsl_complex) );

	work->terms = (gsl_complex**) malloc(4 * sizeof(gsl_complex*) );
	/* Fixme */
	len_terms = 2 * len_freq - 2;
	for (i = 0; i < 4; i++) {
		work->terms[i] = (gsl_complex*) malloc( len_terms * sizeof(gsl_complex) );
	}

	work->fs = (double**) malloc( 4 * sizeof(double*) );
	for (i = 0; i < 4; i++) {
		work->fs[i] = (double*) malloc( 2 * len_terms * sizeof(double) );
	}

	work->temp_ifft = (double*) malloc( len_terms * sizeof(double) );
	s = 2 * len_freq - 2;
	work->fft_wavetable = gsl_fft_complex_wavetable_alloc( s );
	work->fft_workspace = gsl_fft_complex_workspace_alloc( s );

	work->ap_workspace = antenna_patterns_workspace_alloc();
	/* one antenna pattern structure per detector */
	work->ap = (antenna_patterns_t*) malloc (num_detectors * sizeof(antenna_patterns_t) );

	return work;
}

void CN_workspace_free( coherent_network_workspace_t *workspace ) {
	size_t i;

	for (i = 0; i < workspace->num_helpers; i++) {
		CN_helper_free( workspace->helpers[i] );
	}
	free(workspace->helpers);

	SP_free(workspace->sp);

	for (i = 0; i < 4; i++) {
		free(workspace->terms[i]);
		free(workspace->fs[i]);
	}

	free(workspace->terms);
	free(workspace->fs);
	free(workspace->temp_ifft);
	free(workspace->temp_array);
	gsl_fft_complex_workspace_free( workspace->fft_workspace );
	gsl_fft_complex_wavetable_free( workspace->fft_wavetable );

	antenna_patterns_workspace_free(workspace->ap_workspace);
	free(workspace->ap);

	free( workspace );
}

void do_work(gsl_complex *spa, strain_t *regular_strain, gsl_complex *whitened_data, gsl_complex *temp, gsl_complex *out_c) {
	size_t k;
	size_t t_index;
	size_t c_index;

	for (k = 0; k < regular_strain->len; k++) {
		temp[k] = gsl_complex_conjugate(spa[k]);
		temp[k] = gsl_complex_div_real(temp[k], regular_strain->strain[k]);
		temp[k] = gsl_complex_mul( temp[k], whitened_data[k] );

		/*out_c[k] = temp[k];*/
	}

	/* This should extend the array with a flipped conjugated version that has 2 less elements. */
	t_index = regular_strain->len - 2;
	c_index = regular_strain->len;
	for (; t_index > 0; t_index--, c_index++) {
		out_c[c_index] = gsl_complex_conjugate(temp[t_index]);
	}
	/*SS_make_two_sided( regular_strain->len, temp, )*/
}

void CN_save(char* filename, size_t len, double* tmp_ifft) {
	FILE* file;
	size_t i;

	file = fopen(filename, "w");
	for (i = 0; i < len; i++) {
		fprintf(file, "%e\n", tmp_ifft[i]);
	}
	fclose(file);
}

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
		double *out_val)
{
	double UdotU_input;
	double UdotV_input;
	double VdotV_input;
	size_t i;
	double A_input;
	double B_input;
	double C_input;
	double Delta_input;
	double Delta_factor_input;
	double D_input;
	double P1_input;
	double P2_input;
	double P3_input;
	double P4_input;
	double G1_input;
	double G2_input;
	double O11_input;
	double O12_input;
	double O21_input;
	double O22_input;
	size_t s;
	size_t tid;
	size_t did;
	size_t fid;
	size_t j;
	double max;
	/* double std; */

	/* Compute the antenna patterns for each detector */
	for (i = 0; i < net->num_detectors; i++) {
		antenna_patterns(net->detector[i], sky, polarization_angle,
				workspace->ap_workspace, &workspace->ap[i]);
	}

	/* We need to make vectors with the same number of dimensions as the number of detectors in the network */
	UdotU_input = 0.0;
	UdotV_input = 0.0;
	VdotV_input = 0.0;

	/* dot product */
	for (i = 0; i < net->num_detectors; i++) {
		UdotU_input += workspace->ap[i].u * workspace->ap[i].u;
		UdotV_input += workspace->ap[i].u * workspace->ap[i].v;
		VdotV_input += workspace->ap[i].v * workspace->ap[i].v;
	}

	A_input = UdotU_input;
	B_input = UdotV_input;
	C_input = VdotV_input;

	Delta_input = (A_input*C_input) - (B_input*B_input);
	Delta_factor_input = 1.0 / sqrt(2.0*Delta_input);

	D_input = sqrt( gsl_pow_2(A_input-C_input) + 4.0*gsl_pow_2(B_input) ) ;

	P1_input = (C_input-A_input-D_input);
	P2_input = (C_input-A_input+D_input);
	P3_input = sqrt(C_input+A_input+D_input);
	P4_input = sqrt(C_input+A_input-D_input);

	G1_input =  sqrt( gsl_pow_2(P1_input) +  4.0*gsl_pow_2(B_input)) / (2.0*B_input) ;
	G2_input =  sqrt( gsl_pow_2(P2_input) +  4.0*gsl_pow_2(B_input)) / (2.0*B_input) ;

	O11_input = Delta_factor_input * P3_input / G1_input ;
	O12_input = Delta_factor_input * P3_input * P1_input / (2.0*B_input*G1_input) ;
	O21_input = Delta_factor_input * P4_input / G2_input ;
	O22_input  = Delta_factor_input * P4_input * P2_input / (2.0*B_input*G2_input) ;

	/* Loop over each detector to generate a template and do matched filtering */
	for (i = 0; i < net->num_detectors; i++) {
		detector_t* det;
		double coalesce_phase;
		gsl_complex* whitened_data;
		double U_vec_input;
		double V_vec_input;
		double td;

		det = net->detector[i];

		coalesce_phase = 0.0;

		/* Compute time delay */
		time_delay(det, sky, &td);

		SP_compute(coalesce_phase, td,
						chirp, regular_strain,
						f_low, f_high,
						workspace->sp);

		whitened_data = signals[i]->whitened_data;

		do_work(workspace->sp->spa_0, regular_strain, whitened_data, workspace->temp_array, workspace->helpers[i]->c_plus);

		do_work(workspace->sp->spa_90, regular_strain, whitened_data, workspace->temp_array, workspace->helpers[i]->c_minus);

		U_vec_input = workspace->ap[i].u;
		V_vec_input = workspace->ap[i].v;

		workspace->helpers[i]->w_plus_input = (O11_input*U_vec_input +  O12_input*V_vec_input);
		workspace->helpers[i]->w_minus_input = (O21_input*U_vec_input +  O22_input*V_vec_input);
	}

	/* zero the memory */
	s = 2 * regular_strain->len - 2;
	for (tid = 0; tid < 4; tid++) {
		memset( workspace->terms[tid], 0, s * sizeof(gsl_complex) );
		memset( workspace->fs[tid], 0, s * sizeof(gsl_complex) );
	}

	for (did = 0; did < net->num_detectors; did++) {
		for (fid = 0; fid < s; fid++) {
			gsl_complex t;

			t = gsl_complex_mul_real(workspace->helpers[did]->c_plus[fid], workspace->helpers[did]->w_plus_input);
			workspace->terms[0][fid] = gsl_complex_add( workspace->terms[0][fid], t);

			t = gsl_complex_mul_real(workspace->helpers[did]->c_plus[fid], workspace->helpers[did]->w_minus_input);
			workspace->terms[1][fid] = gsl_complex_add( workspace->terms[1][fid], t);

			t = gsl_complex_mul_real(workspace->helpers[did]->c_minus[fid], workspace->helpers[did]->w_plus_input);
			workspace->terms[2][fid] = gsl_complex_add( workspace->terms[2][fid], t);

			t = gsl_complex_mul_real(workspace->helpers[did]->c_minus[fid], workspace->helpers[did]->w_minus_input);
			workspace->terms[3][fid] = gsl_complex_add( workspace->terms[3][fid], t);
		}
	}

	for (i = 0; i < 4; i++) {
		for (j = 0; j < s; j++) {
			workspace->fs[i][2*j + 0] = GSL_REAL( workspace->terms[i][j] );
			workspace->fs[i][2*j + 1] = GSL_IMAG( workspace->terms[i][j] );
		}
		gsl_fft_complex_inverse( workspace->fs[i], 1, s, workspace->fft_wavetable, workspace->fft_workspace );
	}

	memset(workspace->temp_ifft, 0, s * sizeof(double));
	for (i = 0; i < 4; i++) {
		for (j = 0; j < s; j++) {
			/* Take only the real part. The imaginary part should be zero. */
			double x = workspace->fs[i][2*j + 0];
			workspace->temp_ifft[j] += gsl_pow_2(x*s);
		}
	}

	/* CN_save("tmp_ifft.dat", s, workspace->temp_ifft); */

	max = workspace->temp_ifft[0];
	for (i = 1; i < s; i++) {
		double m = workspace->temp_ifft[i];
		if (m > max) {
			max = m;
		}
	}

	/* divide the standard deviation to convert to SNR */
	/* std = gsl_stats_sd(workspace->temp_ifft, 1, s);

	*out_val = sqrt(max) / std; */
	*out_val = sqrt(max);
}

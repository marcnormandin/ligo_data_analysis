#include <assert.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_cblas.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics_double.h>

#include "detector.h"
#include "detector_antenna_patterns.h"
#include "detector_network.h"
#include "detector_time_delay.h"
#include "inspiral_network_statistic.h"
#include "sampling_system.h"
#include "inspiral_chirp.h"
#include "inspiral_chirp_time.h"
#include "hdf5_file.h"

/* this routine was written for the PSO code. */
void CN_template_chirp_time(double f_low, double chirp_time0, double chirp_time1_5, inspiral_chirp_time_t *ct) {
	assert(ct != NULL);

	/* f_low and these are used to compute the required chirp times. */
	ct->chirp_time0 = chirp_time0;
	ct->chirp_time1_5 = chirp_time1_5;

	double calculated_reduced_mass = Chirp_Calc_CalculatedReducedMass(f_low, chirp_time0, chirp_time1_5);
	double calculated_total_mass = Chirp_Calc_CalculatedTotalMass(f_low, chirp_time0, chirp_time1_5);
	double multi_fac_cal = Chirp_Calc_MultiFacCal(f_low, calculated_total_mass);
	double s_mass_ratio_cal = Chirp_Calc_SMassRatioCal(calculated_reduced_mass, calculated_total_mass);

	ct->chirp_time1 =  Chirp_Calc_Time1(f_low, multi_fac_cal, s_mass_ratio_cal);
	ct->chirp_time2 =  Chirp_Calc_Time2(f_low, multi_fac_cal, s_mass_ratio_cal);

	double calc_tchirp = Chirp_Calc_TChirp(ct->chirp_time0, ct->chirp_time1, ct->chirp_time1_5, ct->chirp_time2);

	/* careful that this doesn't use time of arrival because the network statistic wants it as 0 */
	ct->tc = calc_tchirp;
}

coherent_network_helper_t* CN_helper_alloc(size_t num_time_samples) {
	size_t i;
	coherent_network_helper_t *h;

	h = (coherent_network_helper_t*) malloc( sizeof(coherent_network_helper_t) );
	if (h == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory for coherent_network_helper_t. Exiting.\n");
		exit(-1);
	}

	h->c_plus = (gsl_complex*) malloc( num_time_samples * sizeof(gsl_complex) );
	if (h->c_plus == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory for h->c_plus. Exiting.\n");
		exit(-1);
	}

	h->c_minus = (gsl_complex*) malloc( num_time_samples * sizeof(gsl_complex) );
	if (h->c_minus == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory for h->c_minus. Exiting.\n");
		exit(-1);
	}

	for (i = 0; i < num_time_samples; i++) {
		h->c_plus[i] = gsl_complex_rect(0.0, 0.0);
		h->c_minus[i] = gsl_complex_rect(0.0, 0.0);
	}

	return h;
}

void CN_helper_free( coherent_network_helper_t* helper) {
	assert(helper != NULL);

	assert(helper->c_plus != NULL);
	free(helper->c_plus);
	helper->c_plus = NULL;

	assert(helper->c_minus != NULL);
	free(helper->c_minus);
	helper->c_minus = NULL;

	free(helper);
}

coherent_network_workspace_t* CN_workspace_alloc(size_t num_time_samples, detector_network_t *net, size_t num_half_freq,
		double f_low, double f_high) {
	assert(net != NULL);

	coherent_network_workspace_t * work;
	size_t i, j;

	work = (coherent_network_workspace_t*) malloc(sizeof(coherent_network_workspace_t));
	if (work == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory: CN_workspace_malloc(). Exiting.\n");
		exit(-1);
	}

	work->num_time_samples = num_time_samples;
	work->num_helpers = net->num_detectors;

	work->helpers = (coherent_network_helper_t**) malloc( work->num_helpers * sizeof(coherent_network_helper_t*));
	if (work->helpers == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory: CN_workspace_malloc(). Exiting.\n");
		exit(-1);
	}

	for (i = 0; i < work->num_helpers; i++) {
		work->helpers[i] = CN_helper_alloc( num_time_samples );
	}

	/* Note, the asd is only needed to get the frequency values and the number of frequency bins. This should be the
	 * same for every ASD used for a detector network, so any detector from the network can be used.
	 */
	work->sp_lookup = SP_workspace_alloc(f_low, f_high, net->detector[0]->asd->len, net->detector[0]->asd->f);

	work->sp = SP_alloc( num_half_freq );

	work->temp_array = (gsl_complex*) malloc( num_half_freq * sizeof(gsl_complex) );
	if (work->temp_array == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory: CN_workspace_malloc(). Exiting.\n");
		exit(-1);
	}

	for (i = 0; i < num_half_freq; i++) {
		work->temp_array[i] = gsl_complex_rect(0.0, 0.0);
	}

	work->terms = (gsl_complex**) malloc(4 * sizeof(gsl_complex*) );
	if (work->terms == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory: CN_workspace_malloc(). Exiting.\n");
		exit(-1);
	}
	for (i = 0; i < 4; i++) {
		work->terms[i] = (gsl_complex*) malloc( num_time_samples * sizeof(gsl_complex) );
		for (j = 0; j < num_time_samples; j++) {
			work->terms[i][j] = gsl_complex_rect(0.0, 0.0);
		}
	}

	work->fs = (double**) malloc( 4 * sizeof(double*) );
	if (work->fs == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory: CN_workspace_malloc(). Exiting.\n");
		exit(-1);
	}

	for (i = 0; i < 4; i++) {
		work->fs[i] = (double*) malloc( 2 * num_time_samples * sizeof(double) );
		if (work->fs[i] == NULL) {
			fprintf(stderr, "Error. Unable to allocate memory: CN_workspace_malloc(). Exiting.\n");
			exit(-1);
		}
	}

	work->temp_ifft = (double*) malloc( num_time_samples * sizeof(double) );
	if (work->temp_ifft == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory: CN_workspace_malloc(). Exiting.\n");
		exit(-1);
	}

	work->fft_wavetable = gsl_fft_complex_wavetable_alloc( num_time_samples );
	work->fft_workspace = gsl_fft_complex_workspace_alloc( num_time_samples );

	work->ap_workspace = Detector_Antenna_Patterns_workspace_alloc();

	/* one antenna pattern structure per detector */
	work->ap = (detector_antenna_patterns_t*) malloc (net->num_detectors * sizeof(detector_antenna_patterns_t) );
	if (work->ap == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory: CN_workspace_malloc(). Exiting.\n");
		exit(-1);
	}

	/* Each detector has a normalization factor for the stationary phase inner product. */
	work->normalization_factors = (double*) malloc( net->num_detectors * sizeof(double) );
	if (work->normalization_factors == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory: CN_workspace_malloc(). Exiting.\n");
		exit(-1);
	}
	for (i = 0; i < net->num_detectors; i++) {
		work->normalization_factors[i] = SP_normalization_factor(net->detector[i]->asd, work->sp_lookup);
	}

	return work;
}

void CN_workspace_free( coherent_network_workspace_t *workspace ) {
	assert(workspace != NULL);

	size_t i;

	for (i = 0; i < workspace->num_helpers; i++) {
		CN_helper_free( workspace->helpers[i] );
		workspace->helpers[i] = NULL;
	}
	free(workspace->helpers);
	workspace->helpers = NULL;

	SP_workspace_free(workspace->sp_lookup);
	workspace->sp_lookup = NULL;

	SP_free(workspace->sp);
	workspace->sp = NULL;

	for (i = 0; i < 4; i++) {
		free(workspace->terms[i]);
		workspace->terms[i] = NULL;

		free(workspace->fs[i]);
		workspace->fs[i] = NULL;
	}

	free(workspace->terms);
	workspace->terms = NULL;

	free(workspace->fs);
	workspace->fs = NULL;

	free(workspace->temp_ifft);
	workspace->temp_ifft = NULL;

	free(workspace->temp_array);
	workspace->temp_array = NULL;

	gsl_fft_complex_workspace_free( workspace->fft_workspace );
	workspace->fft_workspace = NULL;

	gsl_fft_complex_wavetable_free( workspace->fft_wavetable );
	workspace->fft_wavetable = NULL;

	Detector_Antenna_Patterns_workspace_free(workspace->ap_workspace);
	workspace->ap_workspace = NULL;

	free(workspace->ap);
	workspace->ap = NULL;

	free(workspace->normalization_factors);
	workspace->normalization_factors = NULL;

	free( workspace );
}

void CN_do_work(size_t num_time_samples, size_t f_low_index, size_t f_high_index, gsl_complex *spa, asd_t *asd, gsl_complex *half_fft_data, gsl_complex *temp, gsl_complex *out_c) {
	assert(spa != NULL);
	assert(asd != NULL);
	assert(half_fft_data != NULL);
	assert(temp != NULL);
	assert(out_c != NULL);

	size_t k;
	size_t t_index;
	size_t c_index;

	// faster version for (k = f_low_index; k <= f_high_index; k++) {
	for (k = 0; k < asd->len; k++) {
		temp[k] = gsl_complex_conjugate(spa[k]);
		temp[k] = gsl_complex_div_real(temp[k], asd->asd[k]);
		temp[k] = gsl_complex_mul( temp[k], half_fft_data[k] );
	}

	/* This should extend the array with a flipped conjugated version. */
	SS_make_two_sided( asd->len, temp, num_time_samples, out_c);
}

void CN_save(char* filename, size_t len, double* tmp_ifft) {
	assert(filename != NULL);
	assert(tmp_ifft != NULL);

	FILE* file;
	size_t i;

	file = fopen(filename, "w");
	if (file == NULL) {
		fprintf(stderr, "Error. CN_save: Unable to open the file (%s) for writing the network statistic series. Exiting.\n", filename);
		exit(-1);
	}

	for (i = 0; i < len; i++) {
		fprintf(file, "%e\n", tmp_ifft[i]);
	}
	fclose(file);
}

/* DANGER. This assumes that the coalece phase is 0 */
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
		char *out_network_css_filename)
{
	assert(net);
	assert(chirp);
	assert(sky);
	assert(network_strain);
	assert(workspace);
	assert(out_network_css_value);
	assert(out_network_css_index);

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
	size_t tid;
	size_t did;
	size_t fid;
	size_t j;
	double max_value;
	size_t max_index;

	/* WARNING: This assumes that all of the signals have the same lengths. */
	size_t num_time_samples = network_strain->num_time_samples;

	/* Compute the antenna patterns for each detector */
	for (i = 0; i < net->num_detectors; i++) {
		double polarization_angle = 0.0; // Shihan said only u and v are needed for templates.
		Detector_Antenna_Patterns_compute(net->detector[i], sky, polarization_angle,
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
	O22_input  = Delta_factor_input * P4_input * P2_input / (2.0*B_input*G2_input);

	/* Loop over each detector to generate a template and do matched filtering */
	for (i = 0; i < net->num_detectors; i++) {
		detector_t* det;
		double inspiral_coalesce_phase;
		gsl_complex* whitened_data;
		double U_vec_input;
		double V_vec_input;
		double detector_time_delay;

		det = net->detector[i];

		/* For reconstruction use the phase as 0 */
		inspiral_coalesce_phase = 0.0;

		/* Compute time delay */
		Detector_time_delay(det, sky, &detector_time_delay);

		/*printf("g = %0.21e\n", workspace->normalization_factors[i]);*/

		SP_compute(		detector_time_delay, workspace->normalization_factors[i],
						inspiral_coalesce_phase, chirp,
						workspace->sp_lookup,
						workspace->sp);

		/*
		printf("Detector SPA_0: %s\n", det->name);
		for (j = 0; j < workspace->sp->len; j++) {
			printf("%0.21e \t %0.21e\n", GSL_REAL(workspace->sp->spa_0[j]), GSL_IMAG(workspace->sp->spa_0[j]));
		}

		printf("Detector SPA_90: %s\n", det->name);
		for (j = 0; j < workspace->sp->len; j++) {
			printf("%0.21e \t %0.21e\n", GSL_REAL(workspace->sp->spa_90[j]), GSL_IMAG(workspace->sp->spa_90[j]));
		}*/

		whitened_data = network_strain->strains[i]->half_fft;

		/* compute c_plus */
		CN_do_work(num_time_samples, workspace->sp_lookup->f_low_index, workspace->sp_lookup->f_high_index, workspace->sp->spa_0, det->asd, whitened_data, workspace->temp_array, workspace->helpers[i]->c_plus);

		/* compute c_minus */
		CN_do_work(num_time_samples, workspace->sp_lookup->f_low_index, workspace->sp_lookup->f_high_index, workspace->sp->spa_90, det->asd, whitened_data, workspace->temp_array, workspace->helpers[i]->c_minus);

		U_vec_input = workspace->ap[i].u;
		V_vec_input = workspace->ap[i].v;

		workspace->helpers[i]->w_plus_input = (O11_input*U_vec_input +  O12_input*V_vec_input);
		workspace->helpers[i]->w_minus_input = (O21_input*U_vec_input +  O22_input*V_vec_input);
	}

	/* zero the memory */
	for (tid = 0; tid < 4; tid++) {
		memset( workspace->terms[tid], 0, num_time_samples * sizeof(gsl_complex) );
		memset( workspace->fs[tid], 0, num_time_samples * sizeof(gsl_complex) );
	}

	for (did = 0; did < net->num_detectors; did++) {
		for (fid = 0; fid < num_time_samples; fid++) {
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
		for (j = 0; j < num_time_samples; j++) {
			workspace->fs[i][2*j + 0] = GSL_REAL( workspace->terms[i][j] );
			workspace->fs[i][2*j + 1] = GSL_IMAG( workspace->terms[i][j] );
		}
		gsl_fft_complex_inverse( workspace->fs[i], 1, num_time_samples, workspace->fft_wavetable, workspace->fft_workspace );
	}

	memset(workspace->temp_ifft, 0, num_time_samples * sizeof(double));
	for (i = 0; i < 4; i++) {
		for (j = 0; j < num_time_samples; j++) {
			/* Take only the real part. The imaginary part should be zero. */
			double x = workspace->fs[i][2*j + 0];
			workspace->temp_ifft[j] += gsl_pow_2(x*num_time_samples);
		}
	}

	/*CN_save("tmp_ifft.dat", s, workspace->temp_ifft);*/

	max_index = 0;
	max_value = workspace->temp_ifft[0];

	/* check statistical behavior of this time series */
	for (i = 1; i < num_time_samples; i++) {
		double m = workspace->temp_ifft[i];
		if (m > max_value) {
			max_value = m;
			max_index = i;
		}
	}

	/* check, sqrt sbould behave according to chi */
	/* check, use this with just noise and see if the mean is 4, std should be sqrt(8). Chi-sqre if not sqrt(max). Check 'max' dist.*/
	double old_snr_definition = sqrt(max_value) / sqrt(2.0);

	// temp hack to do this quickly
	//workspace->temp_ifft[max_index] = 0.0;
	//double std = gsl_stats_sd (workspace->temp_ifft, 1, workspace->num_time_samples);
	//workspace->temp_ifft[max_index] = max_value;

	//double new_snr_definition = max_value / std;
	*out_network_css_value = old_snr_definition;
	*out_network_css_index = max_index;

	/*
	if (hdf5_filename != NULL && hdf5_dataset_name != NULL) {
		hdf5_save_array( hdf5_filename, "/", hdf5_dataset_name, workspace->num_time_samples, workspace->temp_ifft);

		for (i = 0; i < workspace->num_time_samples; i++) {
			printf("%f ", workspace->temp_ifft[i]);
		}
	}
	*/
	if (out_network_css_filename != NULL) {
		CN_save( out_network_css_filename, num_time_samples, workspace->temp_ifft);
	}
}

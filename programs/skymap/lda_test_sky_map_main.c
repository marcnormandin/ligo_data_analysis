/*
 * lda_sky_map_main.c
 *
 *  Created on: Mar 16, 2017
 *      Author: marcnormandin
 */

#include "detector_antenna_patterns.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "detector_antenna_patterns.h"
#include "inspiral_chirp.h"
#include "detector.h"
#include "detector_network.h"
#include "inspiral_source.h"
#include "inspiral_stationary_phase.h"
#include "strain.h"
#include "strain_interpolate.h"
#include "signal.h"
#include "inspiral_network_statistic.h"

#include "inspiral_pso_fitness.h"
#include "inspiral_signal.h"

#include "common.h"
#include "random.h"
#include "options.h"

#include "settings.h"

void snr_sky_map(pso_fitness_function_parameters_t *splParams, int num_ra, int num_dec, const char* output_file) {
	size_t ra_i, dec_i;
	const size_t N_ra = num_ra;
	const size_t N_dec = num_dec;

	const double delta_ra = (2.0 * M_PI) / N_ra;
	const double delta_dec = (1.0*M_PI) / N_dec;

	FILE *fid = fopen(output_file, "w");

	for (ra_i = 0; ra_i < N_ra; ra_i++) {
		for (dec_i = 0; dec_i < N_dec; dec_i++) {
			inspiral_chirp_factors_t chirp;
			double recovered_snr;
			double ra = -M_PI + ra_i * delta_ra;
			double dec = -0.5*M_PI + dec_i * delta_dec;

			/* apply the pso particle location */
			splParams->source->sky.ra = ra;
			splParams->source->sky.dec = dec;

			CF_compute(splParams->f_low, splParams->source, &chirp);
			/* The network statistic requires the time of arrival to be zero
			   in order for the matched filtering to work correctly. */
			chirp.ct.tc = chirp.t_chirp;

			coherent_network_statistic(
					splParams->network,
					splParams->strain,
					splParams->f_low,
					splParams->f_high,
					&chirp.ct,
					&splParams->source->sky,
					splParams->source->polarization_angle,
					splParams->signals,
					splParams->workspace,
					&recovered_snr);

			fprintf(fid, "%f\t %f\t %f\n", ra, dec, recovered_snr);
		}
	}

	fclose(fid);

}

int main(int argc, char* argv[]) {
	size_t i;

	source_t source;
	Source_load_testsource(&source);
	gslseed_t seed;
	int last_index;

	seed = 0;

	/* Allow options to override the defaults settings */
	if (process_command_options(argc, argv, &source, &seed, &last_index) != 0) {
		exit(0);
	}

	/* somehow these need to be set */
	if ((argc - last_index) < 4) {
		printf("Error: Usage -> [settings file] [num RA divisions] [num DEC divisions] [output data file]\n");
		exit(-1);
	}

	char* arg_settings_file = argv[last_index++];
	int arg_num_ra = atoi(argv[last_index++]);
	int arg_num_dec = atoi(argv[last_index++]);
	char* arg_output_file = argv[last_index++];



	/* Settings */
	settings_file_t *settings_file = settings_file_open(arg_settings_file);
	if (settings_file == NULL) {
		printf("Error opening the settings file (%s). Aborting.\n", arg_settings_file);
		abort();
	}

	printf("Using the following settings:\n");
	settings_file_print(settings_file);

	const double f_low = atof(settings_file_get_value(settings_file, "f_low"));
	const double f_high = atof(settings_file_get_value(settings_file, "f_high"));
	const double sampling_frequency = atof(settings_file_get_value(settings_file, "sampling_frequency"));
	const double num_time_samples = atoi(settings_file_get_value(settings_file, "num_time_samples"));

	settings_file_close(settings_file);

	printf("num_ra: %d\n", arg_num_ra);
	printf("num_dec: %d\n", arg_num_dec);
	printf("\n\n");

	/*
	const double f_low = 40.0;
	const double f_high = 700.0;
	const double sampling_frequency = 2048.0;
	const size_t num_time_samples = 131072;
	*/


	detector_network_t net;
	Init_Detector_Network(&net);

	strain_t *strain = Strain_simulated(f_low, f_high, sampling_frequency, num_time_samples);

	/* Random number generator */
	gsl_rng *rng = random_alloc(seed);

	/* Simulate data for all the detectors composing the network */
	strain_half_fft_t **signals = simulate_inspiral(rng, f_low, f_high, &net, strain, &source);

	coherent_network_workspace_t *workspace = CN_workspace_alloc( net.num_detectors, Strain_one_sided_length(strain) );

	/* Setup the parameter structure for the pso fitness function */
	pso_fitness_function_parameters_t params;
	params.f_low = f_low;
	params.f_high = f_high;
	params.network = &net;
	params.signals = signals;
	params.source = &source;
	params.strain = strain;
	params.workspace = workspace;

	printf("The true source is at: RA = %f, DEC = %f\n", params.source->sky.ra, params.source->sky.dec);

	snr_sky_map(&params, arg_num_ra, arg_num_dec, arg_output_file);

	printf("The data has been saved to snr_sky_map.dat.\n");
	printf("Plots can be made using lda_test_sky_map.ipynb.\n");

	CN_workspace_free( workspace );

	/* Free the data */
	for (i = 0; i < net.num_detectors; i++) {
		strain_half_fft_free(signals[i]);
	}
	free(signals);

	Detector_Network_free(&net);
	Strain_free(strain);
	random_free(rng);

	return 0;
}




/*
 * matlab_data_serial.c
 *
 *  Created on: Apr 17, 2017
 *      Author: marcnormandin
 */

#ifdef HAVE_CONFIG_H
	#include "config.h"
#endif


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
#include "strain.h"
#include "inspiral_stationary_phase.h"
#include "inspiral_network_statistic.h"

#include "inspiral_pso_fitness.h"
#include "random.h"
#include "settings_file.h"
#include "detector_mapping.h"
#include "hdf5_file.h"
#include "sampling_system.h"


void load_shihan_inspiral_data( const char* hdf_filename, strain_half_fft_t *strain){
	size_t j;

	size_t num_time_samples = hdf5_get_num_time_samples( hdf_filename );
	size_t half_size = SS_half_size( num_time_samples );

	double *real_array = (double*) malloc( half_size * sizeof(double) );
	double *imag_array = (double*) malloc( half_size * sizeof(double) );

	hdf5_load_array( hdf_filename, "/shihan/whitened_data_real", real_array);
	hdf5_load_array( hdf_filename, "/shihan/whitened_data_imag", imag_array);

	for (j = 0; j < half_size; j++) {
		strain->half_fft[j] = gsl_complex_rect(real_array[j], imag_array[j]);
	}

	free(real_array);
	free(imag_array);
}

void pso_result_save(FILE *fid, pso_result_t *result) {
	fprintf(fid, "%20.17g %20.17g %20.17g %20.17g %20.17g %20zu %20zu %20.17g",
			result->ra, result->dec, result->chirp_t0, result->chirp_t1_5, result->snr,
			result->total_iterations, result->total_func_evals, result->computation_time_secs);
}

void pso_result_print(pso_result_t *result) {
	printf("%20.17g %20.17g %20.17g %20.17g %20.17g %20zu %20zu %20.17g",
			result->ra, result->dec, result->chirp_t0, result->chirp_t1_5, result->snr,
			result->total_iterations, result->total_func_evals, result->computation_time_secs);
}

typedef struct callback_function_params_s {
	FILE *fid;
	pso_result_t result;
	pso_ranges_t *ranges;

} callback_function_params_t;

callback_function_params_t* callback_function_params_alloc( const char* pso_settings_filename, const char* callback_save_filename ) {
	callback_function_params_t* p = (callback_function_params_t*) malloc( sizeof(callback_function_params_t) );
	p->ranges = pso_ranges_alloc( pso_settings_filename );
	pso_ranges_init( pso_settings_filename, p->ranges );
	p->fid = fopen( callback_save_filename, "a");
	return p;
}

void callback_function_params_free( callback_function_params_t* p ) {
	fclose(p->fid);
	pso_ranges_free(p->ranges);
	free(p);
}

void callback(void* params, returnData_t *pso_domain_result) {
	callback_function_params_t* p = (callback_function_params_t*) params;
	return_data_to_pso_results( p->ranges, pso_domain_result, &p->result );
	pso_result_save( p->fid, &p->result );
	fprintf(p->fid, "\n");
}

int main(int argc, char* argv[]) {
	size_t i;

	/* somehow these need to be set */
	if (argc != 7) {
		printf("argc = %d\n", argc);
		printf("Error: Usage -> [settings file] [detector mapping file] [rng seed] [input pso settings file] [output final pso results file] [output intermediate results file]\n");
		exit(-1);
	}

	char* arg_settings_file = argv[1];
	char* arg_detector_mapping_file = argv[2];
	const gslseed_t seed = atoi(argv[3]);
	char* arg_pso_settings_file = argv[4];
	char* arg_pso_results_file = argv[5];
	char* arg_pso_record_file = argv[6];

	/* Load the general Settings */
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
	const int arg_pso_record_interval = atoi(settings_file_get_value(settings_file, "pso_callback_interval"));

	settings_file_close(settings_file);

	detector_network_mapping_t *dmap = Detector_Network_Mapping_load( arg_detector_mapping_file );
	size_t num_time_samples = hdf5_get_num_time_samples( dmap->data_filenames[0] );
	detector_network_t *net = Detector_Network_load(
			arg_detector_mapping_file, num_time_samples, sampling_frequency, f_low, f_high );

	network_strain_half_fft_t *network_strain = network_strain_half_fft_alloc(dmap->num_detectors, num_time_samples );
	for (i = 0; i < net->num_detectors; i++) {
		load_shihan_inspiral_data( dmap->data_filenames[i], network_strain->strains[i] );
	}

	pso_fitness_function_parameters_t *fitness_function_params =
				pso_fitness_function_parameters_alloc(f_low, f_high, net, network_strain);

	callback_function_params_t* cbp = callback_function_params_alloc(arg_pso_settings_file, arg_pso_record_file);

	current_result_callback_params_t callback_params;
	callback_params.interval = arg_pso_record_interval;
	callback_params.callback = callback;
	callback_params.callback_params = cbp;

	pso_result_t pso_result;
	pso_estimate_parameters(arg_pso_settings_file, fitness_function_params, &callback_params, seed, &pso_result);

	callback_function_params_free( cbp );

	FILE *fid = fopen(arg_pso_results_file, "a");
	pso_result_save(fid, &pso_result);
	fprintf(fid, "\n");
	fclose(fid);

	pso_result_print(&pso_result);


	pso_fitness_function_parameters_free(fitness_function_params);

	/* Free the data */
	network_strain_half_fft_free(network_strain);

	Detector_Network_free(net);

	return 0;
}

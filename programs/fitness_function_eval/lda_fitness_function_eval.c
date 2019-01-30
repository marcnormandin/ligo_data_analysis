#ifdef HAVE_CONFIG_H
	#include "config.h"
#endif

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <memory.h>
#include <assert.h>

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

#include "settings_file.h"
#include "detector_mapping.h"
#include "hdf5_file.h"
#include "sampling_system.h"

#define MAX_FILENAME_LEN 255

typedef struct compute_workspace_s {
	double f_low;
	double f_high;
	detector_network_t *network;
	network_strain_half_fft_t *network_strain;
	coherent_network_workspace_t *workspace;
} compute_workspace_t;

typedef struct likelihood_result_s {
	// The values to be computed and saved
	inspiral_chirp_time_t chirp_time;
	double css_value;
	int css_index;
	double tc_seconds; // converted from the css_index

} likelihood_result_t;

// Template parameters to evaluate the likelihood at
typedef struct template_parameters_s {
	double ra;
	double dec;
	double chirp_time_0;
	double chirp_time_1_5;
} template_parameters_t;

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

compute_workspace_t* compute_workspace_alloc(
		double f_low, double f_high, detector_network_t* network, network_strain_half_fft_t *network_strain)
{
	assert(network != NULL);
	assert(network_strain != NULL);

	compute_workspace_t *params = (compute_workspace_t*) malloc( sizeof(compute_workspace_t) );
	if (params == 0) {
		fprintf(stderr, "Error. Unable to allocate memory for compute_workspace_t structure. Exiting.");
		exit(-1);
	}

	params->workspace = CN_workspace_alloc(
				network_strain->num_time_samples, network, network->detector[0]->asd->len,
				f_low, f_high);

	/* Setup the parameter structure for the pso fitness function */
	params->f_low = f_low;
	params->f_high = f_high;
	params->network = network;
	params->network_strain = network_strain;

	//fprintf(stderr, "Number of threads: %lu\n", parallel_get_max_threads());

	return params;
}

void compute_workspace_free(compute_workspace_t *params) {
	assert(params != NULL);

	CN_workspace_free(params->workspace);
	params->workspace = NULL;

	free(params);
}

void compute_statistic(compute_workspace_t* workspace, template_parameters_t* params, likelihood_result_t* result, char* css_time_series) {
	assert(workspace != NULL);
	assert(params != NULL);

	CN_template_chirp_time(workspace->f_low, params->chirp_time_0, params->chirp_time_1_5, &result->chirp_time);

	sky_t sky;
	sky.ra = params->ra;
	sky.dec = params->dec;

	coherent_network_statistic(
			workspace->network,
			workspace->f_low,
			workspace->f_high,
			&result->chirp_time,
			&sky,
			workspace->network_strain,
			workspace->workspace,
			&result->css_value,
			&result->css_index,
			css_time_series);
}

int main(int argc, char* argv[]) {
	size_t i;

	/* somehow these need to be set */
	if (argc != 7) {
		printf("argc = %d\n", argc);
		printf("Error: Usage -> [settings file] [map filename] [ra] [dec] [tau 0] [tau 1.5]\n");
		exit(-1);
	}

	// Map the values from the command line
	char* arg_settings_file = argv[1];
	char* arg_dmap_filename = argv[2];
	// ra argv[3]
	// dec argv[4]
	// tau 0 argv[5]
	// tau 1.5 argv[6]
	//char* arg_results_file  = argv[7];


	/* Load the general Settings */
	settings_file_t *settings_file = settings_file_open(arg_settings_file);
	if (settings_file == NULL) {
		printf("Error opening the settings file (%s). Aborting.\n", arg_settings_file);
		abort();
	}

	settings_file_print(settings_file);

	const double f_low = atof(settings_file_get_value(settings_file, "f_low"));
	const double f_high = atof(settings_file_get_value(settings_file, "f_high"));
	const double sampling_frequency = atof(settings_file_get_value(settings_file, "sampling_frequency"));

	settings_file_close(settings_file);

	detector_network_mapping_t *dmap = Detector_Network_Mapping_load( arg_dmap_filename );
	size_t num_time_samples = hdf5_get_num_time_samples( dmap->data_filenames[0] );
	detector_network_t *net = Detector_Network_load(
			arg_dmap_filename, num_time_samples, sampling_frequency, f_low, f_high );

	network_strain_half_fft_t *network_strain = network_strain_half_fft_alloc(dmap->num_detectors, num_time_samples );
	for (i = 0; i < net->num_detectors; i++) {
		load_shihan_inspiral_data( dmap->data_filenames[i], network_strain->strains[i] );
	}


	// Make the template paramater structure holding the template parametes to evaluate
	template_parameters_t params;
	params.ra = atof(argv[3]);
	params.dec = atof(argv[4]);
	params.chirp_time_0 = atof(argv[5]);
	params.chirp_time_1_5 = atof(argv[6]);

	// The result structure that the fitness values will be saved to
	likelihood_result_t result;

	// Change to a filename if you want the time series to be saved to a file
	char* css_time_series = NULL;

	// Perform the fitness function computation
	compute_workspace_t *workspace = compute_workspace_alloc(f_low, f_high, net, network_strain);
	compute_statistic(workspace, &params, &result, css_time_series);	

	// convert the index of the tc into the tc value in seconds
	result.tc_seconds = ((double) result.css_index) / sampling_frequency;

	// print the results to the screen
	printf("css_value at tc_index: %20.17g\n", result.css_value);
	printf("tc_index: %20zu\n", result.css_index);
	printf("tc_seconds: %20.17g\n", result.tc_seconds);
	printf("calculated chirp time (using tau 0 and tau 1.5): %20.17g\n", result.chirp_time);


	// Clean up
	compute_workspace_free(workspace);
	network_strain_half_fft_free(network_strain);
	Detector_Network_free(net);

	return 0;
}



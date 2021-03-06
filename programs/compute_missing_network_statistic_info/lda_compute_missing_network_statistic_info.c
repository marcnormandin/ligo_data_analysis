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

typedef struct est_params_s {
	// The values that were calculated
	double ra, dec, chirp_time_0, chirp_time_1_5;
	double css;
	size_t total_iterations, total_func_evals;
	double computation_time_secs;

	// Init this
	char out_network_css_filename[MAX_FILENAME_LEN];

	// The values to be computed and saved
	inspiral_chirp_time_t chirp_time;
	double new_css_value;
	int new_css_index;
	double new_tc_value; // converted from the new_css_index

} est_params_t;


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

void compute_missing(compute_workspace_t* workspace, est_params_t* params) {
	assert(workspace != NULL);
	assert(params != NULL);

	CN_template_chirp_time(workspace->f_low, params->chirp_time_0, params->chirp_time_1_5, &params->chirp_time);

	sky_t sky;
	sky.ra = params->ra;
	sky.dec = params->dec;

	coherent_network_statistic(
			workspace->network,
			workspace->f_low,
			workspace->f_high,
			&params->chirp_time,
			&sky,
			workspace->network_strain,
			workspace->workspace,
			&params->new_css_value,
			&params->new_css_index,
			params->out_network_css_filename);
}

void parse_pso_filename(char *pso_run_filename, char *dataset_and_rn, int *pso_rn) {
	int i = 0; // index into string we are parsing
	int j = 0; // flag for which token we are writing to
	int k = 0; // index into the token we are writing to

	char pso_rn_str[MAX_FILENAME_LEN];
	memset(pso_rn_str, '\0', sizeof(char)*MAX_FILENAME_LEN);

	for (i = 0; i < strlen(pso_run_filename); i++) {
		char c = pso_run_filename[i];
		//printf("%c ", c);

		if (c == '.') {
			// finish the string we were writing
			if (j == 0) {
				dataset_and_rn[k] = '\0';
			} else if (j == 1) {
				pso_rn_str[k] = '\0';
			}

			j += 1;
			k = 0;
			//printf("\n");
		}
		else {
			if (j == 0) {
				dataset_and_rn[k] = c;
				k++;
			} else if (j == 1) {
				pso_rn_str[k] = c;
				k++;
			} else if (j==2) {
				// We don't care about the extension
				break;
			}
		}
	}

	*pso_rn = atoi( pso_rn_str );

	//printf("done!\n");
}

void get_dmap_filename(char *pso_run_filename, char *dmap_filename) {
	char dataset_and_rn[MAX_FILENAME_LEN];
	memset(dataset_and_rn, 0, sizeof(char)*MAX_FILENAME_LEN);

	int pso_rn;

	parse_pso_filename(pso_run_filename, dataset_and_rn, &pso_rn);

	sprintf(dmap_filename, "%s.map\0", dataset_and_rn);
	//printf("dmap: %s\n", dmap_filename);
}

void get_network_statistic_series_filename(char *pso_run_filename, char *series_filename) {
	char dataset_and_rn[MAX_FILENAME_LEN];
	memset(dataset_and_rn, 0, sizeof(char)*MAX_FILENAME_LEN);

	int pso_rn;

	parse_pso_filename(pso_run_filename, dataset_and_rn, &pso_rn);

	sprintf(series_filename, "%s.%d.series", dataset_and_rn, pso_rn);
	//printf("network statistic inside get_network_statistic: %s\n", series_filename);
}

void get_rerun_filename(char *pso_run_filename, char *rerun_filename) {
	char dataset_and_rn[MAX_FILENAME_LEN];
	memset(dataset_and_rn, 0, sizeof(char)*MAX_FILENAME_LEN);

	int pso_rn;

	parse_pso_filename(pso_run_filename, dataset_and_rn, &pso_rn);
	sprintf(rerun_filename, "%s.%d.rerun\0", dataset_and_rn, pso_rn);
	//printf("%s\n", rerun_filename);
}

void init_est_params(char *pso_run_file, est_params_t* est_params) {

	memset(est_params->out_network_css_filename, '\0', MAX_FILENAME_LEN);

	FILE *fid = fopen(pso_run_file, "r");
	if (fid == NULL) {
		fprintf(stderr, "Error. Unable to open the run file (%s) for reading. Exiting.\n", pso_run_file);
		exit(-1);
	}

	//memset(est_params, 0, sizeof(est_params_t));

	int ret = fscanf(fid, "%lf %lf %lf %lf %lf %zu %zu %lf",
			&est_params->ra, &est_params->dec, &est_params->chirp_time_0, &est_params->chirp_time_1_5, &est_params->css,
			&est_params->total_iterations, &est_params->total_func_evals, &est_params->computation_time_secs);

	fclose(fid);

	if (ret != 8) {
            fprintf(stderr, "Error. Only (%d) values read from (%s), but expected 8. Exiting.\n", ret, pso_run_file);
            exit(-1);
    }

	get_network_statistic_series_filename(pso_run_file, est_params->out_network_css_filename);
	//printf("network statistic inside init_est_params: %s\n", est_params->out_network_css_filename);
}



void est_params_save_to_file(est_params_t *est_params, char *rerun_filename) {
	FILE *fid = fopen(rerun_filename, "w");
	fprintf(fid, "%20.17g %20.17g %20.17g %20.17g %20.17g %20.17g %20.17g %20.17g %20.17g %20zu %20.17g %20zu %20zu %20.17g %s",
			est_params->ra, est_params->dec, 
			est_params->chirp_time_0, est_params->chirp_time.chirp_time1, est_params->chirp_time_1_5, est_params->chirp_time.chirp_time2, 
			est_params->chirp_time.tc, est_params->css, est_params->new_css_value, est_params->new_css_index, est_params->new_tc_value,
			est_params->total_iterations, est_params->total_func_evals, est_params->computation_time_secs,
			est_params->out_network_css_filename
	);
	//fprintf(fid, "%s", est_params->out_network_css_filename);
	fclose(fid);
}


void est_params_print(est_params_t *est_params) {
	printf("%20.17g %20.17g %20.17g %20.17g %20.17g %20.17g %20.17g %20.17g %20.17g %d %20.17g %20zu %20zu %20.17g %s\n",
			est_params->ra, est_params->dec, 
			est_params->chirp_time_0, est_params->chirp_time.chirp_time1, est_params->chirp_time_1_5, est_params->chirp_time.chirp_time2, 
			est_params->chirp_time.tc, est_params->css, est_params->new_css_value, est_params->new_css_index, est_params->new_tc_value,
			est_params->total_iterations, est_params->total_func_evals, est_params->computation_time_secs,
			est_params->out_network_css_filename
	);
	//printf("%s\n", est_params->out_network_css_filename);
	//printf("filename length = %d\n", strlen(est_params->out_network_css_filename));
}

int main(int argc, char* argv[]) {
	size_t i;

	/* somehow these need to be set */
	if (argc != 3) {
		printf("argc = %d\n", argc);
		printf("Error: Usage -> [settings file] [input final pso est_params file]\n");
		exit(-1);
	}

	char* arg_settings_file = argv[1];
	//char* arg_detector_mapping_file = argv[2];
	char* arg_pso_run_file = argv[2];

	char dmap_filename[MAX_FILENAME_LEN];
	get_dmap_filename(arg_pso_run_file, dmap_filename);

	/* Load the general Settings */
	settings_file_t *settings_file = settings_file_open(arg_settings_file);
	if (settings_file == NULL) {
		printf("Error opening the settings file (%s). Aborting.\n", arg_settings_file);
		abort();
	}

	//printf("Using the following settings:\n");

	settings_file_print(settings_file);

	const double f_low = atof(settings_file_get_value(settings_file, "f_low"));
	const double f_high = atof(settings_file_get_value(settings_file, "f_high"));
	const double sampling_frequency = atof(settings_file_get_value(settings_file, "sampling_frequency"));

	settings_file_close(settings_file);

	detector_network_mapping_t *dmap = Detector_Network_Mapping_load( dmap_filename );
	size_t num_time_samples = hdf5_get_num_time_samples( dmap->data_filenames[0] );
	detector_network_t *net = Detector_Network_load(
			dmap_filename, num_time_samples, sampling_frequency, f_low, f_high );

	network_strain_half_fft_t *network_strain = network_strain_half_fft_alloc(dmap->num_detectors, num_time_samples );
	for (i = 0; i < net->num_detectors; i++) {
		load_shihan_inspiral_data( dmap->data_filenames[i], network_strain->strains[i] );
	}

	// COMPUTE
	char rerun_filename[MAX_FILENAME_LEN];
	get_rerun_filename(arg_pso_run_file, rerun_filename);

	est_params_t est_params;
	init_est_params(arg_pso_run_file, &est_params);
	//printf("TEST: %s\n", est_params.out_network_css_filename);
	//est_params.out_network_css_filename[0] = 'M';
	//est_params.out_network_css_filename[1] = '\0';
	//est_params_print(&est_params);

	compute_workspace_t *workspace = compute_workspace_alloc(f_low, f_high, net, network_strain);
	compute_missing(workspace, &est_params);	

	// convert the tc_index into the tc value
	est_params.new_tc_value = ((double) est_params.new_css_index) / sampling_frequency;

	est_params_print(&est_params);


	// SAVE
	est_params_save_to_file(&est_params, rerun_filename);

	compute_workspace_free(workspace);
	network_strain_half_fft_free(network_strain);
	Detector_Network_free(net);

	return 0;
}



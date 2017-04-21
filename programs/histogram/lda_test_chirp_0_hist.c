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
#include "inspiral_pso_fitness.h"

void load_shihan_inspiral_data( const char* hdf_filename, strain_half_fft_t *strain){
	size_t i, j;

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

int main(int argc, char* argv[]) {
	size_t i;

	/* somehow these need to be set */
	if (argc != 4) {
		printf("argc = %d\n", argc);
		printf("Error: Must supply [settings file] [detector mapping file] [results file]!\n");
		exit(-1);
	}

	char* arg_settings_file = argv[1];
	char* arg_detector_mapping_file = argv[2];
	char* arg_results_file = argv[3];

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
	settings_file_close(settings_file);

	/* Load the data */
	detector_network_t *net = Detector_Network_load( arg_detector_mapping_file, f_low, f_high );
	detector_network_mapping_t *dmap = Detector_Network_Mapping_load( arg_detector_mapping_file );
	size_t num_time_samples = hdf5_get_num_time_samples( dmap->data_filenames[0] );
	network_strain_half_fft_t *network_strain = network_strain_half_fft_alloc(dmap->num_detectors, num_time_samples );
	for (i = 0; i < net->num_detectors; i++) {
		load_shihan_inspiral_data( dmap->data_filenames[i], network_strain->strains[i] );
	}

	size_t num_half_freq = network_strain->strains[0]->half_fft_len;
	coherent_network_workspace_t *work = CN_workspace_alloc(num_time_samples, net, num_half_freq, f_low, f_high);

	/* Load the true values */
	sky_t sky;
	sky.ra = -2.14;
	sky.dec = 0.72;
	double chirp_time_1_5 = 0.7284;

	/* Run analysis */
	size_t num_points = 150;
	double min_time = 5.0;
	double max_time = 15.0;
	double delta = (max_time - min_time) / (1.0*(num_points-1.0));

	double *snr = (double*) malloc( num_points * sizeof(double) );

	for (i = 0; i < num_points; i++) {
		double chirp_time_0 = min_time + i*delta;
		inspiral_chirp_time_t ct;
		CN_template_chirp_time(f_low, chirp_time_0, chirp_time_1_5, &ct);

		coherent_network_statistic(
				net,
				f_low,
				f_high,
				&ct,
				&sky,
				network_strain,
				work,
				&snr[i]);
	}

	FILE *fid;
	fid = fopen(arg_results_file, "w+");
	for (i = 0; i < num_points; i++) {
		fprintf(fid, "%0.21e %0.21e",  min_time + i*delta, snr[i]);
		if (i < num_points-1) {
			fprintf(fid, "\n");
		}
	}
	fclose(fid);

	/* Clean up */
	free(snr);
	CN_workspace_free(work);
	network_strain_half_fft_free(network_strain);
	Detector_Network_free(net);

	return 0;
}

#ifdef HAVE_CONFIG_H
	#include "config.h"
#endif

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
#include "inspiral_signal.h"
#include "inspiral_chirp_factors.h"
#include "inspiral_stationary_phase.h"
#include "signal.h"
#include "inspiral_network_statistic.h"
#include "sampling_system.h"

#include "inspiral_pso_fitness.h"
#include "random.h"
#include "options.h"
#include "settings_file.h"

#include "spectral_density.h"
#include "detector_mapping.h"

#include "hdf5_file.h"

//#include <hdf5/hdf5.h>
//#include <hdf5/hdf5_hl.h>

/* Compute the chirp factors so that we know the true chirp times, and save them to file. */
void save_true_signal(double f_low, source_t *source, char *filename) {
	inspiral_chirp_factors_t temp_chirp;
	CF_compute(f_low, source, &temp_chirp);
	FILE *true_parameters = fopen(filename, "w");
	fprintf(true_parameters, "RA DEC CHIRP_TIME_0 CHIRP_TIME_1_5 NETWORK_SNR\n");
	fprintf(true_parameters, "%0.21f %0.21f %0.21f %0.21f %0.21f",
			source->sky.ra, source->sky.dec, temp_chirp.ct.chirp_time0, temp_chirp.ct.chirp_time1_5, source->snr);
	fclose(true_parameters);
}

void append_index_to_prefix(char* buff, size_t buff_len, const char *prefix, size_t index) {
	memset(buff, '\0', buff_len * sizeof(char));
	sprintf(buff, "%s%lu", prefix, index+1);
}

void hdf5_save_simulated_data( size_t len_template, double *template, const char *hdf5_noise_filename, const char *hdf5_output_filename ) {
	//hid_t file_id, dataset_id, dspace_id, output_file_id;
	//herr_t status;
	char dset_name[255];

	size_t strain_len = hdf5_get_dataset_array_length( hdf5_noise_filename, "/strain/Strain_1" );
	printf("The strain length is %lu\n", strain_len);

	size_t num_strains = hdf5_get_num_strains( hdf5_noise_filename );
	printf("There are %zu strains in the file (%s).\n", num_strains, hdf5_noise_filename);

	printf("Checking that the noise in the file matches the length of the template... ");
	if (len_template != strain_len) {
		fprintf(stderr, "Error: The template length (%lu) doesn't equal the noise strain length (%lu). Aborting.\n",
				len_template, strain_len);
		abort();
	}
	printf("yes.\n");

	/* Create the groups in the output file */
	hdf5_create_group( hdf5_output_filename, "/signal");
	hdf5_create_group( hdf5_output_filename, "/noise");
	hdf5_create_group( hdf5_output_filename, "/strain");

	size_t i, j;
	for (i = 0; i < num_strains; i++) {
		/* Read the input LIGO noise strain */
		append_index_to_prefix(dset_name, 255, "/strain/Strain_", i);
		size_t strain_len = hdf5_get_dataset_array_length( hdf5_noise_filename, "/strain/Strain_1" );
		double *noise = (double*) malloc( strain_len * sizeof(double) );
		if (noise == NULL) {
			fprintf(stderr, "Error. Unable to allocate memory to read the noise. Aborting.\n");
			abort();
		}
		hdf5_load_array( hdf5_noise_filename, dset_name, noise);

		double *template_plus_noise = (double*) malloc( len_template * sizeof(double));
		if (template_plus_noise == NULL) {
			fprintf(stderr, "Error. Unable to allocate memory for the template + noise array. Aborting.\n");
			abort();
		}

		/* Add the signal to the strain */
		for (j = 0; j < len_template; j++) {
			template_plus_noise[j] = template[j] + noise[j];
		}

		/* Save the template */
		append_index_to_prefix(dset_name, 255, "Signal_", i);
		hdf5_save_array( hdf5_output_filename, "/signal", dset_name, len_template, template );

		/* Save the noise */
		append_index_to_prefix(dset_name, 255, "Noise_", i);
		hdf5_save_array( hdf5_output_filename, "/noise", dset_name, len_template, noise );

		/* Save the signal = template + noise */
		append_index_to_prefix(dset_name, 255, "Strain_", i);
		hdf5_save_array( hdf5_output_filename, "/strain", dset_name, len_template, template_plus_noise);

		free(noise);
		free(template_plus_noise);
	}

	printf("Finished writing to the HDF5 file.\n\n");
}

int main(int argc, char* argv[]) {
	size_t i, j;

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
	if ((argc - last_index) < 2) {
		printf("last_index = %d\n", last_index);
		printf("argc = %d\n", argc);
		printf("Error: Must supply [input settings file] [detector realizations mapping file]\n");
		exit(-1);
	}

	char* arg_settings_file = argv[last_index++];
	char* arg_detector_mappings_file = argv[last_index++];

	/* Load the general Settings */
	settings_file_t *settings_file = settings_file_open(arg_settings_file);
	if (settings_file == NULL) {
		printf("Error opening the settings file (%s). Aborting.\n", arg_settings_file);
		abort();
	}

	const double f_low = atof(settings_file_get_value(settings_file, "f_low"));
	const double f_high = atof(settings_file_get_value(settings_file, "f_high"));
	/*const double sampling_frequency = atof(settings_file_get_value(settings_file, "sampling_frequency"));*/
	//const double num_time_samples = atoi(settings_file_get_value(settings_file, "num_time_samples"));

	settings_file_close(settings_file);


	/* Initialize the detector network based on the detector mappings file. */
	detector_network_t *net = Detector_Network_load(arg_detector_mappings_file, f_low, f_high);

	/* For each detector, read in the PSD from file, and compute interpolated ASD for the analysis. */
	for (i = 0; i < net->num_detectors; i++) {
		/* Create the output file */
		char output_filename[255];
		memset(output_filename, '\0', 255 * sizeof(char));
		sprintf(output_filename, "%s.hdf5", net->detector[i]->name);
		hdf5_create_file( output_filename );

		/* Save the PSD to file */
		PSD_save( output_filename, net->detector[i]->psd );
	}

	/* Generate the templates for all the detectors composing the network */
	detector_network_mapping_t *dmap = Detector_Network_Mapping_load( arg_detector_mappings_file );
	size_t num_time_samples = hdf5_get_num_time_samples( dmap->data_filenames[0] );

	network_strain_half_fft_t *network_strain = inspiral_template(f_low, f_high, num_time_samples, net, &source);
	for (i = 0; i < net->num_detectors; i++) {
		strain_half_fft_t *one_sided = network_strain->strains[i];

		/* multiply by the ASD before performing the IFFT */
		for (j = 0; j < net->detector[i]->asd->len; j++) {
			one_sided->half_fft[j] = gsl_complex_mul_real(one_sided->half_fft[j], net->detector[i]->asd->asd[j]);
		}

		/* Form the two-sided template so we can take the inverse FFT */
		printf("Forming the two-sided template from the one-sided version.\n");
		gsl_complex *two_sided = (gsl_complex*) malloc( one_sided->full_len * sizeof(gsl_complex) );
		if (two_sided == NULL) {
			fprintf(stderr, "Error. Unable to allocate memory for the two-sided array. Aborting.\n");
			abort();
		}
		fprintf(stderr, "one-sided->full_len = %lu\n", one_sided->full_len);
		fprintf(stderr, "one-sided->half_fft_len = %lu\n", one_sided->half_fft_len);


		SS_make_two_sided (one_sided->half_fft_len, one_sided->half_fft, one_sided->full_len, two_sided);

		gsl_fft_complex_wavetable *fft_wavetable;
		gsl_fft_complex_workspace *fft_workspace;
		fft_wavetable = gsl_fft_complex_wavetable_alloc( one_sided->full_len );
		fft_workspace = gsl_fft_complex_workspace_alloc( one_sided->full_len );

		printf("Creating strided array for the ifft.\n");
		double *template_ifft = (double*) malloc( 2 * one_sided->full_len * sizeof(double) );
		if (template_ifft == NULL) {
			fprintf(stderr, "Error. Unable to allocate memory for the template_ifft. Aborting.\n");
			abort();
		}
		for (j = 0; j < one_sided->full_len; j++) {
			template_ifft[2*j + 0] = GSL_REAL( two_sided[j] );
			template_ifft[2*j + 1] = GSL_IMAG( two_sided[j] );
		}

		printf("Computing the IFFT of the template...\n");
		gsl_fft_complex_inverse( template_ifft, 1, one_sided->full_len, fft_wavetable, fft_workspace );

		double *template = (double*) malloc( one_sided->full_len * sizeof(double) );
		if (template == NULL) {
			fprintf(stderr, "Error. Unable to allocate memory for the template. Aborting.\n");
			abort();
		}
		for (j = 0; j < one_sided->full_len; j++) {
			/* Only save the real part */
			template[j] = template_ifft[2*j + 0];
		}

		gsl_fft_complex_workspace_free( fft_workspace );
		gsl_fft_complex_wavetable_free( fft_wavetable );

		fprintf(stderr, "Creating the output filename... ");
		const char *detector_name = net->detector[i]->name;
		char hdf5_output_filename[255];
		memset(hdf5_output_filename, '\0', 255 * sizeof(char));
		sprintf(hdf5_output_filename, "%s.hdf5", detector_name);
		fprintf(stderr, " (%s) done.\n", hdf5_output_filename);

		const char *hdf5_noise_filename = dmap->data_filenames[i];

		fprintf(stderr, "Saving the simulated data to file (%s)\n", hdf5_output_filename);
		hdf5_save_simulated_data( one_sided->full_len, template, hdf5_noise_filename, hdf5_output_filename );

		free(two_sided);
		free(template_ifft);
		free(template);
	}


	network_strain_half_fft_free(network_strain);

	Detector_Network_Mapping_close(dmap);

	Detector_Network_free(net);

	printf("program ended successfully.\n");

	return 0;
}

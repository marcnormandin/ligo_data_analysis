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

#include "common.h"
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
#include "sampling_system.h"

#include "inspiral_pso_fitness.h"
#include "inspiral_signal.h"
#include "random.h"
#include "options.h"
#include "settings.h"

#include <hdf5.h>
#include <hdf5_hl.h>


/* Compute the chirp factors so that we know the true chirp times, and save them to file. */
void save_true_signal(double f_low, source_t *source, char *filename) {
	chirp_factors_t temp_chirp;
	CF_compute(f_low, source, &temp_chirp);
	FILE *true_parameters = fopen(filename, "w");
	fprintf(true_parameters, "RA DEC CHIRP_TIME_0 CHIRP_TIME_1_5 NETWORK_SNR\n");
	fprintf(true_parameters, "%0.21f %0.21f %0.21f %0.21f %0.21f",
			source->sky.ra, source->sky.dec, temp_chirp.ct.chirp_time0, temp_chirp.ct.chirp_time1_5, source->snr);
	fclose(true_parameters);
}

size_t hdf5_get_strain_length( const char *hdf_filename ) {
	hid_t file_id, dataset_id, dspace_id;
	herr_t status;

	/* Open the file */
	file_id = H5Fopen( hdf_filename, H5F_ACC_RDWR, H5P_DEFAULT);
	if (file_id < 0) {
		fprintf(stderr, "Error opening hdf5 file (%s) for reading. Aborting.\n", hdf_filename);
	}

	/* Read the dataset */
	const char* dataset_name = "/strain/Strain_1";
	dataset_id = H5Dopen2(file_id, "/strain/Strain_1", H5P_DEFAULT);
	if (dataset_id < 0) {
		fprintf(stderr, "Error opening the dataset (%s) from the file (%s). Aborting.\n",
				dataset_name, hdf_filename);
		abort();
	}

	/* Get the length of the dataset */
	dspace_id = H5Dget_space(dataset_id);
	int ndims = H5Sget_simple_extent_ndims(dspace_id);
	hssize_t len = H5Sget_simple_extent_npoints(dspace_id);
	H5Fclose(file_id);

	return len;
}

size_t hdf5_get_num_strains( const char* hdf_filename ) {
	hid_t file_id, dataset_id, dspace_id;
	herr_t status;

	/* Open the file */
	file_id = H5Fopen( hdf_filename, H5F_ACC_RDWR, H5P_DEFAULT);
	if (file_id < 0) {
		fprintf(stderr, "Error opening the hdf5 file (%s). Aborting.\n",
				hdf_filename);
		abort();
	}

	/* Read the attribute */
	double num;
	const char *attribute_name = "num_strains";
	status = H5LTget_attribute_double( file_id, "/", attribute_name, &num);
	if (status < 0) {
		fprintf(stderr, "Error reading the attribute (%s) from the hdf5 file (%s). Aborting.\n",
				attribute_name, hdf_filename);
		abort();
	}

	H5Fclose(file_id);

	return (int)num;
}

void hdf5_save_simulated_data( size_t len_template, double *template, char *hdf5_noise_filename, char *hdf5_output_filename ) {
	hid_t file_id, dataset_id, dspace_id, output_file_id;
	herr_t status;

	size_t strain_len = hdf5_get_strain_length( hdf5_noise_filename );
	printf("The strain length is %lu\n", strain_len);

	size_t num_strains = hdf5_get_num_strains( hdf5_noise_filename );
	printf("There are %zu strains in the file (%s).\n", num_strains, hdf5_noise_filename);

	/* TEMP HACK */
	len_template--;

	printf("Checking that the noise in the file matches the length of the template... ");
	if (len_template != strain_len) {
		fprintf(stderr, "Error: The template length (%lu) doesn't equal the noise strain length (%lu). Aborting.\n",
				len_template, strain_len);
		abort();
	}
	printf("yes.\n");

	/* Open the input file */
	printf("Opening %s for reading.\n", hdf5_noise_filename);
	file_id = H5Fopen( hdf5_noise_filename, H5F_ACC_RDWR, H5P_DEFAULT);
	if (file_id < 0) {
		fprintf(stderr, "Error: Unable to open the HDF5 file (%s). Aborting.\n", hdf5_noise_filename);
		abort();
	}

	/* Open the output file */
	/*output_file_id = H5Fcreate( hdf5_output_filename, H5F_ACC_RDWR, H5P_DEFAULT);*/
	printf("Creating %s HDF5 file for writing.\n", hdf5_output_filename);
	output_file_id = H5Fopen( hdf5_output_filename, H5F_ACC_RDWR, H5P_DEFAULT);
	/*output_file_id = H5Fcreate( hdf5_output_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);*/
	if (output_file_id < 0) {
		fprintf(stderr, "Error opening the output hdf5 file (%s). Aborting.\n",
				hdf5_output_filename);
		abort();
	}

	/* Create the groups in the output file */
	hid_t template_group_id = H5Gcreate(output_file_id, "/template", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t noise_group_id = H5Gcreate(output_file_id, "/noise", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	hid_t signal_group_id = H5Gcreate(output_file_id, "/strain", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	char *dataset_prefix = "/strain/Strain_";

	size_t i, j;
	for (i = 0; i < num_strains; i++) {
		/* Form the strain name */
		char dataset_name[255];
		memset(dataset_name, '\0', 255 *sizeof(char));
		sprintf(dataset_name, "%s%lu", dataset_prefix, i+1);

		/* Read the noise strain */
		double *noise = (double*) malloc( strain_len * sizeof(double) );
		if (noise == NULL) {
			fprintf(stderr, "Error. Unable to allocate memory to read the noise. Aborting.\n");
			abort();
		}
		status = H5LTread_dataset_double( file_id, dataset_name, noise );
		if (status < 0) {
			fprintf(stderr, "Error reading the dataset (%s) from the file (%s). Aborting.\n",
					dataset_name, hdf5_noise_filename);
			abort();
		}

		double *template_plus_noise = (double*) malloc( len_template * sizeof(double));
		if (template_plus_noise == NULL) {
			fprintf(stderr, "Error. Unable to allocate memory for the template + noise array. Aborting.\n");
			abort();
		}

		/* Add the signal to the strain */
		for (j = 0; j < strain_len; j++) {
			template_plus_noise[j] = template[j] + noise[j];
		}

		/* Save the template */
		const char *template_dset_name_prefix = "Template_";
		char template_dset_name[255];
		memset(template_dset_name, '\0', 255 * sizeof(char));
		sprintf(template_dset_name, "%s%lu", template_dset_name_prefix, i+1);
		hsize_t dims[1];
		dims[0] = len_template;
		status = H5LTmake_dataset_double ( template_group_id, template_dset_name, 1, dims, template );
		if (status < 0) {
			fprintf(stderr, "Error saving the dataset (%s) to the hdf5 file (%s). Aborting.\n",
					template_dset_name, hdf5_output_filename);
			abort();
		}

		/* Save the noise */
		const char *noise_dset_name_prefix = "Noise_";
		char noise_dset_name[255];
		memset(noise_dset_name, '\0', 255 * sizeof(char));
		sprintf(noise_dset_name, "%s%lu", noise_dset_name_prefix, i+1);
		status = H5LTmake_dataset_double ( noise_group_id, noise_dset_name, 1, dims, noise );
		if (status < 0) {
			fprintf(stderr, "Error saving the dataset (%s) to the hdf5 file (%s). Aborting.\n",
					noise_dset_name, hdf5_output_filename);
			abort();
		}

		/* Save the signal = template + noise */
		const char *signal_dset_name_prefix = "Strain_";
		char signal_dset_name[255];
		memset(signal_dset_name, '\0', 255 * sizeof(char));
		sprintf(signal_dset_name, "%s%lu", signal_dset_name_prefix, i+1);
		status = H5LTmake_dataset_double ( signal_group_id, signal_dset_name, 1, dims, template_plus_noise );
		if (status < 0) {
			fprintf(stderr, "Error saving the dataset (%s) to the hdf5 file (%s). Aborting.\n",
					signal_dset_name, hdf5_output_filename);
			abort();
		}

		free(noise);
		free(template_plus_noise);
	}

	status = H5Gclose (template_group_id);
	status = H5Gclose (noise_group_id);
	status = H5Gclose (signal_group_id);

	H5Fclose(output_file_id);
	H5Fclose(file_id);

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
	const double sampling_frequency = atof(settings_file_get_value(settings_file, "sampling_frequency"));
	//const double num_time_samples = atoi(settings_file_get_value(settings_file, "num_time_samples"));

	settings_file_close(settings_file);

	/* open the mapping file */
	settings_file = settings_file_open(arg_detector_mappings_file);
	if (settings_file == NULL) {
		printf("Error opening the detector mappings file (%s). Aborting.\n", arg_detector_mappings_file);
		abort();
	}

	/* Initialize the detector network based on the detector mappings file. */
	int num_detectors = settings_file_num_settings(settings_file);
	printf("Initializing a detector network composed of %d detectors.\n", num_detectors);
	detector_network_t net;
	Alloc_Detector_Network(num_detectors, &net);
	for (i = 0; i < num_detectors; i++) {
		detector_t *det = net.detector[i];
		char *dname = settings_file_get_key_by_index(settings_file, i);
		Detector_init_name( dname, det);
	}


	char *hdf_filename = settings_file_get_value(settings_file, settings_file_get_key_by_index(settings_file, 0));
	int num_time_samples = hdf5_get_strain_length( hdf_filename );
	printf("Simulated data will contain %d time samples.\n", num_time_samples);
	/*int num_strains = hdf5_get_num_strains ( hdf_filename );
	printf("There are %d strains in the file.\n", num_strains);*/

	/*
	strain_t *irregular_strain = Strain_readFromFile("strain.txt");
	strain_t *hack = Strain_simulated(irregular_strain, f_low, f_high, sampling_frequency, num_time_samples);*/

	strain_t **strain = (strain_t**) malloc( num_detectors * sizeof(strain_t*) );

	/* HUGE TEMP HACK */
	for (i = 0; i < num_detectors; i++) {
		/*strain[i] = Strain_simulated(f_low, f_high, sampling_frequency, num_time_samples);*/
		//strain[i] = Strain_malloc(ST_ONE_SIDED, hack->len);
		//memcpy(strain[i]->freq, hack->freq, hack->len * sizeof(double));

		hid_t file_id;

		/* Open the input file to get the PSD */
		char *psd_filename = settings_file_get_value(settings_file, settings_file_get_key_by_index(settings_file, i));
		printf("Opening %s to get the PSD.\n", psd_filename);
		file_id = H5Fopen( psd_filename, H5F_ACC_RDWR, H5P_DEFAULT);
		if (file_id < 0) {
			fprintf(stderr, "Error: Unable to open the HDF5 file (%s). Aborting.\n", psd_filename);
			abort();
		}

		/* Read the PSD, but it has two many samples */
		printf("Reading the PSD..\n");
		strain_t* strain_two_sided = Strain_malloc(ST_TWO_SIDED, 2097152);
		herr_t status = H5LTread_dataset_double( file_id, "/psd/PSD", strain_two_sided->strain );
		if (status < 0) {
			fprintf(stderr, "Error reading the dataset (%s) from the file (%s). Aborting.\n",
					"/psd/PSD", psd_filename);
			abort();
		}
		H5Fclose( file_id );

		SS_frequency_array(4096.0, 2097152, strain_two_sided->freq);

		/* Downsample the PSD */
		printf("Downsampling the PSD for use in the network statistic that is used to scale the signals by network SNR... ");
		strain[i] = Strain_simulated(strain_two_sided, f_low, f_high, sampling_frequency, num_time_samples);
		printf("done.\n");





		/* Save the PSD to the output file. */
		char output_filename[255];
		memset(output_filename, '\0', 255 * sizeof(char));
		sprintf(output_filename, "%s.hdf5", settings_file_get_key_by_index(settings_file, i));
		printf("Creating %s HDF5 file for writing the PSD.\n", output_filename);
		hid_t output_file_id = H5Fcreate( output_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		if (output_file_id < 0) {
			fprintf(stderr, "Error opening the output hdf5 file (%s). Aborting.\n",
					output_filename);
			abort();
		}

		/* Create the groups in the output file */
		hid_t psd_group_id = H5Gcreate(output_file_id, "/psd", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

		hsize_t dims[1];
		dims[0] = strain[i]->len;
		status = H5LTmake_dataset_double ( psd_group_id, "PSD", 1, dims, strain[i]->strain );
		if (status < 0) {
			fprintf(stderr, "Error saving the dataset (%s) to the hdf5 file (%s). Aborting.\n",
					"PSD", output_filename);
			abort();
		}
		printf("Downsampled PSD saved to HDF5 file (%s).\n", output_filename);
		H5Gclose( psd_group_id );
		H5Fclose(output_file_id);
	}

	/* Generate the templates for all the detectors composing the network */
	inspiral_signal_half_fft_t **templates = inspiral_template(f_low, f_high, &net, strain, &source);
	for (i = 0; i < num_detectors; i++) {
		inspiral_signal_half_fft_t *one_sided = templates[i];

		/* Form the two-sided template so we can take the inverse FFT */
		printf("Forming the two-sided template from the one-sided version.\n");
		gsl_complex *two_sided = (gsl_complex*) malloc( one_sided->full_len * sizeof(gsl_complex) );
		if (two_sided == NULL) {
			fprintf(stderr, "Error. Unable to allocate memory for the two-sided array. Aborting.\n");
			abort();
		}
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

		char *detector_name = settings_file_get_key_by_index(settings_file, i);
		char *hdf5_noise_filename = settings_file_get_value(settings_file, detector_name);
		char hdf5_output_filename[255];
		memset(hdf5_output_filename, '\0', 255 * sizeof(char));
		sprintf(hdf5_output_filename, "%s.hdf5", detector_name);

		printf("Saving the simulated data to file (%s)\n", hdf5_output_filename);
		hdf5_save_simulated_data( one_sided->full_len, template, hdf5_noise_filename, hdf5_output_filename );

		free(two_sided);
		free(template_ifft);
		free(template);
	}


	/* Free the data */
	for (i = 0; i < net.num_detectors; i++) {
		inspiral_signal_half_fft_free(templates[i]);
	}
	free(templates);

	Free_Detector_Network(&net);

	for (i = 0; i < num_detectors; i++) {
		char strain_filename[255];
		memset(strain_filename, '\0', 255 *sizeof(char));
		sprintf(strain_filename, "%s.strain", settings_file_get_key_by_index(settings_file, i));
		Strain_saveToFile(strain_filename, strain[i]);
		Strain_free(strain[i]);
	}
	free(strain);

	settings_file_close(settings_file);

	/*Strain_free(hack);*/

	printf("program ended successfully.\n");

	return 0;
}

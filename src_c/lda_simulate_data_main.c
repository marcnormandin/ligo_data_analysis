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
#include <gsl/gsl_interp.h>


#include "common.h"
#include "detector_antenna_patterns.h"
#include "inspiral_chirp.h"
#include "detector.h"
#include "detector_network.h"
#include "inspiral_source.h"
#include "inspiral_stationary_phase.h"
#include "signal.h"
#include "inspiral_network_statistic.h"
#include "sampling_system.h"

#include "inspiral_pso_fitness.h"
#include "inspiral_template.h"
#include "random.h"
#include "options.h"
#include "settings.h"

#include "spectral_density.h"

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

size_t hdf5_get_dataset_array_length( const char *hdf_filename, const char* dataset_name ) {
	hid_t file_id, dataset_id, dspace_id;
	herr_t status;

	/* Open the file */
	file_id = H5Fopen( hdf_filename, H5F_ACC_RDWR, H5P_DEFAULT);
	if (file_id < 0) {
		fprintf(stderr, "Error opening hdf5 file (%s) for reading. Aborting.\n", hdf_filename);
	}

	/* Read the dataset */
	dataset_id = H5Dopen2(file_id, dataset_name, H5P_DEFAULT);
	if (dataset_id < 0) {
		fprintf(stderr, "Error opening the dataset (%s) from the file (%s). Aborting.\n",
				dataset_name, hdf_filename);
		abort();
	}

	/* Get the length of the dataset */
	dspace_id = H5Dget_space(dataset_id);
	int ndims = H5Sget_simple_extent_ndims(dspace_id);
	hssize_t len = H5Sget_simple_extent_npoints(dspace_id);

	H5Dclose(dataset_id);
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

void hdf5_load_array( const char *hdf_filename, const char *dataset_name, double *data) {
	hid_t file_id;
	herr_t status;

	file_id = H5Fopen( hdf_filename, H5F_ACC_RDWR, H5P_DEFAULT );
	if (file_id < 0) {
		fprintf(stderr, "Error: Unable to open the HDF5 file (%s). Aborting.\n", hdf_filename);
		abort();
	}

	status = H5LTread_dataset_double( file_id, dataset_name, data );
	if (status < 0) {
		fprintf(stderr, "Error reading the dataset (%s) from the file (%s). Aborting.\n",
				dataset_name, hdf_filename);
		abort();
	}

	H5Fclose( file_id );
}

psd_t* hdf5_load_psd( const char *hdf_filename ) {
	size_t len_psd;

	len_psd = hdf5_get_dataset_array_length( hdf_filename, "/psd/PSD" );

	psd_t* psd = PSD_malloc ( len_psd );

	hdf5_load_array( hdf_filename, "/psd/PSD", psd->psd );
	hdf5_load_array( hdf_filename, "/psd/Freq", psd->f );

	psd->type = PSD_ONE_SIDED;

	return psd;
}

void hdf5_save_psd( const char *hdf_filename, psd_t *psd ) {
	/* Save the PSD to the output file. */
	hid_t output_file_id = H5Fopen( hdf_filename, H5F_ACC_RDWR, H5P_DEFAULT );
	if (output_file_id < 0) {
		fprintf(stderr, "Error opening the output hdf5 file (%s). Aborting.\n",
				hdf_filename);
		abort();
	}

	/* Create the group in the output file */
	hid_t psd_group_id = H5Gcreate(output_file_id, "/psd", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	if (psd_group_id < 0) {
		fprintf(stderr, "Error creating the group (%s) in the file (%s). Aborting.\n",
				"/psd", hdf_filename);
		abort();
	}

	hsize_t dims[1];
	dims[0] = psd->len;
	herr_t status;

	/* Save the PSD */
	status = H5LTmake_dataset_double ( psd_group_id, "PSD", 1, dims, psd->psd );
	if (status < 0) {
		fprintf(stderr, "Error saving the dataset (%s) to the hdf5 file (%s). Aborting.\n",
				"PSD", hdf_filename);
		abort();
	}

	/* Save the frequencies */
	status = H5LTmake_dataset_double ( psd_group_id, "Freq", 1, dims, psd->f );
	if (status < 0) {
		fprintf(stderr, "Error saving the dataset (%s) to the hdf5 file (%s). Aborting.\n",
				"Freq", hdf_filename);
		abort();
	}

	H5Gclose( psd_group_id );
	H5Fclose(output_file_id);
}

void hdf5_save_simulated_data( size_t len_template, double *template, char *hdf5_noise_filename, char *hdf5_output_filename ) {
	hid_t file_id, dataset_id, dspace_id, output_file_id;
	herr_t status;

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

void hdf5_create_file( const char* hdf_filename ) {
	hid_t file_id = H5Fcreate( hdf_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (file_id < 0) {
		fprintf(stderr, "Error. Unable to create the HDF5 file (%s). Aborting.\n",
				hdf_filename);
		abort();
	}
	H5Fclose(file_id);
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
	int num_time_samples = hdf5_get_dataset_array_length( hdf_filename, "/strain/Strain_1" );
	printf("Simulated data will contain %d time samples.\n", num_time_samples);

	size_t half_size = SS_half_size(num_time_samples);
	printf("The half size is: %lu\n", half_size);

	asd_t **net_asd = (asd_t**) malloc( num_detectors * sizeof(asd_t*) );

	/* For each detector, read in the PSD from file, and compute interpolated ASD for the analysis. */
	for (i = 0; i < num_detectors; i++) {
		/* Get the input filename */
		char *hdf_filename = settings_file_get_value(settings_file, settings_file_get_key_by_index(settings_file, i));

		/* Load the ASD from file */
		psd_t *psd_unprocessed = hdf5_load_psd( hdf_filename );
		asd_t *asd_unprocessed = ASD_malloc( psd_unprocessed->len );
		ASD_init_from_psd( psd_unprocessed, asd_unprocessed );
		PSD_free(psd_unprocessed);

		/* Set the frequencies we want to use. */
		net_asd[i] = ASD_malloc( half_size );
		asd_t *asd = net_asd[i];
		asd->type = ASD_ONE_SIDED;
		SS_frequency_array(sampling_frequency, num_time_samples, asd->len, asd->f);

		/* Interpolate the ASD to the desired frequencies. */
		gsl_interp* interp = gsl_interp_alloc(gsl_interp_linear, asd_unprocessed->len);
		gsl_interp_accel* acc = gsl_interp_accel_alloc();
		gsl_interp_init(interp, asd_unprocessed->f, asd_unprocessed->asd, asd_unprocessed->len);
		for (j = 0; j < asd->len; j++) {
			double asd_interpolated = gsl_interp_eval(interp, asd_unprocessed->f, asd_unprocessed->asd, asd->f[j], acc);
			asd->asd[j] = asd_interpolated;
		}
		gsl_interp_accel_free(acc);
		gsl_interp_free(interp);

		/* Convert the ASD to a PSD */
		psd_t *psd = PSD_malloc( asd->len );
		PSD_init_from_asd( asd, psd );

		/* Create the output file */
		char output_filename[255];
		memset(output_filename, '\0', 255 * sizeof(char));
		sprintf(output_filename, "%s.hdf5", settings_file_get_key_by_index(settings_file, i));
		hdf5_create_file( output_filename );

		/* Save the PSD to file */
		hdf5_save_psd( output_filename, psd );

		/* Clean up */
		PSD_free(psd);
	}

	/* Generate the templates for all the detectors composing the network */
	inspiral_template_half_fft_t **templates = inspiral_template(f_low, f_high, num_time_samples, &net, net_asd, &source);
	for (i = 0; i < num_detectors; i++) {
		inspiral_template_half_fft_t *one_sided = templates[i];

		/* Form the two-sided template so we can take the inverse FFT */
		printf("Forming the two-sided template from the one-sided version.\n");
		gsl_complex *two_sided = (gsl_complex*) malloc( one_sided->full_len * sizeof(gsl_complex) );
		if (two_sided == NULL) {
			fprintf(stderr, "Error. Unable to allocate memory for the two-sided array. Aborting.\n");
			abort();
		}
		fprintf(stderr, "one-sided->full_len = %lu\n", one_sided->full_len);

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

	/* Free the ASD */
	for (i = 0; i < num_detectors; i++) {
		ASD_free( net_asd[i] );
	}
	free(net_asd);

	settings_file_close(settings_file);

	printf("program ended successfully.\n");

	return 0;
}

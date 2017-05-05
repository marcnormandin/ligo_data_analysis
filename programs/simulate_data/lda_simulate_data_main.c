#ifdef HAVE_CONFIG_H
	#include "config.h"
#endif

#include <assert.h>
#include <stddef.h>
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
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>


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

#include "strain.h"

#include "hdf5_file.h"

//#include <hdf5/hdf5.h>
//#include <hdf5/hdf5_hl.h>

#define FILENAME_MAX_SIZE 255

typedef struct program_settings_s {
	source_t source;
	gslseed_t alpha_seed;
	double f_low, f_high;
	double sampling_frequency;
	size_t num_time_samples;
	char detector_mapping_filename[FILENAME_MAX_SIZE];
	char output_filename[FILENAME_MAX_SIZE];
	size_t num_realizations;

} program_settings_t;

/* Compute the chirp factors so that we know the true chirp times, and save them to file. */
void save_true_signal(double f_low, source_t *source, char *filename) {
	inspiral_chirp_factors_t temp_chirp;
	CF_compute_for_signal(f_low, source->m1, source->m2, source->time_of_arrival, &temp_chirp);
	FILE *true_parameters = fopen(filename, "w");
	fprintf(true_parameters, "RA DEC CHIRP_TIME_0 CHIRP_TIME_1_5 NETWORK_SNR\n");
	fprintf(true_parameters, "%0.21f %0.21f %0.21f %0.21f %0.21f",
			source->sky.ra, source->sky.dec, temp_chirp.ct.chirp_time0, temp_chirp.ct.chirp_time1_5, source->snr);
	fclose(true_parameters);
}

void append_index_to_prefix(char* buff, size_t buff_len, const char *prefix, size_t index) {
	memset(buff, '\0', buff_len * sizeof(char));
	sprintf(buff, "%s%lu", prefix, index);
}

void simulated_strain_file_create( const char *filename ) {
	hdf5_create_file( filename );
}

void simulated_strain_file_save_source( const char *output_filename, const source_t *source ) {
	hdf5_create_group( output_filename, "/source_parameters");

	hdf5_save_attribute_double( output_filename, "/source_parameters", "m1", 1, &source->m1 );
	hdf5_save_attribute_double( output_filename, "/source_parameters", "m2", 1, &source->m2 );
	hdf5_save_attribute_double( output_filename, "/source_parameters", "time_of_arrival", 1, &source->time_of_arrival );
	hdf5_save_attribute_double( output_filename, "/source_parameters", "right_ascension", 1, &source->sky.ra );
	hdf5_save_attribute_double( output_filename, "/source_parameters", "declination", 1, &source->sky.dec );
	hdf5_save_attribute_double( output_filename, "/source_parameters", "polarization_angle", 1, &source->polarization_angle );
	hdf5_save_attribute_double( output_filename, "/source_parameters", "coalescence_phase", 1, &source->coalescence_phase );
	hdf5_save_attribute_double( output_filename, "/source_parameters", "inclination_angle", 1, &source->inclination_angle );
	hdf5_save_attribute_double( output_filename, "/source_parameters", "snr", 1, &source->snr );
}

void simulated_strain_file_save_settings( const char *output_filename, const program_settings_t *ps) {
	hdf5_save_attribute_ulong( output_filename, "/", "alpha_seed", 1, &ps->alpha_seed );
	hdf5_save_attribute_ulong( output_filename, "/", "num_time_samples", 1, &ps->num_time_samples );
	hdf5_save_attribute_double( output_filename, "/", "sampling_frequency", 1, &ps->sampling_frequency );
	hdf5_save_attribute_double( output_filename, "/", "f_low", 1, &ps->f_low );
	hdf5_save_attribute_double( output_filename, "/", "f_high", 1, &ps->f_high );
	hdf5_save_attribute_string( output_filename, "/", "detector_mapping_filename", ps->detector_mapping_filename );
	hdf5_save_attribute_ulong( output_filename, "/", "num_realizations", 1, &ps->num_realizations );

	simulated_strain_file_save_source( output_filename, &ps->source );
}

void simulated_strain_file_save_detector( const char *output_filename, const detector_t* detector, size_t detector_num ) {
	char group_name[255];
	append_index_to_prefix(group_name, 255, "/detector_", detector_num );
	hdf5_create_group( output_filename, group_name );
	//hdf5_save_attribute_ulong( output_filename, group_name, "id", 1, &detector->id );
	hdf5_save_attribute_string( output_filename, group_name, "name", detector->name );
	hdf5_save_attribute_gsl_vector( output_filename, group_name, "location", detector->location );
	hdf5_save_attribute_gsl_vector( output_filename, group_name, "arm_x", detector->arm_x );
	hdf5_save_attribute_gsl_vector( output_filename, group_name, "arm_y", detector->arm_y );

	// save the PSD
	char psd_group[255];
	memset(psd_group, '\0', 255 *sizeof(char));
	sprintf(psd_group, "%s%s", group_name, "/psd");
	hdf5_create_group( output_filename, psd_group);
	hdf5_save_array( output_filename, psd_group, "PSD", detector->psd->len, detector->psd->psd );
	hdf5_save_array( output_filename, psd_group, "Freq", detector->psd->len, detector->psd->f );

/* Todo
	gsl_matrix *detector_tensor;

	psd_t *psd;

	asd_t *asd; */

}

void simulated_strain_file_save_detector_network( const char *output_filename, const detector_network_t *dnet ) {
	size_t i;
	for (i = 0; i < dnet->num_detectors; i++) {
		simulated_strain_file_save_detector( output_filename, dnet->detector[i], i+1 );
	}
	hdf5_save_attribute_ulong( output_filename, "/", "num_detectors", 1, &dnet->num_detectors );
}

void simulated_strain_file_save_detector_signal( const char *output_filename, strain_t *signal, size_t detector_num ) {
	char buff[255];
	memset(buff, '\0', 255 * sizeof(char));

	sprintf(buff, "detector_%d/signal", detector_num);
	hdf5_create_group( output_filename, buff );
	hdf5_save_array( output_filename, buff, "Signal", signal->num_time_samples, signal->samples );
}

void hdf5_save_simulated_data( size_t len_template, double *signal, const char *hdf5_noise_filename, const char *hdf5_output_filename ) {
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
			template_plus_noise[j] = signal[j] + noise[j];
		}

		/* Save the template */
		append_index_to_prefix(dset_name, 255, "Signal_", i);
		hdf5_save_array( hdf5_output_filename, "/signal", dset_name, len_template, signal );

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

void program_settings_init(int argc, char *argv[], program_settings_t *ps) {
	int last_index;

	/* Allow options to override the defaults settings */
	Source_load_testsource(&ps->source);
	ps->alpha_seed = 0;
	if (process_command_options(argc, argv, &ps->source, &ps->alpha_seed, &last_index) != 0) {
		exit(0);
	}

	/* somehow these need to be set */
	if ((argc - last_index) < 3) {
		printf("last_index = %d\n", last_index);
		printf("argc = %d\n", argc);
		printf("Error: Must supply [settings file] [detector mapping file] [output filename]\n");
		exit(-1);
	}
	char* arg_settings_file = argv[last_index++];
	char* arg_detector_mappings_file = argv[last_index++];
	char* arg_output_filename = argv[last_index++];

	/* Load the general Settings */
	settings_file_t *settings_file = settings_file_open(arg_settings_file);
	if (settings_file == NULL) {
		printf("Error opening the settings file (%s). Aborting.\n", arg_settings_file);
		abort();
	}
	ps->f_low = atof(settings_file_get_value(settings_file, "f_low"));
	ps->f_high = atof(settings_file_get_value(settings_file, "f_high"));
	ps->num_time_samples = atoi(settings_file_get_value(settings_file, "num_time_samples"));
	ps->sampling_frequency = atof(settings_file_get_value(settings_file, "sampling_frequency"));
	ps->num_realizations = atoi(settings_file_get_value(settings_file, "num_realizations"));

	memset( ps->detector_mapping_filename, '\0', FILENAME_MAX_SIZE * sizeof(char) );
	strncpy( ps->detector_mapping_filename, arg_detector_mappings_file, FILENAME_MAX_SIZE );

	memset( ps->output_filename, '\0', FILENAME_MAX_SIZE * sizeof(char) );
	strncpy( ps->output_filename, arg_output_filename, FILENAME_MAX_SIZE );

	settings_file_close(settings_file);
}

void simulate( program_settings_t *ps, detector_network_t *net) {
	size_t i;
	size_t j;
	size_t k;
	size_t l;

	/* Simulate the signals */
	network_strain_half_fft_t *network_strain_half_fft = inspiral_template(ps->f_low, ps->f_high, ps->num_time_samples, net, &ps->source);

	strain_t **signals = (strain_t**) malloc ( net->num_detectors * sizeof(strain_t*) );
	if (signals == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory for detector strains. Exiting.\n");
		exit(-1);
	}

	/* Compute the coloured signal strain series */
	for (i = 0; i < net->num_detectors; i++) {
		strain_half_fft_t *one_sided = network_strain_half_fft->strains[i];

		/* Colour the template. multiply by the ASD before performing the IFFT */
		for (j = 0; j < net->detector[i]->asd->len; j++) {
			one_sided->half_fft[j] = gsl_complex_mul_real(one_sided->half_fft[j], net->detector[i]->asd->asd[j]);
		}

		/* Make the two-sided fft template */
		strain_full_fft_t *two_sided = strain_half_to_full( one_sided );

		/* Compute the strain */
		signals[i] = strain_full_fft_to_strain( two_sided );

		/* Normalize by the number of samples */
		for (j = 0; j < ps->num_time_samples; j++) {
			signals[i]->samples[j] *= ps->num_time_samples;
		}

		strain_full_fft_free(two_sided);
	}

	/* random number generator */
	gsl_rng *rng = random_alloc( ps->alpha_seed );

	/* Compute the coloured noise strain series */
	gsl_fft_real_wavetable *fft_real_wavetable = gsl_fft_real_wavetable_alloc( ps->num_time_samples );
	gsl_fft_halfcomplex_wavetable *fft_complex_wavetable = gsl_fft_halfcomplex_wavetable_alloc( ps->num_time_samples );
	gsl_fft_real_workspace *fft_workspace = gsl_fft_real_workspace_alloc( ps->num_time_samples );
	double *noise = (double*) malloc( ps->num_time_samples * sizeof(double) );
	double *strain = (double*) malloc( ps->num_time_samples * sizeof(double) );
	for (i = 0; i < net->num_detectors; i++) {
		asd_t *asd_one_sided = net->detector[i]->asd;
		//double *asd_two_sided = (double*) malloc( ps->num_time_samples * sizeof(double) );
		//SS_make_two_sided_real( asd_one_sided->len, asd_one_sided->asd, ps->num_time_samples, asd_two_sided);

		char buff[255];
		memset(buff, '\0', 255 * sizeof(char));
		sprintf(buff, "detector_%d/noise", i+1);
		hdf5_create_group( ps->output_filename, buff );

		char buff3[255];
		memset(buff3, '\0', 255 * sizeof(char));
		sprintf(buff3, "detector_%d/strain", i+1);
		hdf5_create_group( ps->output_filename, buff3 );

		for (j = 0; j < ps->num_realizations; j++) {
			for (k = 0; k < ps->num_time_samples; k++) {
				noise[k] = gsl_ran_gaussian( rng, 1.0 );
			}

			gsl_fft_real_transform(noise, 1, ps->num_time_samples, fft_real_wavetable, fft_workspace);

			// DC term doesn't have an imaginary component
			noise[0] *= asd_one_sided->asd[0];

			size_t lu = SS_last_unique_index( ps->num_time_samples );
			if (SS_has_nyquist_term(ps->num_time_samples)) {
				lu--;
			}

			//fprintf(stderr, "last unique index = %d\n", lu);
			for (k = 1, l = 1; l <= lu; k+=2, l++) {
				double s = asd_one_sided->asd[l];
				noise[k+0] *= s;
				noise[k+1] *= s;
			}

			// If nyquist term is present, it doesn't have an imaginary component
			if (SS_has_nyquist_term(ps->num_time_samples)) {
				noise[ps->num_time_samples-1] *= asd_one_sided->asd[asd_one_sided->len-1];
			}

			gsl_fft_halfcomplex_inverse( noise, 1, ps->num_time_samples, fft_complex_wavetable, fft_workspace );

			char buff2[255];
			memset(buff2, '\0', 255 * sizeof(char));
			sprintf(buff2, "Noise_%d", j+1);
			hdf5_save_array( ps->output_filename, buff, buff2, ps->num_time_samples, noise );

			for (k = 0; k < ps->num_time_samples; k++) {
				strain[k] = signals[i]->samples[k] + noise[k];
			}

			char buff4[255];
			memset(buff4, '\0', 255 * sizeof(char));
			sprintf(buff4, "Strain_%d", j+1);
			hdf5_save_array( ps->output_filename, buff3, buff4, ps->num_time_samples, strain );
		}
	}
	gsl_fft_real_workspace_free( fft_workspace );
	gsl_fft_halfcomplex_wavetable_free( fft_complex_wavetable );
	gsl_fft_real_wavetable_free( fft_real_wavetable );
	free(strain);
	free(noise);

	random_free( rng );

	/* Save the signals */
	for (i = 0; i < net->num_detectors; i++) {
		simulated_strain_file_save_detector_signal( ps->output_filename, signals[i], i+1);
	}



	/* Free memory */
	for (i = 0; i < net->num_detectors; i++) {
		strain_free(signals[i]);
	}
	free(signals);

	network_strain_half_fft_free(network_strain_half_fft);
}

int main(int argc, char* argv[])
{
	size_t i, j;
	program_settings_t ps;
	program_settings_init( argc, argv, &ps );

	simulated_strain_file_create( ps.output_filename );
	simulated_strain_file_save_settings( ps.output_filename, &ps );

	detector_network_t *net = Detector_Network_load( ps.detector_mapping_filename,
			ps.num_time_samples, ps.sampling_frequency, ps.f_low, ps.f_high );
	simulated_strain_file_save_detector_network( ps.output_filename, net );

	simulate( &ps, net );

	Detector_Network_free( net );

	return 0;
}

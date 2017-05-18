/*
 * simulation_file.c
 *
 *  Created on: May 8, 2017
 *      Author: marcnormandin
 */

#include <stddef.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "simulation_file.h"
#include "settings_file.h"
#include "strain.h"
#include "spectral_density.h"
#include "sampling_system.h"
#include "hdf5_file.h"
#include "detector.h"
#include "detector_network.h"
#include "random.h"
#include "inspiral_chirp_factors.h"

void append_index_to_prefix(char* buff, size_t buff_len, const char *prefix, size_t index) {
	memset(buff, '\0', buff_len * sizeof(char));
	sprintf(buff, "%s%lu", prefix, index);
}


/* Compute the chirp factors so that we know the true chirp times, and save them to file. */
void simulated_strain_file_save_chirp_factors(const char *output_filename, const double f_low, const source_t *source ) {
	inspiral_chirp_factors_t temp_chirp;
	CF_compute_for_signal(f_low, source->m1, source->m2, source->time_of_arrival, &temp_chirp);

	hdf5_create_group( output_filename, "/chirp_factors" );

	hdf5_save_attribute_double( output_filename, "/chirp_factors", "total_mass", 1, &temp_chirp.total_mass );
	hdf5_save_attribute_double( output_filename, "/chirp_factors", "reduced_mass", 1, &temp_chirp.reduced_mass );
	hdf5_save_attribute_double( output_filename, "/chirp_factors", "chirp_mass", 1, &temp_chirp.chirp_mass );
	hdf5_save_attribute_double( output_filename, "/chirp_factors", "s_mass_ratio", 1, &temp_chirp.s_mass_ratio );
	hdf5_save_attribute_double( output_filename, "/chirp_factors", "multi_fac", 1, &temp_chirp.multi_fac );
	hdf5_save_attribute_double( output_filename, "/chirp_factors", "calculated_reduced_mass", 1, &temp_chirp.calculated_reduced_mass );
	hdf5_save_attribute_double( output_filename, "/chirp_factors", "calculated_total_mass", 1, &temp_chirp.calculated_total_mass );
	hdf5_save_attribute_double( output_filename, "/chirp_factors", "t_chirp", 1, &temp_chirp.t_chirp );
	hdf5_save_attribute_double( output_filename, "/chirp_factors", "s_mass_ratio_cal", 1, &temp_chirp.s_mass_ratio_cal );
	hdf5_save_attribute_double( output_filename, "/chirp_factors", "multi_fac_call", 1, &temp_chirp.multi_fac_cal );

	hdf5_save_attribute_double( output_filename, "/chirp_factors", "chirp_time_0", 1, &temp_chirp.ct.chirp_time0 );
	hdf5_save_attribute_double( output_filename, "/chirp_factors", "chirp_time_1", 1, &temp_chirp.ct.chirp_time1 );
	hdf5_save_attribute_double( output_filename, "/chirp_factors", "chirp_time_1_5", 1, &temp_chirp.ct.chirp_time1_5 );
	hdf5_save_attribute_double( output_filename, "/chirp_factors", "chirp_time_2", 1, &temp_chirp.ct.chirp_time2 );
	hdf5_save_attribute_double( output_filename, "/chirp_factors", "time_of_coalescence", 1, &temp_chirp.ct.tc );
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

void simulated_strain_file_save_settings( const char *output_filename, const simulation_settings_t *ps) {
	hdf5_save_attribute_ulong( output_filename, "/", "alpha_seed", 1, &ps->alpha_seed );
	hdf5_save_attribute_ulong( output_filename, "/", "num_time_samples", 1, &ps->num_time_samples );
	hdf5_save_attribute_double( output_filename, "/", "sampling_frequency", 1, &ps->sampling_frequency );
	hdf5_save_attribute_double( output_filename, "/", "f_low", 1, &ps->f_low );
	hdf5_save_attribute_double( output_filename, "/", "f_high", 1, &ps->f_high );
	hdf5_save_attribute_string( output_filename, "/", "detector_mapping_filename", ps->detector_mapping_filename );
	hdf5_save_attribute_ulong( output_filename, "/", "num_realizations", 1, &ps->num_realizations );

	simulated_strain_file_save_source( output_filename, &ps->source );
	simulated_strain_file_save_chirp_factors( output_filename, ps->f_low, &ps->source );
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

	sprintf(buff, "detector_%zu/signal", detector_num);
	hdf5_create_group( output_filename, buff );
	hdf5_save_array( output_filename, buff, "Signal", signal->num_time_samples, signal->samples );
}

void simulate( simulation_settings_t *ps, detector_network_t *net) {
	size_t i;
	size_t j;
	size_t k;
	size_t l;

	/* Simulate the signals */
	network_strain_half_fft_t *network_strain_half_fft = inspiral_template(ps->f_low, ps->f_high, ps->num_time_samples, net, &ps->source, ps->output_filename);

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

	double *noise = (double*) malloc( ps->num_time_samples * sizeof(double) );
	double *strain = (double*) malloc( ps->num_time_samples * sizeof(double) );

	for (i = 0; i < net->num_detectors; i++) {
		asd_t *asd_one_sided = net->detector[i]->asd;
		//double *asd_two_sided = (double*) malloc( ps->num_time_samples * sizeof(double) );
		//SS_make_two_sided_real( asd_one_sided->len, asd_one_sided->asd, ps->num_time_samples, asd_two_sided);

		char buff[255];
		memset(buff, '\0', 255 * sizeof(char));
		sprintf(buff, "detector_%zu/noise", i+1);
		hdf5_create_group( ps->output_filename, buff );

		char buff3[255];
		memset(buff3, '\0', 255 * sizeof(char));
		sprintf(buff3, "detector_%zu/strain", i+1);
		hdf5_create_group( ps->output_filename, buff3 );

		for (j = 0; j < ps->num_realizations; j++) {
			for (k = 0; k < ps->num_time_samples; k++) {
				noise[k] = gsl_ran_gaussian( rng, 1.0 );
			}

			/* Compute the coloured noise strain series */
			SS_colour_timeseries( net->detector[i]->psd, ps->num_time_samples, noise);

			char buff2[255];
			memset(buff2, '\0', 255 * sizeof(char));
			sprintf(buff2, "Noise_%zu", j+1);
			hdf5_save_array( ps->output_filename, buff, buff2, ps->num_time_samples, noise );

			for (k = 0; k < ps->num_time_samples; k++) {
				strain[k] = signals[i]->samples[k] + noise[k];
			}

			char buff4[255];
			memset(buff4, '\0', 255 * sizeof(char));
			sprintf(buff4, "Strain_%zu", j+1);
			hdf5_save_array( ps->output_filename, buff3, buff4, ps->num_time_samples, strain );
		}
	}

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


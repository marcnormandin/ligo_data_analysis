/*
 * simulation_settings.c
 *
 *  Created on: May 8, 2017
 *      Author: marcnormandin
 */

#include <stddef.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "settings_file.h"
#include "simulation_settings.h"
#include "simulation_options.h"

void simulation_settings_init(int argc, char *argv[], simulation_settings_t *ps) {
	int last_index;

	/* Allow options to override the defaults settings */
	Source_load_testsource(&ps->source);
	ps->alpha_seed = 0;
	if (process_command_options(argc, argv, &ps->source, &ps->alpha_seed, &last_index) != 0) {
		exit(0);
	}

	/* somehow these need to be set */
	if ((argc - last_index) < 3) {
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

	memset( ps->ns_timeseries_filename, '\0', FILENAME_MAX_SIZE * sizeof(char) );
	strncpy( ps->ns_timeseries_filename, "network_timeseries.hdf5", FILENAME_MAX_SIZE);

	memset( ps->output_filename, '\0', FILENAME_MAX_SIZE * sizeof(char) );
	strncpy( ps->output_filename, arg_output_filename, FILENAME_MAX_SIZE );

	settings_file_close(settings_file);
}


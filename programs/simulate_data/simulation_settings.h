/*
 * simulation_settings.h
 *
 *  Created on: May 8, 2017
 *      Author: marcnormandin
 */

#ifndef PROGRAMS_SIMULATE_DATA_SIMULATION_SETTINGS_H_
#define PROGRAMS_SIMULATE_DATA_SIMULATION_SETTINGS_H_

#include <stddef.h>

#include "inspiral_signal.h"
#include "random.h"

#define FILENAME_MAX_SIZE 255

typedef struct simulation_settings_s {
	source_t source;
	gslseed_t alpha_seed;
	double f_low, f_high;
	double sampling_frequency;
	size_t num_time_samples;
	char detector_mapping_filename[FILENAME_MAX_SIZE];
	char output_filename[FILENAME_MAX_SIZE];
	size_t num_realizations;

} simulation_settings_t;

void simulation_settings_init(int argc, char *argv[], simulation_settings_t *ps);



#endif /* PROGRAMS_SIMULATE_DATA_SIMULATION_SETTINGS_H_ */

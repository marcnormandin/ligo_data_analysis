/*
 * simulation_file.h
 *
 *  Created on: May 8, 2017
 *      Author: marcnormandin
 */

#ifndef PROGRAMS_SIMULATE_DATA_SIMULATION_FILE_H_
#define PROGRAMS_SIMULATE_DATA_SIMULATION_FILE_H_

#include <stddef.h>

#include "simulation_settings.h"
#include "detector.h"
#include "strain.h"
#include "detector_network.h"
#include "inspiral_signal.h"

void append_index_to_prefix(char* buff, size_t buff_len, const char *prefix, size_t index);

void simulated_strain_file_create( const char *filename );

void simulated_strain_file_save_source( const char *output_filename, const source_t *source );

void simulated_strain_file_save_settings( const char *output_filename, const simulation_settings_t *ps);

void simulated_strain_file_save_detector( const char *output_filename, const detector_t* detector, size_t detector_num );

void simulated_strain_file_save_detector_network( const char *output_filename, const detector_network_t *dnet );

void simulated_strain_file_save_detector_signal( const char *output_filename, strain_t *signal, size_t detector_num );

void simulate( simulation_settings_t *ps, detector_network_t *net);


#endif /* PROGRAMS_SIMULATE_DATA_SIMULATION_FILE_H_ */

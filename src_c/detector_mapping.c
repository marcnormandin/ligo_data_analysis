/*
 * detector_mapping.c
 *
 *  Created on: Apr 10, 2017
 *      Author: marcnormandin
 */

#include <string.h>
#include <stdlib.h>
#include "detector_mapping.h"
#include "settings.h"

detector_mapping_t* detector_mapping_load( const char* filename ) {
	size_t i;

	detector_mapping_t *dm = (detector_mapping_t*) malloc(sizeof(detector_mapping_t));

	/* open the mapping file */
	settings_file_t *settings_file = settings_file_open( filename );
	if (settings_file == NULL) {
		printf("Error opening the detector mappings file (%s). Aborting.\n", filename );
		abort();
	}

	/* This assumes that the number of settings is the number of detectors */
	int num_detectors = settings_file_num_settings(settings_file);
	dm->num_detectors = num_detectors;

	dm->detector_names = (char**) malloc( num_detectors * sizeof(char*));
	dm->data_filenames = (char**) malloc( num_detectors * sizeof(char*));

	/* Store each detectors information */
	for (i = 0; i < num_detectors; i++) {
		/* Get the detector name */
		const char *detector_name = settings_file_get_key_by_index(settings_file, i);

		/* Get the input filename */
		const char *hdf_filename = settings_file_get_value(settings_file, detector_name);

		dm->detector_names[i] = (char*) malloc( (strlen(detector_name)+1) * sizeof(char) );
		strcpy( dm->detector_names[i], detector_name);

		dm->data_filenames[i] = (char*) malloc( (strlen(hdf_filename) +1) * sizeof(char) );
		strcpy( dm->data_filenames[i], hdf_filename);
	}

	settings_file_close(settings_file);

	return dm;
}

void detector_mapping_close( detector_mapping_t *dm) {
	size_t i;
	for (i = 0; i < dm->num_detectors; i++) {
		free(dm->detector_names[i]);
		free(dm->data_filenames[i]);
	}
	free(dm->data_filenames);
	free(dm->detector_names);
	free(dm);
	dm = NULL;
}


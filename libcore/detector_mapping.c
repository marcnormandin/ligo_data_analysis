#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "detector_mapping.h"
#include "settings_file.h"

detector_network_mapping_t* Detector_Network_Mapping_load( const char* filename ) {
	assert(filename != NULL);

	size_t i;

	detector_network_mapping_t *dm = (detector_network_mapping_t*) malloc(sizeof(detector_network_mapping_t));
	if (dm == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory for detector_mappting_t. Aborting.\n");
		abort();
	}

	/* open the mapping file */
	settings_file_t *settings_file = settings_file_open( filename );
	if (settings_file == NULL) {
		fprintf(stderr, "Error opening the detector mappings file (%s). Aborting.\n", filename );
		abort();
	}

	/* Careful: This assumes that the number of settings is the number of detectors */
	int num_detectors = settings_file_num_settings(settings_file);
	dm->num_detectors = num_detectors;

	dm->detector_names = (char**) malloc( num_detectors * sizeof(char*));
	if (dm->detector_names == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory for dm->detector_names. Aborting.\n");
		abort();
	}

	dm->data_filenames = (char**) malloc( num_detectors * sizeof(char*));
	if (dm->data_filenames == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory for dm->data_filenames. Aborting.\n");
		abort();
	}

	/* Store each detectors information */
	for (i = 0; i < num_detectors; i++) {
		/* Get the detector name */
		const char *detector_name = settings_file_get_key_by_index(settings_file, i);

		/* Get the input filename */
		const char *hdf_filename = settings_file_get_value(settings_file, detector_name);

		dm->detector_names[i] = (char*) malloc( (strlen(detector_name)+1) * sizeof(char) );
		if (dm->detector_names[i] == NULL) {
			fprintf(stderr, "Error. Unable to allocate memory for dm->detector_names[i]. Aborting.\n");
			abort();
		}
		strcpy( dm->detector_names[i], detector_name);

		dm->data_filenames[i] = (char*) malloc( (strlen(hdf_filename) +1) * sizeof(char) );
		if (dm->data_filenames[i] == NULL) {
			fprintf(stderr, "Error. Unable to allocate memory for dm->data_filnames[i]. Aborting.\n");
			abort();
		}
		strcpy( dm->data_filenames[i], hdf_filename);
	}

	settings_file_close(settings_file);

	return dm;
}

void Detector_Network_Mapping_close( detector_network_mapping_t *dm) {
	assert(dm != NULL);
	assert(dm->detector_names != NULL);
	assert(dm->data_filenames != NULL);

	size_t i;
	for (i = 0; i < dm->num_detectors; i++) {
		assert(dm->detector_names[i] != NULL);
		free(dm->detector_names[i]);
		dm->detector_names[i] = NULL;

		assert(dm->data_filenames[i] != NULL);
		free(dm->data_filenames[i]);
		dm->data_filenames[i] = NULL;
	}

	free(dm->detector_names);
	dm->detector_names = NULL;

	free(dm->data_filenames);
	dm->data_filenames = NULL;

	free(dm);
	dm = NULL;
}


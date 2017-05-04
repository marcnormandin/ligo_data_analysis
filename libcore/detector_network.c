#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_interp.h>

#include "detector.h"
#include "detector_mapping.h"
#include "detector_network.h"
#include "hdf5_file.h"
#include "sampling_system.h"

detector_network_t* Detector_Network_alloc(size_t num_detectors) {
	size_t i;

	detector_network_t *net = (detector_network_t*) malloc( sizeof(detector_network_t) );
	if (net == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory. Exiting.\n");
		exit(-1);
	}

	net->num_detectors = num_detectors;
	net->detector = (detector_t**) malloc(net->num_detectors * sizeof(detector_t*));
	if (net->detector == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory. Exiting.\n");
		exit(-1);
	}

	for (i = 0; i < net->num_detectors; i++) {
		net->detector[i] = Detector_alloc();
	}

	return net;
}

void Detector_Network_free(detector_network_t* net) {
	assert(net != NULL);

	size_t i;
	for (i = 0; i < net->num_detectors; i++) {
		Detector_free(net->detector[i]);
	}
	free(net->detector);
	net->detector = NULL;
}

void Detector_Network_print(detector_network_t* net) {
	assert(net != NULL);

	size_t i;

	printf("DETECTOR NEWTORK: %d detectors\n", net->num_detectors);
	for (i = 0; i < net->num_detectors; i++) {
		printf("%s ", net->detector[i]->name);
		printf("\n");
	}
}

detector_network_t* Detector_Network_load( const char* detector_mapping_file,
		size_t num_time_samples, double sampling_frequency, double f_low, double f_high ) {
	assert(detector_mapping_file != NULL);

	size_t i, j;

	detector_network_mapping_t *dmap = Detector_Network_Mapping_load( detector_mapping_file );

	/* sampling frequency */
	//double fs = hdf5_get_sampling_frequency( dmap->data_filenames[0] );
	//printf("The sampling frequency is %f.\n", fs);

	//int num_time_samples = hdf5_get_num_time_samples( dmap->data_filenames[0] );
	//printf("Simulated data will contain %d time samples.\n", num_time_samples);

	fprintf(stderr, "Allocating a (%lu) detector network... ", dmap->num_detectors);
	detector_network_t* net = Detector_Network_alloc ( dmap->num_detectors );
	fprintf(stderr, "done.\n");

	for (i = 0; i < net->num_detectors; i++) {
		/* Load the PSD from file */
		fprintf(stderr, "Loading the PSD for detector %lu from file (%s).\n", i, dmap->data_filenames[i]);
		psd_t *psd_unprocessed = PSD_load( dmap->data_filenames[i] );

		psd_t *psd = PSD_make_suitable_for_network_analysis(psd_unprocessed, num_time_samples, sampling_frequency, f_low, f_high);

		detector_t *det = net->detector[i];
		Detector_init_name( dmap->detector_names[i], psd, det);
	}

	printf("GW Detector network created: ");
	for (i = 0; i < net->num_detectors; i++) {
		printf("%s ", net->detector[i]->name);
	}
	printf("\n");

	Detector_Network_Mapping_close(dmap);

	return net;
}


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

	net->num_detectors = num_detectors;
	net->detector = (detector_t**) malloc(net->num_detectors * sizeof(detector_t*));

	for (i = 0; i < net->num_detectors; i++) {
		net->detector[i] = Detector_alloc();
	}

	return net;
}

void Detector_Network_free(detector_network_t* net) {
	size_t i;
	for (i = 0; i < net->num_detectors; i++) {
		Detector_free(net->detector[i]);
	}
	free(net->detector);
}

void Detector_Network_print(detector_network_t* net) {
	size_t i;

	printf("DETECTOR NEWTORK: %d detectors\n", net->num_detectors);
	for (i = 0; i < net->num_detectors; i++) {
		Detector_print(net->detector[i]);
		printf("\n");
	}
}

detector_network_t* Detector_Network_load( const char* detector_mapping_file ) {
	size_t i, j;

	detector_network_mapping_t *dmap = Detector_Network_Mapping_load( detector_mapping_file );

	/* sampling frequency */
	double fs = hdf5_get_sampling_frequency( dmap->data_filenames[0] );
	printf("The sampling frequency is %f.\n", fs);

	int num_time_samples = hdf5_get_num_time_samples( dmap->data_filenames[0] );
	printf("Simulated data will contain %d time samples.\n", num_time_samples);

	size_t half_size = SS_half_size(num_time_samples);
	printf("The half size is: %lu\n", half_size);

	fprintf(stderr, "Allocating a (%lu) detector network... ", dmap->num_detectors);
	detector_network_t* net = Detector_Network_alloc ( dmap->num_detectors );
	fprintf(stderr, "done.\n");

	for (i = 0; i < net->num_detectors; i++) {
		/* Load the PSD from file */
		fprintf(stderr, "Loading the PSD for detector %lu from file (%s).\n", i, dmap->data_filenames[i]);
		psd_t *psd_unprocessed = PSD_load( dmap->data_filenames[i] );

		/* Set the frequencies we want to use. */
		psd_t *psd = PSD_malloc( half_size );
		psd->type = PSD_ONE_SIDED;
		SS_frequency_array(fs, num_time_samples, psd->len, psd->f);

		/* Interpolate the PSD to the desired frequencies. */
		gsl_interp* interp = gsl_interp_alloc(gsl_interp_linear, psd_unprocessed->len);
		gsl_interp_accel* acc = gsl_interp_accel_alloc();
		gsl_interp_init(interp, psd_unprocessed->f, psd_unprocessed->psd, psd_unprocessed->len);
		for (j = 0; j < psd->len; j++) {
			double asd_interpolated = gsl_interp_eval(interp, psd_unprocessed->f, psd_unprocessed->psd, psd->f[j], acc);
			psd->psd[j] = asd_interpolated;
		}
		gsl_interp_accel_free(acc);
		gsl_interp_free(interp);

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


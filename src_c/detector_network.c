/*
 * detector_network.c
 *
 *  Created on: Mar 2, 2017
 *      Author: marcnormandin
 */

#include <stdlib.h>
#include <stdio.h>

#include "antenna_patterns.h"
#include "detector.h"
#include "detector_network.h"
#include "sky.h"

void Alloc_Detector_Network(int num, detector_network_t* net) {
	net->num_detectors = num;
	net->detector = (detector_t*) malloc(
			net->num_detectors * sizeof(detector_t));
}

void Free_Detector_Network(detector_network_t* net) {
	free(net->detector);
}

void Init_Detector_Network(detector_network_t* net) {
	Alloc_Detector_Network(4, net);

	strncpy(net->detector[0].id, "H1", 2);
	strncpy(net->detector[1].id, "L1", 2);
	strncpy(net->detector[2].id, "V1", 2);
	strncpy(net->detector[3].id, "K1", 2);
}

void Compute_Detector_Network_Antenna_Patterns(
		sky_t* sky,
		double polarization_angle,
		detector_network_t* net)
{
	for (size_t i = 0; i < net->num_detectors; i++) {
		detector_t* det = &net->detector[i];

		// Check return value for errors
		antenna_patterns( det->id, sky, polarization_angle, &det->ant);

		// Time delay
		time_delay(det->id, sky, &det->timedelay);
	}
}

void Print_Detector_Network(detector_network_t* net) {
	printf("DETECTOR NEWTORK: %d detectors\n", net->num_detectors);
	for (size_t i = 0; i < net->num_detectors; i++) {
		Print_Detector(&net->detector[i]);
		printf("\n");
	}
}



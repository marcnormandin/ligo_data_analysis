/*
 * detector_network.c
 *
 *  Created on: Mar 2, 2017
 *      Author: marcnormandin
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "antenna_patterns.h"
#include "detector.h"
#include "detector_network.h"
#include "sky.h"
#include "time_delay.h"

void Alloc_Detector_Network(int num, detector_network_t* net) {
	size_t i;

	net->num_detectors = num;
	net->detector = (detector_t**) malloc(net->num_detectors * sizeof(detector_t*));

	for (i = 0; i < net->num_detectors; i++) {
		net->detector[i] = Detector_alloc();
	}
}

void Free_Detector_Network(detector_network_t* net) {
	size_t i;
	for (i = 0; i < net->num_detectors; i++) {
		Detector_free(net->detector[i]);
	}
	free(net->detector);
}

void Init_Detector_Network(detector_network_t* net) {
	Alloc_Detector_Network(4, net);

	Detector_init(H1, net->detector[0]);
	Detector_init(L1, net->detector[1]);
	Detector_init(V1, net->detector[2]);
	Detector_init(K1, net->detector[3]);
}

void Print_Detector_Network(detector_network_t* net) {
	size_t i;

	printf("DETECTOR NEWTORK: %d detectors\n", net->num_detectors);
	for (i = 0; i < net->num_detectors; i++) {
		Print_Detector(net->detector[i]);
		printf("\n");
	}
}



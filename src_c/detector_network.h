/*
 * detector_network.h
 *
 *  Created on: Mar 2, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_DETECTOR_NETWORK_H_
#define SRC_C_DETECTOR_NETWORK_H_

#include "detector.h"
#include "sky.h"

/* The detector network used in the study. */
typedef struct detector_network_s {
	unsigned int num_detectors;
	detector_t** detector;

} detector_network_t;

/* Detector functions */
void Alloc_Detector_Network(int num, detector_network_t* net);
void Free_Detector_Network(detector_network_t* net);
void Init_Detector_Network(detector_network_t* net);
void Print_Detector_Network(detector_network_t* net);

#endif /* SRC_C_DETECTOR_NETWORK_H_ */

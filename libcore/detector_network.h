/*
 * detector_network.h
 *
 *  Created on: Mar 2, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_DETECTOR_NETWORK_H_
#define SRC_C_DETECTOR_NETWORK_H_

#include "../libcore/detector.h"
#include "../libcore/sky.h"

/* The detector network used in the study. */
typedef struct detector_network_s {
	unsigned int num_detectors;
	detector_t** detector;

} detector_network_t;

/* Detector functions */
detector_network_t* Alloc_Detector_Network(size_t num_detectors);

void Free_Detector_Network(detector_network_t* net);

/*
void Init_Detector_Network(detector_network_t* net);
*/

void Print_Detector_Network(detector_network_t* net);

detector_network_t* detector_network_load( const char* detector_mapping_file );

#endif /* SRC_C_DETECTOR_NETWORK_H_ */

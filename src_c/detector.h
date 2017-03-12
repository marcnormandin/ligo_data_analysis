/*
 * detector.h
 *
 *  Created on: Mar 2, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_DETECTOR_H_
#define SRC_C_DETECTOR_H_

#include "antenna_patterns.h"

/* For a given source, this records the values for a given detector. */
typedef struct {
	char id[3];

	antenna_patterns_t ant;

	double timedelay;
} detector_t;

void Print_Detector(detector_t* det);


#endif /* SRC_C_DETECTOR_H_ */

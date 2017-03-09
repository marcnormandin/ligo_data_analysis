/*
 * detector.c
 *
 *  Created on: Mar 2, 2017
 *      Author: marcnormandin
 */

#include "detector.h"

void Print_Detector(detector_t* det) {
	printf("DETECTOR:\n");
	printf("id: %s\n", det->id);
	printf("u: %f\n", det->ant.u);
	printf("v: %f\n", det->ant.v);
	printf("f_plus: %f\n", det->ant.f_plus);
	printf("f_cross: %f\n", det->ant.f_cross);
}

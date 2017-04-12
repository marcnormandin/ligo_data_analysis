/*
 * detector.h
 *
 *  Created on: Mar 2, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_DETECTOR_H_
#define SRC_C_DETECTOR_H_

#include <stddef.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#include "../libcore/spectral_density.h"

typedef enum {
	L1,
	H1,
	H2,
	V1,
	K1,
	G1,
	T1
} DETECTOR_ID;

const char* detector_id_to_name(DETECTOR_ID id);

DETECTOR_ID detector_name_to_id(const char* name);

/* For a given source, this records the values for a given detector. */
typedef struct detector_s {
	DETECTOR_ID id;

	char name[255];

	/* position vector for location on the Earth */
	gsl_vector *location;

	/* arm direction vectors */
	gsl_vector *arm_x;
	gsl_vector *arm_y;

	/* detector tensor */
	gsl_matrix *detector_tensor;

	/* detector PSD */
	psd_t *psd;

	/* detector ASD */
	asd_t *asd;

} detector_t;

void Print_Detector(detector_t* det);

detector_t* Detector_alloc();
void Detector_free(detector_t *d);

void Detector_init(DETECTOR_ID name, psd_t *psd, detector_t *d);
void Detector_init_name( char *name, psd_t *psd, detector_t *d);

void Detector_init_L1(asd_t *asd, psd_t *psd, detector_t *d);
void Detector_init_H1(asd_t *asd, psd_t *psd, detector_t *d);
void Detector_init_H2(asd_t *asd, psd_t *psd, detector_t *d);
void Detector_init_V1(asd_t *asd, psd_t *psd, detector_t *d);
void Detector_init_G1(asd_t *asd, psd_t *psd, detector_t *d);
void Detector_init_K1(asd_t *asd, psd_t *psd, detector_t *d);
void Detector_init_T1(asd_t *asd, psd_t *psd, detector_t *d);

void Detector_tensor(detector_t *d);


#endif /* SRC_C_DETECTOR_H_ */

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

typedef enum {
	L1,
	H1,
	H2,
	V1,
	K1,
	G1,
	T1
} DETECTOR_NAME;

/* For a given source, this records the values for a given detector. */
typedef struct {
	DETECTOR_NAME name;

	/* position vector for location on the Earth */
	gsl_vector *location;

	/* arm direction vectors */
	gsl_vector *arm_x;
	gsl_vector *arm_y;

	/* detector tensor */
	gsl_matrix *detector_tensor;

} detector_t;

typedef struct detector_data_s {
	DETECTOR_NAME name;
	double sampling_frequency;

	size_t num_samples;
	double *sample_times;
	double *recordings;
	gsl_complex *fft_recordings;

} detector_data_t;


void Print_Detector(detector_t* det);

detector_t* Detector_alloc();
void Detector_free(detector_t *d);

void Detector_init(DETECTOR_NAME name, detector_t *d);
void Detector_init_name( char *name, detector_t *d);

void Detector_init_L1(detector_t *d);
void Detector_init_H1(detector_t *d);
void Detector_init_H2(detector_t *d);
void Detector_init_V1(detector_t *d);
void Detector_init_G1(detector_t *d);
void Detector_init_K1(detector_t *d);
void Detector_init_T1(detector_t *d);
void Detector_tensor(detector_t *d);


#endif /* SRC_C_DETECTOR_H_ */

/*
 * detector.c
 *
 *  Created on: Mar 2, 2017
 *      Author: marcnormandin
 */

#include "detector.h"

#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>


void Print_Detector(detector_t* det) {
	printf("DETECTOR:\n");
}

detector_t* Detector_alloc() {
	detector_t *d = (detector_t*) malloc( sizeof(detector_t) );
	d->location = gsl_vector_alloc(3);
	d->arm_x = gsl_vector_alloc(3);
	d->arm_y = gsl_vector_alloc(3);
	d->detector_tensor = gsl_matrix_alloc(3,3);
	return d;
}

void Detector_free(detector_t *d) {
	gsl_vector_free(d->location);
	gsl_vector_free(d->arm_x);
	gsl_vector_free(d->arm_y);
	gsl_matrix_free(d->detector_tensor);
	free(d);
	d = NULL;
}

void Detector_init(DETECTOR_NAME name, detector_t *d) {
	switch (name) {
	case L1: Detector_init_L1(d); break;
	case H1: Detector_init_H1(d); break;
	case H2: Detector_init_H2(d); break;
	case V1: Detector_init_V1(d); break;
	case G1: Detector_init_G1(d); break;
	case K1: Detector_init_K1(d); break;
	case T1: Detector_init_T1(d); break;
	default:
		printf("ERROR: Invalid detector name. Can not initialize.\n");
		abort();
	}
}

/* Local helper function */
static void init_gsl_vector(double x, double y, double z, gsl_vector *v) {
	gsl_vector_set(v, 0, x);
	gsl_vector_set(v, 1, y);
	gsl_vector_set(v, 2, z);
}

void Detector_init_L1(detector_t *d) {
	/* Generate coordinates for LLO: */
	init_gsl_vector(-0.954574615, -0.141579994, -0.262187738, d->arm_x);
	init_gsl_vector(0.297740169, -0.487910627, -0.820544948, d->arm_y);
	/* X-arm unit vector, L1-|cordnt.X| = 3995.15;
	   Y-arm unit vector, L1-|cordnt.Y| = 3995.15; */

	/* WGS-84 coordinates of "L1" beam splitter */
	init_gsl_vector(-7.427604192e+4, -5.496283721e+6, 3.224257016e+6, d->location);

	Detector_tensor(d);
}

void Detector_init_H1(detector_t *d) {
	/* All LIGO coordinates are from LIGO-P000006-D-E: */

	/* Generate coordinates for LHO:*/

	/* X-arm unit vector relative to earth centered frame. */
	init_gsl_vector(-0.223891216, 0.799830697, 0.556905359, d->arm_x);
	/* Y-arm unit vector relative to earth centered frame. */
	init_gsl_vector(-0.913978490, 0.0260953206, -0.404922650, d->arm_y);

	init_gsl_vector(2.161414928e+6, -3.834695183e+6, 4.600350224e+6, d->location); /* In meters */

	Detector_tensor(d);
}

void Detector_init_H2(detector_t *d) {
	/* Generate coordinates for LHO:*/

	/* X-arm unit vector relative to earth centered frame. */
	init_gsl_vector(-0.223891216, 0.799830697, 0.556905359, d->arm_x);
	/* Y-arm unit vector relative to earth centered frame. */
	init_gsl_vector(-0.913978490, 0.0260953206, -0.404922650, d->arm_y);

	init_gsl_vector(2.161414928e+6, -3.834695183e+6, 4.600350224e+6, d->location); /* In meters */

	Detector_tensor(d);
}

void Detector_init_V1(detector_t *d) {
	init_gsl_vector(-0.70045821479, 0.20848948619, 0.68256166277, d->arm_x);
	init_gsl_vector(-0.05379255368, -0.96908180549, 0.24080451708, d->arm_y);

	/* WGS-84 coordinates of "V1" beam splitter */
	init_gsl_vector(4.5463741e+6, 8.429897e+5, 4.378577e+6, d->location);

	Detector_tensor(d);
}

void Detector_init_G1(detector_t *d) {
	init_gsl_vector(-0.445184239, 0.866534205, 0.225675575, d->arm_x);
	init_gsl_vector(-0.626000687, -0.552167273, 0.550667271, d->arm_y);

	/* WGS-84 coordinates of "G1" beam splitter */
	init_gsl_vector(3.8563112e+6,  6.665978e+5, 5.0196406e+6, d->location);

	Detector_tensor(d);
}

void Detector_init_K1(detector_t *d) {
	/* Source : www.ligo.org/scientists/GW100916/detectors.txt */
	init_gsl_vector(-0.4300, -0.8363, 0.3400, d->arm_x);
	init_gsl_vector(0.6821, -0.0542, 0.7292, d->arm_y);

	/* WGS-84 coordinates of "K1" beam splitter */
	init_gsl_vector(-3.776899062e+6,  3.483900163e+6, 3766657.585, d->location);

	Detector_tensor(d);
}

void Detector_init_T1(detector_t *d) {
	init_gsl_vector(0.648969405, 0.760814505, 0, d->arm_x);
	init_gsl_vector(-0.443713769, 0.378484715, -0.812322234, d->arm_y);

	/* WGS-84 coordinates of "T1" beam splitter */
	init_gsl_vector(-3.946409e6, 3.366259e6, 3.6991507e6, d->location);

	Detector_tensor(d);
}

void Detector_tensor(detector_t *d) {
	/* Detector Tensor
	  D = 0.5*(Xhat"*Xhat - Yhat"*Yhat) where Xhat is the unit vector for the
	  "X" arm and Yhat is the unit vector for the "Y" arm
	 */
	gsl_matrix_set_zero(d->detector_tensor);
	gsl_blas_dger(0.5, d->arm_x, d->arm_x, d->detector_tensor);
	gsl_blas_dger(-0.5, d->arm_y, d->arm_y, d->detector_tensor);
}

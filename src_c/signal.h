/*
 * signal.h
 *
 *  Created on: Mar 4, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_SIGNAL_H_
#define SRC_C_SIGNAL_H_

#include <stddef.h>
#include <gsl/gsl_complex.h>

typedef struct signal_s {
	double* tempval;
	gsl_complex* whitened_sf;
	gsl_complex* h_0;
	gsl_complex* h_90;
	gsl_complex* whitened_signal;
	gsl_complex* whitened_data;

} signal_t;

signal_t* Signal_malloc(size_t size);

void Signal_free(signal_t *s);

#endif /* SRC_C_SIGNAL_H_ */

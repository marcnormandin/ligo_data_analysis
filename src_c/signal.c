/*
 * signal.c
 *
 *  Created on: Mar 9, 2017
 *      Author: marcnormandin
 */

#include <stdlib.h>

#include <gsl/gsl_complex.h>

#include "signal.h"

signal_t* Signal_malloc(size_t size) {
	signal_t *s = (signal_t*) malloc( sizeof(signal_t) );
	s->whitened_sf = (gsl_complex*) malloc( size * sizeof(gsl_complex) );
	s->h_0 = (gsl_complex*) malloc( size * sizeof(gsl_complex) );
	s->h_90 = (gsl_complex*) malloc( size * sizeof(gsl_complex) );
	s->whitened_signal = (gsl_complex*) malloc( size * sizeof(gsl_complex) );
	s->whitened_data = (gsl_complex*) malloc( size * sizeof(gsl_complex) );

	return s;
}

void Signal_free(signal_t *s) {
	free(s->whitened_sf);
	free(s->h_0);
	free(s->h_90);
	free(s->whitened_signal);
	free(s->whitened_data);
	free(s);
}

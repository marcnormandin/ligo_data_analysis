/*
 * datagen.h
 *
 *  Created on: Feb 22, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_DATAGEN_H_
#define SRC_C_DATAGEN_H_

#include <gsl/gsl_rng.h>

#include "inspiral_source.h"
#include "signal.h"
#include "strain.h"
#include "detector_network.h"

/* This is the main routine that simulates the data
   The memory must be freed. */
signal_t** simulate_inspiral(gsl_rng *rng, double f_low, double f_high, detector_network_t *net, strain_t *strain, source_t *source);

signal_t** simulate_inspiral_unscaled(double f_low, double f_high, detector_network_t *net, strain_t *strain, source_t *source);

#endif /* SRC_C_DATAGEN_H_ */

/*
 * simulate_noise.h
 *
 *  Created on: Mar 17, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_SIMULATE_NOISE_H_
#define SRC_C_SIMULATE_NOISE_H_

#include <gsl/gsl_complex.h>
#include <gsl/gsl_rng.h>

/* Simulate white gaussian noise in the complex domain */
gsl_complex SN_whitenoise_frequency_domain (gsl_rng *rng);

#endif /* SRC_C_SIMULATE_NOISE_H_ */

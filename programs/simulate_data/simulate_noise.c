/*
 * simulate_noise.c
 *
 *  Created on: Mar 17, 2017
 *      Author: marcnormandin
 */

#include "simulate_noise.h"

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

gsl_complex SN_whitenoise_frequency_domain (gsl_rng *rng) {
	double noise_f_real = 0.5 * gsl_ran_gaussian(rng, 1.0);
	double noise_f_imag = 0.5 * gsl_ran_gaussian(rng, 1.0);

	return gsl_complex_rect(noise_f_real, noise_f_imag);
}

/*
 * random.c
 *
 *  Created on: Mar 22, 2017
 *      Author: marcnormandin
 */

#include "../libcore/random.h"

#include <time.h>
#include <gsl/gsl_math.h> // M_PI
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>

gsl_rng* random_alloc(gslseed_t seed) {
	const gsl_rng_type *rng_type;
	gsl_rng *rng;
	gsl_rng_env_setup();
	rng_type = gsl_rng_default;
	rng = gsl_rng_alloc(rng_type);
	if (seed == 0) {
		gsl_rng_set(rng, time(0));
	} else {
		gsl_rng_set(rng, seed);
	}

	return rng;
}

void random_free(gsl_rng* rng) {
	gsl_rng_free(rng);
}

gslseed_t random_seed (gsl_rng* rng) {
	/* I wrote the following after looking into the GSL code itself. */
	const gslseed_t offset = rng->type->min;
	const gslseed_t range = rng->type->max - offset;
	gslseed_t rseed = 1 + gsl_rng_uniform_int(rng, range-2);
	return rseed;
}


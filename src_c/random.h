/*
 * random.h
 *
 *  Created on: Mar 22, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_RANDOM_H_
#define SRC_C_RANDOM_H_

#include <gsl/gsl_rng.h>

typedef unsigned long int gslseed_t;

gsl_rng* random_alloc(gslseed_t seed);

void random_free(gsl_rng* rng);

gslseed_t random_seed (gsl_rng* rng);



#endif /* SRC_C_RANDOM_H_ */

#ifndef SRC_C_RANDOM_H_
#define SRC_C_RANDOM_H_

#include <gsl/gsl_rng.h>

#if defined (__cplusplus)
extern "C" {
#endif

typedef unsigned long int gslseed_t;

gsl_rng* random_alloc(gslseed_t seed);

void random_free(gsl_rng* rng);

gslseed_t random_seed (gsl_rng* rng);

#if defined (__cplusplus)
}
#endif

#endif /* SRC_C_RANDOM_H_ */

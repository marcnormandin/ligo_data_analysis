/*
 * ptapso_estimate.h
 *
 *  Created on: Mar 10, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_PTAPSO_ESTIMATE_H_
#define SRC_C_PTAPSO_ESTIMATE_H_

#include <stddef.h>
#include "random.h"
#include "ptapso_func.h"

typedef struct pso_result_s {
	double ra;
	double dec;
	double snr;

} pso_result_t;

int ptapso_estimate(ptapso_fun_params_t *splParams, gslseed_t seed, size_t max_steps, pso_result_t* result);

#endif /* SRC_C_PTAPSO_ESTIMATE_H_ */

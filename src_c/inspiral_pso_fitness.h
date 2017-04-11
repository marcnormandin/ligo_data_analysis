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
#include "inspiral_source.h"
#include "signal.h"
#include "spectral_density.h"
#include "detector_network.h"
#include "inspiral_network_statistic.h"

typedef struct pso_result_s {
	double ra;
	double dec;
	double chirp_t0;
	double chirp_t1_5;
	double snr;

} pso_result_t;

typedef struct ptapso_func_params {
	double f_low;
	double f_high;
	detector_network_t *network;
	inspiral_template_half_fft_t **signals;
	coherent_network_workspace_t *workspace;

} ptapso_fun_params_t;

double ptapso_func(gsl_vector *xVec, void  *inParamsPointer);

int ptapso_estimate_ra_dec_t0_t1_5(char *pso_settings_file, ptapso_fun_params_t *splParams, gslseed_t seed, pso_result_t* result);

#endif /* SRC_C_PTAPSO_ESTIMATE_H_ */

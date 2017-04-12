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
#include "strain.h"
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

typedef struct pso_fitness_function_parameters_s {
	double f_low;
	double f_high;
	detector_network_t *network;
	network_strain_half_fft_t *network_strain;
	coherent_network_workspace_t *workspace;

} pso_fitness_function_parameters_t;

pso_fitness_function_parameters_t* pso_fitness_function_parameters_alloc(
		double f_low, double f_high, detector_network_t* network, network_strain_half_fft_t *network_strain);

void pso_fitness_function_parameters_free(pso_fitness_function_parameters_t *params);

double pso_fitness_function(gsl_vector *xVec, void  *inParamsPointer);

int pso_estimate_parameters(char *pso_settings_file, pso_fitness_function_parameters_t *splParams, gslseed_t seed, pso_result_t* result);

#endif /* SRC_C_PTAPSO_ESTIMATE_H_ */

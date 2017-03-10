/*
 * ptapso_fitnessfunc.h
 *
 *  Created on: Mar 10, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_PTAPSO_FUNC_H_
#define SRC_C_PTAPSO_FUNC_H_

#include <gsl/gsl_vector.h>

#include "detector_network.h"
#include "signal.h"
#include "source.h"
#include "strain.h"
#include "network_analysis.h"

typedef struct ptapso_func_params {
	double f_low;
	double f_high;
	strain_t *strain;
	detector_network_t *network;
	source_t *source;
	signal_t **signals;
	coherent_network_workspace_t *workspace;

} ptapso_fun_params_t;

double ptapso_func(gsl_vector *xVec, void  *inParamsPointer);

#endif /* SRC_C_PTAPSO_FUNC_H_ */

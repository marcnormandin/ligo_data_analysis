#ifndef SRC_C_PTAPSO_ESTIMATE_H_
#define SRC_C_PTAPSO_ESTIMATE_H_

#include <stddef.h>
#include "random.h"
#include "strain.h"
#include "spectral_density.h"
#include "detector_network.h"
#include "inspiral_network_statistic.h"

#include "parallel.h"

#if defined (__cplusplus)
extern "C" {
#endif

/*
typedef enum {
	PSO_LBEST = 0,
	PSO_GBEST
} PSO_VERSION;
*/

typedef struct pso_result_s {
	double ra;
	double dec;
	double chirp_t0;
	double chirp_t1_5;
	double snr;

	/* diagnostics */
	size_t total_iterations;
	size_t total_func_evals;
	double computation_time_secs;

} pso_result_t;

typedef struct pso_fitness_function_parameters_s {
	double f_low;
	double f_high;
	detector_network_t *network;
	network_strain_half_fft_t *network_strain;
	coherent_network_workspace_t **workspace;
} pso_fitness_function_parameters_t;

pso_fitness_function_parameters_t* pso_fitness_function_parameters_alloc(
		double f_low, double f_high, detector_network_t* network, network_strain_half_fft_t *network_strain);

void pso_fitness_function_parameters_free(pso_fitness_function_parameters_t *params);

double pso_fitness_function(gsl_vector *xVec, void  *inParamsPointer);

int pso_estimate_parameters(char *pso_settings_file, pso_fitness_function_parameters_t *splParams, gslseed_t seed, pso_result_t* result);

void CN_template_chirp_time(double f_low, double chirp_time0, double chirp_time1_5, inspiral_chirp_time_t *ct);

#if defined (__cplusplus)
}
#endif

#endif /* SRC_C_PTAPSO_ESTIMATE_H_ */

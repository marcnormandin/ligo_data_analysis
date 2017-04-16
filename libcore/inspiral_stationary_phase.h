#ifndef SRC_C_STATIONARY_PHASE_H_
#define SRC_C_STATIONARY_PHASE_H_

#include <stddef.h>

#include <gsl/gsl_complex.h>

#include "inspiral_chirp_time.h"
#include "spectral_density.h"

typedef struct stationary_phase_lookup_s {
	double f_low, f_high;
	size_t f_low_index, f_high_index;

	size_t len;
	double *g_coeff;
	double *chirp_tc_coeff;
	double *coalesce_phase_coeff;
	double *chirp_time_0_coeff;
	double *chirp_time_1_coeff;
	double *chirp_time1_5_coeff;
	double *chirp_time2_coeff;

} stationary_phase_workspace_t;

typedef struct stationary_phase_s {
	size_t 			len;
	gsl_complex		*spa_0;
	gsl_complex		*spa_90;

} stationary_phase_t;




/* Stationary Phase Workspace functions */
stationary_phase_workspace_t* SP_workspace_alloc(double f_low, double f_high, size_t len_f_array, double *f_array);

void SP_workspace_free( stationary_phase_workspace_t *lookup);




/* Stationary Phase functions */
stationary_phase_t* SP_malloc(size_t size);

void SP_free(stationary_phase_t *sp);

double SP_normalization_factor(asd_t *asd, stationary_phase_workspace_t *lookup);

void SP_compute(
		double detector_time_delay, double detector_normalization_factor,
		double inspiral_coalesce_phase, inspiral_chirp_time_t *chirp,
		stationary_phase_workspace_t *lookup,
		stationary_phase_t *out_sp);

void SP_save(char *filename, asd_t *asd, stationary_phase_t *sp);



#endif /* SRC_C_STATIONARY_PHASE_H_ */

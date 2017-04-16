#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>

#include "inspiral_stationary_phase.h"


int find_index_low(double f_low, size_t len, double *f_array) {
	size_t i;
	for (i = 0; i < len; i++) {
		if (f_array[i] >= f_low) {
			return i;
		}
	}

	/* error. we didn't find the index. */
	return -1;
}

int find_index_high(double f_high, size_t len, double *f_array) {
	size_t i;
	for (i = 0; i < len; i++) {
		if (f_array[i] > f_high) {
			return i-1;
		} else if (f_array[i] == f_high) {
			return i;
		}
	}

	/* The index wasn't found. */
	return -1;
}

/* This is called by SP_workspace_alloc and shouldn't be called otherwise. */
void SP_workspace_init(size_t len_f_array, double *f_array, stationary_phase_workspace_t *lookup) {
	assert(f_array != NULL);
	assert(lookup != NULL);

	size_t i, j;

	for (i = lookup->f_low_index, j = 0; i <= lookup->f_high_index; i++, j++) {
		double f = f_array[i];
		double f_fac = f / lookup->f_low;

		lookup->g_coeff[j] = pow(f, -7.0 / 6.0);
		lookup->chirp_tc_coeff[j] = 2.0 * M_PI * f;
		lookup->constant_coeff[j] = -M_PI / 4.0;
		lookup->chirp_time_0_coeff[j] = 2.0 * M_PI * lookup->f_low * 3.0 * pow(f_fac, -5.0 / 3.0) / 5.0;
		lookup->chirp_time_1_coeff[j] =  2.0 * M_PI * lookup->f_low * pow(f_fac, -5.0 / 3.0);
		lookup->chirp_time1_5_coeff[j] = -2.0 * M_PI * lookup->f_low * 3.0 * pow(f_fac,-2.0 / 3.0) / 2.0;
		lookup->chirp_time2_coeff[j] = 2.0 * M_PI * lookup->f_low * 3.0 * pow(f_fac, -1.0 / 3.0);
	}
}

stationary_phase_workspace_t* SP_workspace_alloc(double f_low, double f_high, size_t len_f_array, double *f_array) {
	assert(f_array != NULL);

	size_t i;
	int index_found;

	stationary_phase_workspace_t *lookup = (stationary_phase_workspace_t*) malloc( sizeof(stationary_phase_workspace_t) );
	if (lookup == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory in SP_workspace_alloc(). Exiting.\n");
		exit(-1);
	}

	lookup->f_low = f_low;
	lookup->f_high = f_high;

	if (lookup->f_low >= lookup->f_high) {
		fprintf(stderr, "Error. f_low (%f) >= f_high(%f). Exiting.\n", lookup->f_low, lookup->f_high);
		exit(-1);
	}

	index_found = find_index_low(lookup->f_low, len_f_array, f_array);
	if (index_found != -1) {
		lookup->f_low_index = index_found;
	} else {
		fprintf(stderr, "Error. Unable to find the index for f_low. Aborting.\n");
		exit(-1);
	}

	index_found = find_index_high(lookup->f_high, len_f_array, f_array);
	if (index_found != -1) {
		lookup->f_high_index = index_found;
	} else {
		fprintf(stderr, "Error. Unable to find the index for f_high. Aborting.\n");
		exit(-1);
	}

	/* Make sure the indices make sense. */
	if (lookup->f_low_index >= lookup->f_high_index) {
		fprintf(stderr, "Error. f_low_index (%lu) >= f_high_index (%lu). Exiting.\n",
				lookup->f_low_index, lookup->f_high_index);
		exit(-1);
	}

	lookup->len = lookup->f_high_index - lookup->f_low_index + 1;

	lookup->g_coeff = (double*) malloc (lookup->len * sizeof(double));
	if (lookup->g_coeff == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory in SP_workspace_alloc(). Exiting.\n");
		exit(-1);
	}

	lookup->chirp_tc_coeff = (double*) malloc (lookup->len * sizeof(double));
	if (lookup->chirp_tc_coeff == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory in SP_workspace_alloc(). Exiting.\n");
		exit(-1);
	}

	lookup->constant_coeff = (double*) malloc (lookup->len * sizeof(double));
	if (lookup->constant_coeff == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory in SP_workspace_alloc(). Exiting.\n");
		exit(-1);
	}

	lookup->chirp_time_0_coeff = (double*) malloc (lookup->len * sizeof(double));
	if (lookup->chirp_time_0_coeff == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory in SP_workspace_alloc(). Exiting.\n");
		exit(-1);
	}

	lookup->chirp_time_1_coeff = (double*) malloc (lookup->len * sizeof(double));
	if (lookup->chirp_time_1_coeff == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory in SP_workspace_alloc(). Exiting.\n");
		exit(-1);
	}

	lookup->chirp_time1_5_coeff = (double*) malloc (lookup->len * sizeof(double));
	if (lookup->chirp_time1_5_coeff == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory in SP_workspace_alloc(). Exiting.\n");
		exit(-1);
	}

	lookup->chirp_time2_coeff = (double*) malloc (lookup->len * sizeof(double));
	if (lookup->chirp_time2_coeff == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory in SP_workspace_alloc(). Exiting.\n");
		exit(-1);
	}

	/* Lookup is now set up and can be initialized with the coefficients. */
	SP_workspace_init(len_f_array, f_array, lookup);

	return lookup;
}

void SP_workspace_free( stationary_phase_workspace_t *lookup) {
	assert(lookup != NULL);

	assert(lookup->g_coeff != NULL);
	free(lookup->g_coeff);
	lookup->g_coeff = NULL;

	assert(lookup->chirp_tc_coeff != NULL);
	free(lookup->chirp_tc_coeff);
	lookup->chirp_tc_coeff = NULL;

	assert(lookup->constant_coeff != NULL);
	free(lookup->constant_coeff);
	lookup->constant_coeff = NULL;

	assert(lookup->chirp_time_0_coeff != NULL);
	free(lookup->chirp_time_0_coeff);
	lookup->chirp_time_0_coeff = NULL;

	assert(lookup->chirp_time_1_coeff != NULL);
	free(lookup->chirp_time_1_coeff);
	lookup->chirp_time_1_coeff = NULL;

	assert(lookup->chirp_time1_5_coeff != NULL);
	free(lookup->chirp_time1_5_coeff);
	lookup->chirp_time1_5_coeff = NULL;

	assert(lookup->chirp_time2_coeff != NULL);
	free(lookup->chirp_time2_coeff);
	lookup->chirp_time2_coeff = NULL;

	free(lookup);
}


stationary_phase_t* SP_alloc(size_t num_half_frequencies) {
	size_t i;

	stationary_phase_t* sp = (stationary_phase_t*) malloc( sizeof(stationary_phase_t) );
	if (sp == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory in SP_malloc. Exiting.\n");
		exit(-1);
	}

	/* common length of the arrays. */
	sp->len = num_half_frequencies;

	sp->spa_0 = (gsl_complex*) malloc(sp->len * sizeof(gsl_complex));
	if (sp->spa_0 == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory in SP_malloc. Exiting.\n");
		exit(-1);
	}

	sp->spa_90 = (gsl_complex*) malloc(sp->len * sizeof(gsl_complex));
	if (sp->spa_90 == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory in SP_malloc. Exiting.\n");
		exit(-1);
	}

	for (i = 0; i < sp->len; i++) {
		sp->spa_0[i] = gsl_complex_rect(0.0, 0.0);
		sp->spa_90[i] = gsl_complex_rect(0.0, 0.0);
	}

	return sp;
}

void SP_free(stationary_phase_t* sp) {
	assert(sp != NULL);

	assert(sp->spa_0);
	free(sp->spa_0);
	sp->spa_0 = NULL;

	assert(sp->spa_90);
	free(sp->spa_90);
	sp->spa_90 = NULL;

	free(sp);
}

void SP_save(char* filename, asd_t* asd, stationary_phase_t* sp) {
	assert(filename != NULL);
	assert(asd != NULL);
	assert(sp != NULL);

	FILE* file;
	size_t i;

	file = fopen(filename, "w");
	if (file == NULL) {
		fprintf(stderr, "Error. Unable to open file (%s) for writing in SP_save. Exiting.\n", filename);
		exit(-1);
	}

	for (i = 0; i < sp->len; i++) {
		fprintf(file, "%e %e %e\n", asd->asd[i], GSL_REAL(sp->spa_0[i]), GSL_REAL(sp->spa_90[i]));
	}

	fclose(file);
}

/* whitening normalization factor for inner product */
double SP_normalization_factor(asd_t* asd, stationary_phase_workspace_t *lookup) {
	assert(asd != NULL);
	assert(lookup != NULL);

	size_t i;
	double sum;

	/* All values that are not used below are assigned a value of 1, so simply
	 * determine how many there would be, and use that instead of iterating and sum += 1.
	 *
	 * It should be checked if the 1's are even needed! But this is what the Matlab version does.
	 */
	sum = lookup->f_low_index + (asd->len - lookup->f_high_index + 1);

	for (i = lookup->f_low_index; i <= lookup->f_high_index; i++) {
		double f = asd->f[i];
		double s = asd->asd[i];

		sum  += pow(f, -7.0 / 3.0) / gsl_pow_2(s);
	}

	return sqrt(sum);
}

void SP_compute(
		double detector_time_delay, double detector_normalization_factor,
		double inspiral_coalesce_phase, inspiral_chirp_time_t *chirp,
		stationary_phase_workspace_t *lookup,
		stationary_phase_t *out_sp)
{
	assert(chirp != NULL);
	assert(lookup != NULL);
	assert(out_sp != NULL);

	size_t i;

	for (i = 0; i < lookup->len; i++) {
		double amp_2pn = (1.0 / detector_normalization_factor) * lookup->g_coeff[i];

		double phase_2pn =
				lookup->chirp_tc_coeff[i] * (chirp->tc - detector_time_delay)
				- 2.0 * inspiral_coalesce_phase
				+ lookup->constant_coeff[i]
				+ lookup->chirp_time_0_coeff[i] * chirp->chirp_time0
				+ lookup->chirp_time_1_coeff[i] * chirp->chirp_time1
				+ lookup->chirp_time1_5_coeff[i] * chirp->chirp_time1_5
				+ lookup->chirp_time2_coeff[i] * chirp->chirp_time2;

		gsl_complex exp_phase = gsl_complex_exp(gsl_complex_rect(0.0, -1.0*phase_2pn));
		out_sp->spa_0[lookup->f_low_index + i] = gsl_complex_mul_real(exp_phase, amp_2pn);
		out_sp->spa_90[lookup->f_low_index + i] = gsl_complex_mul_imag(out_sp->spa_0[lookup->f_low_index + i], -1.0);
	}
}

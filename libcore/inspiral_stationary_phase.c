#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>

#include "inspiral_stationary_phase.h"


stationary_phase_t* SP_malloc(size_t size) {
	size_t i;

	stationary_phase_t* sp = (stationary_phase_t*) malloc( sizeof(stationary_phase_t) );
	sp->len = size;
	sp->spa_0 = (gsl_complex*) malloc(size * sizeof(gsl_complex));
	sp->spa_90 = (gsl_complex*) malloc(size * sizeof(gsl_complex));

	for (i = 0; i < size; i++) {
		sp->spa_0[i] = gsl_complex_rect(0.0, 0.0);
		sp->spa_90[i] = gsl_complex_rect(0.0, 0.0);
	}

	return sp;
}

void SP_free(stationary_phase_t* sp) {
	free(sp->spa_0);
	free(sp->spa_90);
	free(sp);
}

void SP_save(char* filename, asd_t* asd, stationary_phase_t* sp) {
	FILE* file;
	size_t i;

	file = fopen(filename, "w");
	for (i = 0; i < sp->len; i++) {
		fprintf(file, "%e %e %e\n", asd->asd[i], GSL_REAL(sp->spa_0[i]), GSL_REAL(sp->spa_90[i]));
	}
	fclose(file);
}

/* whitening normalization factor for inner product */
double SP_normalization_factor(double f_low, double f_high, asd_t* asd, stationary_phase_workspace_t *lookup) {
	size_t i;
	double sum;

	/* All values that are not used below are assigned a value of 1, so simply
	 * determine how many there would be, and use that instead of iterating and sum += 1.
	 */
	sum = lookup->f_low_index + (asd->len - lookup->f_high_index + 1);

	for (i = lookup->f_low_index; i <= lookup->f_high_index; i++) {
		double f = asd->f[i];
		double s = asd->asd[i];

		sum  += pow(f, -7.0 / 3.0) / gsl_pow_2(s);
	}

	fprintf(stderr, "g = %0.21e\n", sqrt(sum));

	return sqrt(sum);
}

stationary_phase_workspace_t* SP_workspace_alloc(double f_low, double f_high, asd_t *asd) {
	size_t i;

	stationary_phase_workspace_t *lookup = (stationary_phase_workspace_t*) malloc( sizeof(stationary_phase_workspace_t) );

	lookup->f_low = f_low;
	lookup->f_high = f_high;

	int index_found = 0;
	for (i = 0; i < asd->len; i++) {
		if (asd->f[i] > f_low) {
			lookup->f_low_index = i;
			index_found = 1;
			break;
		}
	}
	if (!index_found) {
		fprintf(stderr, "Error. Unable to find the index for f_low. Aborting.\n");
		abort();
	}

	index_found = 0;
	for (i = 0; i < asd->len; i++) {
		if (asd->f[i] > f_high) {
			lookup->f_high_index = i-1;
			index_found = 1;
			break;
		}
	}
	if (!index_found) {
		fprintf(stderr, "Error. Unable to find the index for f_high. Aborting.\n");
		abort();
	}

	lookup->len = lookup->f_high_index - lookup->f_low_index + 1;

	lookup->g_coeff = (double*) malloc (lookup->len * sizeof(double));
	lookup->chirp_tc_coeff = (double*) malloc (lookup->len * sizeof(double));
	lookup->coalesce_phase_coeff = (double*) malloc (lookup->len * sizeof(double));
	lookup->chirp_time_0_coeff = (double*) malloc (lookup->len * sizeof(double));
	lookup->chirp_time_1_coeff = (double*) malloc (lookup->len * sizeof(double));
	lookup->chirp_time1_5_coeff = (double*) malloc (lookup->len * sizeof(double));
	lookup->chirp_time2_coeff = (double*) malloc (lookup->len * sizeof(double));

	return lookup;
}

void SP_workspace_free( stationary_phase_workspace_t *lookup) {
	free(lookup->g_coeff);
	free(lookup->chirp_tc_coeff);
	free(lookup->coalesce_phase_coeff);
	free(lookup->chirp_time_0_coeff);
	free(lookup->chirp_time_1_coeff);
	free(lookup->chirp_time1_5_coeff);
	free(lookup->chirp_time2_coeff);

	free(lookup);
	lookup = NULL;
}

void SP_workspace_init(double f_low, double f_high, asd_t *asd, stationary_phase_workspace_t *lookup) {
	size_t i, j;

	for (i = lookup->f_low_index, j = 0; i <= lookup->f_high_index; i++, j++) {
		double f = asd->f[i];
		double f_fac = f / f_low;

		lookup->g_coeff[j] = pow(f, -7.0 / 6.0);
		lookup->chirp_tc_coeff[j] = 2.0 * M_PI * f;
		lookup->coalesce_phase_coeff[j] = -M_PI / 4.0;
		lookup->chirp_time_0_coeff[j] = 2.0 * M_PI * f_low * 3.0 * pow(f_fac, -5.0 / 3.0) / 5.0;
		lookup->chirp_time_1_coeff[j] =  pow(f_fac, -5.0 / 3.0);
		lookup->chirp_time1_5_coeff[j] = -3.0 * pow(f_fac,-2.0 / 3.0) / 2.0;
		lookup->chirp_time2_coeff[j] = 3.0 * pow(f_fac, -1.0 / 3.0);
	}
}

void SP_compute(double coalesce_phase, double time_delay,
		inspiral_chirp_time_t *chirp, asd_t *asd,
		double f_low, double f_high,
		double g,
		stationary_phase_workspace_t *lookup,
		stationary_phase_t *out_sp)
{
	size_t i;

	for (i = 0; i < lookup->len; i++) {
		double amp_2pn = (1.0 / g) * lookup->g_coeff[i];

		double phase_2pn =
				lookup->chirp_tc_coeff[i] * (chirp->tc - time_delay)
				- 2.0 * coalesce_phase - lookup->coalesce_phase_coeff[i]
				+ lookup->chirp_time_0_coeff[i] * chirp->chirp_time0
				+ lookup->chirp_time_1_coeff[i] * chirp->chirp_time1
				+ lookup->chirp_time1_5_coeff[i] * chirp->chirp_time1_5
				+ lookup->chirp_time2_coeff[i] * chirp->chirp_time2;

		gsl_complex exp_phase = gsl_complex_exp(gsl_complex_rect(0.0, -phase_2pn));
		out_sp->spa_0[lookup->f_low_index + i] = gsl_complex_mul_real(exp_phase, amp_2pn);
		out_sp->spa_90[lookup->f_low_index + i] = gsl_complex_mul_imag(out_sp->spa_0[lookup->f_low_index + i], -1.0);
	}
}

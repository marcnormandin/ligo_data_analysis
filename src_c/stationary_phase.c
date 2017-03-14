/*
 * stationary_phase.c
 *
 *  Created on: Mar 2, 2017
 *      Author: marcnormandin
 */

#include "stationary_phase.h"

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

stationary_phase_t* SP_malloc(size_t size) {
	stationary_phase_t* sp = (stationary_phase_t*) malloc( sizeof(stationary_phase_t) );
	sp->len = size;
	sp->spa_0 = (gsl_complex*) malloc(size * sizeof(gsl_complex));
	sp->spa_90 = (gsl_complex*) malloc(size * sizeof(gsl_complex));
	return sp;
}

void SP_free(stationary_phase_t* sp) {
	free(sp->spa_0);
	free(sp->spa_90);
	free(sp);
}

void SP_save(char* filename, strain_t* interp_strain, stationary_phase_t* sp) {
	FILE* file;
	size_t i;

	file = fopen(filename, "w");
	for (i = 0; i < sp->len; i++) {
		fprintf(file, "%e %e %e\n", interp_strain->freq[i], GSL_REAL(sp->spa_0[i]), GSL_REAL(sp->spa_90[i]));
	}
	fclose(file);
}

double SP_g(double f_low, double f_high, chirp_time_t* chirp, strain_t* interp_strain) {
	size_t i;
	double sum_tempval;

	/* whitening normalization factor */
	double *tempval = (double*) malloc(interp_strain->len * sizeof(double));
	for (i = 0; i < interp_strain->len; i++) {
		double f = interp_strain->freq[i];
		double s = interp_strain->strain[i];

		if (f > f_low && f < f_high) {
			tempval[i] = pow(f, -7.0 / 3.0) / gsl_pow_2(s);
		} else {
			// Matlab version allocates tempval as "ones(1,length(f))"
			tempval[i] = 1.0;
		}
	}

	sum_tempval = 0.0;
	for (i = 0; i < interp_strain->len; i++) {
		sum_tempval += tempval[i];
	}

	free(tempval);

	return sqrt(gsl_pow_2(chirp->amp_fact_1) * gsl_pow_2(chirp->amp_fact_2) * sum_tempval);
}

void SP_compute(double coalesce_phase, double time_delay,
		chirp_time_t *chirp, strain_t *interp_strain,
		double f_low, double f_high,
		stationary_phase_t *out_sp)
{
	size_t j;

	/* This doesn't change unless the strain changes */
	double g = SP_g(f_low, f_high, chirp, interp_strain);

	/* This is not efficient! */
	for (j = 0; j < interp_strain->len; j++) {
		double f = interp_strain->freq[j];
		double s = interp_strain->strain[j];

		if (f > f_low && f < f_high) {
			double amp_2pn = (1.0 / g) * chirp->amp_fact_1 * chirp->amp_fact_2 * pow(f, -7.0 / 6.0);

			double f_fac = f / f_low;

			double phase_2pn =
					2.0 * M_PI * f * (chirp->tc - time_delay)
					- 2.0 * coalesce_phase - M_PI / 4.0
					+ 2.0 * M_PI * f_low * 3.0 * chirp->chirp_time0 * pow(f_fac, -5.0 / 3.0) / 5.0
					+ chirp->chirp_time1 * pow(f_fac, -5.0 / 3.0)
					- 3.0 * chirp->chirp_time1_5 * pow(f_fac,-2.0 / 3.0) / 2.0
					+ 3.0 * chirp->chirp_time2 * pow(f_fac, -1.0 / 3.0);

			gsl_complex exp_phase = gsl_complex_exp(gsl_complex_rect(0.0, -phase_2pn));
			out_sp->spa_0[j] = gsl_complex_mul_real(exp_phase, amp_2pn);

			out_sp->spa_90[j] = gsl_complex_mul_imag(out_sp->spa_0[0], -1.0);
		} else {
			out_sp->spa_0[j] = gsl_complex_rect(0.0, 0.0);
			out_sp->spa_90[j] = gsl_complex_rect(0.0, 0.0);
		}
	}
}

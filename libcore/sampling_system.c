#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>

#include "sampling_system.h"

int SS_has_nyquist_term(size_t N) {
	return GSL_IS_EVEN(N);
}

/* Returns the Nyquist array index for a C indexed-array */
/* WARNING: Nyquist term is only present if N is EVEN */
size_t SS_nyquist_array_index (size_t N) {
	assert(GSL_IS_EVEN(N));
	return N / 2;
}

/* Returns the last unique index. Either the Nyquist index or last regular term before they are mirrored. */
size_t SS_last_unique_index (size_t N) {
	if (GSL_IS_ODD(N)) {
		return (N+1) / 2;
	} else {
		return N / 2;
	}
}

/* Returns the half-size (side with low frequencies including the DC term) */
size_t SS_half_size(size_t N_full) {
	return SS_last_unique_index( N_full ) + 1;
}

/* Takes a one_sided complex array and adds the corresponding mirrored side. */
void SS_make_two_sided (size_t M, gsl_complex *one_sided, size_t N, gsl_complex *two_sided) {
	size_t m, n;

	/* Check that the dimensions make sense */
	if (GSL_IS_ODD(N) && M != (N+1)/2) {
		/* error */
		printf("SS_make_two_sided failed. One-sided length is (%lu) and two-sided length is (%lu).\n",
				M, N);
		abort();
	} else if (GSL_IS_EVEN(N) && M != (N/2 + 1)) {
		/* error */
		printf("SS_make_two_sided failed. One-sided length is (%lu) and two-sided length is (%lu).\n",
						M, N);
		abort();
	}

	/* copy the left side */
	for (m = 0; m < M; m++) {
		two_sided[m] = one_sided[m];
	}

	/* if a Nyquist term is present, then don't mirror it. */
	size_t c = M - 1;
	if (SS_has_nyquist_term(N)) {
		c--;
	}

	/* make the mirror of the left side */
	for (n = M, m = c; n < N; n++, m--) {
		two_sided[n] = gsl_complex_conjugate( one_sided[m] );
	}
}

void SS_frequency_array(double samplingFrequency, size_t num_total_samples, size_t num_desired_freq_samples, double *frequencies)
{
	size_t i;
	double delta_f = samplingFrequency / (1.0*num_total_samples);
	for (i = 0; i < num_desired_freq_samples; i++) {
		frequencies[i] = delta_f * i;
	}
}

void SS_time_array(double samplingFrequency, size_t num_desired_time_samples, double *times)
{
	size_t i;
	double Ts = 1.0 / samplingFrequency;
	for (i = 0; i < num_desired_time_samples; i++) {
		times[i] = i * Ts;
	}
}

#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_halfcomplex.h>

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
	assert(one_sided != NULL);
	assert(two_sided != NULL);
	assert(M <= N);
	assert(one_sided != two_sided);

	size_t m, n;

	/* Check that the dimensions make sense */
	if (GSL_IS_ODD(N) && M != (N+1)/2) {
		/* error */
		fprintf(stderr, "Error. SS_make_two_sided failed. One-sided length is (%lu) and two-sided length is (%lu). Exiting.\n",
				M, N);
		exit(-1);
	} else if (GSL_IS_EVEN(N) && M != (N/2 + 1)) {
		/* error */
		fprintf(stderr, "Error. SS_make_two_sided failed. One-sided length is (%lu) and two-sided length is (%lu). Exiting.\n",
						M, N);
		exit(-1);
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

/* Takes a one_sided double array and adds the corresponding mirrored side. */
void SS_make_two_sided_real (size_t M, double *one_sided, size_t N, double *two_sided) {
	assert(one_sided != NULL);
	assert(two_sided != NULL);
	assert(M <= N);
	assert(one_sided != two_sided);

	size_t m, n;

	/* Check that the dimensions make sense */
	if (GSL_IS_ODD(N) && M != (N+1)/2) {
		/* error */
		fprintf(stderr, "Error. SS_make_two_sided failed. One-sided length is (%lu) and two-sided length is (%lu). Exiting.\n",
				M, N);
		exit(-1);
	} else if (GSL_IS_EVEN(N) && M != (N/2 + 1)) {
		/* error */
		fprintf(stderr, "Error. SS_make_two_sided failed. One-sided length is (%lu) and two-sided length is (%lu). Exiting.\n",
						M, N);
		exit(-1);
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
		two_sided[n] = one_sided[m];
	}
}

void SS_frequency_array(double samplingFrequency, size_t num_total_samples, size_t num_desired_freq_samples, double *frequencies)
{
	assert(frequencies != NULL);

	size_t i;
	double delta_f = samplingFrequency / (1.0*num_total_samples);
	for (i = 0; i < num_desired_freq_samples; i++) {
		frequencies[i] = delta_f * i;
	}
}

void SS_time_array(double samplingFrequency, size_t num_desired_time_samples, double *times)
{
	assert(times != NULL);

	size_t i;
	double Ts = 1.0 / samplingFrequency;
	for (i = 0; i < num_desired_time_samples; i++) {
		times[i] = i * Ts;
	}
}

/* This function colours a time series by multiplying by the ASD in the frequency domain */
void SS_colour_timeseries( psd_t *psd_one_sided, size_t num_time_samples, double *timeseries) {
	assert(psd_one_sided != NULL);
	assert(timeseries != NULL);

	size_t k, l;

	if (psd_one_sided->type != PSD_ONE_SIDED) {
		fprintf(stderr, "SS_colour_timeseries: Error. Algorithm only works with one_sided PSD. Exiting.\n");
		exit(-1);
	}

	if (SS_half_size(num_time_samples) != psd_one_sided->len) {
		fprintf(stderr, "SS_colour_timeseries: Error. The PSD length doesn't match the half-size of the timeseries. Exiting.");
		exit(-1);
	}

	asd_t * asd_one_sided = ASD_alloc( psd_one_sided->len );
	ASD_init_from_psd( psd_one_sided, asd_one_sided );

	gsl_fft_real_wavetable *fft_real_wavetable = gsl_fft_real_wavetable_alloc( num_time_samples );
	gsl_fft_halfcomplex_wavetable *fft_complex_wavetable = gsl_fft_halfcomplex_wavetable_alloc( num_time_samples );
	gsl_fft_real_workspace *fft_workspace = gsl_fft_real_workspace_alloc( num_time_samples );

	gsl_fft_real_transform(timeseries, 1, num_time_samples, fft_real_wavetable, fft_workspace);

	// DC term doesn't have an imaginary component
	timeseries[0] *= asd_one_sided->asd[0];

	size_t lu = SS_last_unique_index( num_time_samples );
	if (SS_has_nyquist_term(num_time_samples)) {
		lu--;
	}

	/* The GSL fft arranges the data as REAL, IMAG pairs except for the DC and Nyquist Terms. */
	for (k = 1, l = 1; l <= lu; k+=2, l++) {
		double s = asd_one_sided->asd[l] / sqrt(2.0);
		timeseries[k+0] *= s;
		timeseries[k+1] *= s;
	}

	// If nyquist term is present, it doesn't have an imaginary component
	if (SS_has_nyquist_term(num_time_samples)) {
		timeseries[num_time_samples-1] *= asd_one_sided->asd[asd_one_sided->len-1];
	}

	gsl_fft_halfcomplex_inverse( timeseries, 1, num_time_samples, fft_complex_wavetable, fft_workspace );

	gsl_fft_real_workspace_free( fft_workspace );
	gsl_fft_halfcomplex_wavetable_free( fft_complex_wavetable );
	gsl_fft_real_wavetable_free( fft_real_wavetable );
}

/* This function whitens a time series by dividing by the ASD in the frequency domain */
void SS_whiten_timeseries( psd_t *psd_one_sided, size_t num_time_samples, double *timeseries) {
	assert(psd_one_sided != NULL);
	assert(timeseries != NULL);

	size_t k, l;

	if (psd_one_sided->type != PSD_ONE_SIDED) {
		fprintf(stderr, "SS_whiten_timeseries: Error. Algorithm only works with one_sided PSD. Exiting.\n");
		exit(-1);
	}

	if (SS_half_size(num_time_samples) != psd_one_sided->len) {
		fprintf(stderr, "SS_whiten_timeseries: Error. The PSD length doesn't match the half-size of the timeseries. Exiting.");
		exit(-1);
	}

	asd_t * asd_one_sided = ASD_alloc( psd_one_sided->len );
	ASD_init_from_psd( psd_one_sided, asd_one_sided );

	gsl_fft_real_wavetable *fft_real_wavetable = gsl_fft_real_wavetable_alloc( num_time_samples );
	gsl_fft_halfcomplex_wavetable *fft_complex_wavetable = gsl_fft_halfcomplex_wavetable_alloc( num_time_samples );
	gsl_fft_real_workspace *fft_workspace = gsl_fft_real_workspace_alloc( num_time_samples );

	gsl_fft_real_transform(timeseries, 1, num_time_samples, fft_real_wavetable, fft_workspace);

	// DC term doesn't have an imaginary component
	timeseries[0] /= asd_one_sided->asd[0];

	size_t lu = SS_last_unique_index( num_time_samples );
	if (SS_has_nyquist_term(num_time_samples)) {
		lu--;
	}

	/* The GSL fft arranges the data as REAL, IMAG pairs except for the DC and Nyquist Terms. */
	for (k = 1, l = 1; l <= lu; k+=2, l++) {
		double s = asd_one_sided->asd[l] / sqrt(2.0);
		timeseries[k+0] /= s;
		timeseries[k+1] /= s;
	}

	// If nyquist term is present, it doesn't have an imaginary component
	if (SS_has_nyquist_term(num_time_samples)) {
		timeseries[num_time_samples-1] /= asd_one_sided->asd[asd_one_sided->len-1];
	}

	gsl_fft_halfcomplex_inverse( timeseries, 1, num_time_samples, fft_complex_wavetable, fft_workspace );

	gsl_fft_real_workspace_free( fft_workspace );
	gsl_fft_halfcomplex_wavetable_free( fft_complex_wavetable );
	gsl_fft_real_wavetable_free( fft_real_wavetable );
}

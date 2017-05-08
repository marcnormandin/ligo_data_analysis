#ifndef SRC_C_SAMPLING_SYSTEM_H_
#define SRC_C_SAMPLING_SYSTEM_H_

#include <stddef.h>

#include <gsl/gsl_complex.h>

#include "spectral_density.h"

int SS_has_nyquist_term(size_t N);

/* Returns the Nyquist array index for a C indexed-array */
/* WARNING: Nyquist term is only present if N is EVEN */
size_t SS_nyquist_array_index (size_t N);

/* Returns the last unique index. Either the Nyquist index or last regular term before they are mirrored. */
size_t SS_last_unique_index (size_t N);

/* Returns the half-size (side with low frequencies including the DC term) */
size_t SS_half_size(size_t N_full);

/* Takes a one_sided complex array and adds the corresponding mirrored side. */
void SS_make_two_sided (size_t M, gsl_complex *one_sided, size_t N, gsl_complex *two_sided);

/* Takes a one_sided real array and adds the corresponding mirrored side. */
void SS_make_two_sided_real (size_t M, double *one_sided, size_t N, double *two_sided);

/* Write the fft frequencies */
void SS_frequency_array(double samplingFrequency, size_t num_total_samples, size_t num_desired_freq_samples, double *frequencies);

void SS_time_array(double samplingFrequency, size_t num_desired_time_samples, double *times);

/* This function colours a time series by multiplying by the ASD in the frequency domain */
void SS_colour_timeseries( psd_t *psd_one_sided, size_t num_time_samples, double *timeseries);

/* This function whitens a time series by dividing by the ASD in the frequency domain */
void SS_whiten_timeseries( psd_t *psd_one_sided, size_t num_time_samples, double *timeseries);

#endif /* SRC_C_SAMPLING_SYSTEM_H_ */

/*
 * sampling_system.h
 *
 *  Created on: Mar 17, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_SAMPLING_SYSTEM_H_
#define SRC_C_SAMPLING_SYSTEM_H_

#include <stddef.h>
#include <gsl/gsl_complex.h>

int SS_has_nyquist_term(size_t N);

/* Returns the Nyquist array index for a C indexed-array */
/* WARNING: Nyquist term is only present if N is EVEN */
size_t SS_nyquist_array_index (size_t N);

/* Returns the last unique index. Either the Nyquist index or last regular term before they are mirrored. */
size_t SS_last_unique_index (size_t N);

/* Takes a one_sided complex array and adds the corresponding mirrored side. */
void SS_make_two_sided (size_t M, gsl_complex *one_sided, size_t N, gsl_complex *two_sided);

/* Write the fft frequencies */
void SS_frequency_array(double samplingFrequency, size_t N, double *frequencies);

void SS_time_array(double samplingFrequency, size_t N, double *times);

#endif /* SRC_C_SAMPLING_SYSTEM_H_ */

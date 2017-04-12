/*
 * strain_interpolate.c
 *
 *  Created on: Mar 4, 2017
 *      Author: marcnormandin
 */

#include <stdio.h>
#include <gsl/gsl_interp.h>

#include "sampling_system.h"
#include "strain.h"
#include "strain_interpolate.h"

strain_t* InterpStrain_malloc_and_compute(strain_t* strain, double sampling_frequency, size_t num_time_samples)
{
	size_t i, N;
	strain_t* out_interpolated;
	gsl_interp* interp;
	gsl_interp_accel* acc;

	interp = gsl_interp_alloc(gsl_interp_linear, strain->len);
	gsl_interp_init(interp, strain->freq, strain->strain, strain->len);
	acc = gsl_interp_accel_alloc();

	N = SS_half_size(num_time_samples);
	out_interpolated = Strain_malloc(ST_ONE_SIDED, N);

	for (i = 0; i < N; i++) {
		/* The frequency we want to have a strain value for */
		double f = i * (0.5 * sampling_frequency / N);

		/* The interpolated strain value */
		double s = gsl_interp_eval(interp, strain->freq, strain->strain, f,
				acc);

		out_interpolated->freq[i] = f;
		out_interpolated->strain[i] = s;
	}

	gsl_interp_accel_free(acc);
	gsl_interp_free(interp);

	return out_interpolated;
}

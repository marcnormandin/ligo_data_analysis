/*
 * strain_interpolate.c
 *
 *  Created on: Mar 4, 2017
 *      Author: marcnormandin
 */

#include <gsl/gsl_interp.h>

#include "strain.h"
#include "strain_interpolate.h"

strain_t* InterpStrain_malloc_and_compute(strain_t* strain)
{
	size_t i, N;
	strain_t* out_interpolated;
	int samplingFrequency;
	gsl_interp* interp;
	gsl_interp_accel* acc;

	interp = gsl_interp_alloc(gsl_interp_linear, strain->len);
	gsl_interp_init(interp, strain->freq, strain->strain, strain->len);
	acc = gsl_interp_accel_alloc();

	N = 65537;
	out_interpolated = Strain_malloc(N);

	samplingFrequency = 2048;
	for (i = 0; i < N; i++) {
		/* The frequency we want to have a strain value for */
		double f = i * (0.5 * samplingFrequency / N);

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

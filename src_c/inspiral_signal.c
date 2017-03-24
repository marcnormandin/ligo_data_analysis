#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include <gsl/gsl_randist.h>

#include "detector_antenna_patterns.h"
#include "inspiral_chirp.h"
#include "detector.h"
#include "detector_network.h"
#include "inspiral_source.h"
#include "inspiral_stationary_phase.h"
#include "strain.h"
#include "strain_interpolate.h"
#include "signal.h"
#include "inspiral_network_statistic.h"
#include "inspiral_signal.h"

#include "detector_time_delay.h"
#include "simulate_noise.h"

signal_t** simulate_inspiral(gsl_rng *rng, double f_low, double f_high, detector_network_t *net, strain_t *strain, source_t *source) {
	size_t i, j;

	/* generate the unscaled signals */
	signal_t **signals = simulate_inspiral_unscaled(f_low, f_high, net, strain, source);

	/* temp hack */
	for (i = 0; i < net->num_detectors; i++) {
		signal_t* signal = signals[i];
		for (j = 0; j < strain->len; ++j) {
			signal->whitened_data[j] = signal->whitened_signal[j];
		}
	}


	/* compute the network statistic so that we can determine the scale parameter */
	chirp_factors_t chirp;
	CF_compute(f_low, source, &chirp);

	coherent_network_workspace_t *cnw = CN_workspace_malloc(net->num_detectors, strain->len);

	double scale = 0;
	coherent_network_statistic(net, strain, f_low, f_high, &chirp.ct, &source->sky, source->polarization_angle,
			signals, cnw, &scale);

	CN_workspace_free(cnw);

	double scale_factor = source->snr / scale;

	printf("scale = %g, scale_factor = %g\n", scale, scale_factor);

	/* now scale the signals so that they have the desired network snr */
	for (i = 0; i < net->num_detectors; i++) {
		signal_t* signal = signals[i];
		for (j = 0; j < strain->len; ++j) {
			signal->whitened_signal[j] = gsl_complex_mul_real(signal->whitened_signal[j], scale_factor);

			/* temp hack */
			signal->whitened_data[j] = signal->whitened_signal[j];
		}
	}

	return signals;
}

signal_t** simulate_inspiral_unscaled(double f_low, double f_high, detector_network_t *net, strain_t *strain, source_t *source) {
	chirp_factors_t chirp;
	signal_t** signals;
	size_t i;
	stationary_phase_t* sp;
	size_t j;

	CF_compute(f_low, source, &chirp);

	/* one antenna pattern structure per detector */
	antenna_patterns_workspace_t *ap_ws = antenna_patterns_workspace_alloc();
	antenna_patterns_t *ap = (antenna_patterns_t*) malloc( net->num_detectors * sizeof(antenna_patterns_t));
	for (i = 0; i < net->num_detectors; i++) {
		antenna_patterns(net->detector[i], &source->sky, source->polarization_angle, ap_ws, &ap[i]);
	}

	/* One measured signal per detector */
	signals = (signal_t**) malloc( net->num_detectors * sizeof(signal_t*) );
	for (i = 0; i < net->num_detectors; i++) {
		signals[i] = Signal_malloc(strain->len);
	}

	sp = SP_malloc(strain->len);

	/* For each detector determine the stationary phase of the signal */
	for (i = 0; i < net->num_detectors; i++) {
		detector_t* det = net->detector[i];
		signal_t *signal = signals[i];
		double td;
		time_delay(det, &source->sky, &td);

		SP_compute(source->coalesce_phase, td,
				&chirp.ct, strain,
				f_low, f_high,
				sp);

		/* This computes the unscaled (in terms of network snr) signal */
		for (j = 0; j < strain->len; ++j) {
			gsl_complex v;
			gsl_complex A;
			gsl_complex B;

			signal->whitened_sf[j] = gsl_complex_div_real(sp->spa_0[j], strain->strain[j]);
			signal->h_0[j] = signal->whitened_sf[j];

			v = gsl_complex_rect(0.0, -1.0);
			signal->h_90[j] = gsl_complex_mul(signal->h_0[j], v);

			A = gsl_complex_mul_real( signal->h_0[j], ap[i].f_plus );
			B = gsl_complex_mul_real( signal->h_90[j], ap[i].f_cross );
			signal->whitened_signal[j] = gsl_complex_add( A, B );
		}
	}

	antenna_patterns_workspace_free(ap_ws);
	free(ap);

	SP_free(sp);

	return signals;
}

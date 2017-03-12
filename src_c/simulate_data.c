#include "simulate_data.h"

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

#include "antenna_patterns.h"
#include "chirp.h"
#include "detector.h"
#include "detector_network.h"
#include "source.h"
#include "stationary_phase.h"
#include "strain.h"
#include "strain_interpolate.h"
#include "signal.h"
#include "network_analysis.h"

signal_t** simulate_data(gsl_rng *rng, double f_low, double f_high, detector_network_t *net, strain_t *strain, source_t *source) {
	chirp_factors_t chirp;
	CF_compute(f_low, source, &chirp);

	Compute_Detector_Network_Antenna_Patterns(&source->sky, source->polarization_angle, net);

	/* One measured signal per detector */
	signal_t** signals;
	signals = (signal_t**) malloc( net->num_detectors * sizeof(signal_t*) );
	for (int i = 0; i < net->num_detectors; i++) {
		signals[i] = Signal_malloc(strain->len);
	}

	/* !Fixme
	   The hardcoded value should be computed using sqrt( sum (F+^2 + Fx^2) ) * SNR
	   old value: double multi_factor = (1.0 / 2.8580) * 20.0;
	 */
	double snr = source->snr;
	double multi_factor = 0.0;
	for (size_t i = 0; i < net->num_detectors; i++) {
		multi_factor += gsl_pow_2(net->detector[i].ant.f_plus);
		multi_factor += gsl_pow_2(net->detector[i].ant.f_cross);
	}
	multi_factor = 0.5 * sqrt(multi_factor);
	multi_factor *= snr;

	stationary_phase_t* sp = SP_malloc(strain->len);

	/* For each detector determine the stationary phase of the signal */
	for (size_t i = 0; i < net->num_detectors; i++) {
		detector_t* det = &net->detector[i];

		SP_compute(source->coalesce_phase, det->timedelay,
				&chirp.ct, strain,
				f_low, f_high,
				sp);

		signal_t *signal = signals[i];

		/* This computes the signal and signal with noise */
		for (size_t j = 0; j < strain->len; ++j) {
				signal->whitened_sf[j] = gsl_complex_div_real(sp->spa_0[j], strain->strain[j]);
				signal->h_0[j] = signal->whitened_sf[j];

				gsl_complex v = gsl_complex_rect(0.0, -1.0);
				signal->h_90[j] = gsl_complex_mul(signal->h_0[j], v);

				gsl_complex A = gsl_complex_mul_real( signal->h_0[j], det->ant.f_plus );
				gsl_complex B = gsl_complex_mul_real( signal->h_90[j], det->ant.f_cross );
				gsl_complex C = gsl_complex_add( A, B );
				signal->whitened_signal[j] = gsl_complex_mul_real(C, multi_factor );

				/* random noise */
				double noise_f_real = gsl_ran_gaussian(rng, 1.0);
				double noise_f_imag = gsl_ran_gaussian(rng, 1.0);
				gsl_complex noise_f = gsl_complex_rect(noise_f_real, noise_f_imag);

				/* signal + noise */
				signal->whitened_data[j] = gsl_complex_add(signal->whitened_signal[j], noise_f);
		}
	}

	SP_free(sp);

	return signals;
}

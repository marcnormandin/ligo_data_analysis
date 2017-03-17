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
#include "simulate_noise.h"


signal_t** simulate_data(gsl_rng *rng, double f_low, double f_high, detector_network_t *net, strain_t *strain, source_t *source) {
	chirp_factors_t chirp;
	signal_t** signals;
	size_t i;
	double multi_factor;
	stationary_phase_t* sp;
	size_t j;

	CF_compute(f_low, source, &chirp);

	Compute_Detector_Network_Antenna_Patterns(&source->sky, source->polarization_angle, net);

	/* One measured signal per detector */
	signals = (signal_t**) malloc( net->num_detectors * sizeof(signal_t*) );
	for (i = 0; i < net->num_detectors; i++) {
		signals[i] = Signal_malloc(strain->len);
	}

	/* !Fixme
	   The hardcoded value should be computed using sqrt( sum (F+^2 + Fx^2) ) * SNR
	   old value: double multi_factor = (1.0 / 2.8580) * 20.0;
	 */
	multi_factor = (1.0 / 2.8580) * source->snr;
	/*
	for (i = 0; i < net->num_detectors; i++) {
		multi_factor += gsl_pow_2(net->detector[i].ant.f_plus);
		multi_factor += gsl_pow_2(net->detector[i].ant.f_cross);
	}
	multi_factor = 0.5 * sqrt(multi_factor);
	multi_factor *= source->snr;
	*/

	sp = SP_malloc(strain->len);

	/* For each detector determine the stationary phase of the signal */
	for (i = 0; i < net->num_detectors; i++) {
		detector_t* det = &net->detector[i];
		signal_t *signal = signals[i];

		SP_compute(source->coalesce_phase, det->timedelay,
				&chirp.ct, strain,
				f_low, f_high,
				sp);



		/* This computes the signal and signal with noise */
		for (j = 0; j < strain->len; ++j) {
			gsl_complex v;
			gsl_complex A;
			gsl_complex B;
			gsl_complex C;
			double noise_f_real;
			double noise_f_imag;
			gsl_complex noise_f;

			signal->whitened_sf[j] = gsl_complex_div_real(sp->spa_0[j], strain->strain[j]);
			signal->h_0[j] = signal->whitened_sf[j];

			v = gsl_complex_rect(0.0, -1.0);
			signal->h_90[j] = gsl_complex_mul(signal->h_0[j], v);

			A = gsl_complex_mul_real( signal->h_0[j], det->ant.f_plus );
			B = gsl_complex_mul_real( signal->h_90[j], det->ant.f_cross );
			C = gsl_complex_add( A, B );
			signal->whitened_signal[j] = gsl_complex_mul_real(C, multi_factor );

			/* simulated noise */
			noise_f = SN_wn_fd(rng);

			/* signal + noise */
			signal->whitened_data[j] = gsl_complex_add(signal->whitened_signal[j], noise_f);
		}
	}

	SP_free(sp);

	return signals;
}

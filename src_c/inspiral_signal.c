#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <stddef.h>

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
#include "sampling_system.h"

inspiral_signal_half_fft_t* inspiral_signal_half_fft_alloc(size_t full_len) {
	inspiral_signal_half_fft_t *signal = (inspiral_signal_half_fft_t*) malloc( sizeof(inspiral_signal_half_fft_t) );
	signal->full_len = full_len;
	signal->half_fft_len = SS_half_size(signal->full_len);
	signal->half_fft = (gsl_complex*) malloc( signal->half_fft_len * sizeof(gsl_complex) );
	return signal;
}

void inspiral_signal_half_fft_free(inspiral_signal_half_fft_t *signal) {
	free(signal->half_fft);
	free(signal);
	signal = NULL;
}

inspiral_signal_half_fft_t* inspiral_signal_half_fft(double f_low, double f_high, detector_t *det, strain_t *half_strain, source_t *source) {
	chirp_factors_t chirp;
	size_t i;
	stationary_phase_t* sp;
	size_t j;

	CF_compute(f_low, source, &chirp);

	/* one antenna pattern structure per detector */
	antenna_patterns_workspace_t *ap_ws = antenna_patterns_workspace_alloc();
	antenna_patterns_t ap;

	antenna_patterns(det, &source->sky, source->polarization_angle, ap_ws, &ap);

	inspiral_signal_half_fft_t *signal = inspiral_signal_half_fft_alloc(Strain_two_sided_length(half_strain));

	sp = SP_malloc(half_strain->len);
	double td;
	time_delay(det, &source->sky, &td);

	SP_compute(source->coalesce_phase, td,
			&chirp.ct, half_strain,
			f_low, f_high,
			sp);

	/* This computes the unscaled (in terms of network snr) inspiral signal */
	for (j = 0; j < half_strain->len; ++j) {
		gsl_complex v;
		gsl_complex A;
		gsl_complex B;

		gsl_complex whitened_sf = gsl_complex_div_real(sp->spa_0[j], half_strain->strain[j]);
		gsl_complex h_0 = whitened_sf;

		v = gsl_complex_rect(0.0, -1.0);
		gsl_complex h_90 = gsl_complex_mul(h_0, v);

		A = gsl_complex_mul_real( h_0, ap.f_plus );
		B = gsl_complex_mul_real( h_90, ap.f_cross );
		signal->half_fft[j] = gsl_complex_add( A, B );
	}
	antenna_patterns_workspace_free(ap_ws);

	SP_free(sp);

	return signal;
}

inspiral_signal_half_fft_t** simulate_inspiral(gsl_rng *rng, double f_low, double f_high, detector_network_t *net, strain_t *half_strain, source_t *source) {
	size_t i, j;

	/* generate the unscaled signals */
	inspiral_signal_half_fft_t **signals = simulate_inspiral_unscaled(f_low, f_high, net, half_strain, source);

	/* compute the network statistic so that we can determine the scale parameter */
	chirp_factors_t chirp;
	CF_compute(f_low, source, &chirp);
	coherent_network_workspace_t *cnw = CN_workspace_malloc(net->num_detectors, half_strain->len);
	double scale = 0;
	coherent_network_statistic(net, half_strain, f_low, f_high, &chirp.ct, &source->sky, source->polarization_angle,
			signals, cnw, &scale);
	CN_workspace_free(cnw);

	/* now determine the scale factor we need to use in order to get the desired network snr */
	double scale_factor = source->snr / scale;
	printf("scale = %g, scale_factor = %g\n", scale, scale_factor);


	/* now apply scale factor to the signals so that they have the desired network snr */
	for (i = 0; i < net->num_detectors; i++) {
		inspiral_signal_half_fft_t* signal = signals[i];
		for (j = 0; j < half_strain->len; ++j) {
			signal->half_fft[j] = gsl_complex_mul_real(signal->half_fft[j], scale_factor);
		}
	}

	/* add noise to each signal */
	/*
	for(i = 0; i < net->num_detectors; i++) {
		inspiral_signal_half_fft_t *signal = signals[i];
		for (j = 0; j < half_strain->len; ++j) {

			gsl_complex noise = gsl_complex_mul_real(SN_wn_fd(rng), 0.5);

			signal->half_fft[j] = gsl_complex_add(signal->half_fft[j], noise);
		}

	}*/

	return signals;
}

inspiral_signal_half_fft_t** simulate_inspiral_unscaled(double f_low, double f_high, detector_network_t *net, strain_t *half_strain, source_t *source) {
	size_t i;

	/* One measured signal per detector */
	inspiral_signal_half_fft_t **signals = (inspiral_signal_half_fft_t**) malloc( net->num_detectors * sizeof(inspiral_signal_half_fft_t*) );
	for (i = 0; i < net->num_detectors; i++) {
		signals[i] = inspiral_signal_half_fft( f_low, f_high, net->detector[i], half_strain, source);
	}

	return signals;
}

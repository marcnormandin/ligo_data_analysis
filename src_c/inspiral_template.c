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
#include "inspiral_template.h"

#include "detector_time_delay.h"
#include "simulate_noise.h"
#include "sampling_system.h"

inspiral_template_half_fft_t* inspiral_template_half_fft_alloc(size_t num_time_samples) {
	inspiral_template_half_fft_t *signal = (inspiral_template_half_fft_t*) malloc( sizeof(inspiral_template_half_fft_t) );
	signal->full_len = num_time_samples;
	signal->half_fft_len = SS_half_size(signal->full_len);
	signal->half_fft = (gsl_complex*) malloc( signal->half_fft_len * sizeof(gsl_complex) );
	return signal;
}

void inspiral_signal_half_fft_free(inspiral_template_half_fft_t *signal) {
	free(signal->half_fft);
	free(signal);
	signal = NULL;
}

inspiral_template_half_fft_t* inspiral_template_half_fft(double f_low, double f_high, size_t num_time_samples, detector_t *det, asd_t *asd, source_t *source) {
	chirp_factors_t chirp;
	size_t i;
	stationary_phase_t* sp;
	size_t j;

	CF_compute(f_low, source, &chirp);

	/* one antenna pattern structure per detector */
	antenna_patterns_workspace_t *ap_ws = antenna_patterns_workspace_alloc();
	antenna_patterns_t ap;

	antenna_patterns(det, &source->sky, source->polarization_angle, ap_ws, &ap);

	inspiral_template_half_fft_t *signal = inspiral_template_half_fft_alloc( num_time_samples );

	sp = SP_malloc(asd->len);
	double td;
	time_delay(det, &source->sky, &td);

	SP_compute(source->coalesce_phase, td,
			&chirp.ct, asd,
			f_low, f_high,
			sp);

	/* This computes the unscaled (in terms of network snr) inspiral template */
	for (j = 0; j < asd->len; ++j) {
		gsl_complex v;
		gsl_complex A;
		gsl_complex B;

		gsl_complex whitened_sf = gsl_complex_div_real(sp->spa_0[j], asd->asd[j]);
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

inspiral_template_half_fft_t** inspiral_template_unscaled(double f_low, double f_high, size_t num_time_samples, detector_network_t *net, asd_t **net_asd, source_t *source) {
	size_t i;

	/* One measured signal per detector */
	inspiral_template_half_fft_t **signals = (inspiral_template_half_fft_t**) malloc( net->num_detectors * sizeof(inspiral_template_half_fft_t*) );
	for (i = 0; i < net->num_detectors; i++) {
		signals[i] = inspiral_template_half_fft( f_low, f_high, num_time_samples, net->detector[i], net_asd[i], source);
	}

	return signals;
}

inspiral_template_half_fft_t** inspiral_template(double f_low, double f_high, size_t num_time_samples, detector_network_t *net, asd_t **net_asd, source_t *source) {
	size_t i, j;

	/* generate the unscaled signals */
	inspiral_template_half_fft_t **signals = inspiral_template_unscaled(f_low, f_high, num_time_samples, net, net_asd, source);

	/* compute the network statistic so that we can determine the scale parameter */
	chirp_factors_t chirp;
	CF_compute(f_low, source, &chirp);
	/* Assume that each strain has the same length */
	coherent_network_workspace_t *cnw = CN_workspace_malloc(num_time_samples, net->num_detectors, net_asd[0]->len);
	double scale = 0;
	coherent_network_statistic(net, net_asd, f_low, f_high, &chirp.ct, &source->sky, source->polarization_angle,
			signals, cnw, &scale);
	CN_workspace_free(cnw);

	/* now determine the scale factor we need to use in order to get the desired network snr */
	double scale_factor = source->snr / scale;
	/*printf("scale = %g, scale_factor = %g\n", scale, scale_factor);*/


	/* now apply scale factor to the signals so that they have the desired network snr */
	for (i = 0; i < net->num_detectors; i++) {
		inspiral_template_half_fft_t* signal = signals[i];
		for (j = 0; j < net_asd[i]->len; ++j) {
			signal->half_fft[j] = gsl_complex_mul_real(signal->half_fft[j], scale_factor);
		}
	}

	return signals;
}

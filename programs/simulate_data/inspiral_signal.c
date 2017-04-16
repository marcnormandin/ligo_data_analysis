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
#include <gsl/gsl_const_mksa.h>

#include <gsl/gsl_randist.h>

#include "detector_antenna_patterns.h"
#include "inspiral_chirp.h"
#include "inspiral_chirp_time.h"
#include "inspiral_chirp_factors.h"
#include "detector.h"
#include "detector_network.h"
#include "inspiral_stationary_phase.h"
#include "inspiral_signal.h"
#include "inspiral_network_statistic.h"
#include "inspiral_signal.h"

#include "detector_time_delay.h"
//#include "simulate_noise.h"
#include "sampling_system.h"



strain_half_fft_t* inspiral_template_half_fft(double f_low, double f_high, size_t num_time_samples, detector_t *det, source_t *source) {
	inspiral_chirp_factors_t chirp;
	size_t i;
	stationary_phase_workspace_t *sp_lookup;
	stationary_phase_t* sp;
	size_t j;

	CF_compute(f_low, source, &chirp);
	//CF_CT_compute(double f_low, double chirp_time0, double chirp_time1_5, inspiral_chirp_time_t *ct);

	/* one antenna pattern structure per detector */
	detector_antenna_patterns_workspace_t *ap_ws = Detector_Antenna_Patterns_workspace_alloc();
	detector_antenna_patterns_t ap;

	Detector_Antenna_Patterns_compute(det, &source->sky, source->polarization_angle, ap_ws, &ap);

	strain_half_fft_t *signal = strain_half_fft_alloc( num_time_samples );

	sp_lookup = SP_workspace_alloc(f_low, f_high, det->asd->len, det->asd->f);

	sp = SP_malloc(det->asd->len);
	double td;
	Detector_time_delay(det, &source->sky, &td);

	double normalization_factor = SP_normalization_factor(det->asd, sp_lookup);

	SP_compute(td, normalization_factor,
			source->coalesce_phase, &chirp.ct,
			sp_lookup,
			sp);

	/* This computes the unscaled (in terms of network snr) inspiral template */
	for (j = 0; j < det->asd->len; ++j) {
		gsl_complex v;
		gsl_complex A;
		gsl_complex B;

		gsl_complex whitened_sf = gsl_complex_div_real(sp->spa_0[j], det->asd->asd[j]);
		gsl_complex h_0 = whitened_sf;

		v = gsl_complex_rect(0.0, -1.0);
		gsl_complex h_90 = gsl_complex_mul(h_0, v);

		A = gsl_complex_mul_real( h_0, ap.f_plus );
		B = gsl_complex_mul_real( h_90, ap.f_cross );
		signal->half_fft[j] = gsl_complex_add( A, B );
	}
	Detector_Antenna_Patterns_workspace_free(ap_ws);

	SP_workspace_free(sp_lookup);
	SP_free(sp);

	return signal;
}

network_strain_half_fft_t* inspiral_template_unscaled(double f_low, double f_high, size_t num_time_samples, detector_network_t *net, source_t *source) {
	size_t i;

	/* One measured signal per detector */
	network_strain_half_fft_t *network_strain = network_strain_half_fft_alloc( net->num_detectors, num_time_samples );
	for (i = 0; i < net->num_detectors; i++) {
		network_strain->strains[i] = inspiral_template_half_fft( f_low, f_high, num_time_samples, net->detector[i], source);
	}

	return network_strain;
}

network_strain_half_fft_t* inspiral_template(double f_low, double f_high, size_t num_time_samples, detector_network_t *net, source_t *source) {
	size_t i, j;

	/* generate the unscaled signals */
	network_strain_half_fft_t *network_strain = inspiral_template_unscaled(f_low, f_high, num_time_samples, net, source);

	/* compute the network statistic so that we can determine the scale parameter */
	inspiral_chirp_factors_t chirp;
	CF_compute(f_low, source, &chirp);

	/* Assume that each strain has the same length */
	size_t len_asd = net->detector[0]->asd->len;
	coherent_network_workspace_t *cnw = CN_workspace_malloc(num_time_samples, net, len_asd, f_low, f_high);
	double scale = 0;
	coherent_network_statistic(net, f_low, f_high, &chirp.ct, &source->sky, source->polarization_angle,
			network_strain, cnw, &scale);
	CN_workspace_free(cnw);

	/* now determine the scale factor we need to use in order to get the desired network snr */
	double scale_factor = source->snr / scale;
	/*printf("scale = %g, scale_factor = %g\n", scale, scale_factor);*/


	/* now apply scale factor to the signals so that they have the desired network snr */
	for (i = 0; i < net->num_detectors; i++) {
		strain_half_fft_t* signal = network_strain->strains[i];
		for (j = 0; j < net->detector[i]->asd->len; ++j) {
			signal->half_fft[j] = gsl_complex_mul_real(signal->half_fft[j], scale_factor);
		}
	}

	return network_strain;
}


void Source_print(source_t* source) {
	printf("right ascension: %f\n", source->sky.ra);
	printf("declination: %f\n", source->sky.dec);
	printf("polarization angle: %f\n", source->polarization_angle);
	printf("coalesce phase: %f\n", source->coalesce_phase);
	printf("inclination: %f\n", source->inclination_angle);
	printf("binary mass 1: %f\n", source->m1);
	printf("binary mass 2: %f\n", source->m2);
	printf("time of arrival: %f\n", source->time_of_arrival);
	printf("snr: %f\n", source->snr);
}

void Source_load_testsource(source_t* source) {
	source->sky.ra = -2.14;
	source->sky.dec = 0.72;
	source->polarization_angle = 0.0;
	source->coalesce_phase = 0.0;
	source->inclination_angle = 0.0;
	source->m1 = 1.4 * GSL_CONST_MKSA_SOLAR_MASS; /* binary mass 1 */
	source->m2 = 1.4 * GSL_CONST_MKSA_SOLAR_MASS; /* binary mass 2 */
	source->time_of_arrival = 16.0;
	source->snr = 9.0;
}


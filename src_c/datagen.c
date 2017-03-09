#include "datagen.h"

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include <gsl/gsl_rng.h>
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


void Print_Source(source_t* source) {
	printf("right ascension: %f\n", source->sky.ra);
	printf("declination: %f\n", source->sky.dec);
	printf("polarization angle: %f\n", source->polarization_angle);
	printf("coalesce phase: %f\n", source->coalesce_phase);
	printf("inclination: %f\n", source->inclination_angle);
	printf("binary mass 1: %f\n", source->m1);
	printf("binary mass 2: %f\n", source->m2);
	printf("time of arrival: %f\n", source->time_of_arrival);
}

void Load_Source(source_t* source) {
	//source->sky.ra = -2.14;
	//source->sky.dec = 0.72;
	source->sky.ra = -0.5;
	source->sky.dec = 0.1;
	source->polarization_angle = M_PI / 6.0;
	source->coalesce_phase = M_PI / 6.0;
	source->inclination_angle = 0.0;
	source->m1 = 1.4 * GSL_CONST_MKSA_SOLAR_MASS; // binary mass 1
	source->m2 = 4.6 * GSL_CONST_MKSA_SOLAR_MASS; // binary mass 2
	source->time_of_arrival = 32.0;
}

signal_t* Signal_malloc(size_t size) {
	signal_t *s = (signal_t*) malloc( sizeof(signal_t) );
	s->whitened_sf = (gsl_complex*) malloc( size * sizeof(gsl_complex) );
	s->h_0 = (gsl_complex*) malloc( size * sizeof(gsl_complex) );
	s->h_90 = (gsl_complex*) malloc( size * sizeof(gsl_complex) );
	s->whitened_signal = (gsl_complex*) malloc( size * sizeof(gsl_complex) );
	s->whitened_data = (gsl_complex*) malloc( size * sizeof(gsl_complex) );

	return s;
}

void Signal_free(signal_t *s) {
	free(s->whitened_sf);
	free(s->h_0);
	free(s->h_90);
	free(s->whitened_signal);
	free(s->whitened_data);
	free(s);
}

int ComplexFreqArray_save(char* filename, strain_t *strain, gsl_complex *array) {
	FILE* file;
	file = fopen(filename, "w");
	if (file) {
		for (size_t i = 0; i < strain->len; i++) {
			fprintf(file, "%e\t %e\t %e\n", strain->freq[i], GSL_REAL(array[i]), GSL_IMAG(array[i]));
		}
		fclose(file);
		return 0;
	} else {
		printf("Error: Unable to complex frequency array to file (%s).\n", filename);
		return -1;
	}
}

// Must free the memory after using this function
char* concat(const char *s1, const char *s2)
{
    char *result = malloc(strlen(s1)+strlen(s2)+1);
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}

void DataGen() {
	// Random number generator
	const gsl_rng_type *rng_type;
	gsl_rng *rng;
	gsl_rng_env_setup();
	rng_type = gsl_rng_default;
	rng = gsl_rng_alloc(rng_type);
	gsl_rng_set(rng, time(0));

	// Settings
	const double f_low = 40.0; // seismic cutoff. All freqs low set to zero
	const double f_high = 700.0; // most stable inner orbit (last stable orbit related)

	source_t source;
	printf("Inspiral source parameters:\n");
	Load_Source(&source);
	Print_Source(&source);

	printf("\n\nChirp Factors:\n");
	chirp_factors_t chirp;
	CF_compute(f_low, &source, &chirp);
	Print_Chirp_Factors(&chirp);

	detector_network_t net;
	Init_Detector_Network(&net);
	Compute_Detector_Network_Antenna_Patterns(&source.sky, source.polarization_angle, &net);
	//Print_Detector_Network(&net);

	strain_t* irregular_strain = Strain_readFromFile("strain.txt");
	//Strain_print(irregular_strain);

	// find the strains to use at the ends
	double strain_f_low;
	for (int i = 0; i < irregular_strain->len; i++) {
		double f = irregular_strain->freq[i];
		double s = irregular_strain->strain[i];
		if (f >= f_low) {
			strain_f_low = s;
			break;
		}
	}

	double strain_f_high;
	for (int i = 0; i < irregular_strain->len; i++) {
		double f = irregular_strain->freq[i];
		double s = irregular_strain->strain[i];
		if (f >= f_high) {
			strain_f_high = s;
			break;
		}
	}

	// fix the strains
	for (int i = 0; i < irregular_strain->len; i++) {
		double f = irregular_strain->freq[i];
		double s = irregular_strain->strain[i];
		if (f < f_low) {
			irregular_strain->strain[i] = strain_f_low;
		} else if (f > f_high) {
			irregular_strain->strain[i] = strain_f_high;
		}
	}

	strain_t* regular_strain = InterpStrain_malloc_and_compute(irregular_strain);
	//Strain_print(regular_strain);
	Strain_saveToFile("interp.txt", regular_strain);

	// Signal
	signal_t* signals[4];
	for (int i = 0; i < 4; i++) {
		signals[i] = Signal_malloc(regular_strain->len);
	}

	// !Fixme
	// The hardcoded value should be computed using sqrt( sum (F+^2 + Fx^2) ) * SNR
	//double multi_factor = (1.0 / 2.8580) * 20.0;
	double snr = 20.0;
	double multi_factor = 0.0;
	for (size_t i = 0; i < net.num_detectors; i++) {
		multi_factor += gsl_pow_2(net.detector[i].ant.f_plus);
		multi_factor += gsl_pow_2(net.detector[i].ant.f_cross);
	}
	multi_factor = 0.5/sqrt(multi_factor);
	printf("multi_factor = %e\n", multi_factor);
	multi_factor *= snr;

	// For each detector determine the stationary phase of the signal
	for (size_t i = 0; i < net.num_detectors; i++) {
		detector_t* det = &net.detector[i];

		stationary_phase_t* sp = SP_malloc(regular_strain->len);

		SP_compute(source.coalesce_phase, det->timedelay,
				&chirp.ct, regular_strain,
				f_low, f_high,
				sp);

		signal_t *signal = signals[i];

		/* This computes the signal and signal with noise */
		for (size_t j = 0; j < regular_strain->len; ++j) {
			//double f = regular_strain->freq[j];
			//if (f > f_low && f < f_high) {
				signal->whitened_sf[j] = gsl_complex_div_real(sp->spa_0[j], regular_strain->strain[j]);
				signal->h_0[j] = signal->whitened_sf[j];

				gsl_complex v = gsl_complex_rect(0.0, -1.0);
				signal->h_90[j] = gsl_complex_mul(signal->h_0[j], v);

				gsl_complex A = gsl_complex_mul_real( signal->h_0[j], det->ant.f_plus );
				gsl_complex B = gsl_complex_mul_real( signal->h_90[j], det->ant.f_cross );
				gsl_complex C = gsl_complex_add( A, B );
				signal->whitened_signal[j] = gsl_complex_mul_real(C, multi_factor );

				// random noise
				double noise_f_real = gsl_ran_gaussian(rng, 1.0);
				double noise_f_imag = gsl_ran_gaussian(rng, 1.0);
				gsl_complex noise_f = gsl_complex_rect(noise_f_real, noise_f_imag);

				// signal + noise
				signal->whitened_data[j] = gsl_complex_add(signal->whitened_signal[j], noise_f);
		}

		char* fn = concat(det->id, ".whitened_signal");
		ComplexFreqArray_save(fn, regular_strain, signal->whitened_signal);
		free(fn);

		fn = concat(det->id, ".whitened_data");
		ComplexFreqArray_save(fn, regular_strain, signal->whitened_data);

		fn = concat(det->id, ".h_0");
		ComplexFreqArray_save(fn, regular_strain, signal->h_0);
		free(fn);

		fn = concat(det->id, ".h_90");
		ComplexFreqArray_save(fn, regular_strain, signal->h_90);
		free(fn);

		fn = concat(det->id, ".whitened_sf");
		ComplexFreqArray_save(fn, regular_strain, signal->whitened_sf);
		free(fn);

		fn = concat(det->id, ".sp");
		SP_save(fn, regular_strain, sp);
		free(fn);

		SP_free(sp);
	}

	// For the template matching, use time_of_arrival = 0, so tc = t_chirp.
	chirp.ct.tc = chirp.t_chirp;

	double out_val = -1.0;
	coherent_network_statistic(
			&net,
			regular_strain,
			f_low,
			f_high,
			&chirp.ct,
			&source.sky,
			source.polarization_angle,
			signals,
			&out_val);

	printf("out_val = %e\n", out_val);

	for (int i = 0; i < 4; i++) {
		Signal_free(signals[i]);
	}


	Free_Detector_Network(&net);
	Strain_free(irregular_strain);
	Strain_free(regular_strain);

	gsl_rng_free(rng);
}

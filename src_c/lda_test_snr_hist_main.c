/*
 * lda_snr_hist_main.c
 *
 *  Created on: Mar 16, 2017
 *      Author: marcnormandin
 */

#include "antenna_patterns.h"
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
#include "simulate_inspiral.h"

int main(int argc, char* argv[]) {
	size_t n;
	/* Settings */
	const double f_low = 40.0; /* seismic cutoff. */
	const double f_high = 700.0; /* most stable inner orbit (last stable orbit related) */

	source_t source;
	Source_load_testsource(&source);

	detector_network_t net;
	Init_Detector_Network(&net);

	strain_t *strain = Strain_simulated(f_low, f_high);

	/* Random number generator */
	const gsl_rng_type *rng_type;
	gsl_rng *rng;
	gsl_rng_env_setup();
	rng_type = gsl_rng_default;
	rng = gsl_rng_alloc(rng_type);
	gsl_rng_set(rng, time(0));

	size_t num_realizations = 100000;
	double *results_snr = (double*) malloc( num_realizations * sizeof(double) );
	coherent_network_workspace_t *workspace = CN_workspace_malloc( net.num_detectors, strain->len );
	for (n = 0; n < num_realizations; n++) {
		signal_t **signals = simulate_inspiral(rng, f_low, f_high, &net, strain, &source);

		/* For the template matching, use time_of_arrival = 0, so tc = t_chirp. */
		chirp_factors_t chirp;
		CF_compute(f_low, &source, &chirp);
		chirp.ct.tc = chirp.t_chirp;

		double out_val = -1.0;

		coherent_network_statistic(
				&net,
				strain,
				f_low,
				f_high,
				&chirp.ct,
				&source.sky,
				source.polarization_angle,
				signals,
				workspace,
				&out_val);

		printf("%e\n", out_val);
		results_snr[n] = out_val;

		size_t i;
		for (i = 0; i < net.num_detectors; i++) {
			Signal_free(signals[i]);
		}
		free(signals);
	}
	CN_workspace_free( workspace );

	FILE* fid = fopen("snr_hist.dat", "w");
	for (n = 0; n < num_realizations; n++) {
		fprintf(fid, "%e\n", results_snr[n]);
	}
	fclose(fid);

	Free_Detector_Network(&net);
	Strain_free(strain);
	free(results_snr);

	gsl_rng_free(rng);

	return 0;
}

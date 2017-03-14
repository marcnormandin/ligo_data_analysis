#include "simulate_data.h"
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

#include "ptapso_estimate.h"
#include "ptapso_func.h"

#include "snr_sky_map.h"

/*
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
}*/

int test_antenna_patterns() {
	sky_t sky;
	int res;
	double polarization_angle;
	antenna_patterns_t ant;

    sky.dec = 0.5;
    sky.ra = 0.5;
    polarization_angle = 0.5;

    res = antenna_patterns("H1", &sky, polarization_angle, &ant);
    if (res != 0) {
        printf("An error occured!");
        return -1;
    }

    printf("%f %f %f %f\n", ant.u, ant.v, ant.f_plus, ant.f_cross);

    return 0;
}

int old_main(int argc, char* argv[]) {
	size_t n;
	/* Settings */
	const double f_low = 40.0; /* seismic cutoff. */
	const double f_high = 700.0; /* most stable inner orbit (last stable orbit related) */

	source_t source;
	Load_Source(&source);

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
		signal_t **signals = simulate_data(rng, f_low, f_high, &net, strain, &source);

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

		for (int i = 0; i < net.num_detectors; i++) {
			Signal_free(signals[i]);
		}
		free(signals);
	}
	CN_workspace_free( workspace );

	FILE* fid = fopen("results.snr", "w");
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

int main(int argc, char* argv[]) {
	size_t i;

	/* Settings */
	const double f_low = 40.0; /* seismic cutoff */
	const double f_high = 700.0; /* most stable inner orbit (last stable orbit related) */

	source_t source;
	Load_Source(&source);

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

	/* Simulate data for all the detectors composing the network */
	signal_t **signals = simulate_data(rng, f_low, f_high, &net, strain, &source);
	gsl_rng_free(rng);

	coherent_network_workspace_t *workspace = CN_workspace_malloc( net.num_detectors, strain->len );

	/* Setup the parameter structure for the pso fitness function */
	ptapso_fun_params_t params;
	params.f_low = f_low;
	params.f_high = f_high;
	params.network = &net;
	params.signals = signals;
	params.source = &source;
	params.strain = strain;
	params.workspace = workspace;

	printf("The real values are: RA = %f, DEC = %f\n", params.source->sky.ra, params.source->sky.dec);

	snr_sky_map(&params, "snr_sky_map.dat");

	CN_workspace_free( workspace );

	/* Free the data */
	for (i = 0; i < net.num_detectors; i++) {
		Signal_free(signals[i]);
	}
	free(signals);

	Free_Detector_Network(&net);
	Strain_free(strain);


	return 0;
}

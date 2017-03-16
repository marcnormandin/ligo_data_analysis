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

int main(int argc, char* argv[]) {
	size_t i;

	/* Settings */
	const double f_low = 40.0; /* seismic cutoff */
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

	ptapso_estimate(&params);

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

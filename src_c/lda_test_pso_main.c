#ifdef HAVE_CONFIG_H
	#include "config.h"
#endif

#ifdef HAVE_MPI
	#include <mpi.h>
#endif

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

#include "common.h"
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
#include "simulate_inspiral.h"
#include "random.h"
#include "options.h"

void pso_result_save(FILE *fid, pso_result_t *result) {
	fprintf(fid, "%f\t %f\t %f", result->ra, result->dec, result->snr);
}

int main(int argc, char* argv[]) {
	size_t i;

	source_t source;
	Source_load_testsource(&source);
	gslseed_t seed;
	int last_index;

#ifdef HAVE_MPI
	MPI_Init(&argc, &argv);
#endif

	seed = 0;

	/* Allow options to override the defaults settings */
	if (process_command_options(argc, argv, &source, &seed, &last_index) != 0) {
		exit(0);
	}

	/* somehow these need to be set */
	if ((argc - last_index) < 2) {
		printf("last_index = %d\n", last_index);
		printf("argc = %d\n", argc);
		printf("Error: Must supply [num_pso_evaluations] [pso max_steps]!\n");
		exit(-1);
	}

	int arg_num_pso_evaluations = atoi(argv[last_index++]);
	size_t arg_pso_max_steps = atoi(argv[last_index++]);



	/* Settings */
	const double f_low = 40.0; /* seismic cutoff */
	const double f_high = 700.0; /* most stable inner orbit (last stable orbit related) */

	detector_network_t net;
	Init_Detector_Network(&net);

	strain_t *strain = Strain_simulated(f_low, f_high);


	/* Random number generator */
	gsl_rng *rng = random_alloc(seed);

	/* Simulate data for all the detectors composing the network */
	signal_t **signals = simulate_inspiral(rng, f_low, f_high, &net, strain, &source);

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

	gslseed_t *seeds = (gslseed_t*) malloc ( arg_num_pso_evaluations * sizeof(gslseed_t) );
	for (i = 0; i < arg_num_pso_evaluations; i++) {
		seeds[i] = random_seed(rng);
	}

#ifndef HAVE_MPI
	FILE *fid = fopen("pso_results.dat", "w");
	for (i = 0; i < arg_num_pso_evaluations; i++) {
		printf("EVALUATING PSO ESTIMATE #(%lu)...\n", i+1);
		pso_result_t pso_result;
		ptapso_estimate(&params, seeds[i], arg_pso_max_steps, &pso_result);
		pso_result_save(fid, &pso_result);
		if (i < arg_num_pso_evaluations-1) {
			fprintf(fid, "\n");
		}
		printf("ESTIMATE RECORDED.\n\n");
	}
	fclose(fid);
#else
	double buff[3];
	int num_workers;
	int num_jobs = arg_num_pso_evaluations;
	int rank;
	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &num_workers);
	num_workers--; /* Rank 0 doesn't do any pso evaluations */
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int tag = 0;
	if (rank == 0) {
		int num_jobs_done = 0;

		/* Rank 0 will accept the results and write them to file. */
		FILE* fid = fopen("pso_results.dat", "w");
		while (num_jobs_done != num_jobs) {
			MPI_Recv(buff, 3, MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
			pso_result_t pso_result;
			pso_result.ra = buff[0];
			pso_result.dec = buff[1];
			pso_result.snr = buff[2];
			pso_result_save(fid, &pso_result);
			num_jobs_done++;
			if (num_jobs_done != num_jobs) {
				fprintf(fid, "\n");
			}
		}
		fclose(fid);
	} else {
		/* All other ranks are workers. */
		int num_jobs_completed = 0;

		while(1) {
			/* get seed number to process */
			int r = (rank-1) + (num_jobs_completed * num_workers);
			if (r >= num_jobs) {
				break;
			}
			pso_result_t pso_result;
			ptapso_estimate(&params, seeds[r], arg_pso_max_steps, &pso_result);
			buff[0] = pso_result.ra;
			buff[1] = pso_result.dec;
			buff[2] = pso_result.snr;
			MPI_Send(buff, 3, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
			num_jobs_completed++;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

#endif
	free(seeds);

	CN_workspace_free( workspace );

	/* Free the data */
	for (i = 0; i < net.num_detectors; i++) {
		Signal_free(signals[i]);
	}
	free(signals);

	Free_Detector_Network(&net);
	Strain_free(strain);
	random_free(rng);

#ifdef HAVE_MPI
	MPI_Finalize();
#endif

	return 0;
}

#ifdef HAVE_CONFIG_H
	#include "config.h"
#endif

#ifdef HAVE_MPI
	#include <mpi.h>
#endif

#include "detector_antenna_patterns.h"
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

#include "inspiral_pso_fitness.h"
#include "inspiral_signal.h"
#include "random.h"
#include "options.h"
#include "settings.h"

void pso_result_save(FILE *fid, pso_result_t *result) {
	fprintf(fid, "%20.17g %20.17g %20.17g %20.17g %20.17g",
			result->ra, result->dec, result->chirp_t0, result->chirp_t1_5, result->snr);
}

void pso_result_print(pso_result_t *result) {
	printf("%20.17g %20.17g %20.17g %20.17g %20.17g",
			result->ra, result->dec, result->chirp_t0, result->chirp_t1_5, result->snr);
}

int i_am_master() {
#ifdef HAVE_MPI
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	return rank == 0;
#else
	return 1;
#endif
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
	if ((argc - last_index) < 4) {
		printf("last_index = %d\n", last_index);
		printf("argc = %d\n", argc);
		printf("Error: Must supply [input settings file] [input pso settings file] [num pso trials] [pso results file]!\n");
		exit(-1);
	}

	char* arg_settings_file = argv[last_index++];
	char* arg_pso_settings_file = argv[last_index++];
	int arg_num_pso_evaluations = atoi(argv[last_index++]);
	char* arg_pso_results_file = argv[last_index++];

	/* Load the general Settings */
	settings_file_t *settings_file = settings_file_open(arg_settings_file);
	if (settings_file == NULL) {
		printf("Error opening the settings file (%s). Aborting.\n", arg_settings_file);
		abort();
	}

#ifdef HAVE_MPI
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
#endif
	printf("Using the following settings:\n");
	settings_file_print(settings_file);
#ifdef HAVE_MPI
	}
#endif

	const double f_low = atof(settings_file_get_value(settings_file, "f_low"));
	const double f_high = atof(settings_file_get_value(settings_file, "f_high"));
	const double sampling_frequency = atof(settings_file_get_value(settings_file, "sampling_frequency"));
	const double num_time_samples = atoi(settings_file_get_value(settings_file, "num_time_samples"));

	settings_file_close(settings_file);

	/* Compute the chirp factors so that we know the true chirp times, and save them to file. */
#ifdef HAVE_MPI
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
#endif

	inspiral_chirp_factors_t temp_chirp;
	CF_compute(f_low, &source, &temp_chirp);
	FILE *true_parameters = fopen("true_parameters.dat", "w");
	fprintf(true_parameters, "RA DEC CHIRP_TIME_0 CHIRP_TIME_1_5 NETWORK_SNR\n");
	fprintf(true_parameters, "%0.21f %0.21f %0.21f %0.21f %0.21f",
			source.sky.ra, source.sky.dec, temp_chirp.ct.chirp_time0, temp_chirp.ct.chirp_time1_5, source.snr);
	fclose(true_parameters);

#ifdef HAVE_MPI
	}
#endif

	detector_network_t net;
	Init_Detector_Network(&net);

	strain_t *strain = Strain_simulated(f_low, f_high, sampling_frequency, num_time_samples);


	/* Random number generator */
	gsl_rng *rng = random_alloc(seed);

	/* Simulate data for all the detectors composing the network */
	strain_half_fft_t **signals = simulate_inspiral(rng, f_low, f_high, &net, strain, &source);

	coherent_network_workspace_t *workspace = CN_workspace_alloc( net.num_detectors, Strain_one_sided_length(strain) );

	/* Setup the parameter structure for the pso fitness function */
	pso_fitness_function_parameters_t params;
	params.f_low = f_low;
	params.f_high = f_high;
	params.network = &net;
	params.signals = signals;
	params.strain = strain;
	params.workspace = workspace;


	gslseed_t *seeds = (gslseed_t*) malloc ( arg_num_pso_evaluations * sizeof(gslseed_t) );
	for (i = 0; i < arg_num_pso_evaluations; i++) {
		seeds[i] = random_seed(rng);
	}

#ifndef HAVE_MPI

	for (i = 0; i < arg_num_pso_evaluations; i++) {
		printf("EVALUATING PSO ESTIMATE #(%lu)...\n", i+1);

		pso_result_t pso_result;
		pso_estimate_parameters(arg_pso_settings_file, &params, seeds[i], &pso_result);

		FILE *fid = fopen(arg_pso_results_file, "a");
		pso_result_save(fid, &pso_result);
		if (i < arg_num_pso_evaluations) {
			fprintf(fid, "\n");
		}
		fclose(fid);

		pso_result_print(&pso_result);

		printf("\nESTIMATE RECORDED.\n\n");
	}

#else
	double buff[5];
	int num_workers;
	int num_jobs = arg_num_pso_evaluations;
	MPI_Status status;
	MPI_Comm_size(MPI_COMM_WORLD, &num_workers);
	num_workers--; /* Rank 0 doesn't do any pso evaluations */
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int tag = 0;
	if (rank == 0) {
		int num_jobs_done = 0;

		/* Rank 0 will accept the results and write them to file. */
		while (num_jobs_done != num_jobs) {
			MPI_Recv(buff, 3, MPI_DOUBLE, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);

			pso_result_t pso_result;
			pso_result.ra = buff[0];
			pso_result.dec = buff[1];
			pso_result.chirp_t0 = buff[2];
			pso_result.chirp_t1_5 = buff[3];
			pso_result.snr = buff[4];

			FILE *fid = fopen(arg_pso_results_file, "a");
			pso_result_save(fid, &pso_result);
			fclose(fid);

			pso_result_print(&pso_result);
			num_jobs_done++;
			if (num_jobs_done != num_jobs) {
				fprintf(fid, "\n");
			}
		}
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
			pso_estimate_parameters(arg_pso_settings_file, &params, seeds[r], &pso_result);
			buff[0] = pso_result.ra;
			buff[1] = pso_result.dec;
			buff[2] = pso_result.chirp_t0;
			buff[3] = pso_result.chirp_t1_5;
			buff[4] = pso_result.snr;

			MPI_Send(buff, 5, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
			num_jobs_completed++;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

#endif
	free(seeds);

	CN_workspace_free( workspace );

	/* Free the data */
	for (i = 0; i < net.num_detectors; i++) {
		strain_half_fft_free(signals[i]);
	}
	free(signals);

	Detector_Network_free(&net);
	Strain_free(strain);
	random_free(rng);

#ifdef HAVE_MPI
	MPI_Finalize();
	if (rank == 0) {
#else
	printf("program ended successfully.\n");
#endif

#ifdef HAVE_MPI
	}
#endif

	return 0;
}

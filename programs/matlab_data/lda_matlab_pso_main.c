#ifdef HAVE_CONFIG_H
	#include "config.h"
#endif

#ifdef HAVE_MPI
	#include <mpi.h>
#endif

/*#include "gperftools/profiler.h"*/

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

#include "detector_antenna_patterns.h"
#include "inspiral_chirp.h"
#include "detector.h"
#include "detector_network.h"
#include "strain.h"
#include "inspiral_stationary_phase.h"
#include "inspiral_network_statistic.h"

#include "inspiral_pso_fitness.h"
#include "random.h"
#include "settings_file.h"
#include "detector_mapping.h"
#include "hdf5_file.h"
#include "sampling_system.h"

void load_shihan_inspiral_data( const char* hdf_filename, strain_half_fft_t *strain){
	size_t i, j;

	size_t num_time_samples = hdf5_get_num_time_samples( hdf_filename );
	size_t half_size = SS_half_size( num_time_samples );

	double *real_array = (double*) malloc( half_size * sizeof(double) );
	double *imag_array = (double*) malloc( half_size * sizeof(double) );

	hdf5_load_array( hdf_filename, "/shihan/whitened_data_real", real_array);
	hdf5_load_array( hdf_filename, "/shihan/whitened_data_imag", imag_array);

	for (j = 0; j < half_size; j++) {
		strain->half_fft[j] = gsl_complex_rect(real_array[j], imag_array[j]);
	}

	free(real_array);
	free(imag_array);
}

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

	gslseed_t seed;

#ifdef HAVE_MPI
	MPI_Init(&argc, &argv);
#endif

	seed = 1;

	/* somehow these need to be set */
	if (argc != 6) {
		printf("argc = %d\n", argc);
		printf("Error: Must supply [settings file] [detector mapping file] [input pso settings file] [num pso trials] [pso results file]!\n");
		exit(-1);
	}

	char* arg_settings_file = argv[1];
	char* arg_detector_mapping_file = argv[2];
	char* arg_pso_settings_file = argv[3];
	int arg_num_pso_evaluations = atoi(argv[4]);
	char* arg_pso_results_file = argv[5];

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

	settings_file_close(settings_file);

	detector_network_t *net = Detector_Network_load( arg_detector_mapping_file, f_low, f_high );
	detector_network_mapping_t *dmap = Detector_Network_Mapping_load( arg_detector_mapping_file );

	size_t num_time_samples = hdf5_get_num_time_samples( dmap->data_filenames[0] );
	network_strain_half_fft_t *network_strain = network_strain_half_fft_alloc(dmap->num_detectors, num_time_samples );
	for (i = 0; i < net->num_detectors; i++) {
		load_shihan_inspiral_data( dmap->data_filenames[i], network_strain->strains[i] );
	}

	/* Random number generator */
	gsl_rng *rng = random_alloc(seed);

	pso_fitness_function_parameters_t *fitness_function_params =
			pso_fitness_function_parameters_alloc(f_low, f_high, net, network_strain);

	gslseed_t *seeds = (gslseed_t*) malloc ( arg_num_pso_evaluations * sizeof(gslseed_t) );
	for (i = 0; i < arg_num_pso_evaluations; i++) {
		seeds[i] = random_seed(rng);
	}

#ifndef HAVE_MPI

	for (i = 0; i < arg_num_pso_evaluations; i++) {
		printf("EVALUATING PSO ESTIMATE #(%lu)...\n", i+1);

		pso_result_t pso_result;
		pso_estimate_parameters(arg_pso_settings_file, fitness_function_params, seeds[i], &pso_result);

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

	pso_fitness_function_parameters_free(fitness_function_params);

	/* Free the data */
	network_strain_half_fft_free(network_strain);
	/*free(signals);*/

	Detector_Network_free(net);
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

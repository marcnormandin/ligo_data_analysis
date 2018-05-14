#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <gsl/gsl_interp.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>


#include "detector.h"
#include "detector_mapping.h"
#include "detector_network.h"
#include "detector_antenna_patterns.h"
#include "hdf5_file.h"
#include "sampling_system.h"

detector_network_t* Detector_Network_alloc(size_t num_detectors) {
	size_t i;

	detector_network_t *net = (detector_network_t*) malloc( sizeof(detector_network_t) );
	if (net == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory. Exiting.\n");
		exit(-1);
	}

	net->num_detectors = num_detectors;
	net->detector = (detector_t**) malloc(net->num_detectors * sizeof(detector_t*));
	if (net->detector == NULL) {
		fprintf(stderr, "Error. Unable to allocate memory. Exiting.\n");
		exit(-1);
	}

	for (i = 0; i < net->num_detectors; i++) {
		net->detector[i] = Detector_alloc();
	}

	return net;
}

void Detector_Network_free(detector_network_t* net) {
	assert(net != NULL);

	size_t i;
	for (i = 0; i < net->num_detectors; i++) {
		Detector_free(net->detector[i]);
	}
	free(net->detector);
	net->detector = NULL;
}

void Detector_Network_print(detector_network_t* net) {
	assert(net != NULL);

	size_t i;

	printf("DETECTOR NEWTORK: %zu detectors\n", net->num_detectors);
	for (i = 0; i < net->num_detectors; i++) {
		printf("%s ", net->detector[i]->name);
		printf("\n");
	}
}

detector_network_t* Detector_Network_load( const char* detector_mapping_file,
		size_t num_time_samples, double sampling_frequency, double f_low, double f_high ) {
	assert(detector_mapping_file != NULL);

	size_t i, j;

	detector_network_mapping_t *dmap = Detector_Network_Mapping_load( detector_mapping_file );

	/* sampling frequency */
	//double fs = hdf5_get_sampling_frequency( dmap->data_filenames[0] );
	//printf("The sampling frequency is %f.\n", fs);

	//int num_time_samples = hdf5_get_num_time_samples( dmap->data_filenames[0] );
	//printf("Simulated data will contain %d time samples.\n", num_time_samples);

	fprintf(stderr, "Allocating a (%lu) detector network... ", dmap->num_detectors);
	detector_network_t* net = Detector_Network_alloc ( dmap->num_detectors );
	fprintf(stderr, "done.\n");

	for (i = 0; i < net->num_detectors; i++) {
		/* Load the PSD from file */
		fprintf(stderr, "Loading the PSD for detector %lu from file (%s).\n", i, dmap->data_filenames[i]);
		psd_t *psd_unprocessed = PSD_load( dmap->data_filenames[i] );

		psd_t *psd = PSD_make_suitable_for_network_analysis(psd_unprocessed, num_time_samples, sampling_frequency, f_low, f_high);

		detector_t *det = net->detector[i];
		Detector_init_name( dmap->detector_names[i], psd, det);
	}

	printf("GW Detector network created: ");
	for (i = 0; i < net->num_detectors; i++) {
		printf("%s ", net->detector[i]->name);
	}
	printf("\n");

	Detector_Network_Mapping_close(dmap);

	return net;
}

double Detector_Network_condition_number_M(detector_network_t* net, sky_t* sky, double polarization_angle) {
	size_t i;

	detector_antenna_patterns_t ant;
	detector_antenna_patterns_workspace_t* work = Detector_Antenna_Patterns_workspace_alloc();

	gsl_vector* u_vec = gsl_vector_alloc( net->num_detectors );
	gsl_vector* v_vec = gsl_vector_alloc( net->num_detectors );

	gsl_matrix* M = gsl_matrix_alloc(4, 4);
	gsl_matrix_set_zero(M);


	for (i = 0; i < net->num_detectors; i++) {
		Detector_Antenna_Patterns_compute(net->detector[i], sky, polarization_angle, work, &ant);
		gsl_vector_set(u_vec, i, ant.u);
		gsl_vector_set(v_vec, i, ant.v);
	}

	double A, B, C;
	gsl_blas_ddot(u_vec, u_vec, &A);
	gsl_blas_ddot(u_vec, v_vec, &B);
	gsl_blas_ddot(v_vec, v_vec, &C);

	gsl_matrix_set(M,0,0,A);
	gsl_matrix_set(M,0,1,B);
	gsl_matrix_set(M,1,0,B);
	gsl_matrix_set(M,1,1,C);
	gsl_matrix_set(M,2,2,A);
	gsl_matrix_set(M,2,3,B);
	gsl_matrix_set(M,3,2,B);
	gsl_matrix_set(M,3,3,C);


	gsl_matrix* V = gsl_matrix_alloc(4,4);
	gsl_vector* S = gsl_vector_alloc(4);
	gsl_vector* W = gsl_vector_alloc(4);

	gsl_linalg_SV_decomp(M, V, S, W);

	double singular_min = gsl_vector_get(S, 0);
	double singular_max = gsl_vector_get(S, 0);
	for (i = 1; i < 4; i++) {
		double x = gsl_vector_get(S, i);
		if (x < singular_min) {
			singular_min = x;
		}
		if (x > singular_max) {
			singular_max = x;
		}
	}

	gsl_vector_free(W);
	gsl_vector_free(S);
	gsl_matrix_free(V);

	double condition_number = singular_max / singular_min;

	// Clean-up

	gsl_matrix_free(M);

	gsl_vector_free(u_vec);
	gsl_vector_free(v_vec);

	Detector_Antenna_Patterns_workspace_free(work);

	return condition_number;
}

double Detector_Network_condition_number_F(detector_network_t* net, sky_t* sky, double polarization_angle) {
	size_t i;

	detector_antenna_patterns_t ant;
	detector_antenna_patterns_workspace_t* work = Detector_Antenna_Patterns_workspace_alloc();

	const int M = net->num_detectors;
	const int N = 2;

	// Create the antenna pattern matrix
	gsl_matrix* F = gsl_matrix_alloc(M, N);
	gsl_matrix_set_zero(F);
	for (i = 0; i < M; i++) {
		Detector_Antenna_Patterns_compute(net->detector[i], sky, polarization_angle, work, &ant);
		gsl_matrix_set(F,i,0, ant.f_plus);
		gsl_matrix_set(F,i,1, ant.f_cross);
	}

	gsl_matrix* V = gsl_matrix_alloc(N,N);
	gsl_vector* S = gsl_vector_alloc(N);
	gsl_vector* W = gsl_vector_alloc(N);

	gsl_linalg_SV_decomp(F, V, S, W);

	double singular_min = gsl_vector_get(S, 0);
	double singular_max = gsl_vector_get(S, 0);
	for (i = 1; i < N; i++) {
		double x = gsl_vector_get(S, i);
		if (x < singular_min) {
			singular_min = x;
		}
		if (x > singular_max) {
			singular_max = x;
		}
	}

	gsl_vector_free(W);
	gsl_vector_free(S);
	gsl_matrix_free(V);

	double condition_number = singular_max / singular_min;

	// Clean-up
	gsl_matrix_free(F);
	Detector_Antenna_Patterns_workspace_free(work);

	return condition_number;
}



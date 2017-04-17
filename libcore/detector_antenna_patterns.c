#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "detector.h"
#include "detector_antenna_patterns.h"
#include "sky.h"

detector_antenna_patterns_workspace_t* Detector_Antenna_Patterns_workspace_alloc() {
	detector_antenna_patterns_workspace_t *ws = (detector_antenna_patterns_workspace_t*)
			malloc( sizeof(detector_antenna_patterns_workspace_t) );
	if (ws == NULL) {
		fprintf(stderr, "Error. Unable to allocate detector_antenna_patterns_workspace_t. Aborting.\n");
		abort();
	}

	ws->n_hat = gsl_vector_alloc(3);
	ws->ex_i = gsl_vector_alloc(3);
	ws->ey_j = gsl_vector_alloc(3);
	ws->epsilon_plus = gsl_matrix_alloc(3, 3);
	ws->epsilon_cross = gsl_matrix_alloc(3, 3);
	ws->temp = gsl_matrix_alloc(3, 3);
	ws->wp_matrix = gsl_matrix_alloc(2, 2);

	return ws;
}

void Detector_Antenna_Patterns_workspace_free( detector_antenna_patterns_workspace_t* ws ) {
	assert(ws != NULL);
	gsl_vector_free(ws->n_hat);
	gsl_vector_free(ws->ex_i);
	gsl_vector_free(ws->ey_j);
	gsl_matrix_free(ws->epsilon_plus);
	gsl_matrix_free(ws->epsilon_cross);
	gsl_matrix_free(ws->temp);
	gsl_matrix_free(ws->wp_matrix);
	free(ws);
	ws = NULL;
}

static int trace(gsl_matrix* A, double* r) {
	assert(A != NULL);
	assert(r != NULL);

	int i;
	double sum = 0.0;
	size_t n = A->size1;

	if (A->size2 != n) {
		GSL_ERROR("Trace is not defined for non-square matrices.", EDOM);
		return -1;
	}

	for (i = 0; i < n; i++) {
		sum += gsl_matrix_get(A, i, i);
	}
	*r = sum;

	return 0;
}

int Detector_Antenna_Patterns_compute(detector_t *d, sky_t *sky, double polarization_angle,
		detector_antenna_patterns_workspace_t *ws, detector_antenna_patterns_t *ant)
{
	assert(d != NULL);
	assert(sky != NULL);
	assert(ws != NULL);
	assert(ant != NULL);

	/* Define X and Y unit vectors in gw frame. (exEarth)
	   (This is found from rotating about z counterclockwise angle
	   sky->ra and then rotating counterclockwise about y angle
	   sky->dec) */


	init_gsl_vector(
			gsl_sf_cos(sky->ra) * gsl_sf_cos(sky->dec),
			gsl_sf_sin(sky->ra) * gsl_sf_cos(sky->dec),
			gsl_sf_sin(sky->dec),
			ws->n_hat);


	init_gsl_vector(
			gsl_sf_sin(sky->ra),
			-gsl_sf_cos(sky->ra),
			0,
			ws->ex_i);


	init_gsl_vector(
			-gsl_sf_cos(sky->ra) * gsl_sf_sin(sky->dec),
			-gsl_sf_sin(sky->ra) * gsl_sf_sin(sky->dec),
			gsl_sf_cos(sky->dec),
			ws->ey_j);

	/* Define epsilon_plus and epsilon_cross. (Polarization basis tensors)
	   ePlusEarth = ex"*ex - ey"*ey
	   eCrossEarth  = ex"*ey + ey"*ex
	 */

	/* Form a matrix through the outerproduct of two vectors
	   See: https://www.math.utah.edu/software/lapack/lapack-blas/dger.html */

	gsl_matrix_set_zero(ws->epsilon_plus);
	gsl_blas_dger(1.0, ws->ex_i, ws->ex_i, ws->epsilon_plus);
	gsl_blas_dger(-1.0, ws->ey_j, ws->ey_j, ws->epsilon_plus);


	gsl_matrix_set_zero(ws->epsilon_cross);
	gsl_blas_dger(1.0, ws->ex_i, ws->ey_j, ws->epsilon_cross);
	gsl_blas_dger(1.0, ws->ey_j, ws->ex_i, ws->epsilon_cross);


	/* Tensor contraction.(Polarization Angel Independent) */
	gsl_matrix_set_zero(ws->temp);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, d->detector_tensor,
			ws->epsilon_plus, 0.0, ws->temp);
	/* cblas uses: C = alpha * A * B + beta * C */

	trace(ws->temp, &ant->u);

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, d->detector_tensor,
			ws->epsilon_cross, 0.0, ws->temp);
	trace(ws->temp, &ant->v);

	gsl_matrix_set_zero(ws->wp_matrix);
	gsl_matrix_set(ws->wp_matrix, 0, 0, gsl_sf_cos(2 * polarization_angle));
	gsl_matrix_set(ws->wp_matrix, 0, 1, gsl_sf_sin(2 * polarization_angle));
	gsl_matrix_set(ws->wp_matrix, 1, 0, -gsl_sf_sin(2 * polarization_angle));
	gsl_matrix_set(ws->wp_matrix, 1, 1, gsl_sf_cos(2 * polarization_angle));

	ant->f_plus = gsl_matrix_get(ws->wp_matrix, 0, 0) * ant->u
			    + gsl_matrix_get(ws->wp_matrix, 0, 1) * ant->v;
	ant->f_cross = gsl_matrix_get(ws->wp_matrix, 1, 0) * ant->u
			     + gsl_matrix_get(ws->wp_matrix, 1, 1) * ant->v;

	/* Everything is fine */
	return 0;
}

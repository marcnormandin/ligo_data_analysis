/*
 * GW Interferometer Antenna Patterns
 * Written by Marc Normandin, 2017.
 */

#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "antenna_patterns.h"
#include "sky.h"

/* Local helper function */
static void init_vector(gsl_vector* v, double x, double y, double z) {
	gsl_vector_set(v, 0, x);
	gsl_vector_set(v, 1, y);
	gsl_vector_set(v, 2, z);
}

static int trace(gsl_matrix* A, double* r) {
	double sum = 0.0;
	size_t n = A->size1;
	if (A->size2 != n) {
		GSL_ERROR("Trace is not defined for non-square matrices.", EDOM);
		return -1;
	}
	for (int i = 0; i < n; i++) {
		sum += gsl_matrix_get(A, i, i);
	}
	*r = sum;
	return 0;
}

int antenna_patterns(const char *iid, sky_t *sky, double polarization_angle, antenna_patterns_t *ant)
{
	gsl_vector* cordnt_x = gsl_vector_alloc(3);
	gsl_vector* cordnt_y = gsl_vector_alloc(3);

	/* All LIGO coordinates are from LIGO-P000006-D-E: */
	if (strcmp(iid, "H1") == 0 || strcmp(iid, "H2") == 0) {
		/* Generate coordinates for LHO:*/

		/* X-arm unit vector relative to earth centered frame. */
		init_vector(cordnt_x, -0.223891216, 0.799830697, 0.556905359);
		/* Y-arm unit vector relative to earth centered frame. */
		init_vector(cordnt_y, -0.913978490, 0.0260953206, -0.404922650);

	} else if (strcmp(iid, "L1") == 0) {
		/* Generate coordinates for LLO: */
		init_vector(cordnt_x, -0.954574615, -0.141579994, -0.262187738);
		init_vector(cordnt_y, 0.297740169, -0.487910627, -0.820544948);
		/* X-arm unit vector, L1-|cordnt.X| = 3995.15;
		   Y-arm unit vector, L1-|cordnt.Y| = 3995.15; */
	}
	/* Source : www.ligo.org/scientists/GW100916/detectors.txt */
	else if (strcmp(iid, "V1") == 0) {
		init_vector(cordnt_x, -0.70045821479, 0.20848948619, 0.68256166277);
		init_vector(cordnt_y, -0.05379255368, -0.96908180549, 0.24080451708);
	}
	/* Source : www.ligo.org/scientists/GW100916/detectors.txt */
	else if (strcmp(iid, "G1") == 0) {
		init_vector(cordnt_x, -0.445184239, 0.866534205, 0.225675575);
		init_vector(cordnt_y, -0.626000687, -0.552167273, 0.550667271);
	}
	/* Source : www.ligo.org/scientists/GW100916/detectors.txt */
	else if (strcmp(iid, "K1") == 0) {
		init_vector(cordnt_x, -0.4300, -0.8363, 0.3400);
		init_vector(cordnt_y, 0.6821, -0.0542, 0.7292);
	}
	/* Double check following coordinates. At the moment I am not ugsl_sf_sing this detector */
	else if (strcmp(iid, "T1") == 0) {
		init_vector(cordnt_x, 0.648969405, 0.760814505, 0);
		init_vector(cordnt_y, -0.443713769, 0.378484715, -0.812322234);
	} else {
		GSL_ERROR("getifo: The iid you entered is not currently listed.",
				GSL_EINVAL);
		return -1;
	}

	/* Define X and Y unit vectors in gw frame. (exEarth)
	   (This is found from rotating about z counterclockwise angle
	   sky->ra and then rotating counterclockwise about y angle
	   sky->dec) */

	gsl_vector* n_hat = gsl_vector_alloc(3);
	init_vector(n_hat, gsl_sf_cos(sky->ra) * gsl_sf_cos(sky->dec),
			gsl_sf_sin(sky->ra) * gsl_sf_cos(sky->dec),
			gsl_sf_sin(sky->dec));

	gsl_vector* ex_i = gsl_vector_alloc(3);
	init_vector(ex_i, gsl_sf_sin(sky->ra), -gsl_sf_cos(sky->ra),
			0);

	gsl_vector* ey_j = gsl_vector_alloc(3);
	init_vector(ey_j, -gsl_sf_cos(sky->ra) * gsl_sf_sin(sky->dec),
			-gsl_sf_sin(sky->ra) * gsl_sf_sin(sky->dec),
			gsl_sf_cos(sky->dec));

	/* Define epsilon_plus and epsilon_cross. (Polarization basis tensors)
	   ePlusEarth = ex"*ex - ey"*ey
	   eCrossEarth  = ex"*ey + ey"*ex
	 */

	/* Form a matrix through the outerproduct of two vectors
	   See: https://www.math.utah.edu/software/lapack/lapack-blas/dger.html */
	gsl_matrix* epsilon_plus = gsl_matrix_alloc(3, 3);
	gsl_matrix_set_zero(epsilon_plus);
	gsl_blas_dger(1.0, ex_i, ex_i, epsilon_plus);
	gsl_blas_dger(-1.0, ey_j, ey_j, epsilon_plus);

	gsl_matrix* epsilon_cross = gsl_matrix_alloc(3, 3);
	gsl_matrix_set_zero(epsilon_cross);
	gsl_blas_dger(1.0, ex_i, ey_j, epsilon_cross);
	gsl_blas_dger(1.0, ey_j, ex_i, epsilon_cross);

	/* Detector Tensor
	  D = 0.5*(Xhat"*Xhat - Yhat"*Yhat) where Xhat is the unit vector for the
	  "X" arm and Yhat is the unit vector for the "Y" arm
	 */
	gsl_matrix* D = gsl_matrix_alloc(3, 3);
	gsl_matrix_set_zero(D);
	gsl_blas_dger(0.5, cordnt_x, cordnt_x, D);
	gsl_blas_dger(-0.5, cordnt_y, cordnt_y, D);


	/* Tensor contraction.(Polarization Angel Independent) */
	gsl_matrix* temp = gsl_matrix_alloc(3, 3);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, D, epsilon_plus, 0.0, temp);
	/* cblas uses: C = alpha * A * B + beta * C */

	trace(temp, &ant->u);

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, D, epsilon_cross, 0.0,
			temp);
	trace(temp, &ant->v);

	gsl_matrix* wp_matrix = gsl_matrix_alloc(2, 2);
	gsl_matrix_set(wp_matrix, 0, 0, gsl_sf_cos(2 * polarization_angle));
	gsl_matrix_set(wp_matrix, 0, 1, gsl_sf_sin(2 * polarization_angle));
	gsl_matrix_set(wp_matrix, 1, 0, -gsl_sf_sin(2 * polarization_angle));
	gsl_matrix_set(wp_matrix, 1, 1, gsl_sf_cos(2 * polarization_angle));

	ant->f_plus = gsl_matrix_get(wp_matrix, 0, 0) * ant->u
			    + gsl_matrix_get(wp_matrix, 0, 1) * ant->v;
	ant->f_cross = gsl_matrix_get(wp_matrix, 1, 0) * ant->u
			     + gsl_matrix_get(wp_matrix, 1, 1) * ant->v;

	/* Free the memory for the local vectors */
	gsl_vector_free(cordnt_x);
	gsl_vector_free(cordnt_y);
	gsl_vector_free(n_hat);
	gsl_vector_free(ex_i);
	gsl_vector_free(ey_j);
	gsl_matrix_free(epsilon_plus);
	gsl_matrix_free(epsilon_cross);
	gsl_matrix_free(temp);
	gsl_matrix_free(wp_matrix);

	/* Everything is fine */
	return 0;
}

#ifndef SRC_C_ANTENNA_PATTERNS_H_
#define SRC_C_ANTENNA_PATTERNS_H_

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "detector.h"
#include "sky.h"

/* GW Interferometer Antenna Patterns */
typedef struct detector_antenna_patterns_s {
	/* polarization-independent antenna patterns */
	double u, v;

	/* polarization-dependent antenna patterns */
	double f_plus, f_cross;

} detector_antenna_patterns_t;

typedef struct detector_antenna_patterns_workspace_s {
	gsl_vector* n_hat;
	gsl_vector* ex_i;
	gsl_vector* ey_j;
	gsl_matrix* epsilon_plus;
	gsl_matrix* epsilon_cross;
	gsl_matrix* temp;
	gsl_matrix* wp_matrix;

} detector_antenna_patterns_workspace_t;

detector_antenna_patterns_workspace_t* Detector_Antenna_Patterns_workspace_alloc();
void Detector_Antenna_Patterns_workspace_free( detector_antenna_patterns_workspace_t* workspace );

int Detector_Antenna_Patterns_compute(detector_t *d, sky_t *sky, double polarization_angle,
		detector_antenna_patterns_workspace_t *workspace, detector_antenna_patterns_t *ant);

#endif /* SRC_C_ANTENNA_PATTERNS_H_ */

#ifndef LDA_INC
#define LDA_INC

#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

extern int trace(gsl_matrix* A, double* r);

extern int timedelay(const char* Intefe_ID, double declination,
		double right_ascension, double *delay);

extern int antennapattern(double declination, double right_ascention,
		double polarization_angle, const char* Intefe_ID, double *u, double *v,
		double *F_Plus, double *F_Cross);

extern void conditionnum_cal(const char* detId[], double declination_input,
		double rightascension_input, double polarization_angle_input);

#endif // #ifndef LDA_INC

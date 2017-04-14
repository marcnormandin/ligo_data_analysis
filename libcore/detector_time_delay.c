#include <assert.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_vector.h>

#include "detector.h"
#include "detector_time_delay.h"
#include "sky.h"

int Detector_time_delay(detector_t *d, sky_t *sky, double *td)
{
	assert(d != NULL);
	assert(sky != NULL);
	assert(td != NULL);

    double xifo;
    double yifo;
    double zifo;
    double xgw;
    double ygw;
    double zgw;
    double gunit;
    double C;

    /* Define the X, Y & Z WGS-84 coordinates for iid */
    xifo = gsl_vector_get(d->location, 0); /* m */
    yifo = gsl_vector_get(d->location, 1); /* m */
    zifo = gsl_vector_get(d->location, 2); /* m */

    /* Define the X, Y & Z coordinates of the GW"s sky location */
    xgw = gsl_sf_cos(sky->dec)*gsl_sf_cos(sky->ra);
    ygw = gsl_sf_cos(sky->dec)*gsl_sf_sin(sky->ra);
    zgw = gsl_sf_sin(sky->dec);

    /* xgw, ygw and zgw  are unit vector components. So modulus of the vector
       should be one. (For testing purposes) */

    /* Define the GW"s unit vector */
    gunit = sqrt( gsl_pow_2(xgw) + gsl_pow_2(ygw) + gsl_pow_2(zgw) );

    /* Define the speed of light */
    C = GSL_CONST_MKSA_SPEED_OF_LIGHT;

    /* Double check the negative sigh infront. */
    /* Shihan"s expalnation: You don't have to put minus here. Plus or minus will be decided by the sky location */
    *td = (((xifo*xgw)+(yifo*ygw)+(zifo*zgw))/(C));

    return 0;
}

#include <math.h>
#include <string.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "time_delay.h"
#include "sky.h"

/* local helper function */
static void set_wgs(double x, double y, double z, double* wgs) {
    wgs[0] = x;
    wgs[1] = y;
    wgs[2] = z;
}

int time_delay(const char *iid, sky_t *sky, double *td)
{
    /* WGS-84 coordinates of beam splitter */
    double WGS[3];

    double xifo;
    double yifo;
    double zifo;
    double xgw;
    double ygw;
    double zgw;
    double gunit;
    double C;


    /* One coordinate specified in the document for both H1 and H2) */
    if (strcmp(iid,"H1") == 0 || strcmp(iid,"H2") == 0) {
        set_wgs(2.161414928e+6, -3.834695183e+6, 4.600350224e+6, WGS); /* In meters */
    }
    else if (strcmp(iid, "L1") == 0) {
        /* WGS-84 coordinates of "L1" beam splitter */
        set_wgs(-7.427604192e+4, -5.496283721e+6, 3.224257016e+6, WGS);
    }
    else if (strcmp(iid, "V1") == 0) {
        /* WGS-84 coordinates of "V1" beam splitter */
        set_wgs(4.5463741e+6, 8.429897e+5, 4.378577e+6, WGS);
    }
    else if (strcmp(iid, "G1") == 0) {
        /* WGS-84 coordinates of "G1" beam splitter */
        set_wgs(3.8563112e+6,  6.665978e+5, 5.0196406e+6, WGS);
    }
    else if (strcmp(iid, "K1") == 0) {
        /* WGS-84 coordinates of "K1" beam splitter */
        set_wgs(-3.776899062e+6,  3.483900163e+6, 3766657.585, WGS);
    }
    else if (strcmp(iid, "T1") == 0) {
        /* WGS-84 coordinates of "T1" beam splitter */
        set_wgs(-3.946409e6, 3.366259e6, 3.6991507e6, WGS);
    }
    else {
        GSL_ERROR("getifo: The interferometer id you entered is not currently listed.", GSL_EINVAL);
        return -1;
    }

    /* Define the X, Y & Z WGS-84 coordinates for iid */
    xifo = WGS[0]; /* m */
    yifo = WGS[1]; /* m */
    zifo = WGS[2]; /* m */


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

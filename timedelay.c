//// 5/5/2014 
// Calculating the time delay of gravitational waves (with respect to earth centered detector) at different 
// LIGO Detectors according to the sky locations.

//  Shihan Weerathunga
//// Timedelay New : May-23-2016
#include <math.h>
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

int timedelay(const char* Intefe_ID, double declination, double right_ascension, double *delay) 
{
    // global WGS; // WGS-84 coordinates of beam splitter

    // One coordinate specified in the document for both H1 and H2)
    if (strcmp(Intefe_ID,'H1') == 0 || strcmp(Intefe_ID,'H2') == 0) {
        WGS = [-2.161414928e+6, -3.834695183e+6, 4.600350224e+6]; // In meters
    }
    else if (strcmp(Intefe_ID, 'L1') == 0) {
        // WGS-84 coordinates of 'L1' beam splitter
        WGS = [-7.427604192e+4, -5.496283721e+6, 3.224257016e+6];
    }
    else if (strcmp(Intefe_ID, 'V1') == 0) {
        // WGS-84 coordinates of 'V1' beam splitter
        WGS = [4.5463741e+6, 8.429897e+5, 4.378577e+6];
    }
    else if (strcmp(Intefe_ID, 'G1') == 0) {
        // WGS-84 coordinates of 'G1' beam splitter
        WGS = [3.8563112e+6,  6.665978e+5, 5.0196406e+6];
    }
    else if (strcmp(Intefe_ID, 'K1') == 0) {
        // WGS-84 coordinates of 'K1' beam splitter
        WGS = [-3.776899062e+6,  3.483900163e+6, 3766657.585];
    }
    else if (strcmp(Intefe_ID, 'T1') == 0) {
        // WGS-84 coordinates of 'T1' beam splitter
        WGS = [-3.946409e6, 3.366259e6, 3.6991507e6];
    }
    else {
        GSL_ERROR("getifo: The Intefe_ID you entered is not currently listed.", GSL_EINVAL);
        return -1;
    }
    // -----------------------------------------------------------------------

    // Define the X, Y & Z WGS-84 coordinates for Intefe_ID

    const double xifo = WGS[1]; // m
    const double yifo = WGS[2]; // m
    const double zifo = WGS[3]; // m


    // Define the X, Y & Z coordinates of the GW's sky location

    // Unit vector n_hat(right_ascension,declination) projections along  x,y,z.

    //xgw = sin(declination)*cos(right_ascension);
    //ygw = sin(declination)*sin(right_ascension);
    //zgw = cos(declination);

    const double xgw = gsl_sf_cos(declination)*gsl_sf_cos(right_ascension);
    const double ygw = gsl_sf_cos(declination)*gsl_sf_sin(right_ascension);
    const double zgw = gsl_sf_sin(declination);

    // xgw, ygw and zgw  are unit vector components. So modulus of the vector
    // should be one. (For testing purposes)

    //Define the GW's unit vector
    const double gunit = sqrt( (gsl_pow_2(xgw) + gsl_pow_2(ygw) + gsl_pow_2(zgw) );

    // Define the speed of light
    const double C=299792458; // m/s

    //Double check the negative sigh infront.
    // Shihan's expalnation: You don't have to put minus here. Plus or minus will be decided by the sky location
    *delay = (((xifo*xgw)+(yifo*ygw)+(zifo*zgw))/(C));

    return 0;
}

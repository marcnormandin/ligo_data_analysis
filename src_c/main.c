#include "lda.h"

int main(int argc, char* argv[]) {
    double declination = 0.5;
    double right_ascention = 0.5;
    double polarization_angle = 0.5;
    const char* Intefe_ID = "H1";

    double u, v, F_Plus, F_Cross;
    int res = antennapattern(declination, right_ascention, polarization_angle, Intefe_ID,
                &u, &v, &F_Plus, &F_Cross);
    if (res != 0) {
        printf("An error occured!");
        return -1;
    }

    printf("%f %f %f %f\n", u, v, F_Plus, F_Cross);

    return 0;
}

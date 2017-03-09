#include "datagen.h"
#include "antenna_patterns.h"

int test_antenna_patterns() {
	sky_t sky;
    sky.dec = 0.5;
    sky.ra = 0.5;
    double polarization_angle = 0.5;
    const char* iid = "H1";

    antenna_patterns_t ant;
    int res = antenna_patterns(iid, &sky, polarization_angle, &ant);
    if (res != 0) {
        printf("An error occured!");
        return -1;
    }

    printf("%f %f %f %f\n", ant.u, ant.v, ant.f_plus, ant.f_cross);

    return 0;
}

int main(int argc, char* argv[]) {
	DataGen();

	return 0;
}

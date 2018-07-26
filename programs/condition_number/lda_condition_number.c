#ifdef HAVE_CONFIG_H
	#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "detector_network.h"
#include "sky.h"

int main(int argc, char *argv[]) {
	/* somehow these need to be set */
	if (argc != 5) {
		printf("argc = %d\n", argc);
		printf("Error: Usage -> [detector mapping file] [right-ascension, degrees] [declination, degrees] [polarization_angle, degrees]\n");
		exit(-1);
	}

	sky_t sky;

	char* arg_detector_mapping_file = argv[1];
	sky.ra = atof(argv[2]) * M_PI / 180.0;
	sky.dec = atof(argv[3]) * M_PI / 180.0;
	const double polarization_angle = atof(argv[4]) * M_PI / 180.0;

	detector_network_t *net = Detector_Network_load(
			arg_detector_mapping_file, 1024, 2048.0, 10.0, 1000.0 );

	double condition_number_M = Detector_Network_condition_number_M(net, &sky, polarization_angle);
	double condition_number_F = Detector_Network_condition_number_F(net, &sky, polarization_angle);

	Detector_Network_free(net);

	printf("Condition number (using matrix M): %20.17g\n", condition_number_M);
	printf("Condition number (using matrix F): %20.17g\n", condition_number_F);


	return 0;
}

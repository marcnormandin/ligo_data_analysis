#include <stdio.h>
#include <stdlib.h>
#include "spectral_density.h"

int main(int argc, char *argv[]) {

	if (argc != 3) {
		printf("Error: <input data filename> <output psd filename>\n. Exiting.\n");
		exit(-1);
	}

	char *data_filename = argv[1];
	char *output_filename = argv[2];

	int num_time_samples = 131074;
	double sampling_frequency = 2048.0;
	double f_low = 10.0;
	double f_high = 1000.0;

	psd_t *psd_unprocessed = PSD_load( data_filename );
	psd_t *psd = PSD_make_suitable_for_network_analysis(psd_unprocessed, num_time_samples, sampling_frequency, f_low, f_high);

	FILE* fid = fopen(output_filename, "w");
	int i;
	for (i = 0; i < psd->len; i++) {
		fprintf(fid, "%20.17g\t%20.17g", psd->f[i], psd->psd[i]);
		if (i < psd->len-1) {
			fprintf(fid, "\n");
		}
	}
	fclose(fid);

	return 0;
}

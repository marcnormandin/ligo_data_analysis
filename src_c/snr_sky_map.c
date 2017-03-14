/*
 * snr_sky_map.c
 *
 *  Created on: Mar 10, 2017
 *      Author: marcnormandin
 */

#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics_double.h>

#include "snr_sky_map.h"
#include "network_analysis.h"

void snr_sky_map(ptapso_fun_params_t *splParams, const char* output_file) {
	size_t ra_i, dec_i;
	const size_t N_ra = 4;
	const size_t N_dec = 2;

	const double delta_ra = (2.0 * M_PI) / N_ra;
	const double delta_dec = (1.0*M_PI) / N_dec;

	FILE *fid = fopen(output_file, "w");

	for (ra_i = 0; ra_i < N_ra; ra_i++) {
		for (dec_i = 0; dec_i < N_dec; dec_i++) {
			chirp_factors_t chirp;
			double recovered_snr;
			double ra = -M_PI + ra_i * delta_ra;
			double dec = -0.5*M_PI + dec_i * delta_dec;

			/* apply the pso particle location */
			splParams->source->sky.ra = ra;
			splParams->source->sky.dec = dec;

			CF_compute(splParams->f_low, splParams->source, &chirp);
			/* The network statistic requires the time of arrival to be zero
			   in order for the matched filtering to work correctly. */
			chirp.ct.tc = chirp.t_chirp;

			coherent_network_statistic(
					splParams->network,
					splParams->strain,
					splParams->f_low,
					splParams->f_high,
					&chirp.ct,
					&splParams->source->sky,
					splParams->source->polarization_angle,
					splParams->signals,
					splParams->workspace,
					&recovered_snr);

			fprintf(fid, "%f\t %f\t %f\n", ra, dec, recovered_snr);
		}
	}

	fclose(fid);

}

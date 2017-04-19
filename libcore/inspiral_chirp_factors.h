/*
 * inspiral_chirp_factors.h
 *
 *  Created on: Apr 12, 2017
 *      Author: marcnormandin
 */

#ifndef COMMON_INSPIRAL_CHIRP_FACTORS_H_
#define COMMON_INSPIRAL_CHIRP_FACTORS_H_

#include "inspiral_chirp_time.h"

/* This structure holds the values that are computed using the inspiral source parameters.
 * They are stored in a structure because they are repeatedly used throughout the code.
 */
typedef struct inspiral_chirp_factors_s {
	/* f_low, m1, m2 are needed to compute these */
	double total_mass;
	double reduced_mass;
	double chirp_mass;
	double s_mass_ratio;

	/* These all depend on f_low */
	double multi_fac;
	double calculated_reduced_mass;
	double calculated_total_mass;
	double t_chirp;
	double s_mass_ratio_cal;
	double multi_fac_cal;

	/* Needed by the Stationary Phase Approximation */
	inspiral_chirp_time_t ct;

} inspiral_chirp_factors_t;

void CF_print(inspiral_chirp_factors_t* f);
void CF_compute_for_signal(double f_low, double m1, double m2, double time_of_arrival, inspiral_chirp_factors_t *out_cf);

#endif /* COMMON_INSPIRAL_CHIRP_FACTORS_H_ */

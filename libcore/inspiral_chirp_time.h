/*
 * inspiral_chirp_time.h
 *
 *  Created on: Apr 12, 2017
 *      Author: marcnormandin
 */

#ifndef COMMON_INSPIRAL_CHIRP_TIME_H_
#define COMMON_INSPIRAL_CHIRP_TIME_H_

typedef struct inspiral_chirp_time_s {
	double chirp_time0;
	double chirp_time1_5;

	double chirp_time1;
	double chirp_time2;
	double tc;

} inspiral_chirp_time_t;

void CF_CT_compute(double f_low, double chirp_time0, double chirp_time1_5, inspiral_chirp_time_t *ct);

#endif /* COMMON_INSPIRAL_CHIRP_TIME_H_ */

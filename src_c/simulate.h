/*
 * simulate.h
 *
 *  Created on: Mar 17, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_SIMULATE_H_
#define SRC_C_SIMULATE_H_

#include <stddef.h>

typedef struct simulation_s {
	double f_low, f_high;
	double fs;
	size_t num_samples;

} simulation_t;



#endif /* SRC_C_SIMULATE_H_ */

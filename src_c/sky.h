/*
 * sky.h
 *
 *  Created on: Mar 2, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_SKY_H_
#define SRC_C_SKY_H_

/* Sky location */
typedef struct sky_s {
	/* declination */
	double dec;

	/* right-ascension */
	double ra;

} sky_t;


#endif /* SRC_C_SKY_H_ */

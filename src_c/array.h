/*
 * array.h
 *
 *  Created on: Mar 17, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_ARRAY_H_
#define SRC_C_ARRAY_H_

#include <stddef.h>

typedef struct array_s {
	size_t len;
	double* data;

} array_t;

array_t* Array_alloc(size_t N);

void Array_free(array_t *a);


#endif /* SRC_C_ARRAY_H_ */

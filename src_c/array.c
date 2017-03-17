/*
 * array.c
 *
 *  Created on: Mar 17, 2017
 *      Author: marcnormandin
 */


#include "array.h"

#include <stdlib.h>

array_t* Array_alloc(size_t N) {
	array_t *a = (array_t*) malloc( sizeof(array_t) );
	a->len = N;
	a->data = (double*) malloc( a->len * sizeof(double) );
	return a;
}

void Array_free(array_t *a) {
	free(a->data);
	free(a);
}


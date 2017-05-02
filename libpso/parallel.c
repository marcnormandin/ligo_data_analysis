/*
 * parallel.c
 *
 *  Created on: May 2, 2017
 *      Author: marcnormandin
 */

#include "parallel.h"

#ifdef HAVE_CONFIG_H
	#include "config.h"
#endif

#ifdef HAVE_OPENMP
	#include "omp.h"

size_t parallel_get_thread_num() {
	return omp_get_thread_num();
}

size_t parallel_get_max_threads() {
	return omp_get_max_threads();
}

#else

size_t parallel_get_thread_num() {
	return 0;
}

size_t parallel_get_max_threads() {
	return 1;
}

#endif

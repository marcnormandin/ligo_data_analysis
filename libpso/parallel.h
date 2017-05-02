/*
 * parallel.h
 *
 *  Created on: May 2, 2017
 *      Author: marcnormandin
 */

#ifndef LIBPSO_PARALLEL_H_
#define LIBPSO_PARALLEL_H_

#ifdef HAVE_CONFIG_H
	#include "config.h"
#endif

#include <stddef.h>

size_t parallel_get_thread_num();
size_t parallel_get_max_threads();

#endif /* LIBPSO_PARALLEL_H_ */

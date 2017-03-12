/*
 * strain.h
 *
 *  Created on: Mar 2, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_STRAIN_H_
#define SRC_C_STRAIN_H_

#include <stddef.h>

/* The irregularly-spaced LIGO strain, and the regularly-spaced interpolated
 * strain both use this structure.
 */
typedef struct strain_s {
	size_t 	 len;
	double 	*freq;
	double 	*strain;

} strain_t;

strain_t* Strain_malloc(size_t size);
void Strain_free(strain_t* strain);
void Strain_print(strain_t* strain);

/* This allocates memory using Strain_malloc and must be freed using Strain_free */
strain_t* Strain_readFromFile(char* filename);

int Strain_saveToFile(char* filename, strain_t* strain);

/* Memory must be freed using Strain_free() */
strain_t* Strain_simulated(double f_low, double f_high);

#endif /* SRC_C_STRAIN_H_ */

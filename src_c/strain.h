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
typedef enum {
	ST_ONE_SIDED,
	ST_TWO_SIDED,
	ST_IRREGULAR_ONE_SIDED,

} strain_type_e;

typedef struct strain_s {
	strain_type_e type;
	size_t 	 len;
	double 	*freq;
	double 	*strain;

} strain_t;

strain_t* Strain_malloc(strain_type_e t, size_t len);
void Strain_free(strain_t* strain);
void Strain_print(strain_t* strain);

/* This allocates memory using Strain_malloc and must be freed using Strain_free */
strain_t* Strain_readFromFile(char* filename);

int Strain_saveToFile(char* filename, strain_t* strain);

/* Memory must be freed using Strain_free() */
strain_t* Strain_simulated(strain_t *irregular_strain, double f_low, double f_high, double sampling_frequency, size_t num_time_samples);

size_t Strain_two_sided_length(strain_t *strain);
size_t Strain_one_sided_length(strain_t *strain);

#endif /* SRC_C_STRAIN_H_ */

/*
 * strain.c
 *
 *  Created on: Mar 2, 2017
 *      Author: marcnormandin
 */

#include <stdio.h>
#include <stdlib.h>
#include "strain.h"

#include "../../common/sampling_system.h"
#include "strain_interpolate.h"

strain_t* Strain_malloc(strain_type_e t, size_t len) {
	strain_t* strain = (strain_t*) malloc(sizeof(strain_t));
	strain->len = len;
	strain->freq = (double*) malloc(strain->len * sizeof(double));
	strain->strain = (double*) malloc(strain->len * sizeof(double));
	strain->type = t;
	return strain;
}

void Strain_free(strain_t* strain) {
	free(strain->freq);
	free(strain->strain);
	free(strain);
}

void Strain_print(strain_t* strain) {
	size_t i;

	printf("Strain Curve: %ld samples\n", strain->len);
	for (i = 0; i < strain->len; i++) {
		printf("%.10e %.10e\n", strain->freq[i], strain->strain[i]);
	}
}

/* This is only used by Strain_readFile */
static int Read_Num_Strain_Samples(char* filename) {
	FILE* file;
	size_t len;

	file = fopen(filename, "r");
	len = 0;

	if (file) {
		double freq, strain;
		while (fscanf(file, "%lf %lf", &freq, &strain) == 2) {
			len++;
		}
		fclose(file);

		return len;
	} else {
		printf("Error: Unable to open strain file (%s) for reading.\n",
				filename);
		return -1;
	}
}

strain_t* Strain_readFromFile(char* filename) {
	FILE* file;

	strain_t* strain = Strain_malloc(ST_IRREGULAR_ONE_SIDED, Read_Num_Strain_Samples(filename));

	file = fopen(filename, "r");
	if (file) {
		size_t i = 0;
		double f, s;
		while (fscanf(file, "%lf %lf", &f, &s) == 2) {
			strain->freq[i] = f;
			strain->strain[i] = s;
			i++;
		}
		fclose(file);

		return strain;
	} else {
		printf("Error: Unable to read the strain file into memory.\n");
		return NULL;
	}
}

int Strain_saveToFile(char* filename, strain_t* strain) {
	FILE* file;
	size_t i;

	file = fopen(filename, "w");
	if (file) {
		for (i = 0; i < strain->len; i++) {
			fprintf(file, "%0.10e %0.10e\n", strain->freq[i], strain->strain[i]);
		}
		fclose(file);
		return 0;
	} else {
		printf("Error: Unable to save strain data to file (%s).\n", filename);
		return -1;
	}
}

/* Memory must be freed using Strain_free() */
strain_t* Strain_simulated(strain_t *irregular_strain, double f_low, double f_high, double sampling_frequency, size_t num_time_samples) {
	size_t i;
	strain_t* regular_strain;
	double strain_f_low;
	double strain_f_high;

	/*strain_t* irregular_strain = Strain_readFromFile(strain_file);*/
	/*Strain_saveToFile("irregular_strain.txt", irregular_strain);*/

	/* find the strains to use at the ends */
	for (i = 0; i < irregular_strain->len; i++) {
		double f = irregular_strain->freq[i];
		double s = irregular_strain->strain[i];
		if (f >= f_low) {
			strain_f_low = s;
			break;
		}
	}

	for (i = 0; i < irregular_strain->len; i++) {
		double f = irregular_strain->freq[i];
		double s = irregular_strain->strain[i];
		if (f >= f_high) {
			strain_f_high = s;
			break;
		}
	}

	/* extend the strains */
	for (i = 0; i < irregular_strain->len; i++) {
		double f = irregular_strain->freq[i];
		if (f < f_low) {
			irregular_strain->strain[i] = strain_f_low;
		} else if (f > f_high) {
			irregular_strain->strain[i] = strain_f_high;
		}
	}

	regular_strain = InterpStrain_malloc_and_compute(irregular_strain, sampling_frequency, num_time_samples);

	/*Strain_saveToFile("regular_strain.txt", regular_strain);*/

	Strain_free(irregular_strain);

	return regular_strain;
}

/*
size_t Strain_two_sided_length(strain_t *strain) {
	switch (strain->type) {
	case ST_ONE_SIDED:
		return SS_full_size(strain->len);
	case ST_TWO_SIDED:
		return strain->len;
	case ST_IRREGULAR_ONE_SIDED:
		printf("Error: Two-sided length is being called on irregularly spaced strain. This should never happen.\n");
		abort();
	default:
		printf("Error: Two-sided strain encountered unknown strain type.\n");
		abort();
	}
}

size_t Strain_one_sided_length(strain_t *strain) {
	switch (strain->type) {
	case ST_ONE_SIDED:
		return strain->len;
	case ST_TWO_SIDED:
		return SS_full_size(strain->len);
	case ST_IRREGULAR_ONE_SIDED:
		printf("Error: One-sided length is being called on irregularly spaced strain. This should never happen.\n");
		abort();
	default:
		printf("Error: One-sided strain encountered unknown strain type.\n");
		abort();
	}
}

*/


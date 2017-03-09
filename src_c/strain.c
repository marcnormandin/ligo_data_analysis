/*
 * strain.c
 *
 *  Created on: Mar 2, 2017
 *      Author: marcnormandin
 */

#include <stdio.h>
#include <stdlib.h>
#include "strain.h"

strain_t* Strain_malloc(size_t len) {
	strain_t* strain = (strain_t*) malloc(sizeof(strain_t));
	strain->len = len;
	strain->freq = (double*) malloc(strain->len * sizeof(double));
	strain->strain = (double*) malloc(strain->len * sizeof(double));
	return strain;
}

void Strain_free(strain_t* strain) {
	free(strain->freq);
	free(strain->strain);
	free(strain);
}

void Strain_print(strain_t* strain) {
	printf("Strain Curve: %ld samples\n", strain->len);
	for (size_t i = 0; i < strain->len; i++) {
		printf("%e %e\n", strain->freq[i], strain->strain[i]);
	}
}

/* This is only used by Strain_readFile */
static int Read_Num_Strain_Samples(char* filename) {
	FILE* file;
	file = fopen(filename, "r");
	size_t len = 0;

	if (file) {
		double freq, strain;
		while (fscanf(file, "%lf %lf", &freq, &strain) != EOF) {
			len++;
			//printf("%e %e\n", freq, strain);
		}
		//printf("%lu lines\n", len);
		fclose(file);

		return len;
		// now store the val
	} else {
		printf("Error: Unable to open strain file (%s) for reading.\n",
				filename);
		return -1;
	}
}

strain_t* Strain_readFromFile(char* filename) {
	FILE* file;

	strain_t* strain = Strain_malloc(Read_Num_Strain_Samples(filename));

	file = fopen(filename, "r");
	if (file) {
		int i = 0;
		while (fscanf(file, "%lf %lf", &strain->freq[i], &strain->strain[i])
				!= EOF) {
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
	file = fopen(filename, "w");
	if (file) {
		for (size_t i = 0; i < strain->len; i++) {
			fprintf(file, "%e\t %e\n", strain->freq[i], strain->strain[i]);
		}
		fclose(file);
		return 0;
	} else {
		printf("Error: Unable to save strain data to file (%s).\n", filename);
		return -1;
	}
}

/*
 * detector_mapping.h
 *
 *  Created on: Apr 10, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_DETECTOR_MAPPING_H_
#define SRC_C_DETECTOR_MAPPING_H_

#include <stddef.h>

typedef struct detector_mapping_s {
	size_t num_detectors;
	char **detector_names;
	char **data_filenames;

} detector_mapping_t;

detector_mapping_t* detector_mapping_load( const char* filename );
void detector_mapping_close( detector_mapping_t *dm);

#endif /* SRC_C_DETECTOR_MAPPING_H_ */

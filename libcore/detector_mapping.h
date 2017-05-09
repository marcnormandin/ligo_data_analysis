#ifndef SRC_C_DETECTOR_MAPPING_H_
#define SRC_C_DETECTOR_MAPPING_H_

#include <stddef.h>

#if defined (__cplusplus)
extern "C" {
#endif

typedef struct detector_network_mapping_s {
	size_t num_detectors;
	char **detector_names;
	char **data_filenames;

} detector_network_mapping_t;

detector_network_mapping_t* Detector_Network_Mapping_load( const char* filename );
void Detector_Network_Mapping_close( detector_network_mapping_t *dm);

#if defined (__cplusplus)
}
#endif

#endif /* SRC_C_DETECTOR_MAPPING_H_ */

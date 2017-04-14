#ifndef SRC_C_DETECTOR_H_
#define SRC_C_DETECTOR_H_

#include <stddef.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include "spectral_density.h"

#define DETECTOR_MAX_NAME_LENGTH 255

/**
 * Detector identification.
 *
 * Each ID is mapped to unique detector properties such as location, and detector tensor.
 */
typedef enum {
	L1 = 0, /**< enum LIGO Livingston */
	H1, /**< enum LIGO Hanford */
	H2, /**< enum LIGO Hanford (dismantled) */
	V1, /**< enum Virgo */
	K1, /**< enum KAGRA */
	G1, /**< enum GEO300 */
	T1  /**< enum Dunno */
} DETECTOR_ID;


/**
 *  This structure represents a GW interferometer detector.
 */
typedef struct detector_s {
	DETECTOR_ID id; 				/**< Detector ID */

	char name[DETECTOR_MAX_NAME_LENGTH]; 				/**< Detector name */

	gsl_vector *location; 			/**< Position vector of detector arm vertex */

	gsl_vector *arm_x; 				/**< Detector arm direction vector, x */
	gsl_vector *arm_y; 				/**< Detector arm direction vector, y */

	gsl_matrix *detector_tensor; 	/**< Detector tensor */

	psd_t *psd; 					/**< Detector power spectral density */

	asd_t *asd; 					/**< Detector amplitude spectral density */

} detector_t;


detector_t* Detector_alloc();
void Detector_free(detector_t *d);

void Detector_init(DETECTOR_ID name, psd_t *psd, detector_t *d);
void Detector_init_name( char *name, psd_t *psd, detector_t *d);

const char* Detector_id_to_name(DETECTOR_ID id);

DETECTOR_ID Detector_name_to_id(const char* name);

#endif /* SRC_C_DETECTOR_H_ */

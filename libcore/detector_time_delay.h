#ifndef SRC_C_DETECTOR_TIME_DELAY_H_
#define SRC_C_DETECTOR_TIME_DELAY_H_

#include "detector.h"
#include "sky.h"

#if defined (__cplusplus)
extern "C" {
#endif

int Detector_time_delay(detector_t *d, sky_t *sky, double *td);

#if defined (__cplusplus)
}
#endif

#endif /* SRC_C_DETECTOR_TIME_DELAY_H_ */

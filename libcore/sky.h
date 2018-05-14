#ifndef SRC_C_SKY_H_
#define SRC_C_SKY_H_

#if defined (__cplusplus)
extern "C" {
#endif

/* Sky location */
typedef struct sky_s {
	/* right-ascension */
	double ra;

	/* declination */
	double dec;

} sky_t;

#if defined (__cplusplus)
}
#endif

#endif /* SRC_C_SKY_H_ */

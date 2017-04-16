/*
 * libcore_test.c
 *
 *  Created on: Apr 15, 2017
 *      Author: marcnormandin
 */

#ifdef HAVE_CONFIG_H
	#include "config.h"
#endif

#ifdef HAVE_GTEST
	#include <gtest/gtest.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../libcore/hdf5_file.c"
#include "../libcore/inspiral_stationary_phase.c"
#include "../libcore/spectral_density.c"

#ifdef HAVE_GTEST
TEST(find_index_low, left_end) {
	size_t N = 100;
	double f_array[N];
	for (int i = 0; i < N; i++) {
		f_array[i] = i;
	}

	double f_low = 0.0;
	int index = find_index_low(f_low, N, f_array);

	EXPECT_EQ(index, 0);
}

TEST(find_index_low, right_end) {
	size_t N = 100;
	double f_array[N];
	for (int i = 0; i < N; i++) {
		f_array[i] = i;
	}

	double f_low = N-1;
	int index = find_index_low(f_low, N, f_array);

	EXPECT_EQ(index, N-1);
}

TEST(find_index_low, not_exact) {
	size_t N = 100;
	double f_array[N];
	for (int i = 0; i < N; i++) {
		f_array[i] = i;
	}

	// f_array = [0, 1, 2, 3]
	double f = 1.1;

	// this should return the index representing 2, which is array index 2
	int index = find_index_low(f, N, f_array);

	EXPECT_EQ(index, 2);
}

TEST(find_index_low, not_found) {
	size_t N = 100;
	double f_array[N];
	for (int i = 0; i < N; i++) {
		f_array[i] = i;
	}

	double f_low = N;
	int index = find_index_low(f_low, N, f_array);

	EXPECT_EQ(index, -1);
}

TEST(find_index_high, left_end) {
	size_t N = 100;
	double f_array[N];
	for (int i = 0; i < N; i++) {
		f_array[i] = i;
	}

	double f = 0.0;
	int index = find_index_high(f, N, f_array);

	EXPECT_EQ(index, 0);
}

TEST(find_index_high, right_end) {
	size_t N = 100;
	double f_array[N];
	for (int i = 0; i < N; i++) {
		f_array[i] = i;
	}

	double f = N-1;
	int index = find_index_high(f, N, f_array);

	EXPECT_EQ(index, N-1);
}

TEST(find_index_high, not_exact) {
	size_t N = 100;
	double f_array[N];
	for (int i = 0; i < N; i++) {
		f_array[i] = i;
	}
	// f_array = [0, 1, 2, 3]

	double f = 2.3;
	int index = find_index_high(f, N, f_array);

	EXPECT_EQ(index, 2);
}

TEST(find_index_high, not_found) {
	size_t N = 100;
	double f_array[N];
	for (int i = 0; i < N; i++) {
		f_array[i] = i;
	}

	double f = N;
	int index = find_index_high(f, N, f_array);

	EXPECT_EQ(index, -1);
}

TEST(SP_workspace, frequencies_and_indices) {
	size_t N = 10;
	asd_t *asd = ASD_alloc( N );
	for (int i = 0; i < N; i++) {
		asd->f[i] = i;
		asd->asd[i] = 1.0;
	}

	double f_low = 2.0;
	double f_high = 6.0;

	stationary_phase_workspace_t *sp_workspace = SP_workspace_alloc(f_low, f_high, asd->len, asd->f);

	EXPECT_EQ( sp_workspace->f_low_index, 2);
	EXPECT_EQ( sp_workspace->f_high_index, 6);
	EXPECT_EQ( sp_workspace->f_low, f_low);
	EXPECT_EQ( sp_workspace->f_high, f_high);
	EXPECT_EQ( sp_workspace->len, 5);

	SP_workspace_free(sp_workspace);
	ASD_free(asd);
}

TEST(SP_workspace, coefficients) {
	size_t N = 10;
	asd_t *asd = ASD_alloc( N );
	for (int i = 0; i < N; i++) {
		asd->f[i] = i;
		asd->asd[i] = 1.0;
	}

	double f_low = 1.0;
	double f_high = 6.0;

	stationary_phase_workspace_t *sp_workspace = SP_workspace_alloc(f_low, f_high, asd->len, asd->f);

	/* Manually change the asd to be uniform so that we can calculate the coefficients easily
	 * since all the powers will be 1.0.
	 */
	for (int i = 0; i < asd->len; i++) {
		asd->f[i] = 1.0;
	}

	/* reinit */
	SP_workspace_init( asd->len, asd->f, sp_workspace );

	for (int i = 0; i < sp_workspace->len; i++) {
		EXPECT_EQ( sp_workspace->g_coeff[i],  1.0);
		EXPECT_EQ( sp_workspace->chirp_tc_coeff[i], 2.0 * M_PI );
		EXPECT_EQ( sp_workspace->constant_coeff[i], -M_PI / 4.0 );
		EXPECT_EQ( sp_workspace->chirp_time_0_coeff[i], 2.0 * M_PI * 3.0 / 5.0 );
		EXPECT_EQ( sp_workspace->chirp_time_1_coeff[i],  1.0 );
		EXPECT_EQ( sp_workspace->chirp_time1_5_coeff[i], -3.0 / 2.0 );
		EXPECT_EQ( sp_workspace->chirp_time2_coeff[i], 3.0 );
	}

	EXPECT_EQ( sp_workspace->f_low_index, 1);
	EXPECT_EQ( sp_workspace->f_high_index, 6);
	EXPECT_EQ( sp_workspace->f_low, f_low);
	EXPECT_EQ( sp_workspace->f_high, f_high);
	EXPECT_EQ( sp_workspace->len, 6);

	SP_workspace_free(sp_workspace);
	ASD_free(asd);
}

TEST(SP_compute, value) {
	double detector_time_delay = 1.0;
	double detector_normalization_factor = 2.0;
	double inspiral_coalesce_phase = 3.0;
	inspiral_chirp_time_t ct;
	ct.chirp_time0 = 4.0;
	ct.chirp_time1 = 5.0;
	ct.chirp_time1_5 = 6.0;
	ct.chirp_time2 = 7.0;
	ct.tc = 8.0;
	double f_low = 1.0;
	double f_high = 6.0;

	size_t len_f_array = 10;
	double f_array[10];
	for (int i = 0; i < len_f_array; i++) {
		f_array[i] = i;
	}

	stationary_phase_workspace_t *w = SP_workspace_alloc(f_low, f_high, len_f_array, f_array);
	stationary_phase_t *s = SP_alloc(10);
	SP_compute(detector_time_delay, detector_normalization_factor,
			inspiral_coalesce_phase, &ct,
			w,
			s);

	for (int i = 0; i < s->len; i++) {
		printf("%0.21e \t %0.21e\n", GSL_REAL(s->spa_0[i]), GSL_IMAG(s->spa_0[i]));
	}

	SP_free(s);
	SP_workspace_free(w);

}

#endif


#include "lda.h"

#include <math.h>
#include <string.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

//// Generating 2 nd postnewtonian signal templates for correlation operations in frequency domain.
//// Signal phase is 2nd postnewtonian.

int templategen(const char* Id, const double chirp_time0_input,
		const double chirp_time1_5_input, double rightascension_input,
		double declination_input, const double g, const size_t len_f,
		const double* f, const double f_low, const double f_high, double* SPA_0,
		double* SPA_90) {

	//Velocity of Light(c) ; Gravitational Constant(G) ; Solar mass
	const double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
	const double G = GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT;

	//solar_mass = 1.9891e30; //polarization_angle  = 0;   //solar_mass = 1.9891e30;

	const double chirp_time0 = chirp_time0_input;
	const double chirp_time1_5 = chirp_time1_5_input;
	const double declination = declination_input;         //Declination Angle
	const double right_ascension = rightascension_input; //Right Ascension Angle

	// Time delay calculations relative to earth centered detector.
	double delay;
	timedelay(Id, declination, right_ascension, &delay);

	const double Amp_fact1 = 1.0;
	const double Amp_fact2 = 1.0;

	//Calculated reduced masses using  chirp_time0 and chirp_time1_5
	const double calculated_reduced_mass = (1.0 / (16.0 * gsl_pow_2(f_low)))
			* pow(
					5.0
							/ (4.0 * gsl_pow_4(M_PI) * chirp_time0
									* gsl_pow_2(chirp_time1_5)), 1.0 / 3.0)
			/ (G / gsl_pow_3(c));
	//Calculated total mass  using  chirp_time0 and chirp_time1_5
	const double calculated_total_mass = (5.0 / (32.0 * f_low))
			* (chirp_time1_5 / (gsl_pow_2(M_PI) * chirp_time0))
			/ (G / gsl_pow_3(c));

	const double s_mass_ratio_cal = calculated_reduced_mass
			/ calculated_total_mass;
	const double multi_fac_cal = f_low * M_PI * G * calculated_total_mass
			/ (gsl_pow_3(c));

	const double chirp_time1 = (5.0 / (192.0 * M_PI)) / f_low / multi_fac_cal
			/ s_mass_ratio_cal
			* ((743.0 / 336.0) + ((11.0 / 4.0) * s_mass_ratio_cal));
	const double chirp_time2 = (5.0 / (128.0 * M_PI)) / f_low
			/ pow(multi_fac_cal, 1.0 / 3.0) / s_mass_ratio_cal
			* ((3058673.0 / 1016064.0) + ((5429.0 / 1008.0) * s_mass_ratio_cal)
					+ ((617.0 / 144.0) * (gsl_pow_2(s_mass_ratio_cal))));

	//Total Chirp time (Duration of the chirp signal)
	const double T_chirp = chirp_time0 + chirp_time1 - chirp_time1_5
			+ chirp_time2;

	//// Template Generation.
	const double time_of_arrival = 0.0;     // Time of arrival = 0 for templates
	const double tc_template = time_of_arrival + T_chirp;    //coalescence time.
	const double coalesce_phase = 0.0;

	// Pre allocation of memory for faster computing.
	//amp_SPA = ones(1, len_f);     phase_SPA = ones(1,len_f);
	//SPA_0 = ones(1, len_f);       SPA_90= ones(1,len_f);

	for (int jj = 0; jj < len_f; ++jj) {

		if (f_low < f[jj] && f[jj] <= f_high) {
			////  h_0 template.

			const double amp_SPA = (1.0 / g) * Amp_fact1 * Amp_fact2
					* (pow(f[jj], (-7 / 6)));
			const double f_fac = f[jj] / f_low;

			const double phase_SPA = 2.0 * M_PI * f[jj] * (tc_template - delay)
					- 2.0 * coalesce_phase - M_PI / 4.0
					+ +(2.0 * M_PI * f_low)
							* ((3.0 * chirp_time0 * (pow(f_fac, -5.0 / 3.0))
									/ 5.0)
									+ (chirp_time1 * pow(f_fac, -5.0 / 3.0))
									- (3.0 * chirp_time1_5
											* pow(f_fac, -2.0 / 3.0) / 2.0)
									+ (3.0 * chirp_time2
											* pow(f_fac, -1.0 / 3.0)));

			//
			// gsl_complex gsl_complex_exp (gsl_complex z)
			SPA_0[jj] = amp_SPA * exp(-1i * phase_SPA);

			////  h_M_PI/2 template.
			SPA_90[jj] = -1i * SPA_0[jj];
		} else {
			SPA_0[jj] = 0.0;
			SPA_90[jj] = 0.0;
		}
	}

	return 0;
}


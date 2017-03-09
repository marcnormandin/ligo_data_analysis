#include "notused_templategen.h"

#include "lda.h"
#include <math.h>
#include <string.h>
#include <stddef.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

int templategen(char* iid, template_parameters_t* tp, double g, frequency_t* freq, template_t* tem) {

	//Velocity of Light(c) ; Gravitational Constant(G) ; Solar mass
	const double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
	const double G = GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT;

	// Time delay calculations relative to earth centered detector.
	double time_delay;
	timedelay(iid, &tp->sky, &time_delay);

	const double Amp_fact1 = 1.0;
	const double Amp_fact2 = 1.0;

	//Calculated reduced masses using  chirp_time0 and chirp_time1_5
	const double calculated_reduced_mass = Chirp_Calc_CalculatedReducedMass(freq->f_low, G, tp->chirp_time0, tp->chirp_time1_5);
	const double calculated_total_mass = Chirp_Calc_CalculatedTotalMass(freq->f_low, G, tp->chirp_time0, tp->chirp_time1_5);
	const double s_mass_ratio_cal = Chirp_Calc_SMassRatioCal(calculated_reduced_mass, calculated_total_mass);
	const double multi_fac_cal = Chirp_Calc_MultiFacCal(freq->f_low, G, calculated_total_mass);
	const double chirp_time1 = Chirp_Calc_Time1(freq->f_low, multi_fac_cal, s_mass_ratio_cal);
	const double chirp_time2 = Chirp_Calc_Time2(freq->f_low, multi_fac_cal, s_mass_ratio_cal);

	//Total Chirp time (Duration of the chirp signal)
	const double T_chirp = Chirp_Calc_TChirp(tp->chirp_time0, chirp_time1, tp->chirp_time1_5, chirp_time2);

	const double tc_template = Chirp_Calc_TC(tp->time_of_arrival, T_chirp);

	for (int jj = 0; jj < freq->len; ++jj) {
		if (freq->f_low < freq->strain.freq[jj] && freq->strain.freq[jj] <= freq->f_high) {
			const double amp_SPA = (1.0 / g) * Amp_fact1 * Amp_fact2 * (pow(freq->strain.freq[jj], (-7 / 6)));
			const double f_fac = freq->strain.freq[jj] / freq->f_low;

			const double phase_SPA = 2.0 * M_PI * freq->strain.freq[jj] * (tc_template - time_delay)
					- 2.0 * tp->coalesce_phase - M_PI / 4.0
					+ (2.0 * M_PI * freq->f_low)
							* ((3.0 * tp->chirp_time0 * (pow(f_fac, -5.0 / 3.0))
									/ 5.0)
									+ (chirp_time1 * pow(f_fac, -5.0 / 3.0))
									- (3.0 * tp->chirp_time1_5
											* pow(f_fac, -2.0 / 3.0) / 2.0)
									+ (3.0 * chirp_time2
											* pow(f_fac, -1.0 / 3.0)));

			tem->spa_0[jj] = amp_SPA * exp(-1i * phase_SPA);
			tem->spa_90[jj] = -1i * tem->spa_0[jj];
		} else {
			tem->spa_0[jj] = 0.0;
			tem->spa_90[jj] = 0.0;
		}
	}

	return 0;
}


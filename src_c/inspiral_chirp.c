/*
 * chirp.c
 *
 *  Created on: Mar 2, 2017
 *      Author: marcnormandin
 */

#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_const_mksa.h>

#include "inspiral_chirp.h"

static double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
static double G = GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT;

void Print_Chirp_Factors(chirp_factors_t* f) {
	printf("total mass: %f\n", f->total_mass);
	printf("reduced mass: %f\n", f->reduced_mass);
	printf("chirp mass: %f\n", f->chirp_mass);
	printf("s mass ratio: %f\n", f->s_mass_ratio);
	printf("multi fac: %f\n", f->multi_fac);

	printf("chirp time 0: %f\n", f->ct.chirp_time0);
	printf("chirp time 1_5: %f\n", f->ct.chirp_time1_5);
	printf("calculated reduced mass: %f\n", f->calculated_reduced_mass);
	printf("calculated total mass: %f\n", f->calculated_total_mass);
	printf("s mass ratio cal: %f\n", f->s_mass_ratio_cal);
	printf("multi fac cal: %f\n", f->multi_fac_cal);
	printf("chirp time 1: %f\n", f->ct.chirp_time1);
	printf("chirp time 2: %f\n", f->ct.chirp_time2);
	printf("t chirp: %f\n", f->t_chirp);
	printf("tc: %f\n", f->ct.tc);
}

double Chirp_Calc_TotalMass(double m1, double m2) {
	return m1 + m2;
}

double Chirp_Calc_ReducedMass(double m1, double m2, double total_mass) {
	return m1 * m2 / (total_mass);
}

double Chirp_Calc_ChirpMass(double reduced_mass, double total_mass) {
	return pow(reduced_mass, 3.0 / 5.0)
			* pow(total_mass, 2.0 / 5.0);
}

double Chirp_Calc_SMassRatio(double reduced_mass, double total_mass) {
	return reduced_mass / total_mass;
}

double Chirp_Calc_SMassRatioCal(double calculated_reduced_mass, double calculated_total_mass) {
	return calculated_reduced_mass
			/ calculated_total_mass;
}

double Chirp_Calc_MultiFac(double f_low, double total_mass) {
	return f_low * M_PI * G * total_mass / gsl_pow_3(c);
}

double Chirp_Calc_ChirpTime0(double f_low, double multi_fac, double s_mass_ratio) {
	return (5.0 / (256.0 * M_PI)) / f_low * pow(multi_fac, -5.0 / 3.0) / s_mass_ratio;
}

double Chirp_Calc_ChirpTime1_5(double f_low, double multi_fac, double s_mass_ratio) {
	return (1.0 / 8.0) / f_low * pow(multi_fac, -2.0 / 3.0) / s_mass_ratio;
}

double Chirp_Calc_CalculatedReducedMass(double f_low, double chirp_time0, double chirp_time1_5) {
	return (1.0 / (16.0 * gsl_pow_2(f_low)))
			* pow(
					5.0
							/ (4.0 * gsl_pow_4(M_PI) * chirp_time0
									* gsl_pow_2(chirp_time1_5)), 1.0 / 3.0)
			/ (G / gsl_pow_3(c));
}

double Chirp_Calc_CalculatedTotalMass(double f_low, double chirp_time0, double chirp_time1_5) {
	return (5.0 / (32.0 * f_low))
			* (chirp_time1_5 / (gsl_pow_2(M_PI) * chirp_time0))
			/ (G / gsl_pow_3(c));
}

double Chirp_Calc_MultiFacCal(double f_low, double calculated_total_mass) {
	return f_low * M_PI * G * calculated_total_mass / gsl_pow_3(c);
}

double Chirp_Calc_Time1(double f_low, double multi_fac_cal, double s_mass_ratio_cal) {
	return (5.0 / (192.0 * M_PI)) / f_low / multi_fac_cal
				/ s_mass_ratio_cal
				* ((743.0 / 336.0) + ((11.0 / 4.0) * s_mass_ratio_cal));
}

double Chirp_Calc_Time2(double f_low, double multi_fac_cal, double s_mass_ratio_cal) {
	return (5.0 / (128.0 * M_PI)) / f_low
			* pow(multi_fac_cal, -1.0 / 3.0) / s_mass_ratio_cal
			* ((3058673.0 / 1016064.0)
					+ ((5429.0 / 1008.0) * s_mass_ratio_cal)
					+ ((617.0 / 144.0) * gsl_pow_2(s_mass_ratio_cal)));
}

double Chirp_Calc_TChirp(double chirp_time0, double chirp_time1, double chirp_time1_5, double chirp_time2) {
	return chirp_time0 + chirp_time1 - chirp_time1_5
			+ chirp_time2;
}

double Chirp_Calc_TC(double time_of_arrival, double t_chirp) {
	return time_of_arrival + t_chirp;
}


void CF_compute(double f_low, source_t *source, chirp_factors_t *fac) {
	fac->total_mass = Chirp_Calc_TotalMass(source->m1, source->m2);
	fac->reduced_mass = Chirp_Calc_ReducedMass(source->m1, source->m2, fac->total_mass);
	fac->chirp_mass = Chirp_Calc_ChirpMass(fac->reduced_mass, fac->total_mass);

	/* Change parameters accordingly to estimate chirp times. */
	fac->s_mass_ratio = Chirp_Calc_SMassRatio(fac->reduced_mass, fac->total_mass);
	fac->multi_fac = Chirp_Calc_MultiFac(f_low, fac->total_mass);
	fac->ct.chirp_time0 = Chirp_Calc_ChirpTime0(f_low, fac->multi_fac, fac->s_mass_ratio);
	fac->ct.chirp_time1_5 = Chirp_Calc_ChirpTime1_5(f_low, fac->multi_fac, fac->s_mass_ratio);

	/* Calculated reduced masses using  chirp_time0 and chirp_time1_5 */
	fac->calculated_reduced_mass = Chirp_Calc_CalculatedReducedMass(f_low, fac->ct.chirp_time0, fac->ct.chirp_time1_5);
	fac->calculated_total_mass = Chirp_Calc_CalculatedTotalMass(f_low, fac->ct.chirp_time0, fac->ct.chirp_time1_5);
	fac->s_mass_ratio_cal = Chirp_Calc_SMassRatioCal(fac->calculated_reduced_mass, fac->calculated_total_mass);
	fac->multi_fac_cal = Chirp_Calc_MultiFacCal(f_low, fac->calculated_total_mass);
	fac->ct.chirp_time1 = Chirp_Calc_Time1(f_low, fac->multi_fac_cal, fac->s_mass_ratio_cal);
	fac->ct.chirp_time2 = Chirp_Calc_Time2(f_low, fac->multi_fac_cal, fac->s_mass_ratio_cal);
	fac->t_chirp = Chirp_Calc_TChirp(fac->ct.chirp_time0, fac->ct.chirp_time1, fac->ct.chirp_time1_5, fac->ct.chirp_time2);

	fac->ct.tc = Chirp_Calc_TC(source->time_of_arrival, fac->t_chirp);
}

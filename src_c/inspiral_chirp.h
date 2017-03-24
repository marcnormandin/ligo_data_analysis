/*
 * chirp.h
 *
 *  Created on: Mar 2, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_CHIRP_H_
#define SRC_C_CHIRP_H_

#include <stddef.h>
#include "inspiral_source.h"

typedef struct chrip_time_s {
	double chirp_time0;
	double chirp_time1;
	double chirp_time1_5;
	double chirp_time2;
	double tc;

} chirp_time_t;

/* This structure holds the values that are computed using the inspiral source parameters.
 * They are stored in a structure because they are repeatedly used throughout the code.
 */
typedef struct chirp_factors_s {
	/* m1, m2 are needed to compute these */
	double total_mass;
	double reduced_mass;
	double chirp_mass;
	double s_mass_ratio;

	/* These all depend on f_low */
	double multi_fac;
	double calculated_reduced_mass;
	double calculated_total_mass;
	double t_chirp;
	double s_mass_ratio_cal;
	double multi_fac_cal;

	/* Needed by the Stationary Phase Approximation */
	chirp_time_t ct;

} chirp_factors_t;


double Chirp_Calc_TotalMass(double m1, double m2);

double Chirp_Calc_ReducedMass(double m1, double m2, double total_mass);

double Chirp_Calc_ChirpMass(double reduced_mass, double total_mass);

double Chirp_Calc_SMassRatio(double reduced_mass, double total_mass);

double Chirp_Calc_MultiFac(double f_low, double total_mass);

double Chirp_Calc_ChirpTime0(double f_low, double multi_fac, double s_mass_ratio);

double Chirp_Calc_ChirpTime1_5(double f_low, double multi_fac, double s_mass_ratio);

double Chirp_Calc_CalculatedReducedMass(double f_low, double chirp_time0, double chirp_time1_5);

double Chirp_Calc_CalculatedTotalMass(double f_low, double chirp_time0, double chirp_time1_5);

double Chirp_Calc_SMassRatioCal(double calculated_reduced_mass, double calculated_total_mass);

double Chirp_Calc_MultiFacCal(double f_low, double calculated_total_mass);

double Chirp_Calc_Time1(double f_low, double multi_fac_cal, double s_mass_ratio_cal);

double Chirp_Calc_Time2(double f_low, double multi_fac_cal, double s_mass_ratio_cal);

double Chirp_Calc_TChirp(double chirp_time0, double chirp_time1, double chirp_time1_5, double chirp_time2);

double Chirp_Calc_TC(double time_of_arrival, double t_chirp);

void Print_Chirp_Factors(chirp_factors_t* f);

void CF_compute(double f_low, source_t *source, chirp_factors_t *out_cf);



#endif /* SRC_C_CHIRP_H_ */

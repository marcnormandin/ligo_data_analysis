/*
 * inspiral_chirp_factors.c
 *
 *  Created on: Apr 12, 2017
 *      Author: marcnormandin
 */

#include <stdio.h>

#include "inspiral_chirp.h"
#include "inspiral_chirp_factors.h"

void CF_print(inspiral_chirp_factors_t* f) {
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

void CF_compute_for_signal(double f_low, double m1, double m2, double time_of_arrival, inspiral_chirp_factors_t *fac) {
	fac->total_mass = Chirp_Calc_TotalMass(m1, m2);
	fac->reduced_mass = Chirp_Calc_ReducedMass(m1, m2, fac->total_mass);
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

	fac->ct.tc = Chirp_Calc_TC(time_of_arrival, fac->t_chirp);
}

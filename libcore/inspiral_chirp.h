/*
 * chirp.h
 *
 *  Created on: Mar 2, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_CHIRP_H_
#define SRC_C_CHIRP_H_


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


#endif /* SRC_C_CHIRP_H_ */

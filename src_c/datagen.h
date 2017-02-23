/*
 * datagen.h
 *
 *  Created on: Feb 22, 2017
 *      Author: marcnormandin
 */

#ifndef SRC_C_DATAGEN_H_
#define SRC_C_DATAGEN_H_


// earth-centered frame
typedef struct {
	double right_ascension;
	double declination;
	double polarization_angle;
	double coalesce_phase;
	double inclination;
	double m1; // binary mass 1
	double m2; // binary mass 2
	double time_of_arrival;

} inspiral_source_parameters_t;


typedef struct {
	double total_mass;
	double reduced_mass;
	double chirp_mass;
	double amp_fact_1;
	double amp_fact_2;
	double s_mass_ratio;
	double multi_fac;

	double chirp_time0;
	double chirp_time1_5;
	double calculated_reduced_mass;
	double calculated_total_mass;
	double s_mass_ratio_cal;
	double multi_fac_cal;
	double chirp_time1;
	double chirp_time2;
	double t_chirp;
	double tc;

	size_t len_tchirp;
	double* tchirp;

} inspiral_chirp_factors_t;


typedef struct {
	char id[3];

	double u, v;
	double f_plus, f_cross;

	double timedelay;
} detector_t;


typedef struct {
	int num_detectors;
	detector_t* detector;

} detector_network_t;

typedef struct {
	size_t len;
	double* freq;
	double* strain;

} strain_t;

void DataGen();

int Read_Num_Strain_Samples(char* filename);
int Load_Strain(char* filename, strain_t* strain);
void Print_Strain(strain_t* strain);
void Free_Strain(strain_t* strain);
void Compute_Interpolated_Strain(strain_t* strain, strain_t* interpolated);
void Alloc_Strain(size_t len, strain_t* strain);
void Free_Strain(strain_t* interpolated);
void Save_Strain(char* filename, strain_t* strain);

void Alloc_Detector_Network(int num, detector_network_t* net);
void Free_Detector_Network(detector_network_t* net);
void Init_Detector_Network(detector_network_t* net);
void Compute_Detector_Network_Antenna_Patterns(inspiral_source_parameters_t* source, detector_network_t* net);
void Print_Detector_Network(detector_network_t* net);
void Print_Detector(detector_t* det);

void Print_Source(inspiral_source_parameters_t* source);
void Load_Source(inspiral_source_parameters_t* source);

void Print_Chirp_Factors(inspiral_chirp_factors_t* f);

void Free_Chirp_Factors(inspiral_chirp_factors_t* fac);
void Compute_Chirp_Factors(
		double f_low,
		inspiral_source_parameters_t* source,
		inspiral_chirp_factors_t* fac);


#endif /* SRC_C_DATAGEN_H_ */

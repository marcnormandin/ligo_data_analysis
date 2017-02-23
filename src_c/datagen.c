#include "lda.h"
#include "datagen.h"

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_interp.h>

void Alloc_Detector_Network(int num, detector_network_t* net) {
	net->num_detectors = num;
	net->detector = (detector_t*) malloc(
			net->num_detectors * sizeof(detector_t));
}

void Free_Detector_Network(detector_network_t* net) {
	free(net->detector);
}

void Init_Detector_Network(detector_network_t* net) {
	Alloc_Detector_Network(4, net);

	strncpy(net->detector[0].id, "H1", 2);
	strncpy(net->detector[1].id, "L1", 2);
	strncpy(net->detector[2].id, "V1", 2);
	strncpy(net->detector[3].id, "K1", 2);
}

void Compute_Detector_Network_Antenna_Patterns(
		inspiral_source_parameters_t* source, detector_network_t* net) {
	for (size_t i = 0; i < net->num_detectors; i++) {
		detector_t* det = &net->detector[i];

		// Fixme!
		// Check return value for errors
		antennapattern(source->declination, source->right_ascension,
				source->polarization_angle, det->id, &det->u, &det->v,
				&det->f_plus, &det->f_cross);

		// Time delay
		timedelay(det->id, source->declination, source->right_ascension,
				&det->timedelay);
	}
}

void Print_Detector(detector_t* det) {
	printf("DETECTOR:\n");
	printf("id: %s\n", det->id);
	printf("u: %f\n", det->u);
	printf("v: %f\n", det->v);
	printf("f_plus: %f\n", det->f_plus);
	printf("f_cross: %f\n", det->f_cross);
}

void Print_Detector_Network(detector_network_t* net) {
	printf("DETECTOR NEWTORK: %d detectors\n", net->num_detectors);
	for (size_t i = 0; i < net->num_detectors; i++) {
		Print_Detector(&net->detector[i]);
		printf("\n");
	}
}

void Print_Source(inspiral_source_parameters_t* source) {
	printf("right ascension: %f\n", source->right_ascension);
	printf("declination: %f\n", source->declination);
	printf("polarization angle: %f\n", source->polarization_angle);
	printf("coalesce phase: %f\n", source->coalesce_phase);
	printf("inclination: %f\n", source->inclination);
	printf("binary mass 1: %f\n", source->m1);
	printf("binary mass 2: %f\n", source->m2);
	printf("time of arrival: %f\n", source->time_of_arrival);
}

void Load_Source(inspiral_source_parameters_t* source) {
	source->right_ascension = -2.14;
	source->declination = 0.72;
	source->polarization_angle = M_PI / 6.0;
	source->coalesce_phase = 0.0;
	source->inclination = 0.0;
	source->m1 = 1.4 * GSL_CONST_MKSA_SOLAR_MASS; // binary mass 1
	source->m2 = 4.6 * GSL_CONST_MKSA_SOLAR_MASS; // binary mass 2
	source->time_of_arrival = 10.0;
}

void Print_Chirp_Factors(inspiral_chirp_factors_t* f) {
	printf("total mass: %f\n", f->total_mass);
	printf("reduced mass: %f\n", f->reduced_mass);
	printf("chirp mass: %f\n", f->chirp_mass);
	printf("amp factor 1: %f\n", f->amp_fact_1);
	printf("amp factor 2: %f\n", f->amp_fact_2);
	printf("s mass ratio: %f\n", f->s_mass_ratio);
	printf("multi fac: %f\n", f->multi_fac);

	printf("chirp time 0: %f\n", f->chirp_time0);
	printf("chirp time 1_5: %f\n", f->chirp_time1_5);
	printf("calculated reduced mass: %f\n", f->calculated_reduced_mass);
	printf("calculated total mass: %f\n", f->calculated_total_mass);
	printf("s mass ratio cal: %f\n", f->s_mass_ratio_cal);
	printf("multi fac cal: %f\n", f->multi_fac_cal);
	printf("chirp time 1: %f\n", f->chirp_time1);
	printf("chirp time 2: %f\n", f->chirp_time2);
	printf("t chirp: %f\n", f->t_chirp);
	printf("tc: %f\n", f->tc);

	printf("len time series chirp: %lu\n", f->len_tchirp);
	printf("time series chirp: \n");
	for (size_t i = 0; i < f->len_tchirp; i++) {
		printf(" %f,", f->tchirp[i]);
	}
	printf("\n");
}

void Compute_Chirp_Factors(double f_low, inspiral_source_parameters_t* source,
		inspiral_chirp_factors_t* fac) {

	const double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
	const double G = GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT;

	fac->total_mass = source->m1 + source->m2;
	fac->reduced_mass = source->m1 * source->m2 / (fac->total_mass);
	fac->chirp_mass = pow(fac->reduced_mass, 3.0 / 5.0)
			* pow(fac->total_mass, 2.0 / 5.0);

	fac->amp_fact_1 = 1.0;
	fac->amp_fact_2 = 1.0;

	// Change parameters accordingly to estimate chirp times.
	fac->s_mass_ratio = fac->reduced_mass / fac->total_mass;
	fac->multi_fac = f_low * M_PI * G * fac->total_mass / gsl_pow_3(c);

	fac->chirp_time0 = (5.0 / (256.0 * M_PI)) / f_low
			* pow(fac->multi_fac, -5.0 / 3.0) / fac->s_mass_ratio;

	fac->chirp_time1_5 = (1.0 / 8.0) / f_low * pow(fac->multi_fac, -2.0 / 3.0)
			/ fac->s_mass_ratio;

	//Calculated reduced masses using  chirp_time0 and chirp_time1_5
	fac->calculated_reduced_mass = (1.0 / (16.0 * gsl_pow_2(f_low)))
			* pow(
					5.0
							/ (4.0 * gsl_pow_4(M_PI) * fac->chirp_time0
									* gsl_pow_2(fac->chirp_time1_5)), 1.0 / 3.0)
			/ (G / gsl_pow_3(c));

	//Calculated total mass  using  chirp_time0 and chirp_time1_5
	fac->calculated_total_mass = (5.0 / (32.0 * f_low))
			* (fac->chirp_time1_5 / (gsl_pow_2(M_PI) * fac->chirp_time0))
			/ (G / gsl_pow_3(c));

	//Calculating chirp_time1 and chrip_time2 using calculated reduced_mass
	//and calculated_total_mass.
	fac->s_mass_ratio_cal = fac->calculated_reduced_mass
			/ fac->calculated_total_mass;
	fac->multi_fac_cal = f_low * M_PI * G * fac->calculated_total_mass
			/ gsl_pow_3(c);

	fac->chirp_time1 = (5.0 / (192.0 * M_PI)) / f_low / fac->multi_fac_cal
			/ fac->s_mass_ratio_cal
			* ((743.0 / 336.0) + ((11.0 / 4.0) * fac->s_mass_ratio_cal));
	fac->chirp_time2 = (5.0 / (128.0 * M_PI)) / f_low
			* pow(fac->multi_fac_cal, -1.0 / 3.0) / fac->s_mass_ratio_cal
			* ((3058673.0 / 1016064.0)
					+ ((5429.0 / 1008.0) * fac->s_mass_ratio_cal)
					+ ((617.0 / 144.0) * gsl_pow_2(fac->s_mass_ratio_cal)));

	fac->t_chirp = fac->chirp_time0 + fac->chirp_time1 - fac->chirp_time1_5
			+ fac->chirp_time2;

	fac->tc = source->time_of_arrival + fac->t_chirp;

	double dt = 0.01;
	fac->len_tchirp = ceil((fac->tc - source->time_of_arrival) / dt);
	fac->tchirp = (double*) malloc(fac->len_tchirp); //time_of_arrival:0.01:tc; //Time vector for the chirp.
	for (size_t i = 0; i < fac->len_tchirp; i++) {
		fac->tchirp[i] = source->time_of_arrival + i * dt;
	}
}

void Free_Chirp_Factors(inspiral_chirp_factors_t* fac) {
	free(fac->tchirp);
	fac->len_tchirp = 0;
}

void Compute_Interpolated_Strain(strain_t* strain, strain_t* interpolated) {
	gsl_interp* interp = gsl_interp_alloc(gsl_interp_linear, strain->len);
	gsl_interp_init(interp, strain->freq, strain->strain, strain->len);
	gsl_interp_accel* acc = gsl_interp_accel_alloc();

	int N = 65537;
	Alloc_Strain(N, interpolated);

	int samplingFrequency = 2048;
	for (int i = 0; i < N; i++) {
		// The frequency we want to have a strain value for
		double f = i * (0.5 * samplingFrequency / N);

		// The interpolated strain value
		double s = gsl_interp_eval(interp, strain->freq, strain->strain, f,
				acc);

		// store
		interpolated->freq[i] = f;
		interpolated->strain[i] = s;
	}

	gsl_interp_accel_free(acc);
	gsl_interp_free(interp);
}

void Save_Strain(char* filename, strain_t* strain) {
	FILE* file;
	file = fopen(filename, "w");
	if (file) {
		for (size_t i = 0; i < strain->len; i++) {
			fprintf(file, "%e\t %e\n", strain->freq[i], strain->strain[i]);
		}
		fclose(file);
	} else {
		printf("Error: Unable to save strain data to file (%s).\n", filename);
	}
}

int Read_Num_Strain_Samples(char* filename) {
	FILE* file;
	file = fopen(filename, "r");
	size_t len = 0;

	if (file) {
		double freq, strain;
		while (fscanf(file, "%lf %lf", &freq, &strain) != EOF) {
			len++;
			//printf("%e %e\n", freq, strain);
		}
		//printf("%lu lines\n", len);
		fclose(file);

		return len;
		// now store the val
	} else {
		printf("Error: Unable to open strain file (%s) for reading.\n",
				filename);
		return -1;
	}
}

void Print_Strain(strain_t* strain) {
	printf("Strain Curve: %ld samples\n", strain->len);
	for (size_t i = 0; i < strain->len; i++) {
		printf("%e %e\n", strain->freq[i], strain->strain[i]);
	}
}

int Load_Strain(char* filename, strain_t* strain) {
	FILE* file;

	Alloc_Strain(Read_Num_Strain_Samples(filename), strain);

	file = fopen(filename, "r");
	if (file) {
		int i = 0;
		while (fscanf(file, "%lf %lf", &strain->freq[i], &strain->strain[i])
				!= EOF) {
			i++;
		}
		fclose(file);
	} else {
		printf("Error: Unable to read the strain file into memory.\n");
		return -1;
	}

	return 0;
}

void Alloc_Strain(size_t len, strain_t* strain) {
	strain->len = len;
	strain->freq = (double*) malloc(strain->len * sizeof(double));
	strain->strain = (double*) malloc(strain->len * sizeof(double));
}

void Free_Strain(strain_t* strain) {
	free(strain->freq);
	free(strain->strain);
}

/*
 void Load_Frequencies() {
 typedef struct {

 } sampling_system_t;


 const double T = 64.0;
 const double sampling_frequency = 2048.0;
 const delta_t = 1.0 / sampling_frequency;
 t=(0:(floor(T/delta_t)-1))*delta_t;
 const int nSamples = length(t);
 T = nSamples / sampling_frequency;
 const int j_nyq = floor(nSamples / 2);
 f=(0:j_nyq)/T;
 */

typedef struct {
	size_t len;
	double* tempval;
	gsl_complex* spa_2pn;
	gsl_complex* whitened_sf;
	gsl_complex* h_0;
	gsl_complex* h_90;
	gsl_complex* whitened_signal;
	gsl_complex* whitened_data;
} stationary_phase_t;

double Compute_G(double f_low, double f_high, inspiral_chirp_factors_t* chirp, strain_t* interp_strain) {
	// whitening normalization factor
		double *tempval = (double*) malloc(interp_strain->len * sizeof(double));
		for (size_t i = 0; i < interp_strain->len; i++) {
			double f = interp_strain->freq[i];
			double s = interp_strain->strain[i];

			if (f > f_low && f < f_high) {
				tempval[i] = pow(f, -7.0 / 3.0) / gsl_pow_2(s);
			}
		}

		double sum_tempval = 0.0;
		for (size_t i = 0; i < interp_strain->len; i++) {
			sum_tempval += tempval[i];
		}

		free(tempval);

		double g = sqrt(gsl_pow_2(chirp->amp_fact_1) * gsl_pow_2(chirp->amp_fact_2) * sum_tempval);

		return g;
}

void Alloc_Stationary_Phase(size_t len, stationary_phase_t* sp) {
	sp->len = len;
	sp->spa_2pn = (gsl_complex*) malloc(len * sizeof(gsl_complex));
	sp->whitened_sf = (gsl_complex*) malloc(len * sizeof(gsl_complex));
	sp->h_0 = (gsl_complex*) malloc(len * sizeof(gsl_complex));
	sp->h_90 = (gsl_complex*) malloc(len * sizeof(gsl_complex));
	sp->whitened_signal = (gsl_complex*) malloc( len * sizeof(gsl_complex) );
	sp->whitened_data = (gsl_complex*) malloc( len * sizeof(gsl_complex) );
}

void Free_Stationary_Phase(stationary_phase_t* sp) {
	free(sp->spa_2pn);
	free(sp->whitened_sf);
	free(sp->h_0);
	free(sp->h_90);
	free(sp->whitened_signal);
	free(sp->whitened_data);
}

void Stationary_Phase(inspiral_source_parameters_t* source,
		inspiral_chirp_factors_t* chirp, strain_t* interp_strain,
		detector_t* det, double f_low, double f_high,
		stationary_phase_t* sp)
{
	double g = Compute_G(f_low, f_high, chirp, interp_strain);
	printf("g = %e\n", g);

		for (size_t j = 0; j < interp_strain->len; j++) {
			double f = interp_strain->freq[j];
			double s = interp_strain->strain[j];

			if (f > f_low && f < f_high) {
				double amp_2pn = (1.0 / g) * chirp->amp_fact_1
						* chirp->amp_fact_2 * pow(f, -7.0 / 6.0);
				double f_fac = f / f_low;
				double phase_2pn =
						2.0 * M_PI * f * (chirp->tc - det->timedelay)
								- 2.0 * source->coalesce_phase - M_PI / 4.0
								+ (2.0 * M_PI * f_low)
										* ((3.0 * chirp->chirp_time0
												* (pow(f_fac, -5.0 / 3.0)) / 5.0)
												+ (chirp->chirp_time1
														* (pow(f_fac, -5.0 / 3.0)))
												- (3.0 * chirp->chirp_time1_5
														* (pow(f_fac,
																-2.0 / 3.0))
														/ 2.0)
												+ (3.0 * chirp->chirp_time2
														* (pow(f_fac,
																-1.0 / 3.0))));

				gsl_complex exp_phase = gsl_complex_exp(gsl_complex_rect(0.0, -phase_2pn));
				sp->spa_2pn[j] = gsl_complex_mul_real(exp_phase, amp_2pn);
			} else {
				sp->spa_2pn[j] = gsl_complex_rect(0.0, 0.0);
			}

			sp->whitened_sf[j] = gsl_complex_div_real(sp->spa_2pn[j], s);
			sp->h_0[j] = sp->whitened_sf[j];

			gsl_complex v = gsl_complex_rect(0.0, -1.0);
			sp->h_90[j] = gsl_complex_mul(sp->h_0[j], v);

			double multi_factor = (1.0 / 2.8580) * 9.0;
			gsl_complex A = gsl_complex_mul_real( sp->h_0[j], det->f_plus );
			gsl_complex B = gsl_complex_mul_real( sp->h_90[j], det->f_cross );
			gsl_complex C = gsl_complex_add( A, B );
			sp->whitened_signal[j] = gsl_complex_mul_real(C, multi_factor );
		}
		/*
		 whitened_SF= (SPA_2pn)./ interp_strain;
		 h_0 = whitened_SF; //(Here h_0 should not have any antenna pattern dependency)
		 h_90 = -1i * h_0;

		 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		 //// Here E(noise_f.* conj(noise_f))= 1 ; For non whitened case it should be E(noise_f.* conj(noise_f))= S(f)
		 //// This is two sided spectral density since I am concatanating bothsides(full frequency range)
		 //// Matched filter noise output level is consistent with other existing pipelines (PYCBC,GSTLAL)
		 // noise_f=(1/(sqrt(2)))*randn(1,(j_nyq+1))+(1/(sqrt(2)))*1i*randn(1,(j_nyq+1));
		 //// Edited E(noise_f.* conj(noise_f))= 1/2 (It should be one sided)
		 noise_f = (1 / (sqrt(4))) * randn(1, (j_nyq + 1))
		 + (1 / (sqrt(4))) * 1i * randn(1, (j_nyq + 1));
		 ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		 //noise_f = 0 ;
		 //noise_f=randn(1,(j_nyq+1))+1i*randn(1,(j_nyq+1));
		 //noise_f = 0 ;
		 multi_factor = 1;
		 multi_factor = (1 / 2.8580) * 9;
		 multi_factor = 0;
		 whitened_signal = multi_factor
		 * ((h_0 * F_Plus_vec(id)) + (h_90 * F_Cross_vec(id)));
		 whitened_data
		 {
		 1, id
		 }
		 = whitened_signal + (noise_f);  //(whitened signal+whitened noise)

		 */
}

void Save_Stationary_Phase_Signal(char* filename, stationary_phase_t* sp) {
	FILE* file;
	file = fopen(filename, "w");
	for (size_t i = 0; i < sp->len; i++) {
		fprintf(file, "%e\n", GSL_REAL(sp->whitened_signal[i]));
	}
	fclose(file);
}

void DataGen() {
	const double f_low = 40.0; // seismic cutoff. All freqs low set to zero
	const double f_high = 700.0; // most stable inner orbit (last stable orbit related)

	inspiral_source_parameters_t source;
	printf("Inspiral source parameters:\n");
	Load_Source(&source);
	Print_Source(&source);

	printf("\n\nChirp Factors:\n");
	inspiral_chirp_factors_t chirp;
	Compute_Chirp_Factors(f_low, &source, &chirp);
	Print_Chirp_Factors(&chirp);

	detector_network_t net;
	Init_Detector_Network(&net);
	Compute_Detector_Network_Antenna_Patterns(&source, &net);
	Print_Detector_Network(&net);

	strain_t strain;
	Load_Strain("strain.txt", &strain);
	Print_Strain(&strain);

	strain_t interp_strain;
	Compute_Interpolated_Strain(&strain, &interp_strain);
	Print_Strain(&interp_strain);
	Save_Strain("interp.txt", &interp_strain);

	// For each detector determine the stationary phase of the signal
	for (size_t i = 0; i < net.num_detectors; i++) {
		detector_t* det = &net.detector[i];

		stationary_phase_t sp;

		Alloc_Stationary_Phase(interp_strain.len, &sp);

		Stationary_Phase(&source,
				&chirp, &interp_strain,
				det, f_low, f_high,
				&sp);

		Save_Stationary_Phase_Signal(det->id, &sp);

		Free_Stationary_Phase(&sp);


	}

	Free_Detector_Network(&net);
	Free_Chirp_Factors(&chirp);
	Free_Strain(&strain);
	Free_Strain(&interp_strain);

	const double oneMpc = 3.08567758e22; // 1 Mpc

}


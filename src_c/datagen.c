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

//// 07-October-2015 : Network Analysis (Signal generation for network of detectors)
//  Shihan Weerathunga

//  Ver 4.1
// This script generates gravitational wave signals in different detectors (in frequency domain),
// phase shifting according to time delays  and add all of them. 

//// *02-February-2016: Modification to add 2PN Signals*  

typedef struct {
	double right_ascension;
	double declination;
	double polarization_angle;
	double coalesce_phase;
	double inclination;
	double m1; // binary mass 1
	double m2; // binary mass 2
	double time_of_arrival;
} Signal;

void datagen() {

// Universal Constants //
//Velocity of Light(c) ; Gravitational Constant(G) ; Solar mass
const double c = GSL_CONST_MKSA_SPEED_OF_LIGHT;
const double G = GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT;
const double solar_mass = GSL_CONST_MKSA_SOLAR_MASS; //1.9891e30;

//  Signal Paramaeters //
const Signal s = {-2.14, 0.72, M_M_PI/6.0, 0.0, 0.0, 1.4, 4.6, 10.0};

const double oneMpc = 3.08567758e22; // 1 Mpc


// m1 = binary mass 1;    //m2 = binary mass 2;
     //m1=1.4;              m2=4.6;
// Conversion of binary masses to the solar mass.
     s.m1 *= solar_mass;
     s.m2 *= solar_mass;
  
//Arrival time of the signal in frequency band of the detector.
    //time_of_arrival=10;     //

    //// Time And Frequency Vectors
// time and frequency parameters.
// f_low          (Lower cut off frequency)
// fc          (Upper cut off frequency)
// T           (Total Time (Duration of the time vector))
// sampling_frequency
// delta_t     (Sampling Interval)
// t           (Time Vector (This is spaced by delta_t))
// nSamples    (Number of Samples)
// j_nyq       (Nyquist Frequency Index)
// f           (Frequency vector up to Nyquist frequency)

    const double f_low = 40.0;
    const double f_high = 700.0;
    const double T = 64.0;
    const double sampling_frequency = 2048.0;
    const delta_t = 1.0 / sampling_frequency;
    t=(0:(floor(T/delta_t)-1))*delta_t; 
    const int nSamples = length(t);
    T = nSamples/sampling_frequency; 
    const int j_nyq = floor(nSamples/2);
    f=(0:j_nyq)/T;

    ////  (LIGO Strain Values)
//   Load LIGO strain values(initial LIGO sensitivity curve) and interpolate 
//   it according to defined defined frequency vector.

    str=xlsread('strain.xlsx');
    frequency=str(:,1);strain=str(:,2);//square root value

    interp_fre =f;
    interp_strain=interp1(frequency,strain,interp_fre);

//Set interp_strain to a non-zero value below f_a
    indxfa = find(f<f_low);
    interp_strain(indxfa) = interp_strain(indxfa(end));

//Similarly for f > fc
    indxfc = find(f> f_high);
    interp_strain(indxfc) = interp_strain(indxfc(1));
////  Calculate chirp mass and chirp time for newtonian signal.
// Masses of two compact objects.
// total mass      (m1+m2)
// redu_mass       (m1*m2/(total_mass))
// chirp_mass
// chirp_time      (Duration of the chirp)
// tc              (Time of coalescence)(time of arrival + duration of the chirp)
// tchirp         (Time vector for the chirp)(time_of_arrival:0.01:tc)


    const double total_mass = s.m1 + s.m2;
    const double reduced_mass = s.m1*s.m2/(total_mass);
    const double chirp_mass = pow(reduced_mass, 3.0/5.0) * pow(total_mass, 2.0/5.0) ;
    
    //Amp_fact1 = sqrt(5/24)*(1/(M_PI^(2/3)*oneMpc));
    //Amp_fact2 = chirp_mass^(-5/3);
    
    const double Amp_fact1 = 1.0;
    const double Amp_fact2 = 1.0;
    
    //chirp_time = 34.5*((chirp_mass/solar_mass)^(-5/3))*((f_low/40)^(-8/3));
    //total_mass = m1 + m2 ; //Total mass
    //reduced_mass = (m1*m2)/total_mass;
    
    // Change parameters accordingly to estimate chirp times.
    
    const double s_mass_ratio = reduced_mass/total_mass;
    const double multi_fac = f_low * M_PI * G * total_mass/(c^3);
 
    //multi_fac = (5 / (256*s_mass_ratio))* (G*total_mass/(c^3));
    //v_low = (G*total_mass * M_PI * f_low / (c^3))^(1/3);
    
    // For 1.4xSolar_mass binaries
    // chirp_time0 = 24.8499; chirp_time1 = 1.3861; 
    // chirp_time1_5 = 0.8662; chirp_time2 = 0.0479; 
    
    chirp_time0 = (5/(256*M_PI)) * f_low^-1 * multi_fac^(-5/3) * s_mass_ratio^-1;
    
    //chirp_time1 = (5/(192*M_PI)) * f_low^-1 * multi_fac^(-1) * s_mass_ratio^-1 * ...
    //                                              ((743/336)+((11/4)*s_mass_ratio));
                                              
    chirp_time1_5 = (1/8) * f_low^-1 * multi_fac^(-2/3) * s_mass_ratio^-1;
    
    //chirp_time2 =   (5/(128*M_PI)) * f_low^-1 * multi_fac^(-1/3) * s_mass_ratio^-1 * ...
    //            ((3058673/1016064)+((5429/1008)*s_mass_ratio)+ ((617/144)*(s_mass_ratio^2)));
 
    //Calculated reduced masses using  chirp_time0 and chirp_time1_5
  calculated_reduced_mass = (1/(16*f_low^2))*(5/(4*M_PI^4*chirp_time0*chirp_time1_5^2))^(1/3)...
                                                                              *(G/c^3)^-1;
   //Calculated total mass  using  chirp_time0 and chirp_time1_5                                                                       
  calculated_total_mass = (5/(32*f_low))*(chirp_time1_5/(M_PI^2*chirp_time0)) * (G/c^3)^-1;
    
   
   //Calculating chirp_time1 and chrip_time2 using calculated reduced_mass
   //and calculated_total_mass.
   
   s_mass_ratio_cal = calculated_reduced_mass/calculated_total_mass; 
   multi_fac_cal = f_low * M_PI * G * calculated_total_mass/(c^3);
   
   chirp_time1 = (5/(192*M_PI)) * f_low^-1 * multi_fac_cal^(-1) * s_mass_ratio_cal^-1 * ...
                                                  ((743/336)+((11/4)*s_mass_ratio_cal));
   chirp_time2 =   (5/(128*M_PI)) * f_low^-1 * multi_fac_cal^(-1/3) * s_mass_ratio_cal^-1 * ...
                ((3058673/1016064)+((5429/1008)*s_mass_ratio_cal)+ ((617/144)*(s_mass_ratio_cal^2)));                                          
                                              
   
    T_chirp = chirp_time0 + chirp_time1 - chirp_time1_5 + chirp_time2 ;
    //tc=time_of_arrival+chirp_time;
    
    tc = time_of_arrival+T_chirp;
    tchirp=time_of_arrival:0.01:tc; //Time vector for the chirp.

//fprintf('Phase 1, done');
////  Calculation of  normalization factor g for PN signals.
// If I use different spectral density values here, I will get different 
// normalization constants according to detectors.

    tempval = ones(1,length(f));
    for j=1:1:length(f)
    
        if (f_low< f(j)&& f(j)<=f_high)
            tempval(j)=(interp_fre(j)^(-7/3))/(interp_strain(j)^2);
        end
    end
    //g=sqrt((8/3)*sum(tempval)); //(2 Re convention)
    g = sqrt((Amp_fact1^2)*(Amp_fact2^2)*sum(tempval));

//// Antenna pattern  calculation for the know signal parameters.

    detId = cellstr(['H1';'L1';'V1';'K1']);
    //detId = cellstr(['H1';'L1']);

// data generation
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

U_vec = ones(1,length(detId)); V_vec = ones(1,length(detId)); 
F_Plus_vec = ones(1,length(detId));
F_Cross_vec = ones(1,length(detId)); 


SNR_d = ones(1,length(detId));

for id = 1:1:length(detId)    
    [U_vec(id),V_vec(id),F_Plus_vec(id),F_Cross_vec(id)]= antennapattern(declination,right_ascension ,polarization_angle,detId(id));
end
// //// For checking, assigning the same strength for all the antenna patterns
//  U_vec(1)=1;    U_vec(2)=1; U_vec(3)=1;    U_vec(4)=1;
//  V_vec(1)=1;    V_vec(2)=1; V_vec(3)=1;    V_vec(4)=1; 
//  F_Plus_vec(1)=1;  F_Plus_vec(2)=1; F_Plus_vec(3)=1;  F_Plus_vec(4)=1;
//  F_Cross_vec(1)=1; F_Cross_vec(2)=1;F_Cross_vec(3)=1; F_Cross_vec(4)=1;
// 
// //// Comment above for different antenna patterns
 
for id = 1:1:length(detId)
    
    delay(id)=timedelay(detId(id),declination, right_ascension);
    
    amp_2pn = ones(1,length(f));
    phase_2pn = ones(1,length(f));
    SPA_2pn = ones(1,length(f));
    
    //// 2PN  Waveform
    for jj=1:1:length(f)
        
        if (f_low< f(jj)&& f(jj)<=f_high)
            
            //amp_ge(jj) = (2/g) * sqrt(2/(3*f_low)) *((f(jj)/f_low)^(-7/6));
            //v(jj) = (G*total_mass * M_PI * f(jj) / (c^3))^(1/3);
            amp_2pn(jj)=  (1/g)*Amp_fact1*Amp_fact2*(f(jj)^(-7/6));
            f_fac(jj) = f(jj)/f_low;
            
            phase_2pn(jj) = 2*M_PI*f(jj)*(tc - delay(id)) - 2*coalesce_phase - M_PI/4 + ...
                            + (2*M_PI*f_low) * (...
                                             ( 3*chirp_time0 *( f_fac(jj)^(-5/3) )/5  )...
                                            +( chirp_time1*( f_fac(jj)^(-5/3) ) )...
                                            -( 3*chirp_time1_5* ( f_fac(jj)^(-2/3) )/2)...
                                            +( 3*chirp_time2* ( f_fac(jj)^(-1/3) ) )...
                                         );
            
            SPA_2pn(jj) = amp_2pn(jj)* exp(-1i*phase_2pn(jj));
         else
            SPA_2pn(jj)=0;
        end
     end
       
    whitened_SF= (SPA_2pn)./ interp_strain;
    h_0 = whitened_SF;  //(Here h_0 should not have any antenna pattern dependency)
    h_90 = -1i*h_0;
   
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //// Here E(noise_f.* conj(noise_f))= 1 ; For non whitened case it should be E(noise_f.* conj(noise_f))= S(f)
    //// This is two sided spectral density since I am concatanating bothsides(full frequency range)
    //// Matched filter noise output level is consistent with other existing pipelines (PYCBC,GSTLAL)    
     // noise_f=(1/(sqrt(2)))*randn(1,(j_nyq+1))+(1/(sqrt(2)))*1i*randn(1,(j_nyq+1));
    //// Edited E(noise_f.* conj(noise_f))= 1/2 (It should be one sided)
       noise_f=(1/(sqrt(4)))*randn(1,(j_nyq+1))+(1/(sqrt(4)))*1i*randn(1,(j_nyq+1));
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //noise_f = 0 ;
    //noise_f=randn(1,(j_nyq+1))+1i*randn(1,(j_nyq+1));
    //noise_f = 0 ;
    multi_factor = 1;
    multi_factor = (1/2.8580)*9;
    multi_factor = 0;
     whitened_signal=multi_factor*((h_0 * F_Plus_vec(id))+ (h_90 * F_Cross_vec(id)));
    whitened_data{1,id}= whitened_signal + (noise_f);  //(whitened signal+whitened noise)
   

end
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//// whitened_signal  + white noise = whitened_data is already generated for each detector at this moment.
//// whitened signal at the detector is time shifted accroding to the detector locations relative to the fide detector   ////////
     
 
sec1param.detId=detId ;
sec1param.g = g;
sec1param.f = f;
sec1param.f_low = f_low;
sec1param.f_high = f_high;
sec1param.interp_strain = interp_strain;
sec1param.whitened_data = whitened_data;   
sec1param.nSamples = nSamples;
sec1param.polarization_angle=polarization_angle;

}

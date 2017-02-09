 
%% Generating 2 nd postnewtonian signal templates for correlation operations in frequency domain.
%% Signal phase is 2nd postnewtonian.

function [SPA_0 ,SPA_90] = templategen(Id,chirp_time0_input,chirp_time1_5_input,...
                                                    rightascension_input,declination_input,g,f,f_low,f_high)
    
    %Velocity of Light(c) ; Gravitational Constant(G) ; Solar mass
    c=299792458;            G=6.67408e-11;              
    
    %solar_mass = 1.9891e30; %polarization_angle  = 0;   %solar_mass = 1.9891e30;
    
    chirp_time0 = chirp_time0_input;
    chirp_time1_5 = chirp_time1_5_input;
    declination = declination_input;         %Declination Angle
    right_ascension = rightascension_input;  %Right Ascension Angle
    
    
    % Time delay calculations relative to earth centered detector.
    delay=timedelay(Id, declination, right_ascension);

    
   Amp_fact1=1;     Amp_fact2 = 1;    
  %Calculated reduced masses using  chirp_time0 and chirp_time1_5
  calculated_reduced_mass = (1/(16*f_low^2))*(5/(4*pi^4*chirp_time0*chirp_time1_5^2))^(1/3)...
                                                                              *(G/c^3)^-1;
  %Calculated total mass  using  chirp_time0 and chirp_time1_5                                                                       
  calculated_total_mass = (5/(32*f_low))*(chirp_time1_5/(pi^2*chirp_time0)) * (G/c^3)^-1;
    
  s_mass_ratio_cal = calculated_reduced_mass/calculated_total_mass; 
  multi_fac_cal = f_low * pi * G * calculated_total_mass/(c^3);
   
  chirp_time1 = (5/(192*pi)) * f_low^-1 * multi_fac_cal^(-1) * s_mass_ratio_cal^-1 * ...
                                                  ((743/336)+((11/4)*s_mass_ratio_cal));
  chirp_time2 =   (5/(128*pi)) * f_low^-1 * multi_fac_cal^(-1/3) * s_mass_ratio_cal^-1 * ...
                ((3058673/1016064)+((5429/1008)*s_mass_ratio_cal)+ ((617/144)*(s_mass_ratio_cal^2))); 
    
    
  %Total Chirp time (Duration of the chirp signal)
  T_chirp = chirp_time0 + chirp_time1 -  chirp_time1_5 + chirp_time2 ;
    
   
%% Template Generation.
  time_of_arrival=0;                  % Time of arrival = 0 for templates
  tc_template=time_of_arrival+T_chirp;        %coalescence time.
  coalesce_phase = 0;

% Pre allocation of memory for faster computing.
   amp_SPA = ones(1,length(f));     phase_SPA = ones(1,length(f));
   SPA_0 = ones(1,length(f));       SPA_90= ones(1,length(f));
    
   
    for jj=1:1:length(f)
        
        if (f_low< f(jj)&& f(jj)<=f_high)
            %%  h_0 template.
            
           
            amp_SPA(jj)=  (1/g)*Amp_fact1*Amp_fact2*(f(jj)^(-7/6));
            f_fac(jj) = f(jj)/f_low;
            
            phase_SPA(jj) = 2*pi*f(jj)*(tc_template - delay) - 2*coalesce_phase - pi/4 + ...
                            + (2*pi*f_low) * (...
                                             ( 3*chirp_time0 *( f_fac(jj)^(-5/3) )/5  )...
                                            +( chirp_time1*( f_fac(jj)^(-5/3) ) )...
                                            -( 3*chirp_time1_5* ( f_fac(jj)^(-2/3) )/2)...
                                            +( 3*chirp_time2* ( f_fac(jj)^(-1/3) ) )...
                                         );
                                    
                                                     
%            
            % gsl_complex gsl_complex_exp (gsl_complex z)
             SPA_0(jj) = amp_SPA(jj)* exp(-1i*phase_SPA(jj));
            
            
            %%  h_pi/2 template.

            SPA_90(jj) =  -1i*SPA_0(jj);
        else
            SPA_0(jj)=0;
            SPA_90(jj)=0;
        end
    end
    
    

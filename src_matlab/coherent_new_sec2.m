%%% This section will be evaluated using PSO evaluated variables.
function [out_val] = coherent_new_sec2(sec1param,chirp_time0_input,chirp_time1_5_input,...
                                      rightascension_input,declination_input)


detId   = sec1param.detId;
g       = sec1param.g;
f       = sec1param.f;
f_low     =  sec1param.f_low;
f_high     =  sec1param.f_high;
interp_strain = sec1param.interp_strain;
whitened_data = sec1param.whitened_data;     
nSamples = sec1param.nSamples;
polarization_angle_input = 0;

 for id = 1:1:length(detId)
    [U_vec_input(id),V_vec_input(id),F_Plus_vec_input(id),F_Cross_vec_input(id)]...
        = antennapattern(declination_input,rightascension_input ,polarization_angle_input,detId(id));   
 end
% %% For checking, assigning the same strength for all the antenna patterns
%  U_vec_input(1)=1;    U_vec_input(2)=1; U_vec_input(3)=1;    U_vec_input(4)=1;
%  V_vec_input(1)=1;    V_vec_input(2)=1; V_vec_input(3)=1;    V_vec_input(4)=1; 
%  F_Plus_vec_input(1)=1;  F_Plus_vec_input(2)=1; F_Plus_vec_input(3)=1;  F_Plus_vec_input(4)=1;
%  F_Cross_vec_input(1)=1; F_Cross_vec_input(2)=1;F_Cross_vec_input(3)=1; F_Cross_vec_input(4)=1;
% 
% %% Comment above for different antenna patterns
 
UdotU_input=  dot(U_vec_input,U_vec_input);
UdotV_input=  dot(U_vec_input,V_vec_input);
VdotV_input=  dot(V_vec_input,V_vec_input);
 
A_input = UdotU_input; B_input = UdotV_input; C_input = VdotV_input;

%%  Network vector non-zero elements

M(1,1)= A_input ;         M(1,2)= B_input ;
M(2,1)= B_input ;         M(2,2)= C_input ;
M(3,3)= A_input ;         M(3,4)= B_input ;
M(4,3)= B_input;          M(4,4)= C_input;

Delta_input = (A_input*C_input) - (B_input*B_input);
Delta_factor_input = 1 / sqrt(2*Delta_input) ; 

D_input = sqrt(((A_input-C_input)^2) + 4*(B_input^2)) ;

P1_input = (C_input-A_input-D_input);
P2_input = (C_input-A_input+D_input);
P3_input = sqrt(C_input+A_input+D_input);
P4_input = sqrt(C_input+A_input-D_input);

G1_input =  sqrt((P1_input^2) +  4*(B_input^2))/ (2*B_input) ;
G2_input =  sqrt((P2_input^2) +  4*(B_input^2))/ (2*B_input) ;


O11_input = Delta_factor_input * P3_input /G1_input ;
O12_input = Delta_factor_input * P3_input * P1_input / (2*B_input*G1_input) ;
O21_input = Delta_factor_input * P4_input /G2_input ;
O22_input  = Delta_factor_input * P4_input * P2_input / (2*B_input*G2_input) ;

 
 for id = 1:1:length(detId)
           
            w_plus_input(id) = (O11_input*U_vec_input(id) +  O12_input*V_vec_input(id));
            w_minus_input(id)= (O21_input*U_vec_input(id) +  O22_input*V_vec_input(id));
            
            
           [SPA_0 ,SPA_90] = templategen(detId(id),chirp_time0_input,chirp_time1_5_input,...
                              rightascension_input,declination_input,g,f,f_low,f_high);
           %% Whitenning and Conjugation of Templates
            
            conjSF0=conj(SPA_0); % conj(h_c)
            conjSF90=conj(SPA_90);% conj(h_s)
            
            whitened_conj_SF0= (conjSF0)./ interp_strain;   % conj(h_c)/sqrt (S_n(f))
            whitened_conj_SF90= (conjSF90)./ interp_strain; % conj(h_s)/sqrt (S_n(f))
            
            %% Filtering Operation
            % matchedfilter;
              
            s0=((whitened_conj_SF0.*whitened_data{1,id})); %
            s00=conj(s0(end-1:-1:2));
            s0_s00=(horzcat(s0,s00));   %Prepare the full frequency vector with both positive and negative frequency components
            c_plus{1,id} =(s0_s00);
            
            s1=((whitened_conj_SF90.*whitened_data{1,id}));
            s11=conj(s1(end-1:-1:2));
            s1_s11=(horzcat(s1,s11));   %Prepare the full frequency vector with both positive and negative frequency components
            c_minus{1,id}= (s1_s11);
                   
            
 end

term_1(1,:) = (w_plus_input(1)*c_plus{1,1});
term_2(1,:) = (w_minus_input(1)*c_plus{1,1});
term_3(1,:) = (w_plus_input(1)*c_minus{1,1});
term_4(1,:) = (w_minus_input(1)*c_minus{1,1});

for id = 2:1:length(detId)
    term_1(1,:) = term_1(1,:)+ (w_plus_input(id)*c_plus{1,id});
    term_2(1,:) = term_2(1,:)+ (w_minus_input(id)*c_plus{1,id});
    term_3(1,:) = term_3(1,:)+ (w_plus_input(id)*c_minus{1,id});
    term_4(1,:) = term_4(1,:)+ (w_minus_input(id)*c_minus{1,id});
end

f1_tmp=(nSamples)*real(ifft(term_1(1,:)));
f2_tmp=(nSamples)*real(ifft(term_2(1,:)));
f3_tmp=(nSamples)*real(ifft(term_3(1,:)));
f4_tmp=(nSamples)*real(ifft(term_4(1,:)));
   
              
            
tmp_ifft = (f1_tmp.^2 + f2_tmp.^2+f3_tmp.^2+f4_tmp.^2);
%plot(sqrt(tmp_ifft));
out_val=max(sqrt(tmp_ifft));

end

%% Condition Number calculation for the detector network matrix.
%  Network matrix values depend on the sky locations. Well-coditioned and ill-conditioned
%  nature of the network matrix (value of the condtion number) gives and idea of invertability of 
%  the matrix according the sky location. 
%  Large condition number -> large numerical errors in estimations.


%% Shihan Weerathunga 
function cond_num = conditionnum_cal(detId,declination_input,...
                             rightascension_input,polarization_angle_input)


for id = 1:1:length(detId)
 [U_vec_input(id),V_vec_input(id),~,~] = antennapattern(declination_input,rightascension_input,...
                                                                  polarization_angle_input,detId(id));   
end


UdotU_input=  dot(U_vec_input,U_vec_input);
UdotV_input=  dot(U_vec_input,V_vec_input);
VdotV_input=  dot(V_vec_input,V_vec_input);
 
A_input = UdotU_input; B_input = UdotV_input; C_input = VdotV_input;

%%  Network vector non-zero elements
M(1,1)= A_input ;         M(1,2)= B_input ;
M(2,1)= B_input ;         M(2,2)= C_input ;
M(3,3)= A_input ;         M(3,4)= B_input ;
M(4,3)= B_input;          M(4,4)= C_input;


cond_num = cond(M);
end

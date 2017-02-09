% Generating gravitational wave antenna responses for different detectors according to 
% the sky location and selected basis tensors.

%% Shihan Weerathunga / September-8-2014
%% Antennapattern_New May-23-2016 (Adding extra information according to
%% Earth Centered, Earth Fixed coordinates in relation to latitude and longitude) 


% 5/23/2016 Gravitational Wave frame unit vectors are changed according to Drew Keppel/ John T. Whelan's papers.

function [u,v,F_Plus,F_Cross]= antennapattern(declination,right_ascention,polarization_angle, Intefe_ID)

% All LIGO coordinates are from LIGO-P000006-D-E:
if strcmp(Intefe_ID, 'H1') || strcmp(Intefe_ID, 'H2')
    % Generate coordinates for LHO:
    
    % X-arm unit vector relative to earth centered frame. 
    cordnt.X = [-0.223891216, 0.799830697, 0.556905359];
    % Y-arm unit vector relative to earth centered frame.
    cordnt.Y = [-0.913978490, 0.0260953206, -0.404922650]; 
    
    %Unit vectors in detector frame to generate antenna patterns. (Detector frame antenna patterns)
    %cordnt.X = [1, 0, 0];  
    %cordnt.Y = [0, 1, 0];  
    
elseif strcmp(Intefe_ID, 'L1')
    % Generate coordinates for LLO:
    cordnt.X = [-0.954574615, -0.141579994, -0.262187738]; 
    cordnt.Y = [ 0.297740169, -0.487910627, -0.820544948]; 
    % X-arm unit vector, L1-|cordnt.X| = 3995.15;
    % Y-arm unit vector, L1-|cordnt.Y| = 3995.15;

% Source : www.ligo.org/scientists/GW100916/detectors.txt 
elseif strcmp(Intefe_ID, 'V1') 
    cordnt.X = [-0.70045821479,  0.20848948619, 0.68256166277];
    cordnt.Y = [-0.05379255368, -0.96908180549, 0.24080451708];
    
% Source : www.ligo.org/scientists/GW100916/detectors.txt     
elseif strcmp(Intefe_ID, 'G1') 
    cordnt.X = [ -0.445184239,  0.866534205,  0.225675575];
    cordnt.Y = [ -0.626000687, -0.552167273,  0.550667271];
    
% Source : www.ligo.org/scientists/GW100916/detectors.txt   
elseif strcmp(Intefe_ID, 'K1') 
    cordnt.X = [ -0.4300,   -0.8363,  0.3400];
    cordnt.Y = [  0.6821,   -0.0542,  0.7292];
    
% Double check following coordinates. At the moment I am not using this detector
elseif strcmp(Intefe_ID, 'T1')
    cordnt.X = [ 0.648969405, 0.760814505, 0];
    cordnt.Y = [-0.443713769, 0.378484715, -0.812322234];
   
else
    error('getifo: The Intefe_ID you entered is not currently listed.')
end

%------------------------------------------------------------------------------

% Define X and Y unit vectors in gw frame. (exEarth)
% (This is found from rotating about z counterclockwise angle 
% right_ascention and then rotating counterclockwise about y angle 
% declination)

 
n_hat = [cos(right_ascention)*cos(declination),...
               sin(right_ascention)*cos(declination), (sin(declination))];

ex_i = [(sin(right_ascention)), -cos(right_ascention), 0]; 
ey_j = [-cos(right_ascention)*sin(declination),...
               -sin(right_ascention)*sin(declination), (cos(declination))];
   

% Define epsilon_plus and epsilon_cross. (Polarization basis tensors) 
%ePlusEarth = ex'*ex - ey'*ey   %eCrossEarth  = ex'*ey + ey'*ex

epsilon_plus = (ex_i'*ex_i)-(ey_j'*ey_j);
epsilon_cross = (ex_i'*ey_j)+(ey_j'*ex_i);


%%  Detector Tensor Calculations --------------------------------------------------
% Make sure that the unit vectors are in 1x3 form
Cordnt.U =cordnt.X(:)'; cordnt.V = cordnt.Y(:)';

% D = 0.5*(Xhat'*Xhat - Yhat'*Yhat) where Xhat is the unit vector for the
% 'X' arm and Yhat is the unit vector for the 'Y' arm
D = 0.5*((Cordnt.U'*Cordnt.U)-(cordnt.V'*cordnt.V));
%----------------------------------------------------------------------------------

% Tensor contraction.(Polarization Angel Independent) 
u = trace(D*epsilon_plus);   % u = D:epsilon_plus
v = trace(D*epsilon_cross);  % v = D:epsilon_cross
%----------------------------------------------------------------------------------

% wave polarizarion matrix(Polarization Angel Dependent) 
wp_matrix = [cos(2*polarization_angle) sin(2*polarization_angle);...
                    -sin(2*polarization_angle) cos(2*polarization_angle) ];
                
    F = wp_matrix * [u v]';
    
    F_Plus =F(1);     
    F_Cross=F(2);

end




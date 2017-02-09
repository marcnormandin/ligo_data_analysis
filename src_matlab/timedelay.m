%% 5/5/2014 
% Calculating the time delay of gravitational waves (with respect to earth centered detector) at different 
% LIGO Detectors according to the sky locations.

%  Shihan Weerathunga
%% Timedelay New : May-23-2016

function [delay] = timedelay_new(Intefe_ID,declination,right_ascension)
% global WGS; % WGS-84 coordinates of beam splitter

% One coordinate specified in the document for both H1 and H2)
if strcmp(Intefe_ID,'H1') || strcmp(Intefe_ID,'H2')
    
    WGS = [-2.161414928e+6, -3.834695183e+6, 4.600350224e+6]; % In meters
    
elseif strcmp(Intefe_ID, 'L1')
    % WGS-84 coordinates of 'L1' beam splitter
    WGS = [-7.427604192e+4, -5.496283721e+6, 3.224257016e+6];
    
elseif strcmp(Intefe_ID, 'V1')
    % WGS-84 coordinates of 'V1' beam splitter
    WGS = [4.5463741e+6, 8.429897e+5, 4.378577e+6];
    
elseif strcmp(Intefe_ID, 'G1')
    % WGS-84 coordinates of 'G1' beam splitter
    WGS = [3.8563112e+6,  6.665978e+5, 5.0196406e+6];
    
elseif strcmp(Intefe_ID, 'K1')
    % WGS-84 coordinates of 'K1' beam splitter
    WGS = [-3.776899062e+6,  3.483900163e+6, 3766657.585];
    
elseif strcmp(Intefe_ID, 'T1')
    % WGS-84 coordinates of 'T1' beam splitter
    WGS = [-3.946409e6, 3.366259e6, 3.6991507e6];
    
end
% -----------------------------------------------------------------------

% Define the X, Y & Z WGS-84 coordinates for Intefe_ID

xifo = WGS(1); % m
yifo = WGS(2); % m
zifo = WGS(3); % m


% Define the X, Y & Z coordinates of the GW's sky location

% Unit vector n_hat(right_ascension,declination) projections along  x,y,z.

%xgw = sin(declination)*cos(right_ascension);
%ygw = sin(declination)*sin(right_ascension);
%zgw = cos(declination);

xgw = cos(declination)*cos(right_ascension);
ygw = cos(declination)*sin(right_ascension);
zgw = sin(declination);

% xgw, ygw and zgw  are unit vector components. So modulus of the vector
% should be one. (For testing purposes)

%Define the GW's unit vector
gunit = sqrt((xgw^2)+(ygw^2)+(zgw^2));

% Define the speed of light
C=299792458; % m/s

%Double check the negative sigh infront.
% Shihan's expalnation: You don't have to put minus here. Plus or minus will be decided by the sky location
delay = (((xifo*xgw)+(yifo*ygw)+(zifo*zgw))/(C));

end

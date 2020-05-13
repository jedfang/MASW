% MASWaves Combination
% Version: 06.2018
%%
%   [c_mean,c_min_std,c_plus_std,lambda_mean,nPoints_fin],...
%     = MASWaves_combined_dispersion_curve(no_measurements,c_all,lambda_all,a)
%
%%
%   The function MASWaves_combined_dispersion_curve is used to obtain 
%   a composite experimental dispersion curve, along with upper and lower 
%   boundaries, by combining a number of elementary dispersion curve. 
%   (A dispersion curve obtained from a single multichannel surface wave 
%   record is here referred to as an elementary dispersion curve.)
%
%   The composite dispersion curve is obtained by adding the elementary
%   curves up within 1/a octave wavelength intervals. The upper and lower 
%   boundary curves correspond to plus/minus one standard deviation from 
%   the mean dispersion curve.
%
%   The combined mean dispersion curve can then be used as an input in 
%   the inversion analysis
% 
%   This tool is written based on the assumption that all the measured
%   elementary dispersion curve wavelength values are in the range of 
%   lambda = 0.1 m to lambda = 256 m. If this assumption is not valid, 
%   the lower/upper values for the measured wavelengths can be altered 
%   by changing the values of q_min and q_max in lines 103 and 104 below.
%
%% Input
%  no_measurements  Number of measurements (number of extracted elementary 
%                   dispersion curves)
%  lambda_all       A cell array of length no_measurements, where 
%                   - cell no. 1 contains the wavelength values of
%                     elementary dispersion curve no. 1
%                   - cell no. 2 contains the wavelength values of
%                     elementary dispersion curve no. 2
%                   - etc.
%  c_all            A cell array of length no_measurements, where 
%                   - cell no. 1 contains the Rayleigh wave phase velocity
%                     values of elementary dispersion curve no. 1
%                   - cell no. 2 contains the Rayleigh wave phase velocity
%                     values of elementary dispersion curve no. 2
%                   - etc.
%  a                Combination parameter. 
%                   The elementary dispersion curves stores in c_all and 
%                   lambda_all are added up within log_a distributed
%                   wavelength intervals.
%                   For selection of a-vales refer to:
%                   Ólafsdóttir, E.Á., Bessason, B. & Erlingsson, S. (2018). 
%                   Combination of dispersion curves from MASW measurements.
%                   Soil Dynamics and Earthquake Engineering, 113, 473-487.
%                   [Open access] 
%                   https://doi.org/10.1016/j.soildyn.2018.05.025
%
%% Output
%  Combined dispersion curve
%  c_mean           Rayleigh wave velocity [m/s]
%  lambda_mean      Wavelength [m]
%  c_plus_std       Upper bound Rayleigh wave velocity
%                   [m/s] (Mean value plus one standard deviation)
%  c_min_std        Lower bound Rayleigh wave velocity
%                   [m/s] (Mean value minus one standard deviation)
%
%  nPoints_fin      Number of wavelength intervals (number of data points
%                   included in the combined dispersion curve)
%
%% Subfunctions
%  (None)
%
%% 
function [c_mean,c_min_std,c_plus_std,lambda_mean,nPoints_fin] = MASWaves_combined_dispersion_curve(no_measurements,c_all,lambda_all,a)
%%
% Reformat lambda_all and c_all into
% lambda_vec: Vector containing the wavelength values of all the elementary
%             dispersion curves, i.e. 
%             lambda_vec = [lambda_all{1} lambda_all{2} ... lambda_all{no_measurements}]
% cdata_vec:  Vector containing the phase velocity values of all the  
%             elementary dispersion curves, i.e.
%             cdata_vec = [c_all{1} c_all{2} ... c_all{no_measurements}]

vec = 1:no_measurements;
temp = zeros(no_measurements,1);
for i = 1:no_measurements
    temp(i) = length(lambda_all{i,1});
end

lambda_data = zeros(length(vec),max(temp));
c_data = zeros(length(vec),max(temp));
k = 1;
for i=1:no_measurements
    lambda_data(k,1:length(lambda_all{i,1})) = lambda_all{i,1};
    c_data(k,1:length(c_all{i,1})) = c_all{i,1};
    k = k+1;
end
[lambda_size1,lambda_size2] = size(lambda_data);
lambda_vec = reshape(lambda_data',lambda_size1*lambda_size2,1);
c_vec = reshape(c_data',lambda_size1*lambda_size2,1);
lambda_zeros = find(lambda_vec==0);
lambda_vec(lambda_zeros)=[]; 
c_vec(lambda_zeros)=[]; 

% Define wavelengths intervals
q_min = -9; % The first wavelength interval starts at wavelength 2^(((q_min-1)/3)-1/6) m (approximately 0.1 m)
q_max = 50; % The last wavelength interval ends at wavelength 2^(((q_max-1)/3)+1/6) m (approx 287 m)

Third = 2.^(1./a);
Sixth = 2.^(1./(2*a));
m = 1;
lambda_low(m) = 1./Sixth; % Frequency lower bound 
lambda_mean(m) = 2^(((q_min-1)/a)); % Center 
lambda_up(m) = 1.*Sixth; % Frequency upper bound
for i = (q_min + 1):q_max
    m = m+1;
    lambda_mean(m) = lambda_mean(m-1)*Third;
    lambda_low(m) = lambda_mean(m)/Sixth;
    lambda_up(m) = lambda_mean(m)*Sixth;
end
%%
% Initialize
nPoints_initial = length(q_min:1:q_max);
nPoints_fin = 0; % Number of points (counter)
c_sum  = zeros(nPoints_initial,1);
c_mean = zeros(nPoints_initial,1);
c_sqrt  = zeros(nPoints_initial,1);
N_lambda = zeros(nPoints_initial,1);

% Find elementary dispersion curve data points that fall within each wavelength interval 
for i = 1:length(lambda_vec)
    for j = 1:length(lambda_mean)
        % Data point within wavelength interval j
        if lambda_low(j) <= lambda_vec(i) && lambda_vec(i) < lambda_up(j)
            c_sum(j) = c_sum(j) + c_vec(i);
            c_sqrt(j) = c_sqrt(j) + c_vec(i)^2;
            N_lambda(j) = N_lambda(j) + 1; % Number of data points within wavelength interval j
        end
    end
end

% Compute the mean phase velocity value and standard deviation for each
% wavelength interval
for i = 1:length(N_lambda)
    if N_lambda(i) > 1 
        c_mean(i) = c_sum(i)/N_lambda(i);
        VelStd(i) =  sqrt((c_sqrt(i)- N_lambda(i)*c_mean(i)^2)/(N_lambda(i)-1));
        nPoints_fin = nPoints_fin + 1; 
        c_min_std(i) = c_mean(i) - VelStd(i);
        c_plus_std(i) = c_mean(i) + VelStd(i);
    elseif  N_lambda(i) <= 1  
        c_mean(i) = 0;
        nPoints_fin = nPoints_fin + 1;
    end
end

% Remove zero values
ref = find(c_mean ~= 0);
c_mean = c_mean(ref);
c_plus_std = c_plus_std(ref);
c_min_std = c_min_std(ref);
lambda_mean = lambda_mean(ref);

end

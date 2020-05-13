%   MASWaves Dispersion
%   Version: 07.2017
%%
%  [f,c,A] = MASWaves_dispersion_imaging(u,n,x,fs,cT_min,cT_max,delta_cT)
%% 
%  The function MASWaves_dispersion_imaging carries out the first three  
%  steps of the dispersion analysis of the recorded surface wave data.
%  The analysis is carried out using the phase-shift method.
%
%% Input
%  u          u(x,t) offset-time surface wave data
%  N          Number of receivers
%  x          Location of receivers, distance from seismic source [m]
%  fs         Recording frequency [Hz]
%  cT_min     Testing Rayleigh wave phase velocity (minimum value) [m/s]
%  cT_max     Testing Rayleigh wave phase velocity (maximum value) [m/s]
%  delta_cT   Testing Rayleigh wave phase velocity increment [m/s]
% 
%% Output
%  f          Frequency [Hz]
%  c          Rayleigh wave velocity [m/s]
%  A          Summed (slant-stacked) amplitude corresponding to 
%             different combinations of omega=2*pi*f and cT
%
%% Subfunctions 
%  (None)
%
%%
function [f,c,A] = MASWaves_dispersion_imaging(u,N,x,fs,cT_min,cT_max,delta_cT)

% Converting measuring frequency from Hz to rad/sec
omega_fs = 2*pi*fs; 

% Number of samples in each trace
Lu = length(u(:,1));

% Empty matrices with Lu lines and n columns
U = zeros(Lu,N);
P = zeros(Lu,N);
Unorm = zeros(Lu,N);

% Apply discrete Fourier transform to the time axis of u
% (The MATLAB function fft( ) is used)
for j = 1:N
    U(:,j) = fft(u(:,j));
end

% Number of samples in each transformed trace
LU = length(U(:,1));

% Normalize U in offset and frequency domains
% Compute the phase spectrum of U
i = sqrt(-1);
for j = 1:N
    for k = 1:LU
        Unorm(k,j) = U(k,j)/abs(U(k,j));
    end
    P(:,j) = exp(i.*-angle(U(:,j)));
end

% Frequency range for U
% Lomega = LU
omega = (0:LU-1)*(omega_fs/LU); 

% Compute the slant-stack (summed) amplitude corresponding to each set of
% omega and cT, A(omega,cT).
cT = cT_min:delta_cT:cT_max;
LcT = length(cT);

% Empty matrices with LU lines and LcT columns
c = zeros(LU,LcT);
f = zeros(LU,LcT);
A = zeros(LU,LcT); 

for j = 1:LU % Frequency component j
    for k = 1:LcT % Testing Rayleigh wave phase velocity k
        % Frequency (in [Hz]) corresponding to angular frequency omega
        f(j,k) = omega(j)/(2*pi);
        % Testing phase velocity [m/s]
        c(j,k) = cT(k);
        % Determining the amount of phase shifts required to counterbalance
        % the time delay corresponding to specific offsets for a given set 
        % of omega and cT
        delta = omega(j)/cT(k);
        % Applying the phase shifts (delta) to distinct traces of the 
        % transformed shot gather
        % Obtaining the (normalized) slant-stack (summed) amplitude
        % corresponding to each set of omega and cT
        temp = 0;
        for l = 1:N
            temp = temp + exp(-i*delta*x(l))*P(j,l);     
        end
        % Compute absolute value and normalize with respect to number of 
        % receivers
        A(j,k) = abs(temp)/N; 
    end
end
end
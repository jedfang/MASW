%  MASWaves Dispersion
%  Version: 07.2017
%%
%  [u,T,Tmax,L,x] = MASWaves_read_data(Filename,HeaderLines,fs,N,dx,x1,Direction)
%%
%  The function MASWaves_read_data loads recorded surface wave data into
%  MATLAB and determines the length of the receiver spread, the location
%  of individual receivers and the total recording time.
%
%% Input
%  Filename     Path of file where recorded data is stored [string]
%                - Recorded data should be stored in an ASCII-delimited
%                  text file (in ASCII format).
%                - Each recorded surface wave trace should be stored in
%                  a single column.
%  HeaderLines  Number of header lines
%  fs           Measuring frequency [Hz]
%  N            Number of receivers
%  dx           Receiver spacing [m]
%  x1           Source offset [m]
%  Direction    'forward' or 'backward' [string]
%                - 'forward':  Forward measurement.
%                              Seismic source is applied next to receiver 1
%                - 'backward': Backward measurement.
%                              Seismic source is applied next to receiver N
%
%% Output
%  u           u(x,t) offset-time shot gather
%  T           Time of individual recordings [s]
%  Tmax        Total recording time [s]
%  L           Length of receiver spread [m]
%  x           Location of receivers, distance from seismic source [m]
%
%% Subfunctions
%  (None)
%
%%
function [u,T,Tmax,L,x] = MASWaves_read_data(Filename,HeaderLines,fs,N,dx,x1,Direction)

% Read recorded data into matrix u
if strcmp(Direction,'forward')
    u = dlmread(Filename,'',HeaderLines,0);
elseif strcmp(Direction,'backward')
    u_temp = dlmread(Filename,'',HeaderLines,0); 
    for j = 1:N
        u(:,j) = u_temp(:,N-j+1);
    end
else
    disp('Error')
    u = [];
end

% Total recording time [s]
Tmax = length(u(:,1))/fs - 1/fs; 

% Time of individual recordings [s]
T = linspace(0,Tmax,length(u(:,1))); 

% Length of receiver spread [m]
L = (N-1)*dx; 

% Location of receivers, distance from seismic source [m]
x = (x1):dx:(L+x1);

end

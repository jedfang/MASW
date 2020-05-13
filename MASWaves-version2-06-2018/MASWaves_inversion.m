%  MASWaves Inversion
%  Version: 07.2017
%%
%  [c_t,lambda_t,e] = MASWaves_inversion(c_test,h,alpha,beta,rho,n,...
%    up_low_boundaries,c_curve0,lambda_curve0,c_curve0_up,lambda_curve0_up,...
%    c_curve0_low,lambda_curve0_low)
%
%%
%  The function MASWaves_inversion can be used to carry out the inversion 
%  analysis through a single .m file (manual inversion). 
%  The function: 
%  (1) computes the theoretical fundamental mode dispersion  
%  curve for the layer model defined by h, alpha, beta, rho and n at the same 
%  wavelengths as are included in the experimental curve, 
%  (2) plots the theoretical and experimental curves and 
%  (3) evaluates the misfit between the theoretical and experimental curves. 
%  
%  For each iteration, the function MASWaves_inversion allows the user to 
%  choose between saving the theoretical dispersion curve obtained in the 
%  current iteration (in a text file), to stop without saving or to iterate 
%  again.
%
%% Input
%  c_test            Testing Rayleigh wave phase velocity vector [m/s]
%  h                 Layer thicknesses [m] (vector of length n)
%  alpha             Compressional wave velocity [m/s] (vector of length n+1)
%  beta              Shear wave velocity [m/s] (vector of length n+1)
%  rho               Mass density [kg/m^3] (vector of length n+1)
%  n                 Number of finite thickness layers
%
%  up_low_boundaries - 'yes'      Upper/lower boundaries for the experimental 
%                                 fundamental mode dispersion curve are wanted.
%                    - 'no'       Upper/lower boundaries for the experimental 
%                                 fundamental mode dispersion curve are not wanted.
%
%  Experimental fundamental mode dispersion curve
%  c_curve0          Phase velocity [m/s]
%  lambda_curve0     Wavelength [m]
%  c_curve0_up       Phase velocity, upper bound curve [m/s]
%                    (Can be assigned as 'nan' or [] if upper/lower boundaries
%                    are not wanted.)
%  lambda_curve0_up  Wavelength, upper bound curve [m]
%                    (Can be assigned as 'nan' or [] if upper/lower boundaries
%                    are not wanted.)
%  c_curve0_low      Phase velocity, upper bound curve [m/s]
%                    (Can be assigned as 'nan' or [] if upper/lower boundaries
%                    are not wanted.)
%  lambda_curve0_low Wavelength, upper bound curve [m]
%                    (Can be assigned as 'nan' or [] if upper/lower boundaries
%                    are not wanted.)
%
%% Output
%  c_t               Phase velocity (fundamental mode theoretical
%                    dispersion curve) [m/s]
%  lambda_t          Wavelength (theoretical fundamental mode dispersion 
%                    curve) [m]
%  e                 Estimated misfit between theoretical and experimental
%                    dispersion curves [%]
%
%% Subfunctions
%  MASWaves_theoretical_dispersion_curve
%  MASWaves_plot_theor_exp_dispersion_curves
%  MASWaves_misfit
%
%%
function [c_t,lambda_t,e] = MASWaves_inversion(c_test,h,alpha,beta,rho,n,...
    up_low_boundaries,c_curve0,lambda_curve0,c_curve0_up,lambda_curve0_up,...
    c_curve0_low,lambda_curve0_low)

% Compute theoretical dispersion curve for the current set of model
% parameters
[c_t,lambda_t] = MASWaves_theoretical_dispersion_curve...
    (c_test,lambda_curve0,h,alpha,beta,rho,n);

% Plot theoretical and experimental dispersion curves
FigWidth = 8; % cm
FigHeight = 10; % cm
FigFontSize = 8; % pt
figure
MASWaves_plot_theor_exp_dispersion_curves(c_t,lambda_t,...
    c_curve0,lambda_curve0,c_curve0_up,lambda_curve0_up,...
    c_curve0_low,lambda_curve0_low,up_low_boundaries,...
    FigWidth,FigHeight,FigFontSize)
hold on
% Compute misfit [%] between theoretical and experimental dispersion curves
e = MASWaves_misfit(c_t,c_curve0);
title(['Misfit: \epsilon = ', num2str(e), '%'])

i = input('Iterate (input 1) \nSave the theoretical dispersion curve (in a text file) and stop iterating (input 2) \nStop iterating without saving the theoretical dispersion curve (in a text file) (input 0): ');

   
if i == 1
    
    ii = input('Change layer thicknesses and shear wave velocity (input 1) \nChange shear wave velocity (input 2)\nChange layer thicknesses (input 3): ');
    
    if ii == 1
        disp('Layer thicknesses in last iteration')
        fprintf('%6.2f', h)
        fprintf('\n')
        h = input(['Input a new layer thickness vector (of length ', num2str(n),') ' ]);
        
        fprintf('\n')
        disp('Shear wave velocities in last iteration')
        fprintf('%8.2f', beta)
        fprintf('\n')
        beta = input(['Input a new shear wave velocity vector (of length ', num2str(n+1), ') ' ]);
        
    elseif ii == 2
        
        fprintf('\n')
        disp('Shear wave velocities in last iteration')
        fprintf('%8.2f', beta)
        fprintf('\n')
        beta = input(['Input a new shear wave velocity vector (of length ', num2str(n+1), ') ' ]);
        
    elseif ii == 3
        
        disp('Layer thicknesses in last iteration')
        fprintf('%6.2f', h)
        fprintf('\n')
        h = input(['Input a new layer thickness vector (of length ', num2str(n),') ' ]);
        
    end
    
    [c_t,lambda_t,e] = MASWaves_inversion(c_test,h,alpha,beta,rho,n,...
        up_low_boundaries,c_curve0,lambda_curve0,c_curve0_up,lambda_curve0_up,...
        c_curve0_low,lambda_curve0_low);

elseif i == 2
    
    A = [c_t' ; lambda_t]';
    fileID = fopen('Results.txt','w');
    fprintf(fileID,'%16s %16s\r\n','c (theor.)','lambda (theor.)');
    fprintf(fileID,'%16.4f %16.4f\r\n',A);
    fclose(fileID);
    
end
end


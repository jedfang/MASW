%  MASWaves Inversion
%  Version: 07.2017
%%
%  Ke = MASWaves_theoretical_dispersion_curve(c_test,k,h,alpha,beta,rho,n)
%%
%  The function MASWaves_theoretical_dispersion_curve computes the
%  fundamental mode theortical dispersion curve for the layer model defined 
%  by h, alpha, beta, rho and n at wavelengths lambda.
%
%% Input
%  c_test        Testing Rayleigh wave phase velocity vector [m/s]
%  lambda        Wavelength vector [m]
%  h             Layer thicknesses [m] (vector of length n)
%  alpha         Compressional wave velocity [m/s] (vector of length n+1)
%  beta          Shear wave velocity [m/s] (vector of length n+1)
%  rho           Mass density [kg/m^3] (vector of length n+1)
%  n             Number of finite thickness layers
%
%% Output
%  c_t           Rayleigh wave phase velocity vector (theoretical
%                fundamental mode dispersion curve) [m/s]
%  lambda_t      Rayleigh wave wavelength (theoretical fundamental mode 
%                dispersion curve) [m]
%
%% Subfunctions
%  MASWaves_stiffness_matrix
%
%%
function [c_t,lambda_t] = MASWaves_theoretical_dispersion_curve(c_test,lambda,h,alpha,beta,rho,n)

% Compute wave numbers that correspond to wavelengths lambda
k = (2*pi)./lambda;

D = zeros(length(c_test),length(k));
c_t = zeros(length(k),1);

% For each wave number, recompute the system stiffness matrix using
% different values of c_test until its determinant has a sign change. 
for l = 1:length(k)
    for m = 1:length(c_test)
        D(l,m) = MASWaves_stiffness_matrix(c_test(m),k(l),h,alpha,beta,rho,n);
        if m==1;
            sign_old = sign(D(l,m));
        else
            sign_old = signD;
        end
        signD = sign(D(l,m));
        if sign_old*signD == -1
            c_t(l)=c_test(m);
            lambda_t(l)=2*pi/k(l);
            break
        end
    end
end
end
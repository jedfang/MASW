%  MASWaves Inversion
%  Version: 07.2017
%%
%  e = MASWaves_plot_theor_exp_dispersion_curves(c_t,lamda_t,...
%    c_curve0,lambda_curve0)
%
%%
%  The function MASWaves_misfit is used to evaluate the misfit between the 
%  theoretical and the experimental fundamental mode dispersion curves. 
%  The theoretical and experimental curves are assumed to have been 
%  evaluated at the same wavelengths. 
%
%% Input
%
%  Theoretical fundamental mode dispersion curve
%  c_t               Phase velocity [m/s]
%
%  Experimental fundamental mode dispersion curve
%  c_curve0          Phase velocity [m/s]
%
%% Output
%  e                 Misfit [%]
%
%% Subfunctions
%  (None)
%
%%
function  e = MASWaves_misfit(c_t,c_curve0)

Q = length(c_t); % Number of data points in theoretical and experimental dispersion curves.

temp = 0;
for q = 1:Q
    temp = temp + sqrt((c_curve0(q) - c_t(q))^2)/c_curve0(q);
end
e = 1/Q * temp * 100;

disp(['Misfit: e = ', num2str(e), '%'])

end

%  MASWaves Inversion
%  Version: 07.2017
%%
%  Ke_halfspace = MASWaves_Ke_halfspace(alpha,beta,rho,c_test,k)
%%
%  The function MASWaves_Ke_halfspace computes the element stiffness matrix
%  for the half-space (layer n+1) of the stratified earth model that is
%  used in the inversion analysis.
%
%% Input
%  alpha         Half-space compressional wave velocity [m/s]
%  beta          Half-space shear wave velocity [m/s]
%  rho           Half-space mass density [kg/m^3]
%  c_test        Testing Rayleigh wave phase velocity [m/s]
%  k             Wave number
%
%% Output
%  Ke_halfspace  Half-space element stiffness matrix
%
%% Subfunctions
%  (None)
%
%%
function Ke_halfspace = MASWaves_Ke_halfspace(alpha,beta,rho,c_test,k)

r = sqrt(1-c_test^2/alpha^2);
s = sqrt(1-c_test^2/beta^2);

k_11 = k*rho*beta^2*(r*(1-s^2))/(1-r*s);
k_12 = k*rho*beta^2*(1-s^2)/(1-r*s) - 2*k*rho*beta^2;
k_21 = k_12;
k_22 = k*rho*beta^2*(s*(1-s^2))/(1-r*s);

Ke_halfspace = [k_11 k_12 ; k_21 k_22];

end
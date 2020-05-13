%  MASWaves Inversion
%  Version: 07.2017
%%
%  Ke = MASWaves_Ke_layer(h,alpha,beta,rho,c_test,k)
%%
%  The function MASWaves_Ke_layer computes the element stiffness matrix
%  of the j-th layer (j = 1,...,n) of the stratified earth
%  model that is used in the inversion analysis.
%
%% Input
%  h             Layer thickness (thickness of the j-th finite thickness
%                layer) [m]
%  alpha         Compressional wave velocity of the j-th layer [m/s]
%  beta          Shear wave velocity of the j-th layer [m/s]
%  rho           Mass density of the j-th layer [kg/m^3]
%  c_test        Testing Rayleigh wave phase velocity [m/s]
%  k             Wave number
%
%% Output
%  Ke            Element stiffness matrix of the j-th layer
%
%% Subfunctions
%  (None)
%
%%
function Ke = MASWaves_Ke_layer(h,alpha,beta,rho,c_test,k)

r = sqrt(1-c_test^2/alpha^2);
s = sqrt(1-c_test^2/beta^2);

Cr = cosh(k*r*h);
Sr = sinh(k*r*h);
Cs = cosh(k*s*h);
Ss = sinh(k*s*h);

D = 2*(1-Cr*Cs) + (1/(r*s) + r*s)*Sr*Ss;

k11_e = (k*rho*c_test^2)/D * (s^(-1)*Cr*Ss - r*Sr*Cs);
k12_e = (k*rho*c_test^2)/D * (Cr*Cs - r*s*Sr*Ss - 1) - k*rho*beta^2*(1+s^2);
k13_e = (k*rho*c_test^2)/D * (r*Sr - s^(-1)*Ss);
k14_e = (k*rho*c_test^2)/D * (-Cr + Cs);
k21_e = k12_e;
k22_e = (k*rho*c_test^2)/D * (r^(-1)*Sr*Cs - s*Cr*Ss);
k23_e = -k14_e;
k24_e = (k*rho*c_test^2)/D * (-r^(-1)*Sr + s*Ss);
k31_e = k13_e;
k32_e = k23_e;
k33_e = k11_e;
k34_e = -k12_e;
k41_e = k14_e;
k42_e = k24_e;
k43_e = -k21_e;
k44_e = k22_e;

Ke = [k11_e k12_e k13_e k14_e;
    k21_e k22_e k23_e k24_e;
    k31_e k32_e k33_e k34_e;
    k41_e k42_e k43_e k44_e];

end
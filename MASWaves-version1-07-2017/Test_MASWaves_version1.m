%  MASWaves - Quick Start Guide
%  Version: 07.2017
%%
% Note: 
% This script is written in cell mode, i.e. the code is divided into
% several code sections. Each code section begins with two comment characters (%%).
% Below are instructions on evaluating code sections.
%
% - Run the code in the current section
% (1) Place the cursor in the code section. 
% (2) On the Editor tab, click Run Section OR hit Crtl and Return on the
% keyboard
% 
% - Run the code in the current section, and then move to the next section.
% (1) Place the cursor in the code section. 
% (2) On the Editor tab, click Run and Advance
%
%%
clear all
close all
clc
%%
% ------------- Dispersion analysis ---------------
%%
Filename = 'SampleData.dat';
HeaderLines = 7;
fs = 1000; % Hz
N = 24;
x1 = 10; % m
dx = 1; % m
Direction = 'forward';

[u,T,Tmax,L,x] = MASWaves_read_data(Filename,HeaderLines,fs,N,dx,x1,Direction);

%%
du = 1/75;
FigWidth = 6; % cm
FigHeight = 8; % cm
FigFontSize = 8; % pt

figure
MASWaves_plot_data(u,N,dx,x1,L,T,Tmax,du,FigWidth,FigHeight,FigFontSize)

%%
cT_min = 50; % m/s
cT_max = 400; % m/s
delta_cT = 1; % m/s

[f,c,A] = MASWaves_dispersion_imaging(u,N,x,fs,cT_min,cT_max,delta_cT);

%%
resolution = 100;
fmin = 0; % Hz
fmax = 50; % Hz
FigWidth = 7; % cm
FigHeight = 7; % cm
FigFontSize = 8; % pt
figure
[fplot,cplot,Aplot] = MASWaves_plot_dispersion_image_2D(f,c,A,fmin,fmax,...
    resolution,FigWidth,FigHeight,FigFontSize);

%%
fmin = 1; % Hz
fmax = 50; % Hz
FigWidth = 10; % cm
FigHeight = 10; % cm
FigFontSize = 8; % pt
figure
[fplot,cplot,Aplot] = MASWaves_plot_dispersion_image_3D(f,c,A,fmin,fmax,...
    FigWidth,FigHeight,FigFontSize);

%%
f_receivers = 4.5; % Hz
select = 'numbers';
up_low_boundary = 'yes'; 
p = 95; % Percentage
[f_curve0,c_curve0,lambda_curve0,...
    f_curve0_up,c_curve0_up,lambda_curve0_up,...
    f_curve0_low,c_curve0_low,lambda_curve0_low] = ...
    MASWaves_extract_dispersion_curve(f,c,A,fmin,fmax,f_receivers,...
    select,up_low_boundary,p);

%%
FigWidth = 9; % cm
FigHeight = 6; % cm
FigFontSize = 8; % pt
type = 'f_c';
up_low_boundary = 'yes';
figure
MASWaves_plot_dispersion_curve(f_curve0,c_curve0,lambda_curve0,...
     f_curve0_up,c_curve0_up,lambda_curve0_up,f_curve0_low,c_curve0_low,...
     lambda_curve0_low,type,up_low_boundary,FigWidth,FigHeight,FigFontSize)

FigWidth = 7; % cm
FigHeight = 9; % cm
FigFontSize = 8; % pt
type = 'c_lambda';
up_low_boundary = 'yes';
figure
MASWaves_plot_dispersion_curve(f_curve0,c_curve0,lambda_curve0,...
     f_curve0_up,c_curve0_up,lambda_curve0_up,f_curve0_low,c_curve0_low,...
     lambda_curve0_low,type,up_low_boundary,FigWidth,FigHeight,FigFontSize)
 
%%
% ------------- Inversion analysis ---------------
%%
% Repeated use of MASWaves_theoretical_dispersion_curve.m, MASWaves_misfit.m 
% and MASWaves_plot_theor_exp_dispersion_curves.m
% (For iteration, the layer parameters should be updated and this code section run again). 

c_test_min = 0; % m/s
c_test_max = 500; % m/+s
delta_c_test = 0.5; % m/s
c_test = c_test_min:delta_c_test:c_test_max; % m/s

% Layer parameters
n = 6;
alpha = [1440 1440 1440 1440 1440 1440 1440]; % m/s
h = [1 1 2 2 4 5 Inf]; % m
beta = [75 90 150 180 240 290 290]; % m/s
rho = [1850 1850 1850 1850 1850 1850 1850]; % kg/m^3

up_low_boundary = 'yes';
[c_t,lambda_t] = MASWaves_theoretical_dispersion_curve...
    (c_test,lambda_curve0,h,alpha,beta,rho,n);

up_low_boundary = 'yes';
FigWidth = 8; % cm
FigHeight = 10; % cm
FigFontSize = 8; % pt
figure
MASWaves_plot_theor_exp_dispersion_curves(c_t,lambda_t,...
    c_curve0,lambda_curve0,c_curve0_up,lambda_curve0_up,...
    c_curve0_low,lambda_curve0_low,up_low_boundary,...
    FigWidth,FigHeight,FigFontSize)

e = MASWaves_misfit(c_t,c_curve0);

%%
% Use of MASWaves_inversion

c_test_min = 0; % m/s
c_test_max = 500; % m/+s
delta_c_test = 0.5; % m/s
c_test = c_test_min:delta_c_test:c_test_max; % m/s

% Layer parameters
n = 6;
alpha = [1440 1440 1440 1440 1440 1440 1440]; % m/s
h = [1 1 2 2 4 5 Inf]; % m
beta = [75 90 150 180 240 290 290]; % m/s
rho = [1850 1850 1850 1850 1850 1850 1850]; % kg/m^3

up_low_boundary = 'yes';
[c_t,lambda_t,e] = MASWaves_inversion(c_test,h,alpha,beta,rho,n,...
    up_low_boundary,c_curve0,lambda_curve0,c_curve0_up,lambda_curve0_up,...
    c_curve0_low,lambda_curve0_low);

% View the results
up_low_boundary = 'yes';
FigWidth = 8; % cm
FigHeight = 10; % cm
FigFontSize = 8; % pt
figure
MASWaves_plot_theor_exp_dispersion_curves(c_t,lambda_t,...
    c_curve0,lambda_curve0,c_curve0_up,lambda_curve0_up,...
    c_curve0_low,lambda_curve0_low,up_low_boundary,...
    FigWidth,FigHeight,FigFontSize)

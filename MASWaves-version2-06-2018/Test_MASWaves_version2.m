%  MASWaves (version 2) - Quick Start Guide
%  Version: 06.2018
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
% ------------- Dispersion analysis / MASWaves Dispersion ---------------
%%
% For demonstration purposes, the dispersion analysis is only carried out
% for two separate multichannel surface wave analysis records in this
% example, resulting in two elementary dispersion curves
%
% For later steps of the analysis (the dispersion curve combination and the
% inversion, a new elementary dispersion curve dataset is loaded

No_records = 2;
c_all = cell(No_records,1);
lambda_all = cell(No_records,2);

Filename = cell(No_records,1);
Filename{1} = 'SampleData1_dx-1m_x1-3m.dat';
Filename{2} = 'SampleData2_dx-1m_x1-5m.dat';
HeaderLines = 7;
fs = 1000; % Hz
N = 24;
dx = [1 1]; % m
x1 = [3 5]; % m
Direction = 'forward';

%%
for i = 1:No_records
    
    [u,T,Tmax,L,x] = MASWaves_read_data(Filename{i},HeaderLines,fs,N,dx(i),x1(i),Direction);
    
    du = 1/75;
    FigWidth = 6; % cm
    FigHeight = 8; % cm
    FigFontSize = 8; % pt
    
    figure
    MASWaves_plot_data(u,N,dx(i),x1(i),L,T,Tmax,du,FigWidth,FigHeight,FigFontSize)
    
    cT_min = 50; % m/s
    cT_max = 400; % m/s
    delta_cT = 1; % m/s
    
    [f,c,A] = MASWaves_dispersion_imaging(u,N,x,fs,cT_min,cT_max,delta_cT);
    
    resolution = 100;
    fmin = 0; % Hz
    fmax = 50; % Hz
    FigWidth = 7; % cm
    FigHeight = 7; % cm
    FigFontSize = 8; % pt
    figure
    [fplot,cplot,Aplot] = MASWaves_plot_dispersion_image_2D(f,c,A,fmin,fmax,...
        resolution,FigWidth,FigHeight,FigFontSize);
    
    fmin = 1; % Hz
    fmax = 50; % Hz
    FigWidth = 10; % cm
    FigHeight = 10; % cm
    FigFontSize = 8; % pt
    figure
    [fplot,cplot,Aplot] = MASWaves_plot_dispersion_image_3D(f,c,A,fmin,fmax,...
        FigWidth,FigHeight,FigFontSize);
    
    f_receivers = 4.5; % Hz
    select = 'numbers';
    up_low_boundary = 'no';
    p = 95; % Percentage
    [f_curve0,c_curve0,lambda_curve0] = MASWaves_extract_dispersion_curve(f,c,A,fmin,fmax,f_receivers,...
        select,up_low_boundary,p);
   
    FigWidth = 9; % cm
    FigHeight = 6; % cm
    FigFontSize = 8; % pt
    type = 'f_c';
    up_low_boundary = 'no';
    figure
    MASWaves_plot_dispersion_curve(f_curve0,c_curve0,lambda_curve0,...
        [],[],[],[],[],[],type,up_low_boundary,FigWidth,FigHeight,FigFontSize)
  
    FigWidth = 7; % cm
    FigHeight = 9; % cm
    FigFontSize = 8; % pt
    type = 'c_lambda';
    up_low_boundary = 'no';
    figure
    MASWaves_plot_dispersion_curve(f_curve0,c_curve0,lambda_curve0,...
        [],[],[],[],[],[],type,up_low_boundary,FigWidth,FigHeight,FigFontSize)
    
    c_all{i} = c_curve0;
    lambda_all{i} = lambda_curve0;
end

%%
% ------------- Evaluation of a composite dispersion curve ---------------
clear lambda_all
clear c_all

% Import sample data
load('lambda_all.mat')
load('c_all.mat')
no_measurements = length(c_all);
%%
% Evaluate the composite dispersion curve
a = 3; % Combination parameter
[c_mean,c_min_std,c_plus_std,lambda_mean,nPoints_fin] = MASWaves_combined_dispersion_curve(no_measurements,c_all,lambda_all,a);

%%
% Display the composite dispersion curve
PlotAll = 1; % False (1 if true)
FigWidth = 7; % 
FigHeight = 9; % cm
FigFontSize = 8; % pt

MASWaves_plot_combined_dispersion_curve(no_measurements,c_all,lambda_all,...
    c_mean,c_plus_std,c_min_std,lambda_mean,PlotAll,FigWidth,FigHeight,FigFontSize);

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
h = [1 1 2 4 6 6 Inf]; % m
beta = [52 84 115 175 250 330 330]; % m/s
rho = [1850 1850 1850 1850 1850 1850 1850]; % kg/m^3

[c_t,lambda_t] = MASWaves_theoretical_dispersion_curve...
    (c_test,lambda_mean,h,alpha,beta,rho,n);
 
up_low_boundary = 'yes';
FigWidth = 8; % cm
FigHeight = 10; % cm
FigFontSize = 8; % pt
figure
MASWaves_plot_theor_exp_dispersion_curves(c_t,lambda_t,...
    c_mean,lambda_mean,c_min_std,lambda_mean,...
    c_plus_std,lambda_mean,up_low_boundary,...
    FigWidth,FigHeight,FigFontSize)

e = MASWaves_misfit(c_t,c_mean);

up_low_boundary = 'yes';
FigWidth = 16; % cm
FigHeight = 10; % cm
FigFontSize = 8; % pt
figure
MASWaves_plot_inversion_results_one_iteation...
    (c_t,lambda_t,c_mean,lambda_mean,c_min_std,lambda_mean,...
    c_plus_std,lambda_mean,n,beta,h,e,up_low_boundary,...
    FigWidth,FigHeight,FigFontSize)
%%
% Use of MASWaves_inversion
c_test_min = 0; % m/s
c_test_max = 500; % m/+s
delta_c_test = 0.5; % m/s
c_test = c_test_min:delta_c_test:c_test_max; % m/s

% Layer parameters
n = 6;
alpha = [1440 1440 1440 1440 1440 1440 1440]; % m/s
h = [1 1 2 4 6 6 Inf]; % m
beta = [52 84 115 175 250 330 330]; % m/s
rho = [1850 1850 1850 1850 1850 1850 1850]; % kg/m^3

up_low_boundary = 'yes';
[c_t,lambda_t,e] = MASWaves_inversion(c_test,h,alpha,beta,rho,n,...
    up_low_boundary,c_mean,lambda_mean,c_min_std,lambda_mean,...
    c_plus_std,lambda_mean);

% View the results
up_low_boundary = 'yes';
FigWidth = 16; % cm
FigHeight = 10; % cm
FigFontSize = 8; % pt
figure
MASWaves_plot_inversion_results_one_iteation...
    (c_t,lambda_t,c_mean,lambda_mean,c_min_std,lambda_mean,...
    c_plus_std,lambda_mean,n,beta,h,e,up_low_boundary,...
    FigWidth,FigHeight,FigFontSize)

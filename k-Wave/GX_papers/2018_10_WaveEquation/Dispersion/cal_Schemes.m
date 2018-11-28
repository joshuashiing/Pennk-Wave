% Define different cases (media parameters and frequency range) and 
% calculate the wavenumber dispersion by calling cal_k_disp.m
% Author: Guangchi Xing
% Date: 03/03/2018

clear;
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
casename = 'case_schemes.mat';
fprintf([casename, '\n']);
% Set up the parameters for Pierre shale (Wuenschel, 1965)
f0_m = 1500;     % Reference frequency of the medium (Hz)
c0_m = 2164;    % Reference phase velocity (m/s)
rho  = 2200;    % Density (kg/m^3)
Q  = 10;        % Quality factor
fl = 10;        % Lower bound of frequency range (Hz)
fu = 60;        % Upper bound of frequency range (Hz)
nf = 100;       % Number of frequency sample points
cal_k_schemes(f0_m, c0_m, rho, Q, fl, fu, nf, casename);

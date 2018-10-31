clear;
clc;

load('conversion_bs03j_V1.mat');

% =========================================================================
% Grid Parameters
% =========================================================================

Nx = 434;           % number of grid points in the x direction
Ny = 467;           % number of grid points in the y direction
dx = 0.15;          % grid point spacing in the x direction [m]
dy = 0.15;          % grid point spacing in the y direction [m]

x0 = 1621.45;       % Top-left corner x coordinate [m]
y0 = -19.95;        % Top-left corner y coordinate [m]

% =========================================================================
% Correction from the original coordinate to the k-Wave coordinate
% =========================================================================
kgrid_tmp = kWaveGrid(Nx, dx, Ny, dy);
xc = kgrid_tmp.x_vec(1) - x0;   % Correction for x direction
yc = kgrid_tmp.y_vec(1) - y0;   % Correction for y direction

% =========================================================================
% Source & Receivers
% =========================================================================

x_src = 1658.05;                            % Source location x in original coordinate [m]
y_src = 0;                                  % Source location y in original coordinate [m]
x_rec = [1635.1; 1650.1; 1658.05; 1680.1];	% Receiver location x in original coordinate [m]
y_rec = [30; 30; 30; 30];                   % Receiver location y in original coordinate [m]



f_rw_c = 100;       % Center frequency of ricker wavelet [Hz]



x_src = -975;          % Source location in the x direction [m] 25m
y_src = 12.5;       % Source location in the x direction [m] 2500 m
x_rec = -975 * ones(398, 1);          % Receiver location in the x direction [m]
y_rec = (-2487.5 : 12.5 : 2475)';        % Receiver location in the y direction [m]
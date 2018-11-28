%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate synthetic seismograms for different Q
% Author: Guangchi Xing
% Date: 11/27/2018

clear;
clc

% =========================================================================
% Grid Parameters
% =========================================================================

Nx = 150;           % number of grid points in the x direction
Ny = 600;           % number of grid points in the y direction
dx = 1.0;           % grid point spacing in the x direction [m]
dy = 1.0;           % grid point spacing in the y direction [m]

x0 = 0;             % Top-left corner x coordinate [m]
y0 = 0;             % Top-left corner y coordinate [m]

% =========================================================================
% Correction from the original coordinate to the k-Wave coordinate
% =========================================================================

kgrid_tmp = kWaveGrid(Nx, dx, Ny, dy);
xc = kgrid_tmp.x_vec(1) - x0;   % Correction for x direction
yc = kgrid_tmp.y_vec(1) - y0;   % Correction for y direction

% =========================================================================
% Medium Parameters
% =========================================================================

f0_m = 1500;           % Reference frequency [Hz]
c0 = 2164;            % Phase velocity at reference frequency [m/s]
% Q = 10;             % Quality factor
rho = 2200;           % Density [kg/m^3]

% =========================================================================
% Simulation Parameters
% =========================================================================

dt = 0.1e-3;        % Time interval [s]
t_max = 0.31;        % Simulation end time [s]
f0 = 75;           % Reference frequency for simulation
mod_mech = 'TF111110';

% =========================================================================
% Source & Receivers
% =========================================================================

x_src = Nx * dx / 2;          % Source location in the x direction [m]
y_src = Ny * dy * 1/6;       % Source location in the y direction [m]
x_rec = x_src;          % Receiver location in the x direction [m]
y_rec = Ny * dy * 5/6;        % Receiver location in the y direction [m]

% x_src = Nx * dx / 2;          % Source location in the x direction [m]
% y_src = Ny * dy / 2;       % Source location in the y direction [m]
% x_rec = x_src;          % Receiver location in the x direction [m]
% y_rec = y_src;        % Receiver location in the y direction [m]

f_rw_c = 75;        % Center frequency of ricker wavelet [Hz]
A_rw = 1e3;         % Amplitude of the ricker wavelet
t_rw_c = 1 / f_rw_c * 2;
nt_rw = 1 / f_rw_c * 4 / dt;
stf = rickerwavelet(f_rw_c, dt, nt_rw, t_rw_c, 1) * A_rw;


% =========================================================================
% MISC
% =========================================================================

x_src = x_src + xc;
y_src = y_src + yc;
x_rec = x_rec + xc;
y_rec = y_rec + yc;

f0_m = f0_m * ones(Nx, Ny);
f0 = f0 * ones(Nx, Ny);
c0 = c0 * ones(Nx, Ny);
rho = rho * ones(Nx, Ny);

% =========================================================================
% Forward Modeling for different cases
% (Test for different quality factor)
% =========================================================================

mat_name = 'case01.mat';
Q = ones(Nx, Ny) * 999999;
[d, t_axis, d2] = homog_simu(Nx, Ny, dx, dy, f0_m, f0, c0, Q, rho, ...
                    stf, x_src, y_src, x_rec, y_rec, ...
                    dt, t_max, mod_mech, -1, '');
save(mat_name, 't_axis', 'd', 'd2');

mat_name = 'case02.mat';
Q = ones(Nx, Ny) * 100;
[d, t_axis, d2] = homog_simu(Nx, Ny, dx, dy, f0_m, f0, c0, Q, rho, ...
                    stf, x_src, y_src, x_rec, y_rec, ...
                    dt, t_max, mod_mech, -1, '');
save(mat_name, 't_axis', 'd', 'd2');

mat_name = 'case03.mat';
Q = ones(Nx, Ny) * 32;
[d, t_axis, d2] = homog_simu(Nx, Ny, dx, dy, f0_m, f0, c0, Q, rho, ...
                    stf, x_src, y_src, x_rec, y_rec, ...
                    dt, t_max, mod_mech, -1, '');
save(mat_name, 't_axis', 'd', 'd2');

mat_name = 'case04.mat';
Q = ones(Nx, Ny) * 10;
[d, t_axis, d2] = homog_simu(Nx, Ny, dx, dy, f0_m, f0, c0, Q, rho, ...
                    stf, x_src, y_src, x_rec, y_rec, ...
                    dt, t_max, mod_mech, -1, '');
save(mat_name, 't_axis', 'd', 'd2');

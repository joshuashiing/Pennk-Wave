clear;
clc

% =========================================================================
% Grid Parameters
% =========================================================================

Nx = 256;           % number of grid points in the x direction
Ny = 1024;           % number of grid points in the y direction
dx = 0.5;           % grid point spacing in the x direction [m]
dy = 0.5;           % grid point spacing in the y direction [m]

% =========================================================================
% Medium Parameters
% =========================================================================

f0_m = 100;           % Reference frequency [Hz]
c0_m = 2164;          % Phase velocity at reference frequency [m/s]
% Q = 32;             % Quality factor
density = 2200;     % Density [kg/m^3]

% =========================================================================
% Source & Receivers
% =========================================================================

f_rw_c = 100;       % Center frequency of ricker wavelet [Hz]
x_src = 0.0;          % Source location in the x direction [m]
y_src = -200.0;       % Source location in the y direction [m]
x_rec = 0.0;          % Receiver location in the x direction [m]
y_rec = 200.0;        % Receiver location in the y direction [m]

% =========================================================================
% Simulation Parameters
% =========================================================================

dt = 5e-5;        % Time interval [s]
t_max = 0.40;        % Simulation end time [s]
f0 = 100;           % Reference frequency for simulation
mod_mech = 'TF111110';  % Numerical modeling mechanism

% =========================================================================
% Forward Modeling for different cases
% (Test for different quality factor)
% =========================================================================
casename = 'case01';
Q = 1e6;
[t_axis, d1, d2, d3] = homog_syn_mod(Nx, Ny, dx, dy,...
    f0_m, c0_m, Q, density, f_rw_c, x_src, y_src, x_rec, y_rec, ...
    dt, t_max, f0, mod_mech, 1, 1, casename);

casename = 'case02';
Q = 100;
[t_axis, d1, d2, d3] = homog_syn_mod(Nx, Ny, dx, dy,...
    f0_m, c0_m, Q, density, f_rw_c, x_src, y_src, x_rec, y_rec, ...
    dt, t_max, f0, mod_mech, 1, 1, casename);

casename = 'case03';
Q = 32;
[t_axis, d1, d2, d3] = homog_syn_mod(Nx, Ny, dx, dy,...
    f0_m, c0_m, Q, density, f_rw_c, x_src, y_src, x_rec, y_rec, ...
    dt, t_max, f0, mod_mech, 1, 1, casename);

casename = 'case04';
Q = 10;
[t_axis, d1, d2, d3] = homog_syn_mod(Nx, Ny, dx, dy,...
    f0_m, c0_m, Q, density, f_rw_c, x_src, y_src, x_rec, y_rec, ...
    dt, t_max, f0, mod_mech, 1, 1, casename);

plot(t_axis, d1, 'k', 'linewidth', 2); hold on;
plot(t_axis, d2, 'b--', 'linewidth', 2); hold on;
plot(t_axis, d3, 'r--', 'linewidth', 2);
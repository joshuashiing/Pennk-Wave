%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate snapshot for separate effects (amplitude decay & velocity
% dispersion)

clear;
clc

% =========================================================================
% Grid Parameters
% =========================================================================

Nx = 512;           % number of grid points in the x direction
Ny = 512;           % number of grid points in the y direction
dx = 1.0;           % grid point spacing in the x direction [m]
dy = 1.0;           % grid point spacing in the y direction [m]

% =========================================================================
% Medium Parameters
% =========================================================================

f0_m = 1500;           % Reference frequency [Hz]
c0_m = 2164;          % Phase velocity at reference frequency [m/s]
% Q = 32;             % Quality factor
density = 2200;     % Density [kg/m^3]

% =========================================================================
% Source & Receivers
% =========================================================================

f_rw_c = 100;       % Center frequency of ricker wavelet [Hz]
x_src = 0.0;          % Source location in the x direction [m]
y_src = 0.0;       % Source location in the y direction [m]
x_rec = 0.0;          % Receiver location in the x direction [m]
y_rec = 200.0;        % Receiver location in the y direction [m]

% =========================================================================
% Simulation Parameters
% =========================================================================

dt = 2e-4;        % Time interval [s]
t_max = 0.15;        % Simulation end time [s]
f0 = 100;           % Reference frequency for simulation
% mod_mech = 'TF111110';  % Numerical modeling mechanism

% =========================================================================
% Forward Modeling for different cases
% (Test for different quality factor)
% =========================================================================
casename = 'testS_case01';
mod_mech = 'lossless';
Q = 0.1;  
t_snap = 0.12;
homog_snap(Nx, Ny, dx, dy, f0, c0_m, Q, density, f_rw_c, x_src, ...
    y_src, x_rec, y_rec, dt, t_max, f0, mod_mech, casename, t_snap);

casename = 'testS_case02';
mod_mech = 'TF111110_ld';
Q = 32;
t_snap = 0.12;
homog_snap(Nx, Ny, dx, dy, f0, c0_m, Q, density, f_rw_c, x_src, ...
    y_src, x_rec, y_rec, dt, t_max, f0, mod_mech, casename, t_snap);

casename = 'testS_case03';
mod_mech = 'TF111110_dd';
Q = 32;
t_snap = 0.12;
homog_snap(Nx, Ny, dx, dy, f0_m, c0_m, Q, density, f_rw_c, x_src, ...
    y_src, x_rec, y_rec, dt, t_max, f0, mod_mech, casename, t_snap);

casename = 'testS_case04';
mod_mech = 'TF111110';
Q = 32;
t_snap = 0.12;
homog_snap(Nx, Ny, dx, dy, f0_m, c0_m, Q, density, f_rw_c, x_src, ...
    y_src, x_rec, y_rec, dt, t_max, f0, mod_mech, casename, t_snap);

clear;
clc;

% =========================================================================
% Grid Parameters
% =========================================================================

Nx = 200;           % number of grid points in the x direction
Ny = 200;           % number of grid points in the y direction
dx = 8;           % grid point spacing in the x direction [m]
dy = 8;           % grid point spacing in the y direction [m]

% =========================================================================
% Source & Receivers
% =========================================================================

f_rw_c = 25;       % Center frequency of ricker wavelet [Hz]
x_src = 0;          % Source location in the x direction [m]
y_src = 0;       % Source location in the y direction [m]
x_rec = 400.;          % Receiver location in the x direction [m]
y_rec = 0;        % Receiver location in the y direction [m]

% =========================================================================
% Simulation Parameters
% =========================================================================

dt = 1e-3;        % Time interval [s]
t_max = 0.6;        % Simulation end time [s]
t_snap = 0.30;

% =========================================================================
% Different Cases
% =========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model 1

model_f0_m = ones(Nx, Ny) * 25;
model_f0 = ones(Nx, Ny) * 25; 
model_c0_m = ones(Nx, Ny) * 3600;   model_c0_m(1 : 130, :) = 1800; 
model_Q = ones(Nx, Ny) * 30.;       model_Q(1 : 130, :) = 30;
model_density = ones(Nx, Ny) * 2200;

mat_name = 'case01_a';
mod_mech = 'TF111110';  % Numerical modeling mechanism

heter_snap(Nx, Ny, dx, dy, model_f0_m, model_f0, model_c0_m,...
    model_Q, model_density, f_rw_c, x_src, y_src, x_rec, y_rec, ...
    dt, t_max, mod_mech, t_snap, mat_name)

mat_name = 'case01_b';
mod_mech = 'TZ14';  % Numerical modeling mechanism

heter_snap(Nx, Ny, dx, dy, model_f0_m, model_f0, model_c0_m,...
    model_Q, model_density, f_rw_c, x_src, y_src, x_rec, y_rec, ...
    dt, t_max, mod_mech, t_snap, mat_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model 2

model_f0_m = ones(Nx, Ny) * 25;
model_f0 = ones(Nx, Ny) * 25; 
model_c0_m = ones(Nx, Ny) * 1800;   model_c0_m(1 : 130, :) = 1800; 
model_Q = ones(Nx, Ny) * 100;       model_Q(1 : 130, :) = 30;
model_density = ones(Nx, Ny) * 2200;

mat_name = 'case02_a';
mod_mech = 'TF111110';  % Numerical modeling mechanism

heter_snap(Nx, Ny, dx, dy, model_f0_m, model_f0, model_c0_m,...
    model_Q, model_density, f_rw_c, x_src, y_src, x_rec, y_rec, ...
    dt, t_max, mod_mech, t_snap, mat_name)

mat_name = 'case02_b';
mod_mech = 'TZ14';  % Numerical modeling mechanism

heter_snap(Nx, Ny, dx, dy, model_f0_m, model_f0, model_c0_m,...
    model_Q, model_density, f_rw_c, x_src, y_src, x_rec, y_rec, ...
    dt, t_max, mod_mech, t_snap, mat_name)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model 3

model_f0_m = ones(Nx, Ny) * 25;
model_f0 = ones(Nx, Ny) * 25; 
model_c0_m = ones(Nx, Ny) * 3600;   model_c0_m(1 : 130, :) = 1800; 
model_Q = ones(Nx, Ny) * 100;       model_Q(1 : 130, :) = 30;
model_density = ones(Nx, Ny) * 2200;

mat_name = 'case03_a';
mod_mech = 'TF111110';  % Numerical modeling mechanism

heter_snap(Nx, Ny, dx, dy, model_f0_m, model_f0, model_c0_m,...
    model_Q, model_density, f_rw_c, x_src, y_src, x_rec, y_rec, ...
    dt, t_max, mod_mech, t_snap, mat_name)

mat_name = 'case03_b';
mod_mech = 'TZ14';  % Numerical modeling mechanism

heter_snap(Nx, Ny, dx, dy, model_f0_m, model_f0, model_c0_m,...
    model_Q, model_density, f_rw_c, x_src, y_src, x_rec, y_rec, ...
    dt, t_max, mod_mech, t_snap, mat_name)



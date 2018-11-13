clear;
clc;

load('BP_cut_model_gas_chimney');

% =========================================================================
% Grid Parameters
% =========================================================================

Nx = nz;           % number of grid points in the x direction
Ny = nx;           % number of grid points in the y direction
dx = 12.5;           % grid point spacing in the x direction [m]
dy = 12.5;           % grid point spacing in the y direction [m]
%%% y_vec -2487.5 : 12.5 : 2475
%%% x_vec -1000 : 12.5 : 1000

% =========================================================================
% Source & Receivers
% =========================================================================

f_rw_c = 15;       % Center frequency of ricker wavelet [Hz]
x_src = -975;          % Source location in the x direction [m] 25m
y_src = 12.5;       % Source location in the x direction [m] 2500 m
x_rec = -975 * ones(398, 1);          % Receiver location in the x direction [m]
y_rec = (-2487.5 : 12.5 : 2475)';        % Receiver location in the y direction [m]

% =========================================================================
% Simulation Parameters
% =========================================================================

dt = 1e-3;        % Time interval [s]
t_max = 2;        % Simulation end time [s]
t_snap = 0.8;

% =========================================================================
% Media Parameters
% =========================================================================

model_f0_m = ones(Nx, Ny) * 100;
model_f0 = ones(Nx, Ny) * 15; 
model_c0_m = velp;
model_Q = Qp;
model_density = ones(Nx, Ny) * 2200;

% =========================================================================
% Different Cases
% =========================================================================

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% case 1
mat_name = 'case01';
mod_mech = 'lossless';
[t_axis, d] = heter_simu(Nx, Ny, dx, dy, model_f0, model_f0, model_c0_m,...
    model_Q, model_density, f_rw_c, x_src, y_src, x_rec, y_rec, ...
    dt, t_max, mod_mech, t_snap, mat_name);
save(mat_name, 't_axis', 'd', '-append');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% case 2
mat_name = 'case02';
mod_mech = 'TF111110';
[t_axis, d] = heter_simu(Nx, Ny, dx, dy, model_f0_m, model_f0, model_c0_m,...
    model_Q, model_density, f_rw_c, x_src, y_src, x_rec, y_rec, ...
    dt, t_max, mod_mech, t_snap, mat_name);
save(mat_name, 't_axis', 'd', '-append');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% case 3
mat_name = 'case03';
mod_mech = 'TZ14';
[t_axis, d] = heter_simu(Nx, Ny, dx, dy, model_f0_m, model_f0, model_c0_m,...
    model_Q, model_density, f_rw_c, x_src, y_src, x_rec, y_rec, ...
    dt, t_max, mod_mech, t_snap, mat_name);
save(mat_name, 't_axis', 'd', '-append');


% pmax = max(abs(d(:)));
% plim = pmax / 50;
% imagesc(d', [-plim, plim]);
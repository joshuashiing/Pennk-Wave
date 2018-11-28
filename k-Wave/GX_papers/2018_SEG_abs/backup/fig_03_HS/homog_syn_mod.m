function [t_axis, d1, d2, d3] = homog_syn_mod(Nx, Ny, dx, dy,...
    f0_m, c0_m, Q, density, f_rw_c, x_src, y_src, x_rec, y_rec, ...
    dt, t_max, f0, mod_mech, a_green, v_green, mat_name)

% clear;
% clc
% 
% % =========================================================================
% % Grid Parameters
% % =========================================================================
% 
% Nx = 128;           % number of grid points in the x direction
% Ny = 256;           % number of grid points in the y direction
% dx = 0.5;           % grid point spacing in the x direction [m]
% dy = 0.5;           % grid point spacing in the y direction [m]
% 
% % =========================================================================
% % Medium Parameters
% % =========================================================================
% 
% f0_m = 200;           % Reference frequency [Hz]
% c0_m = 2089;          % Phase velocity at reference frequency [m/s]
% Q = 10;             % Quality factor
% density = 2200;     % Density [kg/m^3]
% 
% % =========================================================================
% % Source & Receivers
% % =========================================================================
% 
% f_rw_c = 200;       % Center frequency of ricker wavelet [Hz]
% x_src = 0.;          % Source location in the x direction [m]
% y_src = -34.5;       % Source location in the y direction [m]
% x_rec = 0.;          % Receiver location in the x direction [m]
% y_rec = 34.5;        % Receiver location in the y direction [m]
% 
% % =========================================================================
% % Simulation Parameters
% % =========================================================================
% 
% dt = 1e-5;        % Time interval [s]
% t_max = 0.06;        % Simulation end time [s]
% f0 = 200;           % Reference frequency for simulation
% mod_mech = 'TF111110';  % Numerical modeling mechanism

% =========================================================================
% Simulation Set-up
% =========================================================================

% Grid
kgrid = kWaveGrid(Nx, dx, Ny, dy);
kgrid.t_array = 0 : dt: t_max;

% Medium
model_f0_m = f0_m * ones(Nx, Ny);
model_c0_m = c0_m * ones(Nx, Ny);
model_f0   = f0   * ones(Nx, Ny);
model_Q    = Q    * ones(Nx, Ny);
model_gamma= atan(1 ./ model_Q) / pi;
model_c0   = model_c0_m .* (model_f0 ./ model_f0_m) .^ model_gamma;

medium.f0 = model_f0;
medium.c0 = model_c0;
medium.Q  = model_Q;
medium.density = density * ones(Nx, Ny);
medium.sound_speed = medium.c0;
medium.mod_mech = mod_mech;

% Source
[nx_src, ny_src] = close_grid_2d(kgrid, x_src, y_src);
source.p_mask = zeros(Nx, Ny);
source.p_mask(nx_src, ny_src) = 1;
t_rw_c = 1 / f_rw_c * 2;
nt_rw = 1 / f_rw_c * 4 / dt;
rw = rickerwavelet(f_rw_c, dt, nt_rw, t_rw_c, 1);
source.p = rw*100;

% Sensor
[nx_rec, ny_rec] = close_grid_2d(kgrid, x_rec, y_rec);
sensor.mask = zeros(Nx, Ny);
sensor.mask(nx_rec, ny_rec) = 1;

% =========================================================================
% Compute the Analytical Solution
% =========================================================================

if a_green
    d1 = lole_analytical_2d(kgrid, medium, source, sensor);
else
    d1 = zeros(size(kgrid.t_array));
end

if v_green
    d2 = kjar_analytical_2d(kgrid, medium, source, sensor);
else
    d2 = zeros(size(kgrid.t_array));
end

if ~strcmp(mod_mech, '')
    d3 = kspaceFirstOrder2D(kgrid, medium, source, sensor);
else
    d3 = zeros(size(kgrid.t_array));
end

t_axis = kgrid.t_array;

if ~strcmp(mat_name, '')
    save(mat_name, 't_axis', 'd1', 'd2', 'd3');
end
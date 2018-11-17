function [d, p_save, t_axis] = rtm_fm_simu_dev(Nx, Ny, dx, dy, model_f0_m, model_f0, ...
    model_c0_m, model_Q, model_density, stf, x_src, y_src, x_rec, y_rec, ...
    dt, t_max, mod_mech, f_cutoff, taper_ratio, t_snap, mat_name, args)

% clear;
% clc
% 
% % =========================================================================
% % Grid Parameters
% % =========================================================================
% 
% Nx = 200;           % number of grid points in the x direction
% Ny = 200;           % number of grid points in the y direction
% dx = 8;           % grid point spacing in the x direction [m]
% dy = 8;           % grid point spacing in the y direction [m]
% 
% % =========================================================================
% % Medium Parameters
% % =========================================================================
% 
% model_f0_m = ones(Nx, Ny) * 25;
% model_f0 = ones(Nx, Ny) * 25; 
% model_c0_m = ones(Nx, Ny) * 3600;   model_c0_m(1 : 130, :) = 1800; 
% model_Q = ones(Nx, Ny) * 100;       model_Q(1 : 130, :) = 30;
% model_density = ones(Nx, Ny) * 2200;
% 
% % =========================================================================
% % Source & Receivers
% % =========================================================================
% 
% f_rw_c = 25;       % Center frequency of ricker wavelet [Hz]
% x_src = 0;          % Source location in the x direction [m]
% y_src = 0;       % Source location in the y direction [m]
% x_rec = 400.;          % Receiver location in the x direction [m]
% y_rec = 0;        % Receiver location in the y direction [m]
% 
% % =========================================================================
% % Simulation Parameters
% % =========================================================================
% 
% dt = 1e-3;        % Time interval [s]
% t_max = 0.6;        % Simulation end time [s]
% f0 = 25;           % Reference frequency for simulation
% mod_mech = 'TF111110';  % Numerical modeling mechanism
% t_snap = 0.30;
% mat_name = 'case01_a';

% =========================================================================
% Simulation Set-up
% =========================================================================

% Grid
kgrid = kWaveGrid(Nx, dx, Ny, dy);
kgrid.t_array = 0 : dt: t_max;

% Medium
model_gamma= atan(1 ./ model_Q) / pi;
model_c0   = model_c0_m .* (model_f0 ./ model_f0_m) .^ model_gamma;

medium.f0 = model_f0;
medium.c0 = model_c0;
medium.Q  = model_Q;
medium.density = model_density;
medium.sound_speed = medium.c0;
medium.mod_mech = mod_mech;

% Source
[nx_src, ny_src] = close_grid_2d(kgrid, x_src, y_src);
source.p_mask = zeros(Nx, Ny);
source.p_mask(nx_src, ny_src) = 1;
source.p = stf;

% Sensor
[nx_rec, ny_rec] = close_grid_2d(kgrid, x_rec, y_rec);
sensor.mask = zeros(Nx, Ny);
sensor.mask(nx_rec, ny_rec) = 1;
if ~strcmp(mat_name, '') && t_snap >= 0
    sensor.t_snap = t_snap;
    sensor.snap_name = mat_name;
end

% =========================================================================
% RTM Q-compensated Forward modeling
% =========================================================================

medium.alpha_sign = [1, -1];
% medium.alpha_filter = getAlphaFilter(kgrid, medium, f_cutoff,taper_ratio);
medium.save_ndt = 1;
medium.save_pressure = true;


if exist('args', 'var')
    [d, p_save] = kspaceFirstOrder2D(kgrid, medium, source, sensor, args{:});
else
    [d, p_save] = kspaceFirstOrder2D(kgrid, medium, source, sensor);
end
% d = kspaceFirstOrder2D(kgrid, medium, source, sensor, 'PMLInside', false, 'PlotSim', false);
t_axis = kgrid.t_array;

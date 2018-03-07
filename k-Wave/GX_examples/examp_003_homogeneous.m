% Author: Guangchi Xing
% Date: 11/09/2017
% High frequency two stations

clear;

% =========================================================================
% Grid Parameters
% =========================================================================

Nx = 128;           % number of grid points in the x direction
Ny = 256;           % number of grid points in the y direction
dx = 0.5;           % grid point spacing in the x direction [m]
dy = 0.5;           % grid point spacing in the y direction [m]

% =========================================================================
% Medium Parameters
% =========================================================================

f0_m = 300;           % Reference frequency [Hz]
c0_m = 2089;          % Phase velocity at reference frequency [m/s]
Q = 5;             % Quality factor
density = 2200;     % Density [kg/m^3]

% =========================================================================
% Source & Receivers
% =========================================================================

f_rw_c = 200;       % Center frequency of ricker wavelet [Hz]
x_src = 0.;          % Source location in the x direction [m]
y_src = -34.5;       % Source location in the y direction [m]
x_rec = [0.; 0.];          % Receiver location in the x direction [m]
y_rec = [-14.5; 35.5];        % Receiver location in the y direction [m]

% =========================================================================
% Simulation Parameters
% =========================================================================

dt = 1e-5;        % Time interval [s]
t_max = 0.06;        % Simulation end time [s]
f0 = 200;           % Reference frequency for simulation

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

% medium.c0 = c0 * ones(Nx, Ny);
% medium.f0 = f0 * ones(Nx, Ny);
% medium.Q = Q * ones(Nx, Ny);
medium.f0 = model_f0;
medium.c0 = model_c0;
medium.Q  = model_Q;
medium.density = density * ones(Nx, Ny);
medium.sound_speed = medium.c0;

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

% source.p = source.p / c0;
d1 = lole_analytical_2d(kgrid, medium, source, sensor);
d2 = kjar_analytical_2d(kgrid, medium, source, sensor);

% =========================================================================
% Compute the Numerical Solution
% =========================================================================
% medium.mod_mech = 'TZ14';
% medium.mod_mech = 'TZ17';
medium.mod_mech = 'FT17';

d3 = kspaceFirstOrder2D(kgrid, medium, source, sensor);
d3 = d3 * 4;

% =========================================================================
% Post-Processing
% =========================================================================

[f_vec, d1_vec] = comp_spec(kgrid.t_array, d1);
[f_vec, d2_vec] = comp_spec(kgrid.t_array, d2);
[f_vec, d3_vec] = comp_spec(kgrid.t_array, d3);

r = sqrt((x_rec - x_src).^2 + (y_rec - y_src).^2);
alpha_meas1 = measure_alpha(d1_vec(1, :), d1_vec(2, :), r(1), r(2), abs(r(2) - r(1)));
alpha_meas2 = measure_alpha(d2_vec(1, :), d2_vec(2, :), r(1), r(2), abs(r(2) - r(1)));
alpha_meas3 = measure_alpha(d3_vec(1, :), d3_vec(2, :), r(1), r(2), abs(r(2) - r(1)));
cp_meas1 = measure_cp(f_vec, d1_vec(1, :), d1_vec(2, :), abs(r(2) - r(1)));
cp_meas2 = measure_cp(f_vec, d2_vec(1, :), d2_vec(2, :), abs(r(2) - r(1)));
cp_meas3 = measure_cp(f_vec, d3_vec(1, :), d2_vec(2, :), abs(r(2) - r(1)));

f_b = 0;
f_e = 500;
figure(1);
subplot(211);
% plot(f_vec, alpha_meas1, 'k', 'linewidth', 3); hold on;
plot(f_vec, alpha_meas2, 'b', 'linewidth', 2); hold on;
plot(f_vec, alpha_meas3, 'r--', 'linewidth', 2);
xlim([f_b f_e]);
subplot(212);
% plot(f_vec, cp_meas1, 'k', 'linewidth', 3); hold on;
plot(f_vec, cp_meas2, 'b', 'linewidth', 2); hold on;
plot(f_vec, cp_meas3, 'r--', 'linewidth', 2);
xlim([f_b f_e]);


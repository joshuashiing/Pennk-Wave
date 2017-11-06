% Author: Guangchi Xing
% Date: 10/19/2017

clear;

% =========================================================================
% Grid Parameters
% =========================================================================

Nx = 256;           % number of grid points in the x direction
Ny = 512;           % number of grid points in the y direction
dx = 1.0;           % grid point spacing in the x direction [m]
dy = 1.0;           % grid point spacing in the y direction [m]

% =========================================================================
% Medium Parameters
% =========================================================================

f0 = 100;           % Reference frequency [Hz]
c0 = 2164;          % Phase velocity at reference frequency [m/s]
Q = 132.5;             % Quality factor
density = 2200;     % Density [kg/m^3]

% =========================================================================
% Source & Receivers
% =========================================================================

f_rw_c = 100;       % Center frequency of ricker wavelet [Hz]
x_src = 0;          % Source location in the x direction [m]
y_src = -128;       % Source location in the y direction [m]
x_rec = 0;          % Receiver location in the x direction [m]
y_rec = 128;        % Receiver location in the y direction [m]

% =========================================================================
% Simulation Parameters
% =========================================================================

dt = 1.0e-4;        % Time interval [s]
t_max = 0.2;        % Simulation end time [s]

% =========================================================================
% Simulation Set-up
% =========================================================================

% Grid
kgrid = kWaveGrid(Nx, dx, Ny, dy);
kgrid.t_array = 0 : dt: t_max;

% Medium
medium.f0 = f0 * ones(Nx, Ny);
medium.c0 = c0 * ones(Nx, Ny);
medium.Q = Q * ones(Nx, Ny);
medium.density = density * ones(Nx, Ny);
medium.sound_speed = medium.c0;

% Source
[nx_src, ny_src] = close_grid_2d(kgrid, x_src, y_src);
source.p_mask = zeros(Nx, Ny);
source.p_mask(nx_src, ny_src) = 1;
t_rw_c = 1 / f_rw_c * 2;
nt_rw = 1 / f_rw_c * 4 / dt;
rw = rickerwavelet(f_rw_c, dt, nt_rw, t_rw_c, 1);
source.p = rw;

% Sensor
[nx_rec, ny_rec] = close_grid_2d(kgrid, x_rec, y_rec);
sensor.mask = zeros(Nx, Ny);
sensor.mask(nx_rec, ny_rec) = 1;


% =========================================================================
% Compute the Analytical Solution
% =========================================================================

d1 = kjar_analytical_2d(kgrid, medium, source, sensor);
d2 = lole_analytical_2d(kgrid, medium, source, sensor);
% medium2.sound_speed = c0 * ones(Nx, Ny);
% medium2.density = density * ones(Nx, Ny);
% d3 = kspaceFirstOrder2D(kgrid, medium, source, sensor);

plot(kgrid.t_array, d1, 'k', 'linewidth', 2); hold on;
plot(kgrid.t_array, d2, 'r--', 'linewidth', 2);
% plot(kgrid.t_array, d3, 'b--', 'linewidth', 2);
% plot(kgrid.t_array, sensor_data_analytical);

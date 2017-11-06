% Author: Guangchi Xing
% Date: 10/19/2017

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

f0 = 100;           % Reference frequency [Hz]
c0 = 2089;          % Phase velocity at reference frequency [m/s]
Q = 32.5;             % Quality factor
density = 2200;     % Density [kg/m^3]

% =========================================================================
% Source & Receivers
% =========================================================================

f_rw_c = 200;       % Center frequency of ricker wavelet [Hz]
x_src = 0.;          % Source location in the x direction [m]
y_src = -34.5;       % Source location in the y direction [m]
x_rec = 0.;          % Receiver location in the x direction [m]
y_rec = 34.5;        % Receiver location in the y direction [m]

% =========================================================================
% Simulation Parameters
% =========================================================================

dt = 1e-5;        % Time interval [s]
t_max = 0.06;        % Simulation end time [s]

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
source.p = rw*100;
% source.p = filterTimeSeries(kgrid, medium, source.p);

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

% source.p = source.p * c0;
% d3 = kspaceFirstOrder2D(kgrid, medium, source, sensor);
% medium.mod_mech = 'TZ14';
% medium.mod_mech = 'TZ17';
medium.mod_mech = 'TF17';

% medium.max_p_vec = zeros(size(kgrid.t_array));
% medium.p_test = zeros(size(kgrid.t_array));
d3 = kspaceFirstOrder2D(kgrid, medium, source, sensor);



figure;
plot(kgrid.t_array, d1, 'k', 'linewidth', 3); hold on;
plot(kgrid.t_array, d2, 'b--', 'linewidth', 3);
plot(kgrid.t_array, d3 * 4, 'r--', 'linewidth', 2);


% figure;
% plot(kgrid.t_array, d1 / max(d1), 'k', 'linewidth', 3); hold on;
% plot(kgrid.t_array, d2 / max(d2), 'b--', 'linewidth', 3);
% plot(kgrid.t_array, d3 / max(d3), 'r--', 'linewidth', 2);
% plot(kgrid.t_array, sensor_data_analytical);
%  xlim([0 0.2]);
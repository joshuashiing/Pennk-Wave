% Author: Guangchi Xing
% Date: 04/16/2018
% Demonstrate the two effects of seismic attenuation

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

f0_m = 200;           % Reference frequency [Hz]
c0_m = 2089;          % Phase velocity at reference frequency [m/s]
Q = 10;             % Quality factor
density = 2200;     % Density [kg/m^3]

% =========================================================================
% Source & Receivers
% =========================================================================

f_rw_c = 200;       % Center frequency of ricker wavelet [Hz]
x_src = 0.;          % Source location in the x direction [m]
y_src = -60;       % Source location in the y direction [m]
x_rec = 0.;          % Receiver location in the x direction [m]
y_rec = 30;        % Receiver location in the y direction [m]

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

n_rec = 10;
y_rec_list = linspace(y_src, y_rec, n_rec);
d_array = zeros(n_rec - 1, length(kgrid.t_array));

for i_rec = 2 : n_rec
    yi_rec = y_rec_list(i_rec);
    dist = yi_rec - y_src;
    [nx_rec, ny_rec] = close_grid_2d(kgrid, x_rec, yi_rec);
    sensor.mask = zeros(Nx, Ny);
    sensor.mask(nx_rec, ny_rec) = 1;
    d = kjar_analytical_2d(kgrid, medium, source, sensor);
    d_ps = -phase_shift(kgrid.t_array, d, pi/4);
    d_ps = d_ps * sqrt(dist);
    d_array(i_rec - 1, :) = d_ps;
%     plot(kgrid.t_array, d_ps + i_rec); hold on;
%     plot(kgrid.t_array, d_ps * sqrt(dist) + i_rec * 5); hold on;

end
save('visc', 'd_array');
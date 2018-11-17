% =========================================================================
% Test the smoothing parameters
% Date 11/13/2018
% Author: Guangchi Xing


clear;
clc

% =========================================================================
% Grid Parameters
% =========================================================================

Nx = 100;           % number of grid points in the x direction
Ny = 300;           % number of grid points in the y direction
dx = 10;            % grid point spacing in the x direction [m]
dy = 10;            % grid point spacing in the y direction [m]

x0 = 0;             % Top-left corner x coordinate [m]
y0 = 0;             % Top-left corner y coordinate [m]

% =========================================================================
% Correction from the original coordinate to the k-Wave coordinate
% =========================================================================
kgrid_tmp = kWaveGrid(Nx, dx, Ny, dy);
xc = kgrid_tmp.x_vec(1) - x0;   % Correction for x direction
yc = kgrid_tmp.y_vec(1) - y0;   % Correction for y direction

% =========================================================================
% Simulation Parameters
% =========================================================================

dt = 1e-3;        % Time interval [s]
t_max = 1.6;        % Simulation end time [s]
f0 = 25;           % Reference frequency for simulation
args = {'PMLInside', false, 'PlotSim', false};

% =========================================================================
% Source & Receivers
% =========================================================================

y_src = (5 : 10 : Ny)' * dy;
x_src = ones(size(y_src)) * 5 * dx;
y_rec = (5 : 2 : Ny)' * dy;
x_rec = ones(size(y_rec)) * 5 * dx;

x_src = x_src + xc;
y_src = y_src + yc;
x_rec = x_rec + xc;
y_rec = y_rec + yc;

f_rw_c = 25;        % Center frequency of ricker wavelet [Hz]
A_rw = 1e2;         % Amplitude of the ricker wavelet
t_rw_c = 1 / f_rw_c * 2;
nt_rw = 1 / f_rw_c * 4 / dt;
stf = rickerwavelet(f_rw_c, dt, nt_rw, t_rw_c, 1) * A_rw;

% =========================================================================
% Medium Parameters
% =========================================================================

vp = ones(Nx, Ny) * 2600;
vp(1 : 50, :) = 1800;
rho = ones(Nx, Ny) * 2200;
rho(1 : 50, :) = 1800;
Q = ones(Nx, Ny) * 100;
Q(1 : 50, :) = 30;

f0_model = ones(Nx, Ny) * 100;
f0 = ones(Nx, Ny) * f0;

sigma_list = [2, 4, 8, 16];
% c0 = imgaussfilt(c0, 4);
% rho = imgaussfilt(rho, 4);
% =========================================================================
% Run!!!
% =========================================================================

mod_mech = 'TF111110';
% mod_mech = 'lossless';

i = 10;

[d0, t_axis] = heter_simu(Nx, Ny, dx, dy, f0_model, f0, vp, Q, rho, ...
                    stf, x_src(i), y_src(i), x_rec, y_rec, ...
                    dt, t_max, mod_mech, -1, '', args);
d = zeros(size(d0, 1), size(d0, 2), length(sigma_list));
for j = 1 : length(sigma_list)
    sigma = sigma_list(j);
    vp_sm = imgaussfilt(vp, sigma);
    Q_sm = imgaussfilt(Q, sigma);
    rho_sm = imgaussfilt(rho, sigma);
    [d(:, :, j), t_axis] = heter_simu(Nx, Ny, dx, dy, f0_model, f0, vp_sm, Q_sm, rho_sm, ...
                    stf, x_src(i), y_src(i), x_rec, y_rec, ...
                    dt, t_max, mod_mech, -1, '', args);
end

r = 100;
plot(t_axis, d0(r, :), 'k', 'linewidth', 2); hold on;
clist = {'b', 'r', 'c', 'g', 'y'};
for j = 1 : length(sigma_list)
    plot(t_axis, d(r, :, j), [clist{j} '--'], 'linewidth', 2);
end

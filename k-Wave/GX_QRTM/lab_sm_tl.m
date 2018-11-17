% =========================================================================
% Example to develop Q-compensated RTM
% Date 11/01/2018
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

% dt = 1e-3;        % Time interval [s]
dt = 1e-3;        % Time interval [s]
t_max = 1.6;        % Simulation end time [s]
f0 = 15;           % Reference frequency for simulation
% args = {'PMLInside', false, 'PlotSim', false};
args = {'PMLInside', false};

% =========================================================================
% Source & Receivers
% =========================================================================

y_src = (5 : 10 : Ny)' * dy;
x_src = ones(size(y_src)) * 5 * dx;
% x_src = ones(size(y_src)) * 25 * dx;
y_rec = (5 : 2 : Ny)' * dy;
x_rec = ones(size(y_rec)) * 5 * dx;

x_src = x_src + xc;
y_src = y_src + yc;
x_rec = x_rec + xc;
y_rec = y_rec + yc;

f_rw_c = 15;        % Center frequency of ricker wavelet [Hz]
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
% Q(:) = 999999;

vp1 = vp;
rho1 = rho;
Q1 = Q;

% vp2 = vp1;
% vp2(51:end, :) = 3600;
% rho2 = rho1;
% Q2 = Q1;

sigma = 8;
vp2 = imgaussfilt(vp, sigma);
rho2 = imgaussfilt(rho, sigma);
Q2 = imgaussfilt(Q, sigma);


vp_tl = mean(mean(vp(1:3, :)));
rho_tl = mean(mean(rho(1:3, :)));
Q_tl = mean(mean(Q(1:3, :)));
vp3(:) = vp_tl;
rho3(:) = rho_tl;
Q3(:) = Q_tl;

f0_model = ones(Nx, Ny) * 100;
f0 = ones(Nx, Ny) * f0;

% =========================================================================
% Run!!!
% =========================================================================

mod_mech = 'TF111110';
% mod_mech = 'lossless';

i = 10;

[d1, t_axis] = heter_simu(Nx, Ny, dx, dy, f0_model, f0, vp1, Q1, rho1, ...
                    stf, x_src(i), y_src(i), x_rec, y_rec, ...
                    dt, t_max, mod_mech, -1, '', args);

[d2, t_axis] = heter_simu(Nx, Ny, dx, dy, f0_model, f0, vp2, Q2, rho2, ...
                    stf, x_src(i), y_src(i), x_rec, y_rec, ...
                    dt, t_max, mod_mech, -1, '', args);
                
[d3, t_axis] = heter_simu(Nx, Ny, dx, dy, f0_model, f0, vp3, Q3, rho3, ...
                    stf, x_src(i), y_src(i), x_rec, y_rec, ...
                    dt, t_max, mod_mech, -1, '', args);

d1 = d1';
d2 = d2';
d3 = d3';

clim = max(d1(:)) * 0.02;

figure(1);
subplot(131); imagesc(d1, [-clim, clim]);
subplot(132); imagesc(d1 - d2, [-clim, clim]);
subplot(133); imagesc(d1 - d3, [-clim, clim]);

figure(2);
i = 80;
plot(t_axis, d1(:, i), 'k', 'linewidth', 2); hold on;
plot(t_axis, d2(:, i), 'b--', 'linewidth', 2);
plot(t_axis, d3(:, i), 'r--', 'linewidth', 2);
% =========================================================================
% Example to develop Q-compensated RTM
% Date 11/01/2018
% Author: Guangchi Xing


clear;
clc

load('GC_model/BP_cut_model_gas_chimney');

% =========================================================================
% Grid Parameters
% =========================================================================

Nx = nz;           % number of grid points in the x direction
Ny = nx;           % number of grid points in the y direction
% dx = 12.5;            % grid point spacing in the x direction [m]
% dy = 12.5;            % grid point spacing in the y direction [m]
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

dt = 0.5e-3;        % Time interval [s]
t_max = 2;        % Simulation end time [s]
f0 = 15;           % Reference frequency for simulation
% args = {'PMLInside', false, 'PlotSim', false};
args = {'PMLInside', false};

% =========================================================================
% Source & Receivers
% =========================================================================

y_src = (3 : 10 : (Ny-2))' * dy;
x_src = ones(size(y_src)) * 2 * dx;
y_rec = (3 : 2 : (Ny-2))' * dy;
x_rec = ones(size(y_rec)) * 2 * dx;

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

c0 = velp;
rho = ones(Nx, Ny) * 2200;
Q = Qp;

f0_model = ones(Nx, Ny) * 100;
f0 = ones(Nx, Ny) * f0;

Q(:) = 999999;

c01 = c0;
rho1 = rho;
Q1 = Q;

sigma = 8;
c02 = imgaussfilt(c0, sigma);
rho2 = imgaussfilt(rho, sigma);
Q2 = imgaussfilt(Q, sigma);


c0_tl = mean(mean(c0(1:3, :)));
rho_tl = mean(mean(rho(1:3, :)));
Q_tl = mean(mean(Q(1:3, :)));
c03(:) = c0_tl;
rho3(:) = rho_tl;
Q3(:) = Q_tl;

% =========================================================================
% Run!!!
% =========================================================================

% mod_mech = 'TF111110';
mod_mech = 'lossless';

i = 10;

[d1, t_axis] = heter_simu(Nx, Ny, dx, dy, f0_model, f0, c01, Q1, rho1, ...
                    stf, x_src(i), y_src(i), x_rec, y_rec, ...
                    dt, t_max, mod_mech, -1, '', args);

% [d2, t_axis] = heter_simu(Nx, Ny, dx, dy, f0_model, f0, c02, Q2, rho2, ...
%                     stf, x_src(i), y_src(i), x_rec, y_rec, ...
%                     dt, t_max, mod_mech, -1, '', args);
                
[d3, t_axis] = heter_simu(Nx, Ny, dx, dy, f0_model, f0, c03, Q3, rho3, ...
                    stf, x_src(i), y_src(i), x_rec, y_rec, ...
                    dt, t_max, mod_mech, -1, '', args);

d1 = d1';
% d2 = d2';
d3 = d3';

clim = max(d1(:)) * 0.02;

figure(1);
subplot(131); imagesc(d1, [-clim, clim]);
% subplot(132); imagesc(d1 - d2, [-clim, clim]);
subplot(133); imagesc(d1 - d3, [-clim, clim]);

figure(2);
i = 80;
plot(t_axis, d1(:, i), 'k', 'linewidth', 2); hold on;
% plot(t_axis, d2(:, i), 'b--', 'linewidth', 2);
plot(t_axis, d3(:, i), 'r--', 'linewidth', 2);
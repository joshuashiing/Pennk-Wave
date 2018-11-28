% =========================================================================
% Homogeneous test
% Date 11/27/2018
% Author: Guangchi Xing


clear;
clc

% load('GC_model/BP_cut_model_gas_chimney');

% =========================================================================
% Grid Parameters
% =========================================================================

Nx = 250;           % number of grid points in the x direction
Ny = 300;           % number of grid points in the y direction
dx = 2;            % grid point spacing in the x direction [m]
dy = 2;            % grid point spacing in the y direction [m]

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

dt = 0.2e-3;        % Time interval [s]
t_max = 0.35;        % Simulation end time [s]
f0 = 35;           % Reference frequency for simulation
% args = {'PMLInside', false, 'PlotSim', false};
args = {'PMLInside', false};

% =========================================================================
% Source & Receivers
% =========================================================================

x_src = 250;
y_src = 100;
x_rec = 250;
y_rec = 500;

x_src = x_src + xc;
y_src = y_src + yc;
x_rec = x_rec + xc;
y_rec = y_rec + yc;

f_rw_c = 35;        % Center frequency of ricker wavelet [Hz]
A_rw = 1e3;         % Amplitude of the ricker wavelet
t_rw_c = 1 / f_rw_c * 2;
nt_rw = 1 / f_rw_c * 4 / dt;
stf = rickerwavelet(f_rw_c, dt, nt_rw, t_rw_c, 1) * A_rw;

% =========================================================================
% Medium Parameters
% =========================================================================

c0 = ones(Nx, Ny) * 2164;
f0_model = ones(Nx, Ny) * 1500;
Q = ones(Nx, Ny) * 10;
rho = ones(Nx, Ny) * 2200;
f0 = ones(Nx, Ny) * f0;

% c0 = velp;
% rho = ones(Nx, Ny) * 2200;
% Q = Qp;
% 
% f0_model = ones(Nx, Ny) * 100;
% f0 = ones(Nx, Ny) * f0;

% =========================================================================
% Run!!!
% =========================================================================

mod_mech = 'TF111110';
% mod_mech = 'lossless';

[d, d2, t_axis] = homo_simu(Nx, Ny, dx, dy, f0_model, f0, c0, Q, rho, ...
                    stf, x_src, y_src, x_rec, y_rec, ...
                    dt, t_max, mod_mech, -1, '', args);
% d2 = kjar_analytical_2d(kgrid, medium, source, sensor);
% figure(3);
plot(t_axis, d2 * 16, 'k', 'linewidth', 2); hold on;
plot(t_axis, d, 'r--', 'linewidth', 2);
xlim([0.15 0.35]);

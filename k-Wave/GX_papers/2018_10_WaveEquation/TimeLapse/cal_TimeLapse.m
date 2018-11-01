clear;
clc;

load('conversion_bs03j_V1.mat');

% =========================================================================
% Grid Parameters
% =========================================================================

Nx = 434;           % number of grid points in the x direction
Ny = 467;           % number of grid points in the y direction
dx = 0.15;          % grid point spacing in the x direction [m]
dy = 0.15;          % grid point spacing in the y direction [m]

x0 = 1621.45;       % Top-left corner x coordinate [m]
y0 = -19.95;        % Top-left corner y coordinate [m]

% =========================================================================
% Correction from the original coordinate to the k-Wave coordinate
% =========================================================================
kgrid_tmp = kWaveGrid(Nx, dx, Ny, dy);
xc = kgrid_tmp.x_vec(1) - x0;   % Correction for x direction
yc = kgrid_tmp.y_vec(1) - y0;   % Correction for y direction

% =========================================================================
% Source & Receivers
% =========================================================================

f_rw_c = 1000;                               % Center frequency of ricker wavelet [Hz]
% x_src = 1658.05;                            % Source location x in original coordinate [m]
x_src = 1650.10;
y_src = 0;                                  % Source location y in original coordinate [m]
% x_rec = [1635.1; 1650.1; 1658.05; 1680.1];	% Receiver location x in original coordinate [m]
% y_rec = [30; 30; 30; 30];                   % Receiver location y in original coordinate [m]
x_rec = (1635.1 : 0.3 : 1680.1)';
y_rec = 30 * ones(size(x_rec));

x_src = x_src + xc;
y_src = y_src + yc;
x_rec = x_rec + xc;
y_rec = y_rec + yc;

% =========================================================================
% Simulation Parameters
% =========================================================================

dt = 2e-5;        % Time interval [s]
t_max = 3e-2;        % Simulation end time [s]
% t_snap = 0.8;

% =========================================================================
% Media Parameters
% =========================================================================

Q_max = 1000;
f0_m = 5000;
t_snap = 1.6e-2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% case 1
% Time slice 01
mat_name = 'case01'
it = 1;

vp  = cell2mat(seisProps.vp(it));
vs  = cell2mat(seisProps.vs(it));
rho = cell2mat(seisProps.rho(it));
Qi  = cell2mat(seisProps.atten(it));
Q = 1 ./ Qi;
Q(isinf(Q)) = Q_max;

f0 = f_rw_c * ones(size(vp));
f0_model = f0_m * ones(size(vp));

mod_mech = 'TF111110';
[t_axis, d] = heter_simu(Nx, Ny, dx, dy, f0_model, f0, vp, Q, rho,...
                        f_rw_c, x_src, y_src, x_rec, y_rec, ...
                        dt, t_max, mod_mech, t_snap, mat_name);
save(mat_name, 'vp', 'rho', 'Q', 't_axis', 'd', '-append');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% case 2
% Time slice 21
mat_name = 'case02'
it = 5;

vp  = cell2mat(seisProps.vp(it));
vs  = cell2mat(seisProps.vs(it));
rho = cell2mat(seisProps.rho(it));
Qi  = cell2mat(seisProps.atten(it));
Q = 1 ./ Qi;
Q(isinf(Q)) = Q_max;

f0 = f_rw_c * ones(size(vp));
f0_model = f0_m * ones(size(vp));

mod_mech = 'TF111110';
[t_axis, d] = heter_simu(Nx, Ny, dx, dy, f0_model, f0, vp, Q, rho,...
                        f_rw_c, x_src, y_src, x_rec, y_rec, ...
                        dt, t_max, mod_mech, t_snap, mat_name);
save(mat_name, 'vp', 'rho', 'Q', 't_axis', 'd', '-append');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% case 3
% Time slice 41
mat_name = 'case03'
it = 41;

vp  = cell2mat(seisProps.vp(it));
vs  = cell2mat(seisProps.vs(it));
rho = cell2mat(seisProps.rho(it));
Qi  = cell2mat(seisProps.atten(it));
Q = 1 ./ Qi;
Q(isinf(Q)) = Q_max;

f0 = f_rw_c * ones(size(vp));
f0_model = f0_m * ones(size(vp));

mod_mech = 'TF111110';
[t_axis, d] = heter_simu(Nx, Ny, dx, dy, f0_model, f0, vp, Q, rho,...
                        f_rw_c, x_src, y_src, x_rec, y_rec, ...
                        dt, t_max, mod_mech, t_snap, mat_name);
save(mat_name, 'vp', 'rho', 'Q', 't_axis', 'd', '-append');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% case 4
% Time slice 41 lossless
mat_name = 'case04'
it = 41;

vp  = cell2mat(seisProps.vp(it));
vs  = cell2mat(seisProps.vs(it));
rho = cell2mat(seisProps.rho(it));
Qi  = cell2mat(seisProps.atten(it));
Q = 1 ./ Qi;
Q(isinf(Q)) = Q_max;

f0 = f_rw_c * ones(size(vp));

mod_mech = 'lossless';
[t_axis, d] = heter_simu(Nx, Ny, dx, dy, f0, f0, vp, Q, rho,...
                        f_rw_c, x_src, y_src, x_rec, y_rec, ...
                        dt, t_max, mod_mech, t_snap, mat_name);
save(mat_name, 'vp', 'rho', 'Q', 't_axis', 'd', '-append');

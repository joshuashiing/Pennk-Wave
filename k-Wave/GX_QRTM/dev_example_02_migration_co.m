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
dx = 12.5;            % grid point spacing in the x direction [m]
dy = 12.5;            % grid point spacing in the y direction [m]

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
t_max = 2;        % Simulation end time [s]
f0 = 15;           % Reference frequency for simulation
args = {'PMLInside', false, 'PlotSim', false};
% args = {'PMLInside', false};

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

c0 = imgaussfilt(c0, 4);
rho = imgaussfilt(rho, 4);

f0_model = ones(Nx, Ny) * 100;
f0 = ones(Nx, Ny) * f0;
% =========================================================================
% Run!!!
% =========================================================================

mod_mech = 'TF111110';
% mod_mech = 'lossless';
f_cutoff = 50;
taper_ratio = 0.5;

    
% [d, p_save_fm] = rtm_fm_simu(Nx, Ny, dx, dy, f0_model, f0, vp, Q, rho, ...
%                     stf, x_src(15), y_src(15), x_rec, y_rec, ...
%                     dt, t_max, mod_mech, f_cutoff, taper_ratio, -1, '', args);
% 
% tmp = load('Data_example_01_HQ/CSG_015.mat');
% p_data = tmp.d;
% clear tmp
% [d, p_save_bp] = rtm_bp_simu(Nx, Ny, dx, dy, f0_model, f0, vp, Q, rho, ...
%                     x_rec, y_rec, ...
%                     dt, t_max, mod_mech, f_cutoff, taper_ratio, -1, '', p_data, args);


% Migration
% mig1 = zeros(Nx, Ny, 30);
mig1 = zeros(Nx, Ny, length(x_src));
mig2 = zeros(Nx, Ny, length(x_src));
for j = 1 : length(x_src)
    fprintf(['Working on shot #', num2str(j, '%.3i'), '\n']);
    wffile = ['Data_example_02/WF_co_' num2str(j, '%.3i') '.mat'];
    wf = load(wffile);
    p_save_fm = wf.p_save_fm;
    p_save_bp = wf.p_save_bp;
    n = size(p_save_fm, 2);
    p_save_fm = reshape(p_save_fm, Nx, Ny, n);
    p_save_bp = reshape(p_save_bp, Nx, Ny, n - 1);
    
    for i = 1 : (n - 1)
        mig1(:, :, j) = p_save_fm(:, :, i) .* p_save_bp(:, :, (n-i)) + mig1(:, :, j);
        mig2(:, :, j) = del2(p_save_fm(:, :, i)) .* p_save_bp(:, :, (n-i)) + mig2(:, :, j);
    end
end

save('Data_example_02/mig_co.mat', 'mig1', 'mig2');

% subplot(121); imagesc(mig1); colorbar;
% subplot(122); imagesc(mig2); colorbar;

% subplot(121); imagesc(mig1(20:end, :)); colorbar;
% subplot(122); imagesc(mig2(20:end, :)); colorbar;
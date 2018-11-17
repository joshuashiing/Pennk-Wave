% =========================================================================
% Example to develop Q-compensated RTM
% Date 11/13/2018
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

dt = 0.25e-3;        % Time interval [s]
t_max = 2;        % Simulation end time [s]
f0 = 10;           % Reference frequency for simulation
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

f_rw_c = 10;        % Center frequency of ricker wavelet [Hz]
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

c0 = imgaussfilt(c0, 8);
rho = imgaussfilt(rho, 8);
Q = imgaussfilt(Q, 8);
Q(:) = 999999;

f0_model = ones(Nx, Ny) * 100;
f0 = ones(Nx, Ny) * f0;

% =========================================================================
% Run!!!
% =========================================================================

mod_mech = 'TF111110';
% mod_mech = 'lossless';
f_cutoff = 1000;
taper_ratio = 0.5;

ex_name = 'Data_example_02';

mig1 = zeros(Nx, Ny, length(x_src));
mig2 = zeros(Nx, Ny, length(x_src));

for i = 1 : length(x_src)
% for i = 15 : 15

    clear p_save_fm p_save_bp
    fprintf(['Working on shot #', num2str(i, '%.3i'), '\n']);
    
    [d, p_save_fm] = rtm_fm_simu(Nx, Ny, dx, dy, f0_model, f0, c0, Q, rho, ...
                    stf, x_src(i), y_src(i), x_rec, y_rec, ...
                    dt, t_max, mod_mech, f_cutoff, taper_ratio, -1, '', args);
    
    datafile = [ex_name '/CSG_ac_' num2str(i, '%.3i') '.mat'];
    tmp1 = load(datafile);
    datafile = [ex_name '/CSG_ac_tl_' num2str(i, '%.3i') '.mat'];
    tmp2 = load(datafile);
    p_data = tmp1.d - tmp2.d;
    clear tmp
    
    [d, p_save_bp] = rtm_bp_simu(Nx, Ny, dx, dy, f0_model, f0, c0, Q, rho, ...
                    x_rec, y_rec, ...
                    dt, t_max, mod_mech, f_cutoff, taper_ratio, -1, '', p_data, args);
                
    n = size(p_save_fm, 2);
    p_save_fm = reshape(p_save_fm, Nx, Ny, n);
    p_save_bp = reshape(p_save_bp, Nx, Ny, n - 1);
    
    for it = 1 : (n - 1)
        mig1(:, :, i) = p_save_fm(:, :, it) .* p_save_bp(:, :, (n-it)) + mig1(:, :, i);
        mig2(:, :, i) = del2(p_save_fm(:, :, it)) .* p_save_bp(:, :, (n-it)) + mig2(:, :, i);
    end
    
end

mig_file = [ex_name '/mig_cd_rf.mat'];
save(mig_file, 'mig1', 'mig2');


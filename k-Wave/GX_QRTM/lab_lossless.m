% =========================================================================
% Example to develop Q-compensated RTM
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
dt = 0.5e-3;        % Time interval [s]
t_max = 1.6;        % Simulation end time [s]
f0 = 15;           % Reference frequency for simulation
% args = {'PMLInside', false, 'PlotSim', false};
args = {'PMLInside', false};

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
Q(:) = 999999;

c0 = imgaussfilt(vp, 8);
rho = imgaussfilt(rho, 8);
Q = imgaussfilt(Q, 8);

f0_model = ones(Nx, Ny) * 100;
f0 = ones(Nx, Ny) * f0;

% =========================================================================
% Run!!!
% =========================================================================

ex_name = 'Data_example_01';

% mig1 = zeros(Nx, Ny, length(x_src));
% mig2 = zeros(Nx, Ny, length(x_src));

i = 10;
scale = 0.01;

f_cutoff = 1e3;
taper_ratio = 0.2;

% mod_mech = 'TF111110';
% [~, p_save_fm1] = rtm_fm_simu(Nx, Ny, dx, dy, f0_model, f0, c0, Q, rho, ...
%                     stf, x_src(i), y_src(i), x_rec, y_rec, ...
%                     dt, t_max, mod_mech, f_cutoff, taper_ratio, -1, '', args);
% 
% mod_mech = 'lossless';
% [~, p_save_fm2] = rtm_fm_simu(Nx, Ny, dx, dy, f0_model, f0, c0, Q, rho, ...
%                     stf, x_src(i), y_src(i), x_rec, y_rec, ...
%                     dt, t_max, mod_mech, f_cutoff, taper_ratio, -1, '', args);
% 
%                 
% % f_cutoff = 1e1;
% % taper_ratio = 0.2;
% % [~, p_save_fm3] = rtm_fm_simu(Nx, Ny, dx, dy, f0_model, f0, c0, Q, rho, ...
% %                     stf, x_src(i), y_src(i), x_rec, y_rec, ...
% %                     dt, t_max, mod_mech, f_cutoff, taper_ratio, -1, '', args);
% 
% 
% n = size(p_save_fm1, 2);
% p_save_fm1 = reshape(p_save_fm1, Nx, Ny, n);
% p_save_fm2 = reshape(p_save_fm2, Nx, Ny, n);
% % p_save_fm3 = reshape(p_save_fm3, Nx, Ny, n);
% 
% 
% ns = 5;
% nc = 2;
% for is = 1 : ns
%     it = ceil(n / (ns + 1) * is);
%     
%     tmp = p_save_fm1(:, :, it);
%     clim = max(abs(tmp(:))) * scale;
%     
%     subplot(ns, nc, (is-1)*nc + 1);
%     imagesc(p_save_fm1(:, :, it), [-clim clim]); colorbar;
%     
%     subplot(ns, nc, (is-1)*nc + 2);
% %     imagesc(p_save_fm2(:, :, it), [-clim clim]); colorbar;
%     imagesc(p_save_fm2(:, :, it) - p_save_fm1(:, :, it), [-clim clim]); colorbar;
%     
% %     subplot(ns, nc, (is-1)*nc + 3);
% %     imagesc(p_save_fm3(:, :, it), [-clim clim]); colorbar;
% % %     imagesc(p_save_fm3(:, :, it) - p_save_fm1(:, :, it), [-clim clim]); colorbar;
% end


datafile = [ex_name '/CSG_' num2str(i, '%.3i') '.mat'];
tmp1 = load(datafile);
datafile = [ex_name '/CSG_tl_' num2str(i, '%.3i') '.mat'];
tmp2 = load(datafile);
p_data = tmp1.d - tmp2.d;
clear tmp

mod_mech = 'TF111110';
[d, p_save_bp1] = rtm_bp_simu(Nx, Ny, dx, dy, f0_model, f0, c0, Q, rho, ...
                    x_rec, y_rec, ...
                    dt, t_max, mod_mech, f_cutoff, taper_ratio, -1, '', p_data, args);

mod_mech = 'lossless';
[d, p_save_bp2] = rtm_bp_simu(Nx, Ny, dx, dy, f0_model, f0, c0, Q, rho, ...
                    x_rec, y_rec, ...
                    dt, t_max, mod_mech, f_cutoff, taper_ratio, -1, '', p_data, args);


n = size(p_save_bp1, 2);
p_save_bp1 = reshape(p_save_bp1, Nx, Ny, n);
p_save_bp2 = reshape(p_save_bp2, Nx, Ny, n);
% p_save_fm3 = reshape(p_save_fm3, Nx, Ny, n);


ns = 5;
nc = 2;
for is = 1 : ns
    it = ceil(n / (ns + 1) * is);
    
    tmp = p_save_bp1(:, :, it);
    clim = max(abs(tmp(:))) * scale;
    
    subplot(ns, nc, (is-1)*nc + 1);
    imagesc(p_save_bp1(:, :, it), [-clim clim]); colorbar;
    
    subplot(ns, nc, (is-1)*nc + 2);
%     imagesc(p_save_fm2(:, :, it), [-clim clim]); colorbar;
    imagesc(p_save_bp2(:, :, it) - p_save_bp1(:, :, it), [-clim clim]); colorbar;
    
%     subplot(ns, nc, (is-1)*nc + 3);
%     imagesc(p_save_fm3(:, :, it), [-clim clim]); colorbar;
% %     imagesc(p_save_fm3(:, :, it) - p_save_fm1(:, :, it), [-clim clim]); colorbar;
end


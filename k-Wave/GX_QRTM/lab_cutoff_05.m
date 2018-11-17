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

dt = 1e-3;        % Time interval [s]
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

c0 = imgaussfilt(c0, 8);
rho = imgaussfilt(rho, 8);
% Q = imgaussfilt(Q, 8);

f0_model = ones(Nx, Ny) * 100;
f0 = ones(Nx, Ny) * f0;

% =========================================================================
% Run!!!
% =========================================================================

mod_mech = 'TF111110';
% mod_mech = 'lossless';

ex_name = 'Data_example_05';

i = 10;


datafile = [ex_name '/CSG_' num2str(i, '%.3i') '.mat'];
tmp1 = load(datafile);
datafile = [ex_name '/CSG_tl_' num2str(i, '%.3i') '.mat'];
tmp2 = load(datafile);
p_data = tmp1.d - tmp2.d;
clear tmp

f_cutoff = 1e3;
taper_ratio = 0.2;
% [~, p_save_fm1] = rtm_fm_simu_dev(Nx, Ny, dx, dy, f0_model, f0, c0, Q, rho, ...
%                     stf, x_src(i), y_src(i), x_rec, y_rec, ...
%                     dt, t_max, mod_mech, f_cutoff, taper_ratio, -1, '', args);
[d, p_save_bp1] = rtm_bp_simu_dev(Nx, Ny, dx, dy, f0_model, f0, c0, Q, rho, ...
                    x_rec, y_rec, ...
                    dt, t_max, mod_mech, f_cutoff, taper_ratio, -1, '', p_data, args);
                

f_cutoff = 1e3;
taper_ratio = 0.2;
[~, p_save_fm1] = rtm_fm_simu(Nx, Ny, dx, dy, f0_model, f0, c0, Q, rho, ...
                    stf, x_src(i), y_src(i), x_rec, y_rec, ...
                    dt, t_max, mod_mech, f_cutoff, taper_ratio, -1, '', args);
[d, p_save_bp1] = rtm_bp_simu(Nx, Ny, dx, dy, f0_model, f0, c0, Q, rho, ...
                    x_rec, y_rec, ...
                    dt, t_max, mod_mech, f_cutoff, taper_ratio, -1, '', p_data, args);

f_cutoff = 1e2;
taper_ratio = 0.2;
[~, p_save_fm2] = rtm_fm_simu(Nx, Ny, dx, dy, f0_model, f0, c0, Q, rho, ...
                    stf, x_src(i), y_src(i), x_rec, y_rec, ...
                    dt, t_max, mod_mech, f_cutoff, taper_ratio, -1, '', args);
[d, p_save_bp2] = rtm_bp_simu(Nx, Ny, dx, dy, f0_model, f0, c0, Q, rho, ...
                    x_rec, y_rec, ...
                    dt, t_max, mod_mech, f_cutoff, taper_ratio, -1, '', p_data, args);
                
f_cutoff = 5e1;
taper_ratio = 0.2;
[~, p_save_fm3] = rtm_fm_simu(Nx, Ny, dx, dy, f0_model, f0, c0, Q, rho, ...
                    stf, x_src(i), y_src(i), x_rec, y_rec, ...
                    dt, t_max, mod_mech, f_cutoff, taper_ratio, -1, '', args);
[d, p_save_bp3] = rtm_bp_simu(Nx, Ny, dx, dy, f0_model, f0, c0, Q, rho, ...
                    x_rec, y_rec, ...
                    dt, t_max, mod_mech, f_cutoff, taper_ratio, -1, '', p_data, args);
                
f_cutoff = 1e1;
taper_ratio = 0.2;
[~, p_save_fm4] = rtm_fm_simu(Nx, Ny, dx, dy, f0_model, f0, c0, Q, rho, ...
                    stf, x_src(i), y_src(i), x_rec, y_rec, ...
                    dt, t_max, mod_mech, f_cutoff, taper_ratio, -1, '', args);
[d, p_save_bp4] = rtm_bp_simu(Nx, Ny, dx, dy, f0_model, f0, c0, Q, rho, ...
                    x_rec, y_rec, ...
                    dt, t_max, mod_mech, f_cutoff, taper_ratio, -1, '', p_data, args);
                
n = size(p_save_fm1, 2);
p_save_fm1 = reshape(p_save_fm1, Nx, Ny, n);
p_save_fm2 = reshape(p_save_fm2, Nx, Ny, n);
p_save_fm3 = reshape(p_save_fm3, Nx, Ny, n);
p_save_fm4 = reshape(p_save_fm4, Nx, Ny, n);

p_save_bp1 = reshape(p_save_bp1, Nx, Ny, n-1);
p_save_bp2 = reshape(p_save_bp2, Nx, Ny, n-1);
p_save_bp3 = reshape(p_save_bp3, Nx, Ny, n-1);
p_save_bp4 = reshape(p_save_bp4, Nx, Ny, n-1);

mig1 = zeros(Nx, Ny);
mig2 = zeros(Nx, Ny);
mig3 = zeros(Nx, Ny);
mig4 = zeros(Nx, Ny);

for it = 1 : (n - 1)
    mig1 = del2(p_save_fm1(:, :, it)) .* p_save_bp1(:, :, (n-it)) + mig1;
    mig2 = del2(p_save_fm2(:, :, it)) .* p_save_bp2(:, :, (n-it)) + mig2;
    mig3 = del2(p_save_fm3(:, :, it)) .* p_save_bp3(:, :, (n-it)) + mig3;
    mig4 = del2(p_save_fm4(:, :, it)) .* p_save_bp4(:, :, (n-it)) + mig4;
end

% plot forward wavefield
figure(1);
scale = 0.1;
ns = 5;
nc = 4;
for is = 1 : ns
    it = ceil(n / (ns + 1) * is);
    
    tmp = p_save_fm1(:, :, it);
    clim = max(abs(tmp(:))) * scale;
    clim = 0.5;
    
    subplot(ns, nc, (is-1)*nc + 1);
    imagesc(p_save_fm1(:, :, it), [-clim clim]); colorbar;
    
    subplot(ns, nc, (is-1)*nc + 2);
    imagesc(p_save_fm2(:, :, it), [-clim clim]); colorbar;
%     imagesc(p_save_fm2(:, :, it) - p_save_fm1(:, :, it), [-clim clim]); colorbar;
    
    subplot(ns, nc, (is-1)*nc + 3);
    imagesc(p_save_fm3(:, :, it), [-clim clim]); colorbar;
%     imagesc(p_save_fm3(:, :, it) - p_save_fm1(:, :, it), [-clim clim]); colorbar;
    
    subplot(ns, nc, (is-1)*nc + 4);
    imagesc(p_save_fm4(:, :, it), [-clim clim]); colorbar;
%     imagesc(p_save_fm4(:, :, it) - p_save_fm1(:, :, it), [-clim clim]); colorbar;
end

% plot backward wavefield
figure(2);
scale = 0.1;
for is = 1 : ns
    it = ceil(n / (ns + 1) * is);
    
    tmp = p_save_bp1(:, :, it);
    clim = max(abs(tmp(:))) * scale;
    clim = 0.5;
    
    subplot(ns, nc, (is-1)*nc + 1);
    imagesc(p_save_bp1(:, :, it), [-clim clim]); colorbar;
    
    subplot(ns, nc, (is-1)*nc + 2);
    imagesc(p_save_bp2(:, :, it), [-clim clim]); colorbar;
%     imagesc(p_save_bp2(:, :, it) - p_save_bp1(:, :, it), [-clim clim]); colorbar;
    
    subplot(ns, nc, (is-1)*nc + 3);
    imagesc(p_save_bp3(:, :, it), [-clim clim]); colorbar;
%     imagesc(p_save_bp3(:, :, it) - p_save_bp1(:, :, it), [-clim clim]); colorbar;
    
    subplot(ns, nc, (is-1)*nc + 4);
    imagesc(p_save_bp4(:, :, it), [-clim clim]); colorbar;
%     imagesc(p_save_bp4(:, :, it) - p_save_bp1(:, :, it), [-clim clim]); colorbar;
end


% plot migration image
scale = 0.1;
cut = 50;
clim = 0.1;
i = 100;

figure(3);
subplot(nc, 1, 1);
img = mig1(cut:end, :);
img = del2(img);
imagesc(img, [-clim, clim]); colorbar;
figure(4);
plot(img(:, i), 'k', 'linewidth', 2); hold on;

figure(3);
subplot(nc, 1, 2);
img = mig2(cut:end, :);
img = del2(img);
imagesc(img, [-clim, clim]); colorbar;
figure(4);
plot(img(:, i), 'r--', 'linewidth', 2); hold on;

figure(3);
subplot(nc, 1, 3);
img = mig3(cut:end, :);
img = del2(img);
imagesc(img, [-clim, clim]); colorbar;
figure(4);
plot(img(:, i), 'b--', 'linewidth', 2); hold on;

figure(3);
subplot(nc, 1, 4);
img = mig4(cut:end, :);
img = del2(img);
imagesc(img, [-clim, clim]); colorbar;
figure(4);
plot(img(:, i), 'y--', 'linewidth', 2); hold on;

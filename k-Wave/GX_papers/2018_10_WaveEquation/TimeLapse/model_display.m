clear;
clc

% Display models at different times


load('conversion_bs03j_V1.mat');

it_list = [1; 21; 41];
ni = length(it_list);

% Grid parameters
dx = 0.15;
dz = 0.15;
x0 = -19.95;
z0 = 1621.45;
Nx = 467;
Nz = 434;

x_src = 0;
z_src = 1658.05;
x_rec = [30; 30; 30; 30];
z_rec = [1635.1; 1650.1; 1658.05; 1680.1];

x_grid = x0 : dx : (x0 + (Nx-1) * dx);
z_grid = z0 : dz : (z0 + (Nz-1) * dz);

Q_max = 1000;

vp0  = cell2mat(seisProps.vp(1));
vs0  = cell2mat(seisProps.vs(1));
rho0 = cell2mat(seisProps.rho(1));
Qi0  = cell2mat(seisProps.atten(1));
Q0 = 1 ./ Qi0;
Q0(isinf(Q0)) = Q_max;

% Plot absolute value
figure(1);
for i = 1 : length(it_list)
    it = it_list(i);
    vp  = cell2mat(seisProps.vp(it));
    vs  = cell2mat(seisProps.vs(it));
    rho = cell2mat(seisProps.rho(it));
    Qi  = cell2mat(seisProps.atten(it));
    Q = 1 ./ Qi;
    Q(isinf(Q)) = Q_max;
    subplot(ni, 4, (i-1)*4+1);
    imagesc(vp); colorbar;
    subplot(ni, 4, (i-1)*4+2);
    imagesc(vs); colorbar;
    subplot(ni, 4, (i-1)*4+3);
    imagesc(rho); colorbar;
    subplot(ni, 4, (i-1)*4+4);
    imagesc(Q); colorbar;
end

% Plot variation
figure(2);
for i = 1 : length(it_list)
    it = it_list(i);
    vp  = cell2mat(seisProps.vp(it));
    vs  = cell2mat(seisProps.vs(it));
    rho = cell2mat(seisProps.rho(it));
    Qi  = cell2mat(seisProps.atten(it));
    Q = 1 ./ Qi;
    Q(isinf(Q)) = Q_max;
    subplot(ni, 4, (i-1)*4+1);
    imagesc(vp - vp0); colorbar;
    subplot(ni, 4, (i-1)*4+2);
    imagesc(vs - vs0); colorbar;
    subplot(ni, 4, (i-1)*4+3);
    imagesc(rho - rho0); colorbar;
    subplot(ni, 4, (i-1)*4+4);
    imagesc(Q - Q0); colorbar;
end

figure(3);
vp  = cell2mat(seisProps.vp(end));
Qi  = cell2mat(seisProps.atten(end));
Q = 1 ./ Qi;
Q(isinf(Q)) = Q_max;
subplot(2, 1, 1);
imagesc(x_grid, z_grid, vp); colorbar; hold on;
plot(x_src, z_src, 'p', 'markersize', 10);
plot(x_rec, z_rec, 'v', 'markersize', 10);
subplot(2, 1, 2);
imagesc(x_grid, z_grid, Q); colorbar; hold on;
plot(x_src, z_src, 'p', 'markersize', 10);
plot(x_rec, z_rec, 'v', 'markersize', 10);



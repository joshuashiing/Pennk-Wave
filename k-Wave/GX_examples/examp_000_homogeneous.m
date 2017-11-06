% GX Test following k-Wave Manual 3.1

clear;
clc

Nx = 128;
Ny = 256;
dx = 0.5;
dy = 0.5;
% kgrid = makeGrid(Nx, dx, Ny, dy); % Old version
kgrid = kWaveGrid(Nx, dx, Ny, dy);

medium.sound_speed = 2200 * ones(Nx, Ny);
% medium.sound_speed(1:50, :) = 1800;
medium.density = 2000 * ones(Nx, Ny);

disc_x_pos = 65;
disc_y_pos = 60;
disc_radius = 1;
disc_mag = 10;
% source.p0 = disc_mag * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

% source.p0 = zeros(Nx, Ny);
% source.p0(disc_x_pos, disc_y_pos) = disc_mag;
source.p_mask = zeros(Nx, Ny);
source.p_mask(disc_x_pos, disc_y_pos) = 1;
source.p = [0 0 0 0 0 1 0 0 0] * 100;
source.p_mode = 'dirichlet';

% sensor_radius = 2.5e-3;
% num_sensor_points = 50;
% sensor.mask = makeCartCircle(sensor_radius, num_sensor_points);

nx_rec = 65 : 1 : 84;
ny_rec = 200;
sensor.mask = zeros(Nx, Ny);
sensor.mask(nx_rec, ny_rec) = 1;



sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor);
figure;
plot(kgrid.t_array, sensor_data(1, :));
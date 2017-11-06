% GX Test following k-Wave Manual 3.1

clear;
clc

Nx = 128;
Ny = 256;
dx = 50e-6;
dy = 50e-6;
% kgrid = makeGrid(Nx, dx, Ny, dy); % Old version
kgrid = kWaveGrid(Nx, dx, Ny, dy);

medium.sound_speed = 1500 * ones(Nx, Ny);
% medium.sound_speed(1:50, :) = 1800;
medium.density = 1040 * ones(Nx, Ny);

disc_x_pos = 75;
disc_y_pos = 120;
disc_radius = 8;
disc_mag = 3;
source.p0 = disc_mag * makeDisc(Nx, Ny, disc_x_pos, disc_y_pos, disc_radius);

sensor_radius = 2.5e-3;
num_sensor_points = 50;
sensor.mask = makeCartCircle(sensor_radius, num_sensor_points);



sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor);
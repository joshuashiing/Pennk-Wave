function [d, t_axis] = heter_simu(Nx, Ny, dx, dy, model_f0_m, model_f0, ...
    model_c0_m, model_Q, model_density, stf, x_src, y_src, x_rec, y_rec, ...
    dt, t_max, mod_mech, t_snap, mat_name, args)

% Grid
kgrid = kWaveGrid(Nx, dx, Ny, dy);
kgrid.t_array = 0 : dt: t_max;

% Medium
model_gamma= atan(1 ./ model_Q) / pi;
model_c0   = model_c0_m .* (model_f0 ./ model_f0_m) .^ model_gamma;

medium.f0 = model_f0;
medium.c0 = model_c0;
medium.Q  = model_Q;
medium.density = model_density;
medium.sound_speed = medium.c0;
medium.mod_mech = mod_mech;

% Source
[nx_src, ny_src] = close_grid_2d(kgrid, x_src, y_src);
source.p_mask = zeros(Nx, Ny);
source.p_mask(nx_src, ny_src) = 1;
source.p = stf;

% Sensor
[nx_rec, ny_rec] = close_grid_2d(kgrid, x_rec, y_rec);
sensor.mask = zeros(Nx, Ny);
sensor.mask(nx_rec, ny_rec) = 1;
if ~strcmp(mat_name, '') && t_snap >= 0
    sensor.t_snap = t_snap;
    sensor.snap_name = mat_name;
end

% =========================================================================
% Forward modeling
% =========================================================================

if exist('args', 'var')
    d = kspaceFirstOrder2D(kgrid, medium, source, sensor, args{:});
else
    d = kspaceFirstOrder2D(kgrid, medium, source, sensor);
end
% d = kspaceFirstOrder2D(kgrid, medium, source, sensor, 'PMLInside', false, 'PlotSim', false);
t_axis = kgrid.t_array;

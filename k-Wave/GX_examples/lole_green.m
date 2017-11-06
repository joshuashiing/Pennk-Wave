function [green_fun_f] = lole_green(kgrid, medium, source, sensor, f_vec)

% Compute the Green's function for lossless model in frequency domain
% Follow Carcione Book P157

w_vec = 2 * pi * f_vec;

% Source-Receiver Geometry
[nx_src, ny_src] = find(source.p_mask);
[nx_rec, ny_rec] = find(sensor.mask);

x_src = kgrid.x_vec(nx_src);
y_src = kgrid.y_vec(ny_src);
x_rec = kgrid.x_vec(nx_rec);
y_rec = kgrid.y_vec(ny_rec);

r = sqrt((x_rec - x_src)^2 + (y_rec - y_src)^2);

% Homogeneous medium
c0 = mean(medium.c0(:));

% % time domain solution
% t = kgrid.t_array;
% gf = (t.^2 - r^2 / c0^2).^(-1/2) .* heaviside(t - r/c0) * 2;
% plot(t, gf);

green_fun_vec = 1i / 4 * besselh(0, 2, w_vec / c0 * r);
green_fun_vec = green_fun_vec .* w_vec.^2 / (c0^2);
green_fun_f = ifftshift(green_fun_vec);
% green_fun_f = ifftshift(-1i / 4 * besselh(0, 2, w_vec * r / c0));
% green_fun_f = ifftshift(-1i / 4 * besselh(0, 2, w_vec * r / c0) * 1i .* w_vec);
% GXNOTE: where this iw comes from?

green_fun_f(1) = 0;   % Fix the NaN problem at f = 0 Hz
green_fun_f = recons_conj(green_fun_f);
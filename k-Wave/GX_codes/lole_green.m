function [green_fun_f] = lole_green(kgrid, medium, source, sensor, f_vec)

% Compute the Green's function for lossless model in frequency domain
% Follow Carcione Book P157

w_vec = 2 * pi * f_vec;
w_vec = reshape(w_vec, 1, length(w_vec));

% Source-Receiver Geometry
[nx_src, ny_src] = find(source.p_mask);
[nx_rec, ny_rec] = find(sensor.mask);

x_src = kgrid.x_vec(nx_src);
y_src = kgrid.y_vec(ny_src);
x_rec = kgrid.x_vec(nx_rec);
y_rec = kgrid.y_vec(ny_rec);

r = sqrt((x_rec - x_src).^2 + (y_rec - y_src).^2);
n_rec = length(r);
green_fun_f = zeros(n_rec, length(w_vec));

% Homogeneous medium
c0 = mean(medium.c0(:));

for i_rec = 1 : n_rec
    ri = r(i_rec);
    green_fun_vec_i = 1i / 4 * besselh(0, 2, w_vec / c0 * ri);
    green_fun_vec_i = green_fun_vec_i .* w_vec.^2 / (c0^2);
    green_fun_f_i = ifftshift(green_fun_vec_i);
    green_fun_f_i(1) = 0;   % Fix the NaN problem at f = 0 Hz
    green_fun_f_i = recons_conj(green_fun_f_i);
    green_fun_f(i_rec, :) = green_fun_f_i;
end
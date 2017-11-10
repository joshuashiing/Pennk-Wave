function [green_fun_f] = kjar_green(kgrid, medium, source, sensor, f_vec)

% Compute the Green's function for Kjartansson's Model in frequency domain
% Follow Carcione Book P158

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
f0 = mean(medium.f0(:));
c0 = mean(medium.c0(:));
rho = mean(medium.density(:));
Q = mean(medium.Q(:));

w0 = 2 * pi * f0;
gamma = atan(1 / Q) / pi;
M0 = rho * c0^2 * (cos(pi * gamma / 2))^2;
beta = 2 - 2 * gamma;
b = M0 / rho * w0^(-2 * gamma);
Omega = -1i * (1i * w_vec).^(beta / 2);

for i_rec = 1 : n_rec
    ri = r(i_rec);
    green_fun_vec_i = 1i / 4 * besselh(0, 2, Omega / sqrt(b) * ri);
    green_fun_vec_i = -(1i*w_vec).^beta / b .* green_fun_vec_i;
    green_fun_f_i = ifftshift(green_fun_vec_i);
    green_fun_f_i(1) = 0;   % Fix the NaN problem at f = 0 Hz
    green_fun_f_i = recons_conj(green_fun_f_i);
    green_fun_f(i_rec, :) = green_fun_f_i;
end
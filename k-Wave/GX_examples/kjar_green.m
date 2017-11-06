function [green_fun_f] = kjar_green(kgrid, medium, source, sensor, f_vec)

% Compute the Green's function for Kjartansson's Model in frequency domain
% Follow Carcione Book P158

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

green_fun_vec = 1i / 4 * besselh(0, 2, Omega / sqrt(b) * r);
green_fun_vec = -(1i*w_vec).^beta / b .* green_fun_vec;
green_fun_f = ifftshift(green_fun_vec);

% green_fun_f = ifftshift(-1i / 4 * besselh(0, 2, Omega * r / sqrt(b)));
% green_fun_f = ifftshift(-1i / 4 * besselh(0, 2, Omega * r / sqrt(b)) * 1i .* w_vec);
% GXNOTE: where this iw comes from?

% GXNOTE: play with the coefficient later
green_fun_f(1) = 0;   % Fix the NaN problem at f = 0 Hz
green_fun_f = recons_conj(green_fun_f);
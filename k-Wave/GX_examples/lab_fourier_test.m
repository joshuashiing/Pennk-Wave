% DESCRIPTION:
%     Trying to understand the Fourier domain manipulation, especially when
%     dealing with complex wavenumber.
%
% Author: Guangchi Xing
% Date: 10/26/2017
clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up the characteristic parameters of exploration seismology
% c0  |  2  - 4 km/s
% f   |  15 - 40 Hz
% rho |  2  - 2.5 g/cm^3
% Q   |  20 - 100

nf  = 100;
c0  = 2500;
% f   = linspace(15, 100, nf);	% Consider a frequency range
f   = 70;                       % Consider a specific frequency
rho = 2000;
Q   = 10;
f0  = 60;   % Reference frequencys

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the relevant parameters (Kjartansson's model)
gamma = 1 / pi * atan(1 / Q);
M0 = rho * c0^2 * (cos(pi*gamma/2))^2;
w0 = 2 * pi * f0;
w = 2 * pi * f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate the dispersion relations
M = M0 * (1i * w / w0) .^ (2 * gamma);
k = sqrt(rho / M0) * (1i * w / w0) .^ (-gamma) .* w;    % Complex wavenum.
k_r = real(k);
k_i = imag(k);
vc = sqrt(M0 / rho) * (1i * w / w0) .^ gamma;           % Complex velocity
cp = c0 * (w / w0) .^ gamma;                            % Phase velocity
alpha = w ./ cp * tan(pi * gamma / 2);                  % Attenuation Fac.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Consider a 1-D problem
Nx = 1024;           % number of grid points in the x direction
dx = 0.5;           % grid point spacing in the x direction [m]

dt = 1e-5;        % Time interval [s]
t_max = 0.06;        % Simulation end time [s]

kgrid = kWaveGrid(Nx, dx);
kgrid.t_array = 0 : dt: t_max;

% Consider the wavefield corresponding to specific frequency
u0 = 1.0;
t = 20;
u = u0 * exp(1i * (-k * kgrid.x_vec + w * t));
max_amp = max(abs(u(:)));

u = real(u);
u_k = fft(u);
u_k_vec = fftshift(u_k);
nabla_exp = -0.7;
nabla_u_k_vec = u_k_vec .* kgrid.k.^nabla_exp;
nabla_u_k = ifftshift(nabla_u_k_vec);
nabla_u_k(1) = 0;
nabla_u = ifft(nabla_u_k);
% nabla_u = ifft(ifftshift(nabla_u_k_vec));
nabla_u_theo = u0 * k.^nabla_exp * exp(1i * (-k * kgrid.x_vec + w * t));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Test plot for complex value
% val = vc;
% plot(f, real(val), 'k', 'linewidth', 2); hold on;
% plot(f, imag(val), 'r', 'linewidth', 2);
% legend('Real', 'Imaginary');

% % Test plot for real value
% val = cp;
% plot(f, val, 'k', 'linewidth', 2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot 1-D wavefield snapshot
val = real(nabla_u);
plot(kgrid.x_vec, val, 'k', 'linewidth', 2);
% xlim([min(kgrid.x_vec) max(kgrid.x_vec)]);
% ylim([-max_amp max_amp]);

% Compare 1-D wavefield snapshot
val1 = real(nabla_u);
val2 = real(nabla_u_theo);
plot(kgrid.x_vec, val1, 'k', 'linewidth', 2); hold on;
plot(kgrid.x_vec, val2, 'r--', 'linewidth', 2);



% % Plot 1-D wavefield animation
% for t_index = 1 : length(kgrid.t_array)
%     t = kgrid.t_array(t_index);
%     u = u0 * exp(1i * (-k * kgrid.x_vec + w * t));
%     plot(kgrid.x_vec, real(u), 'k', 'linewidth', 2);
%     xlim([min(kgrid.x_vec) max(kgrid.x_vec)]);
%     ylim([-max_amp max_amp]);
%     pause(0.001);
% end

% % Plot wavenumber domain
% val = real(u_k_vec);
% plot(kgrid.kx_vec, val, 'k', 'linewidth', 2);

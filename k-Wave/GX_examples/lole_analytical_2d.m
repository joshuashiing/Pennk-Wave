function [sensor_data] = lole_analytical_2d(kgrid, medium, source, sensor)

% Initialize the frequency axis
f_vec = make_dual_dim(length(kgrid.t_array) , kgrid.dt);

% Compute the Green's function in frequency domain
green_fun_f = lole_green(kgrid, medium, source, sensor, f_vec);

% Source wavelet
src_wavelet = zeros(size(kgrid.t_array));
src_wavelet(:, 1 : size(source.p, 2)) = source.p;

% Convolve the source wavelet with the Green's function
src_wavelet_f = fft(src_wavelet);

src_wavelet_f = src_wavelet_f(:);
green_fun_f = green_fun_f(:);

% green_fun_f = ifftshift(fftshift(green_fun_f) .* exp(1i * pi / 4));
% green_fun_f = recons_conj(green_fun_f);

sensor_data_f = src_wavelet_f .* green_fun_f;
sensor_data = ifft(sensor_data_f);
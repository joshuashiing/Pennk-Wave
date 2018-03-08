function [sensor_data] = lole_analytical_2d(kgrid, medium, source, sensor)

% Initialize the frequency axis
f_vec = make_dual_dim(length(kgrid.t_array) , kgrid.dt);

% Compute the Green's function in frequency domain
green_fun_f = lole_green(kgrid, medium, source, sensor, f_vec) / 4; % GXNOTE: figure out coefficient
% green_fun_f = lole_green(kgrid, medium, source, sensor, f_vec);
n_rec = size(green_fun_f, 1);
sensor_data = zeros(n_rec, length(kgrid.t_array));

% Source wavelet
src_wavelet = zeros(size(kgrid.t_array));
src_wavelet(:, 1 : size(source.p, 2)) = source.p;

% Convolve the source wavelet with the Green's function
src_wavelet_f = fft(src_wavelet);
src_wavelet_f = src_wavelet_f(:);

for i_rec = 1 : n_rec
    green_fun_f_i = green_fun_f(i_rec, :);
    green_fun_f_i = green_fun_f_i(:);
    sensor_data_f_i = src_wavelet_f .* green_fun_f_i;
    sensor_data_i = ifft(sensor_data_f_i);
    sensor_data(i_rec, :) = transpose(sensor_data_i(:));
end

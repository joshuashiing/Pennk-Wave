function [f_vec, d_vec] = comp_spec(t_array, d)
% Compute the spectrum of a time series d associated with a time axis
dt = (t_array(end) - t_array(1)) / (length(t_array) - 1);
f_vec = make_dual_dim(length(t_array) , dt);
n = size(d, 1);
d_vec = zeros(n, length(f_vec));
for i = 1 : n
    d_f_i = fft(d(i, :));
    d_vec_i = fftshift(d_f_i);
    d_vec(i, :) = d_vec_i;
end
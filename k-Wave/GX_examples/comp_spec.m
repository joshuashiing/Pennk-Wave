function [f_vec, d_vec] = comp_spec(t_array, d)
% Compute the spectrum of a time series d associated with a time axis
dt = (t_array(end) - t_array(1)) / (length(t_array) - 1);
f_vec = make_dual_dim(length(t_array) , dt);
d_f = fft(d);
d_vec = fftshift(d_f);
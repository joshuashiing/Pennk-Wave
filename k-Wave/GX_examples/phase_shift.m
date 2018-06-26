function d_ps = phase_shift(t_array, d, phi)
dt = (t_array(end) - t_array(1)) / (length(t_array) - 1);
f_vec = make_dual_dim(length(t_array) , dt);
d_f = fft(d);
n_f = length(f_vec);
d_f_p = d_f(1 : ceil(n_f / 2));
d_ps_f_p = d_f_p;
d_ps_f_p(2:end) = d_f_p(2:end) * exp(1i * phi);
d_ps_f_n = transpose(flipud(d_ps_f_p(2:end)'));
d_ps_f_vec = [d_ps_f_n, d_ps_f_p];
d_ps_f = ifftshift(d_ps_f_vec);
d_ps = real(ifft(d_ps_f));

% 
% % d_vec_f = fftshift(fft(d_f));
% 
% 
% 
% n = size(d, 1);
% d_vec = zeros(n, length(f_vec));
% for i = 1 : n
%     d_f_i = fft(d(i, :));
%     d_vec_i = fftshift(d_f_i);
%     d_vec(i, :) = d_vec_i;
% end
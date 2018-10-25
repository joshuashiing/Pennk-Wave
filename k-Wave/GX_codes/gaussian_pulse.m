function [t, f] = gaussian_pulse(A, td, v, dt, nt)

% Function to generate a Gaussian pulse
% A  - amplitude
% td - time delay
% v  - variance
% dt - time step
% nt - number of time steps

t = (0 : (nt - 1)) * dt;
f = A * exp(- (-(t-td) / v) .^ 2);

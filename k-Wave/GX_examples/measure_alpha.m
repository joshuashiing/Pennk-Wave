function [alpha_meas] = measure_alpha(d1_vec, d2_vec, r1, r2, d)
A1_vec = abs(d1_vec) * sqrt(r1);
A2_vec = abs(d2_vec) * sqrt(r2);
alpha_meas = -20 * log(A2_vec ./ A1_vec) ./ log(10) * d;


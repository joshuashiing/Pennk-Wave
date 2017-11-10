function cp_meas = measure_cp(f_vec, d1_vec, d2_vec, d);
w_vec = transpose(f_vec(:)) * 2 * pi;
cp_meas = w_vec * d ./ (angle(d1_vec) - angle(d2_vec));
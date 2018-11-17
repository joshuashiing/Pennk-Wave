% =========================================================================
% Test the seismogram spectrums
% Date 11/13/2018
% Author: Guangchi Xing


clear;
clc

% data_file = 'Data_example_05/CSG_ac_tl_003.mat';
% load(data_file);

data_file1 = 'Data_example_05/CSG_010.mat';
data_file2 = 'Data_example_05/CSG_tl_010.mat';
data1 = load(data_file1);
data2 = load(data_file2);
d = data1.d - data2.d;
t_axis = data1.t_axis;

d_vec = zeros(size(d));

n = size(d, 1);
for i = 1 : n
    [f_vec, d_vec(i, :)] = comp_spec(t_axis, d(i, :));
end
% 
% for i = 1 : n
%     plot(f_vec, abs(d_vec(i, :))); hold on;
% end

spec = sum(abs(d_vec), 1);
% db = 20 * log10(spec);

plot(f_vec, spec);
% plot(f_vec, db);
% xlim([-100, 100]);
% xlim([-500, 500]);
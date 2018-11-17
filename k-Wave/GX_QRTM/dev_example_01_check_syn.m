clear;
clc

% dfile1 = 'Data_example_01/CSG_ac_015.mat';
% dfile2 = 'Data_example_01/CSG_ac_sm_015.mat';
% dfile3 = 'Data_example_01/CSG_ac_tl_015.mat';

dfile1 = 'Data_example_01/CSG_015.mat';
dfile2 = 'Data_example_01/CSG_tl_015.mat';
% dfile3 = 'Data_example_01/CSG_tl_015.mat';

data1 = load(dfile1);
data2 = load(dfile2);
% data3 = load(dfile3);

d1 = data1.d';
d2 = data2.d';
% d3 = data3.d';

clim = max(d1(:)) * 0.01;

figure(1);
subplot(121); imagesc(d1, [-clim, clim]);
subplot(122); imagesc(d1 - d2, [-clim, clim]);
% subplot(133); imagesc(d1 - d3, [-clim, clim]);

figure(2);
i = 100;
plot(data1.t_axis, d1(:, i), 'k', 'linewidth', 2); hold on;
plot(data2.t_axis, d2(:, i), 'b--', 'linewidth', 2);
% plot(data3.t_axis, d3(:, i), 'r--', 'linewidth', 2);
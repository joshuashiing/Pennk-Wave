clear;
clc

xrange = [0.15 0.25];
yrange = [-0.06 0.06];
casename = 'case01';
plot_HS(casename, 1, 1, xrange, [-0.06 0.06]);
casename = 'case02';
plot_HS(casename, 1, 2, xrange, [-0.03 0.03]);
casename = 'case03';
plot_HS(casename, 1, 3, xrange, [-6e-3 6e-3]);
casename = 'case04';
plot_HS(casename, 1, 4, xrange, [-3e-4 3e-4]);

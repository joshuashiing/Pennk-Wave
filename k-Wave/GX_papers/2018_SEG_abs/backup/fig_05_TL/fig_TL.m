clear;
clc

load('case01_a');
p_lim = max(abs(p(:))) / 10;

plot_TL_snap('case01_a', 1, 331, 8, p_lim);
plot_TL_snap('case01_b', 1, 332, 8, p_lim);
plot_TL_snap_dif('case01_a', 'case01_b', 1, 333, 8, p_lim);

plot_TL_snap('case02_a', 1, 334, 8, p_lim);
plot_TL_snap('case02_b', 1, 335, 8, p_lim);
plot_TL_snap_dif('case02_a', 'case02_b', 1, 336, 8, p_lim);

plot_TL_snap('case03_a', 1, 337, 8, p_lim);
plot_TL_snap('case03_b', 1, 338, 8, p_lim);
plot_TL_snap_dif('case03_a', 'case03_b', 1, 339, 8, p_lim);

plot_cs_dif('case01_a', 'case01_b', 100, 2, 311);
plot_cs_dif('case02_a', 'case02_b', 100, 2, 312);
plot_cs_dif('case03_a', 'case03_b', 100, 2, 313);

function plot_cs_dif(casename1, casename2, i_col, figureid, subid)
load(casename1);
cs1 = p(:, i_col);
load(casename2);
cs2 = p(:, i_col);
figure(figureid);
subplot(subid);
plot(cs2, 'b', 'linewidth', 2); hold on;
plot(cs1, 'r--', 'linewidth', 2);
plot(cs1 - cs2, 'k:', 'linewidth', 1);
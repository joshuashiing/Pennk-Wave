function plot_HS(casename, figureid, subfid, xrange, yrange)
load(casename)
figure(figureid);
subplot(2, 2, subfid)
plot(t_axis, d2, 'b', 'linewidth', 2); hold on;
plot(t_axis, d3, 'r--', 'linewidth', 2);
xlim(xrange);
ylim(yrange);
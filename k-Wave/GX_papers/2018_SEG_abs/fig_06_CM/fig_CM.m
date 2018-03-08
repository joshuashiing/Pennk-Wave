clear;

n_PML = 20;
load('BP_cut_model_gas_chimney');
c = velp;
h = 12.5;
scale = 2e2;
clim = [min(c(:)) max(c(:))];



load('case01');
p1 = p; d1 = d;
pmax = max(abs(p(:)));
plim = [-pmax pmax] / 10;
plot_model_snap(c, zeros(size(p)), n_PML, h, 0, clim, 1, 211);
plot_model_snap(Qp, zeros(size(p)), n_PML, h, 0, [min(Qp(:)) max(Qp(:))], 1, 212);
plot_model_snap(c, p1, n_PML, h, scale, clim, 2, 411);
plot_gather(t_axis, d1, h, plim, 3, 221);

load('case02');
p2 = p; d2 = d;
plot_model_snap(c, p2, n_PML, h, scale, clim, 2, 412);
plot_gather(t_axis, d2, h, plim, 3, 222);

load('case03');
p3 = p; d3 = d;
plot_model_snap(c, p3, n_PML, h, scale, clim, 2, 413);
plot_gather(t_axis, d3, h, plim, 3, 223);

dp = p2 - p3;
plot_model_snap(c, dp, n_PML, h, scale, clim, 2, 414);
plot_gather(t_axis, d2 - d3, h, plim, 3, 224);

irec = 241; % 3000 m / 512.5m

figure(4);
subplot(211);
plot(t_axis, d1(irec, :), 'k', 'linewidth', 2); hold on;
plot(t_axis, d2(irec, :), 'r--', 'linewidth', 2);
xlim([0.6 2]);

subplot(212);
plot(t_axis, d3(irec, :), 'k', 'linewidth', 2); hold on;
plot(t_axis, d2(irec, :), 'r--', 'linewidth', 2);
plot(t_axis, d2(irec, :) - d3(irec, :), 'b:', 'linewidth', 1);
xlim([0.6 2]);

function plot_TL_snap(casename, figureid, subid, h, p_lim)
load(casename);

[Nx, Ny] = size(p);
x_axis = (0 : (Nx - 1)) * h;
y_axis = (0 : (Ny - 1)) * h;

figure(figureid)
subplot(subid)
imagesc(y_axis, x_axis, p, [-p_lim p_lim]);
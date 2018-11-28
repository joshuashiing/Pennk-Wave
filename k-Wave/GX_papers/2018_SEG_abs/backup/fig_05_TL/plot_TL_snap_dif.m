function plot_TL_snap_dif(casename1, casename2, figureid, subid, h, p_lim)
load(casename1);
[Nx, Ny] = size(p);
x_axis = (0 : (Nx - 1)) * h;
y_axis = (0 : (Ny - 1)) * h;
p1 = p;
load(casename2);
p2 = p;
figure(figureid)
subplot(subid)
imagesc(y_axis, x_axis, p1 - p2, [-p_lim p_lim]);